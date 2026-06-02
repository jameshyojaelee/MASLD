#!/usr/bin/env python3
"""
13e_multisp_gene_markers.py — gene-level marker producer for MultiSP domains.

WHY THIS EXISTS
---------------
The unified spatial atlas writer ``17a_unified_spatial_integration.py``
(``build_multisp_gene_flags``) expects a *gene-level* MultiSP file
``multisp/domain_marker_genes.csv`` (or ``multisp/validation_gene_flags.csv``),
maps gene -> domain, and flags ``multisp_in_disease_domain``. Neither file was
ever produced: the only MultiSP output, ``multisp/multisp_domains.csv``, is
PER-SPOT (barcode -> multisp_domain / rna_only_domain / sample_id / condition),
with no gene column. So MultiSP contributed 0 genes to the atlas.

The coverage review also found that 13a's "gene activity" is actually 500-bp
genome TILES (imputed chromatin), NOT gene-level activity — so we do NOT use
that. Markers here are derived from the spatial *RNA expression* (the genuine
per-gene signal in the MultiSP-labelled object).

WHAT THIS PRODUCES
------------------
``Analysis/Spatial/results/multisp/domain_marker_genes.csv`` with columns:
  * ``human_symbol``       — gene symbol (the atlas key; identify_gene_column
                              finds this first).
  * ``domain``             — the MultiSP domain the gene most strongly marks.
  * ``is_disease_domain``  — bool; True (every emitted row is restricted to a
                              disease-enriched domain, see below).
  * ``score``              — rank_genes_groups Wilcoxon z-score for that gene
                              in its top domain (higher = stronger marker).
  * ``logfoldchange``      — rank_genes_groups log2FC (descriptive companion).
  * ``steatotic_frac``     — Steatotic-spot fraction of the gene's domain.

Only markers of disease-enriched domains are written. This keeps
``build_multisp_gene_flags``'s current logic correct: it sets
``multisp_in_disease_domain = domain.notna()``, i.e. presence of any row means
"in a disease domain" — true by construction here.

DISEASE-ENRICHED DOMAIN DEFINITION
----------------------------------
A MultiSP domain is "disease-enriched" if its Steatotic-spot fraction exceeds
the cohort baseline by a margin: ``frac_steatotic(domain) >= baseline + MARGIN``
(MARGIN = 0.10). The cohort baseline is the overall Steatotic-spot fraction
(~0.41 in GSE192741). On the current object this selects the domains whose
steatotic fraction is 0.57-0.99 (well separated from the <0.31 healthy-leaning
domains) — the threshold is not knife-edge.

RIGOR CAVEAT (donor n)
----------------------
The spatial cohort is ~5 donors (GSE192741: sample_id JBO014/015/018/019/022;
donor == sample_id). Marker genes are computed by pooling SPOTS across donors,
so these statistics are pseudoreplicated and DESCRIPTIVE — they rank genes
within MultiSP domains, they are NOT donor-level significance claims. The atlas
consumes them only as a binary "in disease domain" flag, for which a descriptive
spot-pooled ranking is adequate. (All 5 donors contribute to nearly every
domain, so the disease-enrichment is spatial structure, not a single-donor
artifact, but the gene p-values themselves remain descriptive.) This mirrors the
spot/cell pseudoreplication caveat documented in spatial_stats.py.

Env: spatial (scanpy, anndata). Reads spatial_with_multisp_domains.h5ad
(falls back to merging multisp_domains.csv onto merged_spatial.h5ad by barcode).
"""

import pathlib
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import RESULTS_DIR, print_header, print_step  # noqa: E402
from spatial_stats import ensure_lognorm  # noqa: E402

# ── Config ──────────────────────────────────────────────────────────────────
MULTISP_DIR = RESULTS_DIR / "multisp"
PREPROC_DIR = RESULTS_DIR / "preprocessed"

DOMAIN_COL = "multisp_domain"
DONOR_COL = "sample_id"        # donor == sample_id for GSE192741
CONDITION_COL = "condition"
DISEASE_LABEL = "Steatotic"    # the disease arm in this cohort

# A domain is disease-enriched if its Steatotic-spot fraction exceeds the cohort
# baseline by at least this margin. 0.10 cleanly separates the 0.57-0.99 disease
# domains from the <0.31 healthy-leaning ones on the current object.
DISEASE_FRAC_MARGIN = 0.10

# rank_genes_groups marker thresholds (descriptive; gate which genes are emitted)
MIN_LOGFC = 0.25               # log2FC floor for a gene to count as a marker
MAX_PADJ = 0.05                # within-domain Wilcoxon BH padj ceiling


# ── Data loading ────────────────────────────────────────────────────────────
def load_domain_labelled_adata():
    """Load the spatial object carrying MultiSP domain labels.

    Primary: ``multisp/spatial_with_multisp_domains.h5ad`` (already has
    ``multisp_domain`` in .obs). Fallback: merge ``multisp/multisp_domains.csv``
    onto ``preprocessed/merged_spatial.h5ad`` by barcode.
    """
    primary = MULTISP_DIR / "spatial_with_multisp_domains.h5ad"
    if primary.exists():
        print_step(f"Loading domain-labelled object: {primary.name}")
        adata = sc.read_h5ad(primary)
        if DOMAIN_COL not in adata.obs.columns:
            raise KeyError(f"{primary} missing '{DOMAIN_COL}' in .obs")
        return adata

    print_step("Primary h5ad absent — merging multisp_domains.csv onto "
               "merged_spatial.h5ad by barcode")
    csv = MULTISP_DIR / "multisp_domains.csv"
    merged = PREPROC_DIR / "merged_spatial.h5ad"
    if not csv.exists() or not merged.exists():
        raise FileNotFoundError(
            f"Need either {primary} or both {csv} and {merged}")
    dom = pd.read_csv(csv, index_col=0)
    adata = sc.read_h5ad(merged)
    common = adata.obs_names.intersection(dom.index)
    if len(common) == 0:
        raise ValueError("No barcode overlap between merged_spatial and "
                         "multisp_domains.csv")
    adata = adata[common].copy()
    dom = dom.loc[common]
    adata.obs[DOMAIN_COL] = dom[DOMAIN_COL].values
    # carry condition / sample_id from the CSV if missing on the object
    for c in (CONDITION_COL, DONOR_COL):
        if c not in adata.obs.columns and c in dom.columns:
            adata.obs[c] = dom[c].values
    print_step(f"Merged labels onto {adata.n_obs} spots "
               f"({adata.n_obs}/{len(dom)} barcodes matched)")
    return adata


# ── Disease-enriched domain identification ──────────────────────────────────
def identify_disease_domains(adata):
    """Return (disease_domains, summary_df).

    A domain is disease-enriched if its Steatotic-spot fraction is at least
    ``DISEASE_FRAC_MARGIN`` above the cohort baseline Steatotic fraction.
    """
    obs = adata.obs
    is_disease_spot = (obs[CONDITION_COL].astype(str) == DISEASE_LABEL)
    baseline = float(is_disease_spot.mean())

    rows = []
    for dom, idx in obs.groupby(DOMAIN_COL, observed=True).groups.items():
        sub = obs.loc[idx]
        frac = float((sub[CONDITION_COL].astype(str) == DISEASE_LABEL).mean())
        rows.append({
            "domain": dom,
            "n_spots": int(len(sub)),
            "steatotic_frac": frac,
            "n_donors": int(sub[DONOR_COL].nunique()) if DONOR_COL in sub else np.nan,
            "is_disease_domain": frac >= baseline + DISEASE_FRAC_MARGIN,
        })
    summary = pd.DataFrame(rows).sort_values("steatotic_frac", ascending=False)
    summary.attrs["baseline"] = baseline

    disease_domains = summary.loc[summary["is_disease_domain"], "domain"].tolist()
    print_step(f"Cohort Steatotic baseline = {baseline:.3f}; "
               f"threshold = baseline + {DISEASE_FRAC_MARGIN} = "
               f"{baseline + DISEASE_FRAC_MARGIN:.3f}")
    print("  Per-domain Steatotic fraction:")
    for _, r in summary.iterrows():
        tag = "  <-- DISEASE" if r["is_disease_domain"] else ""
        print(f"    domain {r['domain']}: frac={r['steatotic_frac']:.3f} "
              f"(n={r['n_spots']}, donors={r['n_donors']}){tag}")
    return disease_domains, summary


# ── Marker computation ──────────────────────────────────────────────────────
def compute_domain_markers(adata, lognorm_layer):
    """rank_genes_groups (Wilcoxon) per MultiSP domain on the lognorm layer.

    Returns a tidy DataFrame: one row per (gene, domain) with score / logFC /
    padj, restricted to genes passing MIN_LOGFC and MAX_PADJ.
    """
    # rank_genes_groups reads .X; point .X at the lognorm layer for the call.
    adata.obs[DOMAIN_COL] = adata.obs[DOMAIN_COL].astype("category")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            adata, groupby=DOMAIN_COL, method="wilcoxon",
            layer=lognorm_layer, use_raw=False, tie_correct=True,
        )
    tidy = sc.get.rank_genes_groups_df(adata, group=None)
    # standardize column names produced by scanpy
    tidy = tidy.rename(columns={
        "group": "domain",
        "names": "human_symbol",
        "scores": "score",
        "logfoldchanges": "logfoldchange",
        "pvals_adj": "padj",
    })
    keep = (tidy["logfoldchange"] >= MIN_LOGFC) & (tidy["padj"] <= MAX_PADJ)
    tidy = tidy.loc[keep].copy()
    print_step(f"rank_genes_groups: {len(tidy)} (gene,domain) markers pass "
               f"logFC>={MIN_LOGFC} & padj<={MAX_PADJ}")
    return tidy


def assign_top_domain(markers, disease_domains, domain_summary):
    """Restrict to disease domains, then assign each gene to its single best
    domain (highest Wilcoxon score). Returns the per-gene marker table.
    """
    markers = markers.copy()
    markers["domain"] = markers["domain"].astype(str)
    disease_domains = [str(d) for d in disease_domains]

    dis = markers[markers["domain"].isin(disease_domains)].copy()
    if dis.empty:
        raise ValueError("No markers fell in a disease-enriched domain — "
                         "check DISEASE_FRAC_MARGIN / condition labels.")

    # One row per gene: the disease domain it most strongly marks (max score).
    dis = dis.sort_values("score", ascending=False)
    top = dis.drop_duplicates(subset="human_symbol", keep="first").copy()

    frac_map = dict(zip(domain_summary["domain"].astype(str),
                        domain_summary["steatotic_frac"]))
    top["is_disease_domain"] = True
    top["steatotic_frac"] = top["domain"].map(frac_map)

    cols = ["human_symbol", "domain", "is_disease_domain",
            "score", "logfoldchange", "steatotic_frac"]
    top = top[cols].sort_values(["domain", "score"],
                                ascending=[True, False]).reset_index(drop=True)
    print_step(f"{len(top)} unique disease-domain marker genes "
               f"across {top['domain'].nunique()} disease domains")
    return top


# ── Main ────────────────────────────────────────────────────────────────────
def main():
    print_header("13e — MultiSP gene-level domain markers")

    adata = load_domain_labelled_adata()
    print_step(f"Object: {adata.n_obs} spots x {adata.n_vars} genes; "
               f"{adata.obs[DONOR_COL].nunique() if DONOR_COL in adata.obs else '?'} donors")

    # log1p-CPM layer (the spatial object's .X may be raw counts in the fallback
    # path; ensure_lognorm builds a normalized layer from 'counts' or .X).
    lognorm_layer = ensure_lognorm(adata)
    print_step(f"Using expression layer '{lognorm_layer}' for marker tests")

    disease_domains, domain_summary = identify_disease_domains(adata)
    if not disease_domains:
        raise ValueError("No disease-enriched domains identified — "
                         f"baseline={domain_summary.attrs['baseline']:.3f}")

    markers = compute_domain_markers(adata, lognorm_layer)
    table = assign_top_domain(markers, disease_domains, domain_summary)

    out_path = MULTISP_DIR / "domain_marker_genes.csv"
    table.to_csv(out_path, index=False)
    print_step(f"WROTE {out_path}  ({len(table)} rows, "
               f"{table['human_symbol'].nunique()} unique genes)")

    # also persist the per-domain disease summary for provenance / review
    summ_path = MULTISP_DIR / "multisp_domain_disease_summary.csv"
    domain_summary.assign(
        cohort_steatotic_baseline=domain_summary.attrs["baseline"]
    ).to_csv(summ_path, index=False)
    print_step(f"WROTE {summ_path} (per-domain Steatotic fraction)")


if __name__ == "__main__":
    main()
