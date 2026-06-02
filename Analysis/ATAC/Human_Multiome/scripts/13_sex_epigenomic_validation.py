#!/usr/bin/env python3
"""
13_sex_epigenomic_validation.py — Targeted validation of sex-dimorphic
transcriptomic effects using scATAC-seq chromatin accessibility.

Tests whether female-biased MASLD DEGs (from RNA-seq Script 26) show
corresponding female-specific open chromatin in hepatocyte promoters.

Modules:
  1. Sex inference from gene activity (XIST vs Y-linked scores, k-means k=2)
  2. Promoter accessibility enrichment (Female_biased vs Male_biased gene sets)
  3. chromVAR sex-differential TF activity (pseudobulk per-donor)
  4. Gene activity concordance (RNA logFC vs ATAC logFC)

Inputs: gene_activity_matrix.h5ad, snapatac2_label_transferred.h5ad,
        sex_deg_classification.csv, chromvar_deviations.h5ad, GENCODE v49 GTF
Output: results/sex_validation/*.csv + validation_summary.json
"""

import argparse
import gzip
import json
import logging
import os
import sys
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.cluster import KMeans
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
ATAC_DIR = os.path.join(PROJECT_ROOT, "Analysis/ATAC/Human_Multiome")
GENE_ACTIVITY_H5AD = os.path.join(
    ATAC_DIR, "results/snapatac2/gene_activity_matrix.h5ad"
)
LABEL_TRANSFERRED_H5AD = os.path.join(
    ATAC_DIR, "results/label_transfer/snapatac2_label_transferred.h5ad"
)
CHROMVAR_H5AD = os.path.join(
    ATAC_DIR, "results/chromvar_v2/chromvar_deviations.h5ad"
)
HEP_PEAKS_BED = os.path.join(
    ATAC_DIR, "results/label_transfer/cell_type_peak_sets_v2/Hepatocytes_peaks.bed"
)
DISEASE_REGULONS = os.path.join(ATAC_DIR, "scenic_plus/disease_regulons.csv")
HEP_REGULONS = os.path.join(ATAC_DIR, "scenic_plus/hepatocyte_regulons.csv")
DONOR_META = os.path.join(
    PROJECT_ROOT, "data/GSE244832/metadata/donor_pairing.csv"
)
SEX_DEG_CSV = os.path.join(
    PROJECT_ROOT,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/"
    "sex_deg_classification.csv",
)
SEX_CAUSAL_CSV = os.path.join(
    PROJECT_ROOT, "RNA-seq/results/stratified_causal/sex_causal_scores.csv"
)
GENCODE_GTF = (
    "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/"
    "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
)

# Sex-linked genes for inference
XIST_GENE = "XIST"
Y_GENES = ["DDX3Y", "UTY", "KDM5D", "EIF1AY", "RPS4Y1", "USP9Y", "ZFY"]

# Promoter window half-width
PROMOTER_HALFWIDTH = 2000


# =========================================================================
# Utility: parse TSS from GENCODE GTF
# (adapted from 08_annotate_peaks_for_l8.py, extended to extract gene_id)
# =========================================================================
def parse_tss_from_gtf(gtf_path):
    """Parse gene TSS positions from GENCODE GTF.

    Returns DataFrame: gene_name, gene_id, chrom, tss, strand, gene_type
    """
    log.info("Parsing TSS from GTF: %s", gtf_path)
    records = []
    opener = gzip.open if gtf_path.endswith(".gz") else open

    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            chrom = fields[0]
            if not chrom.startswith("chr"):
                continue

            strand = fields[6]
            start = int(fields[3])  # 1-based
            end = int(fields[4])

            attrs = fields[8]
            gene_name = None
            gene_type = None
            gene_id = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]
                elif attr.startswith("gene_type"):
                    gene_type = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]
                elif attr.startswith("gene_id"):
                    gene_id = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]

            if gene_name is None:
                continue

            tss = start if strand == "+" else end
            records.append({
                "gene_name": gene_name,
                "gene_id": gene_id or "",
                "chrom": chrom,
                "tss": tss,
                "strand": strand,
                "gene_type": gene_type or "",
            })

    df = pd.DataFrame(records)
    df = df.drop_duplicates(subset="gene_name", keep="first")
    log.info("  Parsed %d gene TSS positions", len(df))
    return df


# =========================================================================
# Module 1: Sex inference from gene activity
# =========================================================================
def module1_sex_inference(ga, out_dir):
    """Infer donor sex from XIST vs Y-linked gene activity scores."""
    log.info("=" * 60)
    log.info("MODULE 1: Sex inference from gene activity")
    log.info("=" * 60)

    var_names = list(ga.var_names)

    # Extract sex-linked gene indices
    sex_genes = [XIST_GENE] + Y_GENES
    found = {g: var_names.index(g) for g in sex_genes if g in var_names}
    log.info("  Found %d/%d sex-linked genes in gene activity matrix", len(found), len(sex_genes))

    if XIST_GENE not in found:
        log.error("XIST not found in gene activity matrix -- cannot infer sex")
        return None

    y_found = [g for g in Y_GENES if g in found]
    if len(y_found) == 0:
        log.error("No Y-linked genes found -- cannot infer sex")
        return None

    # Extract per-cell values
    log.info("  Extracting sex-linked gene activity (loading into memory)...")
    X = ga.X
    if sparse.issparse(X):
        X = X.toarray()

    sex_vals = {g: X[:, idx] for g, idx in found.items()}
    sex_df = pd.DataFrame(sex_vals)
    sex_df["donor_id"] = ga.obs["donor_id"].values

    # Per-donor means
    donor_means = sex_df.groupby("donor_id").mean()
    donor_means["xist_score"] = donor_means[XIST_GENE]
    donor_means["y_score"] = donor_means[y_found].mean(axis=1)

    # k-means (k=2)
    features = donor_means[["xist_score", "y_score"]].values
    km = KMeans(n_clusters=2, random_state=42, n_init=10).fit(features)
    c0_xist = features[km.labels_ == 0, 0].mean()
    c1_xist = features[km.labels_ == 1, 0].mean()
    female_cluster = 0 if c0_xist > c1_xist else 1
    donor_means["inferred_sex"] = [
        "Female" if l == female_cluster else "Male" for l in km.labels_
    ]

    # Merge with donor metadata for condition
    meta = pd.read_csv(DONOR_META)
    result = (
        donor_means[["xist_score", "y_score", "inferred_sex"]]
        .reset_index()
        .merge(meta[["donor_id", "condition"]], on="donor_id", how="left")
    )

    n_f = (result["inferred_sex"] == "Female").sum()
    n_m = (result["inferred_sex"] == "Male").sum()
    log.info("  Inferred sex: %d Female, %d Male", n_f, n_m)
    log.info("  Sex x Condition:")
    ct = pd.crosstab(result["condition"], result["inferred_sex"])
    for cond in ct.index:
        log.info("    %s: %s", cond, dict(ct.loc[cond]))

    out_path = os.path.join(out_dir, "donor_sex_inference.csv")
    result.to_csv(out_path, index=False)
    log.info("  Saved: %s", out_path)

    return result


# =========================================================================
# Module 2: Promoter accessibility enrichment
# =========================================================================
def _load_sex_deg_classification():
    """Load RNA-seq sex DEG classification with gene symbol mapping."""
    log.info("  Loading sex DEG classification...")
    deg = pd.read_csv(SEX_DEG_CSV)

    # Strip Ensembl version to get base ID
    deg["ensembl_id"] = deg["gene"].str.replace(r"\.\d+$", "", regex=True)

    # Load sex causal scores for symbol mapping
    causal = pd.read_csv(SEX_CAUSAL_CSV)
    id2sym = causal.set_index("ensembl_id")["human_symbol"].to_dict()

    deg["symbol"] = deg["ensembl_id"].map(id2sym)
    n_mapped = deg["symbol"].notna().sum()
    log.info("  Mapped %d/%d DEGs to gene symbols (%.1f%%)",
             n_mapped, len(deg), 100 * n_mapped / len(deg))

    return deg


def _intersect_peaks_promoters(peaks_bed, tss_df, halfwidth=PROMOTER_HALFWIDTH):
    """Find genes with at least one peak in their promoter window.

    Returns set of gene_name values with promoter peaks.
    """
    log.info("  Loading hepatocyte peaks from %s...", peaks_bed)
    peaks = pd.read_csv(
        peaks_bed, sep="\t", header=None, names=["chrom", "start", "end"]
    )
    log.info("  Loaded %d peaks", len(peaks))

    # Build promoter windows
    tss_df = tss_df.copy()
    tss_df["prom_start"] = tss_df["tss"] - halfwidth
    tss_df["prom_end"] = tss_df["tss"] + halfwidth

    # Simple interval overlap: sort both by chrom, then check overlaps
    genes_with_peaks = set()
    for chrom in tss_df["chrom"].unique():
        chrom_peaks = peaks[peaks["chrom"] == chrom].sort_values("start")
        chrom_genes = tss_df[tss_df["chrom"] == chrom]
        if len(chrom_peaks) == 0 or len(chrom_genes) == 0:
            continue
        peak_starts = chrom_peaks["start"].values
        peak_ends = chrom_peaks["end"].values
        for _, gene in chrom_genes.iterrows():
            ps, pe = gene["prom_start"], gene["prom_end"]
            # Binary search for candidate peaks
            left = np.searchsorted(peak_ends, ps, side="right")
            right = np.searchsorted(peak_starts, pe, side="left")
            if left < right:
                genes_with_peaks.add(gene["gene_name"])

    log.info("  %d/%d genes have >=1 hepatocyte peak in promoter (TSS +/- %d bp)",
             len(genes_with_peaks), len(tss_df), halfwidth)
    return genes_with_peaks


def _find_hep_label(adata, label="hepatocyte"):
    """Find the hepatocyte cell type label in an AnnData obs."""
    for col in ["cell_type", "cell_type_transferred"]:
        if col not in adata.obs.columns:
            continue
        for val in adata.obs[col].unique():
            if label in str(val).lower():
                return col, val
    return None, None


def _build_pseudobulk(ga, hep_cells_set, donor_sex, disease_only=False):
    """Build CPM-normalized per-donor pseudobulk from gene activity matrix.

    Fixes:
      1. CPM normalization per donor (removes library size confound)
      2. Optional disease-only subset (removes condition composition confound)

    Returns: (donor_means_df, f_donors, m_donors, var_names)
    """
    ga_donor = ga.obs["donor_id"].values
    ga_cell_names = np.array(list(ga.obs_names))

    # Identify hepatocyte cells
    hep_mask = np.array([c in hep_cells_set for c in ga_cell_names])
    log.info("    Hepatocyte cells in gene activity matrix: %d", hep_mask.sum())

    if hep_mask.sum() == 0:
        log.warning("    Cell name matching failed, using all cells")
        hep_mask = np.ones(len(ga_cell_names), dtype=bool)

    ga_hep = ga[hep_mask, :].copy()
    X = ga_hep.X
    if sparse.issparse(X):
        X = X.toarray()

    donors_hep = ga_hep.obs["donor_id"].values
    var_names = list(ga_hep.var_names)

    sex_lookup = donor_sex.set_index("donor_id")["inferred_sex"].to_dict()
    cond_lookup = donor_sex.set_index("donor_id")["condition"].to_dict()

    unique_donors = sorted(set(donors_hep))

    # Optionally restrict to disease donors only
    if disease_only:
        unique_donors = [d for d in unique_donors if cond_lookup.get(d) != "NORMAL"]
        log.info("    Disease-only: %d donors (MASL+MASH)", len(unique_donors))

    # Per-donor mean gene activity with CPM normalization
    donor_means = np.zeros((len(unique_donors), len(var_names)), dtype=np.float64)
    for i, d in enumerate(unique_donors):
        mask = donors_hep == d
        raw_mean = X[mask].mean(axis=0)
        total = raw_mean.sum()
        if total > 0:
            donor_means[i] = raw_mean / total * 1e6  # CPM
        else:
            donor_means[i] = 0.0

    donor_means_df = pd.DataFrame(donor_means, index=unique_donors, columns=var_names)

    f_donors = [d for d in unique_donors if sex_lookup.get(d) == "Female"]
    m_donors = [d for d in unique_donors if sex_lookup.get(d) == "Male"]

    log.info("    Pseudobulk: %d donors x %d genes (CPM-normalized)",
             len(unique_donors), len(var_names))
    log.info("    Female: %d, Male: %d", len(f_donors), len(m_donors))

    # Sanity check: background female-preferential rate
    n_sample = min(5000, len(var_names))
    sample_genes = np.random.choice(var_names, n_sample, replace=False)
    n_fp = sum(
        donor_means_df.loc[f_donors, g].mean() > donor_means_df.loc[m_donors, g].mean()
        for g in sample_genes
    )
    log.info("    Sanity: %.1f%% of random %d genes are female-preferential (expect ~50%%)",
             100 * n_fp / n_sample, n_sample)

    return donor_means_df, f_donors, m_donors, var_names


def _donor_label_permutation_null(
    results_df, donor_means_df, f_donors, m_donors, n_perm=10000, seed=42
):
    """Donor-label permutation null for the sex-class enrichment statistic.

    Fixes review F075 PSEUDOREPLICATION: the gene-set Fisher test (see
    `_run_enrichment_test`) treats ~19k genes as independent observations, but
    every gene's `female_preferential` label is derived from the SAME ~7F/11M
    donor split. Genes are massively non-independent, so the Fisher p is
    anticonservative. The true experimental unit is the DONOR.

    Permutation design
    ------------------
    - Pool = donors with a defined sex present in this subset (f_donors + m_donors).
      Sex labels are shuffled within this pool, holding the observed counts
      (|f_donors| "Female", |m_donors| "Male") fixed. The donor is thus the unit
      of resampling; the gene-level structure is left intact within each draw.
    - For each permutation we recompute every gene's preferential call from the
      per-donor CPM pseudobulk (mean over permuted-Female donors vs mean over
      permuted-Male donors) and, per sex-class gene set, the fraction of genes
      preferential in the tested direction (female_preferential for
      Female_biased/Divergent, male_preferential for Male_biased) -- exactly the
      observed statistic.
    - Empirical p = (1 + #{perm stat >= observed}) / (n_perm + 1), one-sided
      (enrichment), matching the `alternative="greater"` Fisher direction.

    Returns
    -------
    dict[sex_class] -> {"perm_p", "perm_mean_frac", "n_perm"}; plus the special
    key "__neff__" = (n_female_eff, n_male_eff) donor counts used.
    """
    rng = np.random.RandomState(seed)

    pool = list(f_donors) + list(m_donors)
    n_f = len(f_donors)
    n_pool = len(pool)

    # Per-gene per-donor CPM for the pooled (sexed) donors, genes x donors so a
    # boolean donor mask selects columns. Restrict to genes present in the
    # pseudobulk (testable genes are guaranteed to be, but be defensive).
    out = {"__neff__": (n_f, len(m_donors))}
    if n_f == 0 or len(m_donors) == 0 or n_pool < 2:
        log.warning("    Permutation null skipped: need >=1 Female and >=1 Male "
                    "donor (have %dF/%dM)", n_f, len(m_donors))
        for sex_class in ["Female_biased", "Male_biased", "Divergent"]:
            out[sex_class] = {"perm_p": np.nan, "perm_mean_frac": np.nan,
                              "n_perm": 0}
        return out

    # Genes actually tested, grouped by sex_class, in pseudobulk-column order.
    class_genes = {}
    for sex_class in ["Female_biased", "Male_biased", "Divergent"]:
        g = [
            gene for gene in results_df.loc[
                results_df["sex_class"] == sex_class, "gene"
            ]
            if gene in donor_means_df.columns
        ]
        class_genes[sex_class] = g

    all_genes = sorted(
        {g for gs in class_genes.values() for g in gs}
    )
    if len(all_genes) == 0:
        for sex_class in class_genes:
            out[sex_class] = {"perm_p": np.nan, "perm_mean_frac": np.nan,
                              "n_perm": 0}
        return out

    # genes x donors matrix (CPM), restricted to the sexed-donor pool.
    cpm = donor_means_df.loc[pool, all_genes].to_numpy(dtype=np.float64).T  # (G, P)
    gene_pos = {g: i for i, g in enumerate(all_genes)}
    # Precompute per-class gene-row indices once (hoisted out of the perm loop).
    class_idx = {
        sc: np.array([gene_pos[g] for g in genes], dtype=int)
        for sc, genes in class_genes.items() if len(genes) > 0
    }

    # Observed fraction preferential per class (recompute from cpm so the
    # permutation statistic and the observed statistic use identical machinery).
    obs_f_mask = np.zeros(n_pool, dtype=bool)
    obs_f_mask[:n_f] = True  # pool is [f_donors..., m_donors...]
    obs_female_pref = (
        cpm[:, obs_f_mask].mean(axis=1) > cpm[:, ~obs_f_mask].mean(axis=1)
    )  # (G,) per-gene, observed labels

    obs_frac = {}
    perm_ge = {}
    for sex_class, idx in class_idx.items():
        pref = obs_female_pref[idx]
        if sex_class == "Male_biased":
            pref = ~pref  # tested direction is male_preferential
        obs_frac[sex_class] = float(pref.mean())
        perm_ge[sex_class] = 0  # count of perm stat >= observed

    # Permutation loop: shuffle which pool columns are "Female".
    perm_frac_sum = {sc: 0.0 for sc in obs_frac}
    for _ in range(n_perm):
        perm = rng.permutation(n_pool)
        f_cols = perm[:n_f]
        m_cols = perm[n_f:]
        female_pref = (
            cpm[:, f_cols].mean(axis=1) > cpm[:, m_cols].mean(axis=1)
        )  # (G,)
        for sex_class, idx in class_idx.items():
            pref = female_pref[idx]
            if sex_class == "Male_biased":
                pref = ~pref
            frac = float(pref.mean())
            perm_frac_sum[sex_class] += frac
            if frac >= obs_frac[sex_class]:
                perm_ge[sex_class] += 1

    for sex_class in class_genes:
        if sex_class not in obs_frac:
            out[sex_class] = {"perm_p": np.nan, "perm_mean_frac": np.nan,
                              "n_perm": 0}
            continue
        out[sex_class] = {
            "perm_p": (1 + perm_ge[sex_class]) / (n_perm + 1),
            "perm_mean_frac": perm_frac_sum[sex_class] / n_perm,
            "n_perm": n_perm,
        }
    return out


def _run_enrichment_test(results_df, label, donor_means_df=None,
                         f_donors=None, m_donors=None,
                         n_perm=10000, seed=42):
    """Enrichment of sex-class genes for sex-preferential promoter accessibility.

    PRIMARY test: donor-label permutation null (`perm_p`), which makes the
    DONOR the unit of inference and is robust to the gene non-independence
    flagged in review F075.

    The legacy gene-as-unit Fisher's exact test is retained for provenance as
    `fisher_p_gene_level_nominal` (ANTICONSERVATIVE: it treats ~19k genes as
    independent observations though all share the same ~7F/11M donor split --
    do NOT cite as the headline p-value).
    """
    enrichment = []
    concordant_df = results_df[results_df["sex_class"] == "Concordant"]
    concordant_fp = int(concordant_df["female_preferential"].sum())
    concordant_not_fp = len(concordant_df) - concordant_fp

    log.info("  [%s] Concordant background: %d/%d (%.1f%%) female-preferential",
             label, concordant_fp, len(concordant_df),
             100 * concordant_fp / max(len(concordant_df), 1))

    # Donor-label permutation null (primary). Requires the per-donor CPM matrix.
    perm = None
    if donor_means_df is not None and f_donors is not None and m_donors is not None:
        log.info("  [%s] Donor-label permutation null: N=%d perms, seed=%d "
                 "(unit = donor; %dF/%dM)", label, n_perm, seed,
                 len(f_donors), len(m_donors))
        perm = _donor_label_permutation_null(
            results_df, donor_means_df, f_donors, m_donors,
            n_perm=n_perm, seed=seed,
        )
    n_female_eff, n_male_eff = (perm["__neff__"] if perm else (len(f_donors or []),
                                                               len(m_donors or [])))

    for sex_class in ["Female_biased", "Male_biased", "Divergent"]:
        cls_df = results_df[results_df["sex_class"] == sex_class]
        if len(cls_df) == 0:
            continue

        if sex_class == "Male_biased":
            cls_pref = int((~cls_df["female_preferential"]).sum())
            cls_not_pref = len(cls_df) - cls_pref
            ref_pref = concordant_not_fp
            ref_not_pref = concordant_fp
            direction_label = "male_preferential"
        else:
            cls_pref = int(cls_df["female_preferential"].sum())
            cls_not_pref = len(cls_df) - cls_pref
            ref_pref = concordant_fp
            ref_not_pref = concordant_not_fp
            direction_label = "female_preferential"

        table = [[cls_pref, cls_not_pref], [ref_pref, ref_not_pref]]
        # Legacy gene-as-unit Fisher (anticonservative; kept for provenance only).
        odds_ratio, fisher_p_gene = stats.fisher_exact(table, alternative="greater")
        ci = _fisher_ci(table)

        perm_p = perm[sex_class]["perm_p"] if perm else np.nan

        enrichment.append({
            "sex_class": sex_class,
            "subset": label,
            "direction_tested": direction_label,
            "n_genes": len(cls_df),
            "n_preferential": cls_pref,
            "frac_preferential": cls_pref / len(cls_df) if len(cls_df) > 0 else 0,
            "concordant_n_genes": len(concordant_df),
            "concordant_frac_preferential": ref_pref / max(len(concordant_df), 1),
            "odds_ratio": float(odds_ratio),
            # PRIMARY: donor is the unit of inference (review F075 fix).
            "perm_p": float(perm_p) if perm_p is not None and not np.isnan(perm_p) else np.nan,
            "n_female_donors": int(n_female_eff),
            "n_male_donors": int(n_male_eff),
            "n_perm": int(perm[sex_class]["n_perm"]) if perm else 0,
            "perm_null_mean_frac": (
                float(perm[sex_class]["perm_mean_frac"])
                if perm and not np.isnan(perm[sex_class]["perm_mean_frac"]) else np.nan
            ),
            # Legacy gene-as-unit Fisher p -- ANTICONSERVATIVE, provenance only.
            "fisher_p_gene_level_nominal": float(fisher_p_gene),
            "ci_lower": ci[0],
            "ci_upper": ci[1],
        })

        log.info("  [%s] %s: %d/%d %s (%.1f%%), OR=%.2f, perm_p=%s "
                 "(fisher_gene_nominal=%.2e)",
                 label, sex_class, cls_pref, len(cls_df), direction_label,
                 100 * cls_pref / max(len(cls_df), 1), odds_ratio,
                 ("%.4f" % perm_p) if perm_p is not None and not np.isnan(perm_p) else "NA",
                 fisher_p_gene)

    return pd.DataFrame(enrichment)


def module2_promoter_enrichment(ga, donor_sex, tss_df, out_dir):
    """Test enrichment of sex-class genes for sex-preferential promoter accessibility.

    Fixes applied:
      1. CPM normalization per donor (removes library size confound)
      2. Disease-only subset (removes condition composition confound)
      3. Reports both all-donor and disease-only results
    """
    log.info("=" * 60)
    log.info("MODULE 2: Promoter accessibility enrichment (normalized)")
    log.info("=" * 60)

    # Load sex DEG classification
    deg = _load_sex_deg_classification()

    # Identify genes with hepatocyte promoter peaks
    genes_with_peaks = _intersect_peaks_promoters(HEP_PEAKS_BED, tss_df)

    # Load cell type labels to identify hepatocytes
    log.info("  Loading cell types from label-transferred h5ad...")
    lt = ad.read_h5ad(LABEL_TRANSFERRED_H5AD, backed="r")
    ct_col, hep_label = _find_hep_label(lt)
    if hep_label is None:
        log.error("Cannot find hepatocyte label")
        return None, None
    hep_cells = lt.obs[lt.obs[ct_col] == hep_label].index.tolist()
    log.info("  Found %d hepatocyte cells (label='%s')", len(hep_cells), hep_label)

    sex_lookup = donor_sex.set_index("donor_id")["inferred_sex"].to_dict()
    lt_obs = lt.obs.loc[hep_cells, ["donor_id"]].copy()
    lt_obs["inferred_sex"] = lt_obs["donor_id"].map(sex_lookup)
    log.info("  Hepatocyte cells: %d Female, %d Male",
             (lt_obs["inferred_sex"] == "Female").sum(),
             (lt_obs["inferred_sex"] == "Male").sum())
    lt.file.close()

    hep_cells_set = set(hep_cells)

    # Filter to testable genes
    deg_with_symbol = deg.dropna(subset=["symbol"])
    sym2class = {}
    for _, row in deg_with_symbol.iterrows():
        sym2class.setdefault(row["symbol"], row["sex_class"])

    all_results = []
    all_enrichment = []

    # Run for both subsets: all donors and disease-only
    np.random.seed(42)
    for disease_only, label in [(False, "all_donors"), (True, "disease_only")]:
        log.info("  --- Subset: %s ---", label)
        donor_means_df, f_donors, m_donors, var_names = _build_pseudobulk(
            ga, hep_cells_set, donor_sex, disease_only=disease_only
        )

        testable_genes = [
            g for g in sym2class if g in var_names and g in genes_with_peaks
        ]
        log.info("    Testable genes: %d", len(testable_genes))

        results = []
        for gene in testable_genes:
            f_vals = donor_means_df.loc[f_donors, gene].values.astype(float)
            m_vals = donor_means_df.loc[m_donors, gene].values.astype(float)

            if np.std(np.concatenate([f_vals, m_vals])) < 1e-10:
                continue

            try:
                res = stats.mannwhitneyu(f_vals, m_vals, alternative="two-sided")
                pval = float(res.pvalue)
            except ValueError:
                continue

            mean_f = float(f_vals.mean())
            mean_m = float(m_vals.mean())

            results.append({
                "gene": gene,
                "sex_class": sym2class[gene],
                "subset": label,
                "mean_F_cpm": mean_f,
                "mean_M_cpm": mean_m,
                "fold_change": float(np.log2((mean_f + 1e-3) / (mean_m + 1e-3))),
                "wilcoxon_p": pval,
                "female_preferential": mean_f > mean_m,
            })

        results_df = pd.DataFrame(results)
        log.info("    Tested %d genes", len(results_df))

        # Run enrichment. Pass the per-donor CPM pseudobulk + sex split so the
        # primary test is a DONOR-LABEL permutation null (review F075 fix);
        # the gene-as-unit Fisher is kept only as fisher_p_gene_level_nominal.
        enrichment_df = _run_enrichment_test(
            results_df, label,
            donor_means_df=donor_means_df, f_donors=f_donors, m_donors=m_donors,
            n_perm=10000, seed=42,
        )

        all_results.append(results_df)
        all_enrichment.append(enrichment_df)

    # Save combined results
    combined_results = pd.concat(all_results, ignore_index=True)
    out_path = os.path.join(out_dir, "promoter_accessibility_by_sex_class.csv")
    combined_results.to_csv(out_path, index=False)
    log.info("  Saved: %s", out_path)

    combined_enrichment = pd.concat(all_enrichment, ignore_index=True)
    out_path = os.path.join(out_dir, "promoter_enrichment_summary.csv")
    combined_enrichment.to_csv(out_path, index=False)
    log.info("  Saved: %s", out_path)

    # Return disease-only results for downstream (Module 4)
    disease_results = all_results[1] if len(all_results) > 1 else all_results[0]
    return disease_results, combined_enrichment


def _fisher_ci(table, alpha=0.05):
    """Approximate 95% CI for odds ratio from 2x2 table using log-OR SE."""
    a, b = table[0]
    c, d = table[1]
    # Add 0.5 continuity correction if any cell is zero
    if a == 0 or b == 0 or c == 0 or d == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    log_or = np.log((a * d) / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    z = stats.norm.ppf(1 - alpha / 2)
    return (np.exp(log_or - z * se), np.exp(log_or + z * se))


# =========================================================================
# Module 3: chromVAR sex-differential TF activity
# =========================================================================
def module3_chromvar_sex(donor_sex, deg, out_dir):
    """Compare chromVAR TF motif deviations between sexes in hepatocytes."""
    log.info("=" * 60)
    log.info("MODULE 3: chromVAR sex-differential TF activity")
    log.info("=" * 60)

    log.info("  Loading chromVAR deviations...")
    cv = ad.read_h5ad(CHROMVAR_H5AD)

    # Determine hepatocyte label
    ct_col = "cell_type"
    ct_vals = cv.obs[ct_col].unique().tolist()
    hep_label = None
    for candidate in ["Hepatocytes", "Hepatocyte"]:
        if candidate in ct_vals:
            hep_label = candidate
            break
    if "cell_type_transferred" in cv.obs.columns and hep_label is None:
        ct_col = "cell_type_transferred"
        ct_vals = cv.obs[ct_col].unique().tolist()
        for candidate in ["Hepatocytes", "Hepatocyte"]:
            if candidate in ct_vals:
                hep_label = candidate
                break
    if hep_label is None:
        log.error("Cannot find hepatocyte label in chromVAR h5ad")
        return None

    hep_mask = cv.obs[ct_col] == hep_label
    cv_hep = cv[hep_mask, :].copy()
    log.info("  Hepatocyte cells in chromVAR: %d (label='%s' in '%s')",
             cv_hep.n_obs, hep_label, ct_col)

    sex_lookup = donor_sex.set_index("donor_id")["inferred_sex"].to_dict()

    # Per-donor pseudobulk TF deviations
    X = cv_hep.X
    if sparse.issparse(X):
        X = X.toarray()

    donors = cv_hep.obs["donor_id"].values
    tf_names = list(cv_hep.var_names)
    unique_donors = sorted(set(donors))

    donor_tf = np.zeros((len(unique_donors), len(tf_names)), dtype=np.float32)
    for i, d in enumerate(unique_donors):
        mask = donors == d
        if mask.sum() > 0:
            donor_tf[i] = X[mask].mean(axis=0)

    donor_tf_df = pd.DataFrame(donor_tf, index=unique_donors, columns=tf_names)

    # Deduplicate TF columns (chromVAR var_names may not be unique)
    if donor_tf_df.columns.duplicated().any():
        n_dup = donor_tf_df.columns.duplicated().sum()
        log.warning("  %d duplicate TF names in chromVAR; keeping first occurrence", n_dup)
        donor_tf_df = donor_tf_df.loc[:, ~donor_tf_df.columns.duplicated()]
        tf_names = list(donor_tf_df.columns)

    f_donors = [d for d in unique_donors if sex_lookup.get(d) == "Female"]
    m_donors = [d for d in unique_donors if sex_lookup.get(d) == "Male"]

    # Load disease regulon TFs
    disease_tfs = set()
    if os.path.exists(DISEASE_REGULONS):
        dr = pd.read_csv(DISEASE_REGULONS)
        disease_tfs = set(dr["tf_name"].unique())
        log.info("  Disease regulon TFs: %d", len(disease_tfs))

    hep_tfs = set()
    if os.path.exists(HEP_REGULONS):
        hr = pd.read_csv(HEP_REGULONS)
        hep_tfs = set(hr["tf_name"].unique())
        log.info("  Hepatocyte regulon TFs: %d", len(hep_tfs))

    # Build symbol->sex_class lookup from RNA DEG classification
    deg_with_sym = deg.dropna(subset=["symbol"])
    sym2class = deg_with_sym.set_index("symbol")["sex_class"].to_dict()

    # Per-TF Wilcoxon
    results = []
    for tf in tf_names:
        f_vals = donor_tf_df.loc[f_donors, tf].values.astype(float).ravel()
        m_vals = donor_tf_df.loc[m_donors, tf].values.astype(float).ravel()

        if np.std(np.concatenate([f_vals, m_vals])) < 1e-10:
            continue

        try:
            res = stats.mannwhitneyu(f_vals, m_vals, alternative="two-sided")
            pval = float(res.pvalue)
        except ValueError:
            continue

        results.append({
            "tf_name": tf,
            "mean_dev_F": float(f_vals.mean()),
            "mean_dev_M": float(m_vals.mean()),
            "delta_dev": float(f_vals.mean() - m_vals.mean()),
            "wilcoxon_p": pval,
            "in_disease_regulon": tf in disease_tfs,
            "in_hep_regulon": tf in hep_tfs,
            "rna_sex_class": sym2class.get(tf, "NA"),
        })

    results_df = pd.DataFrame(results)

    # BH FDR correction
    if len(results_df) > 0:
        pvals = results_df["wilcoxon_p"].astype(float).values
        _, padj, _, _ = multipletests(pvals, method="fdr_bh")
        results_df["padj"] = padj
    else:
        results_df["padj"] = []

    n_sig = (results_df["padj"] < 0.05).sum() if len(results_df) > 0 else 0
    log.info("  Tested %d TF motifs, %d significant (padj < 0.05)", len(results_df), n_sig)

    # Report disease regulon TFs
    dr_results = results_df[results_df["in_disease_regulon"]]
    if len(dr_results) > 0:
        dr_sig = dr_results[dr_results["padj"] < 0.05]
        log.info("  Disease regulon TFs tested: %d, significant: %d", len(dr_results), len(dr_sig))
        for _, row in dr_sig.iterrows():
            log.info("    %s: delta=%.4f, padj=%.4f", row["tf_name"], row["delta_dev"], row["padj"])

    # Sort by p-value
    results_df = results_df.sort_values("wilcoxon_p")

    out_path = os.path.join(out_dir, "chromvar_sex_differential.csv")
    results_df.to_csv(out_path, index=False)
    log.info("  Saved: %s", out_path)

    return results_df


# =========================================================================
# Module 4: Gene activity concordance with RNA-seq
# =========================================================================
def module4_concordance(promo_results, deg, out_dir):
    """Correlate RNA sex LFC with ATAC sex LFC (CPM-normalized, disease-only)."""
    log.info("=" * 60)
    log.info("MODULE 4: Gene activity concordance (RNA vs ATAC, normalized)")
    log.info("=" * 60)

    if promo_results is None or len(promo_results) == 0:
        log.warning("  No promoter results available, skipping concordance")
        return None

    # RNA sex LFC = logFC_F - logFC_M
    deg_with_sym = deg.dropna(subset=["symbol"])
    rna_lfc = deg_with_sym.set_index("symbol")[["logFC_F", "logFC_M", "sex_class"]].copy()
    rna_lfc["rna_sex_lfc"] = rna_lfc["logFC_F"] - rna_lfc["logFC_M"]

    # ATAC sex LFC from promoter results (CPM-normalized)
    atac_lfc = promo_results.set_index("gene")[["fold_change", "sex_class"]].copy()
    atac_lfc.rename(columns={"fold_change": "atac_sex_lfc"}, inplace=True)

    # Join
    joined = rna_lfc.join(atac_lfc[["atac_sex_lfc"]], how="inner")
    joined = joined.dropna(subset=["rna_sex_lfc", "atac_sex_lfc"])
    joined = joined[np.isfinite(joined["rna_sex_lfc"]) & np.isfinite(joined["atac_sex_lfc"])]
    log.info("  Genes with both RNA and ATAC sex LFC: %d", len(joined))
    log.info("  NOTE: ATAC LFC is CPM-normalized, disease-only (6F vs 7M)")

    # Global Spearman correlation
    if len(joined) > 10:
        rho, pval = stats.spearmanr(joined["rna_sex_lfc"], joined["atac_sex_lfc"])
        log.info("  Global Spearman: rho=%.4f, p=%.2e (n=%d)", rho, pval, len(joined))
    else:
        rho, pval = np.nan, np.nan

    # Stratified correlations
    strat_corrs = []
    for cls in ["Female_biased", "Male_biased", "Divergent", "Concordant"]:
        sub = joined[joined["sex_class"] == cls]
        if len(sub) > 10:
            r, p = stats.spearmanr(sub["rna_sex_lfc"], sub["atac_sex_lfc"])
        else:
            r, p = np.nan, np.nan
        strat_corrs.append({"sex_class": cls, "n": len(sub), "spearman_rho": r, "spearman_p": p})
        if not np.isnan(r):
            log.info("  %s: rho=%.4f, p=%.2e (n=%d)", cls, r, p, len(sub))

    # Save concordance data
    joined_out = joined.reset_index().rename(columns={"index": "gene"})
    out_path = os.path.join(out_dir, "gene_activity_sex_concordance.csv")
    joined_out.to_csv(out_path, index=False)
    log.info("  Saved: %s", out_path)

    return {
        "global_rho": float(rho) if not np.isnan(rho) else None,
        "global_p": float(pval) if not np.isnan(pval) else None,
        "n_genes": len(joined),
        "stratified": strat_corrs,
    }


# =========================================================================
# Main
# =========================================================================
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default=os.path.join(ATAC_DIR, "results/sex_validation"),
        help="Output directory for results",
    )
    args = parser.parse_args()

    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)
    log.info("Output directory: %s", out_dir)

    # Parse TSS
    tss_df = parse_tss_from_gtf(GENCODE_GTF)

    # Load gene activity matrix
    log.info("Loading gene activity matrix (this may take a few minutes)...")
    ga = ad.read_h5ad(GENE_ACTIVITY_H5AD)
    log.info("  Shape: %s", ga.shape)

    # Module 1: Sex inference
    donor_sex = module1_sex_inference(ga, out_dir)
    if donor_sex is None:
        log.error("Sex inference failed, aborting")
        return 1

    # Load DEG classification (needed by modules 2, 3, 4)
    deg = _load_sex_deg_classification()

    # Module 2: Promoter accessibility enrichment
    promo_results, enrichment = module2_promoter_enrichment(ga, donor_sex, tss_df, out_dir)

    # Module 3: chromVAR sex-differential TF activity
    chromvar_results = module3_chromvar_sex(donor_sex, deg, out_dir)

    # Module 4: Gene activity concordance
    concordance = module4_concordance(promo_results, deg, out_dir)

    # Headline enrichment result: the DONOR-LABEL permutation p (review F075).
    # The legacy gene-as-unit Fisher p (fisher_p_gene_level_nominal) is
    # anticonservative and is NOT the headline.
    headline = None
    if enrichment is not None and len(enrichment) > 0:
        # Prefer the female-biased / disease-only cell as the primary statistic.
        prio = enrichment.copy()
        prio["_rank"] = (
            (prio["sex_class"] == "Female_biased").astype(int) * 2
            + (prio["subset"] == "disease_only").astype(int)
        )
        top = prio.sort_values("_rank", ascending=False).iloc[0]
        headline = {
            "test": "donor_label_permutation",
            "unit": "donor",
            "note": ("Primary enrichment p is the donor-label permutation null; "
                     "the gene-level Fisher p is anticonservative "
                     "(pseudoreplication, review F075) and reported only as "
                     "fisher_p_gene_level_nominal."),
            "sex_class": str(top["sex_class"]),
            "subset": str(top["subset"]),
            "direction_tested": str(top["direction_tested"]),
            "frac_preferential": float(top["frac_preferential"]),
            "perm_p": (float(top["perm_p"]) if pd.notna(top["perm_p"]) else None),
            "n_perm": int(top["n_perm"]),
            "n_female_donors": int(top["n_female_donors"]),
            "n_male_donors": int(top["n_male_donors"]),
            "fisher_p_gene_level_nominal": float(top["fisher_p_gene_level_nominal"]),
        }

    # Write summary JSON
    summary = {
        "n_donors": len(donor_sex),
        "n_female": int((donor_sex["inferred_sex"] == "Female").sum()),
        "n_male": int((donor_sex["inferred_sex"] == "Male").sum()),
        "n_testable_genes": int(len(promo_results)) if promo_results is not None else 0,
        "enrichment_headline": headline,
        "enrichment": enrichment.to_dict("records") if enrichment is not None else [],
        "n_tfs_tested": int(len(chromvar_results)) if chromvar_results is not None else 0,
        "n_tfs_significant": int((chromvar_results["padj"] < 0.05).sum()) if chromvar_results is not None else 0,
        "concordance": concordance if concordance else {},
    }

    summary_path = os.path.join(out_dir, "validation_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log.info("Saved summary: %s", summary_path)

    log.info("=" * 60)
    log.info("All modules complete")
    log.info("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
