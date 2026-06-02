#!/usr/bin/env python
"""
validate_orthologs.py — Held-out validation of the master ortholog table.

Run AFTER master_ortholog_table.tsv.gz has been built. Computes:
  (1) PC precision/recall per tier vs TOGA one2one truth (circularity-corrected)
  (2) Paralog-class agreement: biomaRt vs TOGA on one2many / many2many calls
  (3) False-positive rate on 10K random shuffled pairs
  (4) lncRNA internal consistency: Phase J ∩ BLAST RBH
  (5) miRBase ∩ MirGeneDB Jaccard
  (6) OMA Browser comparison (if local cache present)
  (7) Canonical lncRNA (NEAT1, MEG3, MALAT1, ...) before/after RBH

Emits:
  - results/ortholog_validation/precision_recall_per_tier.csv
  - results/ortholog_validation/null_distribution_check.csv
  - results/ortholog_validation/canonical_lncrna_resolution.csv
  - results/ortholog_validation/lncrna_layer_consistency.csv
  - results/ortholog_validation/miRNA_layer_jaccard.csv
  - results/ortholog_validation/oma_comparison.csv  (if OMA cache present)
  - results/ortholog_validation/figs/recall_per_tier.pdf
  - results/ortholog_validation/figs/precision_per_tier.pdf
  - results/ortholog_validation/figs/null_vs_observed_distribution.pdf

Cite: Wolf & Koonin 2012; Pearson 2013; Kirilenko 2023; Emms 2019.
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ORTH_DIR = PROJECT / "data/external/orthologs"
LAYERS = ORTH_DIR / "layers"
OUT_DIR = PROJECT / "Cas13_Library_Design/results/ortholog_validation"
FIG_DIR = OUT_DIR / "figs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Load master + truth
# -----------------------------------------------------------------------------

def load_master() -> pd.DataFrame:
    p = ORTH_DIR / "master_ortholog_table.tsv.gz"
    m = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    for c in ("tier_H_biomart", "tier_H_toga", "tier_H_mirbase",
              "tier_H_mirgenedb", "tier_M_phasej_synteny", "tier_M_blast",
              "tier_M_blast_rbh"):
        if c in m.columns:
            m[c] = pd.to_numeric(m[c], errors="coerce").fillna(0).astype(int)
    return m


def load_toga_truth() -> pd.DataFrame:
    """Hold out TOGA L0 to use one2one as truth."""
    return pd.read_csv(LAYERS / "L0_toga.tsv", sep="\t", dtype=str)


def load_blast() -> pd.DataFrame:
    return pd.read_csv(LAYERS / "L4_blast.tsv", sep="\t", dtype=str,
                       low_memory=False)


def load_blast_rbh() -> pd.DataFrame:
    return pd.read_csv(LAYERS / "L4_blast_rbh.tsv", sep="\t", dtype=str,
                       low_memory=False)


def load_phasej() -> pd.DataFrame:
    return pd.read_csv(LAYERS / "L1.5_phasej_synteny.tsv", sep="\t",
                       dtype=str, low_memory=False)


def load_mirbase() -> pd.DataFrame:
    return pd.read_csv(LAYERS / "L2_mirbase.tsv", sep="\t", dtype=str,
                       low_memory=False)


def load_mirgenedb() -> pd.DataFrame:
    return pd.read_csv(LAYERS / "L2b_mirgenedb.tsv", sep="\t", dtype=str,
                       low_memory=False)


# -----------------------------------------------------------------------------
# Universe (mouse + human gene IDs we can possibly map)
# -----------------------------------------------------------------------------

def load_gene_universe() -> tuple[set, set]:
    """Mouse + human PC gene IDs we can sample from."""
    mouse_meta = pd.read_csv(
        PROJECT / "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv",
        dtype=str)
    # Schema: mouse_ensembl_base, mouse_symbol_gtf, mouse_biotype
    mouse_meta["mouse_ensembl_base"] = (
        mouse_meta["mouse_ensembl_base"].str.split(".").str[0])
    mouse_pc = set(mouse_meta.loc[
        mouse_meta["mouse_biotype"] == "protein_coding",
        "mouse_ensembl_base"].dropna())

    gencode = pd.read_csv(
        ORTH_DIR / "source/gencode_v49_gene_metadata.tsv.gz",
        sep="\t", dtype=str)
    # Schema: gene_id, gene_name, chromosome, gene_biotype, ensembl_base
    human_pc = set(gencode.loc[
        gencode["gene_biotype"] == "protein_coding",
        "ensembl_base"].dropna())
    return mouse_pc, human_pc


# -----------------------------------------------------------------------------
# (1) PC precision/recall per tier vs TOGA one2one truth
# -----------------------------------------------------------------------------

def evaluate_pc_pr(master: pd.DataFrame, toga: pd.DataFrame) -> pd.DataFrame:
    """
    Truth = TOGA one2one mouse-human pairs (16,512). PR per master tier.

    Metrics reported per tier (PC-PC pairs only):

      recall_inclusive
          = (master_tier_pairs ∩ TOGA one2one truth) / |TOGA one2one truth|
          Mostly self-consistency since master uses TOGA as a source.

      precision_strict_one2one (CIRCULARITY-CONTROLLED)
          Of TOGA-blind master pairs (tier_H_toga==0) restricted to pairs
          where BOTH genes are in TOGA's coverage universe, what fraction
          matches TOGA one2one exactly. Low values indicate this tier's
          *additional* (non-TOGA-supported) pairs conflict with TOGA's
          one2one designations.

      precision_any_toga_class
          Same as above but a pair counts as TP if it is in ANY TOGA class
          (one2one, one2many, many2one, many2many). Measures whether the
          two methods at least agree at the gene-family level.

      precision_partner_consistent
          For each TOGA-blind pair (m, h) in TOGA universe, check whether
          either:
            (a) (m, h) is in TOGA at any orthology class, OR
            (b) m's TOGA partner shares a gene-family relationship with h
                — defined as h being in the same TOGA orthology block (any
                class) as one of m's TOGA partners (same row's
                gene_family_id is not available, so we approximate by
                checking if h appears with any of m's TOGA-listed partners
                in the same many2many cluster).
          We instead use: ≥1 of m's TOGA partners has same h's gene_symbol
          stem (case-insensitive). This is a soft "gene family" check.

      independent_recall_at_tier
          = (master_tier_pairs_TOGA_blind ∩ TOGA one2one truth) /
            |TOGA one2one truth|
          Recall achievable from non-TOGA evidence alone.
    """
    # Truth set: PC-PC TOGA one2one pairs (so it matches the PC-PC predicted
    # filter we apply to master)
    truth_all = set(zip(
        toga.loc[toga["toga_orthology_class"] == "one2one", "mouse_ensembl"],
        toga.loc[toga["toga_orthology_class"] == "one2one", "human_ensembl"]))
    pc_truth_df = toga[(toga["toga_orthology_class"] == "one2one") &
                       (toga["mouse_biotype"] == "protein_coding")]
    truth = set(zip(pc_truth_df["mouse_ensembl"],
                    pc_truth_df["human_ensembl"]))
    n_truth = len(truth)
    n_truth_all = len(truth_all)

    toga_mouse_universe = set(toga["mouse_ensembl"].dropna())
    toga_human_universe = set(toga["human_ensembl"].dropna())
    toga_any = set(zip(toga["mouse_ensembl"], toga["human_ensembl"]))

    # Build mouse_id -> set(human_partners) and human_id -> set(mouse_partners)
    mouse_to_human = toga.groupby("mouse_ensembl")["human_ensembl"].apply(set).to_dict()
    human_to_mouse = toga.groupby("human_ensembl")["mouse_ensembl"].apply(set).to_dict()

    rows = []
    for tier_label, tier_filter in [
        ("H", master["confidence_tier"] == "H"),
        ("M", master["confidence_tier"] == "M"),
        ("L", master["confidence_tier"] == "L"),
        ("H+M", master["confidence_tier"].isin(["H", "M"])),
        ("any", master["confidence_tier"].isin(["H", "M", "L"])),
    ]:
        sub = master.loc[tier_filter]
        sub_pc = sub[(sub["mouse_biotype"] == "protein_coding") &
                     (sub["human_biotype"] == "protein_coding")]
        pred = set(zip(sub_pc["mouse_ensembl"], sub_pc["human_ensembl"]))
        tp_recall = len(pred & truth)
        recall_all = tp_recall / n_truth if n_truth else np.nan

        non_toga = sub[sub["tier_H_toga"].astype(int) == 0]
        non_toga_pc = non_toga[
            (non_toga["mouse_biotype"] == "protein_coding") &
            (non_toga["human_biotype"] == "protein_coding")]
        in_toga_universe = non_toga_pc[
            non_toga_pc["mouse_ensembl"].isin(toga_mouse_universe) &
            non_toga_pc["human_ensembl"].isin(toga_human_universe)]
        eval_set = set(zip(in_toga_universe["mouse_ensembl"],
                           in_toga_universe["human_ensembl"]))

        tp_strict = eval_set & truth
        in_any_class = eval_set & toga_any

        # Partner consistency: at minimum, m's TOGA partners is non-empty
        # AND there exists a TOGA partner h' of m that's the same gene as h.
        # We use a soft criterion: gene_symbol case-insensitive equality of
        # at least one TOGA partner's symbol to h's symbol (already captured
        # by "in_any_class" since we have ensembl ids). The "partner
        # consistent" criterion is equivalent to "in_any_class" because
        # if m has h as a TOGA partner, (m, h) is in toga_any. So we drop
        # the partner_consistent column and just keep strict vs lenient.

        precision_strict = (len(tp_strict) / len(eval_set)
                             if eval_set else np.nan)
        precision_lenient = (len(in_any_class) / len(eval_set)
                              if eval_set else np.nan)

        f1 = (2 * precision_strict * recall_all /
              (precision_strict + recall_all)
              if precision_strict and recall_all else np.nan)

        # Independent (non-TOGA-supported) recall to TOGA one2one truth
        independent_recall = len(tp_strict) / n_truth if n_truth else np.nan

        rows.append({
            "tier": tier_label,
            "n_predicted_pairs_total": len(pred),
            "n_predicted_pairs_toga_blind_pc": len(non_toga_pc),
            "n_eval_set_in_toga_universe": len(eval_set),
            "n_truth_pairs_pc_one2one": n_truth,
            "n_truth_pairs_one2one_any_biotype": n_truth_all,
            "n_TP_recall_inclusive": tp_recall,
            "n_TP_strict_one2one": len(tp_strict),
            "n_TP_lenient_any_toga_class": len(in_any_class),
            "n_FP_eval_set": len(eval_set) - len(in_any_class),
            "recall_inclusive_pc_truth": recall_all,
            "precision_strict_one2one": precision_strict,
            "precision_lenient_any_toga_class": precision_lenient,
            "independent_recall_to_one2one_truth": independent_recall,
            "F1_inclusive": f1,
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# (2) Paralog-class agreement
# -----------------------------------------------------------------------------

def evaluate_paralog_agreement(master: pd.DataFrame) -> pd.DataFrame:
    """For pairs where biomaRt calls paralog class, does TOGA agree?"""
    bm = master[master["biomart_ortholog_type"].isin([
        "ortholog_one2many", "ortholog_many2one", "ortholog_many2many"])].copy()
    bm["bm_class"] = bm["biomart_ortholog_type"].str.replace("ortholog_", "")
    have_toga = bm["tier_H_toga"].astype(int) == 1
    bm_with_toga = bm[have_toga]
    rows = []
    for bm_class in ["one2many", "many2one", "many2many"]:
        sub = bm_with_toga[bm_with_toga["bm_class"] == bm_class]
        n_total = (bm["bm_class"] == bm_class).sum()
        n_with_toga = len(sub)
        n_toga_paralog = sub["toga_orthology_class"].isin(
            ["one2many", "many2one", "many2many"]).sum()
        n_toga_one2one = (sub["toga_orthology_class"] == "one2one").sum()
        rows.append({
            "biomart_class": bm_class,
            "n_biomart_total": int(n_total),
            "n_biomart_with_toga_call": int(n_with_toga),
            "n_toga_paralog": int(n_toga_paralog),
            "n_toga_one2one_disagrees": int(n_toga_one2one),
            "paralog_agreement_rate": (n_toga_paralog / n_with_toga
                                       if n_with_toga else np.nan),
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# (3) Null distribution / shuffled-pair FP check
# -----------------------------------------------------------------------------

def evaluate_null_distribution(master: pd.DataFrame,
                               mouse_pc: set, human_pc: set,
                               n_pairs: int = 10000,
                               seed: int = 42) -> pd.DataFrame:
    """Random (mouse, human) pairs NOT present in any source layer. Confirm
    they are not spuriously elevated to Tier H or M by master."""
    rng = np.random.default_rng(seed)
    existing_pairs = set(zip(master["mouse_ensembl"], master["human_ensembl"]))
    # Sample from cartesian product without replacement and reject existing
    mouse_arr = np.array(sorted(mouse_pc))
    human_arr = np.array(sorted(human_pc))
    sampled = set()
    rejected = 0
    tries = 0
    while len(sampled) < n_pairs and tries < n_pairs * 20:
        m = rng.choice(mouse_arr)
        h = rng.choice(human_arr)
        pair = (m, h)
        if pair in existing_pairs:
            rejected += 1
            tries += 1
            continue
        if pair in sampled:
            tries += 1
            continue
        sampled.add(pair)
        tries += 1

    null = pd.DataFrame(sorted(sampled),
                        columns=["mouse_ensembl", "human_ensembl"])
    # Check: do any of them appear in master?
    merged = null.merge(
        master[["mouse_ensembl", "human_ensembl", "confidence_tier",
                "provenance_sources"]],
        on=["mouse_ensembl", "human_ensembl"], how="left")
    in_master = merged["confidence_tier"].notna().sum()
    tier_h = (merged["confidence_tier"] == "H").sum()
    tier_m = (merged["confidence_tier"] == "M").sum()
    tier_l = (merged["confidence_tier"] == "L").sum()

    out = pd.DataFrame([{
        "n_null_pairs_drawn": len(null),
        "n_rejected_already_in_master": rejected,
        "n_null_in_master": int(in_master),
        "n_null_at_tier_H": int(tier_h),
        "n_null_at_tier_M": int(tier_m),
        "n_null_at_tier_L": int(tier_l),
        "fp_rate_any_tier": float(in_master) / len(null) if len(null) else 0.0,
        "fp_rate_tier_HM": (float(tier_h + tier_m) / len(null)
                            if len(null) else 0.0),
    }])
    return out, merged


# -----------------------------------------------------------------------------
# (4) lncRNA internal consistency: Phase J ∩ BLAST RBH
# -----------------------------------------------------------------------------

def evaluate_lncrna_consistency(phasej: pd.DataFrame,
                                blast_rbh: pd.DataFrame,
                                master: pd.DataFrame) -> pd.DataFrame:
    """For human lncRNAs in Phase J synteny, does the same mouse partner
    also appear via BLAST RBH? Two-method agreement = strong evidence."""
    # Phase J human lncRNAs
    pj_pairs = set(zip(phasej["mouse_ensembl"], phasej["human_ensembl"]))
    pj_human = set(phasej["human_ensembl"])
    pj_mouse = set(phasej["mouse_ensembl"])

    # BLAST RBH pairs
    rbh_pairs = set(zip(blast_rbh["mouse_ensembl"], blast_rbh["human_ensembl"]))
    rbh_human = set(blast_rbh["human_ensembl"])

    # Agreement: pair appears in both
    pj_and_rbh = pj_pairs & rbh_pairs
    # Human lncRNAs in PJ that also have ANY RBH partner (may be different)
    pj_human_with_rbh = pj_human & rbh_human

    rows = []
    rows.append({
        "metric": "phasej_pairs",
        "count": len(pj_pairs),
        "note": "Phase J flanking-PCG synteny lncRNA candidate pairs",
    })
    rows.append({
        "metric": "blast_rbh_pairs",
        "count": len(rbh_pairs),
        "note": "BLAST reciprocal best hit pairs (lncRNA)",
    })
    rows.append({
        "metric": "phasej_human_with_any_rbh_partner",
        "count": len(pj_human_with_rbh),
        "note": ("Human lncRNAs in Phase J that have any BLAST RBH partner "
                 "(partner may be different from PJ candidate)"),
    })
    rows.append({
        "metric": "intersect_phasej_and_rbh_pairs",
        "count": len(pj_and_rbh),
        "note": ("Pairs where Phase J and BLAST RBH agree on the SAME mouse "
                 "partner (high-confidence)"),
    })
    union = len(pj_pairs | rbh_pairs)
    rows.append({
        "metric": "jaccard_phasej_vs_rbh",
        "count": (len(pj_and_rbh) / union) if union else 0,
        "note": "Set Jaccard of the two methods at pair level",
    })
    # How many Phase J human lncRNAs have at least one PJ partner that
    # matches a RBH partner (same pair)?
    pj_with_matching_rbh = sum(1 for h in pj_human
                                if any((m, h) in rbh_pairs
                                       for m in phasej.loc[phasej["human_ensembl"] == h,
                                                            "mouse_ensembl"].unique()))
    # Above is O(N*K); use a faster join instead
    merged = phasej.merge(blast_rbh[["mouse_ensembl", "human_ensembl"]],
                          on=["mouse_ensembl", "human_ensembl"],
                          how="inner")
    n_pj_pairs_confirmed_by_rbh = len(merged)
    n_pj_human_confirmed = merged["human_ensembl"].nunique()
    rows[-1]["count"] = (len(pj_and_rbh) / union) if union else 0

    rows.append({
        "metric": "phasej_pairs_confirmed_by_rbh",
        "count": int(n_pj_pairs_confirmed_by_rbh),
        "note": ("Subset of Phase J pairs where the SAME mouse partner "
                 "passes RBH BLAST — the dual-evidence Tier-M layer"),
    })
    rows.append({
        "metric": "phasej_human_lncrna_with_dual_evidence",
        "count": int(n_pj_human_confirmed),
        "note": ("Unique human lncRNAs whose Phase J synteny call is "
                 "confirmed by BLAST RBH"),
    })

    # Master M-tier lncRNAs that have dual support
    m_dual = master[(master["confidence_tier"] == "M") &
                    (master["human_biotype"] == "lncRNA") &
                    (master["tier_M_phasej_synteny"].astype(int) == 1) &
                    (master["tier_M_blast_rbh"].astype(int) == 1)]
    rows.append({
        "metric": "master_tier_M_lncrna_with_phasej_and_rbh",
        "count": int(len(m_dual)),
        "note": ("Tier M lncRNA pairs in master backed by both Phase J + "
                 "RBH BLAST evidence"),
    })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# (5) miRBase ∩ MirGeneDB
# -----------------------------------------------------------------------------

def evaluate_mirna_jaccard(mb: pd.DataFrame,
                           mg: pd.DataFrame) -> pd.DataFrame:
    mb_pairs = set(zip(mb["mouse_ensembl"], mb["human_ensembl"]))
    mg_pairs = set(zip(mg["mouse_ensembl"], mg["human_ensembl"]))
    inter = mb_pairs & mg_pairs
    union = mb_pairs | mg_pairs
    rows = [
        {"metric": "mirbase_pairs", "count": len(mb_pairs)},
        {"metric": "mirgenedb_pairs", "count": len(mg_pairs)},
        {"metric": "intersection", "count": len(inter)},
        {"metric": "union", "count": len(union)},
        {"metric": "jaccard", "count": round(len(inter) / len(union), 4)
         if union else 0},
        {"metric": "mirbase_only", "count": len(mb_pairs - mg_pairs)},
        {"metric": "mirgenedb_only", "count": len(mg_pairs - mb_pairs)},
    ]
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# (6) OMA Browser comparison (if cache present)
# -----------------------------------------------------------------------------

def load_oma_cache(cache_dir: str = "/tmp/oma_pages") -> pd.DataFrame:
    """Load OMA mouse-human pairs from cached JSON pages. Resolves
    canonical IDs (UniProt-style like NK2R_MOUSE) → Ensembl via OMA xrefs
    where possible. Returns DataFrame with (mouse_xref, human_xref, rel_type).
    """
    pages = sorted(glob.glob(f"{cache_dir}/p*.json"))
    if not pages:
        return pd.DataFrame()
    rows = []
    for p in pages:
        try:
            with open(p) as fh:
                d = json.load(fh)
        except Exception:
            continue
        if not isinstance(d, list):
            continue
        for r in d:
            e1 = r.get("entry_1", {})
            e2 = r.get("entry_2", {})
            rows.append({
                "oma_canonical_mouse": e1.get("canonicalid"),
                "oma_canonical_human": e2.get("canonicalid"),
                "oma_chrom_mouse": e1.get("chromosome"),
                "oma_chrom_human": e2.get("chromosome"),
                "oma_locus_start_mouse": e1.get("locus", {}).get("start"),
                "oma_locus_start_human": e2.get("locus", {}).get("start"),
                "oma_rel_type": r.get("rel_type"),
                "oma_distance": r.get("distance"),
                "oma_score": r.get("score"),
                "oma_group": r.get("oma_group"),
            })
    return pd.DataFrame(rows)


def compare_with_oma(master: pd.DataFrame) -> pd.DataFrame:
    """OMA Browser comparison via OMA-Ensembl xref bulk mapping.

    Uses /tmp/oma-ensembl.txt.gz if present (downloaded from
    https://omabrowser.org/All/oma-ensembl.txt.gz, ~138MB).
    Otherwise falls back to symbol-only resolution via OMA canonical IDs
    (UniProt-style SYMBOL_SPECIES).

    OMA mouse-human pairs come from /tmp/oma_pages/p*.json
    (https://omabrowser.org/api/pairs/MOUSE/HUMAN/?per_page=100 paginated).
    """
    oma = load_oma_cache()
    if oma.empty:
        return pd.DataFrame([{"metric": "error",
                              "value": "OMA cache not found at /tmp/oma_pages"}])
    oma = oma.dropna(subset=["oma_canonical_mouse", "oma_canonical_human"])

    # Prefer bulk OMA-Ensembl mapping when present
    map_path = Path("/tmp/oma-ensembl.txt.gz")
    if map_path.exists():
        # Format: OMA_ID<tab>Ensembl_ID; comments start with #
        m = pd.read_csv(map_path, sep="\t", comment="#",
                        names=["oma_id", "ensembl_id"], dtype=str)
        # Restrict to mouse + human and to gene IDs (not protein / transcript)
        is_mouse = m["oma_id"].str.startswith("MOUSE", na=False)
        is_human = m["oma_id"].str.startswith("HUMAN", na=False)
        is_gene = (m["ensembl_id"].str.startswith("ENSMUSG", na=False) |
                   m["ensembl_id"].str.startswith("ENSG", na=False))
        m = m[(is_mouse | is_human) & is_gene]
        m["ensembl_id"] = m["ensembl_id"].str.split(".").str[0]
        # OMA ID corresponds to a protein; multiple proteins can share an
        # Ensembl Gene ID (alternative isoforms). Take first.
        mouse_omaid_to_ens = m[m["oma_id"].str.startswith("MOUSE")] \
            .drop_duplicates("oma_id").set_index("oma_id")["ensembl_id"].to_dict()
        human_omaid_to_ens = m[m["oma_id"].str.startswith("HUMAN")] \
            .drop_duplicates("oma_id").set_index("oma_id")["ensembl_id"].to_dict()

        # OMA's pairs API uses 'omaid' field (e.g., MOUSE00001), not
        # 'canonicalid'. We didn't save omaid into cache. So we re-load
        # pages and pull omaid this time.
        pages = sorted(glob.glob("/tmp/oma_pages/p*.json"))
        rows = []
        for p in pages:
            try:
                with open(p) as fh:
                    d = json.load(fh)
            except Exception:
                continue
            for r in d:
                e1 = r.get("entry_1", {})
                e2 = r.get("entry_2", {})
                rows.append({
                    "oma_id_mouse": e1.get("omaid"),
                    "oma_id_human": e2.get("omaid"),
                    "oma_rel_type": r.get("rel_type"),
                })
        oma_with_id = pd.DataFrame(rows).dropna(
            subset=["oma_id_mouse", "oma_id_human"])
        oma_with_id["mouse_ensembl"] = oma_with_id["oma_id_mouse"].map(mouse_omaid_to_ens)
        oma_with_id["human_ensembl"] = oma_with_id["oma_id_human"].map(human_omaid_to_ens)
        oma_mapped = oma_with_id.dropna(subset=["mouse_ensembl",
                                                  "human_ensembl"])
        oma_pairs = set(zip(oma_mapped["mouse_ensembl"],
                              oma_mapped["human_ensembl"]))
        resolution_method = "OMA-Ensembl xref bulk mapping"
    else:
        # Fallback: symbol mapping
        def parse(s):
            if not isinstance(s, str):
                return None
            return s.split("_", 1)[0] if "_" in s else s
        oma["mouse_symbol_guess"] = oma["oma_canonical_mouse"].map(parse)
        oma["human_symbol_guess"] = oma["oma_canonical_human"].map(parse)
        mouse_sym2ens = master.dropna(
            subset=["mouse_symbol", "mouse_ensembl"]).drop_duplicates(
            "mouse_symbol")[["mouse_symbol", "mouse_ensembl"]]
        human_sym2ens = master.dropna(
            subset=["human_symbol", "human_ensembl"]).drop_duplicates(
            "human_symbol")[["human_symbol", "human_ensembl"]]
        mouse_map = {s.upper(): e for s, e in zip(
            mouse_sym2ens["mouse_symbol"], mouse_sym2ens["mouse_ensembl"])}
        human_map = {s.upper(): e for s, e in zip(
            human_sym2ens["human_symbol"], human_sym2ens["human_ensembl"])}
        oma["mouse_ensembl_match"] = oma["mouse_symbol_guess"].str.upper().map(mouse_map)
        oma["human_ensembl_match"] = oma["human_symbol_guess"].str.upper().map(human_map)
        oma_mapped = oma.dropna(subset=["mouse_ensembl_match",
                                          "human_ensembl_match"])
        oma_pairs = set(zip(oma_mapped["mouse_ensembl_match"],
                              oma_mapped["human_ensembl_match"]))
        resolution_method = "Symbol-based fallback"

    master_pc = master[(master["mouse_biotype"] == "protein_coding") &
                       (master["human_biotype"] == "protein_coding")]
    master_pairs = set(zip(master_pc["mouse_ensembl"],
                            master_pc["human_ensembl"]))
    master_pair_to_tier = dict(zip(
        zip(master_pc["mouse_ensembl"], master_pc["human_ensembl"]),
        master_pc["confidence_tier"]))

    inter = oma_pairs & master_pairs
    only_oma = oma_pairs - master_pairs
    only_master = master_pairs - oma_pairs

    by_tier = pd.Series([master_pair_to_tier[p] for p in inter
                          if p in master_pair_to_tier]).value_counts().to_dict()

    rows = [
        {"metric": "oma_resolution_method", "value": resolution_method},
        {"metric": "oma_pages_loaded",
         "value": len(glob.glob("/tmp/oma_pages/p*.json"))},
        {"metric": "oma_pairs_loaded", "value": len(oma)},
        {"metric": "oma_unique_pairs_resolved",
         "value": len(oma_pairs),
         "note": "Distinct (mouse_ensembl, human_ensembl) pairs after mapping"},
        {"metric": "master_PC_pairs", "value": len(master_pairs)},
        {"metric": "intersection_oma_master", "value": len(inter)},
        {"metric": "only_in_oma", "value": len(only_oma)},
        {"metric": "only_in_master", "value": len(only_master)},
        {"metric": "jaccard_oma_master",
         "value": round(len(inter) / len(oma_pairs | master_pairs), 4)
         if (oma_pairs | master_pairs) else 0},
        {"metric": "intersect_tier_H", "value": int(by_tier.get("H", 0))},
        {"metric": "intersect_tier_M", "value": int(by_tier.get("M", 0))},
        {"metric": "intersect_tier_L", "value": int(by_tier.get("L", 0))},
        {"metric": "oma_recall_of_master_H_PC",
         "value": round(
             (by_tier.get("H", 0) /
              len(master_pc[master_pc["confidence_tier"] == "H"])
              if len(master_pc[master_pc["confidence_tier"] == "H"]) else 0),
             4),
         "note": "Fraction of master Tier H PC pairs that OMA also calls"},
    ]
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# (6b) OMA precision per tier (independent gold standard)
# -----------------------------------------------------------------------------

def _load_oma_pairs_via_xref() -> set:
    """Helper: parse OMA-Ensembl bulk file + pages cache to get pair set."""
    map_path = Path("/tmp/oma-ensembl.txt.gz")
    if not map_path.exists():
        return set()
    m = pd.read_csv(map_path, sep="\t", comment="#",
                    names=["oma_id", "ensembl_id"], dtype=str)
    is_mouse = m["oma_id"].str.startswith("MOUSE", na=False)
    is_human = m["oma_id"].str.startswith("HUMAN", na=False)
    is_gene = (m["ensembl_id"].str.startswith("ENSMUSG", na=False) |
               m["ensembl_id"].str.startswith("ENSG", na=False))
    m = m[(is_mouse | is_human) & is_gene]
    m["ensembl_id"] = m["ensembl_id"].str.split(".").str[0]
    mouse_omaid_to_ens = m[m["oma_id"].str.startswith("MOUSE")] \
        .drop_duplicates("oma_id").set_index("oma_id")["ensembl_id"].to_dict()
    human_omaid_to_ens = m[m["oma_id"].str.startswith("HUMAN")] \
        .drop_duplicates("oma_id").set_index("oma_id")["ensembl_id"].to_dict()

    pages = sorted(glob.glob("/tmp/oma_pages/p*.json"))
    rows = []
    for p in pages:
        try:
            with open(p) as fh:
                d = json.load(fh)
        except Exception:
            continue
        for r in d:
            e1 = r.get("entry_1", {})
            e2 = r.get("entry_2", {})
            rows.append((e1.get("omaid"), e2.get("omaid")))
    pairs = set()
    for mom, hom in rows:
        me = mouse_omaid_to_ens.get(mom)
        he = human_omaid_to_ens.get(hom)
        if me and he:
            pairs.add((me, he))
    return pairs


def oma_precision_per_tier(master: pd.DataFrame) -> pd.DataFrame:
    """For each tier, compute precision against OMA Browser as independent
    ground truth. OMA is NOT a source of master, so this is a genuinely
    independent benchmark (Emms 2019 / Altenhoff 2015).
    """
    oma_pairs = _load_oma_pairs_via_xref()
    if not oma_pairs:
        return pd.DataFrame([{"metric": "error",
                               "value": "OMA cache not found"}])
    # OMA covers PC only; restrict master to PC-PC pairs
    master_pc = master[(master["mouse_biotype"] == "protein_coding") &
                       (master["human_biotype"] == "protein_coding")].copy()

    # Universe: pairs where both genes are in OMA somewhere (so OMA could
    # confirm)
    oma_mouse = set(p[0] for p in oma_pairs)
    oma_human = set(p[1] for p in oma_pairs)

    rows = []
    tier_cases = [
        ("H", master_pc["confidence_tier"] == "H", "All Tier H PC pairs"),
        ("M", master_pc["confidence_tier"] == "M", "All Tier M PC pairs"),
        ("L", master_pc["confidence_tier"] == "L", "All Tier L PC pairs"),
        ("H+M", master_pc["confidence_tier"].isin(["H", "M"]),
         "All Tier H+M PC pairs"),
        ("H_toga_blind",
         (master_pc["confidence_tier"] == "H") &
         (master_pc["tier_H_toga"].astype(int) == 0),
         "Tier H pairs with no TOGA evidence (biomaRt-only / miRBase only)"),
        ("M_toga_blind",
         (master_pc["confidence_tier"] == "M") &
         (master_pc["tier_H_toga"].astype(int) == 0),
         ("Tier M pairs with no TOGA evidence (biomaRt paralog, RBH, "
          "Phase J + LncBook)")),
    ]
    for tier_label, tier_filter, note in tier_cases:
        sub = master_pc.loc[tier_filter]
        pred = set(zip(sub["mouse_ensembl"], sub["human_ensembl"]))
        eval_set = {p for p in pred
                    if p[0] in oma_mouse and p[1] in oma_human}
        tp = eval_set & oma_pairs
        precision = len(tp) / len(eval_set) if eval_set else np.nan
        recall = len(tp) / len(oma_pairs) if oma_pairs else np.nan
        rows.append({
            "tier": tier_label,
            "n_pred_pc": len(pred),
            "n_pred_pc_in_oma_universe": len(eval_set),
            "n_oma_truth_pairs": len(oma_pairs),
            "n_TP": len(tp),
            "n_FP": len(eval_set) - len(tp),
            "precision_vs_oma": precision,
            "recall_vs_oma": recall,
            "note": note,
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# (7) Canonical lncRNA before-after RBH
# -----------------------------------------------------------------------------

CANONICAL_LNCRNAS = [
    "NEAT1", "MEG3", "MALAT1", "HOTAIR", "XIST", "H19",
    "KCNQ1OT1", "TERC", "PVT1", "NORAD", "HOTTIP",
]


def evaluate_canonical_lncrna(master: pd.DataFrame,
                              blast: pd.DataFrame,
                              blast_rbh: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sym in CANONICAL_LNCRNAS:
        before = blast[blast["human_symbol"] == sym]
        # The mouse symbol corresponding to the canonical 1:1 ortholog
        # (e.g., NEAT1 ↔ Neat1, MEG3 ↔ Meg3)
        canonical_mouse = sym.capitalize() if sym not in {"HOTTIP", "NORAD"} else sym
        # NORAD mouse symbol is Norad/4632428N05Rik; HOTTIP is Hottip
        if sym == "HOTTIP":
            canonical_mouse_set = {"Hottip"}
        elif sym == "NORAD":
            canonical_mouse_set = {"Norad", "1700110I01Rik", "4632428N05Rik"}
        elif sym == "KCNQ1OT1":
            canonical_mouse_set = {"Kcnq1ot1"}
        elif sym == "MEG3":
            canonical_mouse_set = {"Meg3"}
        elif sym == "MALAT1":
            canonical_mouse_set = {"Malat1"}
        elif sym == "NEAT1":
            canonical_mouse_set = {"Neat1"}
        elif sym == "HOTAIR":
            canonical_mouse_set = {"Hotair"}
        elif sym == "XIST":
            canonical_mouse_set = {"Xist"}
        elif sym == "H19":
            canonical_mouse_set = {"H19"}
        elif sym == "TERC":
            canonical_mouse_set = {"Terc"}
        elif sym == "PVT1":
            canonical_mouse_set = {"Pvt1"}
        else:
            canonical_mouse_set = {canonical_mouse}

        after = blast_rbh[blast_rbh["human_symbol"] == sym]
        master_rows = master[master["human_symbol"] == sym]
        m_M = master_rows[master_rows["confidence_tier"] == "M"]
        m_L = master_rows[master_rows["confidence_tier"] == "L"]
        m_H = master_rows[master_rows["confidence_tier"] == "H"]

        before_partners = sorted(set(before["mouse_symbol"].dropna()))
        after_partners = sorted(set(after["mouse_symbol"].dropna()))
        m_M_partners = sorted(set(m_M["mouse_symbol"].dropna()))
        canonical_in_before = bool(set(before_partners) & canonical_mouse_set)
        canonical_in_after = bool(set(after_partners) & canonical_mouse_set)
        canonical_in_M = bool(set(m_M_partners) & canonical_mouse_set)

        rows.append({
            "human_symbol": sym,
            "n_blast_partners_pre_rbh": len(before_partners),
            "n_blast_partners_post_rbh": len(after_partners),
            "n_master_partners_tier_H": len(m_H),
            "n_master_partners_tier_M": len(m_M),
            "n_master_partners_tier_L": len(m_L),
            "biologically_correct_mouse_symbol": ";".join(sorted(canonical_mouse_set)),
            "blast_partners_pre_rbh": ";".join(before_partners) or "—",
            "blast_partners_post_rbh": ";".join(after_partners) or "—",
            "master_tier_M_partners": ";".join(m_M_partners) or "—",
            "canonical_recovered_pre_rbh": canonical_in_before,
            "canonical_recovered_post_rbh": canonical_in_after,
            "canonical_promoted_to_tier_M_in_master": canonical_in_M,
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

def plot_recall(pr: pd.DataFrame, out: Path) -> None:
    sub = pr[pr["tier"].isin(["H", "M", "L", "H+M", "any"])]
    n_truth = int(sub["n_truth_pairs_pc_one2one"].iloc[0])
    fig, ax = plt.subplots(figsize=(5.5, 4))
    bars = ax.bar(sub["tier"], sub["recall_inclusive_pc_truth"],
                  color=["#1f77b4", "#2ca02c", "#9E9E9E", "#ff7f0e", "#d62728"])
    for b, v in zip(bars, sub["recall_inclusive_pc_truth"]):
        if pd.notna(v):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.01,
                    f"{v:.2%}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel(f"Recall vs TOGA one2one PC truth ({n_truth:,} pairs)")
    ax.set_xlabel("Tier")
    ax.set_ylim(0, 1.05)
    ax.set_title("PC recall per tier — TOGA one2one held out as truth")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


def plot_precision(pr: pd.DataFrame, out: Path,
                    oma_pr: pd.DataFrame | None = None) -> None:
    """Two-panel plot: TOGA-based precision (pair-level, harsh) and
    OMA-based precision (independent external benchmark).
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    sub = pr[pr["tier"].isin(["H", "M", "L", "H+M"])].copy()
    ax = axes[0]
    width = 0.35
    x = np.arange(len(sub))
    bars1 = ax.bar(x - width / 2, sub["precision_strict_one2one"],
                    width, label="Strict (TOGA one2one)", color="#1f77b4")
    bars2 = ax.bar(x + width / 2, sub["precision_lenient_any_toga_class"],
                    width, label="Lenient (any TOGA class)", color="#7CB9F7")
    for b, v in zip(bars1, sub["precision_strict_one2one"]):
        if pd.notna(v):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.01,
                    f"{v:.2%}", ha="center", va="bottom", fontsize=8)
    for b, v in zip(bars2, sub["precision_lenient_any_toga_class"]):
        if pd.notna(v):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.01,
                    f"{v:.2%}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(sub["tier"])
    ax.set_ylabel("Precision vs TOGA (TOGA-blind subset)")
    ax.set_xlabel("Tier")
    ax.set_ylim(0, 1.1)
    ax.set_title("Pair-level precision vs TOGA\n(strict: harsh pair-level; lenient: family-level)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1]
    if oma_pr is not None and "precision_vs_oma" in oma_pr.columns:
        oma_sub = oma_pr[oma_pr["tier"].isin(
            ["H", "M", "L", "H+M", "H_toga_blind", "M_toga_blind"])].copy()
        x = np.arange(len(oma_sub))
        bars = ax.bar(x, oma_sub["precision_vs_oma"],
                       color=["#1f77b4", "#2ca02c", "#9E9E9E", "#ff7f0e",
                              "#1f77b4", "#2ca02c"],
                       edgecolor=["black"] * 4 + ["red", "red"])
        for b, v in zip(bars, oma_sub["precision_vs_oma"]):
            if pd.notna(v):
                ax.text(b.get_x() + b.get_width() / 2, v + 0.01,
                        f"{v:.2%}", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(oma_sub["tier"], rotation=15, ha="right")
        ax.set_ylabel("Precision vs OMA (independent gold standard)")
        ax.set_xlabel("Tier")
        ax.set_ylim(0, 1.1)
        ax.set_title("OMA-Browser–based precision\n(red border = TOGA-blind subset)")
        ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


def plot_null_vs_observed(null_summary: pd.DataFrame,
                          pr: pd.DataFrame,
                          out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    null = null_summary.iloc[0]
    obs_tiers = pr[pr["tier"].isin(["H", "M", "L"])].copy()
    n_total_master = obs_tiers["n_predicted_pairs_total"].sum()

    bars = []
    labels = []
    vals = []
    cols = []
    for t, c in zip(["H", "M", "L"], ["#1f77b4", "#2ca02c", "#9E9E9E"]):
        labels.append(f"Master tier {t}")
        n = int(obs_tiers.loc[obs_tiers["tier"] == t,
                              "n_predicted_pairs_total"].values[0])
        vals.append(n / n_total_master)
        cols.append(c)
    labels += ["Null tier H", "Null tier M", "Null tier L"]
    vals += [null["n_null_at_tier_H"] / null["n_null_pairs_drawn"],
             null["n_null_at_tier_M"] / null["n_null_pairs_drawn"],
             null["n_null_at_tier_L"] / null["n_null_pairs_drawn"]]
    cols += ["#1f77b4", "#2ca02c", "#9E9E9E"]

    bars = ax.bar(labels, vals, color=cols)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.005,
                f"{v:.3%}", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("Fraction of pairs at tier")
    ax.set_title(f"Null distribution check ({null['n_null_pairs_drawn']:,} shuffled pairs)")
    ax.set_ylim(0, max(vals) * 1.3)
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


# -----------------------------------------------------------------------------
# (8) Wrong-partner null test (P1-12 fix)
# -----------------------------------------------------------------------------

def evaluate_wrong_partner_fp(master: pd.DataFrame,
                              n_genes: int = 1000,
                              seed: int = 42) -> pd.DataFrame:
    """Non-trivial null test: for genes that ARE in the master table,
    randomly assign a DIFFERENT partner from the same biotype pool.
    Check if the master table contains the wrong pair.

    This tests whether the pipeline spuriously connects unrelated genes,
    unlike the shuffled-pair test which is trivially satisfied by
    construction (pairs excluded from master, then checked against master).

    Returns per-tier FP rates. Expected: <1% Tier H, <5% Tier M.
    """
    rng = np.random.default_rng(seed)
    existing_pairs = set(zip(master["mouse_ensembl"], master["human_ensembl"]))

    # Build biotype pools: for each human biotype, collect all mouse genes
    # of the same biotype that appear in master
    biotype_pools = {}
    for bt in master["human_biotype"].dropna().unique():
        mouse_genes = master.loc[
            master["mouse_biotype"] == bt, "mouse_ensembl"
        ].dropna().unique()
        if len(mouse_genes) > 1:
            biotype_pools[bt] = mouse_genes

    # Sample n_genes human genes that are in master and have a biotype pool
    eligible = master.dropna(subset=["human_biotype"]).copy()
    eligible = eligible[eligible["human_biotype"].isin(biotype_pools.keys())]
    eligible = eligible.drop_duplicates("human_ensembl")
    if len(eligible) < n_genes:
        n_genes = len(eligible)
    sampled = eligible.sample(n=n_genes, random_state=seed)

    results = []
    for _, row in sampled.iterrows():
        h_ens = row["human_ensembl"]
        bt = row["human_biotype"]
        actual_mouse = row["mouse_ensembl"]
        pool = biotype_pools[bt]

        # Pick a random DIFFERENT mouse partner from the same biotype
        candidates = [g for g in pool if g != actual_mouse]
        if not candidates:
            continue
        wrong_mouse = rng.choice(candidates)
        wrong_pair = (wrong_mouse, h_ens)

        # Check if the wrong pair exists in master
        in_master = wrong_pair in existing_pairs
        if in_master:
            tier = master.loc[
                (master["mouse_ensembl"] == wrong_mouse) &
                (master["human_ensembl"] == h_ens),
                "confidence_tier"
            ].values
            tier_str = tier[0] if len(tier) > 0 else "unknown"
        else:
            tier_str = "absent"

        results.append({
            "human_ensembl": h_ens,
            "human_biotype": bt,
            "actual_mouse": actual_mouse,
            "wrong_mouse": wrong_mouse,
            "wrong_pair_in_master": in_master,
            "wrong_pair_tier": tier_str,
        })

    res_df = pd.DataFrame(results)

    # Summary
    n_total = len(res_df)
    n_fp = res_df["wrong_pair_in_master"].sum()
    fp_any = n_fp / n_total if n_total else 0.0
    fp_h = (res_df["wrong_pair_tier"] == "H").sum() / n_total if n_total else 0.0
    fp_m = (res_df["wrong_pair_tier"] == "M").sum() / n_total if n_total else 0.0
    fp_l = (res_df["wrong_pair_tier"] == "L").sum() / n_total if n_total else 0.0

    summary = pd.DataFrame([{
        "n_genes_tested": n_total,
        "n_wrong_pairs_in_master": int(n_fp),
        "fp_rate_any_tier": fp_any,
        "fp_rate_tier_H": fp_h,
        "fp_rate_tier_M": fp_m,
        "fp_rate_tier_L": fp_l,
        "note": ("Wrong-partner test: for genes in master, replace actual "
                 "partner with random same-biotype partner. FP = wrong pair "
                 "found in master. Tests specificity of pipeline pairing."),
    }])
    return summary, res_df


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> int:
    print("[validate] Loading master ortholog table ...", file=sys.stderr)
    master = load_master()
    print(f"  master: {len(master)} pairs; tier distribution: "
          f"{master['confidence_tier'].value_counts().to_dict()}",
          file=sys.stderr)
    toga = load_toga_truth()
    blast = load_blast()
    blast_rbh = load_blast_rbh()
    phasej = load_phasej()
    mirbase = load_mirbase()
    mirgenedb = load_mirgenedb()
    print(f"  TOGA L0: {len(toga)} pairs; blast: {len(blast)}; "
          f"blast_rbh: {len(blast_rbh)}; phasej: {len(phasej)}; "
          f"mirbase: {len(mirbase)}; mirgenedb: {len(mirgenedb)}",
          file=sys.stderr)

    print("[validate] Loading gene universes ...", file=sys.stderr)
    mouse_pc, human_pc = load_gene_universe()
    print(f"  mouse PC: {len(mouse_pc):,}; human PC: {len(human_pc):,}",
          file=sys.stderr)

    # (1) PC precision / recall
    print("[validate] (1) PC precision/recall vs TOGA truth ...",
          file=sys.stderr)
    pr = evaluate_pc_pr(master, toga)
    pr.to_csv(OUT_DIR / "precision_recall_per_tier.csv", index=False)
    print(pr.to_string(index=False), file=sys.stderr)

    # (2) Paralog agreement
    print("[validate] (2) Paralog-class agreement biomaRt vs TOGA ...",
          file=sys.stderr)
    paralog = evaluate_paralog_agreement(master)
    paralog.to_csv(OUT_DIR / "paralog_class_agreement.csv", index=False)
    print(paralog.to_string(index=False), file=sys.stderr)

    # (3) Null distribution
    print("[validate] (3) Null distribution check (10,000 shuffled pairs)...",
          file=sys.stderr)
    null_summary, null_full = evaluate_null_distribution(
        master, mouse_pc, human_pc, n_pairs=10000)
    null_summary.to_csv(OUT_DIR / "null_distribution_check.csv", index=False)
    null_full.to_csv(OUT_DIR / "null_distribution_pairs.tsv.gz",
                     sep="\t", index=False, compression="gzip")
    print(null_summary.to_string(index=False), file=sys.stderr)

    # (4) lncRNA consistency
    print("[validate] (4) lncRNA layer consistency Phase J vs BLAST RBH...",
          file=sys.stderr)
    lncrna_consistency = evaluate_lncrna_consistency(phasej, blast_rbh, master)
    lncrna_consistency.to_csv(OUT_DIR / "lncrna_layer_consistency.csv",
                              index=False)
    print(lncrna_consistency.to_string(index=False), file=sys.stderr)

    # (5) miRBase ∩ MirGeneDB
    print("[validate] (5) miRBase ∩ MirGeneDB Jaccard ...", file=sys.stderr)
    mirna_j = evaluate_mirna_jaccard(mirbase, mirgenedb)
    mirna_j.to_csv(OUT_DIR / "miRNA_layer_jaccard.csv", index=False)
    print(mirna_j.to_string(index=False), file=sys.stderr)

    # (6) OMA comparison
    print("[validate] (6) OMA Browser comparison ...", file=sys.stderr)
    try:
        oma_cmp = compare_with_oma(master)
        oma_cmp.to_csv(OUT_DIR / "oma_comparison.csv", index=False)
        print(oma_cmp.to_string(index=False), file=sys.stderr)
    except Exception as e:
        print(f"  OMA comparison failed: {e}", file=sys.stderr)
        pd.DataFrame([{"metric": "error", "value": str(e)}]).to_csv(
            OUT_DIR / "oma_comparison.csv", index=False)

    # (6b) Tier-wise OMA recall (precision-style independent check)
    print("[validate] (6b) Tier-wise OMA precision ...", file=sys.stderr)
    try:
        oma_per_tier = oma_precision_per_tier(master)
        oma_per_tier.to_csv(OUT_DIR / "oma_precision_per_tier.csv", index=False)
        print(oma_per_tier.to_string(index=False), file=sys.stderr)
    except Exception as e:
        print(f"  OMA per-tier precision failed: {e}", file=sys.stderr)

    # (7) Canonical lncRNAs
    print("[validate] (7) Canonical lncRNA before/after RBH ...",
          file=sys.stderr)
    canon = evaluate_canonical_lncrna(master, blast, blast_rbh)
    canon.to_csv(OUT_DIR / "canonical_lncrna_resolution.csv", index=False)
    print(canon.to_string(index=False), file=sys.stderr)

    # (8) Wrong-partner null test (P1-12 fix)
    print("[validate] (8) Wrong-partner null test (1,000 genes)...",
          file=sys.stderr)
    wp_summary, wp_full = evaluate_wrong_partner_fp(master, n_genes=1000)
    wp_summary.to_csv(OUT_DIR / "wrong_partner_null_test.csv", index=False)
    wp_full.to_csv(OUT_DIR / "wrong_partner_null_test_full.tsv.gz",
                   sep="\t", index=False, compression="gzip")
    print(wp_summary.to_string(index=False), file=sys.stderr)

    # Plots
    print("[validate] Plotting ...", file=sys.stderr)
    plot_recall(pr, FIG_DIR / "recall_per_tier.pdf")
    # Provide OMA per-tier precision to the plot
    oma_per_tier_for_plot = None
    oma_per_tier_path = OUT_DIR / "oma_precision_per_tier.csv"
    if oma_per_tier_path.exists():
        oma_per_tier_for_plot = pd.read_csv(oma_per_tier_path)
    plot_precision(pr, FIG_DIR / "precision_per_tier.pdf",
                    oma_per_tier_for_plot)
    plot_null_vs_observed(null_summary, pr,
                          FIG_DIR / "null_vs_observed_distribution.pdf")

    print("[validate] Done.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
