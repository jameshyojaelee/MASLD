#!/usr/bin/env python3
"""
149_patient_target_prioritization.py
Per-Patient Therapeutic Target Prioritization

Personalizes population-level driver scores (Script 126) for each patient using:
  - Patient's dominant transition assignment (Script 147 TAS)
  - Cell-type-specific expression (BayesPrism)
  - CCC L-R pair rewiring relevance
  - Subtype divergence and fate probability

Score formula:
  PersonalizedScore = PopulationScore * 0.50
                    + CelltypeActivity * 0.25
                    + CCC_Rewiring      * 0.15
                    + SubtypeDivergence  * 0.10
  (PathwayCoherence placeholder = 0, weight redistributed to PopulationScore)

Input:
  - results/progression/driver_scores.csv           (238K rows: gene x transition x evidence)
  - results/progression/therapeutic_roadmap.csv      (350 targets with druggability)
  - results/progression/transition_activity_scores.csv (per-patient TAS + dominant transitions)
  - results/progression/transition_subtype_assignments.csv (dominant transition + quiescent flag)
  - results/progression/cibersortx_celltype_expression/bayesprism_proportions.csv (1444 x 13)
  - results/progression/cibersortx_celltype_expression/bayesprism_{Hepatocyte,Stellate,Macrophage,Endothelial}.csv.gz
  - results/progression/ccc_transition_de.csv        (L-R pair rewiring per transition)
  - results/progression/divergence_genes.csv         (S1/S2 divergence)
  - results/progression/fate_probabilities.csv       (S1/S2 + P(F4))
  - results/staging_classifier/modeling_metadata.csv (clinical metadata)
  - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (druggability)

Output (to results/progression/):
  - patient_top10_targets.csv       (~14,440 rows: 1444 patients x 10 targets)
  - patient_pathway_class.csv       (1444 rows: dominant therapeutic modality per patient)
  - patient_drug_recommendations.csv (druggable targets per patient with drug info)
  - personalization_summary.csv     (validation metrics)

SLURM: io partition, 4 CPUs, 16G RAM, 48h
Env:   micromamba activate spatial
"""
import os
import sys
import logging
import numpy as np
import pandas as pd
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# -- Paths -----------------------------------------------------------------
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INTEG = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PROG = os.path.join(INTEG, "results/progression")
STAGING = os.path.join(INTEG, "results/staging_classifier")
RESULTS = os.path.join(BASE, "RNA-seq/results")
OUTDIR = PROG
CELLTYPE_DIR = os.path.join(PROG, "cibersortx_celltype_expression")
os.makedirs(OUTDIR, exist_ok=True)

# -- Weights ---------------------------------------------------------------
W_POP = 0.50       # Population score (0.40 + 0.10 from PathwayCoherence)
W_CELLTYPE = 0.25  # Cell-type activity
W_CCC = 0.15       # CCC rewiring
W_SUBTYPE = 0.10   # Subtype divergence

# -- Transition name mapping -----------------------------------------------
# TAS columns -> driver_scores transition names
TAS_TO_DS = {
    "TAS_F0_to_F1": "F0_to_F1",
    "TAS_F1_to_F2": "F1_to_F2",
    "TAS_F2_to_F3": "F2_to_F3",
    "TAS_F3_to_F4": "F3_to_F4",
    "TAS_NAS01_to_NAS24": "NAS_0_1_to_NAS_2_4",
    "TAS_NAS24_to_NAS5": "NAS_2_4_to_NAS_5",
    "TAS_NAS5_to_NAS68": "NAS_5_to_NAS_6_8",
}

# Fibrosis stage -> default transition for quiescent/unassigned patients
FIB_STAGE_TO_TRANSITION = {
    0: "F0_to_F1",
    1: "F1_to_F2",
    2: "F2_to_F3",
    3: "F3_to_F4",
    4: "F3_to_F4",   # F4 patients assigned to most advanced fibrosis transition
}

# Key cell types for weighted expression (BayesPrism high-confidence)
KEY_CELLTYPES = ["Hepatocyte", "Stellate", "Macrophage", "Endothelial"]

# -- Therapeutic modality gene patterns ------------------------------------
MODALITY_PATTERNS = {
    "anti_fibrotic": [
        "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL6A1",
        "COL6A2", "COL6A3", "COL14A1", "COL15A1", "COL18A1",
        "LOXL1", "LOXL2", "LOXL4", "LOX",
        "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
        "ACTA2", "FAP", "PDGFRA", "PDGFRB", "TIMP1", "TIMP2",
        "FN1", "SPARC", "THBS1", "THBS2", "LTBP2", "MMP2", "MMP14",
    ],
    "anti_inflammatory": [
        "CCL2", "CCL3", "CCL4", "CCL5", "CCL20", "CCL21",
        "CXCL1", "CXCL2", "CXCL5", "CXCL6", "CXCL8", "CXCL9", "CXCL10", "CXCL12",
        "TNF", "TNFRSF1A", "TNFRSF1B",
        "IL1B", "IL6", "IL6R", "IL18", "IL32", "IL33",
        "NFKB1", "NFKB2", "RELA", "RELB",
        "TLR2", "TLR4", "NLRP3", "CASP1",
    ],
    "metabolic": [
        "CYP2E1", "CYP1A2", "CYP3A4", "CYP7A1", "CYP8B1",
        "SLC2A2", "SLC27A5", "SLC10A1", "SLC22A1",
        "PPARA", "PPARG", "PPARGC1A",
        "FABP1", "FABP4", "FABP5",
        "FASN", "SCD", "ACACA", "DGAT2", "PNPLA3",
        "NR1H4", "NR1H3", "THRB", "GLP1R",
    ],
    "vascular": [
        "VEGFA", "VEGFB", "VEGFC", "KDR", "FLT1",
        "PECAM1", "CDH5", "VWF", "PLVAP", "ACKR1",
        "ANGPT1", "ANGPT2", "TEK", "NOS3",
        "ENG", "MCAM", "ICAM1",
    ],
}

# Invert to gene->modality lookup
GENE_TO_MODALITY = {}
for modality, genes in MODALITY_PATTERNS.items():
    for g in genes:
        GENE_TO_MODALITY[g] = modality

# Prefix-based matching for genes not in the explicit list
MODALITY_PREFIXES = {
    "anti_fibrotic": ["COL", "LOXL", "TGFB"],
    "anti_inflammatory": ["CCL", "CXCL", "TNF", "IL", "NFKB"],
    "metabolic": ["CYP", "SLC", "PPAR", "FABP"],
    "vascular": ["VEGF", "PECAM", "CDH5"],
}

TOP_N_CANDIDATES = 200  # Population-level candidates per transition
TOP_K_OUTPUT = 10       # Personalized targets per patient


def assign_modality(gene_symbol):
    """Classify a gene into a therapeutic modality."""
    if gene_symbol in GENE_TO_MODALITY:
        return GENE_TO_MODALITY[gene_symbol]
    for modality, prefixes in MODALITY_PREFIXES.items():
        for prefix in prefixes:
            if gene_symbol.startswith(prefix):
                return modality
    return "other"


def load_driver_scores():
    """Load population-level driver scores and extract top 200 per transition."""
    log.info("Loading driver_scores.csv ...")
    df = pd.read_csv(os.path.join(PROG, "driver_scores.csv"), low_memory=False)
    log.info(f"  Loaded {len(df):,} rows, {df['transition'].nunique()} transitions")

    # Filter to fibrosis transitions only (fib_dominant uses these)
    fib_transitions = ["F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"]
    all_transitions = df["transition"].unique().tolist()
    log.info(f"  Available transitions: {all_transitions}")

    # Top 200 per transition by driver_rank (rank 1 = best)
    candidates = {}
    for trans in all_transitions:
        sub = df[df["transition"] == trans].copy()
        sub = sub.dropna(subset=["driver_rank", "driver_score"])
        sub = sub.sort_values("driver_rank", ascending=True).head(TOP_N_CANDIDATES)
        candidates[trans] = sub
        log.info(f"  {trans}: {len(sub)} candidates (top score={sub['driver_score'].max():.4f})")

    return df, candidates


def load_patient_assignments():
    """Load patient transition assignments and clinical metadata."""
    log.info("Loading patient assignments ...")
    subtype_df = pd.read_csv(os.path.join(PROG, "transition_subtype_assignments.csv"))
    tas_df = pd.read_csv(os.path.join(PROG, "transition_activity_scores.csv"))
    meta_df = pd.read_csv(os.path.join(STAGING, "modeling_metadata.csv"))

    # Merge
    patients = subtype_df.merge(tas_df[["sample_id"]], on="sample_id", how="left")

    log.info(f"  {len(patients)} patients loaded")
    log.info(f"  fib_dominant distribution: {patients['fib_dominant'].value_counts().to_dict()}")
    return patients, meta_df


def assign_patient_transition(row):
    """Determine which transition to use for a patient."""
    fib_dom = row.get("fib_dominant", "")
    is_quiescent = str(row.get("is_quiescent", "FALSE")).upper() == "TRUE"

    # If patient has a valid fibrosis transition assignment, use it
    valid_transitions = {"F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"}
    if pd.notna(fib_dom) and fib_dom in valid_transitions:
        return fib_dom

    # Quiescent or unassigned: use fibrosis stage to pick default
    fib_stage = row.get("fibrosis_stage", np.nan)
    if pd.notna(fib_stage):
        try:
            fib_int = int(float(fib_stage))
            return FIB_STAGE_TO_TRANSITION.get(fib_int, "F0_to_F1")
        except (ValueError, TypeError):
            pass

    # Absolute fallback: F0_to_F1 (earliest transition)
    return "F0_to_F1"


def load_celltype_data():
    """Load BayesPrism proportions and cell-type-specific expression for key types."""
    log.info("Loading cell-type data ...")

    # Proportions (1444 x 13 cell types)
    proportions = pd.read_csv(os.path.join(CELLTYPE_DIR, "bayesprism_proportions.csv"))
    proportions = proportions.set_index("sample_id")
    log.info(f"  Proportions: {proportions.shape}")

    # Cell-type-specific expression for 4 key types
    ct_expr = {}
    for ct in KEY_CELLTYPES:
        fpath = os.path.join(CELLTYPE_DIR, f"bayesprism_{ct}.csv.gz")
        if os.path.exists(fpath):
            df = pd.read_csv(fpath, index_col=0)
            ct_expr[ct] = df
            log.info(f"  {ct}: {df.shape[0]} samples x {df.shape[1]} genes")
        else:
            log.warning(f"  {ct} expression file not found: {fpath}")

    return proportions, ct_expr


def load_ccc_data():
    """Load CCC transition DE and filter for significance."""
    log.info("Loading CCC rewiring data ...")
    ccc = pd.read_csv(os.path.join(PROG, "ccc_transition_de.csv"))
    log.info(f"  Raw CCC: {len(ccc):,} L-R pairs")

    # Filter: padj < 0.1 and |lfc| > 0.25
    ccc = ccc.dropna(subset=["padj", "lfc"])
    ccc_sig = ccc[(ccc["padj"] < 0.1) & (ccc["lfc"].abs() > 0.25)].copy()
    log.info(f"  Significant CCC: {len(ccc_sig):,} L-R pairs")

    # Build lookup: transition -> set of genes involved in rewired pairs
    ccc_genes = defaultdict(lambda: defaultdict(int))
    for _, row in ccc_sig.iterrows():
        trans = row["transition"]
        lig = row["ligand"]
        rec = row["receptor"]
        ccc_genes[trans][lig] += 1
        ccc_genes[trans][rec] += 1

    return ccc_sig, ccc_genes


def load_divergence_data():
    """Load S1/S2 subtype divergence genes."""
    log.info("Loading divergence genes ...")
    div = pd.read_csv(os.path.join(PROG, "divergence_genes.csv"))
    log.info(f"  {len(div)} divergence genes")

    # Use gene_symbol column; significant genes by padj
    div_sig = div[div["padj"] < 0.05].copy()
    div_genes = set(div_sig["gene_symbol"].dropna().unique())
    log.info(f"  {len(div_genes)} significant divergence gene symbols")

    # Build Cohen's d lookup for directionality
    div_d = dict(zip(div_sig["gene_symbol"], div_sig["cohens_d"]))
    return div_genes, div_d


def load_fate_probabilities():
    """Load fate probabilities (P(F4) per patient)."""
    log.info("Loading fate probabilities ...")
    fate = pd.read_csv(os.path.join(PROG, "fate_probabilities.csv"))
    log.info(f"  {len(fate)} patients with fate probabilities")
    fate = fate.set_index("sample_id")

    # Compute median P(F4) for thresholding
    if "fate_prob_F4" in fate.columns:
        median_f4 = fate["fate_prob_F4"].median()
        log.info(f"  Median P(F4) = {median_f4:.4f}")
    else:
        median_f4 = 0.5

    return fate, median_f4


def load_druggability():
    """Load druggability from multi-evidence atlas and therapeutic roadmap."""
    log.info("Loading druggability data ...")

    # Atlas for dgidb_druggable and opentargets_drug
    atlas = pd.read_csv(
        os.path.join(RESULTS, "multi_evidence/multi_evidence_atlas.csv"),
        usecols=["human_symbol", "dgidb_druggable", "opentargets_drug", "lincs_reversal"],
        low_memory=False,
    )
    atlas = atlas.rename(columns={"human_symbol": "gene_symbol"})
    atlas["is_druggable_atlas"] = (
        atlas["dgidb_druggable"].astype(str).str.lower().eq("true")
        | atlas["opentargets_drug"].astype(str).str.lower().eq("true")
    )
    log.info(f"  Atlas druggable: {atlas['is_druggable_atlas'].sum()}")

    # Roadmap has curated is_druggable per gene x transition
    roadmap = pd.read_csv(os.path.join(PROG, "therapeutic_roadmap.csv"))
    roadmap_druggable = set(
        roadmap.loc[roadmap["is_druggable"].astype(str).str.lower() == "true", "gene_symbol"]
    )
    log.info(f"  Roadmap druggable genes: {len(roadmap_druggable)}")

    # Merge: druggable if in either source
    all_druggable = set(atlas.loc[atlas["is_druggable_atlas"], "gene_symbol"]) | roadmap_druggable
    log.info(f"  Combined druggable: {len(all_druggable)}")

    return all_druggable, atlas


def compute_celltype_activity(sample_id, candidate_genes, proportions, ct_expr, driver_sub):
    """
    Compute cell-type weighted expression for a patient's candidate genes.

    For each gene: sum(proportion_ct * expression_ct_gene * sign(expected_direction))
    Expected direction = sign of driver_score's trans_logFC if available, else +1.
    Z-score across all candidates.
    """
    scores = {}

    # Get patient proportions for key cell types
    if sample_id not in proportions.index:
        return pd.Series(0.0, index=candidate_genes)

    props = {}
    for ct in KEY_CELLTYPES:
        if ct in proportions.columns:
            props[ct] = proportions.loc[sample_id, ct]
        else:
            props[ct] = 0.0

    # Get expected direction from trans_logFC in driver_scores
    direction_map = {}
    if "trans_logFC" in driver_sub.columns:
        for _, row in driver_sub.iterrows():
            gene = row["gene_symbol"]
            lfc = row.get("trans_logFC", np.nan)
            if pd.notna(lfc) and lfc != 0:
                direction_map[gene] = np.sign(lfc)
            else:
                direction_map[gene] = 1.0

    for gene in candidate_genes:
        weighted_expr = 0.0
        direction = direction_map.get(gene, 1.0)

        for ct in KEY_CELLTYPES:
            if ct not in ct_expr:
                continue
            ct_df = ct_expr[ct]
            if sample_id not in ct_df.index or gene not in ct_df.columns:
                continue
            expr_val = ct_df.loc[sample_id, gene]
            weighted_expr += props[ct] * expr_val * direction

        scores[gene] = weighted_expr

    result = pd.Series(scores)

    # Z-score across candidates (within this patient)
    if result.std() > 0:
        result = (result - result.mean()) / result.std()
    else:
        result = result * 0.0

    return result


def compute_ccc_rewiring_score(candidate_genes, transition, ccc_genes):
    """
    Score each candidate gene by its involvement in rewired L-R pairs at this transition.
    Normalized by maximum count.
    """
    trans_genes = ccc_genes.get(transition, {})
    scores = {}
    for gene in candidate_genes:
        scores[gene] = trans_genes.get(gene, 0)

    result = pd.Series(scores)
    max_val = result.max()
    if max_val > 0:
        result = result / max_val
    return result


def compute_subtype_divergence_score(
    candidate_genes, sample_id, div_genes, div_d, fate_df, median_f4
):
    """
    Score genes by subtype divergence relevance.
    If patient has high P(F4), upweight divergence genes.
    Uses Cohen's d as magnitude.
    """
    scores = {}

    # Get patient fate probability
    high_risk = False
    if sample_id in fate_df.index and "fate_prob_F4" in fate_df.columns:
        p_f4 = fate_df.loc[sample_id, "fate_prob_F4"]
        high_risk = p_f4 > median_f4

    risk_multiplier = 1.5 if high_risk else 1.0

    for gene in candidate_genes:
        if gene in div_genes:
            # Use absolute Cohen's d as magnitude, scaled by risk
            d_val = abs(div_d.get(gene, 0.5))
            scores[gene] = d_val * risk_multiplier
        else:
            scores[gene] = 0.0

    result = pd.Series(scores)
    max_val = result.max()
    if max_val > 0:
        result = result / max_val
    return result


def classify_patient_modality(top_genes):
    """Classify patient's dominant therapeutic modality from their top targets."""
    modality_counts = defaultdict(int)
    for gene in top_genes:
        mod = assign_modality(gene)
        modality_counts[mod] += 1

    if not modality_counts:
        return "unclassified"

    # Remove 'other' from consideration unless it's the only one
    if len(modality_counts) > 1 and "other" in modality_counts:
        non_other = {k: v for k, v in modality_counts.items() if k != "other"}
        if non_other:
            modality_counts = non_other

    dominant = max(modality_counts, key=modality_counts.get)
    return dominant


def main():
    log.info("=" * 70)
    log.info("Script 149: Per-Patient Therapeutic Target Prioritization")
    log.info("=" * 70)

    # -- Load all data -------------------------------------------------------
    full_drivers, candidates = load_driver_scores()
    patients, meta_df = load_patient_assignments()
    proportions, ct_expr = load_celltype_data()
    ccc_sig, ccc_genes = load_ccc_data()
    div_genes, div_d = load_divergence_data()
    fate_df, median_f4 = load_fate_probabilities()
    all_druggable, atlas_drug = load_druggability()

    # -- Build population score lookup: normalized driver_score per transition --
    pop_score_lookup = {}
    for trans, cands in candidates.items():
        scores = cands.set_index("gene_symbol")["driver_score"]
        # Min-max normalize to [0, 1]
        smin, smax = scores.min(), scores.max()
        if smax > smin:
            norm_scores = (scores - smin) / (smax - smin)
        else:
            norm_scores = scores * 0.0 + 0.5
        pop_score_lookup[trans] = norm_scores

    # -- Assign transitions to patients --------------------------------------
    patients["assigned_transition"] = patients.apply(assign_patient_transition, axis=1)
    trans_counts = patients["assigned_transition"].value_counts()
    log.info(f"Patient transition assignments:\n{trans_counts}")

    # -- Personalize per patient ---------------------------------------------
    all_top_targets = []
    all_pathway_classes = []
    all_drug_recs = []

    n_patients = len(patients)
    log.info(f"Personalizing targets for {n_patients} patients ...")

    for idx, (_, patient) in enumerate(patients.iterrows()):
        sample_id = patient["sample_id"]
        transition = patient["assigned_transition"]

        if idx % 200 == 0:
            log.info(f"  Processing patient {idx+1}/{n_patients} ({sample_id}) ...")

        # Get candidate genes for this transition
        if transition not in candidates:
            log.warning(f"  No candidates for transition {transition}, skipping {sample_id}")
            continue

        cand_df = candidates[transition]
        cand_genes = cand_df["gene_symbol"].tolist()

        if not cand_genes:
            continue

        # 1) Population score (normalized)
        pop_scores = pop_score_lookup.get(transition, pd.Series(dtype=float))
        pop_vec = pd.Series(
            [pop_scores.get(g, 0.0) for g in cand_genes], index=cand_genes
        )

        # 2) Cell-type activity
        ct_vec = compute_celltype_activity(
            sample_id, cand_genes, proportions, ct_expr, cand_df
        )
        # Clip extreme z-scores and rescale to [0,1]
        ct_vec = ct_vec.clip(-3, 3)
        ct_min, ct_max = ct_vec.min(), ct_vec.max()
        if ct_max > ct_min:
            ct_vec = (ct_vec - ct_min) / (ct_max - ct_min)
        else:
            ct_vec = ct_vec * 0.0 + 0.5

        # 3) CCC rewiring
        ccc_vec = compute_ccc_rewiring_score(cand_genes, transition, ccc_genes)

        # 4) Subtype divergence
        div_vec = compute_subtype_divergence_score(
            cand_genes, sample_id, div_genes, div_d, fate_df, median_f4
        )

        # -- Composite personalized score ------------------------------------
        personalized = (
            W_POP * pop_vec
            + W_CELLTYPE * ct_vec.reindex(cand_genes, fill_value=0.0)
            + W_CCC * ccc_vec.reindex(cand_genes, fill_value=0.0)
            + W_SUBTYPE * div_vec.reindex(cand_genes, fill_value=0.0)
        )

        # Top K
        top_k = personalized.nlargest(TOP_K_OUTPUT)

        for rank_i, (gene, score) in enumerate(top_k.items(), 1):
            all_top_targets.append({
                "sample_id": sample_id,
                "gene_symbol": gene,
                "personalized_score": round(score, 6),
                "personalized_rank": rank_i,
                "assigned_transition": transition,
                "pop_score": round(pop_vec.get(gene, 0.0), 6),
                "celltype_score": round(
                    ct_vec.reindex(cand_genes, fill_value=0.0).get(gene, 0.0), 6
                ),
                "ccc_score": round(ccc_vec.get(gene, 0.0), 6),
                "divergence_score": round(div_vec.get(gene, 0.0), 6),
                "is_druggable": gene in all_druggable,
                "modality": assign_modality(gene),
            })

        # -- Pathway classification ------------------------------------------
        top_gene_list = top_k.index.tolist()
        dominant_mod = classify_patient_modality(top_gene_list)

        # Count per modality in top 10
        mod_counts = defaultdict(int)
        for g in top_gene_list:
            mod_counts[assign_modality(g)] += 1

        all_pathway_classes.append({
            "sample_id": sample_id,
            "assigned_transition": transition,
            "dominant_modality": dominant_mod,
            "n_anti_fibrotic": mod_counts.get("anti_fibrotic", 0),
            "n_anti_inflammatory": mod_counts.get("anti_inflammatory", 0),
            "n_metabolic": mod_counts.get("metabolic", 0),
            "n_vascular": mod_counts.get("vascular", 0),
            "n_other": mod_counts.get("other", 0),
        })

        # -- Drug recommendations -------------------------------------------
        for gene in top_gene_list:
            if gene in all_druggable:
                all_drug_recs.append({
                    "sample_id": sample_id,
                    "gene_symbol": gene,
                    "personalized_score": round(personalized.get(gene, 0.0), 6),
                    "personalized_rank": top_gene_list.index(gene) + 1,
                    "assigned_transition": transition,
                    "modality": assign_modality(gene),
                })

    # -- Build output DataFrames ---------------------------------------------
    log.info("Building output tables ...")

    top_targets_df = pd.DataFrame(all_top_targets)
    pathway_class_df = pd.DataFrame(all_pathway_classes)
    drug_recs_df = pd.DataFrame(all_drug_recs)

    log.info(f"  patient_top10_targets: {len(top_targets_df)} rows")
    log.info(f"  patient_pathway_class: {len(pathway_class_df)} rows")
    log.info(f"  patient_drug_recommendations: {len(drug_recs_df)} rows")

    # -- Validation / summary metrics ----------------------------------------
    summary_rows = []

    # 1. Personalization diversity: how many unique genes appear in top-10 across patients?
    unique_genes = top_targets_df["gene_symbol"].nunique()
    total_genes_possible = len(
        set().union(*[set(c["gene_symbol"]) for c in candidates.values()])
    )
    summary_rows.append({
        "metric": "unique_top10_genes",
        "value": unique_genes,
        "description": "Unique genes appearing in any patient's top 10",
    })
    summary_rows.append({
        "metric": "personalization_diversity",
        "value": round(unique_genes / max(total_genes_possible, 1), 4),
        "description": "Fraction of candidate pool appearing in at least one top 10",
    })

    # 2. Mean personalized score
    mean_score = top_targets_df["personalized_score"].mean()
    summary_rows.append({
        "metric": "mean_personalized_score",
        "value": round(mean_score, 4),
        "description": "Average personalized score across all patient-target pairs",
    })

    # 3. Component contribution analysis
    if len(top_targets_df) > 0:
        for col, label in [
            ("pop_score", "population"), ("celltype_score", "celltype"),
            ("ccc_score", "ccc"), ("divergence_score", "divergence"),
        ]:
            mean_val = top_targets_df[col].mean()
            summary_rows.append({
                "metric": f"mean_{label}_component",
                "value": round(mean_val, 4),
                "description": f"Mean {label} component score in top 10 targets",
            })

    # 4. Druggable fraction in top 10
    if len(top_targets_df) > 0:
        druggable_frac = top_targets_df["is_druggable"].mean()
        summary_rows.append({
            "metric": "druggable_fraction_top10",
            "value": round(druggable_frac, 4),
            "description": "Fraction of top 10 targets that are druggable",
        })

    # 5. Modality distribution
    if len(pathway_class_df) > 0:
        mod_dist = pathway_class_df["dominant_modality"].value_counts(normalize=True)
        for mod, frac in mod_dist.items():
            summary_rows.append({
                "metric": f"modality_frac_{mod}",
                "value": round(frac, 4),
                "description": f"Fraction of patients with {mod} as dominant modality",
            })

    # 6. Patients with drug recommendations
    n_patients_with_drugs = drug_recs_df["sample_id"].nunique() if len(drug_recs_df) > 0 else 0
    summary_rows.append({
        "metric": "patients_with_drug_recs",
        "value": n_patients_with_drugs,
        "description": "Number of patients with at least one druggable target in top 10",
    })

    # 7. Per-transition target overlap (Jaccard between patients in same transition)
    for trans in patients["assigned_transition"].unique():
        trans_patients = top_targets_df[top_targets_df["assigned_transition"] == trans]
        if len(trans_patients) == 0:
            continue
        patient_gene_sets = (
            trans_patients.groupby("sample_id")["gene_symbol"].apply(set).tolist()
        )
        if len(patient_gene_sets) >= 2:
            jaccards = []
            # Sample pairwise Jaccard (max 500 pairs to keep it fast)
            import random
            random.seed(42)
            n_pairs = min(500, len(patient_gene_sets) * (len(patient_gene_sets) - 1) // 2)
            pairs_sampled = 0
            for i in range(len(patient_gene_sets)):
                for j in range(i + 1, len(patient_gene_sets)):
                    inter = len(patient_gene_sets[i] & patient_gene_sets[j])
                    union = len(patient_gene_sets[i] | patient_gene_sets[j])
                    if union > 0:
                        jaccards.append(inter / union)
                    pairs_sampled += 1
                    if pairs_sampled >= n_pairs:
                        break
                if pairs_sampled >= n_pairs:
                    break
            mean_jaccard = np.mean(jaccards) if jaccards else 0.0
            summary_rows.append({
                "metric": f"mean_jaccard_{trans}",
                "value": round(mean_jaccard, 4),
                "description": f"Mean pairwise Jaccard of top-10 targets within {trans}",
            })

    summary_df = pd.DataFrame(summary_rows)

    # -- Save outputs --------------------------------------------------------
    out_top10 = os.path.join(OUTDIR, "patient_top10_targets.csv")
    out_class = os.path.join(OUTDIR, "patient_pathway_class.csv")
    out_drugs = os.path.join(OUTDIR, "patient_drug_recommendations.csv")
    out_summary = os.path.join(OUTDIR, "personalization_summary.csv")

    top_targets_df.to_csv(out_top10, index=False)
    pathway_class_df.to_csv(out_class, index=False)
    drug_recs_df.to_csv(out_drugs, index=False)
    summary_df.to_csv(out_summary, index=False)

    log.info(f"Saved: {out_top10}")
    log.info(f"Saved: {out_class}")
    log.info(f"Saved: {out_drugs}")
    log.info(f"Saved: {out_summary}")

    # -- Print key metrics ---------------------------------------------------
    log.info("=" * 70)
    log.info("SUMMARY")
    log.info("=" * 70)
    for _, row in summary_df.iterrows():
        log.info(f"  {row['metric']}: {row['value']}")

    log.info("Script 149 complete.")


if __name__ == "__main__":
    main()
