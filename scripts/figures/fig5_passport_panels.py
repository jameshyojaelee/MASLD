#!/usr/bin/env python
# ─────────────────────────────────────────────────────────────────────────────
# Fig 5 — evidence-passport panels, STEP 1 of 2 (data extraction).
#
# KEY MESSAGE (F5): evidence CLASS, not a score, is what changes interpretation.
# The passport partitions the whole accepted gene universe by (i) what each
# assay actually CALLED and (ii) whether the gene was TESTABLE at all. It does
# not rank genes, does not sum modalities, and produces no probability.
#
# This script only RESHAPES the sealed PASS06 candidate into small tidy CSVs.
# It recomputes NO scientific call: every call_state / testability_state /
# provenance_state is carried through verbatim from the frozen parquet.
#
# Source (candidate, signed PASS06_VALIDATED, canonical_promotion_authorized=false):
#   RNA-seq/results/evidence_passports/candidates/
#     program-context-v2-candidate-2026-08-07/
#
# Env: micromamba activate spatial   (pyarrow)
# ─────────────────────────────────────────────────────────────────────────────
import hashlib
import os
import sys

import pandas as pd

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
PASSPORT_DIR = os.path.join(
    BASE,
    "RNA-seq/results/evidence_passports/candidates",
    "program-context-v2-candidate-2026-08-07",
)
OUT_DIR = os.path.join(BASE, "figures/main/fig5_convergence/data")
os.makedirs(OUT_DIR, exist_ok=True)

# The bundle is sealed; refuse to proceed if the seal is absent.
SEAL = os.path.join(PASSPORT_DIR, "PASS06_VALIDATED")
if not os.path.exists(SEAL):
    sys.exit(f"FATAL: passport seal not found at {SEAL}")
with open(SEAL, "rb") as fh:
    seal_bytes = fh.read()
print(f"[seal] PASS06_VALIDATED sha256 = {hashlib.sha256(seal_bytes).hexdigest()}")

ev = pd.read_parquet(os.path.join(PASSPORT_DIR, "passport_evidence_long.parquet"))
gi = pd.read_parquet(os.path.join(PASSPORT_DIR, "passport_gene_index.parquet"))
nx = pd.read_csv(
    os.path.join(PASSPORT_DIR, "passport_next_experiment.tsv"), sep="\t", dtype=str
)
rb = pd.read_csv(os.path.join(PASSPORT_DIR, "passport_rulebook.tsv"), sep="\t", dtype=str)

n_genes = gi["ensembl_id"].nunique()
print(f"[in] evidence rows {len(ev)}  genes {n_genes}  next-exp rows {len(nx)}")
assert len(ev) == 2 * n_genes, "expected exactly two gene-grain domains per gene"

# ── Guard rails: the passport must not have acquired a score or a negative ───
banned = [c for c in ev.columns if c.lower() in
          ("score", "rank", "priority", "posterior", "modality_count", "tier")]
assert not banned, f"scoring column present in frozen evidence: {banned}"
n_tested_negative = int((ev["call_state"] == "tested_negative").sum())
print(f"[check] tested_negative rows = {n_tested_negative} (must be 0 in this release)")

# ─────────────────────────────────────────────────────────────────────────────
# PANEL 1 — joint evidence-call matrix (genetics call x transcriptomics call).
# One cell per (genetics call, transcriptomics call). The `tested_negative`
# row and column are emitted EXPLICITLY with n = 0 so the absent negative is
# drawn rather than silently omitted.
# ─────────────────────────────────────────────────────────────────────────────
wide = ev.pivot_table(
    index="ensembl_id", columns="evidence_domain", values="call_state", aggfunc="first"
)
wide.columns = ["genetics_call", "transcriptomics_call"]
joint = gi.set_index("ensembl_id").join(wide)

CALL_LEVELS = ["supported", "indeterminate", "tested_negative", "untestable"]

grid = (
    joint.groupby(["genetics_call", "transcriptomics_call"], observed=True)
    .size()
    .rename("n_genes")
    .reset_index()
)
# class label per occupied cell (verbatim from the frozen gene index)
cls = (
    joint.groupby(["genetics_call", "transcriptomics_call"], observed=True)[
        "primary_evidence_class"
    ]
    .agg(lambda s: ";".join(sorted(set(s))))
    .rename("primary_evidence_class")
    .reset_index()
)
grid = grid.merge(cls, on=["genetics_call", "transcriptomics_call"], how="left")

full = pd.MultiIndex.from_product(
    [CALL_LEVELS, CALL_LEVELS], names=["genetics_call", "transcriptomics_call"]
).to_frame(index=False)
grid = full.merge(grid, on=["genetics_call", "transcriptomics_call"], how="left")
grid["n_genes"] = grid["n_genes"].fillna(0).astype(int)
grid["primary_evidence_class"] = grid["primary_evidence_class"].fillna("")

# An inferential result exists only where BOTH layers were testable.
grid["cell_kind"] = "inferential_result"
grid.loc[
    (grid.genetics_call == "untestable") | (grid.transcriptomics_call == "untestable"),
    "cell_kind",
] = "not_measurable"
grid.loc[
    (grid.genetics_call == "tested_negative")
    | (grid.transcriptomics_call == "tested_negative"),
    "cell_kind",
] = "informative_negative_absent"

assert grid.n_genes.sum() == n_genes, "cell counts must partition the gene universe"
grid.to_csv(os.path.join(OUT_DIR, "fig5_passport_call_matrix.csv"), index=False)
print("[out] fig5_passport_call_matrix.csv")
print(grid[grid.n_genes > 0].to_string(index=False))

# ── Per-domain margins: call + the machine-readable testability reason ───────
margin = (
    ev.groupby(["evidence_domain", "call_state", "testability_state", "testability_reason"])
    .size()
    .rename("n_genes")
    .reset_index()
)
margin["assay"] = margin["evidence_domain"].map(
    ev.drop_duplicates("evidence_domain").set_index("evidence_domain")["assay"]
)
margin.to_csv(os.path.join(OUT_DIR, "fig5_passport_domain_margin.csv"), index=False)
print("[out] fig5_passport_domain_margin.csv")
print(margin.to_string(index=False))

# ── Provenance strip: source dependence over all gene-grain evidence rows ────
prov = (
    ev.groupby(["evidence_domain", "provenance_state"]).size().rename("n_rows").reset_index()
)
prov.to_csv(os.path.join(OUT_DIR, "fig5_passport_provenance.csv"), index=False)
print("[out] fig5_passport_provenance.csv")
print(prov.to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# PANEL 2 — the next discriminating experiment.
# One deterministic rule per gene, carried verbatim from the frozen rulebook.
#
# ORDERING: purely descriptive, by gene count descending. An earlier draft
# ordered by evidence-class semantics with `concordant` first; that reads as a
# merit hierarchy, which this paper rejects. Frequency ordering carries no
# merit claim and puts the dominant fact -- how much was never measurable --
# where it belongs. This is a routing table, not a leaderboard.
# ─────────────────────────────────────────────────────────────────────────────
nx_counts = nx.groupby("experiment_rule_id").size().rename("n_genes").reset_index()
rb_exp = rb[rb.rule_type == "next_experiment"][
    [
        "rule_id",
        "biological_model",
        "context",
        "perturbation",
        "primary_readout",
        "falsifying_outcome",
    ]
].rename(columns={"rule_id": "experiment_rule_id"})
nx_tab = nx_counts.merge(rb_exp, on="experiment_rule_id", how="left")

class_of_rule = (
    gi.groupby("next_experiment_rule_id")["primary_evidence_class"]
    .agg(lambda s: ";".join(sorted(set(s))))
    .rename("primary_evidence_class")
    .reset_index()
    .rename(columns={"next_experiment_rule_id": "experiment_rule_id"})
)
nx_tab = nx_tab.merge(class_of_rule, on="experiment_rule_id", how="left")

nx_tab = nx_tab.sort_values(
    ["n_genes", "experiment_rule_id"], ascending=[False, True]
).reset_index(drop=True)
nx_tab["display_order_by_count"] = range(1, len(nx_tab) + 1)
assert nx_tab.biological_model.notna().all(), "unmapped experiment rule"
assert int(nx_tab.n_genes.sum()) == n_genes, "experiment routing must cover every gene"
nx_tab.to_csv(os.path.join(OUT_DIR, "fig5_passport_next_experiment.csv"), index=False)
print("[out] fig5_passport_next_experiment.csv")
print(nx_tab[["experiment_rule_id", "primary_evidence_class", "n_genes"]].to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# Hero-gene passport states — post-results ILLUSTRATIONS of class boundaries.
# Not a discovery set, not a validation set, not ranked.
# ─────────────────────────────────────────────────────────────────────────────
HEROES = ["THRB", "HKDC1", "GLP1R", "MTARC1"]
hero_ev = ev[ev.symbol.isin(HEROES)][
    [
        "symbol",
        "evidence_domain",
        "assay",
        "phenotype",
        "dataset_id",
        "effect_unit",
        "estimate",
        "q_value",
        "direction",
        "call_state",
        "testability_state",
        "testability_reason",
        "provenance_state",
    ]
].copy()
hero_cls = joint[joint.symbol.isin(HEROES)][
    ["symbol", "primary_evidence_class", "role_hypothesis", "next_experiment_rule_id"]
]
hero_ev = hero_ev.merge(hero_cls, on="symbol", how="left")
hero_ev.to_csv(os.path.join(OUT_DIR, "fig5_passport_hero_states.csv"), index=False)
print("[out] fig5_passport_hero_states.csv")
print(hero_ev.to_string(index=False))

# ── Provenance record for the caption ────────────────────────────────────────
with open(os.path.join(OUT_DIR, "fig5_passport_PROVENANCE.txt"), "w") as fh:
    fh.write("source_bundle\t%s\n" % PASSPORT_DIR)
    fh.write("seal_sha256\t%s\n" % hashlib.sha256(seal_bytes).hexdigest())
    fh.write("n_genes\t%d\n" % n_genes)
    fh.write("n_evidence_rows\t%d\n" % len(ev))
    fh.write("n_tested_negative_rows\t%d\n" % n_tested_negative)
    fh.write("gene_grain_domains\tgenetics;transcriptomics\n")
    fh.write("scientific_call_recomputed\tfalse\n")
print("[out] fig5_passport_PROVENANCE.txt")
print("[done]")
