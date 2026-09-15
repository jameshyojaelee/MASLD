#!/usr/bin/env python3
"""Arm loaders for the four-arm aspect-axis decomposition.

Each loader returns an Arm: a participant axis, an aspect table on that axis, a
nuisance table on that axis, and a feature matrix on that axis. Nothing is
pooled here. The arms are four different instruments and each is loaded to its
own native scale.
"""
from __future__ import annotations

import pathlib
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from aspect_axis_common import (BENCH, ROOT, JoinError, read_tsv_strict,
                                require_join, sha256_file)

SC_INPUTS = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs"
PSEUDOBULK = SC_INPUTS / "dialogue_pseudobulk"
SAF_LOCK = (BENCH / "executions/gse202379-saf-recovery-21183702/saf_lock")
DONOR_PAIRING = (BENCH / "drafts/showcase_substrate_dryrun_20260827T150253_ntasks1"
                 / "sources/gse202379_donor_pairing.csv")


@dataclass
class Arm:
    arm_id: str
    cohort: str
    assay: str
    participants: list
    aspects: pd.DataFrame          # participant x aspect, numeric, no NaN
    nuisance: pd.DataFrame         # participant x covariate, numeric, no NaN
    features: np.ndarray           # participant x feature
    feature_ids: list
    composite: str                 # name of the summed activity label
    components: list               # aspects that sum into the composite
    inputs: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    def check(self):
        n = len(self.participants)
        assert self.aspects.shape[0] == n, "aspect table off the participant axis"
        assert self.nuisance.shape[0] == n, "nuisance table off the participant axis"
        assert self.features.shape[0] == n, "feature matrix off the participant axis"
        assert not self.aspects.isna().any().any(), "an aspect carries NaN"
        assert not self.nuisance.isna().any().any(), "a covariate carries NaN"
        assert len(set(self.participants)) == n, "the participant axis repeats"
        assert np.isfinite(self.features).all(), "the feature matrix carries NaN"
        return self


def _realised_universe(X: np.ndarray, ids: list):
    keep = X.std(axis=0) > 0
    return X[:, keep], [g for g, k in zip(ids, keep) if k], int((~keep).sum())


# --------------------------------------------------------------------------
def load_arm_a(*, expected_n: int = 99) -> Arm:
    """GSE267145 -- bulk RNA, NASH-CRN components, fibrosis 0-3 (no F4)."""
    ep_path = BENCH / "executions/model-data-061-21079623/activation/participant_endpoints.tsv"
    mol = BENCH / "executions/model-data-064-21079902/fixture/molecular"
    ep = read_tsv_strict(ep_path, label="arm_a endpoints", dtype=str)
    require_join(len(ep), expected_n, "arm_a endpoints")
    axis = read_tsv_strict(mol / "participant_axis.tsv", label="arm_a participant axis",
                           dtype=str)
    feat = read_tsv_strict(mol / "rna_feature_axis.tsv", label="arm_a feature axis",
                           dtype=str)
    X = np.load(mol / "rna_values.npy")
    require_join(X.shape[0], len(axis), "arm_a molecular participant axis")
    require_join(X.shape[1], len(feat), "arm_a molecular feature axis")

    order = ep.set_index("participant_id")
    keep = [p for p in axis.participant_id if p in order.index]
    require_join(len(keep), expected_n, "arm_a endpoint-to-matrix join")
    idx = {p: i for i, p in enumerate(axis.participant_id)}
    rows = [idx[p] for p in keep]
    ep = order.loc[keep]

    aspects = pd.DataFrame({
        "metabolism_steatosis": ep.steatosis.astype(float).to_numpy(),
        "inflammation_lobular": ep.lobular_inflammation.astype(float).to_numpy(),
        "ballooning": ep.ballooning.astype(float).to_numpy(),
        "fibrosis": ep.fibrosis.astype(float).to_numpy(),
        "composite_nas": ep.nash_crn_component_sum.astype(float).to_numpy(),
    }, index=keep)
    nuisance = pd.DataFrame(
        {"sex_is_male": (ep.recorded_sex == "M").astype(float).to_numpy()},
        index=keep)
    Xa, ids, dropped = _realised_universe(X[rows, :], list(feat.stable_gene_id))
    return Arm(
        arm_id="A", cohort="GSE267145", assay="bulk RNA-seq (liver)",
        participants=keep, aspects=aspects, nuisance=nuisance,
        features=Xa, feature_ids=ids, composite="composite_nas",
        components=["metabolism_steatosis", "ballooning", "inflammation_lobular"],
        inputs={str(p.relative_to(BENCH)): sha256_file(p) for p in
                [ep_path, mol / "participant_axis.tsv", mol / "rna_feature_axis.tsv",
                 mol / "rna_values.npy"]},
        notes=[f"realised feature universe {len(ids)} of {X.shape[1]} "
               f"({dropped} constant across all {expected_n} participants)",
               "fibrosis is graded 0-3 in this cohort; no F4 exists"]).check()


def load_arm_b(*, expected_n: int = 180) -> Arm:
    """GSE135251 -- bulk RNA. Deposits the NAS SUM only: the three components
    are not separate columns, so metabolism and inflammation cannot be asked
    here at all. That is a property of the substrate, not a choice."""
    src = BENCH / "executions/model-data-880-21130257-gse135251-source/source"
    ep = read_tsv_strict(src / "outcomes/participant_endpoints.tsv",
                         label="arm_b endpoints", dtype=str)
    require_join(len(ep), expected_n, "arm_b endpoints")
    axis = read_tsv_strict(src / "molecular/participant_axis.tsv",
                           label="arm_b participant axis", dtype=str)
    feat = read_tsv_strict(src / "molecular/rna_feature_axis.tsv",
                           label="arm_b feature axis", dtype=str)
    X = np.load(src / "molecular/rna_values.npy")
    require_join(X.shape[0], len(axis), "arm_b molecular participant axis")
    require_join(X.shape[1], len(feat), "arm_b molecular feature axis")
    order = ep.set_index("participant_id")
    keep = [p for p in axis.participant_id if p in order.index]
    require_join(len(keep), expected_n, "arm_b endpoint-to-matrix join")
    idx = {p: i for i, p in enumerate(axis.participant_id)}
    rows = [idx[p] for p in keep]
    ep = order.loc[keep]
    aspects = pd.DataFrame({
        "fibrosis": ep.fibrosis_stage.astype(float).to_numpy(),
        "composite_nas": ep.nas_score.astype(float).to_numpy(),
    }, index=keep)
    nuisance = pd.DataFrame(index=keep)
    Xa, ids, dropped = _realised_universe(X[rows, :], list(feat.stable_gene_id))
    return Arm(
        arm_id="B", cohort="GSE135251", assay="bulk RNA-seq (liver)",
        participants=keep, aspects=aspects, nuisance=nuisance,
        features=Xa, feature_ids=ids, composite="composite_nas", components=[],
        inputs={str(p.relative_to(BENCH)): sha256_file(p) for p in
                [src / "outcomes/participant_endpoints.tsv",
                 src / "molecular/participant_axis.tsv",
                 src / "molecular/rna_feature_axis.tsv",
                 src / "molecular/rna_values.npy"]},
        notes=[f"realised feature universe {len(ids)} of {X.shape[1]} "
               f"({dropped} constant across all {expected_n} participants)",
               "no sex column is deposited, so the sex positive control cannot "
               "be run in this arm",
               "NAS components are not deposited: metabolism and inflammation "
               "are untestable here for want of a label, not for want of power"]
        ).check()


def load_arm_c(*, expected_n: int = 58, min_observed: int = 58) -> Arm:
    """PXD051911 -- liver protein. Full NASH-CRN components, Kleiner F0-F3.
    Complete-case proteins only for the primary, matching the 4019-protein
    family the 2026-08-27 protein run used."""
    src = BENCH / "executions/pxd051911-activation-readiness-21109461/source"
    meta = read_tsv_strict(src / "meta_data.txt", label="arm_c metadata", dtype=str)
    liv = meta[meta.liver_proteomics_filename.notna()
               & meta.liver_proteomics_filename.ne("NA")].copy()
    require_join(len(liv), expected_n, "arm_c liver arm")
    require_join(liv.patient_name.nunique(), expected_n, "arm_c patient uniqueness")
    q = read_tsv_strict(src / "liver_protein_quant.txt", label="arm_c quant",
                        low_memory=False)
    cols = [c for c in q.columns if c in set(liv.liver_proteomics_filename)]
    require_join(len(cols), expected_n, "arm_c matrix-to-metadata join")
    liv = liv.set_index("liver_proteomics_filename").loc[cols].reset_index()
    X = q[cols].to_numpy(float).T
    obs = np.isfinite(X).sum(axis=0)
    keepf = obs >= min_observed
    Xc = X[:, keepf]
    ids = [f"{g}|{a}" for g, a in zip(q.Genes.astype(str)[keepf],
                                      q.ProteinAccessions.astype(str)[keepf])]
    Xc, ids, dropped = _realised_universe(Xc, ids)
    F = liv.kleiner_fibrosis_grade.str.extract(r"F(\d)")[0].astype(float)
    aspects = pd.DataFrame({
        "metabolism_steatosis": liv.steatosis_score.astype(float).to_numpy(),
        "inflammation_lobular": liv.lobular_inflammation_score.astype(float).to_numpy(),
        "ballooning": liv.hepatocellular_ballooning_score.astype(float).to_numpy(),
        "fibrosis": F.to_numpy(),
        "composite_nas": liv.nafld_activity_score.astype(float).to_numpy(),
    }, index=list(liv.patient_name))
    nuisance = pd.DataFrame({
        "batch_atlasliver": liv.liver_proteomics_filename.str.contains(
            "ATLASLiver").astype(float).to_numpy(),
        "sex_is_male": (liv.gender == "Male").astype(float).to_numpy(),
    }, index=list(liv.patient_name))
    return Arm(
        arm_id="C", cohort="PXD051911", assay="liver protein (DIA-MS)",
        participants=list(liv.patient_name), aspects=aspects, nuisance=nuisance,
        features=Xc, feature_ids=ids, composite="composite_nas",
        components=["metabolism_steatosis", "ballooning", "inflammation_lobular"],
        inputs={str(p.relative_to(BENCH)): sha256_file(p) for p in
                [src / "meta_data.txt", src / "liver_protein_quant.txt"]},
        notes=[f"complete-case proteins ({min_observed}/{expected_n} observed): "
               f"{len(ids)} of {X.shape[1]} ({dropped} constant)",
               "Kleiner fibrosis is graded F0-F3 in this arm; no F4 exists",
               "two acquisition batches; batch is a covariate in every direction"]
        ).check()


def load_arm_d(*, lineage: str = "Hepatocytes", min_cells: int = 20,
               min_donors: int = 25, expressed_fraction: float = 0.5) -> Arm:
    """GSE202379 -- single-cell, donor-level pseudobulk in one lineage.

    SAF grades a Steatosis, an ACTIVITY and a Fibrosis axis. SAF activity
    conflates ballooning and lobular inflammation into one 0-4 grade, so
    inflammation is not separable in this arm at all.
    """
    saf = read_tsv_strict(SAF_LOCK / "donor_saf_grades.tsv",
                          label="arm_d saf grades", dtype=str)
    graded = saf[saf.saf_recovered == "True"].copy()
    require_join(len(graded), 40, "arm_d SAF-graded donors")
    pair = read_tsv_strict(DONOR_PAIRING, sep=",", label="arm_d donor pairing",
                           dtype=str)
    srr = {r.donor_id: [s for s in str(r.rna_srrs).split(";") if s]
           for _, r in pair.iterrows()}
    counts_path = PSEUDOBULK / f"{lineage.lower().replace(' ', '_')}_counts.tsv.gz"
    cc = read_tsv_strict(PSEUDOBULK / "cell_counts_per_donor_ct.tsv",
                         label="arm_d cell counts")
    cc = cc[cc.cell_type == lineage]
    per_sample = dict(zip(cc["sample"], cc.n_cells))
    mat = read_tsv_strict(counts_path, label=f"arm_d {lineage} counts")
    gene_col = mat.columns[0]
    have = set(mat.columns[1:])

    donors, cells, blocks = [], [], []
    for _, row in graded.iterrows():
        runs = [s for s in srr.get(row.donor_id, []) if s in have]
        n_cells = sum(float(per_sample.get(s, 0)) for s in runs)
        if not runs or n_cells < min_cells:
            continue
        donors.append(row.donor_id)
        cells.append(n_cells)
        blocks.append(mat[runs].to_numpy(float).sum(axis=1))
    if len(donors) < min_donors:
        raise JoinError(f"arm_d {lineage}: only {len(donors)} donors clear "
                        f"{min_cells} cells, floor is {min_donors}")
    C = np.vstack(blocks)                                   # donors x genes
    lib = C.sum(axis=1, keepdims=True)
    if not (lib > 0).all():
        raise JoinError(f"arm_d {lineage}: a donor pseudobulk is empty")
    L = np.log1p(C / lib * 1e6)
    detected = (C > 0).mean(axis=0) >= expressed_fraction
    L = L[:, detected]
    ids = [f"{lineage}|{g}" for g, d in zip(mat[gene_col].astype(str), detected) if d]
    L, ids, dropped = _realised_universe(L, ids)
    g = graded.set_index("donor_id").loc[donors]
    aspects = pd.DataFrame({
        "metabolism_steatosis": g.steatosis_S.astype(float).to_numpy(),
        "fibrosis": g.fibrosis_F.astype(float).to_numpy(),
        "composite_nas": g.activity_A.astype(float).to_numpy(),
    }, index=donors)
    nuisance = pd.DataFrame({
        "sex_is_male": (g.gender == "M").astype(float).to_numpy(),
        "log10_cells": np.log10(np.asarray(cells, float)),
    }, index=donors)
    return Arm(
        arm_id="D", cohort="GSE202379", assay=f"single-cell pseudobulk ({lineage})",
        participants=donors, aspects=aspects, nuisance=nuisance,
        features=L, feature_ids=ids, composite="composite_nas", components=[],
        inputs={str(p): sha256_file(p) for p in
                [SAF_LOCK / "donor_saf_grades.tsv", DONOR_PAIRING, counts_path,
                 PSEUDOBULK / "cell_counts_per_donor_ct.tsv"]},
        notes=[f"{len(donors)} of 40 SAF-graded donors reach {min_cells} "
               f"{lineage} cells and join the pseudobulk",
               f"expressed-gene universe {len(ids)} ({dropped} constant)",
               "SAF activity conflates ballooning and lobular inflammation; "
               "inflammation is untestable in this arm for want of a label",
               "fibrosis is the SAF-recovered F, never the label-map F_343"]
        ).check()
