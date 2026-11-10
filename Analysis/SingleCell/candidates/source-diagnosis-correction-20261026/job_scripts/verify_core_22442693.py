"""Independent checks of source identities, donor inference and testing families."""
import argparse
import json
import math
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy.stats import t as t_dist

ROOT = Path(__file__).resolve().parents[4]
TASK = ROOT / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"


def read(path):
    return pd.read_csv(path, sep="\t")


def bh(values, n=None):
    p = np.asarray(values, dtype=float)
    good = np.flatnonzero(np.isfinite(p))
    order = good[np.argsort(p[good])]
    ans = np.full(len(p), np.nan)
    ans[order] = np.minimum.accumulate((p[order] * (n or len(p)) / np.arange(1, len(order) + 1))[::-1])[::-1].clip(0, 1)
    return ans


def same(a, b):
    assert np.allclose(a, b, atol=1e-10, rtol=1e-8, equal_nan=True)


def core():
    result = {}
    lib = read(TASK / "metadata/library_metadata_source_F.tsv")
    donors = read(TASK / "metadata/donor_metadata.tsv")
    assert len(lib) == 275 and lib["sample"].nunique() == 275 and len(donors) == 104
    assert lib.n_cells.sum() == 1232318 and donors.analyzed_donor.sum() == 102
    counts = lib.groupby("donor").n_cells.sum()
    same(counts.loc[donors.donor], donors.n_cells)
    assert not lib.loc[lib.dataset.eq("Liver_Atlas"), "binary_eligible"].any()
    assert not lib.loc[lib.dataset.eq("GSE189600"), "binary_eligible"].any()
    assert not lib.loc[lib.source_etiology.isin(["Alcohol", "PBC"]), "binary_eligible"].any()
    q = donors.set_index("donor")
    assert q.loc["GSE202379_P98", "disease_stage_coarse"] == "Unresolved"
    assert q.loc["GSE202379_P98", "exclude_stage_analysis"]
    assert q.loc["GSE202379_P30", "binary_eligible"] and not q.loc["GSE202379_P30", "strict_binary_eligible"]
    assert q.loc["GSE244832_D11", "disease_stage_coarse"] == "Steatosis"
    assert q.loc["GSE244832_D14", "disease_stage_coarse"] == "Steatohepatitis"
    f = read(TASK / "metadata/documented_F_source.tsv")
    assert f.F_stage_documented.notna().sum() == 39
    assert f.F_stage_documented.isin([0, 1]).sum() == 12
    assert f.set_index("donor").loc["GSE202379_P98", "F_stage_documented"] == 1
    assert f.loc[f.deposited_SAF.isin(["healthy control", "end stage"]), "F_stage_documented"].isna().all()
    result["metadata"] = dict(libraries=275, source_identifiers=104, analyzed_identifiers=102, cells=1232318, measured_SAF_donors=39, F0_F1_donors=12)
    with h5py.File(TASK / "integrated_atlas_source_qualified_v2.h5ad", "r") as h:
        obs = read_elem(h["obs"])
        assert len(obs) == 1232318 and "preparation_method" in obs
        matched = lib.set_index("sample").loc[obs["sample"].astype(str)]
        assert np.array_equal(obs.donor.astype(str), matched.donor.astype(str))
        assert np.array_equal(obs.condition_harmonized.astype(str), matched.condition_harmonized.astype(str))
        assert np.array_equal(obs.preparation_method.astype(str), matched.preparation_method.astype(str))
        link = h.get("raw", getlink=True)
        assert isinstance(link, h5py.ExternalLink)
        with h5py.File(link.filename, "r") as original:
            assert h["raw/X/data"].id == original["raw/X/data"].id
            assert h["X"].id == original["X"].id
    result["atlas_obs_and_expression_links"] = "passed"
    tests = []
    for folder in ["pseudobulk_de_v2", "pseudobulk_de_strict_control", "pseudobulk_de_coarse"]:
        files = sorted((TASK / folder).glob("*_de.csv"))
        assert len(files) == (27 if folder.endswith("coarse") else 9)
        for p in files:
            d = pd.read_csv(p)
            assert not d.gene.duplicated().any()
            same(bh(d.pvalue), d.padj)
            tests.append(dict(analysis=folder, table=p.name, tested=len(d), BH05=int(d.padj.lt(.05).sum())))
        if not folder.endswith("coarse"):
            for p in (TASK / folder).glob("*_donor_roster.tsv"):
                roster = read(p)
                assert not roster["sample"].duplicated().any()
                allowed = donors[donors.strict_binary_eligible if "strict" in folder else donors.binary_eligible]
                assert set(roster["sample"]).issubset(allowed.donor)
    result["gene_tests"] = tests
    registry = read(TASK / "OLS_ambient/hotspot/program_registry_v2.tsv")
    assert len(registry) == 117
    old = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot"
    assert (old / "program_membership_v2.tsv").read_bytes() == (TASK / "OLS_ambient/hotspot/program_membership_v2.tsv").read_bytes()
    scores = read(TASK / "OLS_ambient/hotspot/donor_program_scores_primary.tsv")
    effects = read(TASK / "OLS_ambient/hotspot/primary_effects.tsv")
    for _, r in effects[effects.estimable].iterrows():
        d = scores[scores.cell_type.eq(r.cell_type) & scores.module.eq(r.module)]
        X = np.column_stack([np.ones(len(d)), d.stage_ordinal, pd.get_dummies(d.dataset, drop_first=True, dtype=float)])
        y = d.score.to_numpy(); bread = np.linalg.inv(X.T @ X); beta = bread @ X.T @ y
        e = y - X @ beta; df = len(d) - np.linalg.matrix_rank(X)
        se = math.sqrt(e @ e / df * bread[1, 1])
        leverage = np.einsum("ij,jk,ik->i", X, bread, X)
        meat = X.T @ ((e / (1 - leverage))[:, None] ** 2 * X)
        hc3 = bread @ meat @ bread
        pv = 2 * t_dist.sf(abs(beta[1] / se), df)
        hp = 2 * t_dist.sf(abs(beta[1] / np.sqrt(hc3[1, 1])), df)
        same([beta[1], se, pv, hp], [r.beta, r.se, r.pvalue, r.hc3_pvalue])
    result["fixed_program_weights_and_independent_OLS_HC3"] = "passed; 61-record sensitivity includes two GSE189600 records of unresolved independence"
    ccc = read(TASK / "CCC/fig4d_program_ccc_associations.tsv")
    assert len(ccc) == 26
    same(bh(ccc.pvalue, 26), ccc.qvalue)
    same(bh(ccc.hc3_pvalue, 26), ccc.hc3_qvalue)
    membership = read(TASK / "CCC/fig4d_program_ccc_matched_donors.tsv")
    assert not membership.donor.eq("GSE202379_P98").any()
    labels = membership.merge(donors[["donor", "disease_stage_coarse"]], on="donor", suffixes=("_fit", "_source"), validate="many_to_one")
    assert labels.disease_stage_coarse_fit.eq(labels.disease_stage_coarse_source).all()
    full = read(TASK / "CCC/full_donor_lr_family_dispositions.tsv")
    assert len(full) == 55402 and full.disposition.eq("fitted").sum() == 374
    result["communication"] = dict(planned=26, estimable=int(ccc.testable.sum()), supported=int(ccc.evidence_state.eq("supported").sum()), full_groups=55402, fitted_groups=374)
    for name in ["lineages", "hepatocyte_subtype", "meta_subtype"]:
        d = read(TASK / f"composition/{name}_results.tsv")
        same(bh(d["P.Value"]), d["adj.P.Val"])
        roster = read(TASK / f"composition/{name}_donor_counts.tsv")
        assert not roster.donor.duplicated().any()
        assert not roster.dataset.isin(["Liver_Atlas", "GSE136103", "GSE189600"]).any()
    comp = pd.read_csv(TASK / "allcell_v2/allcell_vs_bulk_comparison.csv")
    assert len(comp) > 10000 and not comp.gene_key.duplicated().any()
    result["allcell_vs_bulk"] = dict(genes=len(comp), allcell_rho=comp.bulk_lfc.corr(comp.allcell_lfc, method="spearman"), hepatocyte_rho=comp.bulk_lfc.corr(comp.hep_lfc, method="spearman"))
    web = json.loads((TASK / "portal/singlecell_summary.json").read_text())
    assert web["n_cells"] == 1232318 and len(web["composition"]) == 275
    assert web["celltype_disease_activity"] == [] and web["fgsea_transitions"] == []
    assert not json.loads((TASK / "portal/source.json").read_text())["published"]
    assert len(pd.read_parquet(TASK / "portal/gene_pseudobulk_de.parquet")) > 0
    result["local_exports"] = "source-qualified gene and donor-composition data; obsolete label-derived scores withdrawn; not published"
    return result


def models():
    result = {}
    for mode, count in [("programs", 117), ("programs_qc", 117), ("programs_ambient", 106), ("programs_tcell20", 11)]:
        d = read(TASK / mode / "results_primary_family.tsv")
        assert len(d) == count and d.permutation_draws.eq(19999).all()
        assert not d.program.duplicated().any()
        if mode != "programs_tcell20":
            same(bh(d.perm_p, 117), d.perm_q)
        else:
            full = read(TASK / mode / "full_117_primary_p_with_tcell20.tsv")
            same(bh(full.perm_p, 117), full.perm_q)
        result[mode] = dict(programs=count, selected=int(d.selected_H59.sum()), smallest_q=float(d.perm_q.min()))
    assert (TASK / "programs_diagnostics/COMPLETED").exists()
    for tag in ["raw", "corrected"]:
        base = TASK / "N3" / tag / "results_source_F"
        summary = json.loads((base / "summary.json").read_text())
        d = read(base / "target_classes.tsv")
        tested = d[d.testable]
        assert len(d) == 51 and summary["n_donors"] == 69
        assert summary["permutation"]["B"] == 1000 and summary["permutation"]["seed"] == 20260923
        same(bh(tested["P.Value"]), tested.q_bh)
        empirical = read(base / "permutation_target_p.tsv").set_index("symbol")
        same(bh(empirical.loc[tested.symbol, "emp_p"]), tested.q_emp)
        expected = tested.q_bh if summary["permutation"]["calibration_pass"] else tested.q_emp
        same(expected, tested.q_primary)
        f = read(base / "numeric_SAF_donors.tsv")
        assert f.F_documented.notna().all() and not f.donor.isin(["GSE202379_PHL1", "GSE202379_PHL2", "GSE202379_P98"]).any()
        assert f.donor.isin(read(TASK / "metadata/documented_F_source.tsv").dropna(subset=["F_stage_documented"]).donor).all()
        result["N3_" + tag] = dict(donors=69, testable=len(tested), primary_hits=int(tested.q_primary.lt(.05).sum()), numeric_SAF_donors=len(f), calibration_pass=summary["permutation"]["calibration_pass"])
    for marker in ["H17_COMPLETED", "H34_COMPLETED", "BF3_COMPLETED", "FIGURES_COMPLETED", "N3_DOCUMENTED_F_COMPLETED"]:
        assert (TASK / marker).exists(), marker
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["core", "models"], required=True)
    args = parser.parse_args()
    data = core() if args.phase == "core" else models()
    out = TASK / (args.phase + "_verification.json")
    assert not out.exists()
    out.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    print(args.phase, "verification passed")
