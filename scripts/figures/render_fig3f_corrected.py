#!/usr/bin/env python3
"""Render the existing adjacent-stage fits against corrected exact COLOC pairs.

KEY MESSAGE: Cross-sectional expression effects and inherited regulatory support
are different measurements; an unavailable GWAS PIP is not negative evidence.
No models are refitted and no memberships or weights are changed.
"""
import argparse
import hashlib
import json
import os
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bh(p):
    order = np.argsort(p, kind="stable")
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(q, 1)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    src = out / "source_tables"
    src.mkdir()
    genetics = root / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
    # FIG3F_STAGE_ROOT and the three label variables below select another stage
    # release; unset, the 2026-09-25 tables and their labels are used.
    stage_root = Path(os.environ.get("FIG3F_STAGE_ROOT", root / "figures/candidates/bulk-release-corrected-controls-20260925T175215Z/analysis/stage_extensions"))
    paths = {
        "stage": stage_root / "stage_all_gene_results.tsv",
        "roster": stage_root / "stage_extension_sample_manifest.tsv",
        "genes": genetics / "gene_level_coloc_tier12.csv",
        "pairs": genetics / "susie_coloc_signal_pairs.tsv.gz",
        "studies": genetics / "susie_coloc_all_gwas.csv",
        "named": genetics / "summary/named_gene_status.tsv",
        "identity": root / "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/identity/sample_to_individual.tsv",
    }
    assert (genetics / "ASSEMBLY_COMPLETE.json").is_file()
    recorded = dict(line.split(None, 1)[::-1] for line in (genetics / "output_sha256.tsv").read_text().splitlines())
    for name in ("genes", "pairs", "studies"):
        matches = [v for k, v in recorded.items() if Path(k).name == paths[name].name]
        assert matches == [sha(paths[name])], f"Genetics checksum: {name}"
    stage = pd.read_csv(paths["stage"], sep="\t")
    transitions = ["F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"]
    stage = stage[stage.contrast.isin(transitions)].copy()
    assert not stage.duplicated(["gene_id_base", "contrast"]).any()
    assert stage.groupby("contrast").gene_id_base.apply(set).apply(len).nunique() == 1
    summaries = []
    for contrast, rows in stage.groupby("contrast", sort=False):
        assert rows.bh_family_size.nunique() == 1 and len(rows) == rows.bh_family_size.iloc[0]
        assert np.max(np.abs(bh(rows["P.Value"].to_numpy()) - rows.FDR.to_numpy())) < 1e-10
        summaries.append(rows.iloc[0][["contrast", "n_reference", "n_comparison", "n_samples", "n_cohorts", "covariates", "bh_family_size"]].to_dict())
    assert all(set(g.gene_id_base) == set(stage[stage.contrast == transitions[0]].gene_id_base)
               for _, g in stage.groupby("contrast"))
    genes = pd.read_csv(paths["genes"])
    assert genes.ensembl.is_unique
    positive = genes[genes.coloc_best_susie_pp4 > 0.5].copy()
    assert len(positive) == 436
    pairs = pd.read_csv(paths["pairs"], sep="\t")
    assert not pairs.duplicated(["ensembl", "gwas", "idx1", "idx2"]).any()
    # The tier-1/2 strongest driver fixes the study; choose its eligible exact pair.
    chosen = pairs[pairs.eligible.eq(True)].merge(
        positive[["ensembl", "susie_driving_gwas", "coloc_best_susie_pp4"]],
        left_on=["ensembl", "gwas"], right_on=["ensembl", "susie_driving_gwas"],
        validate="many_to_one")
    chosen = chosen.sort_values(["ensembl", "PP.H4", "hit1", "hit2", "idx1", "idx2"],
                                ascending=[True, False, True, True, True, True]).drop_duplicates("ensembl")
    assert len(chosen) == 436
    assert np.max(np.abs(chosen["PP.H4"] - chosen.coloc_best_susie_pp4)) < 1e-8
    assert chosen["PP.H4"].between(0.5, 1, inclusive="right").all()
    # top_snp_PP in the gene-study summary is a coloc SNP posterior, NOT a GWAS PIP.
    # The saved exact-pair exports do not contain the matching GWAS alpha/PIP.
    chosen["pip"] = np.nan
    chosen["pip_evidence_state"] = "untestable_no_saved_exact_signal_gwas_pip"
    chosen["coloc_pair_id"] = chosen.apply(
        lambda r: f"{r.ensembl}|{r.gwas}|GWAS{r.idx1}:{r.hit1}|eQTL{r.idx2}:{r.hit2}", axis=1)
    eligible = set(stage.loc[stage.FDR < 0.05, "gene_id_base"])
    chosen = chosen[chosen.ensembl.isin(eligible)].copy()
    assert len(chosen) > 0
    matrix = stage.merge(chosen, left_on="gene_id_base", right_on="ensembl", validate="many_to_one")
    matrix["evidence_state"] = os.environ.get("FIG3F_EVIDENCE_STATE", "historical_library_roster_conditional_not_donor_qualified")
    matrix["estimand"] = "cross_sectional_higher_minus_lower_adjacent_fibrosis_stage"
    matrix["units"] = "log2_expression_fold_change"
    matrix["uncertainty"] = "pointwise_95_percent_moderated_t_interval"
    matrix["testing_family"] = "BH_all_tested_genes_separately_per_adjacent_contrast"
    # Original row rule: strongest absolute transition, direction, minimum FDR.
    effects = matrix.pivot(index="ensembl", columns="contrast", values="logFC")[transitions]
    fdr = matrix.pivot(index="ensembl", columns="contrast", values="FDR")[transitions]
    strongest = np.abs(effects.to_numpy()).argmax(axis=1)
    ordering = pd.DataFrame({"ensembl": effects.index, "transition_rank": strongest,
                            "direction": np.sign(effects.to_numpy()[np.arange(len(effects)), strongest]),
                            "minimum_FDR": fdr.min(axis=1).to_numpy()})
    ordering = ordering.merge(chosen[["ensembl", "gene"]], on="ensembl", validate="one_to_one")
    ordering = ordering.sort_values(["transition_rank", "direction", "minimum_FDR", "gene"],
                                     ascending=[True, False, True, True]).reset_index(drop=True)
    ordering["row_order"] = np.arange(1, len(ordering)+1)
    # Label every eligible gene, avoiding selection of new biological examples.
    ordering["row_label"] = ordering.gene
    matrix = matrix.merge(ordering[["ensembl", "row_order", "row_label"]], on="ensembl", validate="many_to_one")
    assert len(matrix) == 4*len(chosen)
    matrix = matrix.sort_values(["row_order", "contrast"])
    lane = chosen.merge(ordering[["ensembl", "row_order", "row_label"]], on="ensembl", validate="one_to_one").sort_values("row_order")
    matrix.to_csv(src / "fig3f_stage_genetics_matrix.tsv", sep="\t", index=False)
    lane.to_csv(src / "fig3f_exact_pair_lanes.tsv", sep="\t", index=False)
    ordering.to_csv(src / "fig3f_row_order.tsv", sep="\t", index=False)
    # n is documented from the same fitted roster, without relabeling libraries as people.
    roster = pd.read_csv(paths["roster"], sep="\t")
    identity = pd.read_csv(paths["identity"], sep="\t", dtype={"run": str, "individual_id": str})
    assert identity.run.is_unique
    for item in summaries:
        r = roster[(roster.axis == "fibrosis_adjacent") & (roster.contrast == item["contrast"])].copy()
        assert r.sample_id.is_unique and len(r) == item["n_samples"]
        r = r.merge(identity[["run", "individual_id"]], left_on="sample_id", right_on="run", how="left", validate="one_to_one")
        item["n_genotype_identity_components"] = r.individual_id.nunique()
        item["n_identity_unresolved_records"] = int(r.individual_id.isna().sum())
        item["n_components_on_both_sides"] = int((r.groupby("individual_id").group.nunique() > 1).sum())
        item["inferential_unit_in_existing_fit"] = os.environ.get("FIG3F_INFERENTIAL_UNIT", "library_record_not_verified_independent_donor")
        item["cohorts"] = ";".join(sorted(r.dataset.unique()))
    pd.DataFrame(summaries).to_csv(src / "fig3f_contrast_denominators.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.family": "Helvetica", "font.size": 6, "pdf.use14corefonts": True,
                         "axes.labelsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
                         "legend.fontsize": 6, "axes.unicode_minus": False,
                         "font.weight": "normal", "axes.labelweight": "normal"})
    n = len(ordering)
    height = max(2.6, n*0.105+0.75)
    fig = plt.figure(figsize=(3.6, height))
    grid = fig.add_gridspec(2, 3, height_ratios=[1, .035], width_ratios=[4, .55, .55],
                            left=.22, right=.97, top=.97, bottom=.12, wspace=.22, hspace=.12)
    axes = [fig.add_subplot(grid[0, k]) for k in range(3)]
    ordered = ordering.ensembl.tolist()
    values = effects.loc[ordered].to_numpy()
    qvalues = fdr.loc[ordered].to_numpy()
    lim = float(np.quantile(np.abs(values), .98))
    assert lim > 0
    cmap = LinearSegmentedColormap.from_list("stage", ["#1565C0", "white", "#C9265E"])
    im = axes[0].imshow(values, cmap=cmap, vmin=-lim, vmax=lim, aspect="auto", interpolation="none")
    yy, xx = np.where(qvalues < .05)
    axes[0].scatter(xx, yy, s=1.6, c="black", linewidths=0)
    axes[0].set_xticks(range(4), ["F1-F0", "F2-F1", "F3-F2", "F4-F3"])
    labelled = np.flatnonzero(ordering.row_label.ne("").to_numpy())
    axes[0].set_yticks(labelled, ordering.row_label.iloc[labelled])
    axes[0].set_ylabel(f"{n} genes")
    pp = axes[1].imshow(lane[["PP.H4"]].to_numpy(), cmap=LinearSegmentedColormap.from_list(
        "pph4", ["white", "#ED5565"]), vmin=0, vmax=1, aspect="auto", interpolation="none")
    axes[1].set_xticks([0], ["PP.H4"])
    axes[1].set_yticks([])
    axes[2].imshow(np.ones((n, 1)), cmap=LinearSegmentedColormap.from_list("missing", ["#eeeeee", "#eeeeee"]),
                   vmin=0, vmax=1, aspect="auto", interpolation="none")
    axes[2].scatter(np.zeros(n), np.arange(n), marker="x", s=2.5, linewidths=.35, c="#656D77")
    axes[2].set_xticks([0], ["PIP"])
    axes[2].set_yticks([])
    for ax in axes:
        ax.tick_params(length=0, pad=3)
        for spine in ax.spines.values():
            spine.set_visible(False)
    cb = fig.colorbar(im, cax=fig.add_subplot(grid[1, 0]), orientation="horizontal", ticks=[-lim, 0, lim])
    cb.set_label("log2FC")
    cb.set_ticklabels([f"{-lim:.1f}", "0", f"{lim:.1f}"])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=2, width=.4)
    cb2 = fig.colorbar(pp, cax=fig.add_subplot(grid[1, 1]), orientation="horizontal", ticks=[0, 1])
    cb2.outline.set_visible(False)
    cb2.ax.tick_params(length=2, width=.4)
    fig.savefig(out / "fig3f_stage_deg_genetics_matrix.pdf", bbox_inches="tight", pad_inches=.04)
    plt.close(fig)
    summary = {"genes": n, "rows": len(matrix), "BH_family_size": sorted(stage.bh_family_size.unique().tolist()),
               "exact_PIP_available": 0, "labels": ordering.loc[ordering.row_label.ne(""), "row_label"].tolist(),
               "effect_color_limit": lim, "genetics_positive_family": 436,
               "evidence_state": os.environ.get("FIG3F_SUMMARY_EVIDENCE_STATE", "historical_stage_fits_conditional;September_genetics_candidate"),
               "PIP_reason": "Saved exact-pair table lacks GWAS PIP; top_snp_PP is a different posterior.",
               "python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
               "matplotlib": matplotlib.__version__, "seed": "No stochastic operation or new model fit"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    entries = [{"role": name, "path": str(path), "sha256": sha(path)} for name, path in paths.items()]
    entries.append({"role": "renderer", "path": str(Path(__file__).resolve()), "sha256": sha(__file__)})
    for path in sorted(out.rglob("*")):
        if path.is_file():
            entries.append({"role": "output", "path": str(path), "sha256": sha(path)})
    pd.DataFrame(entries).to_csv(out / "manifest.tsv", sep="\t", index=False)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
