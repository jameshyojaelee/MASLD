#!/usr/bin/env python3
"""KEY MESSAGE: RNA recovers measured chromatin variation and allele models guide assays within their measured scope.

Render six separate vector PDF panels. Values are candidate development results;
this script never edits the adopted Figure 7 or a named release.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT/"Analysis/MASLD_Model_Benchmark"
FIX = BENCH/"executions/model-data-064-21079902/fixture"
OOF = BENCH/"executions/chromatin-stable-rrr-forms-20260908T192024Z/out/stable_oof.npz"
JOIN = BENCH/"executions/gse267145-authoritative-join-21064930/participant_join.tsv"
CONFIRMATION = ROOT/"GWAS/finemapping/results/alphagenome_campaign/fold0-confirmation-20260922T191005Z/results/stage3_contrasts_21846521/summary.json"
PURPLE = "#7B1FA2"
BLUE = "#1565C0"
MAGENTA = "#C9265E"
TEAL = "#00695C"
GRAY = "#9E9E9E"
PEACH = "#F4A674"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6,
                     "axes.titlesize": 6, "axes.labelsize": 6,
                     "xtick.labelsize": 6, "ytick.labelsize": 6,
                     "legend.fontsize": 6, "pdf.fonttype": 42,
                     "pdf.use14corefonts": False, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": .4,
                     "xtick.major.width": .4, "ytick.major.width": .4,
                     "savefig.transparent": False})


def save(fig, path):
    fig.savefig(path, format="pdf", bbox_inches="tight", pad_inches=.045)
    plt.close(fig)


def panel_a(out, panel_file):
    panel = pd.read_csv(panel_file, sep="\t").head(24)
    region = panel.region_index.to_numpy(int)
    data = np.load(OOF, allow_pickle=False)
    axis = pd.read_csv(FIX/"molecular/participant_axis.tsv", sep="\t")
    join = pd.read_csv(JOIN, sep="\t").set_index("participant_id").loc[axis.participant_id]
    order = np.lexsort((join.steatosis.to_numpy(), join.fibrosis.to_numpy()))
    fibrosis = join.fibrosis.to_numpy()[order]
    measured = data["Cres_oof"][:, region].astype(float)[order]
    predicted = data["pred_rrr_offset_cis_wide"][:, region].astype(float)[order]
    # A common region-level scale is a display transform only.
    mean = measured.mean(0)
    sd = measured.std(0)
    sd[sd < 1e-10] = 1
    measured = (measured-mean)/sd
    predicted = (predicted-mean)/sd
    pd.DataFrame({"display_participant": np.repeat(np.arange(1, 100), len(region)),
                  "participant_id": np.repeat(axis.participant_id.to_numpy()[order], len(region)),
                  "fibrosis": np.repeat(fibrosis, len(region)),
                  "steatosis": np.repeat(join.steatosis.to_numpy()[order], len(region)),
                  "region_index": np.tile(region, 99),
                  "region_key": np.tile(panel.region_key, 99),
                  "measured_display_z": measured.ravel(),
                  "predicted_display_z": predicted.ravel()}).to_csv(
                      out/"Figure7A_source.tsv", sep="\t", index=False)
    vmax = np.percentile(np.abs(measured), 98)
    fig, axes = plt.subplots(1, 2, figsize=(4.3, 2.25), sharey=True,
                             layout="constrained")
    stages = np.unique(fibrosis)
    midpoints = [np.flatnonzero(fibrosis == value).mean()+.5 for value in stages]
    boundaries = np.flatnonzero(np.diff(fibrosis))+1
    for ax, array, title in zip(axes, (measured, predicted),
                                ("Measured H3K27ac", "Held-participant RNA prediction")):
        mesh = ax.pcolormesh(np.arange(len(region)+1), np.arange(100), array,
                             cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="flat",
                             edgecolors="none", rasterized=False)
        for boundary in boundaries:
            ax.axhline(boundary, color="white", linewidth=.45)
        ax.set_title(title, fontweight="normal")
        ax.set_xticks([.5, 11.5, 23.5], ["1", "12", "24"])
        ax.set_xlabel("Annotation-fixed regions")
        ax.set_yticks(midpoints, [f"F{value:g}" for value in stages])
        ax.set_xlim(0, 24)
        ax.set_ylim(99, 0)
    axes[0].set_ylabel("99 participants; fibrosis stage")
    colorbar = fig.colorbar(mesh, ax=axes, fraction=.035, pad=.025, shrink=.8,
                           ticks=[-2, 0, 2])
    colorbar.set_label("H3K27ac display z")
    colorbar.solids.set_rasterized(False)
    save(fig, out/"Figure7A_molecular_variation.pdf")


def panel_b(out, allele):
    names = [("native1m", "Native 1 Mb", GRAY), ("adapter", "Adapter 2 kb", BLUE),
             ("combined", "Combined", MAGENTA)]
    subset = allele.assign(_hash=allele.key.map(lambda x: hashlib.sha256(str(x).encode()).hexdigest()))
    subset = subset.sort_values("_hash").head(1500)
    subset[["key", "beta_alt", "native1m", "adapter", "combined", "fold"]].to_csv(
        out/"Figure7B_source.tsv", sep="\t", index=False)
    all_values = np.concatenate([subset.beta_alt.to_numpy()]+[subset[name].to_numpy() for name,_,_ in names])
    limit = max(0.5, float(np.percentile(np.abs(all_values), 99)))
    fig, axes = plt.subplots(1, 3, figsize=(4.8, 1.75), sharex=True, sharey=True,
                             layout="constrained")
    for ax, (column, title, color) in zip(axes, names):
        ax.scatter(subset.beta_alt, subset[column], s=1.3, alpha=.25,
                   color=color, linewidths=0, rasterized=False)
        ax.plot([-limit, limit], [-limit, limit], color=GRAY, linewidth=.5)
        ax.set_title(title, fontweight="normal")
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Measured ALT dosage beta")
    axes[0].set_ylabel("Predicted ALT dosage beta")
    save(fig, out/"Figure7B_allele_effects.pdf")


def panel_c(out, chrom, sequence_path):
    skill = chrom["skill"]
    matched_path = sequence_path.parent/"within_histology_21849320"
    matched = json.loads((matched_path/"results.json").read_text())
    secondary = json.loads((sequence_path.parent/"chromatin_secondary_review/results.json").read_text())
    if matched["n_pairs"] != 22 or matched["n_participants"] != 44:
        raise ValueError("Matched histology participant census differs")
    rows = [("Global − histology/technical", "global_RRR_vs_histology_sex_RNA_technical",
             "global_RRR", "histology_sex_RNA_technical", BLUE),
            ("Global + local − histology/technical", "RRR_local_vs_histology_sex_RNA_technical",
             "RRR_local", "histology_sex_RNA_technical", MAGENTA),
            ("Global + local − global", "RRR_local_vs_global_RRR", "RRR_local", "global_RRR", MAGENTA)]
    source = [{"label": label, "model": key, "endpoint": "96,460-region overall paired skill increment",
               "population": "99 held participants", "n_participants": 99, "n_pairs": 0,
               "reference": "fold-training region mean",
               "estimate": chrom["contrasts"][key]["delta_skill"],
               "ci_low95": chrom["contrasts"][key]["ci95"][0],
               "ci_high95": chrom["contrasts"][key]["ci95"][1],
               "p_holm": chrom["contrasts"][key]["p_holm_primary_family"],
               "absolute_model_skill": skill[model]["skill"],
               "absolute_comparator_skill": skill[comparator]["skill"]}
              for label, key, model, comparator, _ in rows]
    matched_rows = [("Global + local − technical", "RRR_local_vs_technical_nuisance",
                     "RRR_local", "technical_nuisance", MAGENTA),
                    ("Global + local − global", "RRR_local_vs_global_RRR",
                     "RRR_local", "global_RRR", MAGENTA)]
    source.extend({"label": label, "model": "matched_"+key,
                   "endpoint": "96,460-region within-pair paired skill increment",
                   "population": "22 disjoint pairs matched on four histology grades, sex and held fold",
                   "n_participants": 44, "n_pairs": 22, "reference": "zero within-pair difference",
                   "estimate": matched["contrasts"][key]["delta_skill"],
                   "ci_low95": matched["contrasts"][key]["ci95"][0],
                   "ci_high95": matched["contrasts"][key]["ci95"][1],
                   "p_holm": matched["contrasts"][key]["p_holm_two_contrasts"],
                   "absolute_model_skill": matched["skill"][model]["skill"],
                   "absolute_comparator_skill": matched["skill"][comparator]["skill"]}
                  for label, key, model, comparator, _ in matched_rows)
    sequence = secondary["sequence"]
    delta = sequence["interaction_minus_additive_original_skill"]
    interval = sequence["interaction_minus_additive_original_skill_ci95"]
    source.append({"label": "Sequence × RNA residual − additive", "model": "interaction_minus_additive",
                   "endpoint": "1,024-region original-target skill increment",
                   "population": "99 held participants", "n_participants": 99, "n_pairs": 0,
                   "reference": "fold-training region mean of original concentration-residual target",
                   "p_holm": np.nan, "absolute_model_skill": sequence["original_skill"]["interaction"],
                   "absolute_comparator_skill": sequence["original_skill"]["additive"],
                   "sequence_coordinate_semantics": "counted_interval", "sequence_coordinate_residual_bp": 1,
                   "estimate": delta, "ci_low95": interval[0], "ci_high95": interval[1]})
    pd.DataFrame(source).to_csv(out/"Figure7C_source.tsv", sep="\t", index=False)
    fig, (ax, within, residual) = plt.subplots(3, 1, figsize=(4.5, 3.1),
        gridspec_kw={"height_ratios": [2, 1.5, 1]}, layout="constrained")
    for axis, entries, estimates, label in (
            (ax, rows, chrom["contrasts"], "Overall skill increment; 99 participants"),
            (within, matched_rows, matched["contrasts"], "Within-pair skill increment; 22 pairs / 44 participants")):
        for j, (_, key, _, _, color) in enumerate(entries):
            entry = estimates[key]
            low, high = entry["ci95"]
            axis.plot([low, high], [j, j], color=color, linewidth=1.5)
            axis.plot(entry["delta_skill"], j, "o", color=color, markersize=3)
        axis.axvline(0, color=GRAY, linewidth=.5)
        axis.set_yticks(range(len(entries)), [row[0] for row in entries])
        axis.set_ylim(len(entries)-.5, -.5)
        axis.set_xlim(-.02, .18)
        axis.set_xticks([0, .05, .10, .15])
        axis.set_xlabel(label)
    residual.plot([interval[0], interval[1]], [0, 0], color=PURPLE, linewidth=1.5)
    residual.plot(delta, 0, "o", color=PURPLE, markersize=3)
    residual.axvline(0, color=GRAY, linewidth=.5)
    residual.set_yticks([0], ["Sequence × RNA − additive"])
    residual.set_xlim(-.01, .01)
    residual.set_xticks([-.01, 0, .01])
    residual.set_xlabel("Original-target skill increment; 1,024 regions")
    save(fig, out/"Figure7C_information_added.pdf")


def recorded_confirmation():
    """Read the completed aggregate only; no fold-0 row data are opened."""
    raw = CONFIRMATION.read_bytes()
    summary = json.loads(raw)
    if summary["status"] != "complete" or summary["further_fold0_analysis_licensed"]:
        raise ValueError("Historical confirmation summary state differs")
    contrasts = {}
    for comparator in ("native_atac_2048", "frozen_2048_symmetric_shared_mlp64"):
        rows = [entry for entry in summary["contrasts"]
                if entry["arm"] == "adapter_r16_last5_symmetric" and entry["comparator"] == comparator]
        if len(rows) != 1:
            raise ValueError("Historical confirmation contrast is missing or duplicated")
        contrasts[comparator] = {key: rows[0][key] for key in
            ("contrast", "n", "blocks", "mean_squared_error_improvement", "ci_low95", "ci_high95", "p_holm")}
    return {"source": str(CONFIRMATION), "source_sha256": hashlib.sha256(raw).hexdigest(),
            "state": "recorded aggregate reproduced without new fold-0 analysis", "contrasts": contrasts}


def panel_d(out, allele, result):
    confirmation = recorded_confirmation()
    prior_gain = confirmation["contrasts"]["native_atac_2048"]["mean_squared_error_improvement"]
    names = [("native2k", "Native 2 kb", GRAY), ("fixed_head", "Fixed head 2 kb", GRAY),
             ("adapter", "Adapter 2 kb", BLUE), ("native1m", "Native 1 Mb", GRAY),
             ("combined", "Combined", MAGENTA)]
    fig, (ax, note) = plt.subplots(2, 1, figsize=(3.5, 2.15),
                                    gridspec_kw={"height_ratios": [4, 1]},
                                    layout="constrained")
    source = []
    for j, (column, label, color) in enumerate(names):
        fold_mse = allele.groupby("fold").apply(
            lambda x: float(np.mean((x.beta_alt-x[column])**2)), include_groups=False).to_numpy()
        ax.plot(fold_mse, np.full(len(fold_mse), j), "o", color=color,
                markersize=2, alpha=.55)
        ax.plot(result["MSE"][column], j, "D", color=color, markersize=3)
        source.extend({"model": column, "fold": int(f), "MSE": float(v),
                       "four_fold_pooled_MSE": result["MSE"][column]}
                      for f, v in zip(sorted(allele.fold.unique()), fold_mse))
    pd.DataFrame(source).to_csv(out/"Figure7D_source.tsv", sep="\t", index=False)
    ax.set_yticks(range(len(names)), [x[1] for x in names])
    ax.set_ylim(4.5, -.5)
    ax.set_xlabel("Development MSE, Currin beta squared")
    note.axis("off")
    note.text(0, .7, f"Prior fold 0: adapter MSE gain over native 2 kb {prior_gain:+.5f}", fontsize=6)
    note.text(0, .15, "Confirmation only; 1 Mb comparison not tested there", fontsize=6)
    save(fig, out/"Figure7D_sequence_context.pdf")
    return confirmation


def acquisition_source(acquisition_path, histology_path=None):
    summary = pd.read_csv(acquisition_path/"acquisition_summary.tsv", sep="\t")
    people = pd.read_csv(acquisition_path/"acquisition_per_participant.tsv.gz", sep="\t")
    for table in (summary, people):
        table["policy"] = table.policy.replace({"histology": "histology_equal_count_balancing"})
    if histology_path is not None:
        if not (histology_path/"results.json").is_file():
            raise FileNotFoundError(histology_path/"results.json")
        summary = pd.concat([summary, pd.read_csv(histology_path/"acquisition_summary.tsv", sep="\t")],
                            ignore_index=True)
        people = pd.concat([people, pd.read_csv(histology_path/"acquisition_per_participant.tsv.gz", sep="\t")],
                           ignore_index=True)
    averaged = people.groupby(["budget_fraction", "policy", "participant_id"]).agg(
        base_SSE=("base_SSE", "mean"), refit_SSE=("refit_SSE", "mean")).reset_index()
    pooled = averaged.assign(MSE_reduction=(averaged.base_SSE-averaged.refit_SSE)/1000).groupby(
        ["budget_fraction", "policy"]).agg(participant_pooled_MSE_reduction=("MSE_reduction", "mean"),
                                            n_participants=("participant_id", "nunique")).reset_index()
    if not pooled.n_participants.eq(99).all():
        raise ValueError("Acquisition participant-pooled census differs")
    return summary.merge(pooled, on=["budget_fraction", "policy"], validate="one_to_one")


def panel_e(out, acquisition):
    histology = ("histology_proportional_random" if acquisition.policy.eq("histology_proportional_random").any()
                 else "histology_equal_count_balancing")
    colors = {"random": GRAY, histology: PEACH,
              "rna_diversity": TEAL, "global_local_disagreement": PURPLE}
    labels = {"random": "Random", "histology_equal_count_balancing": "Histology balancing",
              "histology_proportional_random": "Histology stratified",
              "rna_diversity": "RNA diversity", "global_local_disagreement": "Model disagreement"}
    markers = {"random": "o", histology: "s", "rna_diversity": "^",
               "global_local_disagreement": "D"}
    acquisition = acquisition.assign(displayed=acquisition.policy.isin(colors))
    acquisition.to_csv(out/"Figure7E_source.tsv", sep="\t", index=False)
    fig, (cap, refit) = plt.subplots(2, 1, figsize=(3.5, 2.8), sharex=True,
                                    layout="constrained")
    for name in colors:
        x = acquisition.loc[acquisition.policy.eq(name)].sort_values("budget_fraction")
        cap.plot(100*x.budget_fraction, x.mean_error_capture, marker=markers[name], color=colors[name],
                 label=labels[name], markersize=2, linewidth=.8)
        refit.plot(100*x.budget_fraction, x.participant_pooled_MSE_reduction, marker=markers[name], color=colors[name],
                   markersize=2, linewidth=.8)
    for ax in (cap, refit):
        ax.axvline(20, color=GRAY, linestyle=":", linewidth=.6)
    refit.axhline(0, color=GRAY, linewidth=.6)
    cap.set_ylabel("Pool prediction error\ncaptured (fraction)")
    refit.set_ylabel("Held-participant MSE\nreduction (log2 CPM squared)")
    refit.set_xlabel("Additional measurements (n out of 99 participants)")
    cap.legend(frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.02),
               borderaxespad=0)
    refit.set_xticks([10, 20, 30])
    save(fig, out/"Figure7E_measurement_value.pdf")


def panel_f(out, allele):
    identity = "6:42911813:C:T"
    row = allele.loc[allele.key.eq(identity)]
    if len(row) != 1:
        raise ValueError("Prespecified example absent from identical-row population")
    row = row.iloc[0]
    pd.DataFrame([{"key": identity, "measured_beta_alt": row.beta_alt,
                   "native1m_beta_alt": row.native1m, "adapter_beta_alt": row.adapter,
                   "combined_beta_alt": row.combined, "peak_id": "peak282340",
                   "peak_interval": "chr6:42911191-42912028",
                   "target_gene_link": "unsupported",
                   "effect_semantics": "association per additional ALT allele, not an isolated editing effect",
                   "next_experiment": "Test the allele substitution at the variant and measure the source peak; assess target-gene links separately"}]).to_csv(
                       out/"Figure7F_source.tsv", sep="\t", index=False)
    fig, (locus, effect) = plt.subplots(2, 1, figsize=(3.5, 2.3),
                                       gridspec_kw={"height_ratios": [1.1, 2]},
                                       layout="constrained")
    locus.hlines(0, 42.905, 42.970, color=GRAY, linewidth=.5)
    locus.plot([42.911813], [0], marker="o", color=BLUE, markersize=3)
    locus.plot([42.960690, 42.963883], [0, 0], color=MAGENTA, linewidth=3)
    locus.text(42.911813, .13, "caQTL peak / lead", ha="left", fontsize=6)
    locus.text(42.960690, -.2, "GNMT", ha="left", fontsize=6, fontstyle="italic")
    locus.set_xlim(42.905, 42.970)
    locus.set_ylim(-.4, .35)
    locus.set_yticks([])
    locus.set_xlabel("GRCh38 chr6 position (Mb)")
    values = [("Measured", row.beta_alt, GRAY),
              ("Native 1 Mb", row.native1m, GRAY),
              ("Adapter 2 kb", row.adapter, BLUE),
              ("Combined", row.combined, MAGENTA)]
    for j, (label, value, color) in enumerate(values):
        effect.plot(value, j, "o", color=color, markersize=3)
    effect.axvline(0, color=GRAY, linewidth=.5)
    effect.set_yticks(range(len(values)), [x[0] for x in values])
    effect.set_ylim(3.5, -.5)
    effect.set_xlabel("ALT-dose association beta; gene link unresolved")
    save(fig, out/"Figure7F_biological_example.pdf")


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if args.out.exists():
        raise FileExistsError(args.out)
    sequence = args.chrom.parent/"sequence_rna_residual_21849236"
    matched_path = args.chrom.parent/"within_histology_21849320"
    secondary_path = args.chrom.parent/"chromatin_secondary_review"
    histology_path = args.histology_acquisition
    for p in (args.chrom, args.acquisition, args.allele, sequence,
              matched_path, secondary_path, histology_path):
        if not p.exists():
            raise FileNotFoundError(p)
        if not (p/"results.json").is_file():
            raise FileNotFoundError(p/"results.json")
    args.out.mkdir(parents=True)
    chrom = json.loads((args.chrom/"results.json").read_text())
    allele_result = json.loads((args.allele/"results.json").read_text())
    allele = pd.read_csv(args.allele/"matched_development_predictions.tsv.gz", sep="\t")
    acquisition = acquisition_source(args.acquisition, histology_path)
    secondary = json.loads((secondary_path/"results.json").read_text())
    matched = json.loads((matched_path/"results.json").read_text())
    panel_a(args.out, args.acquisition/"fixed_region_panel.tsv")
    panel_b(args.out, allele)
    panel_c(args.out, chrom, sequence)
    confirmation = panel_d(args.out, allele, allele_result)
    panel_e(args.out, acquisition)
    panel_f(args.out, allele)
    pdfs = sorted(args.out.glob("Figure7*.pdf"))
    if len(pdfs) != 6:
        raise ValueError("Expected six individual Figure 7 panels")
    primary = acquisition.loc[np.isclose(acquisition.budget_fraction, .2)].set_index("policy")
    acquisition_folds = pd.read_csv(args.acquisition/"acquisition_per_fold.tsv", sep="\t")
    counts = acquisition_folds[["outer_fold", "n_base", "n_pool", "n_evaluation"]].drop_duplicates()
    if len(counts) != 5 or counts.n_base.nunique() != 1:
        raise ValueError("Acquisition participant census differs")
    c = pd.read_csv(args.out/"Figure7C_source.tsv", sep="\t").set_index("model")
    reviewed_q = [row["p_holm_three_policy_family"] for row in secondary["acquisition"]["contrasts"]]
    if not np.allclose(reviewed_q, reviewed_q[0]):
        raise ValueError("Original acquisition Holm summary requires a revised caption")
    prior_native = confirmation["contrasts"]["native_atac_2048"]
    prior_fixed = confirmation["contrasts"]["frozen_2048_symmetric_shared_mlp64"]
    caption = f"""# Candidate Figure 7. Molecular differences recovered and research use

Candidate only; the Resource's working Figure 7 has not changed. The left column
uses 99 paired GSE267145 biological participants; the right column uses the
same matched Currin development variants in held chromosome folds 1–4, except
for the separately identified spent fold-0 confirmation. Predictions are not
measurements. All primary model estimates are source-exposed development results.

**A.** Measured concentration-residual H3K27ac and held-participant RNA
predictions at 24 regions from an annotation-only, hash-fixed 1,000-region
panel, with participants ordered by fibrosis then steatosis. Both matrices use
the measured region mean and standard deviation solely for display. The model
was fitted on other participants, not on unseen regions. White boundaries mark
recorded fibrosis stages. Colors share one display-z scale clipped at the
98th percentile of measured absolute values.

**B.** Measured Currin FastQTL ALT-dosage beta versus native 1-Mb, five-seed
adapter 2-kb and combined predictions. The plot shows 1,500 of
{len(allele):,} matched variants chosen by an outcome-blind hash. Axes retain
the source beta units; points outside the displayed 99th-percentile magnitude
range are clipped for display, and all rows enter Panel D's MSE.
Currin beta is an association per additional ALT allele in the source molecular
phenotype's units; it is not the effect of an isolated allele-editing experiment.

**C.** Paired improvements in held-participant H3K27ac prediction, with each
contrast's own paired interval. Top: overall skill increments on 96,460 regions
in 99 participants; skill uses the fold-training region mean as its reference.
The global RNA, global-plus-local and training-fold-tuned histology, sex and
RNA-technical models have absolute skills {chrom['skill']['global_RRR']['skill']:.5f},
{chrom['skill']['RRR_local']['skill']:.5f} and
{chrom['skill']['histology_sex_RNA_technical']['skill']:.5f}.
Middle: improvements in predicting differences
within 22 disjoint participant pairs (44 participants) exactly matched on
steatosis, ballooning, lobular inflammation, fibrosis, sex and held fold, with
zero within-pair difference as the reference. Pairs were selected from metadata
without outcomes. Twenty pairs are women at F0, one pair is men at F0 and one
pair is men at F1; this sensitivity does not establish performance within
advanced fibrosis. The local model's within-pair skill is
{matched['skill']['RRR_local']['skill']:.5f}; its increment over the technical
comparator is {matched['contrasts']['RRR_local_vs_technical_nuisance']['delta_skill']:+.5f}
[{matched['contrasts']['RRR_local_vs_technical_nuisance']['ci95'][0]:+.5f},
{matched['contrasts']['RRR_local_vs_technical_nuisance']['ci95'][1]:+.5f}].
Intervals resample participants for overall skill increments and disjoint pairs
for within-pair increments, conditional on the fitted predictions. The three
overall and two within-pair contrasts receive Holm correction in separate
families; pair-by-chromosome resampling is a separate
dependence sensitivity and does not test unseen regions.
Bottom: increment of the rank-selected sequence-by-RNA residual over an
additive control at 1,024 fixed regions, using the original concentration-
residual target reference rather than the remaining baseline-error denominator:
{c.loc['interaction_minus_additive', 'estimate']:+.6f}
[{c.loc['interaction_minus_additive', 'ci_low95']:+.6f},
{c.loc['interaction_minus_additive', 'ci_high95']:+.6f}] by paired participant
bootstrap. This falls below the +0.01 original-skill complexity margin.
The frozen GSE267145 sequence features use `counted_interval` semantics with
the recorded 1-bp coordinate residual; correspondence to called-peak boundaries
remains unsupported. Rank and penalty selection reused cross-fitted base-model
residuals generated once per outer training fold. Outer evaluation participants
remain held out; inner selection does not fully nest the complete base-model
fitting procedure.
The three estimands and two region sets remain on separate axes.

**D.** Same-variant development MSE in squared Currin beta units for native
2 kb, matched fixed head, five-seed adapter 2 kb, native 1 Mb and their
combined prediction. Small circles are four chromosome-fold values; diamonds
are pooled values. Native 1-Mb MSE is {allele_result['MSE']['native1m']:.5f},
adapter 2-kb MSE {allele_result['MSE']['adapter']:.5f}, and combination MSE
{allele_result['MSE']['combined']:.5f}. The fold-0 adapter-versus-native-2-kb
confirmation is printed separately and was not used to fit or choose the
combination. Its recorded adapter MSE improvement over native 2 kb is
{prior_native['mean_squared_error_improvement']:+.5f}
[{prior_native['ci_low95']:+.5f}, {prior_native['ci_high95']:+.5f}], and over the
matched fixed head is {prior_fixed['mean_squared_error_improvement']:+.5f}
[{prior_fixed['ci_low95']:+.5f}, {prior_fixed['ci_high95']:+.5f}]. These values
are read from the completed aggregate summary without reopening fold-0 rows
or recomputing its contrasts. Native 1-Mb MSE minus combined MSE is
{allele_result['contrasts']['combined_vs_native1m']['chromosome']['MSE_reduction']:+.5f}
[{allele_result['contrasts']['combined_vs_native1m']['chromosome']['CI95'][0]:+.5f},
{allele_result['contrasts']['combined_vs_native1m']['chromosome']['CI95'][1]:+.5f}];
adapter MSE minus combined MSE is
{allele_result['contrasts']['combined_vs_adapter']['chromosome']['MSE_reduction']:+.5f}
[{allele_result['contrasts']['combined_vs_adapter']['chromosome']['CI95'][0]:+.5f},
{allele_result['contrasts']['combined_vs_adapter']['chromosome']['CI95'][1]:+.5f}].
The paired intervals resample {allele_result['n_chromosomes']} chromosomes,
conditional on saved predictions; they exclude refitting, recipe selection,
fitted-seed variation and runtime numerical variation. The two primary contrasts
receive Holm correction (P={allele_result['contrasts']['combined_vs_native1m']['p_holm_primary_family']:.5f}
and P={allele_result['contrasts']['combined_vs_adapter']['p_holm_primary_family']:.5f}, respectively).
The combination's relative MSE improvement over the better simple model is
{100*allele_result['combined_relative_MSE_gain_over_best_simple']:.2f}%; the 5%
development complexity margin was {'met' if allele_result['development_complexity_margin_met'] else 'not met'}.
These inspected development results do not establish external superiority.

**E.** Fraction of candidate-pool prediction error found, and reduction in
MSE after selected H3K27ac measurements are revealed, on distinct held
participants. Each outer fold starts with {int(counts.n_base.iloc[0])} measured
training participants, a selectable pool of {int(counts.n_pool.min())}–{int(counts.n_pool.max())}
participants, and {int(counts.n_evaluation.min())}–{int(counts.n_evaluation.max())}
separate evaluation participants. The x-axis shows 10, 20 or 30 additional
measurements out of the full cohort of 99, corresponding to the rounded
10%, 20% and 30% budgets. These are not percentages of the selectable pool.
A nominal 20% budget (20 additional measurements)
is primary; at that budget the random, proportional-histology, RNA-diversity
and model-disagreement participant-pooled refit reductions are
{primary.loc['random', 'participant_pooled_MSE_reduction']:.5f},
{primary.loc['histology_proportional_random', 'participant_pooled_MSE_reduction']:.5f},
{primary.loc['rna_diversity', 'participant_pooled_MSE_reduction']:.5f} and
{primary.loc['global_local_disagreement', 'participant_pooled_MSE_reduction']:.5f}.
Global rank-12 loadings are fitted across all 96,460 regions; local correction
and evaluation use the fixed 1,000-region panel. This bounded refit differs
from the complete released recipe. The upper curves average pool-error capture
over five folds; the lower curves weight each of the 99 distinct evaluation
participants equally. Random and proportional-histology sampling average three
fixed-seed repetitions. Proportional joint-histology quotas use largest-remainder
allocation followed by random selection within strata. The original deterministic
equal-count histology balancing remains in the source table, with a primary
participant-pooled reduction of
{primary.loc['histology_equal_count_balancing', 'participant_pooled_MSE_reduction']:.5f}.
The proportional policy was specified after the initial exposed development
results and is a sensitivity, not independent confirmation.
The independently reconstructed 99-participant conditional bootstrap comparisons
for original balancing, RNA diversity and model disagreement versus random all
have intervals including zero and Holm-adjusted approximate P={reviewed_q[0]:.4f}. The
proportional-histology curve is descriptive here. Earlier fold sign flips are
dependence sensitivities, not exact tests: the five training sets overlap.
All such intervals condition on the fitted models and fixed region panel;
retraining, region-population and selection Monte Carlo uncertainty are omitted.

**F.** A development-selected chr6:42,911,813 C>T Currin lead and its source
peak chr6:42,911,191–42,912,028 (GRCh38) near GNMT. The measurements and
predictions have the same signed Currin-beta orientation. The peak-to-GNMT
target link and any MASLD effect are unresolved. Test the allele substitution
at the variant and measure the source peak; assess target-gene links separately.
The observed beta is a source association per additional ALT allele and does
not establish the isolated causal effect of editing that allele.
"""
    (args.out/"Figure7_candidate_caption.md").write_text(caption)
    manifest = {"status": "candidate_not_adopted", "figure_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "panels": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in pdfs},
                "source_tables": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sorted(args.out.glob("Figure7*_source.tsv"))},
                "data_sources": {"chromatin": str(args.chrom), "acquisition": str(args.acquisition),
                                 "allele": str(args.allele), "sequence_residual": str(sequence),
                                 "within_histology": str(matched_path), "secondary_review": str(secondary_path),
                                 "histology_sampling": str(histology_path)},
                "caption_sha256": hashlib.sha256((args.out/"Figure7_candidate_caption.md").read_bytes()).hexdigest(),
                "recorded_fold0_confirmation": confirmation,
                "bounds": ["GSE267145 is one inspected development cohort",
                           "Currin is one inspected source; fold 0 confirmation is separate",
                           "Wenz source not admitted; no external allele-effect transfer shown",
                           "Target assignment for Figure 7F is unresolved"]}
    (args.out/"manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"pdfs": [p.name for p in pdfs], "candidate": True}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", type=Path, required=True)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--allele", type=Path, required=True)
    parser.add_argument("--histology-acquisition", type=Path,
                        default=ROOT/"GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/histology_sampling_21849335")
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
