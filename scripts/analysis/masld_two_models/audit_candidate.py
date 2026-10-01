#!/usr/bin/env python3
"""Second-path numerical checks of the candidate model and figure quantities."""
import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT/"Analysis/MASLD_Model_Benchmark"
FIX = BENCH/"executions/model-data-064-21079902/fixture"
OOF = BENCH/"executions/chromatin-stable-rrr-forms-20260908T192024Z/out/stable_oof.npz"
TECH = BENCH/"executions/chromatin-stable-rrr-forms-20260908T192024Z/out/technical_features.tsv"
JOIN = BENCH/"executions/gse267145-authoritative-join-21064930/participant_join.tsv"


def check_close(actual, expected, label, tol=1e-6):
    if not np.isfinite([actual, expected]).all() or abs(actual-expected) > tol:
        raise ValueError(f"{label}: {actual} differs from {expected} by >{tol}")


def chromatin_check(path):
    result = json.loads((path/"results.json").read_text())
    saved = np.load(OOF, allow_pickle=False)
    y = saved["Cres_oof"].astype(float)
    fold = saved["fold"]
    errors = {}
    denom = np.sum(y*y)
    for label, field in (("global_RRR", "pred_rrr_wide_stable"),
                         ("RRR_local", "pred_rrr_offset_cis_wide")):
        skill = 1-np.sum((y-saved[field])**2)/denom
        check_close(skill, result["skill"][label]["skill"], f"chromatin {label}", 1e-6)
        errors[label] = float(skill)
    axis = pd.read_csv(FIX/"molecular/participant_axis.tsv", sep="\t")
    join = pd.read_csv(JOIN, sep="\t").set_index("participant_id").loc[axis.participant_id]
    hist = join[["steatosis", "ballooning", "lobular_inflammation", "fibrosis"]].to_numpy(float)
    sex = (join.sex.to_numpy() == "F").astype(float)[:, None]
    tech = pd.read_csv(TECH, sep="\t")
    x = np.column_stack([hist, sex, tech[[c for c in tech if not c.startswith("h3k27ac_")]].to_numpy(float)])
    h3 = np.asarray(np.load(FIX/"molecular/h3k27ac_counts.npy", mmap_mode="r"), float)
    total = h3.sum(1)
    ym = np.log2(h3/total[:, None]*1e6+1)
    ratio = h3/total[:, None]
    entropy = -np.sum(np.where(ratio > 0, ratio*np.log(np.maximum(ratio, 1e-300)), 0), axis=1)
    sorted_counts = np.sort(h3, axis=1)
    width = h3.shape[1]
    gini = ((2*np.arange(1, width+1)-width-1)*sorted_counts).sum(1)/(width*total)
    iqr = np.array([np.subtract(*np.percentile(np.log2(v[v > 0]), [75, 25])) for v in h3])
    conc = np.column_stack([gini, entropy, iqr,
                            sorted_counts[:, -int(round(.05*width)):].sum(1)/total])
    predictions = np.zeros_like(y)
    for f in sorted(set(fold)):
        tr, te = fold != f, fold == f
        basis = np.column_stack([np.ones(tr.sum()), conc[tr]])
        coeff = np.linalg.pinv(basis)@ym[tr]
        target = ym[tr]-basis@coeff
        mu, sd = x[tr].mean(0), x[tr].std(0)
        sd[sd < 1e-9] = 1
        a, b = (x[tr]-mu)/sd, (x[te]-mu)/sd
        alpha = np.linalg.eigvalsh(a.T@a).max()
        model = Ridge(alpha=float(alpha), fit_intercept=True, solver="svd")
        model.fit(a, target)
        predictions[te] = model.predict(b)
    hist_skill = 1-np.sum((y-predictions)**2)/denom
    check_close(hist_skill, result["skill"]["histology_sex_RNA_technical"]["skill"],
                "histology/technical independent ridge", 2e-5)
    errors["histology_sex_RNA_technical"] = float(hist_skill)
    delta = errors["RRR_local"]-hist_skill
    check_close(delta, result["contrasts"]["RRR_local_vs_histology_sex_RNA_technical"]["delta_skill"],
                "chromatin paired contrast", 2e-5)
    return errors


def acquisition_check(path):
    folded = pd.read_csv(path/"acquisition_per_fold.tsv", sep="\t")
    donor = pd.read_csv(path/"acquisition_per_participant.tsv.gz", sep="\t")
    summary = pd.read_csv(path/"acquisition_summary.tsv", sep="\t")
    if len(folded) != 90 or folded.outer_fold.nunique() != 5 or donor.participant_id.nunique() != 99:
        raise ValueError("Acquisition census differs from design")
    for keys, part in donor.groupby(["outer_fold", "budget_fraction", "policy", "random_repeat"]):
        f, budget, policy, repeat = keys
        original = folded.loc[folded.outer_fold.eq(f) & folded.budget_fraction.eq(budget)
                              & folded.policy.eq(policy) & folded.random_repeat.eq(repeat)]
        if len(original) != 1 or part.participant_id.duplicated().any():
            raise ValueError("Acquisition fold/policy identity is not unique")
        check_close(part.base_SSE.mean()/1000, original.iloc[0].evaluation_base_MSE, "acquisition base MSE")
        check_close(part.refit_SSE.mean()/1000, original.iloc[0].evaluation_refit_MSE, "acquisition refit MSE")
        check_close((part.base_SSE-part.refit_SSE).mean()/1000,
                    original.iloc[0].evaluation_MSE_reduction, "acquisition improvement")
    derived = folded.groupby(["budget_fraction", "policy"], sort=True).agg(
        mean_error_capture=("captured_pool_error_fraction", "mean"),
        mean_refit_MSE_reduction=("evaluation_MSE_reduction", "mean")).reset_index()
    merged = derived.merge(summary, on=["budget_fraction", "policy"], suffixes=("_audit", "_source"))
    if len(merged) != 12:
        raise ValueError("Missing acquisition budget/policy family")
    for field in ("mean_error_capture", "mean_refit_MSE_reduction"):
        if np.max(np.abs(merged[field+"_audit"]-merged[field+"_source"])) > 1e-10:
            raise ValueError(f"Acquisition summary {field} differs")
    return {"n_fit_comparisons": len(folded), "n_unique_evaluation_participants": donor.participant_id.nunique(),
            "primary_20pct": summary.loc[summary.budget_fraction.eq(.2),
                                          ["policy", "mean_error_capture", "mean_refit_MSE_reduction"]].to_dict("records")}


def allele_check(path):
    table = pd.read_csv(path/"matched_development_predictions.tsv.gz", sep="\t")
    result = json.loads((path/"results.json").read_text())
    if table.fold.eq(0).any() or table.key.duplicated().any() or table.chr.nunique() != result["n_chromosomes"]:
        raise ValueError("Allele fold, identity or chromosome count disagrees")
    mse = {}
    for model in ("native2k", "fixed_head", "adapter", "native1m", "combined"):
        mse[model] = float(np.mean((table.beta_alt-table[model])**2))
        check_close(mse[model], result["MSE"][model], f"allele MSE {model}", 1e-10)
    for weight in result["weights"]:
        part = table.loc[table.fold.eq(weight["held_fold"])]
        predicted = part.native1m*weight["native1m_weight"]+part.adapter*weight["adapter_weight"]
        if np.max(np.abs(predicted-part.combined)) > 1e-9:
            raise ValueError("Combined effect differs from stored fit weights")
    stacking = stacking_check(path, result)
    p_values = []
    for baseline in ("native1m", "adapter"):
        observed = mse[baseline]-mse["combined"]
        declared = result["contrasts"]["combined_vs_"+baseline]["chromosome"]["MSE_reduction"]
        check_close(observed, declared, f"allele contrast combined vs {baseline}", 1e-10)
        losses = table.assign(loss_difference=(table.beta_alt-table[baseline])**2
                              -(table.beta_alt-table.combined)**2)
        for label, grouping in (("chromosome", "chr"), ("historical_bin_sensitivity", "historical_bin")):
            summary = result["contrasts"]["combined_vs_"+baseline][label]
            grouped = losses.groupby(grouping, sort=True).loss_difference.agg(["sum", "size"])
            if len(grouped) != summary["n_groups"]:
                raise ValueError("Allele uncertainty grouping differs")
            # Resample grouped loss sums and row counts rather than raw rows.
            rng = np.random.default_rng(summary["seed"])
            sampled = rng.integers(0, len(grouped), (summary["draws"], len(grouped)))
            draws = grouped["sum"].to_numpy()[sampled].sum(1)/grouped["size"].to_numpy()[sampled].sum(1)
            for value, expected in zip(np.quantile(draws, [.025, .975]), summary["CI95"]):
                check_close(value, expected, f"allele {label} paired interval", 1e-10)
            centered = draws-draws.mean()
            p = (1+np.count_nonzero(np.abs(centered) >= abs(observed)))/(len(draws)+1)
            check_close(p, summary["p_nominal_centred_bootstrap"], f"allele {label} bootstrap P", 1e-10)
            if label == "chromosome":
                p_values.append(p)
    order = np.argsort(p_values)
    adjusted = np.empty(2)
    adjusted[order] = np.minimum(1, np.maximum.accumulate(np.array(p_values)[order]*[2, 1]))
    for baseline, p in zip(("native1m", "adapter"), adjusted):
        check_close(p, result["contrasts"]["combined_vs_"+baseline]["p_holm_primary_family"],
                    "allele primary Holm P", 1e-10)
    best_simple = min(mse["adapter"], mse["native1m"])
    relative_gain = (best_simple-mse["combined"])/best_simple
    check_close(relative_gain, result["combined_relative_MSE_gain_over_best_simple"], "allele relative gain", 1e-10)
    if bool(relative_gain >= .05) != result["development_complexity_margin_met"]:
        raise ValueError("Allele complexity decision differs from declared margin")
    return {"matched_variants": len(table), "MSE": mse, "stacking": stacking}


def stacking_check(path, result):
    """Rebuild the stacking fits with sklearn from the sixty saved inputs."""
    source = ROOT/"GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model/native/native_1048576.tsv"
    membership = pd.read_csv(source, sep="\t", usecols=["heldout_fold"])
    skip = np.flatnonzero(~membership.heldout_fold.isin([1, 2, 3, 4]).to_numpy())+1
    native = pd.read_csv(source, sep="\t", skiprows=skip,
                         usecols=["key", "heldout_fold", "beta_alt", "local_atac_liver"]).set_index("key")
    if native.index.has_duplicates or not native.heldout_fold.isin([1, 2, 3, 4]).all():
        raise ValueError("Native stacking audit population differs")
    completion_paths = [Path(value) for value in result["inputs_sha256"] if Path(value).name == "completion.json"]
    if len(completion_paths) != 1:
        raise ValueError("No unique nested-input completion record")
    inputs = completion_paths[0].parent
    count = 0
    for weight in result["weights"]:
        outer = weight["held_fold"]
        matrices, labels = [], []
        for meta in sorted({1, 2, 3, 4}-{outer}):
            training = sorted({1, 2, 3, 4}-{outer, meta})
            seeds = []
            for seed in range(20260917, 20260922):
                folder = inputs/f"outer{outer}_meta{meta}_seed{seed}"
                split = json.loads((folder/"split.json").read_text())
                if (split["training_folds"] != training or split["validation_fold"] != meta
                        or split["held_fold"] != outer or split["seed"] != seed):
                    raise ValueError("Stacking input fit violates its excluded folds or seed")
                if json.loads((folder/"feasibility.json").read_text())["completed_steps"] != 5000:
                    raise ValueError("Incomplete stacking input fit")
                frame = pd.read_csv(folder/"validation_predictions.tsv", sep="\t")
                frame.index = frame.variant_id.str.removeprefix("chr")
                if frame.index.has_duplicates or not frame.validation_fold.eq(meta).all():
                    raise ValueError("Stacking input validation population differs")
                if seeds and set(frame.index) != set(seeds[0].index):
                    raise ValueError("Stacking input seed populations differ")
                seeds.append(frame)
                count += 1
            keys = seeds[0].index.intersection(native.index[native.heldout_fold.eq(meta)])
            calibration = native.loc[native.heldout_fold.isin(training)]
            x = calibration.local_atac_liver.to_numpy(float)
            slope = np.sum(x*calibration.beta_alt.to_numpy(float))/np.sum(x*x)
            predictions = np.column_stack([frame.loc[keys, "predicted_beta"].to_numpy(float) for frame in seeds])
            y = seeds[0].loc[keys, "observed_beta"].to_numpy(float)
            if any(not np.allclose(y, frame.loc[keys, "observed_beta"], rtol=1e-6, atol=1e-6) for frame in seeds):
                raise ValueError("Stacking input seed outcomes differ")
            if not np.allclose(y, native.loc[keys, "beta_alt"], rtol=1e-6, atol=1e-6):
                raise ValueError("Stacking input labels differ from source")
            matrices.append(np.column_stack([slope*native.loc[keys, "local_atac_liver"], predictions.mean(1)]))
            labels.append(y)
        matrix, target = np.vstack(matrices), np.concatenate(labels)
        scale = np.sqrt(np.mean(matrix**2, axis=0))
        standardized = matrix/scale
        alpha = .01*np.linalg.norm(standardized, ord=2)**2
        model = Ridge(alpha=alpha, fit_intercept=False, solver="svd").fit(standardized, target)
        coefficients = model.coef_/scale
        for actual, name in zip(coefficients, ("native1m_weight", "adapter_weight")):
            check_close(actual, weight[name], "Independent stacking coefficient", 1e-10)
    if count != 60:
        raise ValueError("Stacking fit census differs")
    return {"input_fits_checked": count, "independent_coefficient_fits": 4,
            "fold0_outcome_rows_parsed": 0, "solver": "sklearn Ridge SVD"}


def within_histology_check(path):
    pairs = pd.read_csv(path/"pairs_fixed_from_metadata.tsv", sep="\t")
    source = json.loads((path/"results.json").read_text())
    ids = pd.read_csv(FIX/"molecular/participant_axis.tsv", sep="\t").participant_id.astype(str)
    index = {value: j for j, value in enumerate(ids)}
    a, b = pairs.participant_a.map(index).to_numpy(int), pairs.participant_b.map(index).to_numpy(int)
    if len(pairs) != 22 or len(set(a)|set(b)) != 44:
        raise ValueError("Matched pair replication census differs")
    join = pd.read_csv(JOIN, sep="\t").set_index("participant_id").loc[ids]
    fields = ["outer_fold", "steatosis", "ballooning", "lobular_inflammation", "fibrosis", "sex"]
    if not np.array_equal(join.iloc[a][fields].to_numpy(), join.iloc[b][fields].to_numpy()):
        raise ValueError("Matched pair differs in recorded histology, sex or fold")
    with np.load(OOF, allow_pickle=False) as saved:
        measured = saved["Cres_oof"].astype(float)
        truth = measured[a]-measured[b]
        denominator = np.sum(truth**2)
        reconstructed = {}
        for key, field in (("global_RRR", "pred_rrr_wide_stable"), ("RRR_local", "pred_rrr_offset_cis_wide")):
            pred = saved[field].astype(float)
            skill = float(1-np.sum((truth-(pred[a]-pred[b]))**2)/denominator)
            check_close(skill, source["skill"][key]["skill"], "within-histology skill", 1e-10)
            reconstructed[key] = skill
    return {"n_pairs": 22, "n_distinct_participants": 44,
            "fibrosis_grades": sorted(join.iloc[np.r_[a, b]].fibrosis.unique().astype(int).tolist()),
            "skill": reconstructed}


def figure_check(path, chrom, acquisition, allele):
    manifest = json.loads((path/"manifest.json").read_text())
    if len(manifest["panels"]) != 6 or len(manifest["source_tables"]) != 6:
        raise ValueError("Incomplete six-panel candidate")
    caption = path/"Figure7_candidate_caption.md"
    if not caption.exists() or len(caption.read_text().strip()) < 1000:
        raise ValueError("Figure 7 candidate caption is missing or incomplete")
    if hashlib.sha256(caption.read_bytes()).hexdigest() != manifest.get("caption_sha256"):
        raise ValueError("Figure 7 caption differs from the rendering manifest")
    for pdf in sorted(path.glob("Figure7*.pdf")):
        info = subprocess.run(["pdfinfo", str(pdf)], check=True, text=True, capture_output=True).stdout
        pages = [line.split(":", 1)[1].strip() for line in info.splitlines() if line.startswith("Pages:")]
        if pages != ["1"]:
            raise ValueError(f"Panel must be one PDF page: {pdf}")
        if b"/Subtype /Type3" in pdf.read_bytes():
            raise ValueError(f"Type 3 font in {pdf}")
        text = subprocess.run(["pdftotext", str(pdf), "-"], check=True,
                              text=True, capture_output=True).stdout
        if not text.strip():
            raise ValueError(f"Panel has no editable text: {pdf}")
    c = pd.read_csv(path/"Figure7C_source.tsv", sep="\t").set_index("model")
    published = json.loads((chrom/"results.json").read_text())
    matched = json.loads((chrom.parent/"within_histology_21849320"/"results.json").read_text())
    for key, row in c.iterrows():
        if key in published["contrasts"]:
            check_close(row.estimate, published["contrasts"][key]["delta_skill"], "Figure 7C paired contrast", 1e-10)
            for name, expected in zip(("ci_low95", "ci_high95"), published["contrasts"][key]["ci95"]):
                check_close(row[name], expected, "Figure 7C paired interval", 1e-10)
        elif key == "interaction_minus_additive":
            residual = pd.read_csv(chrom.parent/"sequence_rna_residual_21849236"/
                                   "held_participant_errors.tsv", sep="\t")
            region = pd.read_csv(chrom.parent/"chromatin_secondary_review"/"sequence_regions.tsv", sep="\t").region_index.to_numpy(int)
            with np.load(OOF, allow_pickle=False) as saved:
                denominator = np.sum(saved["Cres_oof"][:, region].astype(float)**2)
            value = (residual.additive_SSE.sum()-residual.interaction_SSE.sum())/denominator
            check_close(row.estimate, value, "Figure 7C original-target skill increment", 1e-8)
        elif key.startswith("matched_"):
            model = key.removeprefix("matched_")
            check_close(row.estimate, matched["contrasts"][model]["delta_skill"], "Figure 7C matched-pair contrast", 1e-10)
            for name, expected in zip(("ci_low95", "ci_high95"), matched["contrasts"][model]["ci95"]):
                check_close(row[name], expected, "Figure 7C matched-pair interval", 1e-10)
        else:
            raise ValueError(f"Unexpected Figure 7C model {key}")
    d = pd.read_csv(path/"Figure7D_source.tsv", sep="\t")
    a = json.loads((allele/"results.json").read_text())
    d_predictions = pd.read_csv(allele/"matched_development_predictions.tsv.gz", sep="\t")
    expected_pairs = {(name, fold) for name in a["MSE"] for fold in (1, 2, 3, 4)}
    if len(d) != 20 or set(zip(d.model, d.fold)) != expected_pairs:
        raise ValueError("Figure 7D model/fold census differs")
    for key, part in d.groupby("model"):
        for _, row in part.iterrows():
            check_close(row.four_fold_pooled_MSE, a["MSE"][key], "Figure 7D pooled MSE", 1e-10)
            population = d_predictions.loc[d_predictions.fold.eq(row.fold)]
            expected = np.mean((population.beta_alt-population[key])**2)
            check_close(row.MSE, expected, "Figure 7D fold MSE", 1e-10)
    historical = manifest["recorded_fold0_confirmation"]
    recorded_bytes = Path(historical["source"]).read_bytes()
    if hashlib.sha256(recorded_bytes).hexdigest() != historical["source_sha256"]:
        raise ValueError("Recorded historical aggregate summary checksum differs")
    recorded = json.loads(recorded_bytes)
    for comparator, values in historical["contrasts"].items():
        matches = [entry for entry in recorded["contrasts"]
                   if entry["arm"] == "adapter_r16_last5_symmetric" and entry["comparator"] == comparator]
        if len(matches) != 1 or any(matches[0][name] != value for name, value in values.items()):
            raise ValueError("Historical annotation differs from completed aggregate summary")
    printed_d = subprocess.run(["pdftotext", str(path/"Figure7D_sequence_context.pdf"), "-"],
                               check=True, text=True, capture_output=True).stdout
    historical_gain = historical["contrasts"]["native_atac_2048"]["mean_squared_error_improvement"]
    if f"{historical_gain:+.5f}" not in printed_d:
        raise ValueError("Figure 7D printed historical number differs from recorded summary")
    e = pd.read_csv(path/"Figure7E_source.tsv", sep="\t")
    original = pd.read_csv(acquisition/"acquisition_summary.tsv", sep="\t")
    original["policy"] = original.policy.replace({"histology": "histology_equal_count_balancing"})
    histology = acquisition.parent/"histology_sampling_21849335"
    original = pd.concat([original, pd.read_csv(histology/"acquisition_summary.tsv", sep="\t")], ignore_index=True)
    comparison = e.merge(original, on=["budget_fraction", "policy"], suffixes=("_figure", "_source"), validate="one_to_one")
    if len(e) != 15 or len(comparison) != 15:
        raise ValueError("Figure 7E policy/budget population differs")
    for name in ("mean_error_capture", "mean_refit_MSE_reduction", "mean_base_MSE"):
        if not np.allclose(comparison[name+"_figure"], comparison[name+"_source"], rtol=0, atol=1e-10):
            raise ValueError("Figure 7E descriptive fold mean differs: "+name)
    people = pd.read_csv(acquisition/"acquisition_per_participant.tsv.gz", sep="\t")
    people["policy"] = people.policy.replace({"histology": "histology_equal_count_balancing"})
    people = pd.concat([people, pd.read_csv(histology/"acquisition_per_participant.tsv.gz", sep="\t")], ignore_index=True)
    for _, entry in e.iterrows():
        part = people.loc[people.policy.eq(entry.policy) & people.budget_fraction.eq(entry.budget_fraction)]
        repeated = part.groupby("participant_id")[["base_SSE", "refit_SSE"]].mean()
        if len(repeated) != 99 or entry.n_participants != 99:
            raise ValueError("Figure 7E participant census differs")
        value = float((repeated.base_SSE-repeated.refit_SSE).sum()/(99*1000))
        check_close(entry.participant_pooled_MSE_reduction, value, "Figure 7E participant-weighted refit value", 1e-10)
    if e.loc[e.displayed, "policy"].nunique() != 4 or e.loc[e.policy.eq("histology_equal_count_balancing"), "displayed"].any():
        raise ValueError("Figure 7E displayed policy set differs")
    matched_predictions = pd.read_csv(allele/"matched_development_predictions.tsv.gz", sep="\t").set_index("key")
    b = pd.read_csv(path/"Figure7B_source.tsv", sep="\t").set_index("key")
    for name in ("beta_alt", "native1m", "adapter", "combined"):
        if len(b) != 1500 or not np.allclose(b[name], matched_predictions.loc[b.index, name], rtol=0, atol=1e-10):
            raise ValueError("Figure 7B effects differ from evaluation table")
    f = pd.read_csv(path/"Figure7F_source.tsv", sep="\t").iloc[0]
    for figure_name, column in (("measured_beta_alt", "beta_alt"), ("native1m_beta_alt", "native1m"),
                                ("adapter_beta_alt", "adapter"), ("combined_beta_alt", "combined")):
        check_close(f[figure_name], matched_predictions.loc[f.key, column], "Figure 7F effect", 1e-10)
    source_a = pd.read_csv(path/"Figure7A_source.tsv", sep="\t")
    axis = pd.read_csv(FIX/"molecular/participant_axis.tsv", sep="\t").participant_id.astype(str)
    indices = {person: i for i, person in enumerate(axis)}
    people = source_a.participant_id.map(indices).to_numpy(int)
    regions = source_a.region_index.to_numpy(int)
    if len(source_a) != 99*24 or len(set(regions)) != 24:
        raise ValueError("Figure 7A census differs")
    with np.load(OOF, allow_pickle=False) as saved:
        truth = saved["Cres_oof"].astype(float)
        mean, sd = truth.mean(0), truth.std(0)
        sd[sd < 1e-10] = 1
        for column, field in (("measured_display_z", "Cres_oof"), ("predicted_display_z", "pred_rrr_offset_cis_wide")):
            actual = (saved[field][people, regions].astype(float)-mean[regions])/sd[regions]
            if not np.allclose(actual, source_a[column], rtol=0, atol=1e-10):
                raise ValueError("Figure 7A display values differ from saved held-participant predictions")
    return {"panels": 6, "vector_pdf_text_checked": True, "type3_fonts": 0}


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if args.out.exists():
        raise FileExistsError(args.out)
    report = {"status": "independent_candidate_reconstruction",
              "chromatin": chromatin_check(args.chrom),
              "within_histology": within_histology_check(args.chrom.parent/"within_histology_21849320"),
              "acquisition": acquisition_check(args.acquisition),
              "allele": allele_check(args.allele),
              "figure": figure_check(args.figure, args.chrom, args.acquisition, args.allele)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for name in ("chrom", "acquisition", "allele", "figure", "out"):
        parser.add_argument("--"+name, type=Path, required=True)
    main(parser.parse_args())
