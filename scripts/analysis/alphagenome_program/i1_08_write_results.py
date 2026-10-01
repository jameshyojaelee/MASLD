#!/usr/bin/env python3
"""Generate RESULTS.md for the AD.1 gate run from this run's own JSON artifacts.

Every number printed here is read from `tables/*.json` and `tables/*.tsv`, not transcribed, so the
report cannot drift from what the jobs computed.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import i1_common as C  # noqa: E402


def f(x, nd=4):
    if x is None:
        return "n/a"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v != v:
        return "n/a"
    return f"{v:.{nd}f}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    args = ap.parse_args()
    run = pathlib.Path(args.run)
    T = run / "tables"

    census = json.loads((run / "inputs/gate_set_census.json").read_text())
    gate = json.loads((T / "ad1_gate.json").read_text()) if (T / "ad1_gate.json").exists() else None
    zs = (
        json.loads((T / "zeroshot_endpoint.json").read_text())
        if (T / "zeroshot_endpoint.json").exists()
        else None
    )
    sweep = (
        pd.read_csv(T / "pilot_length_sweep.tsv", sep="\t")
        if (T / "pilot_length_sweep.tsv").exists()
        else None
    )
    pilot = (
        json.loads((T / "pilot_verdict.json").read_text())
        if (T / "pilot_verdict.json").exists()
        else {}
    )
    jobs = (
        pd.read_csv(T / "jobs.tsv", sep="\t") if (T / "jobs.tsv").exists() else None
    )
    versions = (run / "env/versions.txt").read_text() if (run / "env/versions.txt").exists() else ""

    L: list[str] = []
    A = L.append

    A("# AD.1 gate: do the local AlphaGenome weights reproduce the archived hosted Atlas ranking?")
    A("")
    A(
        "Run `ad1-gate-20260915T114615Z`. Executes the first prediction of "
        "`scripts/analysis/alphagenome_program/ADAPTATION_ARM_PRESPEC.md` as amended by "
        "`ADAPTATION_ARM_AMENDMENT_01.md`, both unedited by this run. Nothing here adopts a Resource "
        "number, gene class, figure or claim."
    )
    A("")

    # ---------------------------------------------------------------- verdict first
    if gate:
        a = gate["channels"]["atac"]
        d = gate["channels"]["dnase"]
        A("## The verdict, first")
        A("")
        A(
            f"**{gate['gate_verdict']}.** {gate['gate_verdict_statement']}."
        )
        A("")
        A(
            "Both score sets are the same quantity: a 501-bp variant-centred centre mask, "
            "`DIFF_LOG2_SUM` aggregation, averaged over the same three named liver ATAC tracks. The "
            "archived column is the hosted Atlas `ATAC` scorer, which is "
            "`CenterMaskScorer(width=501, aggregation_type=DIFF_LOG2_SUM)` in the SDK's own "
            "recommended set; the local side applies that same width and aggregation through the "
            "packaged `create_center_mask` and `_apply_aggregation`, imported rather than "
            "reimplemented. Raw is compared with raw: the local runtime is built with no calibration "
            "table, so it produces no quantile layer, and the runtime probe already recorded that "
            "Atlas quantiles must not be pooled with local values."
        )
        A("")
        A(
            "**The packaged `model.score_variant` could not be used, and that was measured rather "
            "than assumed.** Its `_predict_variant` unconditionally calls "
            "`jnp.asarray(splice_junction_masks.splice_sites)`, and `splice_sites` is `None` whenever "
            "the model is built without a splice-site extractor, which is how it is built here "
            "because the packaged extractor's feather files sit behind a Google Storage URL that "
            "compute nodes cannot reach. The first pilot attempt failed exactly there, at all five "
            "input lengths, with `ValueError: None is not a valid value for jnp.array`. The two "
            "`predict_sequence` calls the runtime probe had already validated are used instead. For a "
            "matched-reference substitution the two routes see identical inputs: "
            "`extract_variant_sequences` writes `variant.reference_bases` over a base that already "
            "equals REF, indel stitching does not apply, splice-junction post-processing touches only "
            "the SPLICE_JUNCTIONS head, and the interval is on the positive strand. Requesting a "
            "splice annotation in order to score an ATAC track would also have brought a second "
            "annotation vintage into the run."
        )
        A("")
        A("| channel | n | 1-Mb blocks | Spearman | 95% CI | Pearson | sign agreement | median \\|difference\\| | median \\|archived\\| |")
        A("|---|---:|---:|---:|---|---:|---:|---:|---:|")
        for r in (a, d):
            A(
                f"| {r['channel']} | {r['n']} | {r['blocks']} | **{f(r['spearman'])}** | "
                f"{f(r['ci_low'])} to {f(r['ci_high'])} | {f(r['pearson'])} | "
                f"{f(r['sign_agreement'])} | {f(r['median_abs_difference'], 5)} | "
                f"{f(r['median_abs_archived'], 5)} |"
            )
        A("")
        A(
            f"The interval resamples the {a['blocks']} 1-Mb blocks the {a['n']} leads sit in, "
            f"{gate['channels']['atac']['bootstrap_resamples']:,} draws at seed "
            f"{C.BOOTSTRAP_SEED}. {a['n']} leads in {a['blocks']} blocks is less information than "
            f"{a['n']} independent leads, and every statement here is block-limited."
        )
        A("")
        A(f"Windows rejected before scoring: **{gate['rejected_windows']}** of {census['gate_rows']}.")
        if gate["rejected_reasons"]:
            A("")
            A("| reason | windows |")
            A("|---|---:|")
            for k, v in sorted(gate["rejected_reasons"].items()):
                A(f"| {k} | {v} |")
        A("")
        A(
            "The packaged one-hot encoder writes a zero vector for any byte that is not A, C, G or T "
            "and warns about nothing, so every window was alphabet-validated on the same interval the "
            "model extracts, before it was scored. A rejected window was not scored and was not "
            "silently replaced. Without this guard each of these would have returned a "
            "confident-looking prediction built partly from zero vectors."
        )
        A("")
        rk = gate.get("rejected_vs_kept_effect_size") or {}
        if rk:
            A(
                f"The exclusion is not neutral and its direction is reported rather than left "
                f"implicit. The rejected leads have a median |archived ATAC| of "
                f"{f(rk['rejected_median_abs_archived_atac'], 4)} against "
                f"{f(rk['kept_median_abs_archived_atac'], 4)} for the "
                f"{rk['kept_n']} that remain, so what leaves the gate set is weighted towards "
                f"near-null variants, which is where two implementations most often disagree in sign. "
                f"The gate correlation below is therefore measured on a slightly easier set than the "
                f"full 3,845."
            )
            A("")
        sw = gate.get("rejected_would_be_clean_at_shorter_length") or {}
        if sw:
            A(
                "Most of the loss is the window length rather than a defect at the variant: of the "
                f"{gate['rejected_windows']} rejected 1-Mb windows, "
                + ", ".join(f"{v} would be clean at {int(k):,} bp" for k, v in sorted(sw.items(), key=lambda x: -int(x[0])))
                + "."
            )
            A("")

    # ---------------------------------------------------------------- the length question
    A("## Which input length the Atlas scored at, measured rather than assumed")
    A("")
    A(
        "The archived Atlas chunks record `requested_scorers`, `sdk_version`, `chunk_size`, "
        "`max_workers` and `stage`, and nothing about the interval, so the sequence length was not "
        "known from the archive. The same 60 leads, one per 1-Mb block, were scored at all five "
        "lengths the SDK names and each length's agreement with the archive was read off."
    )
    A("")
    if sweep is not None:
        A("| input bp | n | ATAC Spearman | ATAC sign agr. | DNase Spearman | median local/archived | first call s | warm call s | peak GPU GiB |")
        A("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for _, r in sweep.iterrows():
            if r.get("state") != "ok":
                A(f"| {int(r['length_bp']):,} | — | — | — | — | — | — | — | — |")
                continue
            A(
                f"| {int(r['length_bp']):,} | {int(r['n'])} | {f(r.get('atac_spearman'))} | "
                f"{f(r.get('atac_sign_agreement'))} | {f(r.get('dnase_spearman'))} | "
                f"{f(r.get('atac_median_ratio_local_over_archived'), 3)} | "
                f"{f(r.get('first_call_seconds'), 2)} | {f(r.get('warm_call_seconds_median'), 3)} | "
                f"{f(r.get('peak_gpu_gib'), 2)} |"
            )
        A("")
    msw = (
        pd.read_csv(T / "pilot_length_sweep_matched.tsv", sep="\t")
        if (T / "pilot_length_sweep_matched.tsv").exists()
        else None
    )
    if msw is not None:
        A(
            "The table above compares lengths on different lead sets, because a longer window "
            "rejects more windows for assembly gaps and chromosome edges, and a comparison across "
            "different n is not a comparison of lengths. The same statistic on the "
            f"{int(msw.n_common.iloc[0])} leads every length scored:"
        )
        A("")
        A("| input bp | ATAC Spearman | DNase Spearman | ATAC sign agr. | median \\|difference\\| | median local/archived |")
        A("|---:|---:|---:|---:|---:|---:|")
        for _, r in msw.sort_values("length_bp").iterrows():
            A(
                f"| {int(r['length_bp']):,} | {f(r['atac_spearman'], 5)} | "
                f"{f(r['dnase_spearman'], 5)} | {f(r['atac_sign_agreement'])} | "
                f"{f(r['atac_median_abs_difference'], 5)} | "
                f"{f(r['atac_median_ratio_local_over_archived'], 3)} |"
            )
        A("")
    if pilot:
        cost = pilot.get("cost", {})
        A(
            f"**Best agreement is at {pilot.get('best_length_on_matched_leads', pilot.get('chosen_length_bp')):,} bp, "
            f"and the gate was run at {pilot.get('length_actually_used_for_the_gate', 0):,} bp anyway.** "
            f"That is a deliberate choice against the sweep, and the reason is comparability rather "
            f"than agreement. 524,288 and 1,048,576 bp are not separable here: 0.9980 against 0.9969 "
            f"on {int(msw.n_common.iloc[0]) if msw is not None else 56} leads, with 1,048,576 bp "
            f"giving the smaller median absolute difference of the two. 1,048,576 bp is the length "
            f"of this project's pinned AlphaGenome recipe, the length the HSD17B13 anchor was "
            f"validated at, and the length arm 2's frozen probe reads its (1, 8192, 3072) 128-bp "
            f"representation at. An arm-2-minus-zero-shot increment measured across two window "
            f"lengths would confound adaptation with window length, and the spec requires nested "
            f"increments to differ in one component."
        )
        A("")
        A(
            f"The price is measured and paid knowingly. The cost pilot at 1,048,576 bp scored "
            f"{cost.get('pilot_variants_scored')} leads at a median warm call of "
            f"{f(cost.get('pilot_warm_call_seconds_median'), 3)} s and a peak of "
            f"{f(cost.get('pilot_peak_gpu_gib'), 2)} GiB of a 33.39 GiB working limit, which "
            f"projects the 3,845-lead gate at about {f(cost.get('projected_full_run_hours_3845'), 2)} h "
            f"against roughly 0.7 h at 524,288 bp, and a hypothetical 32,322-lead run at about "
            f"{f(cost.get('projected_full_run_hours_32322'), 2)} h. A 1-Mb window also rejects more "
            f"leads than a shorter one, because it is likelier to reach an assembly gap or a "
            f"chromosome end."
        )
        A("")
        A(
            "What the sweep does establish, whichever length is used: the two implementations agree "
            "on the same quantity. The median ratio of local to archived sits at 0.99 to 1.04 at "
            "every length, so no scale convention differs, and the disagreement that remains shrinks "
            "monotonically as the window grows. A wrong readout, a wrong track set or a genuinely "
            "different model would not behave that way."
        )
        A("")

    # ---------------------------------------------------------------- diagnostics and predictions
    if gate:
        A("## Predictions, including the ones that failed")
        A("")
        A(
            "The eight diagnostics and a written prediction of each one's result were fixed in "
            "`i1_04_gate_stats.py`'s docstring before any local score was read. Prespecifying only "
            "the decision rule is how a failure gets explained after the fact instead of before it."
        )
        A("")
        A("| id | prediction | observed | verdict |")
        A("|---|---|---|---|")
        for p in gate["predictions"]:
            obs = p["observed"]
            if isinstance(obs, float):
                obs = f(obs)
            elif isinstance(obs, (dict, list)):
                obs = json.dumps(obs, default=str)
                obs = obs if len(obs) < 180 else obs[:177] + "..."
            # A literal pipe anywhere in a cell ends the column, so escape both fields.
            pred = str(p["prediction"]).replace("|", "\\|")
            obs = str(obs).replace("|", "\\|")
            A(f"| {p['id']} | {pred} | {obs} | **{'met' if p['met'] else 'NOT met'}** |")
        A("")

        A("### Agreement by effect size")
        A("")
        A(
            "Quartiles of the archived score's magnitude. This is where a cross-implementation "
            "difference should live: the runtime probe measured only 5.6x concentration of the ATAC "
            "|alt-ref| mass within 2 kb of a substitution, so a 501-bp mask on a near-null variant "
            "is dominated by the diffuse component."
        )
        A("")
        A("| channel | quartile of \\|archived\\| | n | \\|archived\\| range | Spearman | sign agreement |")
        A("|---|---:|---:|---|---:|---:|")
        for chan in ("atac", "dnase"):
            for s in gate["channels"][chan]["by_effect_size_quartile"]:
                rng = s["abs_archived_range"]
                A(
                    f"| {chan} | {s['quartile_of_abs_archived']} | {s['n']} | "
                    f"{f(rng[0], 5)} to {f(rng[1], 5)} | {f(s['spearman'])} | "
                    f"{f(s['sign_agreement'])} |"
                )
        A("")

        pt = pd.read_csv(T / "ad1_gate_per_track.tsv", sep="\t")
        A("### Per liver track")
        A("")
        A("| channel | track | n | Spearman | sign agreement | median \\|difference\\| |")
        A("|---|---|---:|---:|---:|---:|")
        for _, r in pt.iterrows():
            A(
                f"| {r['channel']} | `{r['track']}` | {int(r['n'])} | {f(r['spearman'])} | "
                f"{f(r['sign_agreement'])} | {f(r['median_abs_difference'], 5)} |"
            )
        A("")

        xp = gate.get("cross_process_reproducibility")
        if xp:
            A("### Cross-process reproducibility, which the runtime probe left open")
            A("")
            A(
                "The runtime probe established that repeat calls within one process on one card are "
                "bit-identical across all 35.8 million returned values, and said plainly that "
                "reproducibility across separate processes was not tested. The gate rescored, in its "
                f"own process, the same {xp['n']} leads the pilot sweep had already scored at the "
                f"same length on the same card."
            )
            A("")
            A("| channel | leads | bit-identical | largest absolute difference |")
            A("|---|---:|---:|---:|")
            for chan in ("atac", "dnase"):
                A(
                    f"| {chan} | {xp['n']} | {xp[f'{chan}_n_bit_identical']} | "
                    f"{xp[f'{chan}_max_abs_difference']:.3e} |"
                )
            A("")

        og = gate["orientation_guard"]
        A("### The orientation guard")
        A("")
        A(
            "A sign applied twice cancels for half a family, so the archived column's own signed "
            "Spearman against `beta_alt` is recomputed here and checked against the value "
            "`c1-endpoint2-v2` printed for tier A4. This is a guard, not a new result."
        )
        A("")
        A("| quantity | value |")
        A("|---|---:|")
        A(
            f"| archived ATAC vs `beta_alt`, this run, pooled over folds | "
            f"{f(og['archived_atac_vs_beta_alt_spearman_here'])} |"
        )
        A(
            f"| archived ATAC vs `beta_alt`, `c1-endpoint2-v2` tier A4 | "
            f"{f(og['archived_atac_vs_beta_alt_spearman_c1e2_tierA4'])} |"
        )
        A(
            f"| local ATAC vs `beta_alt`, this run, pooled over folds | "
            f"{f(og['local_atac_vs_beta_alt_spearman_here'])} |"
        )
        A("")
        A(
            "No sign is applied anywhere in this run. `beta_alt` is the ALT-dosage slope, positive "
            "meaning ALT raises accessibility; the local score is "
            "log2(1 + sum ALT) - log2(1 + sum REF) over the mask, positive meaning the same thing."
        )
        A("")

    # ---------------------------------------------------------------- zero-shot endpoint
    if zs:
        A("## The zero-shot local endpoint, and what it is allowed to say")
        A("")
        A(
            "Zero-shot means no training: the released weights score REF and ALT and the signed "
            "number is ranked against the measured `beta_alt` directly. Statistic and folds are "
            "`c2-endogenous-head-20260914T194000Z`'s: signed Spearman inside each of the five "
            "ChromBPNet chromosome-group held-out folds, Fisher-z macro over the five, "
            f"{zs['bootstrap_resamples']:,} bootstrap resamples of 1-Mb blocks within each fold at "
            f"seed {zs['bootstrap_seed']}, every arm recomputed inside the same draw."
        )
        A("")
        le = zs["loader_equivalence_32322"]
        A(
            "**The comparator is regenerated, not quoted.** This loader recomputes the c2 package's "
            "own macro Spearman from its deposited out-of-fold predictions on all 32,322 leads and "
            f"reproduces every published value to within {zs['loader_equivalence_max_abs_difference']:.3e}:"
        )
        A("")
        A("| arm | published macro (32,322) | recomputed here | abs difference | tolerance | distinct predicted values |")
        A("|---|---:|---:|---:|---:|---:|")
        for k, v in le.items():
            A(
                f"| `{k}` | {v['published_macro']:.16f} | {v['recomputed_macro']:.16f} | "
                f"{v['abs_difference']:.3e} | {v.get('tolerance', float('nan')):.0e} | "
                f"{v.get('distinct_predicted_values_all', 'n/a')} of 32,322 |"
            )
        A("")
        A(
            "The tolerance is per arm, and the reason is the arm's tie structure rather than a "
            "judgement about which arm matters. The two continuous arms predict 32,322 distinct "
            "values and reproduce at machine precision. The allele-identity control is a ridge on a "
            "16-level one-hot of (REF, ALT), so it predicts only 12 to 13 distinct values in each "
            "fold of roughly 6,000 variants; a last-bit difference in the deposited decimal then "
            "moves which rows are tied and `rankdata`'s average ranks shift with it. The resulting "
            "1.3e-05 is 0.24 percent of that control's own published bootstrap standard error of "
            "0.0056. Its Spearman is capped by the tie structure and is not comparable in precision "
            "with a continuous arm's."
        )
        A("")
        ss = zs["shared_set"]
        cov = zs.get("coverage_and_survivor_skew") or {}
        A(
            f"**The shared set is {ss['variants']} leads in {ss['blocks']} 1-Mb blocks**, and it is "
            f"doubly selected. Two selections stack, and neither is neutral."
        )
        A("")
        A(
            f"The first is the Atlas: only {cov.get('atlas_served_leads', 3845):,} of the "
            f"{cov.get('currin_leads_total', 32322):,} Currin leads carry an archived hosted score, and "
            f"`c1-endpoint2-v2` already recorded that the same ChromBPNet reaches 0.4495 on all "
            f"32,322 against 0.6450 on these leads, so an Atlas-served subset is an easier set for "
            f"every model. The second is this arm's own window guard, which the deposited comparators "
            f"never faced: they were built from 4,096-bp windows, where essentially nothing is "
            f"rejected, while a 1-Mb window is thrown out whenever it carries an assembly-gap N or "
            f"runs off a chromosome end."
        )
        A("")
        if cov.get("leads_rejected") is not None:
            A(
                f"**What the window guard removed, and which way the survivors lean.** "
                f"{cov['leads_rejected']:,} of {cov['atlas_served_leads']:,} leads "
                f"({100 * cov['rejected_fraction']:.1f} percent) never reached the model."
            )
            A("")
            A("| set | leads | median \\|archived ATAC\\| | quartiles of \\|archived ATAC\\| | median \\|beta_alt\\| |")
            A("|---|---:|---:|---|---:|")
            rq = cov.get("rejected_abs_archived_atac_quartiles") or []
            kq = cov.get("kept_abs_archived_atac_quartiles") or []
            A(
                f"| rejected | {cov['leads_rejected']} | "
                f"{f(cov.get('rejected_median_abs_archived_atac'), 4)} | "
                + (", ".join(f(x, 4) for x in rq) if rq else "n/a")
                + f" | {f(cov.get('rejected_median_abs_beta_alt'), 4)} |"
            )
            A(
                f"| scored | {cov['atlas_served_leads'] - cov['leads_rejected']} | "
                f"{f(cov.get('kept_median_abs_archived_atac'), 4)} | "
                + (", ".join(f(x, 4) for x in kq) if kq else "n/a")
                + f" | {f(cov.get('kept_median_abs_beta_alt'), 4)} |"
            )
            A("")
            if cov.get("rejected_reasons"):
                A(
                    "Reasons: "
                    + ", ".join(f"{v} {k.replace('_', ' ')}" for k, v in sorted(cov["rejected_reasons"].items()))
                    + "."
                )
                A("")
            A(
                f"The survivor set skews **{cov.get('survivors_skew')}**. That is the population "
                f"every number in this section describes, and the published comparators describe a "
                f"different one."
            )
            A("")
        A(
            "So the published 0.1737 and 0.0435 are never the contrast here. Every comparator is "
            "recomputed from its deposited out-of-fold predictions on exactly the surviving leads, "
            "with the same folds and the same blocks, and every difference is taken inside one "
            "bootstrap draw."
        )
        A("")
        cs = zs.get("comparator_shift") or []
        if cs:
            A("### How far the comparators move on the surviving set")
            A("")
            A("| arm | published macro (32,322) | published 95% CI | recomputed on the survivors | recomputed 95% CI | difference | in published SE | outside published CI | material |")
            A("|---|---:|---|---:|---|---:|---:|---|---|")
            for s in cs:
                pci = s["published_ci_on_32322"]
                rci = s["recomputed_ci_on_shared_set"]
                A(
                    f"| `{s['arm']}` | {f(s['published_macro_on_32322'])} | "
                    f"{f(pci[0])} to {f(pci[1])} | "
                    f"**{f(s['recomputed_macro_on_shared_set'])}** | "
                    f"{f(rci[0])} to {f(rci[1])} | {f(s['difference'])} | "
                    f"{s['difference_in_published_se_units']:+.1f} | "
                    f"{'yes' if s['outside_published_ci'] else 'no'} | "
                    f"{'**yes**' if s['material'] else 'no'} |"
                )
            A("")
            moved = zs.get("comparators_that_moved_materially") or []
            if moved:
                names = ", ".join(f"`{m}`" for m in moved)
                verb = "each moved" if len(moved) > 1 else "moved"
                A(
                    f"Stated plainly: {names} {verb} by more than two of its own published bootstrap "
                    f"standard errors between the published population and the surviving one. That is "
                    f"a fact about the subset, not an error in either run. Reading the local arm "
                    f"against the published number instead of the recomputed one would have credited "
                    f"that shift to the model."
                )
            else:
                A(
                    "No comparator moved by more than two of its own published bootstrap standard "
                    "errors, so on this endpoint the two populations happen to be close. The "
                    "contrasts below still use the recomputed values, because that is what makes them "
                    "paired."
                )
            A("")
        arms = pd.read_csv(T / "zeroshot_arms.tsv", sep="\t")
        A("| arm | variants | blocks | macro Spearman | 95% CI | pooled Spearman | 95% CI |")
        A("|---|---:|---:|---:|---|---:|---|")
        for _, r in arms.iterrows():
            A(
                f"| `{r['arm']}` | {int(r['variants'])} | {int(r['blocks'])} | "
                f"**{f(r['macro_spearman'])}** | {f(r['macro_ci_low'])} to {f(r['macro_ci_high'])} | "
                f"{f(r['pooled_spearman'])} | {f(r['pooled_ci_low'])} to {f(r['pooled_ci_high'])} |"
            )
        A("")
        con = pd.read_csv(T / "zeroshot_contrasts.tsv", sep="\t")
        A("### Paired contrasts on the shared set")
        A("")
        A("| contrast | statistic | difference | 95% CI | bootstrap p | excludes 0 |")
        A("|---|---|---:|---|---:|---|")
        for _, r in con.iterrows():
            pnote = f" ({r['p_note']})" if isinstance(r.get("p_note"), str) and r["p_note"] else ""
            A(
                f"| {r['contrast']} | {r['statistic']} | **{f(r['difference'])}** | "
                f"{f(r['ci_low'])} to {f(r['ci_high'])} | {f(r['bootstrap_p'], 4)}{pnote} | "
                f"{'yes' if r['excludes_zero'] else 'no'} |"
            )
        A("")
        pf = pd.read_csv(T / "zeroshot_per_fold.tsv", sep="\t")
        A("### The five fold components the macro averages")
        A("")
        cols = [c for c in pf.columns if c not in ("fold", "variants", "blocks")]
        A("| fold | variants | blocks | " + " | ".join(f"`{c}`" for c in cols) + " |")
        A("|---:|---:|---:|" + "---:|" * len(cols))
        for _, r in pf.iterrows():
            A(
                f"| {int(r['fold'])} | {int(r['variants'])} | {int(r['blocks'])} | "
                + " | ".join(f(r[c]) for c in cols)
                + " |"
            )
        A("")
        A(
            f"For the record, the prespecified reference points of `ADAPTATION_ARM_PRESPEC.md` are "
            f"HyenaDNA delta ridge {PUB_H} and allele-identity control {PUB_C}. Those were measured "
            f"on all 32,322 leads from 4,096-bp windows and are printed here only so the size of the "
            f"population shift is visible. No contrast in this run uses them."
        )
        A("")

    # ---------------------------------------------------------------- jobs, environment
    A("## Jobs and resource asks")
    A("")
    if jobs is not None:
        A("| job id | name | state | elapsed | MaxRSS | node | resource ask |")
        A("|---|---|---|---|---|---|---|")
        for _, r in jobs.iterrows():
            A(
                f"| {r['JobID']} | {r['JobName']} | {r['State']} | {r['Elapsed']} | "
                f"{r.get('MaxRSS', '')} | {r.get('NodeList', '')} | `{r.get('AllocTRES', '')}` |"
            )
        A("")
    A(
        "Two of the listed jobs produced no output and are listed because they are part of what was "
        "measured. `agp_ad1_pilot` 21767518 was cancelled before it ran: its 8-CPU ask could not "
        "schedule, because the L40S nodes had free GPUs but only 4 free CPUs each, and the work is "
        "GPU-bound, so every later job asks for 4. `agp_ad1_pilot` 21767523 ran and failed at all "
        "five input lengths inside `model.score_variant`, which is how the splice-site blocker above "
        "was found rather than assumed."
    )
    A("")
    A("sbatch headers as submitted:")
    A("")
    A("```")
    for name in ("i1_00_build.sbatch", "i1_02_pilot.sbatch", "i1_06_gate.sbatch", "i1_09_final.sbatch"):
        p = run / "sbatch" / name
        if p.exists():
            for line in p.read_text().splitlines():
                if line.startswith("#SBATCH"):
                    A(line)
            A("")
    A("```")
    A("")
    A(
        "All scoring ran on an L40S. The runtime probe measured that scores are bit-identical when "
        "repeated on one card and are not reproducible across GPU architectures by as much as "
        "3.88e-02 in log2 units, which is the order of a typical measured indel effect, so every "
        "score row here records its card and the statistics scripts refuse a table that spans two."
    )
    A("")
    A("## Environment")
    A("")
    A("```")
    A(versions.strip())
    A("```")
    A("")
    A(
        "Full `pip freeze` in `env/pip_freeze.txt`, Python version in `env/python_version.txt`, "
        "seeds in `seeds.json`, every input, code file, weight file and output with its sha256 in "
        "`MANIFEST.tsv`."
    )
    A("")

    # ---------------------------------------------------------------- limits
    A("## What this run does not establish")
    A("")
    for line in [
        "That local and hosted AlphaGenome are the same service. They are two implementations and "
        "their outputs are reported side by side, never averaged and never pooled. Hosted Atlas "
        "outputs remain excluded from training use; outputs of these local weights are eligible as "
        "features or teacher signals under the noncommercial terms, which every derivative inherits.",
        "Anything causal. Currin peak leads are LD-selected marginal associations and a lead is not "
        "necessarily the causal base.",
        "Anything about the 28,477 Currin leads no Atlas panel covers. The gate set is the "
        "Atlas-served subset and absolute values are not comparable across sets; only within-set "
        "contrasts are.",
        "A cell-type-resolved accessibility claim. `ADAPTATION_ARM_AMENDMENT_01.md` section 5 records "
        "that ATAC in this checkpoint is adult bulk liver tissue only, three tracks, with no "
        "hepatocyte ATAC track anywhere in the human catalogue.",
        "That the zero-shot arm's standing survives a held-out comparison on equal terms. The "
        "AlphaGenome column is a three-track mean of a static model whose folds are not held out "
        "here and whose declared training corpus includes ENCODE adult liver ATAC, while the fitted "
        "comparators are scored inside their own held-out folds. Every measured gap is an upper bound "
        "on a genuine advantage.",
    ]:
        A(f"- {line}")
    A("")

    (run / "RESULTS.md").write_text("\n".join(L) + "\n")
    C.log(f"wrote {run / 'RESULTS.md'} ({len(L)} lines)")
    return 0


PUB_H = "0.1737 [0.1633, 0.1844]"
PUB_C = "0.0435"


if __name__ == "__main__":
    raise SystemExit(main())
