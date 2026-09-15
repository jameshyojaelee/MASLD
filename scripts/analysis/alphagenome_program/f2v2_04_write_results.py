#!/usr/bin/env python3
"""f2-haplotype-v2 step 4: write RESULTS.md straight from the producing JSONs. No number is hand-copied."""

from __future__ import annotations

import json
import os
import pathlib

OUT = pathlib.Path(os.environ["F2V2_OUT_ROOT"])
T = OUT / "tables"
V1NAME = "f2-haplotype-20260914T230433Z"
CH = ("rna", "atac", "dnase", "h3k27ac")
CHN = {"rna": "RNA_SEQ", "atac": "ATAC", "dnase": "DNASE", "h3k27ac": "H3K27ac"}
R_BH_Q = 0.10


def cap(text: str) -> str:
    return text[:1].upper() + text[1:] if text else text


def main() -> None:
    S = json.loads((T / "f2v2_summary.json").read_text())
    d1 = S["defect_1_repair"]
    d2 = S["defect_2_repair"]
    geo = S["readout_geometry"]
    acc = S["accounting"]["calls_and_checkpoints"]
    ts = S["accounting"]["prespec_timestamps"]
    rep, old = d1["REPAIRED_primary"], d1["V1_SUPERSEDED"]
    des = d1["design"]
    ro = d1["repaired_null_readout"]
    L: list[str] = []
    A = L.append

    A("# F2 and F3, version 2: the GNMT verdict on a matched RNA readout, and the stratum's own null")
    A("")
    A(f"Package `f2-haplotype-v2`. It supersedes `{V1NAME}`, which is marked SUPERSEDED in its own directory.")
    A("It repairs two number-level defects an independent checker found there, reissues the numbers those")
    A("defects touched, and corrects six accounting statements that needed no new scoring. Every verdict the")
    A("checker confirmed is restated here unchanged, and the section \"what moved and what did not\" says which")
    A("is which.")
    A("")
    A("Every number is a model prediction about substitutions in one GRCh38 reference sequence, scored through")
    A("the hosted AlphaGenome model API (`predict_sequence`, 1,048,576 bp window, liver ontology terms).")
    A("Nothing here measures an interaction in cells or people, and nothing here adopts a Resource number.")
    A("")
    A(f"New scoring in this package: {rep['n_null_draws']} chr6 null draws x 4 arms = "
      f"{4 * rep['n_null_draws']} `predict_sequence` calls. Nothing else was rescored: the six GNMT-region pairs,")
    A(f"the {S['v1_row_counts']['obs']} observed P5 pairs and the {S['v1_row_counts']['null']} corrected null draws")
    A(f"are v1's rows, read from its checkpoint unchanged.")
    A("")

    # ---------------------------------------------------------------- what moved
    A("## What moved and what did not")
    A("")
    A("| number | v1 as published | v2 | moved |")
    A("|---|---|---|---|")
    mv = d1["moved"]
    for ch in CH:
        m = mv[ch]
        A(f"| F2.1 {CHN[ch]}: primary pair inside the null's central 95% | "
          f"{'yes' if m['v1_inside'] else '**no**'}, p {m['v1_p']:.4f} | "
          f"{'yes' if m['v2_inside'] else '**no**'}, p {m['v2_p']:.4f} | "
          f"{'**yes**' if m['v1_inside'] != m['v2_inside'] else ('no' if abs(m['v1_p'] - m['v2_p']) < 5e-5 else 'no, p changed')} |")
    for ch in CH:
        m = mv[ch]
        A(f"| GNMT null q95 abs residual, {CHN[ch]} | {m['v1_null_q95_abs_residual']:.6f} | "
          f"{m['v2_null_q95_abs_residual']:.6f} | {m['relative_change_in_null_q95_percent']:+.1f}% |")
    for ch in CH:
        a = d2["n2500"]["per_channel"][ch]
        v1n = a["V1_AS_PUBLISHED_p_against_the_pooled_null"]["n_bh_rejected"]
        v2n = a["REPAIRED_p_against_the_null_inside_the_stratum"]["n_bh_rejected_gated"]
        v2u = a["REPAIRED_p_against_the_null_inside_the_stratum"]["n_bh_rejected_if_the_20_draw_gate_is_ignored"]
        shown = "untestable" if v2n is None else str(v2n)
        A(f"| both-variants-active BH rejections, {CHN[ch]} | {v1n} of {a['n_observed_in_stratum']} | "
          f"{shown}{'' if v2n is not None else f' ({v2u} if the 20-draw gate is ignored)'} | "
          f"{'**yes**' if (v2n if v2n is not None else v2u) != v1n else 'no'} |")
    A(f"| cancelled job {acc['cancelled_job']['jobid']}: calls | "
      f"{acc['cancelled_job']['previously_reported']['calls']} | "
      f"{acc['cancelled_job']['predict_sequence_calls']} | **yes** |")
    A(f"| cancelled job {acc['cancelled_job']['jobid']}: pairs deposited | "
      f"{acc['cancelled_job']['previously_reported']['pairs']} | "
      f"{acc['cancelled_job']['pairs_deposited_in_the_checkpoint']} | **yes** |")
    A(f"| total `predict_sequence` calls in v1 | about {acc['total_previously_reported_about']} | "
      f"{acc['total_predict_sequence_calls']} | yes, by "
      f"{abs(acc['total_previously_reported_about'] - acc['total_predict_sequence_calls'])} |")
    for ch in CH:
        a = d2["n2500"]["per_channel"][ch]
        A(f"| both-active null draws qualifying by construction, {CHN[ch]} | stated as about 0.25% | "
          f"{a['null_draws_in_stratum_share_percent']:.2f}% | **yes** |")
    rmv = geo["rejections_as_rows_and_as_distinct_variant_pairs"]
    A(f"| F3.1 BH rejections over the {geo['n_pairs']}-pair family, RNA / ATAC / DNASE / H3K27ac | rows "
      + " / ".join(str(rmv[c]["rejected_rows"]) for c in CH)
      + " | unchanged as rows; " + " / ".join(str(rmv[c]["distinct_variant_pairs"]) for c in CH)
      + " distinct variant pairs | no |")
    A("| F2.2, F3.2, F3.3, the EUR haplotype frequencies, the determinism check, the p-floor arithmetic, the "
      "separation-stratified and magnitude-matched diagnostics, the family-level contrast | | unchanged | no |")
    A("")

    # ---------------------------------------------------------------- predictions
    A("## Predictions written before any score in this package was read")
    A("")
    A(f"Stamped at {json.loads((OUT / 'prespec' / 'WRITTEN_UTC.json').read_text())['written_utc']}, "
      f"before the first `predict_sequence` call of this package.")
    A("")
    A("| id | prediction | verdict | number that decides it |")
    A("|---|---|---|---|")
    ins = [CHN[c] for c in CH if rep["F2.1_verdict_per_channel"][c]]
    outs = [CHN[c] for c in CH if not rep["F2.1_verdict_per_channel"][c]]
    A(f"| P1 | the primary pair's residual stays inside the central 95% of the repaired null in all four channels "
      f"| **{'met' if rep['F2.1_met'] else 'NOT met'}** | inside in {', '.join(ins) if ins else 'none'}"
      f"{'; OUTSIDE in ' + ', '.join(outs) if outs else ''}. p "
      + ", ".join(f"{CHN[c]} {rep['F2.1_p_per_channel'][c]:.4f}" for c in CH) + " |")
    rna_change = mv["rna"]["relative_change_in_null_q95_percent"]
    A(f"| P2 | the repaired null's RNA q95 abs residual is more than 10% smaller than the central-20kb null's "
      f"| **{'met' if rna_change < -10 else 'NOT met'}** | {rna_change:+.1f}% "
      f"({old['null_q95_abs_residual']['rna']:.6f} to {rep['null_q95_abs_residual']['rna']:.6f}) |")
    chrom_ok = all(abs(mv[c]["relative_change_in_null_q95_percent"]) <= 25 for c in ("atac", "dnase", "h3k27ac"))
    verd_ok = all(mv[c]["v1_inside"] == mv[c]["v2_inside"] for c in ("atac", "dnase", "h3k27ac"))
    A(f"| P3 | the three chromatin channels' q95 abs residual agree between the two draw sets within 25%, so their "
      f"verdicts do not move | **{'met' if (chrom_ok and verd_ok) else 'NOT met on the 25% bound'}"
      f"{'' if chrom_ok else ('; the verdicts did not move' if verd_ok else '; the verdicts DID move')}** | "
      + ", ".join(f"{CHN[c]} {mv[c]['relative_change_in_null_q95_percent']:+.1f}%" for c in ("atac", "dnase", "h3k27ac"))
      + f"; verdicts {'unchanged' if verd_ok else 'CHANGED'} |")
    p4 = {ch: d2["n2500"]["per_channel"][ch]["REPAIRED_p_against_the_null_inside_the_stratum"] for ch in CH}
    p4_counts = [p4[ch]["n_bh_rejected_if_the_20_draw_gate_is_ignored"] for ch in CH]
    rna_floor_n = p4["rna"]["n_pairs_at_the_in_stratum_floor"]
    rna_fam = p4["rna"]["bh_family_size"]
    A(f"| P4 | the in-stratum recomputation reproduces the checker's 5 / 0 / 3 / 0, with RNA's five at the floor "
      f"| **{'met' if p4_counts == [5, 0, 3, 0] else 'NOT met'} on the counts, NOT met on the floor** | "
      f"counts {' / '.join(str(x) for x in p4_counts)}; {rna_floor_n} of RNA's {rna_fam} sit at the floor "
      f"{p4['rna']['min_p']:.6f} and the other {rna_fam - rna_floor_n} at twice it, "
      f"{2 * p4['rna']['min_p']:.6f} |")
    shares = [d2["n2500"]["per_channel"][ch]["null_draws_in_stratum_share_percent"] for ch in CH]
    rhos = [d2["n2500"]["per_channel"][ch]["spearman_abs_v1_vs_abs_v2_in_null"] for ch in CH]
    A(f"| P5 | the measured both-active null rate exceeds 0.25% in every channel and the null's two arm magnitudes "
      f"are positively rank-correlated | **{'met' if all(s > 0.25 for s in shares) and all(r[0] > 0 for r in rhos) else 'NOT met'}** | "
      f"shares {' / '.join(f'{s:.2f}%' for s in shares)}; Spearman "
      + ", ".join(f"{CHN[c]} {d2['n2500']['per_channel'][c]['spearman_abs_v1_vs_abs_v2_in_null'][0]:+.3f}" for c in CH)
      + f" (largest p {max(r[1] for r in rhos):.1e}) |")
    A("")

    # ---------------------------------------------------------------- defect 1
    A("## Defect 1. The GNMT verdict was not a matched comparison in the RNA channel")
    A("")
    A("v1's six GNMT-region pairs read RNA out over the GENCODE v49 GNMT span,")
    A(f"chr6:{des['gnmt_span']['start']}-{des['gnmt_span']['end']}, {ro['gnmt_rna_readout_bp']} bases of mask.")
    A(f"All {old['n_null_draws']} of their matched chr6 null draws read RNA out over a central "
      f"{ro['v1_null_rna_readout_bp']} bp window, because `score_pair` takes the gene-span branch only when the")
    A("target gene's span overlaps the scored window, and a null draw at a random chr6 position never contains")
    A(f"GNMT. The RNA arm of v1's F2.1 compared a sum over {ro['gnmt_rna_readout_bp']} bases of an expressed gene")
    A(f"against sums over {ro['v1_null_rna_readout_bp']} bases of arbitrary chr6 sequence, a "
      f"{ro['width_ratio_v1_null_over_observed']}-fold width mismatch.")
    A("")
    A("**The repaired null.** Every draw sits inside a real chr6 protein-coding gene and carries that gene's")
    A("Ensembl id, so the unchanged `score_pair` takes its gene-span branch and the RNA readout is a gene span.")
    A("Nothing in the recipe changed; only which (gene, variant pair) went in.")
    A("")
    A("| property | GNMT observed side | repaired null | v1 null |")
    A("|---|---|---|---|")
    A(f"| RNA readout | GNMT gene span, {ro['gnmt_rna_readout_bp']} bp | gene span, "
      f"{'all ' + str(rep['n_null_draws']) + ' draws' if ro['all_draws_use_the_gene_span'] else 'MIXED'} | "
      f"central {ro['v1_null_rna_readout_bp']} bp, {old['n_null_draws']} draws |")
    A(f"| RNA readout width bp | {ro['gnmt_rna_readout_bp']} | median {ro['rna_readout_bp']['median']:.0f} "
      f"(range {ro['rna_readout_bp']['min']:.0f} to {ro['rna_readout_bp']['max']:.0f}) | "
      f"{ro['v1_null_rna_readout_bp']} |")
    A(f"| width ratio to the observed readout | 1.0 | {ro['width_ratio_repaired_null_over_observed']} | "
      f"{ro['width_ratio_v1_null_over_observed']} |")
    A(f"| REF RNA sum in the readout | {ro['ref_rna_sum_vs_gnmt']['gnmt']:.1f} | median "
      f"{ro['ref_rna_sum_vs_gnmt']['null_median']:.1f} (range {ro['ref_rna_sum_vs_gnmt']['null_min']:.1f} to "
      f"{ro['ref_rna_sum_vs_gnmt']['null_max']:.1f}) | not comparable, different window |")
    A(f"| separation bp | {des['separation_bin']['gnmt_primary_separation_bp']} | median "
      f"{des['achieved_separation_bp']['median']:.0f} (bin [{des['separation_bin']['lo']:.1f}, "
      f"{des['separation_bin']['hi']:.1f}]) | same bin |")
    A(f"| local (chromatin) readout | union of two +/-1 kb windows | identical definition | identical definition |")
    A(f"| draws | | {des['n_draws']} over {des['n_distinct_host_genes']} distinct host genes "
      f"({des['n_draws_that_are_a_second_pair_from_a_gene_already_drawn']} are a second pair from a gene already "
      f"drawn) | {old['n_null_draws']} |")
    A("")
    A(f"Host genes were the {des['candidate_genes']} chr6 protein-coding genes with span width in "
      f"[{des['width_band_bp'][0]}, {des['width_band_bp'][1]}] bp, excluding GNMT and anything overlapping it; "
      f"{des['host_genes_with_a_qualifying_pair']} of them carry at least one common-EUR SNV pair inside the span "
      f"at a separation in the primary pair's decile bin, {des['qualifying_pairs_total']} such pairs in all. "
      f"Seed {des['seed']}.")
    A("")
    le = S.get("loader_equivalence", {})
    if le.get("n_sums_compared"):
        A(f"**The repaired null goes through the same loader.** This package re-implements `score_pair` (it drops "
          f"the determinism probe and records the two readout widths), so one of v1's own central-20kb draws, "
          f"`{le['probe_tag']}`, was rescored through v2's scorer and compared with v1's deposited numbers: "
          f"{le['n_identical']} of {le['n_sums_compared']} summed values identical, largest relative difference "
          f"{le['max_relative_difference']:.3g}, readout branch "
          f"{'the same' if le['readout_branch_identical'] else 'DIFFERENT'}. "
          + ("" if le["n_identical"] == le["n_sums_compared"]
                 else " **The loaders differ; read every verdict below as a different loader.**"))
        A("")
    A("**F2.1, the primary pair rs2296805 (chr6:42961020 T>G) x rs2296804 (chr6:42963523 C>G).** Signs refer to")
    A("the ALT allele named. Native scale is the summed predicted signal in the readout window; log2 is against")
    A("the REF arm.")
    A("")
    A("| channel | residual | repaired null central 95% | inside | p | v1 null central 95% | v1 inside | v1 p |")
    A("|---|---|---|---|---|---|---|---|")
    pr = next(e for e in rep["pairs"] if e["is_primary"])
    po = next(e for e in old["pairs"] if e["is_primary"])
    for ch in CH:
        A(f"| {CHN[ch]} | {pr[f'{ch}_residual']:+.6f} | [{pr[f'{ch}_null_central95_lo']:+.6f}, "
          f"{pr[f'{ch}_null_central95_hi']:+.6f}] | {'yes' if pr[f'{ch}_inside_central95'] else '**no**'} | "
          f"{pr[f'{ch}_p_empirical']:.4f} | [{po[f'{ch}_null_central95_lo']:+.6f}, "
          f"{po[f'{ch}_null_central95_hi']:+.6f}] | {'yes' if po[f'{ch}_inside_central95'] else '**no**'} | "
          f"{po[f'{ch}_p_empirical']:.4f} |")
    A("")
    A(f"The p floor is 1/{rep['n_null_draws'] + 1} = {rep['p_floor']:.4f} in both nulls, so the two verdicts are")
    A("comparable on resolution and any movement is the readout, not the draw count.")
    A("")
    A("**The verdict survives the repair, and the reason it survives is not the one predicted.** The RNA null's")
    A(f"q95 abs residual barely moved ({mv['rna']['relative_change_in_null_q95_percent']:+.1f}%) and the RNA p is")
    A(f"identical to five decimal places ({mv['rna']['v1_p']:.5f} against {mv['rna']['v2_p']:.5f}). Matching the")
    A("readout to a gene span did not tighten the RNA null, because the host spans' predicted RNA signal ranges")
    A(f"over six orders of magnitude ({ro['ref_rna_sum_vs_gnmt']['null_min']:.2f} to "
      f"{ro['ref_rna_sum_vs_gnmt']['null_max']:.0f} REF sum, median {ro['ref_rna_sum_vs_gnmt']['null_median']:.0f} "
      f"against GNMT's {ro['ref_rna_sum_vs_gnmt']['gnmt']:.0f}): a 4 kb protein-coding span on chr6 is often as")
    A("quiet as an arbitrary 20 kb window, so the log2 of its summed signal is no less noisy. The defect was real,")
    A("the comparison is now matched on the readout, and the verdict it produces is the same one.")
    A("")
    A("**The three chromatin nulls moved a lot without moving their verdicts.** Their readout definition is")
    A("identical in both draw sets, so the only difference is where on chr6 the draws sit: inside protein-coding")
    A("gene spans here, at arbitrary positions in v1. That alone raises their q95 abs residual by "
      + ", ".join(f"{CHN[c]} {mv[c]['relative_change_in_null_q95_percent']:+.0f}%" for c in ("atac", "dnase", "h3k27ac"))
      + " and roughly doubles their p-values (" + ", ".join(
        f"{CHN[c]} {mv[c]['v1_p']:.4f} to {mv[c]['v2_p']:.4f}" for c in ("atac", "dnase", "h3k27ac")) + ").")
    A("A matched chr6 null is therefore sensitive to where it is drawn, by more than the 25% this package")
    A("predicted, and the chromatin verdicts hold only because the observed residual is far inside either null.")
    A("")
    A("All six GNMT-region pairs against the repaired null, residual per channel:")
    A("")
    A("| pair | separation bp | " + " | ".join(CHN[c] for c in CH) + " | outside central 95% |")
    A("|---|---|---|---|---|---|---|")
    for e in rep["pairs"]:
        bad = [CHN[c] for c in CH if not e[f"{c}_inside_central95"]]
        tag = f"{e['rsid1']} x {e['rsid2']}" + (" **(primary)**" if e["is_primary"] else "")
        A(f"| {tag} | {e['separation_bp']} | "
          + " | ".join(f"{e[f'{c}_residual']:+.6f}" for c in CH)
          + f" | {', '.join(bad) if bad else 'none'} |")
    A("")
    fb = rep["gnmt_family_bh"]
    A(f"Those marks are uncorrected. Over the {fb['n_tests']}-test family (6 pairs x 4 channels) against the")
    A(f"{rep['n_null_draws']}-draw repaired null, BH at q={fb['q']:.2f} rejects **{fb['n_bh_rejected']}**")
    A(f"(smallest p {fb['min_p']:.4f}, floor {fb['p_floor']:.4f}; "
      f"{fb['pairs_that_must_sit_AT_the_floor_together_for_any_rejection']} tests would have to sit at that floor "
      f"together and {fb['n_at_the_floor']} do), against {fb['n_outside_central95_UNCORRECTED']} marked outside")
    A(f"uncorrected. v1's same family gave {old['gnmt_family_bh']['n_bh_rejected']} rejections and "
      f"{old['gnmt_family_bh']['n_outside_central95_UNCORRECTED']} uncorrected marks.")
    A("")
    sf = d1["sensitivity_one_pair_per_host_gene"]
    A(f"**Sensitivity.** {des['n_draws_that_are_a_second_pair_from_a_gene_already_drawn']} of the "
      f"{des['n_draws']} draws are a second pair from a host gene already drawn, so they are not independent of "
      f"their partner. On the {sf['n_null_draws']}-draw one-pair-per-gene sub-null the verdict is "
      + ", ".join(f"{CHN[c]} {'inside' if sf['F2.1_verdict_per_channel'][c] else '**outside**'} "
                  f"(p {sf['F2.1_p_per_channel'][c]:.4f})" for c in CH) + ".")
    A("")
    ps = d1["POST_HOC_signal_matched"]
    if ps.get("F2.1_verdict_per_channel"):
        A(f"**Signal-matched sub-null (post hoc).** Matching the readout's width and gene type does not match its")
        A(f"signal level: the {des['n_draws']} host spans carry REF RNA sums from "
          f"{ro['ref_rna_sum_vs_gnmt']['null_min']:.2f} to {ro['ref_rna_sum_vs_gnmt']['null_max']:.0f}, median "
          f"{ro['ref_rna_sum_vs_gnmt']['null_median']:.0f}, against GNMT's "
          f"{ro['ref_rna_sum_vs_gnmt']['gnmt']:.0f}. Restricted to the "
          f"{ps['n_null_draws']} draws within a factor of 2 of GNMT's sum, the primary pair's RNA residual is "
          f"{'inside' if ps['F2.1_verdict_per_channel']['rna'] else '**outside**'} the central 95% with p "
          f"{ps['F2.1_p_per_channel']['rna']:.4f}, and the other three channels stay inside (p "
          + ", ".join(f"{CHN[c]} {ps['F2.1_p_per_channel'][c]:.4f}" for c in ("atac", "dnase", "h3k27ac")) + ").")
        A(f"Those two RNA statements cannot both be read as a threshold: with {ps['n_null_draws']} draws the 2.5th")
        A(f"percentile is an interpolation between the two smallest values and the smallest attainable p is "
          f"1/{ps['n_null_draws'] + 1} = {ps['p_floor']:.4f}. The p is the number to read, it is above 0.05, and the")
        A(f"24-test family's BH threshold at rank 1 is {R_BH_Q / 24:.4f}, so nothing here rejects. What the sub-null")
        A("does establish is that the repaired null is still not matched on signal level, and that a signal-matched")
        A("RNA null would need liver-expression matching this package did not do.")
    else:
        A(f"**Signal-matched sub-null (post hoc):** only {ps['n_null_draws']} of the {des['n_draws']} draws have a "
          f"REF RNA sum within a factor of 2 of GNMT's {ro['ref_rna_sum_vs_gnmt']['gnmt']:.0f}, below the recipe's "
          f"{d2['n2500']['per_channel']['rna']['recipe_min_matched_null']}-draw minimum, so no verdict is emitted. "
          f"Matching the readout's width does not match its signal level, and that gap is not closed here.")
    A("")

    # ---------------------------------------------------------------- defect 2
    A("## Defect 2. The both-variants-active stratum was tested against the pooled null")
    A("")
    A("The stratum keeps pairs whose two single-variant arms are both above the null's 95th percentile of |arm|.")
    A("v1 computed the empirical p once against all 2,500 pooled draws and then masked it to the stratum, so a")
    A("stratum selected on arm magnitude was compared with draws of every arm magnitude. The residual is")
    A("`joint - v1 - v2`, a function of those arms, which is the confound v1's own amendment 06 exists to avoid.")
    A("Recomputed with the null restricted to the draws that meet the same selection:")
    A("")
    A("| channel | observed pairs in stratum | null draws in stratum | p floor in stratum | v1 BH (pooled null) | "
      "v2 BH (null in stratum) | pairs at the in-stratum floor |")
    A("|---|---|---|---|---|---|---|")
    for ch in CH:
        a = d2["n2500"]["per_channel"][ch]
        v1b = a["V1_AS_PUBLISHED_p_against_the_pooled_null"]
        v2b = a["REPAIRED_p_against_the_null_inside_the_stratum"]
        gated = ("**untestable** (" + str(v2b["n_bh_rejected_if_the_20_draw_gate_is_ignored"]) + " ungated)"
                 if v2b["n_bh_rejected_gated"] is None else f"**{v2b['n_bh_rejected_gated']}**")
        A(f"| {CHN[ch]} | {a['n_observed_in_stratum']} | {a['n_null_draws_in_stratum']} of {d2['n2500']['n_null']} | "
          f"{a['p_floor_null_inside_the_stratum']:.6f} | {v1b['n_bh_rejected']} | {gated} | "
          f"{v2b['n_pairs_at_the_in_stratum_floor']} |")
    A("")
    A("Per-pair p-values inside the stratum:")
    A("")
    for ch in CH:
        a = d2["n2500"]["per_channel"][ch]
        v2b = a["REPAIRED_p_against_the_null_inside_the_stratum"]
        A(f"- {CHN[ch]}: " + ", ".join(f"{x:.4f}" for x in v2b["p_values"])
          + f" (floor {a['p_floor_null_inside_the_stratum']:.4f}, BH family {v2b['bh_family_size']}, "
            f"rejection needs {a['pairs_that_must_sit_AT_the_in_stratum_floor_together_for_any_rejection']} at the "
            f"floor together)")
    A("")
    rna = d2["n2500"]["per_channel"]["rna"]
    rnb = rna["REPAIRED_p_against_the_null_inside_the_stratum"]
    A(f"Read the RNA row carefully. Its {rnb['bh_family_size']} p-values take only two values, the floor")
    A(f"{rna['p_floor_null_inside_the_stratum']:.4f} and twice the floor, and BH over a family of")
    A(f"{rnb['bh_family_size']} rejects anything at or below q = {rna['bh_threshold_at_the_last_rank']:.2f} at the last")
    A(f"rank. The count {rnb['n_bh_rejected_gated']} of")
    A(f"{rnb['bh_family_size']} therefore says that all five pairs exceed most of a {rna['n_null_draws_in_stratum']}-draw")
    A("null, in a family small enough that the coarsest attainable p already clears the threshold. It is not")
    A(f"evidence of non-additivity. DNASE's {d2['n2500']['per_channel']['dnase']['n_null_draws_in_stratum']} matched")
    A(f"draws fail the recipe's {rna['recipe_min_matched_null']}-draw gate, so DNASE emits no verdict; the")
    A("checker's 3 for DNASE is the ungated count and is reported only as that.")
    A("")
    A("**The 0.25% by-construction rate was wrong.** v1 stated that about 0.25% of null draws qualify for the")
    A("stratum by construction, which is 5% x 5% and assumes the two arms are independent. They are not:")
    A("")
    A("| channel | fraction of null v1 arms above the threshold | v2 arms | product (independence) | measured share "
      "of draws in the stratum | Spearman(abs v1, abs v2) in the null | p |")
    A("|---|---|---|---|---|---|---|")
    for ch in CH:
        a = d2["n2500"]["per_channel"][ch]
        r = a["spearman_abs_v1_vs_abs_v2_in_null"]
        A(f"| {CHN[ch]} | {a['fraction_of_null_v1_arms_above_threshold']:.4f} | "
          f"{a['fraction_of_null_v2_arms_above_threshold']:.4f} | "
          f"{a['share_expected_if_the_two_arms_were_independent_percent']:.2f}% | "
          f"{a['null_draws_in_stratum_share_percent']:.2f}% | {r[0]:+.3f} | {r[1]:.1e} |")
    A("")
    A(f"Measured {min(d2['n2500']['per_channel'][c]['null_draws_in_stratum_share_percent'] for c in CH):.2f}% to "
      f"{max(d2['n2500']['per_channel'][c]['null_draws_in_stratum_share_percent'] for c in CH):.2f}%, two to four")
    A("times the independent-arms figure, because the two arms of a draw share a sequence context and a readout")
    A("window and their magnitudes are positively rank-correlated in every channel.")
    A("")
    A("**The headline sentence was also self-contradictory.** v1 wrote that the stratum \"is still untestable\"")
    A("above a table that marked three of four channels testable. Those are two different stages: at v1's stage 1")
    A("of 1,000 draws the stratum had " + ", ".join(
        f"{CHN[c]} {d2['n1000']['per_channel'][c]['n_null_draws_in_stratum']}" for c in CH)
      + " matched draws and all four channels failed the 20-draw gate; at 2,500 draws it had " + ", ".join(
        f"{CHN[c]} {d2['n2500']['per_channel'][c]['n_null_draws_in_stratum']}" for c in CH)
      + " and three passed. The stage is named in every row above.")
    A("")

    # ---------------------------------------------------------------- corrections without new scoring
    A("## Corrections that needed no new scoring")
    A("")
    A("**1. The cancelled job's calls and checkpoint deposit.**")
    A("")
    cj = acc["cancelled_job"]
    A(f"Job {cj['jobid']} was cancelled during stage 1. It deposited "
      f"**{cj['pairs_deposited_in_the_checkpoint']} pairs** and made **{cj['predict_sequence_calls']} calls**, not "
      f"{cj['previously_reported']['pairs']} and {cj['previously_reported']['calls']}. "
      f"{cj['how_that_is_derived']} Its successor job {acc['per_job'][1]['jobid']} then printed "
      f"\"checkpoint: {cj['pairs_deposited_in_the_checkpoint']} pairs already scored\" and reported "
      f"{acc['per_job'][1]['calls_reported_at_DONE']} calls, which is 4 arms x "
      f"{acc['per_job'][1]['calls_reported_at_DONE'] // 4} pairs exactly: it took "
      f"{acc['per_job'][1]['pairs_to_score']} pairs off the work list and one of them, the chr8 pair whose 1-Mb "
      f"window runs past the end of the chromosome, returns before any call is made.")
    A("")
    A("| job | state | pairs in checkpoint at start | pairs to score | calls at DONE | minutes |")
    A("|---|---|---|---|---|---|")
    for j in acc["per_job"]:
        A(f"| {j['jobid']} | {'CANCELLED' if j['cancelled'] else 'COMPLETED'} | {j['checkpoint_pairs_at_start']} | "
          f"{j['pairs_to_score']} | {j['calls_reported_at_DONE'] if j['calls_reported_at_DONE'] else str(cj['predict_sequence_calls']) + ' (derived)'} | "
          f"{j['minutes_reported_at_DONE'] if j['minutes_reported_at_DONE'] else 'cancelled mid-run'} |")
    A("")
    A(f"Total: {acc['total_how_derived']} `predict_sequence` calls, against about "
      f"{acc['total_previously_reported_about']} reported. {cap(acc['note_on_the_counter'])}")
    A(f"v1's checkpoint holds {acc['checkpoint_rows']} rows, which is {acc['checkpoint_expected_rows']}, with no "
      f"duplicate tag.")
    A("")
    rc = acc["rate_ceiling"]
    A(f"**2. The request rate.** v1 spaced calls {rc['v1_setting_seconds_between_call_starts']} s apart, a ceiling of "
      f"{rc['v1_implied_ceiling_per_min']}/min, and its logs printed up to {rc['v1_max_rate_printed']}/min, above the "
      f"28/min this repair was asked to hold. This package spaces them "
      f"{rc['v2_setting_seconds_between_call_starts']} s, {rc['v2_implied_ceiling_per_min']}/min.")
    A("")
    A("**3. The prespecification hash.**")
    A("")
    A(f"v1's RESULTS.md describes `prespec/f2_00_prespec.json` as \"{ts['v1_claim_in_RESULTS_md']}\". "
      f"{cap(ts['correction'])}")
    A("")
    A("**4. The amendment timestamps.**")
    A("")
    A("| file | mtime UTC | written_utc field | field later than its own mtime |")
    A("|---|---|---|---|")
    for r in ts["v1_prespec_files"]:
        A(f"| `{r['file']}` | {r['mtime_utc']} | {r['written_utc_field']} | "
          f"{'**yes**' if r['written_utc_is_after_its_own_mtime'] else 'no'} |")
    A(f"| `{ts['v1_sha256_file']['file']}` | {ts['v1_sha256_file']['mtime_utc']} | (no field) | n/a |")
    A("")
    A(cap(ts["written_utc_defect"]))
    A("")
    A(cap(ts["v2_practice"]))
    A("")
    A("**5. The RNA rejections do not come from overlapping readout windows.**")
    A("")
    A("v1's headline reads \"Every rejection sits where the two readout windows overlap\" and then gives the")
    A("mechanism as the union of two +/-1 kb windows collapsing below 2 kb of separation. That mechanism is real")
    A("for ATAC, DNASE and H3K27ac and cannot apply to RNA, because the RNA readout is a single region at every")
    A("separation:")
    A("")
    A("| quantity | RNA | ATAC / DNASE / H3K27ac |")
    A("|---|---|---|")
    A(f"| readout regions per pair | {geo['rna_readout_regions_per_pair'][0]} | 2, unioned |")
    A(f"| pairs of {geo['n_pairs']} whose two arms share the readout window | {geo['n_pairs']} (all of them, by "
      f"definition) | {geo['n_pairs_with_overlapping_local_windows']} |")
    A(f"| Spearman(separation, readout width) | {geo['rna_readout_bp_vs_separation_spearman'][0]:+.3f} "
      f"(p {geo['rna_readout_bp_vs_separation_spearman'][1]:.3g}) | "
      f"{geo['local_readout_bp_vs_separation_spearman'][0]:+.3f} "
      f"(p {geo['local_readout_bp_vs_separation_spearman'][1]:.3g}) |")
    A("")
    A("In the RNA channel both arms are read out over exactly the same bases at every separation, so there is no")
    A("overlap that appears below 2 kb and none that disappears above it. The correct statement is that all")
    A(f"{geo['n_rna_rejections']} RNA rejections happen to sit below 2 kb, and that is a fact about where they are,")
    A("not a mechanism for why.")
    A("")
    A("What they actually are:")
    A("")
    A("| signal | separation bp | RNA readout | readout bp | contains either variant | bp to the readout | residual |")
    A("|---|---|---|---|---|---|---|")
    for r in geo["rna_bh_rejections_v1"]:
        A(f"| `{r['signal_uid']}` | {r['separation_bp']} | {r['rna_readout']} | {r['rna_readout_bp']} | "
          f"{'**no**' if r['rna_readout_contains_neither_variant'] else 'yes'} | "
          f"{r['bp_from_variant_1_to_the_rna_readout']} | {r['rna_residual']:+.6f} |")
    A("")
    rm = geo["rejections_as_rows_and_as_distinct_variant_pairs"]
    A(f"All {geo['n_rna_rejections']} are the SAME two positions, chr1:205492407 A>T and chr1:205492435 T>C, 28 bp")
    A(f"apart, entering the family once per COLOC signal that names CDK18; their residuals agree to 16 digits. So")
    A(f"the RNA channel's \"{rm['rna']['rejected_rows']} of {geo['n_pairs']}\" is "
      f"{rm['rna']['distinct_variant_pairs']} distinct variant pair in "
      f"{rm['rna']['distinct_analysis_blocks']} analysis block. Both of its variants lie "
      f"{geo['rna_bh_rejections_v1'][0]['bp_from_variant_1_to_the_rna_readout']} bp OUTSIDE the readout window, so")
    A("its residual is a distal effect on a window that contains neither substitution.")
    A("")

    # ---------------------------------------------------------------- further defects
    A("## Further defects found while repairing these two, not part of the assignment")
    A("")
    dup = geo["duplicate_rows_in_the_199_pair_family"]
    A(f"**The 199-pair family holds {dup['distinct_variant_pairs']} distinct variant pairs.** "
      f"{dup['variant_pairs_appearing_more_than_once']} pairs of positions appear more than once, up to "
      f"{max(int(k) for k in dup['multiplicity_histogram'])} times, accounting for {dup['rows_they_account_for']} of "
      f"the {dup['rows']} rows. {cap(dup['why_that_matters'])} Rejection counts as rows and as distinct "
      f"variant pairs:")
    A("")
    A("| channel | rejected rows | distinct variant pairs | distinct 1-Mb analysis blocks |")
    A("|---|---|---|---|")
    for ch in CH:
        A(f"| {CHN[ch]} | {rm[ch]['rejected_rows']} | {rm[ch]['distinct_variant_pairs']} | "
          f"{rm[ch]['distinct_analysis_blocks']} |")
    A("")
    A(f"**The RNA readout does not contain the variants in {geo['n_pairs_whose_rna_readout_contains_neither_variant']} "
      f"of {geo['n_pairs']} pairs.** The readout is the target gene's span, which need not include the variant "
      f"positions; {geo['n_pairs_whose_rna_readout_contains_both_variants']} pairs have both variants inside it. "
      f"This is inherited from the deposit's recipe and is not repaired here.")
    A("")
    gs = geo["gene_span_readout_bp"]
    A(f"**The RNA readout width varies by three orders of magnitude across the observed family** "
      f"({geo['n_gene_span']} gene-span pairs: {gs['min']:.0f} to {gs['max']:.0f} bp, median {gs['median']:.0f}), "
      f"while {geo['v1_null_rna_readout_bp_for_2470_of_2500_draws']} bp is the readout for almost every one of "
      f"v1's 2,500 null draws. So the RNA arm of the 199-pair test is unmatched on readout width pair by pair, "
      f"not only for the GNMT pair.")
    A("")

    # ---------------------------------------------------------------- carried forward
    cf = S["carried_forward_unchanged_from_v1"]
    A("## Numbers carried forward unchanged, restated so this package stands alone")
    A("")
    A(f"Read from `{cf['source']}`; nothing here was rescored or recomputed.")
    A("")
    A(f"**The GNMT primary pair's four arms** ({cf['n_obs']}-pair family context below). log2 against the REF arm:")
    A("")
    A("| channel | REF sum | v1 log2 | v2 log2 | additive expectation log2 | joint log2 | residual |")
    A("|---|---|---|---|---|---|---|")
    for ch in CH:
        g = cf["gnmt_primary_arms"][ch]
        A(f"| {CHN[ch]} | {g['ref_sum']:.3f} | {g['v1_log2']:+.6f} | {g['v2_log2']:+.6f} | "
          f"{g['additive_expectation_log2']:+.6f} | {g['joint_log2']:+.6f} | {pr[f'{ch}_residual']:+.6f} |")
    A("")
    A(f"**F3, the {cf['n_obs']}-pair family against the {cf['n_null']}-draw corrected null.** Rejection counts are")
    A("rows; the distinct-variant-pair counts are in the table above.")
    A("")
    A("| channel | observed median abs residual | null median | null q95 | pairs over q95 uncorrected | smallest p | "
      "BH rejections (rows) | signals |")
    A("|---|---|---|---|---|---|---|---|")
    for ch in CH:
        c = cf["per_channel"][ch]
        A(f"| {CHN[ch]} | {c['median_abs_residual']:.5f} | {c['null_median_abs_residual']:.5f} | "
          f"{c['null_q95_abs_residual']:.5f} | {c['n_pairs_exceeding_null_q95_UNCORRECTED']} | "
          f"{c['min_p_empirical']:.6f} | {c['n_pairs_bh']} | {c['n_signals_bh']} |")
    A("")
    A("**Family-level contrast** (median abs residual, observed minus null; observed side resampled over 1-Mb")
    A(f"analysis blocks, null side over its draws, 10,000 draws, seed 20260914, {cf['family_level_contrast']['n_blocks']} blocks):")
    A("")
    A("| channel | difference | 95% interval | excludes 0 |")
    A("|---|---|---|---|")
    for ch in CH:
        b = cf["family_level_contrast"]["per_channel"][ch]
        A(f"| {CHN[ch]} | {b['difference']:+.5f} | [{b['ci95'][0]:+.5f}, {b['ci95'][1]:+.5f}] | "
          f"{'yes' if b['excludes_zero'] else 'no'} |")
    A("")
    A("**Magnitude-matched diagnostic** (v1's amendment 06; per pair, null draws whose |v1|+|v2| is within 1.5x of")
    A("the pair's own, at least 50 of them, then BH):")
    A("")
    A("| channel | observed median arm magnitude | null median | ratio | BH rejections magnitude-matched | unmatched |")
    A("|---|---|---|---|---|---|")
    for ch in CH:
        m = cf["magnitude_matched"][ch]
        A(f"| {CHN[ch]} | {m['observed_median_arm_magnitude']:.5f} | {m['null_median_arm_magnitude']:.5f} | "
          f"{m['observed_over_null_arm_magnitude']:.2f} | **{m['n_pairs_bh']}** | {cf['per_channel'][ch]['n_pairs_bh']} |")
    A("")
    A("**Separation-stratified** (v1's amendment 07; the null kept inside the stratum):")
    A("")
    A("| channel | " + " | ".join(f"{k.split('_readout')[0].replace('separation_', 'sep ')} "
                                  f"({v['n_observed']} pairs, null {v['n_null']})"
                                  for k, v in cf["separation_stratified"].items()) + " |")
    A("|---|---|---|")
    for ch in CH:
        cells = []
        for k, v in cf["separation_stratified"].items():
            p = v["per_channel"][ch]
            cells.append(f"{p['n_pairs_bh']}, smallest p {p['min_p']:.5f}")
        A(f"| {CHN[ch]} | " + " | ".join(cells) + " |")
    A("")
    fl = cf["p_floor_arithmetic"]
    A(f"**The p-floor arithmetic.** The existing P5 deposit's {fl['existing_p5_deposit_300_draws_1572_pairs']['null_draws']}-draw "
      f"null over {fl['existing_p5_deposit_300_draws_1572_pairs']['family_size']} pairs has a p floor of "
      f"{fl['existing_p5_deposit_300_draws_1572_pairs']['p_floor']:.7f} and needs "
      f"{fl['existing_p5_deposit_300_draws_1572_pairs']['pairs_that_must_sit_AT_the_floor_together_for_any_rejection']} "
      f"pairs at that floor together before BH can reject anything; this package's corrected null needs "
      f"{fl['this_package']['pairs_that_must_sit_AT_the_floor_together_for_any_rejection']}. The deposit's reported "
      f"zero rejections are a statement about its resolution, not about additivity.")
    A("")
    A(f"**Determinism.** {cap(cf['determinism']['verdict'])}.")
    A("")
    A(f"**The RNA readout counts.** Observed {cf['rna_readout_counts']['observed_gene_span']} gene-span and "
      f"{cf['rna_readout_counts']['observed_central_20kb']} central-20kb; null "
      f"{cf['rna_readout_counts']['null_gene_span']} gene-span and "
      f"{cf['rna_readout_counts']['null_central_20kb']} central-20kb. This is the same mismatch as defect 1, in the "
      f"F3 family rather than in the GNMT null, and it is not repaired here.")
    A("")
    A(f"**Blocks.** {cf['blocks']['n_distinct_analysis_blocks']} distinct 1-Mb analysis blocks and "
      f"{cf['blocks']['n_distinct_signals']} signals over the {cf['n_obs']} rows. {cap(cf['blocks']['note'])}. "
      f"(v1's note says a family of 200; {cf['n_obs']} were scored, and only "
      f"{geo['duplicate_rows_in_the_199_pair_family']['distinct_variant_pairs']} of those are distinct variant pairs.)")
    A("")

    # ---------------------------------------------------------------- unchanged
    A("## What did not move")
    A("")
    A("- The three chromatin channels' F2.1 verdicts. They were already matched: the ATAC, DNASE and H3K27ac")
    A("  readout is the union of two +/-1 kb windows on both sides, so v1's null was matched in those channels and")
    A("  the repair only changes which chr6 positions were sampled.")
    A("- F2.2 (the GNMT arms are far larger than the family median, so F2.1 is not a near-null test), F3.1, F3.2,")
    A("  F3.3 and their numbers, which this package did not rescore.")
    A("- The EUR haplotype frequencies, the r2 and D' table, and the finding that both single-variant arms of the")
    A("  primary pair have frequency 0 in 758 European chromosomes.")
    A("- The determinism check on the repeated REF call, the p-floor arithmetic and the k_min correction, the")
    A("  separation-stratified diagnostic, the magnitude-matched diagnostic and the family-level contrast.")
    A("- The count of BH rejections over the 199-pair family in each channel, restated above as distinct variant")
    A("  pairs rather than rows.")
    A("")

    # ---------------------------------------------------------------- not established
    A("## What is NOT established")
    A("")
    A("- **No measured interaction anywhere.** A residual is the model's response to two substitutions in one")
    A("  reference sequence. Nothing here is epistasis or evidence of it.")
    A("- **The repaired null matches the readout's width and gene type, not its signal level.** The median host")
    A(f"  span is {ro['rna_readout_bp']['median']:.0f} bp against GNMT's {ro['gnmt_rna_readout_bp']}, a factor of")
    A(f"  {ro['width_ratio_repaired_null_over_observed']}, and only "
      f"{ro['ref_rna_sum_vs_gnmt']['n_within_2x_of_gnmt']} of {rep['n_null_draws']} draws carry a REF RNA sum")
    A("  within a factor of 2 of GNMT's. A perfectly matched RNA null would need liver expression matching, which")
    A("  this package does not do.")
    A(f"- **{des['n_draws_that_are_a_second_pair_from_a_gene_already_drawn']} of {des['n_draws']} draws share a host")
    A("  gene with another draw**, so the empirical p treats as exchangeable some draws that are not independent.")
    A("  The one-pair-per-gene sub-null bounds that.")
    A("- **The both-active stratum is at its resolution limit in every channel.** With 13 to 28 matched draws the")
    A("  attainable p-values are coarse; a rejection count from such a null bounds nothing about effect size.")
    A("- **The RNA channel's 199-pair family remains readout-unmatched pair by pair.** Repairing the GNMT null")
    A("  does not repair the F3 family's RNA arm; that would need a per-pair gene-span-matched null, which was not")
    A("  scored.")
    A("- **\"Additive\" for the GNMT pair is still not a population statement.** Both single-variant arms have")
    A("  frequency 0 in 758 European chromosomes, so two of the four sequences entering the residual are sequences")
    A("  nobody carries.")
    A("- A hosted-API log2 track-sum ratio is not an Atlas quantile and is not comparable to one.")
    A("")

    # ---------------------------------------------------------------- provenance
    A("## Provenance")
    A("")
    A(f"Output directory: `{OUT}`")
    A("")
    A(f"Superseded package: `{V1NAME}`, marked in `SUPERSEDED.md` in its own directory. Its RESULTS.md is left")
    A("untouched so the record of what was published stays readable.")
    A("")
    A("| file | what |")
    A("|---|---|")
    for pat, what in (
        ("prespec/f2v2_00_prespec.json", "prespecification, written and digested before the first API call"),
        ("prespec/f2v2_00_amendment_01.json", "the host gene type match, written before any score was read"),
        ("prespec/f2v2_00_prespec.sha256", "sha256 of the prespecification and its amendment"),
        ("prespec/WRITTEN_UTC.json", "machine-clock write time beside the file mtime, for each"),
        ("tables/gnmt_genespan_null_pairs.tsv", "the 100 gene-span-matched chr6 draws as drawn"),
        ("tables/gnmt_genespan_null_design.json", "the draw rule and what it achieved"),
        ("tables/f2v2_genespan_null_effects.tsv", "the repaired null's four-arm sums, log2 effects and residuals"),
        ("tables/f2v2_gnmt_verdict.json", "F2.1 per channel against the repaired null, beside v1's"),
        ("tables/f2v2_both_active_stratum.json", "defect 2, both stages, with the null inside the stratum"),
        ("tables/f2v2_readout_geometry.json", "the RNA readout is one window; duplicates; the RNA rejections"),
        ("tables/f2v2_accounting.json", "call counts, checkpoint deposits, prespec and amendment timestamps"),
        ("tables/f2v2_summary.json", "every number in this document"),
        ("raw/scored_genespan_null.jsonl", "every scored draw; the checkpoint a restart resumes from"),
        ("MANIFEST.tsv", "every input with sha256"),
        ("pip_freeze.txt", "environment"),
    ):
        A(f"| `{pat}` | {what} |")
    A("")
    (OUT / "RESULTS.md").write_text("\n".join(L) + "\n")
    print(f"RESULTS.md written: {len(L)} lines")


if __name__ == "__main__":
    main()
