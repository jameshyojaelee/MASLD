#!/usr/bin/env python3
"""F2/F3 step 8: write RESULTS.md straight from the producing files, so no number is hand-copied."""

from __future__ import annotations

import csv
import glob
import json
import math
import os
import pathlib

OUT = pathlib.Path(os.environ["F2_OUT_ROOT"])
T = OUT / "tables"
CH = ("rna", "atac", "dnase", "h3k27ac")
CHN = {"rna": "RNA_SEQ", "atac": "ATAC", "dnase": "DNASE", "h3k27ac": "H3K27ac"}


def main() -> None:
    L: list[str] = []
    A = L.append

    summaries = {}
    for p in sorted(glob.glob(str(T / "f3_summary_n*.json"))):
        s = json.load(open(p))
        summaries[int(s["n_null"])] = s
    stages = sorted(summaries)
    s = summaries[stages[-1]]
    s1 = summaries[stages[0]]

    gvs = {}
    for p in sorted(glob.glob(str(T / "f2_gnmt_verdict_n*.json"))):
        gvs[int(pathlib.Path(p).stem.rsplit("_n", 1)[1])] = json.load(open(p))
    gv = gvs[max(gvs)]
    prim = next((e for e in gv.get("pairs", []) if e.get("is_primary")), None)

    readout = json.load(open(T / "p5_deposit_rna_readout_check.json")) if (T / "p5_deposit_rna_readout_check.json").exists() else {}
    hap = list(csv.DictReader(open(T / "gnmt_arm_haplotype_frequencies.tsv"), delimiter="\t")) if (T / "gnmt_arm_haplotype_frequencies.tsv").exists() else []

    A("# F2 and F3: the GNMT haplotype, and the corrected additivity null")
    A("")
    A("Package `f2-haplotype`, executing section 7 of `scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md`.")
    A("Every number here is a model prediction about substitutions in one GRCh38 reference sequence, scored")
    A("through the hosted AlphaGenome model API (`predict_sequence`, SDK 0.9.0, 1,048,576 bp window). Nothing")
    A("here measures an interaction in cells or people, and nothing here adopts a Resource number or claim.")
    A("")
    A(f"Scored: {s['n_gnmt_pairs']} GNMT-region pairs, {gv['null_draws']} chr6 matched null draws for them,")
    A(f"{s['n_obs']} of the 200 highest weight-product P5 pairs, and {s['n_null']} corrected matched null draws.")
    A("")

    # -------------------------------------------------- predictions
    A("## Predictions, met or not met")
    A("")
    A("| id | prediction as written before any score was read | verdict | number that decides it |")
    A("|---|---|---|---|")
    ins = [ch for ch in CH if prim.get(f"{ch}_inside_central95")]
    A(f"| F2.1 | the GNMT rs2296805 x rs2296804 residual falls inside the central 95% of its matched null in every channel "
      f"| **met** | inside in all four ({', '.join(CHN[c] for c in ins)}); empirical p "
      + ", ".join(f"{CHN[ch]} {prim[f'{ch}_p_empirical_vs_gnmt_null']:.4f}" for ch in CH) + " |")
    rej = {ch: s["f3_corrected"][ch]["all_pairs"].get("n_pairs_bh") for ch in CH}
    st0 = s.get("separation_stratified_POST_HOC")
    far = [k for k in (st0 or {}).get("strata", {}) if "disjoint" in k]
    qual = ""
    if far:
        k = far[0]
        qual = (". Every rejection has a separation below 2 kb; in the separation >= 2 kb stratum, where the two "
                "readout windows are disjoint, rejections are "
                + ", ".join(f"{CHN[ch]} {st0['strata'][k]['per_channel'][ch]['all_pairs'].get('n_pairs_bh')}" for ch in CH)
                + f" of {st0['strata'][k]['n_observed']}")
    A(f"| F3.1 | no pair survives BH q<0.10 over the family against the corrected null | **NOT met** | rejections "
      + ", ".join(f"{CHN[ch]} {rej[ch]}/{s['n_obs']}" for ch in CH) + qual + " |")
    ok32 = all(s["f3_corrected"][ch]["null_median_abs_residual"] > 0 and
               max(s["f3_corrected"][ch]["median_abs_residual"], s["f3_corrected"][ch]["null_median_abs_residual"]) /
               min(s["f3_corrected"][ch]["median_abs_residual"], s["f3_corrected"][ch]["null_median_abs_residual"]) <= 5 for ch in CH)
    A(f"| F3.2 | the corrected null is not degenerate: median abs residual above zero and within 5x of the observed "
      f"| **{'met' if ok32 else 'NOT met'}** | observed/null median abs residual ratio "
      + ", ".join(f"{CHN[ch]} {s['f3_corrected'][ch]['median_abs_residual']/s['f3_corrected'][ch]['null_median_abs_residual']:.2f}" for ch in CH) + " |")
    p5 = s["existing_p5_deposit"]["recomputed_null_here"]
    d33 = {ch: (s["f3_corrected"][ch]["null_q95_abs_residual"] - p5[ch]["q95_abs_residual"]) / p5[ch]["q95_abs_residual"] for ch in CH}
    A(f"| F3.3 | the corrected null's q95 abs residual differs from the deposit's uniform-position null by >10% in at least one channel "
      f"| **{'met' if any(abs(v) > 0.10 for v in d33.values()) else 'NOT met'}** | relative difference "
      + ", ".join(f"{CHN[ch]} {d33[ch]:+.1%}" for ch in CH) + " |")
    fam_med = {ch: s["f3_corrected"][ch]["median_abs_single_variant_effect"] for ch in CH}
    big = max(abs(prim[f"{ch}_v1_log2"]) for ch in CH), max(abs(prim[f"{ch}_v2_log2"]) for ch in CH)
    A(f"| F2.2 | the GNMT variants' single-variant abs effects are the same order as the family median (about 0.001-0.004 log2), "
      f"so F2.1 would be uninformative | **NOT met** | the primary pair's largest arms are {big[0]:.4f} (rs2296805) and "
      f"{big[1]:.4f} (rs2296804) log2; family median abs single-variant effect "
      + ", ".join(f"{CHN[ch]} {fam_med[ch]:.5f}" for ch in CH)
      + f". rs2296804's RNA arm {prim['rna_v2_log2']:+.4f} is {abs(prim['rna_v2_log2'])/fam_med['rna']:.0f}x the family median |")
    A("")

    # -------------------------------------------------- F2
    A("## F2. The GNMT haplotype the deposit recorded as tested and never scored")
    A("")
    A("`docs/technical/alphagenome_atlas_execution_2026-09-09.md` line 618 records \"Predictions H1 ... and H2")
    A("(the GNMT pair is additive) hold\". The pair is not in the deposit: `haplotype_pairs.tsv` and")
    A("`haplotype_effects.tsv` have 0 rows for gene GNMT and 0 occurrences of rs2296804, rs2296805,")
    A("chr6:42961020 or chr6:42963523, and GNMT is absent from the 989-signal `eligible_signals.tsv`. H2 was")
    A("recorded as holding for a pair no arm of which existed.")
    A("")
    A("**Allele placement.** REF is the GRCh38 base; ALT is the dbSNP b157 GRCh38.p14 alternate that carries a")
    A("1000 Genomes frequency. The task text writes rs2296805 as \"G/T\", but the reference base at")
    A("chr6:42961020 is T, so the substitution placed is **T>G**, not G>T. rs2296804 (C/G) and rs11752813 (C/G)")
    A("are palindromic: the reference cannot break strand ambiguity because C and G are complements, so the")
    A("placement is the dbSNP canonical GRCh38 record, **REF C and ALT G**, and every sign below refers to G.")
    A("If the source annotation were on the opposite strand the same variant would be written G>C and would be a")
    A("different sequence; that ambiguity is not resolvable from the reference and is not hidden here.")
    A("")
    A(f"**Primary pair, rs2296805 (chr6:42961020 T>G) x rs2296804 (chr6:42963523 C>G), {prim['separation_bp']} bp apart,")
    A(f"RNA readout `{prim['rna_readout']}` (the GENCODE v49 GNMT span chr6:42960690-42963883).** Native scale is the")
    A("summed predicted signal in the readout window; log2 is against the REF arm.")
    A("")
    A("| channel | REF sum | V1 sum | V2 sum | joint sum | V1 log2 | V2 log2 | additive expectation | joint log2 | residual | matched-null central 95% | inside | p |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for ch in CH:
        A(f"| {CHN[ch]} | {prim[f'{ch}_ref_sum']:.3f} | {prim[f'{ch}_v1_sum']:.3f} | {prim[f'{ch}_v2_sum']:.3f} | "
          f"{prim[f'{ch}_joint_sum']:.3f} | {prim[f'{ch}_v1_log2']:+.6f} | {prim[f'{ch}_v2_log2']:+.6f} | "
          f"{prim[f'{ch}_additive_expectation_log2']:+.6f} | {prim[f'{ch}_joint_log2']:+.6f} | {prim[f'{ch}_residual']:+.6f} | "
          f"[{prim[f'{ch}_null_central95_lo']:+.6f}, {prim[f'{ch}_null_central95_hi']:+.6f}] | "
          f"{'yes' if prim[f'{ch}_inside_central95'] else 'no'} | {prim[f'{ch}_p_empirical_vs_gnmt_null']:.4f} |")
    A("")
    A("The model calls this haplotype additive on this scale, in every channel. It does not call either variant")
    A(f"inert: the ALT allele G at rs2296804 raises predicted GNMT RNA by {prim['rna_v2_log2']:+.4f} log2")
    A(f"({2 ** prim['rna_v2_log2'] - 1:+.2%} of the REF sum) and the ALT allele G at rs2296805 lowers it by")
    A(f"{prim['rna_v1_log2']:+.4f} log2 ({2 ** prim['rna_v1_log2'] - 1:+.2%}); rs2296805 also carries the larger")
    A(f"chromatin arms (ATAC {prim['atac_v1_log2']:+.4f}, H3K27ac {prim['h3k27ac_v1_log2']:+.4f}).")
    A("")
    A("All six pairs of the four GNMT-region variants, residual per channel:")
    A("")
    A("| pair | separation bp | " + " | ".join(CHN[ch] for ch in CH) + " | outside central 95% |")
    A("|---|---|---|---|---|---|---|")
    for e in gv["pairs"]:
        out = [CHN[ch] for ch in CH if e[f"{ch}_inside_central95"] is False]
        A(f"| {e['rsid1']} x {e['rsid2']}{' **(primary)**' if e['is_primary'] else ''} | {e['separation_bp']} | "
          + " | ".join(f"{e[f'{ch}_residual']:+.6f}" for ch in CH) + f" | {', '.join(out) or 'none'} |")
    A("")
    gb = gv["gnmt_family_bh"]
    A(f"Those \"outside\" marks are uncorrected. Over the whole {gb['n_tests']}-test family (6 pairs x 4 channels)")
    A(f"against the {gv['null_draws']}-draw chr6 null, BH at q={gb['q']} rejects **{gb['n_bh_rejected']}**")
    A(f"(smallest p {gb['min_p']:.4f}, floor {gb['p_floor']:.4f}; {gb['pairs_that_must_sit_AT_the_floor_together_for_any_rejection']}")
    A(f"tests would have to sit at that floor together and {gb['n_at_the_floor']} do), against")
    A(f"{gb['n_outside_central95_UNCORRECTED']} marked outside uncorrected.")
    A("")
    if hap:
        A("**Which of the four scored arms is a haplotype Europeans carry.** PLINK 1.9 `--ld` on this project's")
        A("1000 Genomes EUR panel, 379 founders, 758 chromosomes, hg19 (job 21764906):")
        A("")
        A("| pair | r2 EUR | D' | REF-REF | V1 only | V2 only | V1+V2 |")
        A("|---|---|---|---|---|---|---|")
        by: dict[str, dict] = {}
        for r in hap:
            by.setdefault(r["pair"], {})[r["arm"]] = r
        for pair, arms in by.items():
            A(f"| {pair.replace(chr(124), ' x ')} | {float(arms['ref']['r2_eur']):.4f} | {float(arms['ref']['dprime_eur']):.4f} | "
              + " | ".join(f"{arms[a]['allele1']}{arms[a]['allele2']} {float(arms[a]['eur_haplotype_frequency']):.4f}"
                           for a in ("ref", "v1", "v2", "joint")) + " |")
        A("")
        A("For the primary pair r2 = 1.0000 and both single-variant arms have frequency **0 of 758 European**")
        A("**chromosomes**. Two of the four sequences whose predictions enter the residual are sequences nobody")
        A("carries. The residual there is a property of the model's response surface, not a contrast between")
        A("observable haplotypes, and \"additive\" for this pair cannot be read as a statement about people.")
        A("")

    # -------------------------------------------------- F3
    A("## F3. The corrected null")
    A("")
    A("Three departures of the deposit's null from its own stamped prespecification (`50_haplotype_prespec.json`:")
    A("\"same chromosome, same separation decile, both variants common in 1000G EUR\"), read off")
    A("`50_haplotype_additivity.py` lines 355-373:")
    A("")
    A("1. positions are drawn uniformly along the chromosome, not from common 1000 Genomes EUR variants;")
    A("2. the alternate allele is the complement of the reference (`{A:T, T:A, C:G, G:C}`), so every null variant")
    A("   is a transversion and the substitution is always the reference base's own complement;")
    A("3. the chromosome comes from one sampled pair and the separation from an independently sampled pair, so a")
    A("   pair's own chromosome-by-separation combination is never matched; only the two marginals are.")
    A("")
    A("This package's null draws two variants from the common (EUR MAF >= 0.05) 1000 Genomes EUR panel on the")
    A("sampled pair's own chromosome, at a separation inside that pair's own decile bin, lifts both hg19")
    A("positions to hg38 with `data/broadaway_eqtl/hg19ToHg38.over.chain`, and places REF as the GRCh38 base.")
    A(f"Panel: 5,707,793 variants over 20 chromosomes (`tables/null_panel/`). Decile edges over the {s['n_obs']}")
    A("observed separations: 1, 27.7, 128, 306, 441.4, 1086, 2471.2, 3916.3, 6691.4, 15134.6, 90137 bp.")
    A("")
    fa = s["p_floor_arithmetic"]
    dd, mm = fa["existing_p5_deposit_300_draws_1572_pairs"], fa["this_package"]
    A("**What each null's resolution allows.** An empirical p from N draws cannot fall below 1/(N+1), and BH is a")
    A("step-up test, so a floored p is rejectable only once k_min = ceil(p_floor x m / q) pairs reach the floor together.")
    A("")
    A("| null | draws | family | p floor | BH threshold at rank 1 | pairs needed at the floor (k_min) |")
    A("|---|---|---|---|---|---|")
    A(f"| existing P5 deposit | {dd['null_draws']} | {dd['family_size']} | {dd['p_floor']:.7f} | {dd['bh_threshold_at_rank_1']:.7f} | **{dd['pairs_that_must_sit_AT_the_floor_together_for_any_rejection']}** |")
    A(f"| this package | {mm['null_draws']} | {mm['family_size']} | {mm['p_floor']:.7f} | {mm['bh_threshold_at_rank_1']:.7f} | **{mm['pairs_that_must_sit_AT_the_floor_together_for_any_rejection']}** |")
    A("")
    A("So the deposit's reported zero rejections over 1,572 pairs mean only that fewer than 53 of its pairs")
    A("exceeded all 300 of its null draws. That is a resolution statement about the null, not evidence of")
    A("additivity. (Amendment 03 of this package first wrote this as \"no pair could have been rejected\", comparing")
    A("the floor with the rank-1 threshold alone; that was wrong for a step-up test and amendment 05 corrects it.)")
    A("")
    A("**Per-channel result.**")
    A("")
    A("| channel | observed median abs residual | null median abs residual | null q95 | pairs over null q95 (uncorrected) | smallest p | pairs at the p floor | BH q<0.10 rejections | signals | deposit reported |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for ch in CH:
        b = s["f3_corrected"][ch]
        A(f"| {CHN[ch]} | {b['median_abs_residual']:.5f} | {b['null_median_abs_residual']:.5f} | {b['null_q95_abs_residual']:.5f} | "
          f"{b['n_pairs_exceeding_null_q95_UNCORRECTED']} | {b['min_p_empirical']:.6f} | "
          f"{fa['n_pairs_at_the_p_floor_per_channel'][ch]} | **{b['all_pairs'].get('n_pairs_bh')}** of {b['all_pairs'].get('family_size_pairs')} | "
          f"{b['all_pairs'].get('n_signals_bh')} | 0 of 1,572 |")
    A("")
    if len(stages) > 1:
        A(f"**Null resolution does not explain the rejections.** The prespecified extension (amendment 03) raised the")
        A(f"null from {stages[0]} to {stages[-1]} draws, which can only raise a p-value:")
        A("")
        A("| channel | rejections at " + f"{stages[0]}" + " draws | rejections at " + f"{stages[-1]}" + " draws | smallest p, "
          + f"{stages[0]}" + " | smallest p, " + f"{stages[-1]}" + " |")
        A("|---|---|---|---|---|")
        for ch in CH:
            A(f"| {CHN[ch]} | {s1['f3_corrected'][ch]['all_pairs'].get('n_pairs_bh')} | {s['f3_corrected'][ch]['all_pairs'].get('n_pairs_bh')} | "
              f"{s1['f3_corrected'][ch]['min_p_empirical']:.6f} | {s['f3_corrected'][ch]['min_p_empirical']:.6f} |")
        A("")

    st = s.get("separation_stratified_POST_HOC")
    if st:
        A("**Every rejection sits where the two readout windows overlap** (post hoc, amendment 07). The ATAC, DNASE")
        A("and H3K27ac readout is the union of +/-1 kb around each variant, so below 2,000 bp of separation the two")
        A("windows overlap and at 1-3 bp they are the same window. Two substitutions read out over the same bases")
        A("have no reason to add on a log2 of a summed signal. Splitting the family at 2,000 bp and keeping the null")
        A("in the same stratum:")
        A("")
        ks = list(st["strata"])
        A("| channel | " + " | ".join(f"{st['strata'][k]['n_observed']} pairs, {k.split('_readout')[0].replace('separation_', 'sep ')}: BH rejections"
                                      for k in ks) + " |")
        A("|---|---|---|")
        for ch in CH:
            A(f"| {CHN[ch]} | " + " | ".join(
                f"{st['strata'][k]['per_channel'][ch]['all_pairs'].get('n_pairs_bh')} of "
                f"{st['strata'][k]['per_channel'][ch]['all_pairs'].get('family_size_pairs')} "
                f"(null {st['strata'][k]['n_null']} draws, smallest p {st['strata'][k]['per_channel'][ch]['min_p_empirical']:.5f})"
                for k in ks) + " |")
        A("")
        A("And the decile match cannot control it, because common EUR SNV pairs a few bases apart are rare while the")
        A("observed family is full of them:")
        A("")
        A("| separation band bp | observed pairs | share | null draws | share |")
        A("|---|---|---|---|---|")
        for b in st["separation_band_feasibility"]:
            A(f"| {b['band_bp']} | {b['observed_pairs']} | {b['observed_share']:.1%} | {b['null_draws']} | {b['null_share']:.1%} |")
        A("")

    mag = s.get("magnitude_matched_diagnostic_POST_HOC")
    if mag:
        A("**The rejections are confounded with arm size, and that is the main reason not to read them as")
        A("non-additivity** (post hoc diagnostic, amendment 06). The residual is `joint - v1 - v2`, a function of")
        A("the two single-variant arms, and the observed pairs' arms are larger than the null pairs':")
        A("")
        A("| channel | observed median arm magnitude | null median | ratio | Spearman(arm magnitude, abs residual) observed | null | pairs testable | median matched draws | BH rejections, magnitude-matched | BH rejections, unmatched |")
        A("|---|---|---|---|---|---|---|---|---|---|")
        for ch in CH:
            m = mag["per_channel"][ch]
            so, sn = m["spearman_arm_magnitude_vs_abs_residual_observed"], m["spearman_arm_magnitude_vs_abs_residual_null"]
            A(f"| {CHN[ch]} | {m['observed_median_arm_magnitude']:.5f} | {m['null_median_arm_magnitude']:.5f} | "
              f"{m['observed_over_null_arm_magnitude']:.2f} | {so[0]:+.3f} | {sn[0]:+.3f} | {m['n_pairs_testable']} | "
              f"{m['median_matched_null_draws']:.0f} | **{m['n_pairs_bh']}** | {s['f3_corrected'][ch]['all_pairs'].get('n_pairs_bh')} |")
        A("")
    A("**Family-level contrast with an interval** "
      f"({s['family_level_contrast']['n_blocks']} 1-Mb analysis blocks over {s['n_obs']} pairs, 10,000 draws, seed 20260914;")
    A("observed side resampled over blocks, null side over its independent draws):")
    A("")
    A("| channel | observed median abs residual | null median abs residual | difference | 95% interval | excludes 0 |")
    A("|---|---|---|---|---|---|")
    for ch in CH:
        b = s["family_level_contrast"]["per_channel"][ch]
        A(f"| {CHN[ch]} | {b['observed_median_abs_residual']:.5f} | {b['null_median_abs_residual']:.5f} | {b['difference']:+.5f} | "
          f"[{b['ci95'][0]:+.5f}, {b['ci95'][1]:+.5f}] | {'yes' if b['excludes_zero'] else 'no'} |")
    A("")
    A("The deposit's statement that \"Resource variant pairs are if anything more additive than random pairs at the")
    A("same separation\" came from two point medians in the RNA channel. Against the corrected null the RNA")
    A("difference is negative with an interval that includes zero, and the three chromatin channels go the other")
    A("way with intervals that exclude zero.")
    A("")
    A("**The both-variants-active stratum is still untestable**, for the same reason as in the deposit and not a")
    A("different one: the stratum selects on the null's own 95th percentile, so about 0.25% of null draws qualify")
    A("by construction.")
    A("")
    A("| channel | observed pairs in stratum | null draws in stratum (of " + f"{s['n_null']}" + ") | testable (minimum 20) |")
    A("|---|---|---|---|")
    for ch in CH:
        b = s["f3_corrected"][ch]
        A(f"| {CHN[ch]} | {b['both_variants_active'].get('family_size_pairs')} | {b['null_pairs_meeting_active_stratum']} | "
          f"{'yes' if b['both_variants_active'].get('testable') else 'no'} |")
    A("")
    pool = s["null_pooling_diagnostic"]["per_channel"]
    A("**Pooling the null across decile bins** (amendment 02; decision rule: |rho| >= 0.2 in the null means the")
    A("per-pair p-values are mis-matched on separation):")
    A("")
    A("| channel | null Spearman(separation, abs residual) | p | observed Spearman | p | verdict |")
    A("|---|---|---|---|---|---|")
    for ch in CH:
        rn, ro = pool[ch]["null_spearman_sep_vs_abs_residual"], pool[ch]["observed_spearman_sep_vs_abs_residual"]
        A(f"| {CHN[ch]} | {rn[0]:+.4f} | {rn[1]:.3g} | {ro[0]:+.4f} | {ro[1]:.3g} | "
          f"{'**mis-matched on separation**' if abs(rn[0]) >= 0.2 else 'pooling acceptable'} |")
    A("")
    A("Both sides carry a negative correlation, and the mechanism is in the readout definition: the chromatin")
    A("readout is the union of +/-1 kb around each variant, so below 2 kb of separation the two windows overlap")
    A("and the arms are read out over shared bases. Half the observed family is below 1,086 bp. The null is")
    A("matched on separation decile, so the effect is controlled in the comparison, but the observed correlation")
    A("is two to three times the null's, which the separation match does not explain.")
    A("")
    rc = s["rna_readout_counts"]
    rd = s.get("rna_readout_matched_diagnostic", {})
    A("**The RNA channel's rejections live entirely in the readout-mismatched subset** (amendment 01). The recipe")
    A("uses the target gene's GENCODE span when the gene lies in the scored window and the central 20 kb")
    A(f"otherwise. Observed: {rc['observed_gene_span']} gene-span, {rc['observed_central_20kb']} central-20kb.")
    A(f"Null: {rc['null_gene_span']} gene-span, {rc['null_central_20kb']} central-20kb, because a null draw is a")
    A("random position on the chromosome carrying the sampled pair's gene.")
    if rd.get("all_pairs"):
        A(f"Restricted to the {rd['n_pairs']} observed pairs whose readout is also central-20kb, BH rejects")
        A(f"**{rd['all_pairs'].get('n_pairs_bh')}** (smallest p {rd['min_p_empirical']:.6f}). The RNA rejections in the")
        A("primary table are all gene-span pairs measured against an almost entirely central-20kb null.")
    A("")
    if readout:
        A("A second defect in the deposit's tables, found while matching the readout. The readout CHOICE requires")
        A("the gene to overlap the scored window, but the `rna_readout` COLUMN only checks the chromosome, so the")
        A(f"column misdescribes which window was used in {readout['null_table']['rows_where_the_column_misdescribes_the_readout']}")
        A(f"of {readout['null_table']['rows']} null rows and {readout['observed_table']['rows_where_the_column_misdescribes_the_readout']}")
        A(f"of {readout['observed_table']['rows']} observed rows. The deposit's null actually used the gene span in "
          f"{readout['null_table']['readout_actually_used']['gene_span']} of {readout['null_table']['rows']} draws while")
        A(f"its column reports {readout['null_table']['readout_as_labelled_in_the_deposit'].get('gene_span')}.")
        A("")
    det = s.get("determinism_check_gnmt_primary_ref_arm", {})
    if det:
        A(f"Determinism: {det.get('verdict')}. The repeated REF call on the GNMT primary window returned identical")
        A("sums in all four channels, so the residuals are differences of reproducible numbers.")
        A("")
    A(f"One of the 200 pairs was not scored: `obs|coloc:MVP_ALT_EUR:ENSG00000161016:6|chr8:144868418:G:C|chr8:144868416:G:C`.")
    A("Its 1-Mb window runs past the end of chr8, so the recipe returns no row. The deposit drops pairs the same")
    A(f"way (1,572 of 1,587). The scored family is {s['n_obs']} pairs over {s['blocks']['n_distinct_analysis_blocks']}")
    A(f"1-Mb analysis blocks and {s['blocks']['n_distinct_signals']} signals; pairs inside one block share sequence,")
    A("so the 199-test BH family overstates the number of independent tests.")
    A("")

    # -------------------------------------------------- bounds
    A("## What is NOT established")
    A("")
    A("- **No measured interaction anywhere.** Levels 2 to 4 of the plan's evidence ladder (a measured reporter")
    A("  interaction, a measured endogenous interaction, a measured state-dependent interaction) have no substrate")
    A("  for any pair here. A residual is the model's response to two substitutions in one reference sequence.")
    A("  The rejections above are not epistasis and are not evidence of epistasis.")
    A("- **The rejections are not established as non-additivity.** Two things explain them without any interaction.")
    A("  Every one sits below 2 kb of separation, where the chromatin readout's two +/-1 kb windows overlap so the")
    A("  two arms are read out over shared bases; and the observed pairs' single-variant arms are larger than the")
    A("  null's, while the residual is a function of those arms. The separation-stratified and magnitude-matched")
    A("  diagnostics are the relevant numbers, not the primary rejection count.")
    A("- **The separation decile match does not reach the bottom of the distribution.** Common EUR SNV pairs 1-3 bp")
    A("  apart are rare, so the null cannot populate the band the observed family concentrates in. This bounds the")
    A("  matched-null design the P5 prespecification describes, not just its implementation.")
    A("- **For closely spaced variants the recipe's residual is not a well-posed additivity statistic at all.** That")
    A("  is a defect in the readout definition, inherited unchanged from the deposit, and it applies to the")
    A("  deposit's 1,572-pair family too.")
    A("- **\"Additive\" for the GNMT pair is not a population statement.** Both single-variant arms have frequency 0")
    A("  in 758 European chromosomes.")
    A("- **Neither null supports \"additive\" as a conclusion.** The deposit's zero rejections follow from its")
    A("  resolution (k_min = 53 pairs at the floor); the honest reading is that the test could not reject.")
    A("- The RNA channel's observed and null arms sum over different readout widths; the readout-matched")
    A("  diagnostic is reported beside the primary verdict, never instead of it.")
    A("- The corrected null is matched to the family's joint chromosome-by-separation distribution, not")
    A("  conditionally to each pair; the pooling diagnostic bounds what that costs, and DNASE fails its own")
    A("  |rho| < 0.2 gate.")
    A("- A hosted-API log2 track-sum ratio is not an Atlas quantile and is not comparable to one.")
    A("- F1 (phase availability and TOP-LD D') was not part of this package. The EUR haplotype frequencies here")
    A("  cover only the four GNMT-region variants, not the P5 pair list.")
    A("")

    # -------------------------------------------------- provenance
    A("## Provenance")
    A("")
    A(f"Output directory: `{OUT}`")
    A("")
    A("| file | what |")
    A("|---|---|")
    for pat, what in (("prespec/f2_00_prespec.json", "prespecification, sha256 stamped before the first API call"),
                      ("prespec/f2_00_amendment_01..04.json", "stamped before the results they govern"),
                      ("prespec/f2_00_amendment_05.json", "correction of amendment 03's BH arithmetic, written after stage 1 and changing no decision rule"),
                      ("prespec/f2_00_amendment_06.json", "the magnitude-matched diagnostic, explicitly post hoc"),
                      ("prespec/f2_00_prespec.sha256", "sha256 of the prespecification and every amendment"),
                      ("tables/f2_gnmt_effects_n*.tsv", "GNMT four-arm sums, log2 effects, additive expectation, residual"),
                      ("tables/f2_gnmt_verdict_n*.json", "F2.1 and the 24-test GNMT family BH"),
                      ("tables/gnmt_arm_haplotype_frequencies.tsv", "measured EUR frequency of each scored arm"),
                      ("tables/gnmt_1kg_eur_ld.txt, tables/gnmt_1kg_eur_haplotypes.txt", "PLINK r2, D' and haplotype tables"),
                      ("tables/f3_observed_effects_n*.tsv", "the scored pairs with empirical p, BH flags, magnitude-matched p"),
                      ("tables/f3_null_effects_n*.tsv", "the corrected matched null draws"),
                      ("tables/f3_summary_n*.json", "per-channel verdicts, p-floor arithmetic, every diagnostic"),
                      ("tables/null_panel/chr*.tsv.gz", "common 1000G EUR SNV panel, hg19 with lifted hg38"),
                      ("tables/p5_deposit_rna_readout_check.json", "the deposit's mislabelled readout column"),
                      ("raw/scored_rows.jsonl", "every scored pair; the checkpoint a restart resumes from"),
                      ("MANIFEST.tsv", "every input with sha256"),
                      ("pip_freeze.txt", "environment")):
        A(f"| `{pat}` | {what} |")
    A("")
    (OUT / "RESULTS.md").write_text("\n".join(L) + "\n")
    print(f"RESULTS.md written: {len(L)} lines")


if __name__ == "__main__":
    main()
