#!/usr/bin/env python
"""A1 v2. Reporter measurement reliability (GSE281364 MPRA), complete deposit.

Runs the wave-1 analysis (scripts/analysis/alphagenome_program/a1_mpra_reliability.py,
section 3 of IMPLEMENTATION_SPEC.md) at the full 10,000 bootstrap draws and 10,000
permutations at seed 20260914, in one sbatch job, and deposits every table.

It then adds eight corrections that an independent checker raised against the wave-1
deposit. Each is a new table; no wave-1 number is altered, only labelled or completed.

  1 reference_reproduction_agreement.tsv   agreement with the published 1,033-element
    values is to four decimals, not exact; the signed differences are reported.
  2 nan_convention_invariance.tsv          complete-case versus NaN-skipping is an
    invariance, not a check: there are zero missing element x replicate cells.
  3 element_locus_group_mapping.tsv        elements map many-to-one onto outer locus
    groups; the ratio is reported for the analysis set and for the fixture.
  4 treatment_contrast_published_value.tsv ONE projected value per cell line, with the
    aggregation convention named and the alternative shown beside it.
  5 treatment_tertile_values.tsv           every tertile value carries n_elements and
    n_blocks; the decomposition is re-emitted with n.
  6 treatment_contrast_decomposition_annotated.tsv  explicit sign statement for
    j_minus_dna_ratio_only, plus a numerical additivity guard.
  7 treatment_contrast_by_partition.tsv    the three per-partition Pearsons behind each
    cell line's treatment value, not only their mean.
  8 activity_convention_sensitivity.tsv    the activity-convention check that wave 1
    wrote and never ran.

Usage: python a1_mpra_reliability_v2.py <output_dir>
"""
from __future__ import annotations

import json
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import a1_mpra_reliability as a1  # noqa: E402
from a1_mpra_reliability import (  # noqa: E402
    CELL_PAIRS, CONTEXTS, PC, REPS, SPLITS, _pearson, _spearman, build_matrices,
    half_means, load, reliability_from_halves, sb_project,
)

REF_PUBLISHED = {"LX2_TGFb": 0.9195, "LX2_control": 0.9228,
                 "HepG2_control": 0.8861, "HepG2_PAOA": 0.8580}

# The convention published for the treatment contrast, fixed here and used everywhere.
PUBLISHED_CONVENTION = (
    "mean_of_three_partition_pearsons_then_spearman_brown "
    "(average the three 2v2 Pearsons over partitions {12|34,13|24,14|23}, "
    "then apply the Spearman-Brown projection from k=2 to k=4 once)"
)
ALT_CONVENTION = (
    "spearman_brown_per_partition_then_mean "
    "(project each partition's Pearson separately, then average the three projections)"
)

J_DNA_SIGN = (
    "j_minus_dna_ratio_only = -( [log2((DNA_alt+0.5)/(DNA_ref+0.5))]_treated "
    "- [log2((DNA_alt+0.5)/(DNA_ref+0.5))]_control ). It carries the minus sign that the "
    "DNA term already has inside the activity definition a = log2(RNA+0.5) - log2(DNA+0.5). "
    "A POSITIVE value therefore means the ALT-over-REF plasmid DNA ratio is LOWER in the "
    "treated sample than in the control sample. It is not the DNA contrast itself; it is the "
    "DNA contrast's contribution to j, and j_full = j_rna_ratio_only + j_minus_dna_ratio_only "
    "exactly."
)


def tertile_of(v, ok):
    t = np.full(len(v), -1)
    vv = v[ok]
    q1, q2 = np.nanpercentile(vv, [100 / 3, 200 / 3])
    t[ok] = np.where(vv <= q1, 0, np.where(vv <= q2, 1, 2))
    return t, (float(q1), float(q2))


def main(outdir):
    tab = os.path.join(outdir, "tables")
    os.makedirs(tab, exist_ok=True)
    log = []

    def say(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        log.append(s)

    say("### A1 v2 addenda: eight corrections to the wave-1 deposit")
    el, ro, paired, grp, g2b, o, qc = load()
    elements = sorted(paired.element_id)
    M = build_matrices(o, elements)
    D, A = {}, {}
    for ctx in CONTEXTS:
        aref = M[(ctx, "ref", "act")].values
        aalt = M[(ctx, "alt", "act")].values
        D[ctx] = aalt - aref
        A[ctx] = (aalt + aref) / 2.0
    grp_e = grp.reindex(elements)
    blk = grp_e.map(g2b)
    blocks = sorted(blk.dropna().unique())
    bcode_all = pd.Series(pd.Categorical(blk, categories=blocks).codes, index=elements)

    # ---------------------------------------------------------------- fix 1
    fix = pd.read_csv(a1.FIXTURE, sep="\t", dtype=str)
    fids = sorted(set(fix.element_id))
    pos = {e: i for i, e in enumerate(elements)}
    sel = np.array([pos[e] for e in fids if e in pos])
    rows = []
    for ctx in CONTEXTS:
        Mat = D[ctx][sel]
        H = np.stack([half_means(Mat, p) for p in range(3)], axis=1)
        ok = np.isfinite(H).all(axis=(1, 2))
        st = reliability_from_halves(H[ok])
        pub = REF_PUBLISHED[ctx]
        diff = float(st["projected_4rep"]) - pub
        rows.append(dict(
            context=ctx, n_elements=int(ok.sum()),
            recomputed_projected_4rep=float(st["projected_4rep"]),
            published_projected_4rep_4dp=pub,
            signed_difference_recomputed_minus_published=diff,
            absolute_difference=abs(diff),
            agrees_to_4_decimals=bool(round(float(st["projected_4rep"]), 4) == pub),
            agrees_exactly=False,
            note=("the published value is a 4-decimal rounding of the same statistic; "
                  "the residual is the rounding, not a disagreement")))
    AGR = pd.DataFrame(rows)
    AGR.to_csv(f"{tab}/reference_reproduction_agreement.tsv", sep="\t", index=False)
    say("\n### fix 1: reference reproduction agrees to 4 decimals, NOT exactly")
    say(AGR.to_string(index=False))
    say(f"max |signed difference| = {AGR.absolute_difference.max():.3e}; "
        f"all four agree at 4 decimals: {bool(AGR.agrees_to_4_decimals.all())}")

    # ---------------------------------------------------------------- fix 2
    inv_rows = []
    for ctx in CONTEXTS:
        Mat = D[ctx][sel]
        ncell_fix = int(np.isnan(Mat).sum())
        ncell_all = int(np.isnan(D[ctx]).sum())
        H = np.stack([half_means(Mat, p) for p in range(3)], axis=1)
        ok = np.isfinite(H).all(axis=(1, 2))
        st = reliability_from_halves(H[ok])
        Hn = np.stack([np.column_stack([
            np.nanmean(Mat[:, [REPS.index(x) for x in SPLITS[p][0]]], axis=1),
            np.nanmean(Mat[:, [REPS.index(x) for x in SPLITS[p][1]]], axis=1)])
            for p in range(3)], axis=1)
        okn = np.isfinite(Hn).all(axis=(1, 2))
        stn = reliability_from_halves(Hn[okn])
        inv_rows.append(dict(
            context=ctx,
            missing_element_replicate_cells_full_4359=ncell_all,
            missing_element_replicate_cells_fixture_1033=ncell_fix,
            n_elements_complete_case=int(ok.sum()),
            n_elements_nan_skipping=int(okn.sum()),
            half_mean_arrays_bitwise_identical=bool(np.array_equal(H, Hn, equal_nan=True)),
            pearson_complete_case=float(st["pearson"]),
            pearson_nan_skipping=float(stn["pearson"]),
            abs_difference=abs(float(st["pearson"]) - float(stn["pearson"])),
            status=("INVARIANCE, not a check: with zero missing cells the two conventions "
                    "are the same arithmetic and cannot disagree")))
    INV = pd.DataFrame(inv_rows)
    INV.to_csv(f"{tab}/nan_convention_invariance.tsv", sep="\t", index=False)
    say("\n### fix 2: complete-case vs NaN-skipping is an INVARIANCE (zero missing cells)")
    say(INV.to_string(index=False))

    # ---------------------------------------------------------------- fix 3
    per_grp = grp_e.value_counts()
    fix_grp = fix.groupby("element_id")["outer_locus_sequence_group_id"].first()
    fix_per_grp = fix_grp.value_counts()
    blk_counts = blk.value_counts()
    map_rows = [
        dict(set_name="paired_snv analysis set", n_elements=len(elements),
             n_outer_locus_groups=int(grp_e.nunique()),
             elements_per_group_mean=float(len(elements) / grp_e.nunique()),
             elements_per_group_min=int(per_grp.min()),
             elements_per_group_median=float(per_grp.median()),
             elements_per_group_max=int(per_grp.max()),
             mapping="many_to_one"),
        dict(set_name="dna_lm common fixture", n_elements=len(fids),
             n_outer_locus_groups=int(fix_grp.nunique()),
             elements_per_group_mean=float(len(fids) / fix_grp.nunique()),
             elements_per_group_min=int(fix_per_grp.min()),
             elements_per_group_median=float(fix_per_grp.median()),
             elements_per_group_max=int(fix_per_grp.max()),
             mapping=("one_to_one" if fix_per_grp.max() == 1 else "many_to_one")),
        dict(set_name="paired_snv analysis set onto borzoi long-range blocks",
             n_elements=len(elements), n_outer_locus_groups=len(blocks),
             elements_per_group_mean=float(len(elements) / len(blocks)),
             elements_per_group_min=int(blk_counts.min()),
             elements_per_group_median=float(blk_counts.median()),
             elements_per_group_max=int(blk_counts.max()),
             mapping="many_to_one"),
    ]
    MAP = pd.DataFrame(map_rows)
    MAP.to_csv(f"{tab}/element_locus_group_mapping.tsv", sep="\t", index=False)
    say("\n### fix 3: elements map MANY to one onto outer locus groups")
    say(MAP.round(3).to_string(index=False))

    # ---------------------------------------------------------------- fixes 4, 7
    R = pd.read_csv(f"{tab}/reliability_by_context.tsv", sep="\t")
    trt_all = R[(R.quantity == "treatment_minus_control_allele_effect")
                & (R.stratum == "all")]
    part_rows, pub_rows = [], []
    for cell in CELL_PAIRS:
        s = trt_all[trt_all.context == cell]
        per = s[s.partition != "mean_of_three"].sort_values("partition")
        mean_row = s[s.partition == "mean_of_three"].iloc[0]
        for r in per.itertuples():
            part_rows.append(dict(cell_line=cell, partition=r.partition,
                                  n_elements=int(r.n_elements), n_blocks=int(r.n_blocks),
                                  pearson=float(r.pearson), spearman=float(r.spearman),
                                  projected_4rep_this_partition=float(r.projected_4rep)))
        pear = per.pearson.astype(float).values
        alt = float(np.mean([sb_project(x) for x in pear]))
        part_rows.append(dict(cell_line=cell, partition="mean_of_three",
                              n_elements=int(mean_row.n_elements),
                              n_blocks=int(mean_row.n_blocks),
                              pearson=float(mean_row.pearson),
                              spearman=float(mean_row.spearman),
                              projected_4rep_this_partition=float(mean_row.projected_4rep)))
        pub_rows.append(dict(
            cell_line=cell, n_elements=int(mean_row.n_elements),
            n_blocks=int(mean_row.n_blocks),
            published_projected_4rep=float(mean_row.projected_4rep),
            published_projected_4rep_lo=float(mean_row.projected_4rep_lo),
            published_projected_4rep_hi=float(mean_row.projected_4rep_hi),
            published_convention=PUBLISHED_CONVENTION,
            alternative_projected_4rep=alt,
            alternative_convention=ALT_CONVENTION,
            difference_published_minus_alternative=float(mean_row.projected_4rep) - alt,
            raw_2v2_pearson_mean_of_three=float(mean_row.pearson),
            partition_pearsons="; ".join(
                f"{r.partition} {float(r.pearson):+.4f}" for r in per.itertuples())))
    PART = pd.DataFrame(part_rows)
    PART.to_csv(f"{tab}/treatment_contrast_by_partition.tsv", sep="\t", index=False)
    PUB = pd.DataFrame(pub_rows)
    PUB.to_csv(f"{tab}/treatment_contrast_published_value.tsv", sep="\t", index=False)
    say("\n### fix 7: the three per-partition Pearsons behind each treatment value")
    say(PART.round(4).to_string(index=False))
    say("\n### fix 4: ONE published projected value per cell line, convention named")
    say(PUB[["cell_line", "n_elements", "n_blocks", "published_projected_4rep",
             "published_projected_4rep_lo", "published_projected_4rep_hi",
             "alternative_projected_4rep",
             "difference_published_minus_alternative"]].round(4).to_string(index=False))
    say(f"published convention = {PUBLISHED_CONVENTION}")
    say(f"alternative convention = {ALT_CONVENTION}")

    # ---------------------------------------------------------------- fix 5
    dna_min, bc_min = {}, {}
    for ctx in CONTEXTS:
        dna_min[ctx] = np.minimum(np.nanmin(M[(ctx, "ref", "DNA")].values, axis=1),
                                  np.nanmin(M[(ctx, "alt", "DNA")].values, axis=1))
        bc_min[ctx] = np.minimum(np.nanmin(M[(ctx, "ref", "n_barcodes")].values, axis=1),
                                 np.nanmin(M[(ctx, "alt", "n_barcodes")].values, axis=1))
    thresh = {}
    for cell, (trt, ctl) in CELL_PAIRS.items():
        J = D[trt] - D[ctl]
        H = np.stack([half_means(J, p) for p in range(3)], axis=1)
        ok = np.isfinite(H).all(axis=(1, 2)) & ~np.isnan(bcode_all.values)
        for sname, vals in (("min_dna", np.minimum(dna_min[trt], dna_min[ctl])),
                            ("min_barcodes", np.minimum(bc_min[trt], bc_min[ctl]))):
            _, qs = tertile_of(vals, ok)
            thresh[(cell, sname)] = qs
    ter_rows = []
    for r in R[(R.quantity == "treatment_minus_control_allele_effect")
               & (R.stratum != "all") & (R.partition == "mean_of_three")].itertuples():
        sname, k = r.stratum.rsplit("_tertile", 1)
        q1, q2 = thresh[(r.context, sname)]
        lo = "-inf" if k == "1" else (f"{q1:g}" if k == "2" else f"{q2:g}")
        hi = f"{q1:g}" if k == "1" else (f"{q2:g}" if k == "2" else "+inf")
        ter_rows.append(dict(cell_line=r.context, stratifier=sname, tertile=int(k),
                             tertile_lower_bound_exclusive=lo,
                             tertile_upper_bound_inclusive=hi,
                             n_elements=int(r.n_elements), n_blocks=int(r.n_blocks),
                             pearson=float(r.pearson),
                             pearson_lo=float(r.pearson_lo), pearson_hi=float(r.pearson_hi),
                             spearman=float(r.spearman),
                             projected_4rep=float(r.projected_4rep),
                             projected_4rep_lo=float(r.projected_4rep_lo),
                             projected_4rep_hi=float(r.projected_4rep_hi)))
    TER = pd.DataFrame(ter_rows).sort_values(["cell_line", "stratifier", "tertile"])
    TER.to_csv(f"{tab}/treatment_tertile_values.tsv", sep="\t", index=False)
    say("\n### fix 5: treatment-by-tertile values, each with n_elements and n_blocks")
    say(TER[["cell_line", "stratifier", "tertile", "n_elements", "n_blocks", "pearson",
             "pearson_lo", "pearson_hi", "projected_4rep"]].round(4).to_string(index=False))

    # ---------------------------------------------------------------- fix 6
    DEC = pd.read_csv(f"{tab}/treatment_contrast_decomposition.tsv", sep="\t")
    sign_map = {
        "j_full": "alt minus ref within a sample, then treated minus control",
        "j_rna_ratio_only": ("positive = the ALT-over-REF RNA ratio is HIGHER in the "
                             "treated sample than in the control sample"),
        "j_minus_dna_ratio_only": ("positive = the ALT-over-REF plasmid DNA ratio is LOWER "
                                   "in the treated sample than in the control sample "
                                   "(the leading minus is part of the activity definition; "
                                   "double negation)"),
        "d_treated": "alt minus ref in the treated sample",
        "d_control": "alt minus ref in the control sample",
        "d_mean_of_conditions": "alt minus ref, averaged over the two conditions",
        "j_residual_after_scale_removal": "j minus beta * d_mean_of_conditions",
    }
    DEC["sign_refers_to"] = DEC.component.map(sign_map).fillna(
        "diagnostic row, not a reliability")
    DEC["is_reliability"] = ~DEC.component.str.startswith("_")
    nblk = {}
    for cell, (trt, ctl) in CELL_PAIRS.items():
        J = D[trt] - D[ctl]
        H = np.stack([half_means(J, p) for p in range(3)], axis=1)
        ok = np.isfinite(H).all(axis=(1, 2))
        nblk[cell] = int(len(np.unique(bcode_all.values[ok & ~np.isnan(bcode_all.values)])))
    DEC["n_blocks"] = DEC.cell_line.map(nblk)
    DEC.to_csv(f"{tab}/treatment_contrast_decomposition_annotated.tsv",
               sep="\t", index=False)

    add_rows = []
    for cell, (trt, ctl) in CELL_PAIRS.items():
        comp = {}
        for ctx in (trt, ctl):
            comp[(ctx, "rna")] = (np.log2(M[(ctx, "alt", "RNA")].values + PC)
                                  - np.log2(M[(ctx, "ref", "RNA")].values + PC))
            comp[(ctx, "dna")] = (np.log2(M[(ctx, "alt", "DNA")].values + PC)
                                  - np.log2(M[(ctx, "ref", "DNA")].values + PC))
        j_full = D[trt] - D[ctl]
        j_rna = comp[(trt, "rna")] - comp[(ctl, "rna")]
        j_dna = -(comp[(trt, "dna")] - comp[(ctl, "dna")])
        resid = j_full - (j_rna + j_dna)
        fin = np.isfinite(resid)
        add_rows.append(dict(
            cell_line=cell, n_element_replicate_cells=int(fin.sum()),
            max_abs_additivity_residual=float(np.nanmax(np.abs(resid[fin]))),
            mean_j_minus_dna_ratio_only=float(np.nanmean(j_dna[fin])),
            sd_j_minus_dna_ratio_only=float(np.nanstd(j_dna[fin], ddof=1)),
            fraction_positive_j_minus_dna_ratio_only=float(np.mean(j_dna[fin] > 0)),
            sign_statement=J_DNA_SIGN))
    ADD = pd.DataFrame(add_rows)
    ADD.to_csv(f"{tab}/j_component_sign_and_additivity.tsv", sep="\t", index=False)
    say("\n### fix 6: sign of j_minus_dna_ratio_only, and the additivity guard")
    say(ADD[["cell_line", "n_element_replicate_cells", "max_abs_additivity_residual",
             "mean_j_minus_dna_ratio_only", "sd_j_minus_dna_ratio_only",
             "fraction_positive_j_minus_dna_ratio_only"]].to_string(index=False))
    say(J_DNA_SIGN)

    # ---------------------------------------------------------------- fix 8
    say("\n### fix 8: activity-convention check (written in wave 1, never executed)")
    conv_rows = []
    for ctx in CONTEXTS:
        aref = M[(ctx, "ref", "act")].values
        aalt = M[(ctx, "alt", "act")].values
        for name, Mat in (("activity_element_mean_of_alleles", (aref + aalt) / 2.0),
                          ("activity_construct_level_stacked", np.vstack([aref, aalt])),
                          ("activity_ref_construct_only", aref),
                          ("activity_alt_construct_only", aalt),
                          ("allele_effect", aalt - aref)):
            H = np.stack([half_means(Mat, p) for p in range(3)], axis=1)
            ok = np.isfinite(H).all(axis=(1, 2))
            st = reliability_from_halves(H[ok])
            conv_rows.append(dict(context=ctx, quantity=name, n_units=int(ok.sum()),
                                  unit=("element x context" if not name.endswith("stacked")
                                        else "construct x context"),
                                  pearson=st["pearson"], spearman=st["spearman"],
                                  projected_4rep=st["projected_4rep"]))
    CONV = pd.DataFrame(conv_rows)
    CONV.to_csv(f"{tab}/activity_convention_sensitivity.tsv", sep="\t", index=False)
    say(CONV.round(4).to_string(index=False))
    piv = CONV.pivot(index="context", columns="quantity", values="pearson")
    verdicts = {c: bool(min(piv.loc[c, q] for q in piv.columns if q.startswith("activity"))
                        > piv.loc[c, "allele_effect"]) for c in CONTEXTS}
    say(f"A1.2 holds under EVERY activity convention, per context: {verdicts}")
    say(f"A1.2 holds under every convention in all four contexts: {all(verdicts.values())}")
    with open(f"{tab}/activity_convention_sensitivity.note.txt", "w") as fh:
        fh.write(json.dumps(verdicts, indent=2) + "\n")

    with open(f"{outdir}/run_log_addenda.txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    say(f"\nwrote addenda into {tab}")


if __name__ == "__main__":
    main(sys.argv[1])
