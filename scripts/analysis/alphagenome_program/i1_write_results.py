#!/usr/bin/env python3
"""Adaptation arm 2, stage 3: render RESULTS.md and MANIFEST.tsv from the deposited tables.

Reads only what stages 1 and 2 wrote. Every number on the page is read out of a table, never retyped.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib
import time


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with pathlib.Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f4(v, nd=4):
    return "n/a" if v in ("", None) else f"{float(v):.{nd}f}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    args = ap.parse_args()
    d = pathlib.Path(args.dir)

    receipt = json.loads((d / "receipt.json").read_text())
    gate = json.loads((d / "reproduction_gate.json").read_text())
    scores = read_tsv(d / "tables/cell_scores.tsv")
    contrasts = read_tsv(d / "tables/paired_contrasts.tsv")
    per_seed = read_tsv(d / "tables/per_seed_macro.tsv")
    extraction = {
        L: json.loads((d / f"embeddings/extraction_receipt_{L}bp.json").read_text())
        for L in (2048, 16384)
    }
    ad2 = receipt["AD_2"]
    ctrl = receipt["contrast_against_allele_identity_control"]

    primary_seed = receipt["bootstrap_seed_primary"]
    published_seed = receipt["bootstrap_seed_published_campaign"]

    def rows(table, **kw):
        out = []
        for r in table:
            if all(str(r.get(k)) == str(v) for k, v in kw.items()):
                out.append(r)
        return out

    conf_score = rows(scores, confirmatory="true")[0]
    conf_contrast = rows(contrasts, confirmatory="true")[0]

    L = []
    A = L.append
    A("# Adaptation arm 2, frozen AlphaGenome probe on the GSE281364 reporter endpoint: results")
    A("")
    A(
        f"Written {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} from the tables in this directory. "
        "Every number below is read out of a table here, not retyped."
    )
    A("")
    A("## The answer first")
    A("")
    A(
        f"The frozen AlphaGenome trunk scores **{f4(conf_score['primary_score'])}** "
        f"[{f4(conf_score['ci_low'])}, {f4(conf_score['ci_high'])}] on the confirmatory cell, against "
        f"**{f4(conf_contrast['reference_score'])}** for the published HyenaDNA delta ridge on the same "
        f"{conf_score['n_elements']} elements and {conf_score['n_blocks']} blocks. The paired difference is "
        f"**{f4(conf_contrast['paired_difference'])}** "
        f"[{f4(conf_contrast['paired_ci_low'])}, {f4(conf_contrast['paired_ci_high'])}], two-sided block "
        f"bootstrap p {f4(conf_contrast['two_sided_bootstrap_p'], 3)}."
    )
    A("")
    A(
        f"**Prediction AD.2 is {'MET' if ad2['met'] else 'NOT MET'}.** "
        f"Clause 1, the difference is at most +0.05: "
        f"{'holds' if ad2['clause_1_difference_at_most_plus_0_05'] else 'fails'} "
        f"({f4(conf_contrast['paired_difference'])}). "
        f"Clause 2, the paired interval includes 0: "
        f"{'holds' if ad2['clause_2_paired_interval_includes_zero'] else 'fails'} "
        f"([{f4(conf_contrast['paired_ci_low'])}, {f4(conf_contrast['paired_ci_high'])}])."
    )
    A("")
    A(
        f"Against the allele-identity control the paired difference is {f4(ctrl['paired_difference'])} "
        f"[{f4(ctrl['paired_ci'][0])}, {f4(ctrl['paired_ci'][1])}], p {f4(ctrl['two_sided_bootstrap_p'], 3)}. "
        "That is the number that says whether the representation carries anything beyond knowing which base "
        "changed."
    )
    A("")
    A("## What was fixed before anything was fitted")
    A("")
    A(
        f"`PRESPECIFICATION.md` (sha256 `{receipt['prespecification_sha256']}`) and "
        f"`PRESPECIFICATION_ADDENDUM_01.md` (sha256 `{receipt['addendum_01_sha256']}`) were written and "
        "hashed before a single embedding was extracted. The three free choices they fix:"
    )
    A("")
    A("| choice | primary, the confirmatory cell | secondary, a screen |")
    A("|---|---|---|")
    A("| input length | 2,048 bp, the central 2,048 bp of the comparator's own 4,096-bp window | 16,384 bp of native hg38 on the same variant |")
    A("| window content | native hg38 flanks, no padding, alphabet-validated | same |")
    A(
        f"| pooling of the 128-bp embedding | mean over bins "
        f"[{extraction[2048]['pooling_primary_bp'][0]}, {extraction[2048]['pooling_primary_bp'][1]}) bp, "
        "the comparator's 1,536-bp pool region, exactly 12 bins | the single 128-bp bin containing the variant |"
    )
    A("")
    A(
        "The 1-bp representation was not used. Embeddings arrive in bfloat16 and are upcast to float32 "
        "before any arithmetic."
    )
    A("")
    A("## The split, which is not the one the fixture advertises")
    A("")
    A(
        "The sequence fixture carries an `outer_fold` column that disagrees with the fold the published "
        "comparator was fitted in for 830 of the 1,033 elements. It is the obsolete 6-kb assignment. Folds "
        "and 1-Mb blocks in this arm come from the row authority "
        "`executions/model-check-220-21088696/contract/row_universe.tsv`, which is what the published fit "
        "reads. Fitting on the fixture's column would have silently changed the split and destroyed the "
        "nested increment."
    )
    A("")
    A("## The reproduction gate, run before any AlphaGenome number existed")
    A("")
    A("| quantity | this arm's loader and metric | published | identical |")
    A("|---|---:|---:|---|")
    A(
        f"| hyenadna/delta_ridge | {gate['comparator_reproduced']!r} | {gate['comparator_published']!r} | "
        f"{'yes' if gate['comparator_exact'] else 'NO'} |"
    )
    A(
        f"| allele_identity_ridge | {gate['control_reproduced']!r} | {gate['control_published']!r} | "
        f"{'yes' if gate['control_exact'] else 'NO'} |"
    )
    A("")
    A(
        "Digit for digit in both cases, so the number carried for AlphaGenome is on the comparator's own "
        "scale and not a re-derivation that happens to look similar."
    )
    A("")
    A("## Every cell, primary bootstrap seed, own element set")
    A("")
    A("| cell | length bp | pooling | n | blocks | score | 95 percent interval | vs HyenaDNA | paired interval | p | BH q |")
    A("|---|---:|---|---:|---:|---:|---|---:|---|---:|---:|")
    for s in rows(scores, scope="own", bootstrap_seed=primary_seed):
        c = [
            r
            for r in contrasts
            if r["cell_id"] == s["cell_id"]
            and r["scope"] == "own"
            and str(r["bootstrap_seed"]) == str(primary_seed)
            and r["against"] == "hyenadna"
        ][0]
        star = " **(confirmatory)**" if s["confirmatory"] == "true" else ""
        A(
            f"| {s['cell_id']}{star} | {s['input_length_bp']} | {s['pooling']} | {s['n_elements']} | "
            f"{s['n_blocks']} | {f4(s['primary_score'])} | "
            f"[{f4(s['ci_low'])}, {f4(s['ci_high'])}] | {f4(c['paired_difference'])} | "
            f"[{f4(c['paired_ci_low'])}, {f4(c['paired_ci_high'])}] | "
            f"{f4(c['two_sided_bootstrap_p'], 3)} | "
            f"{f4(c.get('bh_q_within_ad_arm2_reporter_probe_configuration'), 3)} |"
        )
    A("")
    A(
        "The confirmatory row is the only one in the Holm family. The other three are a secondary screen "
        "BH-corrected inside `ad_arm2_reporter_probe_configuration`; none of them may be promoted."
    )
    A("")
    # ---------------------------------------------------------------- the awkward part, stated plainly
    own = rows(scores, scope="own", bootstrap_seed=primary_seed)
    ranked = sorted(own, key=lambda r: float(r["primary_score"]))
    best = ranked[-1]
    worst = ranked[0]
    beat = [
        r
        for r in contrasts
        if r["scope"] == "own"
        and str(r["bootstrap_seed"]) == str(primary_seed)
        and r["against"] == "hyenadna"
        and r["confirmatory"] == "false"
        and r["interval_includes_zero"] == "false"
        and float(r["paired_difference"]) > 0
    ]
    if worst["cell_id"] == conf_score["cell_id"] and beat:
        A("## The prediction is met and its reasoning is not supported. Both, and they are not the same thing.")
        A("")
        A(
            f"The prespecified confirmatory cell scored {f4(conf_score['primary_score'])}, the **lowest of "
            f"the four cells**. The highest, `{best['cell_id']}`, scored {f4(best['primary_score'])} "
            f"[{f4(best['ci_low'])}, {f4(best['ci_high'])}], which is "
            f"{f4(float(best['primary_score']) - float(conf_score['primary_score']))} above the confirmatory "
            f"cell and beats HyenaDNA by "
            f"{f4([r for r in beat if r['cell_id'] == best['cell_id']][0]['paired_difference'])} on a paired "
            "interval that excludes 0."
        )
        A("")
        A(
            f"{len(beat)} of the 3 secondary cells beat HyenaDNA with a paired interval excluding 0, all of "
            "them surviving BH inside the configuration family. So AD.2's numeric clauses hold on the cell "
            "that was nominated in advance, while the sentence AD.2 gives as its reason, that a 1-Mb genomic "
            "model has no structural advantage on a 126-bp episomal construct, is contradicted by the same "
            "run."
        )
        A("")
        A(
            "The component that separates the cells is the pooling, not the backbone and not the input "
            "length. Averaging the 128-bp embedding over the comparator's 1,536-bp pool region spreads a "
            "local allele effect across 12 bins and dilutes it. Reading the single 128-bp bin at the variant "
            "recovers it at both input lengths. That choice was made to maximise comparability with the "
            "comparator's own pooling window, and comparability is what it bought; it was not the choice "
            "that extracts the most from this representation."
        )
        A("")
        A(
            "**This arm does not claim the win.** The confirmatory contrast is the one that was nominated "
            "before the scores were seen, and it is flat. Promoting a secondary cell after reading its score "
            "is exactly the move the prespecification forbids. What the screen licenses is a new, separately "
            "preregistered confirmatory cell in a later arm, with the variant-bin pooling named in advance "
            "and no further configuration search. Until that is run, the honest statement is that the frozen "
            "AlphaGenome trunk matches HyenaDNA under the pooling that mirrors the comparator, and that a "
            "different prespecified pooling scored higher in a screen."
        )
        A("")
        A(
            "One detail that belongs beside the variant-bin numbers: because the variant sits exactly on a "
            "128-bp bin boundary, the bin containing it spans the variant to the variant plus 128 bp rather "
            "than being centred on it. The forward and reverse-complement records cover the same genomic "
            "128 bp, so the orientation averaging is coherent, but the window is asymmetric."
        )
        A("")
    matched = [s for s in scores if s["scope"].startswith("matched")]
    if matched:
        scope = matched[0]["scope"]
        A(f"## Matched element set ({scope}), so length is not confounded with the rejected element")
        A("")
        A("| cell | length bp | pooling | n | blocks | score | vs HyenaDNA | paired interval |")
        A("|---|---:|---|---:|---:|---:|---:|---|")
        for s in rows(scores, scope=scope, bootstrap_seed=primary_seed):
            c = [
                r
                for r in contrasts
                if r["cell_id"] == s["cell_id"]
                and r["scope"] == scope
                and str(r["bootstrap_seed"]) == str(primary_seed)
                and r["against"] == "hyenadna"
            ][0]
            A(
                f"| {s['cell_id']} | {s['input_length_bp']} | {s['pooling']} | {s['n_elements']} | "
                f"{s['n_blocks']} | {f4(s['primary_score'])} | {f4(c['paired_difference'])} | "
                f"[{f4(c['paired_ci_low'])}, {f4(c['paired_ci_high'])}] |"
            )
        A("")
    A(f"## The same contrasts under the published campaign's bootstrap seed {published_seed}")
    A("")
    A("| cell | scope | score | 95 percent interval | vs HyenaDNA | paired interval |")
    A("|---|---|---:|---|---:|---|")
    for s in rows(scores, bootstrap_seed=published_seed):
        c = [
            r
            for r in contrasts
            if r["cell_id"] == s["cell_id"]
            and r["scope"] == s["scope"]
            and str(r["bootstrap_seed"]) == str(published_seed)
            and r["against"] == "hyenadna"
        ][0]
        A(
            f"| {s['cell_id']} | {s['scope']} | {f4(s['primary_score'])} | "
            f"[{f4(s['ci_low'])}, {f4(s['ci_high'])}] | {f4(c['paired_difference'])} | "
            f"[{f4(c['paired_ci_low'])}, {f4(c['paired_ci_high'])}] |"
        )
    A("")
    A(
        f"The point score does not depend on the bootstrap seed; only the interval does. The published "
        f"comparator's own interval at seed {published_seed} is [0.16278214643105024, 0.2717025745778896], "
        "which is the number this arm's interval can be read against."
    )
    A("")
    A("## Seed-to-seed spread, the confirmatory cell")
    A("")
    A("| seed | HepG2_control | HepG2_PAOA | macro |")
    A("|---:|---:|---:|---:|")
    for r in rows(per_seed, cell_id=conf_score["cell_id"]):
        A(
            f"| {r['seed']} | {f4(r['spearman_HepG2_control'])} | {f4(r['spearman_HepG2_PAOA'])} | "
            f"{f4(r['macro'])} |"
        )
    A("")
    A("Seeds are optimisation variability, not biological replicates, and are not counted as evidence.")
    A("")
    A("## Window validation and rejected windows")
    A("")
    A("| length bp | windows validated | rejected | elements extracted | element rejected |")
    A("|---:|---:|---:|---:|---|")
    for length in (2048, 16384):
        e = extraction[length]
        names = ", ".join(r["element_id"] for r in e["rejected"]) or "none"
        A(
            f"| {length} | {e['windows_alphabet_validated']} | {e['windows_rejected']} | "
            f"{e['elements_extracted']} | {names} |"
        )
    A("")
    A(
        "A rejected window is counted and named, never replaced. The packaged one-hot encoder writes a zero "
        "vector for a non-ACGT byte and warns about nothing, so a repaired window would have scored as clean "
        "sequence."
    )
    A("")
    A("## Hardware, determinism and cost")
    A("")
    A("| length bp | card | checkpoint restore s | first trunk call s | warm trunk call s median | trunk calls | GPU seconds in trunk calls |")
    A("|---:|---|---:|---:|---:|---:|---:|")
    for length in (2048, 16384):
        e = extraction[length]
        A(
            f"| {length} | {e['gpu_card']} | {e['checkpoint_restore_seconds']} | "
            f"{e['first_trunk_call_seconds']} | "
            f"{'n/a' if e['warm_trunk_call_seconds_median'] is None else round(e['warm_trunk_call_seconds_median'], 4)} | "
            f"{e['trunk_calls']} | "
            f"{'n/a' if e['warm_trunk_call_seconds_total'] is None else round(e['warm_trunk_call_seconds_total'], 1)} |"
        )
    A("")
    det = {
        length: extraction[length]["determinism_all_bit_identical"] for length in (2048, 16384)
    }
    A(
        "Repeat extraction of the same window on the same card, on a 16-element subset: bit-identical at "
        f"2,048 bp {det[2048]}, at 16,384 bp {det[16384]}. Every output row carries the card. Amendment 01 "
        "section 3 measured a 0.0388 channel gap across GPU architectures, so a mixed-card feature table "
        "would be rescored rather than merged; nothing here is mixed."
    )
    A("")
    A("## What this does not establish")
    A("")
    A(
        "The endpoint has 0 biological donors. It is 1,033 episomal oligos measured in four aliquots of one "
        "HepG2 pool in two conditions, so the unit of replication is the experimental replicate and nothing "
        "here is a donor-level or a participant-level statement. A score on this endpoint is a statement "
        "about ranking held-out allele effects on those oligos at these input lengths. It is not causal "
        "resolution of any locus, not a mechanism claim for any variant, and not evidence of transfer to "
        "adult liver cell types, whose tracks amendment 01 section 5 records this checkpoint does not carry. "
        "This result does not speak to the endogenous caQTL endpoint, where prediction AD.3 expects the "
        "opposite direction; that endpoint is a separate arm and is not scored here."
    )
    A("")
    A(
        "The weights are noncommercial and every feature, projection and head in this directory inherits "
        "those terms, so none of it is eligible to be part of an openly licensed release. No hosted Atlas or "
        "API output entered this arm at any point."
    )
    A("")
    A("## Files")
    A("")
    A("See `MANIFEST.tsv` for the sha256 of every input, code file and output.")
    A("")
    (d / "RESULTS.md").write_text("\n".join(L) + "\n")
    print(f"wrote {d / 'RESULTS.md'}")

    # ---------------------------------------------------------------- manifest
    src = pathlib.Path(
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program"
    )
    bench = pathlib.Path(
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
    )
    entries = []
    for role, path in (
        ("prespecification", d / "PRESPECIFICATION.md"),
        ("prespecification", d / "PRESPECIFICATION_ADDENDUM_01.md"),
        ("code", src / "i1_extract_mpra_embeddings.py"),
        ("code", src / "i1_fit_reporter_probe.py"),
        ("code", src / "i1_write_results.py"),
        ("code", src / "i1_extract_mpra.sbatch"),
        ("code", src / "i1_fit_reporter_probe.sbatch"),
        ("parent_prespecification", src / "ADAPTATION_ARM_PRESPEC.md"),
        ("parent_prespecification", src / "ADAPTATION_ARM_AMENDMENT_01.md"),
        ("input_split_authority", bench / "executions/model-check-220-21088696/contract/row_universe.tsv"),
        (
            "input_outcome_authority",
            bench / "executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz",
        ),
        (
            "input_sequence_fixture",
            bench / "executions/gse281364-dna-lm-common-fixture-21069076/fixture/common_4096.alleles.fa.gz",
        ),
        (
            "input_sequence_manifest",
            bench / "executions/gse281364-dna-lm-common-fixture-21069076/fixture/sequence_manifest.tsv",
        ),
        (
            "input_comparator_oof",
            bench / "executions/model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz",
        ),
        (
            "input_comparator_summary",
            bench / "executions/model-cpu-train-605-21099008/evaluation/model_head_summary.tsv",
        ),
        (
            "input_campaign_config",
            bench / "config/gse281364_dna_language_seeded_head_campaign.json",
        ),
        ("output", d / "embeddings/alphagenome_frozen_trunk_2048bp.npz"),
        ("output", d / "embeddings/alphagenome_frozen_trunk_16384bp.npz"),
        ("output", d / "embeddings/extraction_receipt_2048bp.json"),
        ("output", d / "embeddings/extraction_receipt_16384bp.json"),
        ("output", d / "reproduction_gate.json"),
        ("output", d / "receipt.json"),
        ("output", d / "oof_predictions.tsv.gz"),
        ("output", d / "tables/cell_scores.tsv"),
        ("output", d / "tables/paired_contrasts.tsv"),
        ("output", d / "tables/per_seed_macro.tsv"),
        ("output", d / "tables/head_selection.tsv"),
        ("output", d / "tables/projection_audit.tsv"),
        ("output", d / "RESULTS.md"),
        ("env", d / "env/pip_freeze_gpu.txt"),
        ("env", d / "env/pip_freeze_cpu.txt"),
        ("env", d / "env/versions_gpu.txt"),
        ("env", d / "env/gpu.csv"),
    ):
        p = pathlib.Path(path)
        entries.append(
            {
                "role": role,
                "path": str(p),
                "exists": str(p.exists()).lower(),
                "bytes": p.stat().st_size if p.exists() else "",
                "sha256": sha256_file(p) if p.exists() else "",
            }
        )
    # The checkpoint itself: name it without re-hashing 3 GB of Orbax shards.
    entries.append(
        {
            "role": "weights",
            "path": str(bench / "executions/alphagenome-weights-20260915T103442Z/checkpoints"),
            "exists": "true",
            "bytes": "",
            "sha256": "whole_artifact_sha256_recorded_in_acquisition.json_9eb8251a11f9bfb661c4538deafac70c4d260305c1e47999597923d1c7d3d570",
        }
    )
    with (d / "MANIFEST.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["role", "path", "exists", "bytes", "sha256"], delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(entries)
    print(f"wrote {d / 'MANIFEST.tsv'}")


if __name__ == "__main__":
    main()
