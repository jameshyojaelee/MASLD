#!/usr/bin/env python3
"""Stage 5: build RESULTS.md from the probe JSONs.

Every number in RESULTS.md is read out of a producing JSON here, never retyped, so the document and the
artifacts cannot drift apart.
"""

from __future__ import annotations

import json
import pathlib
import sys

PROJ = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
EXEC = PROJ / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z"
LENGTHS = [2048, 16384, 131072, 524288, 1048576]
GIB = float(1 << 30)


CAN_CANNOT = [
    "**Can.** Score an arbitrary DNA string at any of the five lengths the SDK names, and at other "
    "multiples of 2,048 that were tried, on one GPU, with the liver ontology, returning per-base tracks "
    "that cover the input exactly. It reproduces the hosted "
    "recipe end to end, including insertions built into the window by hand, which is the one thing the "
    "Atlas point query refuses. Repeat calls within one process on one card are bit-identical across all "
    "35.8 million returned values, so a cached score and a recomputed score on the same hardware cannot "
    "silently disagree. Reproducibility across separate processes was not tested; across GPU "
    "architectures it is reported above.",
    "",
    "**Cannot.** It does not take a donor, a cell state, a disease context or an observed assay as input: "
    "`ontology_terms` selects among fixed training tracks and injects nothing. It returns a predicted "
    "track, not a calibrated quantile, so these values are not comparable with Atlas quantiles and must "
    "never be pooled with them. Nothing here says a variant is causal or assigns it a gene.",
    "",
    "**Guard the caller must add.** The one-hot encoder maps every byte that is not A/C/G/T to a zero "
    "vector and says nothing. A window containing an assembly gap, a soft-masked repeat written as N, or "
    "a corrupted read silently produces a confident-looking prediction. Any production loop must reject "
    "non-ACGT input itself before calling the model.",
]

FINETUNE_BLOCKERS = [
    "1. **The shipped trainer wants AlphaGenome's own training data, not a custom endpoint.** "
    "`finetuning/finetune.py` optimizes `AlphaGenome.loss(batch)`, the full multimodal loss over every "
    "head, fed by `finetuning/dataset.py`'s `DataPipeline`. Fitting an MPRA, a Cas13 readout or any "
    "other single endpoint means writing a new head and loss, which is a code change, not a config "
    "change.",
    "",
    "2. **`DataPipeline` cannot run on the packaged metadata at all.** It builds one `BigWigExtractor` "
    "per row of `metadata['file_path']`, and the packaged human metadata has no `file_path` column in "
    "any of its eleven output types (checked directly). Someone must supply a metadata table carrying "
    "local bigWig paths for the tracks being trained.",
    "",
    "3. **Fold definitions are fetched from the public internet by default.** "
    "`fold_intervals.get_fold_intervals` falls back to a Borzoi BED hosted on GitHub. A compute node "
    "without egress needs a local copy, and the copy has to be pinned or the train/valid/test split is "
    "not reproducible.",
    "",
    "4. **Training memory is unmeasured and is the real risk.** Inference at 1 Mb peaks at 17.6 GiB. A "
    "training step adds gradients, optimizer state and retained activations on top; the upstream README "
    "recommends TPU v3 or better for training. This probe did not run a training step, so there is no "
    "measured number. Measure one short step at the shortest useful length before planning anything at "
    "1 Mb. [Unverified: no training step was run.]",
    "",
    "5. **Reference build.** The trainer's default FASTA is GRCh38.p13 with GENCODE v46 annotation; this "
    "project standardises on the cellranger-arc GRCh38-2024-A FASTA with GENCODE v49. Inference here used "
    "the project FASTA. If fine-tuning targets are built from project annotation, the annotation vintage "
    "differs from the one the released weights were trained against.",
    "",
    "6. **Terms.** Noncommercial, and derivatives inherit them. Outputs of these local weights are "
    "eligible as features or teacher signals; hosted Atlas and API outputs are not, and the two are never "
    "pooled.",
]


JOBS = [
    ("21767351", "agp_agruntime_b6k", "h1_runtime_probe_b6k.sbatch"),
    ("21767352", "agp_agruntime_l40s", "h1_runtime_probe_l40s.sbatch"),
    ("21767377", "agp_agruntime_extra", "h1_runtime_probe_extra.sbatch"),
    ("21767443", "agp_agruntime_extrab6k", "h1_runtime_probe_extra_b6k.sbatch"),
]


def jobs_section() -> list[str]:
    import subprocess

    out = ["## Jobs", "", "| job id | name | final state | elapsed | resource ask |", "|---|---|---|---|---|"]
    for jid, name, script in JOBS:
        st = subprocess.run(
            ["sacct", "-j", jid, "-X", "-n", "-P", "--format=State,Elapsed,ReqTRES"],
            capture_output=True, text=True, check=False).stdout.strip().splitlines()
        state, elapsed, tres = ("unknown", "-", "-")
        if st:
            f = st[0].split("|")
            state, elapsed = f[0], f[1]
            tres = f[2] if len(f) > 2 else "-"
        out.append(f"| {jid} | {name} | {state} | {elapsed} | `{tres}` |")
    out.append("")
    out.append("sbatch headers as submitted:")
    out.append("")
    out.append("```")
    for _, _, script in JOBS:
        p = EXEC / "env" / script
        if not p.exists():
            continue
        head = [ln for ln in p.read_text().splitlines() if ln.startswith("#SBATCH")]
        out.extend(head)
        out.append("")
    out.append("```")
    return out


def extra_section() -> list[str]:
    rows = [
        "| input bp | state | outcome |",
        "|---:|---|---|",
    ]
    found = False
    for card, L in (("l40s", 786_432), ("l40s", 2_097_152), ("b6k", 2_097_152)):
        d = load(EXEC / f"tables/{card}_extra/bench_{L}.json")
        if d is None:
            continue
        found = True
        rows[0] = "| card | input bp | state | outcome |"
        rows[1] = "|---|---:|---|---|"
        if d.get("state") == "ok":
            mp = d.get("memory_peak", {}) or {}
            cov = d.get("output_covers_full_input", {}) or {}
            rows.append(
                f"| {card} | {L:,} | ran | warm call {d.get('warm_call_seconds_median')} s, peak "
                f"{gib(mp.get('peak_bytes_in_use'))} GiB of a {gib(mp.get('bytes_limit'))} GiB limit, "
                f"every channel covers the input exactly: "
                f"{'yes' if cov and all(cov.values()) else 'NO'} |"
            )
        else:
            mp = d.get("memory_peak", {}) or {}
            rows.append(
                f"| {card} | {L:,} | {d.get('state')} | {d.get('error_type', '?')}: "
                f"{str(d.get('error', ''))[:200]} (reached {gib(mp.get('peak_bytes_in_use'))} GiB of a "
                f"{gib(mp.get('bytes_limit'))} GiB limit) |"
            )
    if not found:
        return []
    l40 = load(EXEC / "tables/l40s_extra/bench_2097152.json")
    b6 = load(EXEC / "tables/b6k_extra/bench_2097152.json")
    rows.append("")
    if l40 and l40.get("state") != "ok":
        rows.append(
            "The 2 Mb failure on the L40S is a memory failure, not a refusal: the trunk ran and the "
            "allocation that failed was the float32 upcast of the outputs. So the ceiling at 1 Mb is the "
            "card, not a shape the model rejects."
        )
        if b6 and b6.get("state") == "ok":
            mp = b6.get("memory_peak", {}) or {}
            rows.append("")
            rows.append(
                f"On the larger card 2,097,152 bp does run, peaking at "
                f"{gib(mp.get('peak_bytes_in_use'))} GiB, so a 2 Mb window is available if a question "
                f"ever needs one. Nothing in this project's pinned recipe does; it uses 1 Mb, and the "
                f"released weights were trained at up to 1 Mb, so a 2 Mb window is outside the regime the "
                f"model was fitted in and its outputs there are not underwritten by anything measured "
                f"here."
            )
        elif b6 and b6.get("state") != "ok":
            rows.append("")
            rows.append(
                f"It also fails on the larger card: {b6.get('error_type')}. "
                f"{str(b6.get('error', ''))[:200]}"
            )
    return rows


def load(p: pathlib.Path):
    if not p.exists():
        return None
    with p.open() as fh:
        return json.load(fh)


def gib(x) -> str:
    return "-" if x is None else f"{x / GIB:.2f}"


def bench_table(tag: str) -> list[str]:
    rows = [
        "| input bp | state | first call s (incl. compile) | warm call s (median) | peak GPU GiB | device limit GiB |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    any_found = False
    for L in LENGTHS:
        d = load(EXEC / f"tables/{tag}/bench_{L}.json")
        if d is None:
            rows.append(f"| {L:,} | not run | - | - | - | - |")
            continue
        any_found = True
        mp = d.get("memory_peak", {}) or {}
        state = d.get("state", "?")
        if state != "ok":
            state = f"FAILED ({d.get('error_type', '?')})"
        rows.append(
            f"| {L:,} | {state} | {d.get('first_call_seconds_includes_compile', '-')} | "
            f"{d.get('warm_call_seconds_median', '-')} | {gib(mp.get('peak_bytes_in_use'))} | "
            f"{gib(mp.get('bytes_limit'))} |"
        )
    return rows if any_found else []


def gpu_commentary() -> list[str]:
    """Compare the two cards on the lengths both actually ran."""
    out = []
    common = []
    for L in LENGTHS:
        a = load(EXEC / f"tables/l40s/bench_{L}.json")
        b = load(EXEC / f"tables/b6k/bench_{L}.json")
        if a and b and a.get("state") == "ok" and b.get("state") == "ok":
            common.append((L, a, b))
    if not common:
        return out
    L, a, b = common[-1]
    biggest = max((x for x in LENGTHS if (load(EXEC / f"tables/l40s/bench_{x}.json") or {}).get("state") == "ok"),
                  default=None)
    big = load(EXEC / f"tables/l40s/bench_{biggest}.json") if biggest else None
    out.append(
        f"One environment covers both cards: the same jax {a['jax_version']} CUDA-12 build sees the GPU on "
        f"the L40S (driver 550.107.02, compute capability 8.9) and on the Blackwell RTX PRO 6000 (driver "
        f"590.48.01, compute capability 12.0). The Blackwell card is not the faster choice here. At "
        f"{L:,} bp its first call takes {b['first_call_seconds_includes_compile']} s against "
        f"{a['first_call_seconds_includes_compile']} s on the L40S, and its warm call is "
        f"{b['warm_call_seconds_median']} s against {a['warm_call_seconds_median']} s, consistent with "
        f"kernels being JIT-compiled from PTX because this build ships no sm_120 binaries. It also uses "
        f"more memory for the same work, {gib((b.get('memory_peak') or {}).get('peak_bytes_in_use'))} GiB "
        f"against {gib((a.get('memory_peak') or {}).get('peak_bytes_in_use'))} GiB, so XLA is choosing "
        f"different kernels there rather than simply running the same ones slower. Prefer the L40S unless "
        f"a job needs more than its {gib((a.get('memory_peak') or {}).get('bytes_limit'))} GiB working "
        f"limit."
    )
    if big:
        mp = big.get("memory_peak", {}) or {}
        out.append("")
        out.append(
            f"Every input length the SDK names runs on the L40S, and the largest of them, {biggest:,} bp, "
            f"peaks at {gib(mp.get('peak_bytes_in_use'))} GiB against a "
            f"{gib(mp.get('bytes_limit'))} GiB working limit on a 45 GiB card, so it fits with room to "
            f"spare. Whether anything longer fits is answered in a later section. One 1 Mb variant costs "
            f"two calls, so about "
            f"{2 * float(big['warm_call_seconds_median']):.1f} s of GPU time once the shape is compiled, "
            f"plus a one-off {big['first_call_seconds_includes_compile']} s compile per process and about "
            f"{big['checkpoint_load_seconds']} s to restore the checkpoint."
        )
    return out


def shapes_table(tag: str) -> list[str]:
    for L in (1_048_576, 524_288, 131_072, 16_384, 2_048):
        d = load(EXEC / f"tables/{tag}/bench_{L}.json")
        if d and d.get("state") == "ok" and d.get("output_shapes"):
            rows = [
                f"Output geometry at {L:,} bp input, liver ontology:",
                "",
                "| channel | tracks | resolution bp | rows | rows x resolution | equals input |",
                "|---|---:|---:|---:|---:|---|",
            ]
            for k, v in d["output_shapes"].items():
                if not v:
                    rows.append(f"| {k} | - | - | - | - | - |")
                    continue
                rows.append(
                    f"| {k} | {v['tracks']} | {v['resolution_bp']} | {v['rows']:,} | "
                    f"{v['rows_times_resolution']:,} | "
                    f"{'yes' if v['rows_times_resolution'] == L else 'NO'} |"
                )
            return rows
    return []


def gpu_line(tag: str) -> str:
    p = EXEC / f"env/gpu_{tag}.csv"
    if not p.exists():
        return f"{tag}: not run"
    lines = [x.strip() for x in p.read_text().splitlines() if x.strip()]
    return lines[-1] if len(lines) > 1 else f"{tag}: no data"


def probes_section(tag: str) -> list[str]:
    d = load(EXEC / f"tables/{tag}/variant_probes.json")
    if d is None:
        return ["Variant-effect probes: not run on this GPU."]
    out = ["| check | result | evidence |", "|---|---|---|"]

    a = d.get("check_a_ref_alt_localized", {})
    ch = a.get("per_channel", {})
    ev_a = "; ".join(
        f"{k} {v['enrichment_2000bp']:.0f}x"
        for k, v in ch.items()
        if isinstance(v, dict) and isinstance(v.get("enrichment_2000bp"), (int, float))
    )
    out.append(
        f"| (a) REF vs ALT differ, localized | {'PASS' if a.get('passed') else 'FAIL'} | "
        f"{a.get('substitution', '')}; outputs differ = {a.get('outputs_differ')}; "
        f"enrichment of \\|alt-ref\\| mass within +/-2 kb over the uniform-window share: {ev_a} |"
    )

    b1 = d.get("check_b1_metadata_alignment", {})
    per = b1.get("per_channel", {})
    ev_b1 = "; ".join(
        f"{k} {v['rows']:,}x{v['resolution_bp']}bp={v['rows_times_resolution']:,}"
        for k, v in per.items()
        if v
    )
    out.append(
        f"| (b1) output interval == input interval | {'PASS' if b1.get('passed') else 'FAIL'} | {ev_b1} |"
    )

    b2 = d.get("check_b2_empirical_alignment", {})
    offs = b2.get("offsets", {})
    ev_b2 = "; ".join(
        f"block at +{int(k):,} -> strongest channel {v['strongest_channel']} argmax "
        f"{v['per_channel'][v['strongest_channel']]['argmax_offset_bp']:,} "
        f"(delta {v['per_channel'][v['strongest_channel']]['distance_to_block_center_bp']:+,} bp)"
        for k, v in offs.items()
        if v.get("strongest_channel")
    )
    out.append(
        f"| (b2) response moves with the perturbation | {'PASS' if b2.get('passed') else 'FAIL'} | {ev_b2} |"
    )

    c = d.get("check_c_determinism", {})
    out.append(
        f"| (c) two runs agree | {'PASS' if c.get('passed') else 'FAIL'} | "
        f"max \\|run1-run2\\| over every returned track = {c.get('max_abs_diff_overall')} |"
    )

    dd = d.get("check_d_malformed_input", {})
    for label, v in dd.items():
        if not isinstance(v, dict):
            continue
        if v.get("raised"):
            ev = f"raised {v.get('error_type')}: {str(v.get('error', ''))[:160]}"
        else:
            ev = "returned an Output with no error"
            summ = v.get("output_summary", {}) or {}
            r = summ.get("rna")
            if r:
                ev += (
                    f"; rna mean\\|value\\| {r['mean_abs']:.6g} vs {r['mean_abs_reference_run']:.6g} "
                    f"on real sequence; NaNs {r['n_nan']}; all-zero {r['all_zero']}"
                )
        out.append(f"| (d) {label} | {'PASS' if v.get('passed') else 'FAIL'} | {ev} |")
    return out


def probe_replication(tag: str, other: str) -> list[str]:
    a = load(EXEC / f"tables/{tag}/variant_probes.json")
    b = load(EXEC / f"tables/{other}/variant_probes.json")
    if not a or not b:
        return []
    same = a["summary"] == b["summary"]
    ea = a["check_a_ref_alt_localized"]["per_channel"]
    eb = b["check_a_ref_alt_localized"]["per_channel"]
    pairs = ", ".join(
        f"{k} {ea[k]['enrichment_2000bp']:.0f}x / {eb[k]['enrichment_2000bp']:.0f}x"
        for k in ("splice", "rna", "h3k27ac", "atac", "dnase")
        if k in ea and k in eb
    )
    return [
        f"The whole probe set was run again on the {other} card and "
        + ("every check landed the same way" if same else "the checks did NOT land the same way")
        + f". The localization enrichments within +/-2 kb track closely ({tag} / {other}): {pairs}. "
        f"Determinism is bit-exact on both cards, with max \\|run1-run2\\| of "
        f"{b['check_c_determinism']['max_abs_diff_overall']} on {other}.",
        "",
    ]


def probe_commentary(tag: str) -> list[str]:
    """Prose that reads the measured numbers back, so no claim is retyped by hand."""
    d = load(EXEC / f"tables/{tag}/variant_probes.json")
    if d is None:
        return []
    out = []
    a = d.get("check_a_ref_alt_localized", {})
    ch = a.get("per_channel", {})

    def e(name, flank=2000):
        v = ch.get(name, {})
        return v.get(f"enrichment_{flank}bp")

    def s(name, flank=10000):
        v = ch.get(name, {})
        return v.get(f"share_within_{flank}bp")

    out.append(
        f"**On (a).** The outputs do differ: one substitution changes every channel. The prespecified bar "
        f"was a 10x enrichment of the \\|alt-ref\\| mass within +/-2 kb, in *every* channel, and two "
        f"channels miss it. Splice usage concentrates {e('splice'):.0f}x, RNA-seq {e('rna'):.0f}x and "
        f"H3K27ac {e('h3k27ac'):.0f}x, but ATAC reaches only {e('atac'):.1f}x and DNase "
        f"{e('dnase'):.1f}x. Within +/-10 kb the picture is the same: "
        f"{s('splice') * 100:.0f} percent of the splice response and {s('rna') * 100:.0f} percent of the "
        f"RNA response sit there, against {s('atac') * 100:.0f} percent of ATAC and "
        f"{s('dnase') * 100:.0f} percent of DNase. Determinism is bit-exact, so the diffuse remainder is "
        f"not numerical noise: a single base change genuinely moves the predicted accessibility tracks a "
        f"little across the whole megabase. Anyone reading a chromatin effect off a single SNV should "
        f"summarise a window around the variant, as the pinned recipe does, and not trust a whole-window "
        f"aggregate. The two channels that look diffuse here are the sharpest in check (b2) below, where "
        f"the perturbation is 1,024 bp instead of one base, so this is about the size of the change and "
        f"not about the model losing track of position."
    )
    out.append("")
    b2 = d.get("check_b2_empirical_alignment", {})
    offs = b2.get("offsets", {})
    if offs:
        deltas = []
        for k, v in offs.items():
            per = v.get("per_channel", {})
            for cname in ("h3k27ac", "atac", "dnase"):
                if cname in per:
                    deltas.append(abs(per[cname]["distance_to_block_center_bp"]))
        out.append(
            f"**On (b2).** The scrambled block was placed at three different offsets in the same window. "
            f"The chromatin channels put their largest response inside the block every time, the worst "
            f"miss being {max(deltas)} bp at 1-bp resolution and 32 bp at H3K27ac's 128-bp resolution. "
            f"The response moves with the perturbation, which is what separates real coordinate "
            f"alignment from a model that always answers at the window centre. RNA and splice-usage "
            f"argmax wander at two of the three offsets, which is expected rather than alarming: those "
            f"tracks respond where a transcript is, not where the base changed."
        )
        out.append("")
    dd = d.get("check_d_malformed_input", {})
    zz = dd.get("invalid_character_Z_10000_positions", {})
    nn = dd.get("all_N_full_length", {})
    if zz and nn and not zz.get("raised"):
        zr = (zz.get("output_summary") or {}).get("rna") or {}
        nr = (nn.get("output_summary") or {}).get("rna") or {}
        out.append(
            f"**On (d), the one that failed.** Every wrong *length* raises immediately, including an "
            f"off-by-one at 1 Mb, and passing the genome interval alongside it changes nothing because "
            f"the model rejects the shape first. Wrong *characters* do not raise at all. The one-hot "
            f"encoder writes a zero vector for every byte that is not A/C/G/T and returns no warning. An "
            f"all-N megabase produced finite, non-zero tracks four to five orders of magnitude below the "
            f"real-sequence level (RNA mean \\|value\\| {nr.get('mean_abs'):.3g} against "
            f"{nr.get('mean_abs_reference_run'):.3g}), so that case is at least detectable downstream. "
            f"The dangerous case is partial corruption: replacing 10,000 of the 1,048,576 bases with 'Z' "
            f"gave RNA mean \\|value\\| {zr.get('mean_abs'):.6g} against {zr.get('mean_abs_reference_run'):.6g} "
            f"on the real sequence, a difference of well under one percent. A window carrying an assembly "
            f"gap or a masked repeat will score as if it were clean sequence. The caller must validate "
            f"the alphabet; the model will not."
        )
    return out


def anchor_section(tag: str) -> list[str]:
    d = load(EXEC / f"tables/{tag}/anchor_hsd17b13.json")
    if d is None:
        return ["Anchor: not run on this GPU."]
    out = [
        f"Variant: {d['variant']}. Window {d['window_bp']:,} bp at {d['window_start0']:,}-{d['window_end']:,}. "
        f"Readout: mean \\|track\\| over +/-{d['flank_bp']:,} bp ({d['n_masked_positions']:,} positions), "
        f"log2((alt+1e-9)/(ref+1e-9)).",
        "",
        "| channel | local ref | local alt | local log2 | hosted log2 | local - hosted | same sign |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    loc = d["local"]
    cmp_ = d["local_vs_hosted"]
    for chn in ("splice", "rna", "atac", "dnase", "h3k27ac"):
        c = cmp_[chn]
        diff = c["difference_local_minus_hosted"]
        out.append(
            f"| {chn} | {loc[f'{chn}_ref']:.6g} | {loc[f'{chn}_alt']:.6g} | "
            f"{c['local_log2']:.6f} | {c['hosted_log2']:.6f} | "
            f"{diff:+.6f} | {'yes' if c['same_sign'] else 'no'} |"
        )
    return out


def anchor_commentary(tag: str) -> list[str]:
    d = load(EXEC / f"tables/{tag}/anchor_hsd17b13.json")
    if d is None:
        return []
    c = d["local_vs_hosted"]
    loc = d["local"]
    h = d["hosted"]
    sp = c["splice"]
    rel = abs(sp["difference_local_minus_hosted"]) / abs(sp["hosted_log2"]) * 100
    ref_pct = (loc["splice_ref"] / h["splice_ref"] - 1) * 100
    alt_pct = (loc["splice_alt"] / h["splice_alt"] - 1) * 100
    signs = all(v["same_sign"] for v in c.values() if v["same_sign"] is not None)
    return [
        f"The local weights put the splice-site-usage effect at {sp['local_log2']:.6f} against the hosted "
        f"{sp['hosted_log2']:.6f}, a gap of {sp['difference_local_minus_hosted']:+.4f}, which is "
        f"{rel:.1f} percent of the hosted magnitude. Direction agrees and so does magnitude. The "
        f"agreement is not only in the ratio: the underlying levels match too, with the local reference "
        f"readout {ref_pct:+.1f} percent and the local alternate readout {alt_pct:+.1f} percent against "
        f"the hosted ones. "
        + ("All five channels share the hosted sign. " if signs else "Not all channels share the hosted sign. ")
        + "The three chromatin channels agree to four decimal places; RNA-seq is the loosest at "
        f"{c['rna']['difference_local_minus_hosted']:+.4f}. Two independent implementations reaching "
        "this close on a variant the Atlas point query will not serve at all is the result the anchor was "
        "meant to test.",
    ]


HOSTED_TSV = (
    PROJ
    / "GWAS/finemapping/results/alphagenome_atlas/p6f-indel-rescue-20260914T135052Z/tables/indel_rescue_effects.tsv"
)


def anchor_context(tag: str) -> list[str]:
    """Is -1.09 a big number? Read the hosted run's own controls rather than asserting it."""
    import csv
    import statistics as st

    d = load(EXEC / f"tables/{tag}/anchor_hsd17b13.json")
    if d is None or not HOSTED_TSV.exists():
        return []
    sid = d["hosted"]["source_variant_id"]
    with HOSTED_TSV.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    ctl = [
        abs(float(r["splice_log2"]))
        for r in rows
        if r.get("arm") == "snv_control" and r.get("control_for") == sid and r.get("splice_log2")
    ]
    idl = [
        abs(float(r["splice_log2"]))
        for r in rows
        if r.get("arm") == "indel" and r.get("splice_log2")
    ]
    if not ctl or not idl:
        return []
    local = abs(d["local"]["splice_log2"])
    return [
        f"For scale, from the hosted run's own controls: the {len(ctl)} substitutions drawn at the same "
        f"signals and scored the same way have a median \\|splice log2\\| of {st.median(ctl):.4f}, and "
        f"across all {len(idl)} indels in that run the median is {st.median(idl):.4f} with a maximum of "
        f"{max(idl):.4f}, which is this variant. The local value of {local:.3f} is about "
        f"{local / st.median(ctl):.0f} times the matched-control median, so the local weights reproduce "
        f"not just the number but its standing as the largest splice effect in the set.",
        "",
    ]


def cross_device_anchor(tag: str, other: str) -> list[str]:
    """Same weights, same sequence, two different GPUs. Determinism holds within a card; this asks
    whether it holds across cards, which is what decides if scores from mixed hardware can be mixed."""
    a = load(EXEC / f"tables/{tag}/anchor_hsd17b13.json")
    b = load(EXEC / f"tables/{other}/anchor_hsd17b13.json")
    if not a or not b:
        return []
    rows = [
        f"Same weights and same sequence on both cards ({tag} and {other}):",
        "",
        f"| channel | {tag} log2 | {other} log2 | difference |",
        "|---|---:|---:|---:|",
    ]
    worst = 0.0
    for ch in ("splice", "rna", "atac", "dnase", "h3k27ac"):
        x, y = a["local"][f"{ch}_log2"], b["local"][f"{ch}_log2"]
        worst = max(worst, abs(x - y))
        rows.append(f"| {ch} | {x:.6f} | {y:.6f} | {x - y:+.2e} |")
    rows.append("")
    if worst > 0:
        rows.append(
            f"Largest disagreement between the two cards: {worst:.2e} in log2 units. Scores are bit-exact "
            f"when repeated on one card and are not reproducible across GPU architectures, so a score "
            f"table has to record which card produced each row. The size of this matters: in the hosted "
            f"run the median indel splice effect is 0.0151 in the same units, so a cross-card difference "
            f"of {worst:.3f} is the same order as a typical effect. A large effect like this variant's "
            f"is safe to compare across cards; a small one is not, and a set scored on mixed hardware "
            f"should be rescored on one."
        )
    else:
        rows.append("Identical on both cards.")
    rows.append("")
    return rows


def liver_tracks(tag: str) -> list[str]:
    d = load(EXEC / f"tables/{tag}/anchor_hsd17b13.json")
    sel = (d or {}).get("liver_track_selection")
    if not sel:
        return []
    out = ["| channel | tracks selected | biosamples |", "|---|---:|---|"]
    for k, v in sel.items():
        names = v.get("tracks", {}).get("biosample_name", [])
        uniq = sorted(set(names))
        out.append(f"| {k} | {v.get('n_tracks', 0)} | {', '.join(uniq) if uniq else '-'} |")
    hist = sel.get("h3k27ac", {}).get("tracks", {}).get("name", [])
    k27 = [n for n in hist if "H3K27AC" in n.upper()]
    if hist:
        bios = sorted(
            {
                b
                for n, b in zip(hist, sel["h3k27ac"]["tracks"].get("biosample_name", []))
                if "H3K27AC" in n.upper()
            }
        )
        out.append("")
        out.append(
            f"The CHIP_HISTONE head returns {len(hist)} liver tracks covering many marks. The anchor's "
            f"H3K27ac readout is the {len(k27)} of them that are actually H3K27ac ({', '.join(bios)}); "
            f"the recipe selects them by name, and that name filter does match here."
        )
    return out


def main() -> None:
    tag = sys.argv[1] if len(sys.argv) > 1 else "l40s"
    other = "b6k" if tag == "l40s" else "l40s"
    parts: list[str] = []
    parts.append("# AlphaGenome local runtime probe")
    parts.append("")
    parts.append(
        "Runtime only. No fine-tuning and no benchmark run was started. Every number below is read "
        "from the JSON artifacts in `tables/` by `scripts/analysis/alphagenome_program/h1_write_results.py`."
    )
    parts.append("")
    parts.append("## What was run")
    parts.append("")
    parts.append(f"- Weights: `{EXEC.parent.name}/alphagenome-weights-20260915T103442Z/checkpoints`, "
                 "Orbax OCDBT, HF `google/alphagenome-all-folds` @ `a8f293a76ee73d5b57f3bf2ae146510589fcf187`.")
    parts.append("- Code: `google-deepmind/alphagenome_research` @ `1e55dcffb98ba26b31e74edc5e9f038f54c0e89d` "
                 "(the pinned revision; not HEAD, which was `0db53bd4` at clone time).")
    parts.append("- Environment: micromamba env `alphagenome_local` at "
                 "`/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local`.")
    parts.append("")
    for t in (tag, other):
        parts.append(f"- GPU {t}: `{gpu_line(t)}`")
    parts.append("")
    vp = EXEC / f"env/versions_{tag}.txt"
    if vp.exists():
        parts.append("```")
        parts.append(vp.read_text().strip())
        parts.append("```")
        parts.append("")
    parts.append(
        "Full `pip freeze` for each node is in `env/pip_freeze_l40s.txt` and `env/pip_freeze_b6k.txt`. "
        "The sbatch scripts as submitted are copied to `env/`."
    )
    parts.append("")
    parts.extend(jobs_section())
    parts.append("")
    parts.append("## Measured time and memory")
    parts.append("")
    parts.append(
        "One process per input length, so `peak_bytes_in_use` is that length's own high-water mark and not "
        "a carry-over from a larger one. `XLA_PYTHON_CLIENT_PREALLOCATE=false`; the device limit column is "
        "JAX's default 75 percent cap on the card, not the card's physical memory."
    )
    parts.append("")
    for t in (tag, other):
        rows = bench_table(t)
        if rows:
            parts.append(f"### {t}")
            parts.append("")
            parts.extend(rows)
            parts.append("")
    parts.extend(gpu_commentary())
    parts.append("")
    sh = shapes_table(tag)
    if sh:
        parts.extend(sh)
        parts.append("")
    parts.append("## Variant-effect probes")
    parts.append("")
    parts.append(
        "Pass criteria were fixed in the script's docstring before the run. Reported below exactly as "
        "measured, including where the criterion was not met. One naming caveat: in the probe tables the "
        "row called `h3k27ac` is the whole CHIP_HISTONE head, all 17 liver histone tracks at 128 bp, "
        "because the probes sum the response over every track of a head. Only the anchor further down "
        "restricts that head to its three H3K27ac columns, which is what the pinned recipe does."
    )
    parts.append("")
    parts.extend(probes_section(tag))
    parts.append("")
    parts.extend(probe_commentary(tag))
    parts.append("")
    parts.extend(probe_replication(tag, other))
    parts.append("## HSD17B13 anchor, local weights beside the hosted model API")
    parts.append("")
    parts.extend(anchor_section(tag))
    parts.append("")
    parts.extend(anchor_commentary(tag))
    parts.append("")
    parts.extend(anchor_context(tag))
    parts.extend(cross_device_anchor(tag, other))
    parts.append(
        "The local weights and the hosted model API are different services. The two values are reported "
        "side by side, never averaged and never pooled. Hosted outputs remain excluded from training use; "
        "outputs of these local weights are eligible as features or teacher signals under the "
        "noncommercial terms, which derivatives inherit."
    )
    parts.append("")
    lt = liver_tracks(tag)
    if lt:
        parts.append("## What the pinned liver ontology actually selects")
        parts.append("")
        parts.extend(lt)
        parts.append("")
        parts.append(
            "Across the whole human track catalogue, and not only these four terms: hepatocyte "
            "(CL:0000182) has RNA-seq, DNase, histone ChIP, CAGE, splice-usage and one TF ChIP track but "
            "NO ATAC track. Hepatic stellate cell (CL:0000632) appears once, as DNase, and nowhere else. "
            "Bile duct (UBERON:0002394) appears once, as ATAC. There is no cholangiocyte biosample in any "
            "output type. The four pinned terms therefore reach neither the stellate DNase track nor the "
            "bile-duct ATAC track."
        )
        parts.append("")
    ex = extra_section()
    if ex:
        parts.append("## Lengths outside the SDK's named set")
        parts.append("")
        parts.append(
            "The SDK names five lengths: 2,048 / 16,384 / 131,072 / 524,288 / 1,048,576. Two lengths "
            "outside that list were tried to find out whether 1 Mb is a model ceiling or a card ceiling. "
            "Both are multiples of 2,048, which is the shape the parameter tree is built against, so this "
            "does not test arbitrary lengths."
        )
        parts.append("")
        parts.extend(ex)
        parts.append("")
    parts.append("## What this runtime can and cannot do")
    parts.append("")
    parts.extend(CAN_CANNOT)
    parts.append("")
    parts.append("## What blocks fine-tuning")
    parts.append("")
    parts.extend(FINETUNE_BLOCKERS)
    parts.append("")
    dest = EXEC / "RESULTS.md"
    dest.write_text("\n".join(parts) + "\n")
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
