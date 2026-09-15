#!/usr/bin/env python3
"""Step 68b: check step 68's indel orientations by routes that share none of its code.

Step 68 orients 112,863 deferred indels, 67,377 of them only through an allele-frequency tie-break. Before
anything is rescored on those orientations they are checked four ways (prespecification
68b_orientation_check_prespec.json, written before this script ran):

C1  Truth from the studies. Some GWAS files store alleles in a fixed order: measured on their SNVs, UKBB
    lists ALT first and MVP and BBJ list REF first, while FinnGen and most single-cohort files are mixed.
    That order is applied to the same file's indels, and the assumption that indels follow it is itself
    tested on the indels the reference alone can orient. Accuracy is also split by consortium, so a truth
    set that one harmonisation pipeline might have oriented against dbSNP can be seen on its own.
C2  A seeded sample re-derived with bcftools (query, then norm against the repository FASTA) and a matcher
    and FREQ parser written separately from steps 68 and 69.
C3  Cross-tabs: the frequency tie-break against the hg19 exact-record arm, the DISAGREE rows, flip rates by
    rule, and the like-for-like disagreement count in the run before the tie-break existed.
C4  What the final table changes for the 152 rescue targets.

Outputs (tables/): orientation_check.json, study_allele_conventions.tsv, bcftools_rederivation.tsv,
truth_set_disagreements.tsv
"""

from __future__ import annotations

import csv
import hashlib
import json
import pathlib
import random
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lib_atlas as la

HERE = pathlib.Path(__file__).resolve().parent
PRESPEC = HERE / "68b_orientation_check_prespec.json"
RESULTS = la.out_root().parent
TRACK0 = la.track0_root() / "tables"
VCF_HG38 = pathlib.Path("/gpfs/commons/home/jameslee/reference_genome/dbsnp/GCF_000001405.40.gz")
FASTA = pathlib.Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
                     "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")
PRE_FREQUENCY_RUN = RESULTS / "p0-dbsnp-orientation-20260914T174542Z" / "tables" / "dbsnp_indel_orientation.tsv"
PREVIOUS_RESCUE = RESULTS / "p6f-indel-rescue-20260914T135052Z" / "tables" / "indel_rescue_effects.tsv"
SEED = 20260914
PER_STRATUM = 300
FREQ_FLOOR = 0.01
CONVENTION_FLOOR = 0.01
MIN_SNVS = 200
WINDOW = 60


def study_convention(n_allele1_is_ref: int, n_snvs: int) -> str:
    """How a study orders its alleles, measured on SNVs where the reference base decides."""
    if n_snvs < MIN_SNVS:
        return "too_few"
    share = n_allele1_is_ref / n_snvs
    if share >= 1 - CONVENTION_FLOOR:
        return "allele1_is_ref"
    if share <= CONVENTION_FLOOR:
        return "allele1_is_alt"
    return "mixed"


def truth_orientation(a1: str, a2: str, conventions: set) -> tuple:
    """(ref, alt) implied by the fixed-order studies containing this indel, or None and why not."""
    fixed = {c for c in conventions if c in ("allele1_is_ref", "allele1_is_alt")}
    if len(fixed) > 1:
        return None, "conflicting_conventions"
    if not fixed:
        return None, "no_convention"
    return ((a1, a2) if fixed.pop() == "allele1_is_ref" else (a2, a1)), "truth"


def freq_per_alt(field: str, n_alts: int) -> list:
    """Max frequency per ALT across FREQ studies ('STUDY:ref,alt1,...|...'); '.' is unreported."""
    best = [None] * n_alts
    if not field or field == ".":
        return best
    for block in field.split("|"):
        name, _, values = block.partition(":")
        if not values:
            continue
        for k, text in enumerate(values.split(",")[1:]):
            if k >= n_alts or text in ("", "."):
                continue
            try:
                value = float(text)
            except ValueError:
                continue
            if best[k] is None or value > best[k]:
                best[k] = value
    return best


def classify_site(ours: dict, db: dict) -> tuple:
    """Orientation from normalised keys: `ours` key -> (ref, alt); `db` key -> list of per-ALT frequencies."""
    hits = {ours[k]: max((f for f in db[k] if f is not None), default=None) for k in ours if k in db}
    if len(hits) == 1:
        (ref, alt), = hits
        return "unique_record", ref, alt
    if not hits:
        return "absent", None, None
    common = [o for o, f in hits.items() if f is not None and f >= FREQ_FLOOR]
    if len(common) == 1:
        return "frequency", common[0][0], common[0][1]
    return ("both_common" if common else "neither_common"), None, None


def step68_state(row: dict) -> tuple:
    """Step 68's hg38 answer in the vocabulary of classify_site."""
    state, fb = row["dbsnp_hg38_state"], row.get("frequency_tie_break", "")
    if state == "resolved" and fb == "":
        return "unique_record", row["dbsnp_hg38_ref"], row["dbsnp_hg38_alt"]
    if state == "resolved" and fb == "resolved_by_frequency":
        return "frequency", row["dbsnp_hg38_ref"], row["dbsnp_hg38_alt"]
    if state == "both_in_dbsnp" and fb in ("both_common", "neither_common"):
        return fb, None, None
    if state == "absent" and fb == "":
        return "absent", None, None
    raise la.ContractError(f"unexpected step-68 state {state!r} with tie-break {fb!r}")


def latest(pattern: str) -> pathlib.Path:
    found = [p for p in sorted(RESULTS.glob(pattern)) if not (p.parents[1] / "SUPERSEDED.txt").exists()]
    if not found:
        raise la.ContractError(f"no current table matches {pattern}")
    return found[-1]


def uid_alleles(uid: str) -> tuple:
    parts = str(uid or "").split(":")
    return (parts[2], parts[3]) if len(parts) == 4 else None


# ---------------------------------------------------------------------------------------------------- C1
def study_conventions_and_indel_studies(indel_ids: set) -> tuple:
    gwas = {r["signal_uid"]: r["gwas_name"] for r in la.read_tsv(TRACK0 / "eligible_signals.tsv")}
    snv_ref = {}
    with la.open_text(TRACK0 / "variant_crosswalk.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["mapping_status"] == "mapped" and r["is_snv"] == "True" and r["palindromic"] == "False":
                snv_ref[r["source_variant_id"]] = r["hg38_ref"].upper()
    tally = defaultdict(lambda: [0, 0])
    seen = set()
    indel_studies = defaultdict(set)
    with la.open_text(TRACK0 / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            study = gwas.get(r["signal_uid"])
            if study is None:
                raise la.ContractError(f"signal {r['signal_uid']} has no study")
            sid = r["source_variant_id"]
            if sid in indel_ids:
                indel_studies[sid].add(study)
                continue
            ref = snv_ref.get(sid)
            if ref is None or (study, sid) in seen:
                continue
            seen.add((study, sid))
            m = la.SOURCE_ID_RE.match(sid)
            if m is None:
                raise la.ContractError(f"unparsable source id {sid!r}")
            tally[study][0] += m["ref"].upper() == ref
            tally[study][1] += 1
    conv = {s: study_convention(a, n) for s, (a, n) in tally.items()}
    rows = [{"study": s, "n_snvs": n, "n_allele1_is_ref": a, "share_allele1_is_ref": f"{a / n:.5f}",
             "convention": conv[s]} for s, (a, n) in sorted(tally.items())]
    return conv, indel_studies, rows


def study_family(study: str) -> str:
    """The consortium a study's files came from, so a truth set built by one pipeline can be separated."""
    for prefix in ("PanUKBB", "UKBB", "MVP", "BBJ", "FinnGen"):
        if study.startswith(prefix):
            return prefix
    return "other"


def accuracy(pairs) -> dict:
    n = len(pairs)
    k = sum(1 for got, want in pairs if got == want)
    return {"n": n, "correct": k, "accuracy": (k / n) if n else None,
            "determinate": n >= 200}


# ---------------------------------------------------------------------------------------------------- C2
def run(cmd, stdin=None) -> str:
    return subprocess.run(cmd, input=stdin, capture_output=True, check=True, text=True).stdout


def accession_map() -> dict:
    out = {}
    for line in run(["bcftools", "view", "-h", str(VCF_HG38)]).splitlines():
        if line.startswith("##contig=<ID="):
            acc = line.split("ID=", 1)[1].split(",")[0].rstrip(">")
            head, _, _ = acc.partition(".")
            if head.startswith("NC_0000") and head[7:].isdigit():
                n = int(head[7:])
                name = {23: "chrX", 24: "chrY"}.get(n, f"chr{n}" if 1 <= n <= 22 else None)
                if name:
                    out[name] = acc
    return out


def rederive(sample: list) -> list:
    acc = accession_map()
    lengths = {}
    order = {}
    for i, line in enumerate((FASTA.parent / (FASTA.name + ".fai")).read_text().splitlines()):
        name, length = line.split("\t")[:2]
        lengths[name], order[name] = int(length), i
    records = []
    for idx, row in enumerate(sample):
        chrom, pos = row["hg38_chrom"], int(row["hg38_position_1based"])
        region = f"{acc[chrom]}:{max(1, pos - WINDOW + 1)}-{pos + WINDOW}"
        text = run(["bcftools", "query", "-r", region, "-f", "%POS\t%REF\t%ALT\t%INFO/FREQ\n", str(VCF_HG38)])
        for line in text.splitlines():
            p, ref, alts, freq = line.split("\t")
            alt_list = alts.split(",")
            freqs = freq_per_alt(freq, len(alt_list))
            for alt, f in zip(alt_list, freqs):
                if len(alt) == len(ref) or not set(alt.upper()) <= set("ACGTN"):
                    continue
                info = f"TAG=s{idx}|d" + (f";AFI={f!r}" if f is not None else "")
                records.append((chrom, int(p), ref.upper(), alt.upper(), info))
        for ref, alt in ((row["source_allele1"], row["source_allele2"]),
                         (row["source_allele2"], row["source_allele1"])):
            records.append((chrom, pos, ref, alt, f"TAG=s{idx}|o|{ref}|{alt}"))
    records.sort(key=lambda r: (order[r[0]], r[1]))
    header = ["##fileformat=VCFv4.2",
              '##INFO=<ID=TAG,Number=1,Type=String,Description="site and kind">',
              '##INFO=<ID=AFI,Number=1,Type=Float,Description="max FREQ of this ALT">']
    header += [f"##contig=<ID={c},length={lengths[c]}>" for c in sorted({r[0] for r in records}, key=order.get)]
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    vcf = "\n".join(header + [f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t{i}" for c, p, r, a, i in records]) + "\n"
    normed = run(["bcftools", "norm", "-f", str(FASTA), "-c", "x", "-Ov", "-"], stdin=vcf)
    table = run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/TAG\t%INFO/AFI\n", "-"], stdin=normed)
    ours, db = defaultdict(dict), defaultdict(lambda: defaultdict(list))
    for line in table.splitlines():
        chrom, p, ref, alt, tag, afi = line.split("\t")
        key = (chrom, int(p), ref, alt)
        site, kind, *rest = tag.split("|")
        if kind == "o":
            ours[site][key] = (rest[0], rest[1])
        else:
            db[site][key].append(None if afi == "." else float(afi))
    out = []
    for idx, row in enumerate(sample):
        site = f"s{idx}"
        got = classify_site(ours[site], db[site])
        want = step68_state(row)
        out.append({"source_variant_id": row["source_variant_id"], "stratum": want[0],
                    "step68_state": want[0], "step68_ref": want[1] or "", "step68_alt": want[2] or "",
                    "bcftools_state": got[0], "bcftools_ref": got[1] or "", "bcftools_alt": got[2] or "",
                    "agree": str(got == want)})
    return out


# ---------------------------------------------------------------------------------------------------- main
def main() -> None:
    tables = la.out_root() / "tables"
    prespec_sha = hashlib.sha256(PRESPEC.read_bytes()).hexdigest()
    path68 = latest("p0-dbsnp-orientation-*/tables/dbsnp_indel_orientation.tsv")
    rows = la.read_tsv(path68)
    by_id = {r["source_variant_id"]: r for r in rows}
    la.log(f"step 68b: {len(rows)} rows from {path68}; prespec sha256 {prespec_sha}")
    out = {"prespec_sha256": prespec_sha, "step68_table": str(path68),
           "bcftools": run(["bcftools", "--version"]).splitlines()[0]}

    def rule(r):
        if not r["resolved_variant_uid"]:
            return "unresolved"
        return "frequency" if r["frequency_tie_break"] == "resolved_by_frequency" else "unique_record"

    # C3 cross-tabs
    xt = Counter((rule(r), r["dbsnp_hg19_state"], r["assemblies_agree"]) for r in rows)
    freq_with_hg19 = [r for r in rows if r["frequency_tie_break"] == "resolved_by_frequency"
                      and r["dbsnp_hg19_state"] == "resolved"]
    flips = defaultdict(Counter)
    for r in rows:
        if r["resolved_variant_uid"]:
            flips[(r["current_orientation_ambiguous"], rule(r))][r["verdict_vs_current"]] += 1
    pre = la.read_tsv(PRE_FREQUENCY_RUN)
    pre_conflicts = [r for r in pre if r["dbsnp_hg19_state"] == "resolved" and r["dbsnp_hg38_state"] == "resolved"
                     and (r["dbsnp_hg19_ref"], r["dbsnp_hg19_alt"]) != (r["dbsnp_hg38_ref"], r["dbsnp_hg38_alt"])]
    out["C3"] = {
        "rule_x_hg19_state_x_assemblies": {" / ".join(k): v for k, v in sorted(xt.items())},
        "frequency_tiebreak_where_hg19_exact_also_answers": {
            "n": len(freq_with_hg19),
            "agree": sum(r["assemblies_agree"] == "agree" for r in freq_with_hg19),
            "disagree": sum(r["assemblies_agree"] == "DISAGREE" for r in freq_with_hg19)},
        "disagree_rows": [{k: r[k] for k in ("source_variant_id", "dbsnp_hg19_ref", "dbsnp_hg19_alt",
                                             "dbsnp_hg38_ref", "dbsnp_hg38_alt", "frequency_tie_break",
                                             "resolved_alt_frequency", "dbsnp_hg19_rsid", "dbsnp_hg38_rsid")}
                          for r in rows if r["assemblies_agree"] == "DISAGREE"],
        "verdict_vs_step67_by_reference_ambiguity_and_rule": {f"ambiguous={a} / {b}": dict(c)
                                                               for (a, b), c in sorted(flips.items())},
        "pre_frequency_run_two_resolved_arms_with_different_alleles": len(pre_conflicts),
        "pre_frequency_run_disagree_label_count": sum(r["assemblies_agree"] == "DISAGREE" for r in pre)}

    # C1 truth from study conventions
    conv, indel_studies, conv_rows = study_conventions_and_indel_studies(set(by_id))
    la.write_tsv_once(tables / "study_allele_conventions.tsv", conv_rows,
                      ["study", "n_snvs", "n_allele1_is_ref", "share_allele1_is_ref", "convention"])
    truth_state = Counter()
    ref_check, by_rule, by_family, disagreements = [], defaultdict(list), defaultdict(list), []
    for sid, r in by_id.items():
        truth, why = truth_orientation(r["source_allele1"], r["source_allele2"],
                                       {conv.get(s, "too_few") for s in indel_studies.get(sid, set())})
        truth_state[why] += 1
        if truth is None:
            continue
        current = uid_alleles(r["current_variant_uid"])
        if r["current_orientation_ambiguous"] == "False" and current:
            ref_check.append((current, truth))
            continue
        if r["current_orientation_ambiguous"] != "True":
            continue
        got = uid_alleles(r["resolved_variant_uid"])
        by_rule[rule(r)].append((got, truth))
        fams = {study_family(s) for s in indel_studies.get(sid, set())
                if conv.get(s) in ("allele1_is_ref", "allele1_is_alt")}
        by_family[(fams.pop() if len(fams) == 1 else "several", rule(r))].append((got, truth))
        if got and got != truth and len(disagreements) < 500:
            disagreements.append({"source_variant_id": sid, "rule": rule(r), "truth_ref": truth[0],
                                  "truth_alt": truth[1], "step68_ref": got[0], "step68_alt": got[1],
                                  "resolved_alt_frequency": r["resolved_alt_frequency"],
                                  "studies": ";".join(sorted(indel_studies.get(sid, set())))})
        by_rule["source_order_first"].append(((r["source_allele1"], r["source_allele2"]), truth))
    la.write_tsv_once(tables / "truth_set_disagreements.tsv", disagreements,
                      ["source_variant_id", "rule", "truth_ref", "truth_alt", "step68_ref", "step68_alt",
                       "resolved_alt_frequency", "studies"])
    out["C1"] = {"study_conventions": Counter(conv.values()),
                 "truth_availability": dict(truth_state),
                 "assumption_test_reference_unambiguous_indels": accuracy(ref_check),
                 "reference_ambiguous_accuracy_by_rule": {
                     k: accuracy([(g, t) for g, t in v if g is not None]) | {"n_unresolved": sum(g is None for g, _ in v)}
                     for k, v in sorted(by_rule.items())},
                 "reference_ambiguous_accuracy_by_study_family_and_rule": {
                     f"{fam} / {k}": accuracy([(g, t) for g, t in v if g is not None])
                     for (fam, k), v in sorted(by_family.items()) if k != "unresolved"},
                 "disagree_rows_against_truth": {
                     r["source_variant_id"]: truth_orientation(
                         r["source_allele1"], r["source_allele2"],
                         {conv.get(s, "too_few") for s in indel_studies.get(r["source_variant_id"], set())})
                     for r in rows if r["assemblies_agree"] == "DISAGREE"}}

    # C2 bcftools re-derivation on a seeded sample
    rng = random.Random(SEED)
    strata = defaultdict(list)
    for r in rows:
        if r["hg38_position_1based"]:
            strata[step68_state(r)[0]].append(r)
    sample = []
    for name in sorted(strata):
        pool = sorted(strata[name], key=lambda r: r["source_variant_id"])
        sample += pool if len(pool) <= PER_STRATUM else rng.sample(pool, PER_STRATUM)
    ids = {r["source_variant_id"] for r in sample}
    sample += [r for r in rows if r["assemblies_agree"] == "DISAGREE" and r["source_variant_id"] not in ids]
    red = rederive(sample)
    la.write_tsv_once(tables / "bcftools_rederivation.tsv", red, list(red[0]))
    agree = defaultdict(Counter)
    for x in red:
        agree[x["stratum"]][x["agree"]] += 1
    out["C2"] = {k: {"n": sum(c.values()), "agree": c["True"], "share": c["True"] / sum(c.values())}
                 for k, c in sorted(agree.items())}

    # C4 rescue targets
    prev = [r for r in la.read_tsv(PREVIOUS_RESCUE) if r.get("arm") == "indel" and r.get("state") == "scored"]
    c4 = Counter()
    mass = Counter()
    for p in prev:
        r = by_id.get(p["source_variant_id"])
        new = uid_alleles(r["resolved_variant_uid"]) if r else None
        k = rule(r) if r else "not_in_table"
        c4[k] += 1
        if not new and p.get("orientation_source") == "dbsnp":
            c4["previously_dbsnp_now_unresolved"] += 1
        if new and new != (p["ref"], p["alt"]):
            c4[f"{k}_changes_scored_alleles"] += 1
            mass[k] += float(p["posterior_mass"])
    sites = {(p["chrom"], p["pos_hg38"], frozenset((p["ref"], p["alt"]))) for p in prev}
    out["C4"] = {"n_previous_targets": len(prev), "n_distinct_variants": len(sites), "counts": dict(c4),
                 "posterior_mass_changing": {k: round(v, 4) for k, v in mass.items()},
                 "total_posterior_mass": round(sum(float(p["posterior_mass"]) for p in prev), 4)}

    json.dump(out, (tables / "orientation_check.json").open("w"), indent=1, default=str)
    la.log("step 68b: " + json.dumps({"C1": out["C1"]["reference_ambiguous_accuracy_by_rule"],
                                      "C2": out["C2"]}, default=str))


if __name__ == "__main__":
    main()
