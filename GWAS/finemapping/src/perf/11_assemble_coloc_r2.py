#!/usr/bin/env python3
"""
Assemble the corrected COLOC release r2 (review item 2, Phase 2 step 5) and
summarize it against the adopted release.

WHAT CHANGES: only the 1,307 SuSiE-positive (study, chr, gene) rows planned in
results/susie_coloc_eligibility_<ts>/plan/rerun_plan.tsv. For those rows the
SuSiE fields (PP.H3.susie, PP.H4.susie, method, n_cs_pairs) are taken from the
eligibility rerun; every other field, including all ABF fields, stays as the
adopted master wrote it (byte-for-byte for unplanned rows). A new column
n_cs_pairs_eligible holds the number of signal pairs passing coloc's shared-
posterior check, and eligibility_source says why it is empty elsewhere:
  rerun_shared_posterior_check  planned row, value from the rerun
  not_rerun_negative            adopted SuSiE PP.H4 <= 0.5 (cannot become positive)
  not_rerun_no_susie_result     adopted row has no SuSiE result (abf_fallback/abf_only)

CHECKS (assemble fails and writes no ASSEMBLY_COMPLETE.json otherwise):
  - every planned row appears exactly once in the rerun, with method susie or
    susie_untestable_insufficient_shared_posterior, and the same n_cs_pairs;
  - rerun ABF PP.H0-H4 equal the adopted values within 1e-10 (reproduction);
  - the best eligible pair's PP.H4 does not exceed the adopted PP.H4.susie by
    more than 1e-4 (the documented non-EUR replay tolerance). coloc returns NA
    PP.H4 for ineligible pairs, so the adopted value of a pair that is now
    ineligible is reproduced only through the offline check below;
  - the 816 tier-1/2 gene-studies agree with the offline checker
    (10_offline_eqtl_overlap_check.R) on status, n_pairs, n_eligible and
    corrected PP.H4 (within 1e-4), and every signal pair agrees on eligibility
    (eligible pairs also on hit1/hit2 and PP.H4).

BEFORE ANY OF THAT, the five batch logs are scanned. 06 stops only AFTER writing
its outputs when its eligibility call disagrees with coloc's, and 09b skips a
task whose outputs exist on resubmission, so a mismatch can hide behind files
that look complete. Any ELIGIBILITY MISMATCH / stop message / R error / unscanned
resubmission log refuses the assembly.

ALSO WRITTEN (qc/, reviewer requests):
  eligibility_rule_comparison.tsv  every planned row under coloc 5.2.3's rule
      (last lbf column = null) and under the uniform-prior share; the call is
      re-derived under each rule, using replay PP.H4 for pairs the rerun dropped
  eligibility_rule_pair_differences.tsv  every signal pair the two rules disagree on
  call_changes.tsv  each planned row whose call changed: adopted PP.H4, the rerun's
      max over pairs with a finite PP.H4, and the replay's max over ALL pairs
      (tier 1/2 only), which must reproduce the adopted value within 1e-4
A study whose every pair fails the check (method UNTESTABLE, PP.H4.susie NA) is
its own class in every summary, never a negative.

ld_contamination_clusters.csv is copied from the adopted release: 07b builds it
from PP.H4.abf and the ABF top_snp only, and both are unchanged here.

Usage (compute node; see 11b_assemble_coloc_r2.sbatch for the full chain):
  python3 11_assemble_coloc_r2.py assemble  --release-root results/susie_coloc_r2_<UTC ts>
  python3 11_assemble_coloc_r2.py summarize --release-root results/susie_coloc_r2_<UTC ts>
"""
import argparse, collections, csv, glob, gzip, hashlib, io, json, math, os, re, shutil, sys
from datetime import datetime, timezone

FM = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
ADOPTED = os.path.join(FM, "results/susie_coloc")
RERUN = os.path.join(FM, "results/susie_coloc_eligibility_20260923T211525Z")
OFFLINE = os.path.join(FM, "results/coloc_eligibility_offline_20260923T211525Z")
REPLAY = os.path.join(os.path.dirname(os.path.dirname(FM)), "Analysis/Multimodal_Program_Projection",
                      "candidates/atac-context-v3-candidate-2026-08-11-r1/genetics/replay_execution")
LOG_DIR = os.path.join(FM, "logs/coloc_eligibility")
BATCH_JOBS = "21862335,21862336,21862337,21862338,21862339"
LOG_REFUSE = ("ELIGIBILITY MISMATCH", "disagrees with coloc.bf_bf", "Execution halted",
              "SuSiE-COLOC failed", "CANCELLED", "DUE TO TIME LIMIT", "oom-kill", "Killed")
LOG_ERROR = re.compile(r"^Error\b")
OVERLAP_MIN = 0.5   # coloc 5.2.3 overlap.min (src/coloc_signal_eligibility.R)

SUSIE_FIELDS = ["PP.H3.susie", "PP.H4.susie", "method", "n_cs_pairs"]
ABF_FIELDS = ["PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"]
NEW_FIELDS = ["n_cs_pairs_eligible", "eligibility_source"]
UNTESTABLE = "susie_untestable_insufficient_shared_posterior"
ABF_TOL = 1e-10
PP4_REPRO_TOL = 1e-4   # 22_validate_replay_batch.py rel_tol; STATUS.md replay-tolerance row
NAMED_GENES = ["GNMT", "MAT1A", "CYP2C19", "IGFBP7", "BICC1", "RORA", "FABP1"]
WORKED_LOCI = {"RORA": "MVP_NAFLD_EUR", "FABP1": "UKBB_ALT"}   # Fig 2G / 2H


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def absdiff(x, y):
    """|x - y| for two text fields; 0 when both are NA, inf when only one is."""
    a, b = num(x), num(y)
    if a is None and b is None:
        return 0.0
    return math.inf if a is None or b is None else abs(a - b)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def release_path(value):
    root = value if os.path.isabs(value) else os.path.join(FM, value)
    root = os.path.realpath(root)
    if not os.path.basename(root).startswith("susie_coloc_r2_"):
        sys.exit(f"release root must be named susie_coloc_r2_<UTC ts>: {root}")
    if root.startswith(os.path.realpath(ADOPTED) + os.sep):
        sys.exit("release root must not sit inside the adopted results/susie_coloc/")
    return root


def rerun_status(method, pp4):
    if method == UNTESTABLE:
        return "untestable_insufficient_shared_posterior"
    return "retained_pp4_gt_0.5" if pp4 is not None and pp4 > 0.5 else "eligible_pp4_le_0.5"


# ------------------------------------------------------------------ assemble
def load_rerun(rerun_root, plan):
    """Rerun rows keyed (gwas, chr, ensembl) and signal pairs keyed (gwas, ensembl)."""
    tasks = sorted({(r["gwas_name"], r["chr"]) for r in plan}, key=lambda t: (t[0], int(t[1])))
    missing, rows, pairs, pair_header = [], {}, collections.defaultdict(list), None
    for gwas, chrom in tasks:
        csv_path = os.path.join(rerun_root, "coloc", gwas, f"susie_coloc_chr{chrom}.csv")
        pair_path = os.path.join(rerun_root, "coloc", gwas, f"susie_coloc_pairs_chr{chrom}.tsv")
        if not (os.path.isfile(csv_path) and os.path.isfile(pair_path)):
            missing.append(f"{gwas} chr{chrom}")
            continue
        with open(csv_path, newline="") as fh:
            for r in csv.DictReader(fh):
                key = (r["gwas_name"], r["chr"], r["ensembl"])
                if key in rows:
                    sys.exit(f"rerun row duplicated: {key}")
                rows[key] = r
        with open(pair_path, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if pair_header is None:
                pair_header = rd.fieldnames
            elif rd.fieldnames != pair_header:
                sys.exit(f"signal-pair header differs in {pair_path}")
            for r in rd:
                pairs[(r["gwas"], r["ensembl"])].append(r)
    if missing:
        sys.exit(f"{len(missing)} of {len(tasks)} rerun tasks have no output yet, e.g. {missing[:5]}")
    return rows, pairs, pair_header


def validate_plan_rows(plan, rows, pairs, adopted):
    """Per planned row: completeness, method, n_cs_pairs and PP.H4 reproduction."""
    fail, log = [], []
    planned = {(r["gwas_name"], r["chr"], r["ensembl"]) for r in plan}
    extra = sorted(set(rows) - planned)
    if extra:
        fail.append(f"{len(extra)} rerun rows were not planned, e.g. {extra[:3]}")
    for p in plan:
        key = (p["gwas_name"], p["chr"], p["ensembl"])
        r, a = rows.get(key), adopted.get(key)
        if r is None:
            fail.append(f"planned row absent from rerun: {key}")
            continue
        if a is None:
            fail.append(f"planned row absent from adopted master: {key}")
            continue
        pr = pairs.get((key[0], key[2]), [])
        elig_pp4 = [num(x["PP.H4"]) for x in pr if x["eligible"] == "TRUE"]
        if any(v is None for v in elig_pp4) or any(num(x["PP.H4"]) is not None for x in pr
                                                    if x["eligible"] != "TRUE"):
            problems_pp4 = ["PP.H4 is NA for an eligible pair or finite for an ineligible one"]
        else:
            problems_pp4 = []
        max_elig = max(elig_pp4) if elig_pp4 and None not in elig_pp4 else None
        new_pp4, old_pp4 = num(r["PP.H4.susie"]), num(a["PP.H4.susie"])
        delta = new_pp4 - old_pp4 if new_pp4 is not None and old_pp4 is not None else None
        problems = problems_pp4
        if r["method"] not in ("susie", UNTESTABLE):
            problems.append(f"method {r['method']}")
        if num(r["n_cs_pairs"]) != num(a["n_cs_pairs"]) or len(pr) != num(r["n_cs_pairs"]):
            problems.append(f"n_cs_pairs adopted {a['n_cs_pairs']} rerun {r['n_cs_pairs']} pairs {len(pr)}")
        if str(len(elig_pp4)) != r["n_cs_pairs_eligible"]:
            problems.append(f"n_cs_pairs_eligible {r['n_cs_pairs_eligible']} vs {len(elig_pp4)} eligible pairs")
        if delta is not None and delta > PP4_REPRO_TOL:
            problems.append(f"PP.H4.susie rose from {old_pp4} to {new_pp4}")
        if r["method"] == "susie" and (max_elig is None or new_pp4 is None or abs(max_elig - new_pp4) > 1e-12):
            problems.append(f"PP.H4.susie {new_pp4} is not the best eligible pair {max_elig}")
        if r["method"] == UNTESTABLE and (elig_pp4 or new_pp4 is not None):
            problems.append("untestable row carries an eligible pair or a PP.H4")
        if problems:
            fail.append(f"{key}: " + "; ".join(problems))
        status = rerun_status(r["method"], new_pp4)
        log.append({"gwas_name": key[0], "chr": key[1], "gene": p["gene"], "ensembl": key[2],
                    "tier": p["tier"], "tier_label": p["tier_label"],
                    "adopted_pp_h4_susie": a["PP.H4.susie"], "r2_pp_h4_susie": r["PP.H4.susie"],
                    "adopted_method": a["method"], "r2_method": r["method"],
                    "n_cs_pairs": r["n_cs_pairs"], "n_cs_pairs_eligible": r["n_cs_pairs_eligible"],
                    "r2_minus_adopted_pp4": "" if delta is None else delta,
                    "best_pair_unchanged": str(delta is not None and abs(delta) <= PP4_REPRO_TOL).upper(),
                    "r2_status": status})
    return fail, log


def offline_crosscheck(offline_root, rows, pairs):
    """816 tier-1/2 gene-studies and their signal pairs against the offline checker."""
    by_study_gene = {(k[0], k[2]): v for k, v in rows.items()}
    bad, label_diffs, max_pp4 = [], [], 0.0
    gene_rows = read_tsv(os.path.join(offline_root, "eqtl_overlap_check_gene_study.tsv"))
    for o in gene_rows:
        r = by_study_gene.get((o["gwas_name"], o["ensembl"]))
        if r is None:
            bad.append({"level": "gene_study", "gwas_name": o["gwas_name"], "ensembl": o["ensembl"],
                        "field": "presence", "offline": "present", "rerun": "absent"})
            continue
        new_pp4, off_pp4 = num(r["PP.H4.susie"]), num(o["corrected_pp4"])
        checks = [("status", o["status"], rerun_status(r["method"], new_pp4)),
                  ("n_pairs", o["n_pairs"], r["n_cs_pairs"]),
                  ("n_eligible", o["n_eligible"], r["n_cs_pairs_eligible"])]
        for field, a, b in checks:
            if a != b:
                bad.append({"level": "gene_study", "gwas_name": o["gwas_name"], "ensembl": o["ensembl"],
                            "field": field, "offline": a, "rerun": b})
        if (new_pp4 is None) != (off_pp4 is None) or \
                (new_pp4 is not None and abs(new_pp4 - off_pp4) > PP4_REPRO_TOL):
            bad.append({"level": "gene_study", "gwas_name": o["gwas_name"], "ensembl": o["ensembl"],
                        "field": "corrected_pp4", "offline": o["corrected_pp4"], "rerun": r["PP.H4.susie"]})
        elif new_pp4 is not None:
            max_pp4 = max(max_pp4, abs(new_pp4 - off_pp4))
    rerun_pairs = {}
    for (gwas, ens), prs in pairs.items():
        for p in prs:
            rerun_pairs[(gwas, ens, p["idx1"], p["idx2"])] = p
    off_pairs = read_tsv(os.path.join(offline_root, "eqtl_overlap_check.tsv"))
    for o in off_pairs:
        key = (o["gwas_name"], o["ensembl"], o["idx1"], o["idx2"])
        p = rerun_pairs.get(key)
        if p is None:
            bad.append({"level": "signal_pair", "gwas_name": key[0], "ensembl": key[1],
                        "field": f"presence idx {key[2]}/{key[3]}", "offline": "present", "rerun": "absent"})
            continue
        # hit1/hit2 of an ineligible pair are not comparable: the replay cut both
        # fits to the shared SNPs before naming the top SNP, the rerun does not.
        checks = [("eligible", o["eligible"], p["eligible"])]
        pp4_agrees = True
        if o["eligible"] == "TRUE":
            if absdiff(o["replay_PP.H4"], p["PP.H4"]) > PP4_REPRO_TOL:
                checks.append(("PP.H4", o["replay_PP.H4"], p["PP.H4"]))
                pp4_agrees = False
            # A differing top-SNP label with the same eligibility and PP.H4 is a
            # near-tie in the refitted SuSiE signal (non-EUR Gram LD panels are not
            # positive semidefinite, so these fits are not bit-reproducible). It is
            # recorded, not treated as a disagreement; with a PP.H4 difference it is.
            for field in ("hit1", "hit2"):
                if o[field] != p[field]:
                    entry = {"level": "signal_pair", "gwas_name": key[0], "ensembl": key[1],
                             "field": f"{field} idx {key[2]}/{key[3]}", "offline": o[field],
                             "rerun": p[field]}
                    if pp4_agrees and o["eligible"] == p["eligible"]:
                        label_diffs.append(entry)
                    else:
                        bad.append(entry)
        for field, a, b in checks:
            if a != b:
                bad.append({"level": "signal_pair", "gwas_name": key[0], "ensembl": key[1],
                            "field": f"{field} idx {key[2]}/{key[3]}", "offline": a, "rerun": b})
    return bad, label_diffs, len(gene_rows), len(off_pairs), max_pp4


def check_batch_logs(log_dir, jobs):
    """Refusal reasons from the eligibility batch logs (empty list = clean)."""
    problems, n_overlap_warnings = [], 0
    for job in jobs:
        for ext in ("out", "err"):
            path = os.path.join(log_dir, f"coloc_elig_{job}.{ext}")
            if not os.path.isfile(path):
                problems.append(f"missing log {path}")
                continue
            with open(path, errors="replace") as fh:
                text = fh.read()
            for n, line in enumerate(text.splitlines(), 1):
                if LOG_ERROR.match(line) or any(s in line for s in LOG_REFUSE):
                    problems.append(f"{os.path.basename(path)}:{n}: {line.strip()[:160]}")
            n_overlap_warnings += text.count("snp overlap too small between datasets")
            if ext == "out" and not re.search(r"\] batch \d+ done$", text, re.M):
                problems.append(f"{os.path.basename(path)} has no 'batch <k> done' line")
    listed = {f"coloc_elig_{job}.out" for job in jobs}
    others = sorted(os.path.basename(p) for p in glob.glob(os.path.join(log_dir, "coloc_elig_*.out"))
                    if os.path.basename(p) not in listed)
    if others:
        problems.append(f"batch logs not scanned (add their job ids to --batch-jobs): {others}")
    return problems, n_overlap_warnings


def load_replay(replay_root):
    """Replay PP.H4 of EVERY signal pair (tier 1/2 only), keyed (study, gene) -> {(idx1, idx2): pp4}."""
    replay = collections.defaultdict(dict)
    for path in glob.glob(os.path.join(replay_root, "batch_*/exports/*/chr*/*/signal_pairs.tsv")):
        for r in read_tsv(path):
            replay[(r["gwas_name"], r["ensembl"])][(r["idx1"], r["idx2"])] = num(r["PP.H4.abf"])
    return replay


def uniform_eligible(p):
    a, b = num(p["prop_gwas_uniform"]), num(p["prop_eqtl_uniform"])
    return a is not None and b is not None and a >= OVERLAP_MIN and b >= OVERLAP_MIN


def call_from(eligible_pp4):
    """Status from the PP.H4 of the eligible pairs (None = PP.H4 not computed)."""
    if not eligible_pp4:
        return "untestable_insufficient_shared_posterior", None
    known = [v for v in eligible_pp4 if v is not None]
    best = max(known) if known else None
    if best is not None and best > 0.5:
        return "retained_pp4_gt_0.5", best
    if len(known) < len(eligible_pp4):
        return "undetermined_pp4_not_computed", best
    return "eligible_pp4_le_0.5", best


def eligibility_rules(plan, pairs, log, replay):
    """Per planned row: the call under coloc's rule and under the uniform-prior share."""
    status = {(r["gwas_name"], r["ensembl"]): r for r in log}
    rows, pair_diffs = [], []
    for p in plan:
        key = (p["gwas_name"], p["ensembl"])
        if key not in status:          # already a validation failure
            continue
        prs, rp = pairs.get(key, []), replay.get(key, {})
        pp4 = {}
        for x in prs:
            v = num(x["PP.H4"])
            pp4[(x["idx1"], x["idx2"])] = (v, "rerun") if v is not None else \
                (rp.get((x["idx1"], x["idx2"])), "replay" if (x["idx1"], x["idx2"]) in rp else "not_computed")
        coloc_rule = [pp4[(x["idx1"], x["idx2"])][0] for x in prs if x["eligible"] == "TRUE"]
        unif_rule = [pp4[(x["idx1"], x["idx2"])][0] for x in prs if uniform_eligible(x)]
        s_coloc, pp_coloc = call_from(coloc_rule)
        s_unif, pp_unif = call_from(unif_rule)
        if s_coloc != status[key]["r2_status"]:
            raise SystemExit(f"rule re-derivation disagrees with the rerun for {key}")
        differing = [x for x in prs if (x["eligible"] == "TRUE") != uniform_eligible(x)]
        for x in differing:
            v, src = pp4[(x["idx1"], x["idx2"])]
            pair_diffs.append({"gwas_name": key[0], "gene": p["gene"], "ensembl": key[1], "tier": p["tier"],
                               "idx1": x["idx1"], "idx2": x["idx2"], "hit1": x["hit1"], "hit2": x["hit2"],
                               "prop_gwas": x["prop_gwas"], "prop_eqtl": x["prop_eqtl"],
                               "prop_gwas_uniform": x["prop_gwas_uniform"],
                               "prop_eqtl_uniform": x["prop_eqtl_uniform"],
                               "eligible_coloc_rule": x["eligible"],
                               "eligible_uniform_rule": str(uniform_eligible(x)).upper(),
                               "pp4": "" if v is None else v, "pp4_source": src})
        rows.append({"gwas_name": key[0], "chr": p["chr"], "gene": p["gene"], "ensembl": key[1],
                     "tier": p["tier"], "n_pairs": len(prs), "n_eligible_coloc_rule": len(coloc_rule),
                     "n_eligible_uniform_rule": len(unif_rule), "n_pairs_rules_differ": len(differing),
                     "status_coloc_rule": s_coloc, "pp4_coloc_rule": "" if pp_coloc is None else pp_coloc,
                     "status_uniform_rule": s_unif, "pp4_uniform_rule": "" if pp_unif is None else pp_unif,
                     "call_differs": str(s_coloc != s_unif).upper()})
    return rows, pair_diffs


def call_changes(log, pairs, replay):
    """Planned rows whose call changed, with the adopted PP.H4 reproduced before the drop."""
    out, fail = [], []
    for r in log:
        if r["r2_status"] == "retained_pp4_gt_0.5":
            continue
        key = (r["gwas_name"], r["ensembl"])
        finite = [num(x["PP.H4"]) for x in pairs.get(key, []) if num(x["PP.H4"]) is not None]
        rp = [v for v in replay.get(key, {}).values() if v is not None]
        adopted = num(r["adopted_pp_h4_susie"])
        diff = abs(max(rp) - adopted) if rp else None
        if diff is not None and diff > PP4_REPRO_TOL:
            fail.append(f"replay max PP.H4 {max(rp)} does not reproduce adopted {adopted} for {key}")
        out.append({"gwas_name": key[0], "chr": r["chr"], "gene": r["gene"], "ensembl": key[1],
                    "tier": r["tier"], "adopted_pp_h4_susie": r["adopted_pp_h4_susie"],
                    "r2_status": r["r2_status"], "r2_pp_h4_susie": r["r2_pp_h4_susie"],
                    "n_pairs": len(pairs.get(key, [])), "n_pairs_pp4_na": len(pairs.get(key, [])) - len(finite),
                    "rerun_max_pp4_finite_pairs": max(finite) if finite else "",
                    "replay_max_pp4_all_pairs": max(rp) if rp else "",
                    "replay_minus_adopted": "" if diff is None else max(rp) - adopted,
                    "adopted_reproduced_before_drop": "not_in_replay (tier 3/4)" if diff is None
                    else str(diff <= PP4_REPRO_TOL).upper()})
    return out, fail


def assemble(args):
    root = release_path(args.release_root)
    if os.path.lexists(root):
        sys.exit(f"refusing to reuse {root}")
    log_problems, n_overlap_warnings = check_batch_logs(args.log_dir, args.batch_jobs.split(","))
    if log_problems:
        print("*** refusing to assemble: the eligibility batch logs report problems")
        for line in log_problems[:40]:
            print(f"    - {line}")
        sys.exit(2)
    plan = read_tsv(os.path.join(args.rerun_root, "plan/rerun_plan.tsv"))
    keys = [(r["gwas_name"], r["chr"], r["ensembl"]) for r in plan]
    if len(keys) != len(set(keys)):
        sys.exit("rerun plan has duplicated (study, chr, gene) rows")
    planned = set(keys)
    rows, pairs, pair_header = load_rerun(args.rerun_root, plan)
    print(f"rerun: {len(rows):,} rows, {sum(len(v) for v in pairs.values()):,} signal pairs, "
          f"{len(planned):,} planned")

    # One pass over the adopted master: keep the planned rows for validation.
    adopted_planned = {}
    with open(args.master, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        col = {c: i for i, c in enumerate(header)}
        for f in SUSIE_FIELDS + ABF_FIELDS + ["gwas_name", "chr", "ensembl", "n_snps", "top_snp", "top_snp_PP"]:
            if f not in col:
                sys.exit(f"adopted master lacks column {f}")
        for f in NEW_FIELDS:
            if f in col:
                sys.exit(f"adopted master already has column {f}")
        for row in rd:
            key = (row[col["gwas_name"]], row[col["chr"]], row[col["ensembl"]])
            if key in planned:
                if key in adopted_planned:
                    sys.exit(f"adopted master has duplicated planned row {key}")
                adopted_planned[key] = dict(zip(header, row))

    fail, log = validate_plan_rows(plan, rows, pairs, adopted_planned)
    bad, label_diffs, n_gene_off, n_pair_off, max_off_pp4 = offline_crosscheck(args.offline_root, rows, pairs)
    replay = load_replay(args.replay_root)
    rule_rows, rule_pairs = eligibility_rules(plan, pairs, log, replay)
    changes, change_fail = call_changes(log, pairs, replay)
    fail += change_fail

    os.makedirs(os.path.join(root, "qc"))
    log_cols = list(log[0]) if log else ["gwas_name"]
    write_tsv(os.path.join(root, "qc/replacement_log.tsv"), log_cols, log)
    write_tsv(os.path.join(root, "qc/eligibility_rule_comparison.tsv"), list(rule_rows[0]), rule_rows)
    write_tsv(os.path.join(root, "qc/eligibility_rule_pair_differences.tsv"),
              ["gwas_name", "gene", "ensembl", "tier", "idx1", "idx2", "hit1", "hit2", "prop_gwas",
               "prop_eqtl", "prop_gwas_uniform", "prop_eqtl_uniform", "eligible_coloc_rule",
               "eligible_uniform_rule", "pp4", "pp4_source"], rule_pairs)
    write_tsv(os.path.join(root, "qc/call_changes.tsv"),
              ["gwas_name", "chr", "gene", "ensembl", "tier", "adopted_pp_h4_susie", "r2_status",
               "r2_pp_h4_susie", "n_pairs", "n_pairs_pp4_na", "rerun_max_pp4_finite_pairs",
               "replay_max_pp4_all_pairs", "replay_minus_adopted", "adopted_reproduced_before_drop"],
              changes)
    write_tsv(os.path.join(root, "qc/offline_crosscheck_disagreements.tsv"),
              ["level", "gwas_name", "ensembl", "field", "offline", "rerun"], bad)
    write_tsv(os.path.join(root, "qc/offline_hit_label_differences.tsv"),
              ["level", "gwas_name", "ensembl", "field", "offline", "rerun"], label_diffs)
    if fail or bad:
        with open(os.path.join(root, "ASSEMBLY_FAILED.txt"), "w") as fh:
            fh.write("\n".join(fail + [f"offline disagreements: {len(bad)} (qc/offline_crosscheck_disagreements.tsv)"]) + "\n")
        print(f"*** {len(fail)} rerun validation failures, {len(bad)} offline disagreements; no master written")
        for f in fail[:20]:
            print(f"    - {f}")
        sys.exit(2)

    # Second pass: write the r2 master. Unplanned lines are copied verbatim.
    abf_max = {f: 0.0 for f in ABF_FIELDS}
    anchor_diff = collections.Counter()
    lam_max = 0.0
    n_reason = collections.Counter()
    out_master = os.path.join(root, "susie_coloc_all_gwas.csv")
    with open(args.master, newline="") as fin, open(out_master, "w", newline="") as fout:
        header_line = fin.readline()
        fout.write(header_line.rstrip("\r\n") + "," + ",".join(NEW_FIELDS) + "\n")
        for line in fin:
            row = next(csv.reader([line]))
            if len(row) != len(header):
                sys.exit(f"adopted master line has {len(row)} fields, expected {len(header)}")
            key = (row[col["gwas_name"]], row[col["chr"]], row[col["ensembl"]])
            if key not in planned:
                method, pp4 = row[col["method"]], num(row[col["PP.H4.susie"]])
                if method == "susie" and pp4 is not None and pp4 > 0.5:
                    sys.exit(f"SuSiE-positive row missing from the rerun plan: {key}")
                reason = "not_rerun_negative" if method == "susie" else "not_rerun_no_susie_result"
                n_reason[reason] += 1
                fout.write(line.rstrip("\r\n") + "," + "," + reason + "\n")
                continue
            r = rows[key]
            for f in ABF_FIELDS:
                abf_max[f] = max(abf_max[f], absdiff(r[f], row[col[f]]))
            for f in ("n_snps", "top_snp"):
                if r[f] != row[col[f]]:
                    anchor_diff[f] += 1
            if absdiff(r["top_snp_PP"], row[col["top_snp_PP"]]) > ABF_TOL:
                anchor_diff["top_snp_PP"] += 1
            a_lam, r_lam = num(row[col["lambda_s_locus"]]), num(r["lambda_s_locus"])
            if a_lam is not None and r_lam is not None:
                lam_max = max(lam_max, abs(a_lam - r_lam))
            for f in SUSIE_FIELDS:
                row[col[f]] = r[f]
            n_reason["rerun_shared_posterior_check"] += 1
            buf = io.StringIO()
            csv.writer(buf, lineterminator="\n").writerow(row + [r["n_cs_pairs_eligible"],
                                                                "rerun_shared_posterior_check"])
            fout.write(buf.getvalue())

    repro = [{"field": f, "max_abs_diff": abf_max[f], "tolerance": ABF_TOL,
              "pass": str(abf_max[f] <= ABF_TOL).upper()} for f in ABF_FIELDS]
    repro += [{"field": f, "max_abs_diff": f"{anchor_diff[f]} rows differ", "tolerance": "0 rows",
               "pass": str(anchor_diff[f] == 0).upper()} for f in ("n_snps", "top_snp", "top_snp_PP")]
    repro.append({"field": "lambda_s_locus (adopted value kept)", "max_abs_diff": lam_max,
                  "tolerance": "report only", "pass": "NA"})
    write_tsv(os.path.join(root, "qc/reproduction_check.tsv"),
              ["field", "max_abs_diff", "tolerance", "pass"], repro)
    if any(r["pass"] == "FALSE" for r in repro):
        os.replace(out_master, out_master + ".FAILED")
        sys.exit("ABF or anchor fields do not reproduce the adopted master; see qc/reproduction_check.tsv")

    with gzip.open(os.path.join(root, "susie_coloc_signal_pairs.tsv.gz"), "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=pair_header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for key in sorted(pairs):
            w.writerows(pairs[key])

    ld_src = os.path.join(os.path.dirname(args.master), "ld_contamination_clusters.csv")
    ld_dst = os.path.join(root, "ld_contamination_clusters.csv")
    shutil.copyfile(ld_src, ld_dst)
    if sha256(ld_src) != sha256(ld_dst):
        sys.exit("ld_contamination_clusters.csv copy differs from its source")

    status = collections.Counter(r["r2_status"] for r in log)
    status_t12 = collections.Counter(r["r2_status"] for r in log if r["tier"] in ("1", "2"))
    manifest = [{"role": role, "path": p, "sha256": sha256(p)} for role, p in (
        ("adopted_master", args.master),
        ("rerun_plan", os.path.join(args.rerun_root, "plan/rerun_plan.tsv")),
        ("offline_gene_study", os.path.join(args.offline_root, "eqtl_overlap_check_gene_study.tsv")),
        ("offline_signal_pairs", os.path.join(args.offline_root, "eqtl_overlap_check.tsv")),
        ("ld_contamination_clusters", ld_src),
        ("assembler", os.path.abspath(__file__)),
        ("r2_master", out_master))]
    write_tsv(os.path.join(root, "input_manifest.tsv"), ["role", "path", "sha256"], manifest)
    done = {"created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "rerun_root": args.rerun_root, "planned_rows": len(planned),
            "status_all_planned": dict(status), "status_tier12": dict(status_t12),
            "master_rows_by_eligibility_source": dict(n_reason),
            "abf_max_abs_diff": max(abf_max.values()),
            "offline_gene_studies_checked": n_gene_off, "offline_signal_pairs_checked": n_pair_off,
            "offline_corrected_pp4_max_abs_diff": max_off_pp4,
            "batch_logs_scanned": args.batch_jobs.split(","),
            "coloc_overlap_too_small_warnings": n_overlap_warnings,
            "eligibility_rules": {
                "signal_pairs_rules_differ": len(rule_pairs),
                "gene_studies_call_differs": [f"{r['gwas_name']} {r['gene']}: {r['status_coloc_rule']} "
                                              f"(coloc rule) vs {r['status_uniform_rule']} (uniform)"
                                              for r in rule_rows if r["call_differs"] == "TRUE"]},
            "call_changes": {"n": len(changes),
                             "replay_reproduced_adopted": sum(c["adopted_reproduced_before_drop"] == "TRUE"
                                                              for c in changes),
                             "not_in_replay": sum(c["adopted_reproduced_before_drop"].startswith("not_in")
                                                  for c in changes)}}
    with open(os.path.join(root, "ASSEMBLY_COMPLETE.json"), "w") as fh:
        json.dump(done, fh, indent=2)
    print(json.dumps(done, indent=2))


# ----------------------------------------------------------------- summarize
def cascade_sets(master_path, main):
    """Figure 2A definitions (fig1_gwas_creative_options.py): gene SYMBOL, MAIN strata, > 0.5."""
    bs, ba, untestable = {}, {}, set()
    with open(master_path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["gwas_name"] not in main or not r["gene"].strip():
                continue
            if r["method"] == UNTESTABLE:
                untestable.add(r["gene"])
            s, a = num(r["PP.H4.susie"]), num(r["PP.H4.abf"])
            if s is not None and s > bs.get(r["gene"], -1.0):
                bs[r["gene"]] = s
            if a is not None and a > ba.get(r["gene"], -1.0):
                ba[r["gene"]] = a
    multi = {g for g, v in bs.items() if v > 0.5}
    abf = {g for g, v in ba.items() if v > 0.5}
    # A gene whose MAIN SuSiE test is untestable in some study and positive in none
    # is its own class; Fig 2A still draws it as single-signal-only when ABF > 0.5.
    return {"multi_signal": multi, "single_signal_only": abf - multi, "union": multi | abf,
            "susie_untestable": untestable - multi}, bs, ba


def full_portfolio(gene_level_path):
    """All-study SuSiE state per named gene (gene_level_coloc.csv); positive = RESULTS' 745.

    The adopted table predates coloc_susie_state; with no untestable rows its state
    is positive / tested_le_0.5 / no_susie_result from the same columns."""
    state = {}
    with open(gene_level_path, newline="") as fh:
        for r in csv.DictReader(fh):
            if not r["ensembl"].startswith("ENSG"):
                continue
            if "coloc_susie_state" in r:
                state[r["ensembl"]] = r["coloc_susie_state"]
            elif (num(r["coloc_n_gwas_susie_h4_05"]) or 0) > 0:
                state[r["ensembl"]] = "positive"
            else:
                state[r["ensembl"]] = "tested_le_0.5" if num(r["coloc_best_susie_pp4"]) is not None \
                    else "no_susie_result"
    return state


GENE_LEVEL_ADDED = {"coloc_n_gwas_susie_untestable", "coloc_susie_untestable_gwas", "coloc_susie_state"}


def compare_gene_level(a_path, b_path, tol=1e-12):
    """Row-by-row equality of an adopted gene_level_coloc.csv (a) and a rebuild (b), keyed
    by (gene, ensembl), on a's columns; b may add only the untestable-class columns."""
    def load(p):
        with open(p, newline="") as fh:
            rd = csv.DictReader(fh)
            return rd.fieldnames, {(r["gene"], r["ensembl"]): r for r in rd}
    ha, A = load(a_path)
    hb, B = load(b_path)
    diffs = []
    if set(ha) - set(hb) or set(hb) - set(ha) - GENE_LEVEL_ADDED:
        diffs.append(f"columns differ: {sorted(set(ha) ^ set(hb))}")
    if set(A) != set(B):
        diffs.append(f"{len(set(A) ^ set(B))} (gene, ensembl) keys differ")
    for k in set(A) & set(B):
        for c in ha:
            x, y = A[k][c], B[k][c]
            if x == y:
                continue
            fx, fy = num(x), num(y)
            if fx is None or fy is None or abs(fx - fy) > tol:
                diffs.append(f"{k} {c}: {x} vs {y}")
                break
    return diffs


def summarize(args):
    root = release_path(args.release_root)
    if not os.path.isfile(os.path.join(root, "ASSEMBLY_COMPLETE.json")):
        sys.exit("assemble has not completed for this release root")
    out = os.path.join(root, "summary")
    if os.path.lexists(out):
        sys.exit(f"refusing to reuse {out}")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import rebuild_tier12 as rt

    old_master = os.path.join(ADOPTED, "susie_coloc_all_gwas.csv")
    new_master = os.path.join(root, "susie_coloc_all_gwas.csv")
    fail = []

    # (1) The master-in combine path reproduces the adopted gene-level table.
    diffs = compare_gene_level(os.path.join(ADOPTED, "gene_level_coloc.csv"),
                               os.path.join(root, "qc/adopted_gene_level_rerun/gene_level_coloc.csv"))
    print(f"(1) combine from adopted master vs adopted gene_level_coloc.csv: {len(diffs)} differences")
    if diffs:
        fail.append(f"combine path does not reproduce the adopted gene-level table: {diffs[:3]}")

    tier, trait, placement, label = rt.load_tiers()
    main = {s for s, p in placement.items() if p == "main"}
    old_best, _ = rt.build(old_master, tier, trait, placement)
    new_best, _ = rt.build(new_master, tier, trait, placement)

    # (2) The written r2 tier-1/2 table equals the in-memory rebuild.
    with open(os.path.join(root, "gene_level_coloc_tier12.csv"), newline="") as fh:
        t12 = {r["ensembl"]: r for r in csv.DictReader(fh)}
    mism = [g for g, e in new_best.items() if g not in t12
            or num(t12[g]["coloc_best_susie_pp4"]) != e["s"] or num(t12[g]["coloc_best_abf_pp4"]) != e["a"]]
    print(f"(2) r2 gene_level_coloc_tier12.csv vs rebuild: {len(mism)} differences")
    if mism or len(t12) != len(new_best):
        fail.append(f"r2 tier-1/2 table differs from the rebuild for {len(mism)} genes")

    def t12_sets(best):
        s = {g for g, e in best.items() if e["s"] is not None and e["s"] > 0.5}
        a = {g for g, e in best.items() if e["a"] is not None and e["a"] > 0.5}
        return {"tier12_susie_gt_0.5": s, "tier12_abf_gt_0.5": a, "tier12_union_gt_0.5": s | a}

    old_t12, new_t12 = t12_sets(old_best), t12_sets(new_best)
    old_cas, old_bs, old_ba = cascade_sets(old_master, main)
    new_cas, new_bs, new_ba = cascade_sets(new_master, main)
    old_state_fp = full_portfolio(os.path.join(ADOPTED, "gene_level_coloc.csv"))
    new_state_fp = full_portfolio(os.path.join(root, "gene_level_coloc.csv"))
    old_fp = {g for g, st in old_state_fp.items() if st == "positive"}
    new_fp = {g for g, st in new_state_fp.items() if st == "positive"}

    # (3) The definitions reproduce the adopted published counts.
    expected_old = {"tier12_susie_gt_0.5": 462, "tier12_abf_gt_0.5": 822, "tier12_union_gt_0.5": 1013,
                    "multi_signal": 462, "single_signal_only": 551, "union": 1013,
                    "full_portfolio_susie_gt_0.5": 745}
    observed_old = {**{k: len(v) for k, v in old_t12.items()},
                    **{k: len(v) for k, v in old_cas.items()},
                    "full_portfolio_susie_gt_0.5": len(old_fp)}
    for k, v in expected_old.items():
        if observed_old[k] != v:
            fail.append(f"adopted {k} = {observed_old[k]}, published {v}: definition not reproduced")
    if old_cas["susie_untestable"] or any(e["u"] for e in old_best.values()):
        fail.append("the adopted release already carries untestable SuSiE rows")
    n_s_old, c_s_old = rt.class_split(old_best, label, "s", "s_g")
    if (n_s_old, c_s_old.get("liver_enzyme"), c_s_old.get("direct_MASLD")) != (462, 448, 14):
        fail.append(f"adopted SuSiE trait split {n_s_old} {dict(c_s_old)} is not 462 = 448 + 14")
    if fail:
        os.makedirs(out)
        with open(os.path.join(out, "SUMMARY_FAILED.txt"), "w") as fh:
            fh.write("\n".join(fail) + "\n")
        print("*** summary checks failed:\n    - " + "\n    - ".join(fail))
        sys.exit(2)

    os.makedirs(out)
    n_s_new, c_s_new = rt.class_split(new_best, label, "s", "s_g")
    n_a_old, c_a_old = rt.class_split(old_best, label, "a", "a_g")
    n_a_new, c_a_new = rt.class_split(new_best, label, "a", "a_g")

    counts = [("tier12_susie_gt_0.5", old_t12, new_t12), ("tier12_abf_gt_0.5", old_t12, new_t12),
              ("tier12_union_gt_0.5", old_t12, new_t12), ("multi_signal", old_cas, new_cas),
              ("single_signal_only", old_cas, new_cas), ("union", old_cas, new_cas),
              ("susie_untestable", old_cas, new_cas)]
    table = [{"metric": k, "unit": "ensembl gene" if k.startswith("tier12") else "gene symbol (Fig 2A)",
              "adopted": len(o[k]), "r2": len(n[k]), "delta": len(n[k]) - len(o[k]),
              "lost": len(o[k] - n[k]), "gained": len(n[k] - o[k])} for k, o, n in counts]
    table.append({"metric": "full_portfolio_susie_gt_0.5", "unit": "ensembl gene, all 50 studies",
                  "adopted": len(old_fp), "r2": len(new_fp), "delta": len(new_fp) - len(old_fp),
                  "lost": len(old_fp - new_fp), "gained": len(new_fp - old_fp)})
    for name, n_old, c_old, n_new, c_new in (("susie_ge_0.5_by_susie_driver", n_s_old, c_s_old, n_s_new, c_s_new),
                                             ("abf_ge_0.5_by_abf_driver", n_a_old, c_a_old, n_a_new, c_a_new)):
        for cls in ("liver_enzyme", "direct_MASLD"):
            table.append({"metric": f"{name}:{cls}", "unit": "ensembl gene",
                          "adopted": c_old.get(cls, 0), "r2": c_new.get(cls, 0),
                          "delta": c_new.get(cls, 0) - c_old.get(cls, 0), "lost": "", "gained": ""})
        table.append({"metric": f"{name}:total", "unit": "ensembl gene", "adopted": n_old, "r2": n_new,
                      "delta": n_new - n_old, "lost": "", "gained": ""})
    # Untestable is its own class, never folded into the negatives.
    old_st12 = collections.Counter(rt.susie_state(e) for e in old_best.values())
    new_st12 = collections.Counter(rt.susie_state(e) for e in new_best.values())
    old_stfp, new_stfp = collections.Counter(old_state_fp.values()), collections.Counter(new_state_fp.values())
    for name, o, n, unit in (("tier12_susie_state", old_st12, new_st12, "ensembl gene, tier-1/2 studies"),
                             ("full_portfolio_susie_state", old_stfp, new_stfp, "ensembl gene, all 50 studies")):
        for st in ("positive", "untestable", "tested_le_0.5", "no_susie_result"):
            table.append({"metric": f"{name}:{st}", "unit": unit, "adopted": o.get(st, 0),
                          "r2": n.get(st, 0), "delta": n.get(st, 0) - o.get(st, 0), "lost": "", "gained": ""})
    write_tsv(os.path.join(out, "old_vs_new_counts.tsv"),
              ["metric", "unit", "adopted", "r2", "delta", "lost", "gained"], table)

    # Per-gene membership changes.
    changes = []
    for k, o, n in counts:
        is_t12 = k.startswith("tier12")
        for g in sorted(o[k] ^ n[k]):
            ob, nb = (old_best.get(g), new_best.get(g)) if is_t12 else (None, None)
            changes.append({"set": k, "gene_or_ensembl": g, "change": "lost" if g in o[k] else "gained",
                            "adopted_susie_pp4": ob["s"] if ob else old_bs.get(g, ""),
                            "r2_susie_pp4": nb["s"] if nb else new_bs.get(g, ""),
                            "adopted_abf_pp4": ob["a"] if ob else old_ba.get(g, ""),
                            "r2_abf_pp4": nb["a"] if nb else new_ba.get(g, ""),
                            "adopted_susie_driver": ob["s_g"] if ob else "",
                            "r2_susie_driver": nb["s_g"] if nb else "",
                            "r2_susie_state": rt.susie_state(nb) if nb else
                            ("untestable" if g in new_cas["susie_untestable"] else
                             "positive" if g in new_cas["multi_signal"] else "not_positive")})
    for g in sorted(old_fp ^ new_fp):
        changes.append({"set": "full_portfolio_susie_gt_0.5", "gene_or_ensembl": g,
                        "change": "lost" if g in old_fp else "gained",
                        "r2_susie_state": new_state_fp.get(g, "")})
    write_tsv(os.path.join(out, "gene_membership_changes.tsv"),
              ["set", "gene_or_ensembl", "change", "adopted_susie_pp4", "r2_susie_pp4",
               "adopted_abf_pp4", "r2_abf_pp4", "adopted_susie_driver", "r2_susie_driver",
               "r2_susie_state"], changes)

    # Named genes: tier-1/2 bests by method, plus every rerun gene-study with its SuSiE anchor.
    log = read_tsv(os.path.join(root, "qc/replacement_log.tsv"))
    with gzip.open(os.path.join(root, "susie_coloc_signal_pairs.tsv.gz"), "rt", newline="") as fh:
        pairs = collections.defaultdict(list)
        for p in csv.DictReader(fh, delimiter="\t"):
            if p["gene"] in NAMED_GENES:
                pairs[(p["gwas"], p["ensembl"])].append(p)
    # Adopted anchor = the replay pair with the highest PP.H4 over ALL pairs (the
    # rerun reports NA for ineligible pairs); tier-1/2 gene-studies only.
    replay = collections.defaultdict(list)
    for p in read_tsv(os.path.join(OFFLINE, "eqtl_overlap_check.tsv")):
        if p["gene"] in NAMED_GENES:
            replay[(p["gwas_name"], p["ensembl"])].append(p)
    ens_of = collections.defaultdict(set)
    with open(new_master, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["gene"] in NAMED_GENES and r["ensembl"].startswith("ENSG"):
                ens_of[r["gene"]].add(r["ensembl"])
    named = []
    for gene in NAMED_GENES:
        for ens in sorted(ens_of.get(gene, {""})):
            ob, nb = old_best.get(ens), new_best.get(ens)
            named.append({"gene": gene, "ensembl": ens, "row": "tier12_best",
                          "gwas_name": "", "adopted_susie": f"{ob['s_g']} {ob['s']}" if ob and ob["s"] is not None else "",
                          "r2_susie": f"{nb['s_g']} {nb['s']}" if nb and nb["s"] is not None else "",
                          "adopted_abf": f"{ob['a_g']} {ob['a']}" if ob and ob["a"] is not None else "",
                          "r2_abf": f"{nb['a_g']} {nb['a']}" if nb and nb["a"] is not None else "",
                          "r2_status": "", "n_cs_pairs": "", "n_cs_pairs_eligible": "",
                          "adopted_anchor_gwas_eqtl": "", "r2_anchor_gwas_eqtl": ""})
            for r in log:
                if r["ensembl"] != ens:
                    continue
                prs = pairs.get((r["gwas_name"], ens), [])
                best_all = max(replay.get((r["gwas_name"], ens), []),
                               key=lambda p: num(p["replay_PP.H4"]), default=None)
                elig = [p for p in prs if p["eligible"] == "TRUE"]
                best_el = max(elig, key=lambda p: num(p["PP.H4"]), default=None)
                named.append({"gene": gene, "ensembl": ens, "row": "rerun_gene_study",
                              "gwas_name": r["gwas_name"], "adopted_susie": r["adopted_pp_h4_susie"],
                              "r2_susie": r["r2_pp_h4_susie"], "adopted_abf": "", "r2_abf": "",
                              "r2_status": r["r2_status"], "n_cs_pairs": r["n_cs_pairs"],
                              "n_cs_pairs_eligible": r["n_cs_pairs_eligible"],
                              "adopted_anchor_gwas_eqtl": f"{best_all['hit1']} {best_all['hit2']}" if best_all else "",
                              "r2_anchor_gwas_eqtl": f"{best_el['hit1']} {best_el['hit2']}" if best_el else ""})
    write_tsv(os.path.join(out, "named_gene_status.tsv"), list(named[0]), named)

    worked = {}
    with open(new_master, newline="") as fh:
        for r in csv.DictReader(fh):
            if WORKED_LOCI.get(r["gene"]) == r["gwas_name"]:
                worked[r["gene"]] = {"gwas_name": r["gwas_name"], "ensembl": r["ensembl"],
                                     "PP.H4.susie": num(r["PP.H4.susie"]), "method": r["method"],
                                     "top_snp_PP": num(r["top_snp_PP"])}
    status = json.load(open(os.path.join(root, "ASSEMBLY_COMPLETE.json")))
    summary = {
        "release_root": root, "adopted_root": ADOPTED,
        "fig2": {"multi_signal": len(new_cas["multi_signal"]),
                 "single_signal_only": len(new_cas["single_signal_only"]),
                 "multi_or_single_signal_coloc_gene_union": len(new_cas["union"])},
        "fig2_adopted": {"multi_signal": len(old_cas["multi_signal"]),
                         "single_signal_only": len(old_cas["single_signal_only"]),
                         "multi_or_single_signal_coloc_gene_union": len(old_cas["union"])},
        "fig2_susie_untestable": {
            "n": len(new_cas["susie_untestable"]),
            "drawn_as_single_signal_only": len(new_cas["susie_untestable"] & new_cas["single_signal_only"]),
            "genes": sorted(new_cas["susie_untestable"])},
        "tier12": {k: len(v) for k, v in new_t12.items()},
        "tier12_susie_state": dict(new_st12), "tier12_susie_state_adopted": dict(old_st12),
        "full_portfolio_susie_gt_0.5": len(new_fp),
        "full_portfolio_susie_state": dict(new_stfp), "full_portfolio_susie_state_adopted": dict(old_stfp),
        "eligibility_rules": status["eligibility_rules"], "call_changes": status["call_changes"],
        "trait_split": {"susie_ge_0.5_by_susie_driver": {"total": n_s_new, **dict(c_s_new)},
                        "abf_ge_0.5_by_abf_driver": {"total": n_a_new, **dict(c_a_new)}},
        "trait_split_adopted": {"susie_ge_0.5_by_susie_driver": {"total": n_s_old, **dict(c_s_old)},
                                "abf_ge_0.5_by_abf_driver": {"total": n_a_old, **dict(c_a_old)}},
        "worked_loci": worked,
        "gene_study_status_all_planned": status["status_all_planned"],
        "gene_study_status_tier12": status["status_tier12"],
    }
    with open(os.path.join(root, "r2_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    for row in table:
        print(f"  {row['metric']:<42} {row['adopted']:>6} -> {row['r2']:>6}")
    print(json.dumps(summary["worked_loci"], indent=2))
    print(f"wrote {out} and {os.path.join(root, 'r2_summary.json')}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="step", required=True)
    a = sub.add_parser("assemble")
    a.add_argument("--release-root", required=True, help="new results/susie_coloc_r2_<UTC ts>")
    a.add_argument("--rerun-root", default=RERUN)
    a.add_argument("--offline-root", default=OFFLINE)
    a.add_argument("--master", default=os.path.join(ADOPTED, "susie_coloc_all_gwas.csv"))
    a.add_argument("--replay-root", default=REPLAY)
    a.add_argument("--log-dir", default=LOG_DIR)
    a.add_argument("--batch-jobs", default=BATCH_JOBS,
                   help="comma-separated 09b job ids whose logs are scanned; include resubmissions")
    s = sub.add_parser("summarize")
    s.add_argument("--release-root", required=True)
    args = ap.parse_args()
    assemble(args) if args.step == "assemble" else summarize(args)
