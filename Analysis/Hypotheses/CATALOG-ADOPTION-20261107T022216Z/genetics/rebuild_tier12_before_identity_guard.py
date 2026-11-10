#!/usr/bin/env python3
"""
Rebuild gene_level_coloc_tier12.csv -- the Fig2 tier-1/2 headline file.

WHY THIS EXISTS: promoting the rerun to canonical replaced gene_level_coloc.csv
and susie_coloc_all_gwas.csv, but gene_level_coloc_tier12.csv is dated
2026-07-06 and NO script in the repo produces it. It is consumed by
scripts/figures/fig4a_overview_candidates.py,
scripts/figures/fig1_gwas_creative_options.py and
scripts/portal/generate_ancestry_coloc_data.py -- so left alone it would feed a
stale headline into figures and the portal with nothing to signal the mismatch.

METHOD, AND THE GUARD ON IT: the file's semantics are inferred, not documented.
So this script FIRST rebuilds the OLD file from the OLD master table and diffs
it against the real old file. Only if that reproduces does it rebuild from the
requested master. If the reconstruction does not reproduce, it writes nothing
and says so -- shipping a guessed headline file is worse than shipping none.

Tier-1/2 studies come from config/gwas_trait_tier.tsv (15 tier-1 + 20 tier-2,
all placement=main; the 15 tier-3/4 are supp).

METHOD-SPECIFIC DRIVERS (review item 5, 2026-09-23): `driving_gwas` is the best
SuSiE study whenever any tier-1/2 SuSiE PP.H4 is finite, else the best ABF
study. It is kept unchanged for existing consumers, but it is NOT the ABF
driver: GNMT's ABF best is 2021_34841290_NAFLD_EUR (0.705832) while its
driving_gwas is UKBB_ALT from a SuSiE PP.H4 of 0.128219. The best SuSiE and
best ABF studies are now tracked separately (`susie_driving_gwas`,
`abf_driving_gwas` + trait + tier), and an ABF trait-class split must use the
abf_* columns. Tie rule, same for both methods: higher PP.H4 wins; equal PP.H4
-> the lexicographically smaller study id. Blank or non-Ensembl gene ids are
left out of the named table and written to a separate diagnostics file.

UNTESTABLE SuSiE (review item 2): a tier-1/2 row whose every signal pair failed
coloc's shared-posterior check has PP.H4.susie = NA and method
susie_untestable_insufficient_shared_posterior. It is counted in
susie_untestable_t12_studies, and susie_state_t12 is `untestable` unless another
tier-1/2 study is SuSiE-positive; such a gene is never a SuSiE negative.
susie_state_t12: positive (> 0.5) > untestable > tested_le_0.5 > no_susie_result.

Usage:
  python3 rebuild_tier12.py --out-dir results/coloc_tier12_<UTC timestamp> \
      [--master <susie_coloc_all_gwas.csv>]
"""
import argparse, csv, sys, os, re, collections

FM = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
TIERCFG = os.path.join(FM, "config/gwas_trait_tier.tsv")
ADOPTED = os.path.join(FM, "results/susie_coloc")
OLD = os.path.join(FM, "results/susie_coloc_preRerun_backup_2026-08-17")
ENSG = re.compile(r"^ENSG\d{11}$")
UNTESTABLE = "susie_untestable_insufficient_shared_posterior"
HEADER = ["ensembl", "coloc_best_susie_pp4", "coloc_best_abf_pp4", "driving_gwas",
          "driving_trait", "driving_tier", "any_main", "any_supp", "tier34_only",
          "susie_driving_gwas", "abf_driving_gwas", "abf_driving_trait", "abf_driving_tier",
          "susie_untestable_t12_studies", "susie_state_t12"]


def load_tiers():
    tier, trait, placement, label = {}, {}, {}, {}
    with open(TIERCFG) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tier[r["study_name"]] = int(r["tier"])
            trait[r["study_name"]] = r["trait"]
            placement[r["study_name"]] = r["placement"]
            label[r["study_name"]] = r["tier_label"]
    return tier, trait, placement, label


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def better(pp, study, best_pp, best_study):
    """Declared tie rule: higher PP.H4; equal PP.H4 -> smaller study id."""
    return best_pp is None or pp > best_pp or (pp == best_pp and study < best_study)


def build(master_path, tier, trait, placement):
    """Per gene: best SuSiE and best ABF over TIER-1/2 studies only, tracked separately."""
    best = collections.defaultdict(lambda: {"s": None, "s_g": "", "a": None, "a_g": "",
                                            "any_main": False, "any_supp": False,
                                            "any_t12": False, "u": 0})
    unmapped = []
    with open(master_path) as fh:
        rd = csv.DictReader(fh)
        gcol = "gwas_name" if "gwas_name" in rd.fieldnames else "gwas"
        scol = "PP.H4.susie" if "PP.H4.susie" in rd.fieldnames else None
        acol = "PP.H4.abf" if "PP.H4.abf" in rd.fieldnames else None
        for r in rd:
            g = r.get("ensembl", "")
            st = r.get(gcol, "")
            if not ENSG.match(g):
                unmapped.append({"gwas_name": st, "gene": r.get("gene", ""), "ensembl": g,
                                 "chr": r.get("chr", ""), "tier": tier.get(st, ""),
                                 "PP.H4.abf": r.get(acol, "") if acol else "",
                                 "PP.H4.susie": r.get(scol, "") if scol else "",
                                 "method": r.get("method", "")})
                continue
            t = tier.get(st)
            e = best[g]
            if placement.get(st) == "main":
                e["any_main"] = True
            if placement.get(st) == "supp":
                e["any_supp"] = True
            if t in (1, 2):
                e["any_t12"] = True
                s = num(r.get(scol)) if scol else None
                a = num(r.get(acol)) if acol else None
                if s is not None and better(s, st, e["s"], e["s_g"]):
                    e["s"], e["s_g"] = s, st
                if a is not None and better(a, st, e["a"], e["a_g"]):
                    e["a"], e["a_g"] = a, st
                if r.get("method") == UNTESTABLE:
                    e["u"] += 1
    return best, unmapped


def susie_state(e):
    """Tier-1/2 SuSiE evidence state; untestable is never a negative."""
    if e["s"] is not None and e["s"] > 0.5:
        return "positive"
    if e["u"]:
        return "untestable"
    return "tested_le_0.5" if e["s"] is not None else "no_susie_result"


def driving(e):
    """Legacy driver: the SuSiE study if any SuSiE PP.H4 is finite, else the ABF study."""
    return e["s_g"] if e["s"] is not None else e["a_g"]


def write(best, unmapped, tier, trait, out, diag):
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        for g, e in best.items():
            d = driving(e)
            w.writerow([g,
                        "" if e["s"] is None else e["s"],
                        "" if e["a"] is None else e["a"],
                        d, trait.get(d, ""), tier.get(d, ""),
                        str(e["any_main"]).upper(), str(e["any_supp"]).upper(),
                        str(not e["any_t12"]).upper(),
                        e["s_g"], e["a_g"], trait.get(e["a_g"], ""), tier.get(e["a_g"], ""),
                        e["u"], susie_state(e)])
    with open(diag, "w", newline="") as fh:
        cols = ["gwas_name", "gene", "ensembl", "chr", "tier", "PP.H4.abf", "PP.H4.susie", "method"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(unmapped)


def compare(a_path, b_path):
    """Compare on the numeric headline only -- the thing Fig2 actually uses."""
    def load(p):
        d = {}
        for r in csv.DictReader(open(p)):
            d[r["ensembl"]] = (num(r.get("coloc_best_susie_pp4")), num(r.get("coloc_best_abf_pp4")))
        return d
    A, B = load(a_path), load(b_path)
    common = set(A) & set(B)
    sus_gt = lambda d, k: (d[k][0] or 0) > 0.5
    same_s = sum(1 for k in common if sus_gt(A, k) == sus_gt(B, k))
    n5A = sum(1 for k in A if (A[k][0] or 0) > 0.5)
    n5B = sum(1 for k in B if (B[k][0] or 0) > 0.5)
    print(f"    rows: real {len(A):,} | rebuilt {len(B):,} | common {len(common):,}")
    print(f"    SuSiE>0.5: real {n5A} | rebuilt {n5B}")
    print(f"    agreement on the >0.5 call: {same_s:,}/{len(common):,} "
          f"({100*same_s/max(len(common),1):.2f}%)")
    return n5A, n5B, same_s == len(common)


def class_split(best, label, key_pp, key_g):
    """Genes with PP.H4 >= 0.5 by the trait class of the method's own driver."""
    c = collections.Counter(label.get(e[key_g], "") for e in best.values()
                            if e[key_pp] is not None and e[key_pp] >= 0.5)
    return sum(c.values()), c


def regression_checks(tier, trait, placement, label):
    """Review item 5 on the adopted (pre-correction) release. Returns failures."""
    fail = []
    best, _ = build(os.path.join(ADOPTED, "susie_coloc_all_gwas.csv"), tier, trait, placement)
    real = {r["ensembl"]: r for r in csv.DictReader(open(os.path.join(ADOPTED, "gene_level_coloc_tier12.csv")))}
    gl = {r["ensembl"]: r for r in csv.DictReader(open(os.path.join(ADOPTED, "gene_level_coloc.csv")))}

    # (a) The legacy columns still reproduce the adopted file for every named gene.
    named = [g for g in real if ENSG.match(g)]
    diff = [g for g in named if g not in best
            or num(real[g]["coloc_best_susie_pp4"]) != best[g]["s"]
            or num(real[g]["coloc_best_abf_pp4"]) != best[g]["a"]
            or real[g]["driving_gwas"] != driving(best[g])]
    print(f"  (a) legacy columns vs adopted file: {len(named) - len(diff):,}/{len(named):,} named genes identical")
    if diff:
        fail.append(f"legacy columns differ for {len(diff)} genes, e.g. {diff[:3]}")

    # (b) The review's reconstruction: ABF PP.H4 >= 0.5 records whose all-study best ABF
    # study (gene_level_coloc.csv coloc_best_gwas) is itself tier 1/2, so it is that
    # record's tier-1/2 ABF driver. Classes of the stored driving_gwas vs that study.
    rec = [g for g, r in real.items() if (num(r["coloc_best_abf_pp4"]) or 0) >= 0.5
           and g in gl and tier.get(gl[g]["coloc_best_gwas"]) in (1, 2)]
    flips = collections.Counter((label.get(real[g]["driving_gwas"]), label.get(gl[g]["coloc_best_gwas"]))
                                for g in rec
                                if label.get(real[g]["driving_gwas"]) != label.get(gl[g]["coloc_best_gwas"]))
    e2d = flips[("liver_enzyme", "direct_MASLD")]
    d2e = flips[("direct_MASLD", "liver_enzyme")]
    print(f"  (b) recoverable ABF-positive records: {len(rec)} | class disagreements: "
          f"{sum(flips.values())} (enzyme->direct {e2d}, direct->enzyme {d2e})")
    if (len(rec), e2d, d2e, sum(flips.values())) != (772, 52, 10, 62):
        fail.append("review reconstruction is not 772 records / 62 = 52 + 10 disagreements")

    # (c) The new ABF driver equals that recovered study on every recoverable named gene.
    bad = [g for g in rec if ENSG.match(g) and best[g]["a_g"] != gl[g]["coloc_best_gwas"]]
    print(f"  (c) abf_driving_gwas == recovered ABF study: "
          f"{sum(1 for g in rec if ENSG.match(g)) - len(bad)}/{sum(1 for g in rec if ENSG.match(g))}")
    if bad:
        fail.append(f"abf_driving_gwas differs from the recovered ABF study for {len(bad)} genes")

    # (d) High ABF, low finite SuSiE: GNMT.
    e = best["ENSG00000124713"]
    print(f"  (d) GNMT: abf {e['a_g']} {e['a']:.6f} | susie {e['s_g']} {e['s']:.6f} | "
          f"driving_gwas {driving(e)}")
    if (e["a_g"], round(e["a"], 6), e["s_g"], round(e["s"], 6)) != \
            ("2021_34841290_NAFLD_EUR", 0.705832, "UKBB_ALT", 0.128219):
        fail.append("GNMT regression case does not reproduce")

    # (e) Complete adopted denominators, each method split by its own driver.
    n_a, c_a = class_split(best, label, "a", "a_g")
    n_s, c_s = class_split(best, label, "s", "s_g")
    legacy = collections.Counter((label.get(driving(e)), label.get(e["a_g"])) for e in best.values()
                                 if e["a"] is not None and e["a"] >= 0.5
                                 and label.get(driving(e)) != label.get(e["a_g"]))
    print(f"  (e) adopted, named genes: ABF >= 0.5 {n_a} by ABF driver {dict(c_a)} | "
          f"SuSiE >= 0.5 {n_s} by SuSiE driver {dict(c_s)}")
    print(f"      legacy driving class vs ABF driver class, all ABF >= 0.5 named genes: "
          f"{sum(legacy.values())} disagree {dict(legacy)}")
    return fail


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="new directory; refuses to reuse one")
    ap.add_argument("--master", default=os.path.join(ADOPTED, "susie_coloc_all_gwas.csv"))
    args = ap.parse_args()
    out_dir = args.out_dir if os.path.isabs(args.out_dir) else os.path.join(FM, args.out_dir)
    if os.path.lexists(out_dir):
        sys.exit(f"refusing to reuse {out_dir}")
    os.makedirs(out_dir)

    tier, trait, placement, label = load_tiers()
    print(f"tier-1/2 studies: {sum(1 for v in tier.values() if v in (1,2))} of {len(tier)}")

    OLD_MASTER = os.path.join(OLD, "susie_coloc_all_gwas.csv")
    OLD_TIER12 = os.path.join(OLD, "gene_level_coloc_tier12.csv")

    print("\n=== VALIDATION: rebuild the OLD file from the OLD master ===")
    old_rebuilt = os.path.join(out_dir, "tier12_old_rebuilt.csv")
    old_best, old_unmapped = build(OLD_MASTER, tier, trait, placement)
    write(old_best, old_unmapped, tier, trait, old_rebuilt,
          os.path.join(out_dir, "tier12_old_rebuilt_unmapped_ids.csv"))
    n_real, n_reb, exact = compare(OLD_TIER12, old_rebuilt)

    if not exact:
        print("\n*** RECONSTRUCTION DOES NOT REPRODUCE THE OLD FILE ***")
        print("    Writing NOTHING. The semantics of this file are not recoverable")
        print("    from the repo, so it must be rebuilt by whoever owns its producer.")
        sys.exit(2)

    print("\n=== REGRESSION: method-specific drivers on the adopted release ===")
    failures = regression_checks(tier, trait, placement, label)
    if failures:
        print("\n*** REGRESSION CHECKS FAILED; writing no tier-1/2 table ***")
        for f in failures:
            print(f"    - {f}")
        sys.exit(3)

    print(f"\n=== rebuilding from {args.master} ===")
    best, unmapped = build(args.master, tier, trait, placement)
    out = os.path.join(out_dir, "gene_level_coloc_tier12.csv")
    diag = os.path.join(out_dir, "gene_level_coloc_tier12_unmapped_ids.csv")
    write(best, unmapped, tier, trait, out, diag)
    n_a, c_a = class_split(best, label, "a", "a_g")
    n_s, c_s = class_split(best, label, "s", "s_g")
    print(f"    wrote {out} ({len(best):,} named genes)")
    print(f"    wrote {diag} ({len(unmapped):,} master rows without an Ensembl id)")
    print(f"    ABF PP.H4 >= 0.5: {n_a} genes by ABF driver {dict(c_a)}")
    print(f"    SuSiE PP.H4 >= 0.5: {n_s} genes by SuSiE driver {dict(c_s)}")
    print(f"    tier-1/2 SuSiE state: {dict(collections.Counter(susie_state(e) for e in best.values()))}")
    compare(os.path.join(ADOPTED, "gene_level_coloc_tier12.csv"), out)
