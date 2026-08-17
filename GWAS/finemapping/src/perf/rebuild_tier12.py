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
new one. If the reconstruction does not reproduce, it writes nothing and says
so -- shipping a guessed headline file is worse than shipping none.

Tier-1/2 studies come from config/gwas_trait_tier.tsv (15 tier-1 + 20 tier-2,
all placement=main; the 15 tier-3/4 are supp).
"""
import csv, sys, os, collections

FM = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
TIERCFG = os.path.join(FM, "config/gwas_trait_tier.tsv")


def load_tiers():
    tier, trait, placement = {}, {}, {}
    with open(TIERCFG) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tier[r["study_name"]] = int(r["tier"])
            trait[r["study_name"]] = r["trait"]
            placement[r["study_name"]] = r["placement"]
    return tier, trait, placement


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def build(master_path, tier, trait, placement):
    """Per gene: best SuSiE and best ABF over TIER-1/2 studies only."""
    best = collections.defaultdict(lambda: {"s": None, "a": None, "g": "", "any_main": False,
                                            "any_supp": False, "any_t12": False})
    with open(master_path) as fh:
        rd = csv.DictReader(fh)
        gcol = "gwas_name" if "gwas_name" in rd.fieldnames else "gwas"
        scol = "PP.H4.susie" if "PP.H4.susie" in rd.fieldnames else None
        acol = "PP.H4.abf" if "PP.H4.abf" in rd.fieldnames else None
        for r in rd:
            g = r.get("ensembl", "")
            st = r.get(gcol, "")
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
                if s is not None and (e["s"] is None or s > e["s"]):
                    e["s"] = s; e["g"] = st
                if a is not None and (e["a"] is None or a > e["a"]):
                    e["a"] = a
                    if e["s"] is None:
                        e["g"] = st
    return best


def write(best, tier, trait, out):
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ensembl", "coloc_best_susie_pp4", "coloc_best_abf_pp4", "driving_gwas",
                    "driving_trait", "driving_tier", "any_main", "any_supp", "tier34_only"])
        for g, e in best.items():
            w.writerow([g,
                        "" if e["s"] is None else e["s"],
                        "" if e["a"] is None else e["a"],
                        e["g"], trait.get(e["g"], ""), tier.get(e["g"], ""),
                        str(e["any_main"]).upper(), str(e["any_supp"]).upper(),
                        str(not e["any_t12"]).upper()])


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


if __name__ == "__main__":
    tier, trait, placement = load_tiers()
    print(f"tier-1/2 studies: {sum(1 for v in tier.values() if v in (1,2))} of {len(tier)}")

    OLD_MASTER = os.path.join(FM, "results/susie_coloc_preRerun_backup_2026-08-17/susie_coloc_all_gwas.csv")
    OLD_TIER12 = os.path.join(FM, "results/susie_coloc_preRerun_backup_2026-08-17/gene_level_coloc_tier12.csv")
    NEW_MASTER = os.path.join(FM, "results/susie_coloc/susie_coloc_all_gwas.csv")
    NEW_TIER12 = os.path.join(FM, "results/susie_coloc/gene_level_coloc_tier12.csv")

    print("\n=== VALIDATION: rebuild the OLD file from the OLD master ===")
    old_rebuilt = "/tmp/tier12_old_rebuilt.csv"
    write(build(OLD_MASTER, tier, trait, placement), tier, trait, old_rebuilt)
    n_real, n_reb, exact = compare(OLD_TIER12, old_rebuilt)

    if not exact:
        print("\n*** RECONSTRUCTION DOES NOT REPRODUCE THE OLD FILE ***")
        print("    Writing NOTHING. The semantics of this file are not recoverable")
        print("    from the repo, so it must be rebuilt by whoever owns its producer.")
        sys.exit(2)

    print("\n=== reconstruction validated; rebuilding from the PROMOTED master ===")
    write(build(NEW_MASTER, tier, trait, placement), tier, trait, NEW_TIER12)
    print(f"    wrote {NEW_TIER12}")
    compare(OLD_TIER12, NEW_TIER12)
