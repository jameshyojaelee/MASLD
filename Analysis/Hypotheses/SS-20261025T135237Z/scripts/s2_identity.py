#!/usr/bin/env python3
"""SS job S2, identity part (PRESPEC_SS section 3, Q2 allele rule and Q3). Reads no stage and no expression.

Inputs: the per-library SteatoSITE VCFs written by s1_one.sh (DIG1 job 62 flags, all S0d sites), the S0d annotated
site list (role "identity" = A1 identity panel minus C>T and G>A), and the 10 cohort VCFs of the A1 identity run
(1,466 libraries, a1-full-20260924T025858Z).
Steps:
  1. SteatoSITE: keep identity sites, set GT missing where DP < 10, merge libraries.
  2. Existing cohorts: the same sites, merged with step 1; biallelic SNVs only; REF/ALT must match the panel.
  3. plink2 KING-robust kinship; a pair is called only with >= 1,000 sites called in both (nsnp).
  4. Q2 allele rule: per SteatoSITE library, pooled other-allele read fraction at sites called homozygous with
     DP >= 20 (reads not of the called allele / all reads).
Outputs (in --out): kinship_pairs.tsv (kinship >= 0.0884, nsnp >= 1,000), ss_called_sites.tsv,
ss_hom_other_allele.tsv, identity_summary.json. VCFs stay in --out/vcf (genotypes; not for sharing).
Usage: python -I s2_identity.py --ss-vcf-dir D --targets-annotated F --a1-vcf-dir D --libraries F --out NEW_DIR
"""
import argparse
import json
import os
import subprocess
from pathlib import Path

import pandas as pd

BCF = "/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools"
COHORTS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE167523", "GSE174478", "GSE193066",
           "GSE213621", "GSE240729", "PRJNA512027"]
MIN_DP_CALL = 10
MIN_DP_HOM = 20
MIN_SHARED = 1000


def run(cmd, **kw):
    print("+", " ".join(map(str, cmd)), flush=True)
    return subprocess.run(list(map(str, cmd)), check=True, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ss-vcf-dir", required=True)
    ap.add_argument("--targets-annotated", required=True)
    ap.add_argument("--a1-vcf-dir", required=True)
    ap.add_argument("--libraries", required=True, help="one sample_alias per line")
    ap.add_argument("--out", required=True)
    ap.add_argument("--plink2", default=os.environ.get("PLINK2", "plink2"))
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"refusing to overwrite {out}")
    (out / "vcf").mkdir(parents=True)

    t = pd.read_csv(a.targets_annotated, sep="\t", dtype=str)
    ident = t[t["role"] == "identity"][["chrom", "pos", "ref", "alt"]].drop_duplicates()
    bad = ident[((ident.ref == "C") & (ident.alt == "T")) | ((ident.ref == "G") & (ident.alt == "A"))]
    if len(bad):
        raise SystemExit("identity panel contains C>T or G>A sites")
    ident = ident.assign(pos=ident.pos.astype(int)).sort_values(["chrom", "pos"])
    regions = out / "identity_sites.tsv"
    ident[["chrom", "pos"]].to_csv(regions, sep="\t", header=False, index=False)
    panel = {(r.chrom, r.pos): (r.ref, r.alt) for r in ident.itertuples()}

    libs = [x.strip() for x in open(a.libraries) if x.strip()]
    missing = [s for s in libs if not (Path(a.ss_vcf_dir) / f"{s}.vcf.gz").exists()]
    if missing:
        raise SystemExit(f"{len(missing)} SteatoSITE VCFs missing, e.g. {missing[:3]}")

    # Step 1 and the Q2 allele rule, one library at a time (each VCF is small).
    tmp = out / "vcf" / "ss_tmp"
    tmp.mkdir()
    called_rows, hom_rows = [], []
    for s in libs:
        src = Path(a.ss_vcf_dir) / f"{s}.vcf.gz"
        q = run([BCF, "query", "-R", regions, "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\t[%DP]\n", src],
                capture_output=True, text=True).stdout
        n_called = n_hom = other = depth = 0
        for line in q.splitlines():
            chrom, pos, ref, alt, gt, ad, dp = line.split("\t")
            if (chrom, int(pos)) not in panel or dp in (".", "") or int(dp) < MIN_DP_CALL or "." in gt:
                continue
            n_called += 1
            if int(dp) < MIN_DP_HOM:
                continue
            al = gt.replace("|", "/").split("/")
            if len(al) != 2 or al[0] != al[1]:
                continue
            ads = [int(x) for x in ad.split(",") if x != "."]
            k = int(al[0])
            if k >= len(ads) or sum(ads) == 0:
                continue
            n_hom += 1
            other += sum(ads) - ads[k]
            depth += sum(ads)
        called_rows.append({"sample_alias": s, "identity_sites_called_dp10": n_called})
        hom_rows.append({"sample_alias": s, "hom_sites_dp20": n_hom, "other_allele_reads": other, "reads": depth,
                         "other_allele_fraction": other / depth if depth else float("nan")})
        f = tmp / f"{s}.vcf.gz"
        p1 = subprocess.Popen([BCF, "view", "-R", str(regions), "-Ou", str(src)], stdout=subprocess.PIPE)
        run([BCF, "+setGT", "-", "-Oz", "-o", f, "--", "-t", "q", "-n", ".", "-i", f"FMT/DP<{MIN_DP_CALL}"], stdin=p1.stdout)
        p1.wait()
        run([BCF, "index", "-t", f])
    pd.DataFrame(called_rows).to_csv(out / "ss_called_sites.tsv", sep="\t", index=False)
    pd.DataFrame(hom_rows).to_csv(out / "ss_hom_other_allele.tsv", sep="\t", index=False)

    lst = out / "vcf" / "ss_list.txt"
    lst.write_text("".join(f"{tmp / (s + '.vcf.gz')}\n" for s in libs))
    ss = out / "vcf" / "steatosite.vcf.gz"
    run([BCF, "merge", "-m", "snps", "--threads", a.threads, "-l", lst, "-Oz", "-o", ss])
    run([BCF, "index", "-t", ss])

    # Step 2: existing cohorts at the same sites.
    parts = [ss]
    for c in COHORTS:
        f = out / "vcf" / f"{c}.sites.vcf.gz"
        run([BCF, "view", "-R", regions, "-Oz", "-o", f, Path(a.a1_vcf_dir) / f"{c}.vcf.gz"])
        run([BCF, "index", "-t", f])
        parts.append(f)
    allv = out / "vcf" / "all.vcf.gz"
    run([BCF, "merge", "-m", "snps", "--threads", a.threads, "-Oz", "-o", allv] + parts)
    keep = out / "vcf" / "all.panel.vcf.gz"
    p1 = subprocess.Popen([BCF, "view", "-m2", "-M2", "-v", "snps", "-Ou", str(allv)], stdout=subprocess.PIPE)
    run([BCF, "annotate", "-x", "INFO", "--set-id", "%CHROM:%POS", "-Oz", "-o", keep], stdin=p1.stdout)
    p1.wait()
    q = run([BCF, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", keep], capture_output=True, text=True).stdout
    ok = [f"{c}:{p}" for c, p, r, al in (l.split("\t") for l in q.splitlines()) if panel.get((c, int(p))) == (r, al)]
    (out / "vcf" / "panel_match_ids.txt").write_text("\n".join(ok) + "\n")
    for p in parts[1:] + [allv]:
        for f in (str(p), str(p) + ".tbi"):
            if os.path.exists(f):
                os.remove(f)

    # Step 3: KING-robust.
    pref = out / "vcf" / "king"
    run([a.plink2, "--vcf", keep, "--extract", out / "vcf" / "panel_match_ids.txt", "--allow-extra-chr",
         "--chr", ",".join([f"chr{i}" for i in range(1, 23)]), "--make-king-table", "--king-table-filter", "0.0884",
         "--threads", a.threads, "--out", pref])
    k = pd.read_csv(f"{pref}.kin0", sep=r"\s+")
    k.columns = [c.lstrip("#") for c in k.columns]
    k = k.rename(columns={"IID1": "id1", "IID2": "id2", "KINSHIP": "kinship", "NSNP": "nsnp"})
    ssset = set(libs)
    k["pair_type"] = ["ss_ss" if (x in ssset and y in ssset) else ("ss_existing" if (x in ssset or y in ssset)
                      else "existing_existing") for x, y in zip(k.id1, k.id2)]
    k = k[(k.nsnp >= MIN_SHARED) & (k.pair_type != "existing_existing")]
    k["relation"] = pd.cut(k.kinship, [0.0884, 0.177, 0.354, 1.0], labels=["second_degree", "first_degree", "same_individual"],
                           right=False)
    k[["id1", "id2", "pair_type", "nsnp", "kinship", "relation"]].to_csv(out / "kinship_pairs.tsv", sep="\t", index=False)
    # Coverage check (D-S2-5): every pair involving a SteatoSITE library, unfiltered, so pairs below the 1,000-site
    # floor are counted and any same-individual-level kinship among 500-999 shared sites is reported (not excluded).
    (out / "vcf" / "ss_ids.txt").write_text("\n".join(libs) + "\n")
    pa = out / "vcf" / "king_ss_all"
    run([a.plink2, "--vcf", keep, "--extract", out / "vcf" / "panel_match_ids.txt", "--allow-extra-chr",
         "--chr", ",".join([f"chr{i}" for i in range(1, 23)]), "--make-king-table",
         "--king-table-require", out / "vcf" / "ss_ids.txt", "--threads", a.threads, "--out", pa])
    ka = pd.read_csv(f"{pa}.kin0", sep=r"\s+")
    ka.columns = [c.lstrip("#") for c in ka.columns]
    ka = ka.rename(columns={"IID1": "id1", "IID2": "id2", "KINSHIP": "kinship", "NSNP": "nsnp"})
    ka["pair_type"] = ["ss_ss" if (x in ssset and y in ssset) else "ss_existing" for x, y in zip(ka.id1, ka.id2)]
    ka[["id1", "id2", "pair_type", "nsnp", "kinship"]].to_csv(out / "kinship_ss_all_pairs.tsv.gz", sep="\t", index=False)
    os.remove(f"{pa}.kin0")
    low = ka[ka.nsnp < MIN_SHARED]
    sens = ka[(ka.nsnp >= 500) & (ka.nsnp < MIN_SHARED) & (ka.kinship >= 0.354)]
    sens.to_csv(out / "kinship_same_individual_500_999_sites.tsv", sep="\t", index=False)
    coverage = {"ss_pairs_total": int(len(ka)), "ss_pairs_below_min_shared": int(len(low)),
                "ss_pairs_below_500": int((ka.nsnp < 500).sum()),
                "same_individual_at_500_999_sites": {t_: int((sens.pair_type == t_).sum()) for t_ in ["ss_ss", "ss_existing"]},
                "max_kinship_ss_pairs_any_nsnp": float(ka.kinship.max()) if len(ka) else None}
    summ = {"n_libraries": len(libs), "n_panel_sites": len(panel), "n_sites_matched_after_merge": len(ok),
            "min_dp_call": MIN_DP_CALL, "min_dp_hom": MIN_DP_HOM, "min_shared_sites": MIN_SHARED,
            "same_individual_pairs": {t_: int(((k.relation == "same_individual") & (k.pair_type == t_)).sum())
                                      for t_ in ["ss_ss", "ss_existing"]},
            "pair_coverage": coverage, "stage_read": False, "expression_read": False}
    (out / "identity_summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ))


if __name__ == "__main__":
    main()
