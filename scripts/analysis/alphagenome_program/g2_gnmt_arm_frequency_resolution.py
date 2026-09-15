#!/usr/bin/env python3
"""G2. Settle the GNMT arm-frequency disagreement between F1 and F2.

Two deposits report different four-haplotype counts for the same pair on the same panel:

  F1 (`f1-phase-20260914T193608Z/tables/gnmt_pair.json`)  ALT1-REF2 = 2/758, REF1-ALT2 = 2/758
  F2 (`f2-haplotype-20260914T230433Z/RESULTS.md` line 83) both arms = 0.0000, r2 = 1.0000

This recomputes the four two-locus haplotype counts for rs2296805 x rs2296804 directly from
`data/1kg_eur` by three routes that do not share a code path, and tests the three candidate
explanations the assignment names: founders-only counting, dropped missing calls, printing rounding.

Routes
  R1  plink2 --export vcf from the PHASED .pgen, haplotypes counted from the '|' genotype fields
  R2  plink2 --export haps  from the same .pgen, an independent exporter
  R3  plink 1.9 --ld on the UNPHASED .bed, which reports EM-estimated haplotype frequencies

Writes `tables/gnmt_arm_frequency_resolution.tsv` into the g-integration package.
"""

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys

PROJ = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
KG = PROJ / "data" / "1kg_eur"
PLINK2 = "/nfs/sw/easybuild/software/PLINK/2.0a5.13/bin/plink"
PLINK1 = "/nfs/sw/easybuild/software/PLINK/1.9b_6.21-x86_64/plink"

OUT = pathlib.Path(os.environ["G2_OUT_ROOT"])
TAB = OUT / "tables"
WORK = pathlib.Path(os.environ["G2_WORK"])
TAB.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)

CHROM = "chr6"
# hg19 panel ids; GRCh38 positions and alleles from f1-phase tables/gnmt_pair.json
V1 = dict(pvar_id="6:42928758:T:G", rsid="rs2296805", pos38=42961020, ref="T", alt="G")
V2 = dict(pvar_id="6:42931261:C:G", rsid="rs2296804", pos38=42963523, ref="C", alt="G")


def log(msg: str) -> None:
    print(msg, flush=True)


def run(cmd: list[str]) -> subprocess.CompletedProcess:
    log("$ " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log(r.stdout[-3000:])
        log(r.stderr[-3000:])
        raise SystemExit(f"command failed rc={r.returncode}")
    return r


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def counts_from_dosage_pairs(pairs: list[tuple[int, int]]) -> dict:
    """pairs are (alt dosage on v1, alt dosage on v2) per HAPLOTYPE, 0/1."""
    n11 = sum(1 for a, b in pairs if a == 1 and b == 1)
    n10 = sum(1 for a, b in pairs if a == 1 and b == 0)
    n01 = sum(1 for a, b in pairs if a == 0 and b == 1)
    n00 = sum(1 for a, b in pairs if a == 0 and b == 0)
    n = len(pairs)
    p1, p2 = (n11 + n10) / n, (n11 + n01) / n
    D = n11 / n - p1 * p2
    denom = p1 * (1 - p1) * p2 * (1 - p2)
    r = D / denom ** 0.5 if denom > 0 else float("nan")
    dmax = min(p1 * (1 - p2), p2 * (1 - p1)) if D >= 0 else min(p1 * p2, (1 - p1) * (1 - p2))
    return dict(n_haplotypes=n, n_ref_ref=n00, n_alt1_ref2=n10, n_ref1_alt2=n01, n_alt_alt=n11,
                alt1_freq=p1, alt2_freq=p2, D=D, r=r, r2=r * r,
                Dprime=(D / dmax if dmax > 0 else float("nan")))


# --------------------------------------------------------------------------- R1: phased VCF
def route_vcf() -> tuple[dict, dict]:
    ex = WORK / "extract.txt"
    ex.write_text(f"{V1['pvar_id']}\n{V2['pvar_id']}\n")
    pfx = WORK / "vcf_gnmt"
    run([PLINK2, "--pfile", str(KG / f"{CHROM}_eur"), "--extract", str(ex),
         "--export", "vcf", "--out", str(pfx)])
    vcf = pathlib.Path(str(pfx) + ".vcf")
    rows = {}
    diag = dict(n_samples=None, n_phased_calls=0, n_unphased_calls=0, n_missing_calls=0)
    for line in vcf.read_text().splitlines():
        if line.startswith("##"):
            continue
        f = line.split("\t")
        if line.startswith("#CHROM"):
            diag["n_samples"] = len(f) - 9
            continue
        gts = f[9:]
        hap_a, hap_b = [], []
        for g in gts:
            g0 = g.split(":")[0]
            if "|" in g0:
                x, y = g0.split("|")
                diag["n_phased_calls"] += 1
            elif "/" in g0:
                x, y = g0.split("/")
                diag["n_unphased_calls"] += 1
            else:
                x = y = "."
            if x == "." or y == ".":
                diag["n_missing_calls"] += 1
                hap_a.append(-1)
                hap_b.append(-1)
            else:
                hap_a.append(int(x))
                hap_b.append(int(y))
        rows[f[2]] = dict(ref=f[3], alt=f[4], a=hap_a, b=hap_b)
    a1, a2 = rows[V1["pvar_id"]], rows[V2["pvar_id"]]
    assert (a1["ref"], a1["alt"]) == (V1["ref"], V1["alt"]), (a1["ref"], a1["alt"])
    assert (a2["ref"], a2["alt"]) == (V2["ref"], V2["alt"]), (a2["ref"], a2["alt"])
    all_pairs = list(zip(a1["a"] + a1["b"], a2["a"] + a2["b"]))
    kept = [(x, y) for x, y in all_pairs if x >= 0 and y >= 0]
    diag["n_haplotype_slots"] = len(all_pairs)
    diag["n_haplotypes_dropped_for_missing"] = len(all_pairs) - len(kept)
    # per-individual phase pattern, for the EM explanation
    per_ind = list(zip(a1["a"], a1["b"], a2["a"], a2["b"]))
    dbl_het = [i for i, (pa, pb, qa, qb) in enumerate(per_ind)
               if pa + pb == 1 and qa + qb == 1]
    diag["n_double_het_individuals"] = len(dbl_het)
    diag["n_double_het_in_coupling"] = sum(
        1 for i in dbl_het if per_ind[i][0] == per_ind[i][2])
    diag["n_double_het_in_repulsion"] = sum(
        1 for i in dbl_het if per_ind[i][0] != per_ind[i][2])
    return counts_from_dosage_pairs(kept), diag


# --------------------------------------------------------------------------- R2: .haps export
def route_haps() -> dict:
    ex = WORK / "extract.txt"
    pfx = WORK / "haps_gnmt"
    run([PLINK2, "--pfile", str(KG / f"{CHROM}_eur"), "--extract", str(ex),
         "--export", "haps", "--out", str(pfx)])
    haps = pathlib.Path(str(pfx) + ".haps")
    got = {}
    for line in haps.read_text().splitlines():
        f = line.split()
        # chrom id pos a1 a2 then 2 columns per sample
        got[f[1]] = dict(a1=f[3], a2=f[4], h=[int(x) for x in f[5:]])
    h1, h2 = got[V1["pvar_id"]], got[V2["pvar_id"]]
    # .haps codes 0 = the file's A1 allele, 1 = A2. Map to ALT dosage on GRCh38.
    d1 = [(x if h1["a2"] == V1["alt"] else 1 - x) for x in h1["h"]]
    d2 = [(x if h2["a2"] == V2["alt"] else 1 - x) for x in h2["h"]]
    out = counts_from_dosage_pairs(list(zip(d1, d2)))
    out["haps_allele_columns"] = f"{V1['pvar_id']} A1={h1['a1']} A2={h1['a2']}; " \
                                 f"{V2['pvar_id']} A1={h2['a1']} A2={h2['a2']}"
    return out


# --------------------------------------------------------------------------- R3: plink 1.9 --ld
def route_plink1_ld(founders_only: bool) -> dict:
    pfx = WORK / ("ld_founders" if founders_only else "ld_all")
    cmd = [PLINK1, "--bfile", str(KG / f"{CHROM}_eur"), "--ld", V1["pvar_id"], V2["pvar_id"],
           "--out", str(pfx)]
    if not founders_only:
        cmd.append("--nonfounders")
    run(cmd)
    text = pathlib.Path(str(pfx) + ".log").read_text()
    st = re.search(r"R-sq = (\S+)\s+D' = (\S+)", text)
    haps = {}
    for line in text.splitlines():
        mm = re.match(r"\s{5,}([ACGT]{2})\s+(-?[\d.eE+-]+)\s+[\d.eE+-]+\s*$", line)
        if mm:
            haps[mm.group(1)] = float(mm.group(2))
    founders = re.search(r"(\d+) founders and (\d+) nonfounders present", text)
    return dict(r2=float(st.group(1)), Dprime=float(st.group(2)),
                f_ref_ref=haps.get(V1["ref"] + V2["ref"]),
                f_alt1_ref2=haps.get(V1["alt"] + V2["ref"]),
                f_ref1_alt2=haps.get(V1["ref"] + V2["alt"]),
                f_alt_alt=haps.get(V1["alt"] + V2["alt"]),
                printed_lines=[l.strip() for l in text.splitlines()
                               if re.match(r"\s{5,}[ACGT]{2}\s", l)],
                founders=int(founders.group(1)) if founders else None,
                nonfounders=int(founders.group(2)) if founders else None,
                raw_log=str(pfx) + ".log")


# --------------------------------------------------------------------------- missingness
def route_missing() -> dict:
    pfx = WORK / "miss"
    run([PLINK1, "--bfile", str(KG / f"{CHROM}_eur"), "--extract", str(WORK / "extract.txt"),
         "--missing", "--out", str(pfx)])
    lmiss = pathlib.Path(str(pfx) + ".lmiss").read_text().splitlines()
    return dict(lmiss=[" ".join(l.split()) for l in lmiss])


def main() -> None:
    log("== R1 phased VCF from the .pgen")
    vcf_counts, vcf_diag = route_vcf()
    log(json.dumps(vcf_counts, indent=1))
    log(json.dumps(vcf_diag, indent=1))

    log("== R2 .haps export from the same .pgen")
    haps_counts = route_haps()
    log(json.dumps(haps_counts, indent=1))

    log("== R3 plink 1.9 --ld on the .bed, founders only (default)")
    ld_f = route_plink1_ld(founders_only=True)
    log(json.dumps({k: v for k, v in ld_f.items() if k != "printed_lines"}, indent=1))
    log("   printed haplotype block: " + " | ".join(ld_f["printed_lines"]))

    log("== R3b plink 1.9 --ld with --nonfounders")
    ld_a = route_plink1_ld(founders_only=False)
    log(json.dumps({k: v for k, v in ld_a.items() if k != "printed_lines"}, indent=1))

    log("== missingness at the two variants")
    miss = route_missing()
    log("\n".join(miss["lmiss"]))

    fam = (KG / f"{CHROM}_eur.fam").read_text().splitlines()
    n_nonfounder = sum(1 for l in fam if l.split()[2] != "0" or l.split()[3] != "0")

    n = vcf_counts["n_haplotypes"]
    rows = []

    def add(route, arm, count, denom, freq, note):
        rows.append(dict(route=route, arm=arm, haplotype_count=count, denominator=denom,
                         frequency=freq, note=note))

    lbl = {"ref_ref": f"REF1-REF2 ({V1['ref']}{V2['ref']})",
           "alt1_ref2": f"ALT1-REF2 ({V1['alt']}{V2['ref']})",
           "ref1_alt2": f"REF1-ALT2 ({V1['ref']}{V2['alt']})",
           "alt_alt": f"ALT1-ALT2 ({V1['alt']}{V2['alt']})"}
    for key in ("ref_ref", "alt1_ref2", "ref1_alt2", "alt_alt"):
        c = vcf_counts["n_" + key]
        add("R1_plink2_export_vcf_phased_hardcalls", lbl[key], c, n, c / n,
            "direct count of phased hardcalls; this is the number the resolution adopts")
    for key in ("ref_ref", "alt1_ref2", "ref1_alt2", "alt_alt"):
        c = haps_counts["n_" + key]
        add("R2_plink2_export_haps_phased_hardcalls", lbl[key], c, haps_counts["n_haplotypes"],
            c / haps_counts["n_haplotypes"], "independent exporter, same .pgen")
    for key, fk in (("ref_ref", "f_ref_ref"), ("alt1_ref2", "f_alt1_ref2"),
                    ("ref1_alt2", "f_ref1_alt2"), ("alt_alt", "f_alt_alt")):
        f = ld_f[fk]
        add("R3_plink1.9_ld_EM_on_unphased_bed", lbl[key],
            "" if f is None else round(f * n, 4), n, f,
            "EM estimate from genotypes; the .bed format carries no phase")

    with open(TAB / "gnmt_arm_frequency_resolution.tsv", "w") as h:
        h.write("pair\tbuild\troute\tarm\thaplotype_count\tdenominator\tfrequency\tnote\n")
        pair = f"{V1['rsid']} x {V2['rsid']}"
        build = f"GRCh38 chr6:{V1['pos38']} / chr6:{V2['pos38']}"
        for r in rows:
            freq = "" if r["frequency"] is None else format(r["frequency"], ".6f")
            h.write("\t".join([pair, build, r["route"], r["arm"], str(r["haplotype_count"]),
                               str(r["denominator"]), freq, r["note"]]) + "\n")

    resolution = dict(
        pair=f"{V1['rsid']} x {V2['rsid']}",
        panel=str(KG / f"{CHROM}_eur"),
        n_samples=vcf_diag["n_samples"], n_haplotypes=n,
        fam_nonfounders=n_nonfounder,
        plink1_founders=ld_f["founders"], plink1_nonfounders=ld_f["nonfounders"],
        missing_calls=vcf_diag["n_missing_calls"],
        haplotypes_dropped_for_missing=vcf_diag["n_haplotypes_dropped_for_missing"],
        phased_calls=vcf_diag["n_phased_calls"], unphased_calls=vcf_diag["n_unphased_calls"],
        R1=vcf_counts, R2=haps_counts, R3_founders=ld_f, R3_nonfounders=ld_a,
        double_het_individuals=vcf_diag["n_double_het_individuals"],
        double_het_in_coupling=vcf_diag["n_double_het_in_coupling"],
        double_het_in_repulsion=vcf_diag["n_double_het_in_repulsion"],
        em_arithmetic=dict(
            phased_ref_ref=vcf_counts["n_ref_ref"], phased_alt_alt=vcf_counts["n_alt_alt"],
            em_ref_ref=None if ld_f["f_ref_ref"] is None else round(ld_f["f_ref_ref"] * n, 6),
            em_alt_alt=None if ld_f["f_alt_alt"] is None else round(ld_f["f_alt_alt"] * n, 6)),
        inputs=[dict(path=str(KG / f"{CHROM}_eur{ext}"), sha256=sha256(KG / f"{CHROM}_eur{ext}"))
                for ext in (".psam", ".fam", ".bim")],
        plink2_version=subprocess.run([PLINK2, "--version"], capture_output=True,
                                      text=True).stdout.strip(),
        plink1_version=subprocess.run([PLINK1, "--version"], capture_output=True,
                                      text=True).stdout.strip(),
    )
    (TAB / "gnmt_arm_frequency_resolution.json").write_text(json.dumps(resolution, indent=1))
    log("== resolution written")
    log(json.dumps({k: v for k, v in resolution.items()
                    if k not in ("R1", "R2", "R3_founders", "R3_nonfounders", "inputs")}, indent=1))
    log("G2_DONE")


if __name__ == "__main__":
    sys.exit(main())
