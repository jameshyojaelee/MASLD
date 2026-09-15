#!/usr/bin/env python
"""F1: what phase substrate exists for the P5 local-haplotype pairs.

Executes scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md section 7, F1 only.
Decision rules are fixed in the package PRESPECIFICATION.md, written before this ran.

Stages
  A  probe every data/1kg_eur/chr*_eur.pgen for retained phase (plink2 --pgen-info)
  B  lift the 1,587 P5 pairs (1,380 distinct GRCh38 SNVs) plus the GNMT pair to GRCh37,
     validate the lift two independent ways
  C  export phased haplotypes for the lifted variants and count two-locus haplotypes
  D  read TOP-LD EUR GRCh38 R2, Dprime and correlation sign for every covered pair
  E  merge, decide arm plausibility, emit the F1.1 verdict
  F  count donors heterozygous at both sites of a P5 pair in the two allelic deposits
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import hashlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT = Path(sys.argv[1]).resolve()
TAB = OUT / "tables"
PROV = OUT / "provenance"
TMP = Path(os.environ.get("F1_TMP", str(OUT / "_tmp")))
for d in (TAB, PROV, TMP):
    d.mkdir(parents=True, exist_ok=True)

PLINK2 = "/nfs/sw/easybuild/software/PLINK/2.0a5.13/bin/plink"
CHAIN = Path("/gpfs/commons/groups/sanjana_lab/MYCScreen/TCGA/hg38ToHg19.over.chain")
KG = PROJ / "data/1kg_eur"
TOPLD = PROJ / "GWAS/finemapping/data/ld_ref/topld_raw/EUR"
TOPLD_EUR_BLOCKS = PROJ / "GWAS/finemapping/data/ld_ref/topld_eur"
P5 = PROJ / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables"
ASE = PROJ / "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables"

SEED = 20260914
ARM_UNOBSERVED = 0.01
NPROC = int(os.environ.get("F1_NPROC", "8"))

# GNMT pair named in the spec (GRCh38). REF/ALT verified against the TOP-LD GRCh38 info file
# (EUR_chr6_..._info.csv.gz: 42961020,rs2296805,0.4228,T,G and 42963523,rs2296804,0.4232,C,G).
GNMT = [("chr6", 42961020, "rs2296805", "T", "G"), ("chr6", 42963523, "rs2296804", "C", "G")]
GNMT_SIGNAL = "GNMT_spec_pair_not_in_p5_universe"

CHROMS = [f"chr{i}" for i in range(1, 23)]
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def log(msg: str) -> None:
    print(f"[f1] {msg}", flush=True)


def sha256(path: Path, cap: int = 200 * 1024 * 1024) -> str:
    if path.stat().st_size > cap:
        return "not_hashed_over_200MB"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rc(s: str) -> str:
    return "".join(COMP.get(b, "N") for b in reversed(s))


# --------------------------------------------------------------------------- A
def stage_a() -> pd.DataFrame:
    """Does the local 1000 Genomes EUR pgen retain phase?"""
    rows = []
    for c in CHROMS:
        pfx = KG / f"{c}_eur"
        cmd = [PLINK2, "--pfile", str(pfx), "--pgen-info", "--out", str(TMP / f"pgeninfo_{c}")]
        res = subprocess.run(cmd, capture_output=True, text=True)
        txt = res.stdout + res.stderr
        rows.append(
            dict(
                chrom=c,
                command=" ".join(cmd),
                returncode=res.returncode,
                explicitly_phased_hardcalls=("Explicitly phased hardcalls present" in txt),
                no_phase_line=("Explicitly phased hardcalls absent" in txt
                               or "hardcalls are unphased" in txt),
                n_variants=next((int(l.split()[1].replace(",", ""))
                                 for l in txt.splitlines() if l.strip().startswith("Variants:")), -1),
                n_samples=next((int(l.split()[1].replace(",", ""))
                                for l in txt.splitlines() if l.strip().startswith("Samples:")), -1),
                pgen_info_block="; ".join(
                    l.strip() for l in txt.splitlines()
                    if l.startswith("  ") and ("Variants:" in l or "Samples:" in l
                                               or "phased" in l or "dosages" in l
                                               or "REF alleles" in l or "allele count" in l)
                ),
            )
        )
        log(f"A {c}: phased={rows[-1]['explicitly_phased_hardcalls']} n={rows[-1]['n_variants']}")
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "pgen_phase_probe.tsv", sep="\t", index=False)
    return df


# --------------------------------------------------------------------------- B
def load_chain(chain: Path):
    """Parse a UCSC chain file into per-source-chrom block lists.

    Returns {tName: [(score, tStart, tEnd, blocks, qName, qSize, qStrand)]} where blocks is a
    list of (t_block_start, q_block_start, size) in chain-internal coordinates.
    """
    chains: dict[str, list] = {}
    cur = None
    with open(chain) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                cur = None
                continue
            if line.startswith("chain"):
                f = line.split()
                score = int(f[1])
                tName, tSize, tStrand, tStart, tEnd = f[2], int(f[3]), f[4], int(f[5]), int(f[6])
                qName, qSize, qStrand, qStart, qEnd = f[7], int(f[8]), f[9], int(f[10]), int(f[11])
                assert tStrand == "+", "source strand is always + in UCSC chains"
                cur = dict(score=score, tName=tName, tStart=tStart, tEnd=tEnd, qName=qName,
                           qSize=qSize, qStrand=qStrand, blocks=[], _tp=tStart, _qp=qStart)
                chains.setdefault(tName, []).append(cur)
                continue
            if cur is None:
                continue
            f = line.split()
            size = int(f[0])
            cur["blocks"].append((cur["_tp"], cur["_qp"], size))
            if len(f) == 3:
                cur["_tp"] += size + int(f[1])
                cur["_qp"] += size + int(f[2])
    for v in chains.values():
        v.sort(key=lambda d: -d["score"])
    return chains


def lift_point(chains, chrom: str, pos1: int):
    """Lift a 1-based GRCh38 SNV position. Returns (qChrom, qPos1, strand) or None.

    Highest-scoring chain that contains the base wins, matching UCSC liftOver behaviour.
    """
    p0 = pos1 - 1
    for ch in chains.get(chrom, []):
        if not (ch["tStart"] <= p0 < ch["tEnd"]):
            continue
        for tb, qb, size in ch["blocks"]:
            if tb <= p0 < tb + size:
                off = p0 - tb
                q0 = qb + off
                if ch["qStrand"] == "+":
                    return ch["qName"], q0 + 1, "+"
                return ch["qName"], ch["qSize"] - q0, "-"
    return None


def read_pvar(chrom: str) -> pd.DataFrame:
    path = KG / f"{chrom}_eur.pvar"
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            f = line.rstrip("\n").split("\t")
            info = f[7]
            rsq = np.nan
            eur_af = np.nan
            for kv in info.split(";"):
                if kv.startswith("RSQ="):
                    rsq = float(kv[4:])
                elif kv.startswith("EUR_AF="):
                    eur_af = float(kv[7:])
            rows.append((int(f[1]), f[2], f[3], f[4], rsq, eur_af))
    return pd.DataFrame(rows, columns=["pos37", "pvar_id", "pvar_ref", "pvar_alt",
                                       "pvar_rsq", "pvar_eur_af"])


def stage_b(pairs: pd.DataFrame) -> pd.DataFrame:
    log("B: parsing chain")
    chains = load_chain(CHAIN)

    vs = pd.unique(pd.concat([pairs.v1, pairs.v2]))
    gnmt_ids = {f"{c}:{p}:{r}:{a}": rs for c, p, rs, r, a in GNMT}
    recs = []
    for v in vs:
        c, p, r, a = v.split(":")
        recs.append(dict(variant=v, chrom=c, pos38=int(p), ref38=r, alt38=a,
                         source=gnmt_ids.get(v, "p5")))
    lv = pd.DataFrame(recs)

    got = [lift_point(chains, r.chrom, r.pos38) for r in lv.itertuples()]
    lv["chrom37"] = [g[0] if g else None for g in got]
    lv["pos37"] = [g[1] if g else -1 for g in got]
    lv["lift_strand"] = [g[2] if g else None for g in got]
    lv["lifted"] = lv.pos37 > 0
    log(f"B: lifted {int(lv.lifted.sum())} / {len(lv)}")

    # a lift that lands on a different chromosome cannot be looked up in that chromosome's panel
    lv.loc[lv.lifted & (lv.chrom37 != lv.chrom), "pos37"] = -1
    n_cross = int((lv.lifted & (lv.chrom37 != lv.chrom)).sum())
    if n_cross:
        log(f"B: {n_cross} variants lifted to a different chromosome; dropped from the panel branch")

    out = []
    for chrom, sub in lv.groupby("chrom"):
        pv = read_pvar(chrom)
        m = sub.merge(pv, on="pos37", how="left", suffixes=("", "_pv"))
        out.append(m)
    lv = pd.concat(out, ignore_index=True)

    # allele validation, allowing a REF/ALT swap and reverse complement on minus-strand lifts
    def allele_state(r):
        if not isinstance(r.pvar_ref, str):
            return "no_panel_variant_at_lifted_position"
        a38 = {r.ref38, r.alt38}
        if r.lift_strand == "-":
            a38 = {COMP.get(x, "N") for x in a38}
        a37 = {r.pvar_ref, r.pvar_alt}
        if a38 != a37:
            return "allele_mismatch"
        want_ref = COMP.get(r.ref38, "N") if r.lift_strand == "-" else r.ref38
        return "match_same_orientation" if r.pvar_ref == want_ref else "match_ref_alt_swapped"

    lv["allele_state"] = [allele_state(r) for r in lv.itertuples()]
    # dedupe: one panel row per variant, prefer an allele match
    rank = {"match_same_orientation": 0, "match_ref_alt_swapped": 1,
            "allele_mismatch": 2, "no_panel_variant_at_lifted_position": 3}
    lv["_r"] = lv.allele_state.map(rank)
    lv = lv.sort_values(["variant", "_r", "pvar_id"]).drop_duplicates("variant").drop(columns="_r")
    lv["panel_usable"] = lv.allele_state.isin(["match_same_orientation", "match_ref_alt_swapped"])
    log("B allele_state:\n" + lv.allele_state.value_counts().to_string())

    # independent check: rsID crosswalk GRCh38 TOP-LD info -> GRCh37 topld_eur bim
    lv = crosswalk_check(lv)
    lv.to_csv(TAB / "variant_liftover.tsv", sep="\t", index=False)
    return lv


def crosswalk_check(lv: pd.DataFrame) -> pd.DataFrame:
    """Second, chain-independent GRCh38->GRCh37 mapping via rsID."""
    rs38 = {}
    for chrom in sorted(lv.chrom.unique()):
        want = set(lv.loc[lv.chrom == chrom, "pos38"].tolist())
        f = TOPLD / f"EUR_{chrom}_Dprime_extended_0.2_1000000_info.csv.gz"
        if not f.exists():
            continue
        with gzip.open(f, "rt") as fh:
            next(fh)
            for line in fh:
                p, rsid, maf, ref, alt = line.rstrip("\n").split(",")
                if int(p) in want:
                    rs38[(chrom, int(p))] = (rsid, float(maf), ref, alt)
    lv["topld_rsid"] = [rs38.get((r.chrom, r.pos38), (None,) * 4)[0] for r in lv.itertuples()]
    lv["topld_maf"] = [rs38.get((r.chrom, r.pos38), (None, np.nan, None, None))[1]
                       for r in lv.itertuples()]
    lv["topld_ref"] = [rs38.get((r.chrom, r.pos38), (None,) * 4)[2] for r in lv.itertuples()]
    lv["topld_alt"] = [rs38.get((r.chrom, r.pos38), (None,) * 4)[3] for r in lv.itertuples()]

    want_rs = set(x for x in lv.topld_rsid.dropna().tolist() if str(x).startswith("rs"))
    rs37 = {}
    for chrom in sorted(lv.chrom.unique()):
        d = TOPLD_EUR_BLOCKS / chrom
        if not d.is_dir():
            continue
        for blk in sorted(d.iterdir()):
            bim = blk / f"{blk.name}.bim"
            if not bim.exists():
                continue
            with open(bim) as fh:
                for line in fh:
                    f = line.split()
                    if f[1] in want_rs:
                        rs37.setdefault(f[1], (f"chr{f[0]}", int(f[3])))
    lv["crosswalk_pos37"] = [rs37.get(r.topld_rsid, (None, -1))[1] for r in lv.itertuples()]
    ok = lv.crosswalk_pos37 > 0
    lv["crosswalk_state"] = np.where(
        ~ok, "no_rsid_crosswalk",
        np.where(lv.crosswalk_pos37 == lv.pos37, "agrees_with_chain", "disagrees_with_chain"))
    log("B crosswalk:\n" + lv.crosswalk_state.value_counts().to_string())
    return lv


# --------------------------------------------------------------------------- C
def export_haps(chrom: str, ids: list[str]) -> pd.DataFrame | None:
    """Export phased genotypes as VCF for the requested panel variant ids."""
    if not ids:
        return None
    ex = TMP / f"extract_{chrom}.txt"
    ex.write_text("\n".join(ids) + "\n")
    pfx = TMP / f"vcf_{chrom}"
    cmd = [PLINK2, "--pfile", str(KG / f"{chrom}_eur"), "--extract", str(ex),
           "--export", "vcf", "--out", str(pfx)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    vcf = Path(str(pfx) + ".vcf")
    if res.returncode != 0 or not vcf.exists():
        log(f"C {chrom}: plink2 export failed rc={res.returncode}\n{res.stderr[-1500:]}")
        return None
    rows = []
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                continue
            gts = f[9:]
            a1 = np.empty(len(gts), dtype=np.int8)
            a2 = np.empty(len(gts), dtype=np.int8)
            phased = 0
            for i, g in enumerate(gts):
                g = g.split(":")[0]
                if "|" in g:
                    x, y = g.split("|")
                    phased += 1
                elif "/" in g:
                    x, y = g.split("/")
                else:
                    x = y = "."
                a1[i] = -1 if x == "." else int(x)
                a2[i] = -1 if y == "." else int(y)
            rows.append(dict(pvar_id=f[2], pos37=int(f[1]), vcf_ref=f[3], vcf_alt=f[4],
                             n_phased_calls=phased, n_calls=len(gts),
                             h=np.concatenate([a1, a2])))
    return pd.DataFrame(rows)


def stage_c(pairs: pd.DataFrame, lv: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    usable = lv[lv.panel_usable].set_index("variant")
    hap = {}
    meta = {}
    for chrom in sorted(usable.chrom.unique()):
        sub = usable[usable.chrom == chrom]
        df = export_haps(chrom, sorted(sub.pvar_id.unique().tolist()))
        if df is None:
            continue
        for r in df.itertuples():
            hap[(chrom, r.pvar_id)] = r.h
            meta[(chrom, r.pvar_id)] = (r.vcf_ref, r.vcf_alt, r.n_phased_calls, r.n_calls)
        log(f"C {chrom}: {len(df)} variants exported, "
            f"phased calls {df.n_phased_calls.sum()}/{df.n_calls.sum()}")

    qc = pd.DataFrame(
        [dict(chrom=k[0], pvar_id=k[1], vcf_ref=v[0], vcf_alt=v[1],
              n_phased_calls=v[2], n_calls=v[3]) for k, v in meta.items()])
    qc.to_csv(TAB / "panel_variant_export_qc.tsv", sep="\t", index=False)

    rows = []
    for pr in pairs.itertuples():
        rec = dict(signal_uid=pr.signal_uid, v1=pr.v1, v2=pr.v2)
        if pr.v1 not in usable.index or pr.v2 not in usable.index:
            rec["panel_state"] = "variant_not_in_panel"
            rows.append(rec)
            continue
        m1, m2 = usable.loc[pr.v1], usable.loc[pr.v2]
        k1, k2 = (m1.chrom, m1.pvar_id), (m2.chrom, m2.pvar_id)
        if k1 not in hap or k2 not in hap:
            rec["panel_state"] = "variant_not_exported"
            rows.append(rec)
            continue
        h1, h2 = hap[k1], hap[k2]
        ok = (h1 >= 0) & (h2 >= 0)
        h1, h2 = h1[ok], h2[ok]
        # panel coding 0/1 == pvar REF/ALT; map to GRCh38 ALT dosage
        d1 = h1 if m1.allele_state == "match_same_orientation" else 1 - h1
        d2 = h2 if m2.allele_state == "match_same_orientation" else 1 - h2
        n = len(d1)
        n11 = int(np.sum((d1 == 1) & (d2 == 1)))
        n10 = int(np.sum((d1 == 1) & (d2 == 0)))
        n01 = int(np.sum((d1 == 0) & (d2 == 1)))
        n00 = int(np.sum((d1 == 0) & (d2 == 0)))
        p1, p2 = (n11 + n10) / n, (n11 + n01) / n
        h11 = n11 / n
        D = h11 - p1 * p2
        if D >= 0:
            dmax = min(p1 * (1 - p2), p2 * (1 - p1))
        else:
            dmax = min(p1 * p2, (1 - p1) * (1 - p2))
        denom = p1 * (1 - p1) * p2 * (1 - p2)
        rec.update(
            panel_state="ok", panel_n_haplotypes=n,
            panel_alt1_freq=p1, panel_alt2_freq=p2,
            panel_n_ref_ref=n00, panel_n_alt1_ref2=n10,
            panel_n_ref1_alt2=n01, panel_n_alt_alt=n11,
            panel_f_ref_ref=n00 / n, panel_f_alt1_ref2=n10 / n,
            panel_f_ref1_alt2=n01 / n, panel_f_alt_alt=n11 / n,
            panel_D=D,
            panel_Dprime=(abs(D) / dmax) if dmax > 0 else np.nan,
            panel_r=(D / np.sqrt(denom)) if denom > 0 else np.nan,
            panel_r2=(D * D / denom) if denom > 0 else np.nan,
            panel_rsq_v1=m1.pvar_rsq, panel_rsq_v2=m2.pvar_rsq,
        )
        rows.append(rec)
    df = pd.DataFrame(rows)
    log("C panel_state:\n" + df.panel_state.value_counts().to_string())
    return df, qc


# --------------------------------------------------------------------------- D
def topld_one_chrom(args):
    chrom, positions, pair_keys = args
    info = {}
    f_info = TOPLD / f"EUR_{chrom}_Dprime_extended_0.2_1000000_info.csv.gz"
    if f_info.exists():
        with gzip.open(f_info, "rt") as fh:
            next(fh)
            for line in fh:
                p, rsid, maf, ref, alt = line.rstrip("\n").split(",")
                ip = int(p)
                if ip in positions:
                    info[ip] = (rsid, float(maf), ref, alt)
    f_ld = TOPLD / f"EUR_{chrom}_Dprime_extended_0.2_1000000.csv.gz"
    hits = {}
    nrows = 0
    if f_ld.exists():
        posarr = np.fromiter(positions, dtype=np.int64)
        posarr.sort()
        reader = pd.read_csv(f_ld, chunksize=4_000_000,
                             dtype={"SNP1": np.int64, "SNP2": np.int64,
                                    "R2": np.float64, "Dprime": np.float64,
                                    "+/-corr": "category"})
        for ch in reader:
            ch = ch.rename(columns={"+/-corr": "corr_sign"})
            nrows += len(ch)
            i1 = np.searchsorted(posarr, ch["SNP1"].to_numpy())
            i1 = np.clip(i1, 0, len(posarr) - 1)
            keep = posarr[i1] == ch["SNP1"].to_numpy()
            if not keep.any():
                continue
            ch = ch[keep]
            i2 = np.searchsorted(posarr, ch["SNP2"].to_numpy())
            i2 = np.clip(i2, 0, len(posarr) - 1)
            keep = posarr[i2] == ch["SNP2"].to_numpy()
            ch = ch[keep]
            for r in ch.itertuples(index=False):
                a, b = (r.SNP1, r.SNP2) if r.SNP1 <= r.SNP2 else (r.SNP2, r.SNP1)
                if (a, b) in pair_keys:
                    hits[(a, b)] = (float(r.R2), float(r.Dprime), str(r.corr_sign))
    log(f"D {chrom}: scanned {nrows:,} LD rows, {len(hits)} requested pairs found, "
        f"{len(info)} of {len(positions)} positions in panel")
    return chrom, info, hits, nrows


def stage_d(pairs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    jobs = []
    for chrom, sub in pairs.groupby("chrom"):
        pos = set(sub.pos1.tolist()) | set(sub.pos2.tolist())
        keys = set()
        for r in sub.itertuples():
            keys.add((min(r.pos1, r.pos2), max(r.pos1, r.pos2)))
        jobs.append((chrom, pos, keys))
    jobs.sort(key=lambda j: -len(j[1]))
    with Pool(min(NPROC, len(jobs))) as pool:
        res = pool.map(topld_one_chrom, jobs)

    info_rows, hit_rows, scan_rows = [], [], []
    for chrom, info, hits, nrows in res:
        scan_rows.append(dict(chrom=chrom, topld_rows_scanned=nrows))
        for p, (rsid, maf, ref, alt) in info.items():
            info_rows.append(dict(chrom=chrom, pos38=p, topld_rsid=rsid, topld_maf=maf,
                                  topld_ref=ref, topld_alt=alt))
        for (a, b), (r2, dp, sign) in hits.items():
            hit_rows.append(dict(chrom=chrom, posA=a, posB=b, topld_r2=r2,
                                 topld_dprime=dp, topld_corr_sign=sign))
    pd.DataFrame(scan_rows).to_csv(TAB / "topld_scan_log.tsv", sep="\t", index=False)
    inf = pd.DataFrame(info_rows)
    inf.to_csv(TAB / "topld_variant_info.tsv", sep="\t", index=False)
    return inf, pd.DataFrame(hit_rows)


# --------------------------------------------------------------------------- F
def stage_f(pairs: pd.DataFrame) -> pd.DataFrame:
    out = []
    for tag, fn in [("GSE281367", "allelic_donor_counts.tsv.gz"),
                    ("GSE244832", "gse244832_allelic_donor_counts.tsv.gz")]:
        d = pd.read_csv(ASE / fn, sep="\t")
        d["key"] = d.chrom + ":" + d.pos.astype(str)
        het = d[d.is_het].groupby("key").donor.apply(set).to_dict()
        measured = set(d.key)
        for pr in pairs.itertuples():
            k1 = f"{pr.chrom}:{pr.pos1}"
            k2 = f"{pr.chrom}:{pr.pos2}"
            m1, m2 = k1 in measured, k2 in measured
            s1, s2 = het.get(k1, set()), het.get(k2, set())
            out.append(dict(dataset=tag, signal_uid=pr.signal_uid, gene=pr.gene,
                            v1=pr.v1, v2=pr.v2, separation_bp=pr.separation_bp,
                            site1_measured=m1, site2_measured=m2,
                            both_sites_measured=m1 and m2,
                            n_het_site1=len(s1), n_het_site2=len(s2),
                            n_donors_het_at_both=len(s1 & s2),
                            donors_het_at_both=",".join(sorted(s1 & s2))))
    df = pd.DataFrame(out)
    df.to_csv(TAB / "ase_pair_het_donors.tsv", sep="\t", index=False)
    return df


# --------------------------------------------------------------------------- main
def main() -> None:
    np.random.seed(SEED)
    pairs = pd.read_csv(P5 / "haplotype_pairs.tsv", sep="\t")
    pairs["is_gnmt"] = False
    # the spec's named GNMT pair is carried through every stage but never enters the F1.1
    # denominator, because it is absent from the 989-signal universe and from P5.
    g = pd.DataFrame([dict(
        signal_uid=GNMT_SIGNAL, universe="not_in_989_signal_universe", gene="GNMT",
        ensembl="ENSG00000124713", analysis_block="chr6:not_in_p5",
        v1=f"{GNMT[0][0]}:{GNMT[0][1]}:{GNMT[0][3]}:{GNMT[0][4]}",
        v2=f"{GNMT[1][0]}:{GNMT[1][1]}:{GNMT[1][3]}:{GNMT[1][4]}",
        w1=np.nan, w2=np.nan, weight_product=np.nan,
        separation_bp=abs(GNMT[1][1] - GNMT[0][1]), chrom="chr6", is_gnmt=True)])
    pairs = pd.concat([pairs, g], ignore_index=True)
    pairs["pos1"] = pairs.v1.str.split(":").str[1].astype(int)
    pairs["pos2"] = pairs.v2.str.split(":").str[1].astype(int)
    pairs["ref1"] = pairs.v1.str.split(":").str[2]
    pairs["alt1"] = pairs.v1.str.split(":").str[3]
    pairs["ref2"] = pairs.v2.str.split(":").str[2]
    pairs["alt2"] = pairs.v2.str.split(":").str[3]
    log(f"pairs {len(pairs)} distinct variants "
        f"{len(pd.unique(pd.concat([pairs.v1, pairs.v2])))}")

    eff = pd.read_csv(P5 / "haplotype_effects.tsv", sep="\t",
                      usecols=["signal_uid", "v1", "v2"])
    scored = set(zip(eff.signal_uid, eff.v1, eff.v2))
    pairs["p5_scored"] = [(r.signal_uid, r.v1, r.v2) in scored for r in pairs.itertuples()]

    phase = stage_a()
    lv = stage_b(pairs)
    panel, _ = stage_c(pairs, lv)
    inf, hits = stage_d(pairs)

    # ---- merge TOP-LD onto pairs
    pairs["posA"] = np.minimum(pairs.pos1, pairs.pos2)
    pairs["posB"] = np.maximum(pairs.pos1, pairs.pos2)
    m = pairs.merge(hits, on=["chrom", "posA", "posB"], how="left")
    ipos = inf.set_index(["chrom", "pos38"]) if len(inf) else pd.DataFrame()
    def in_topld(c, p):
        return (c, p) in ipos.index if len(ipos) else False
    m["v1_in_topld"] = [in_topld(r.chrom, r.pos1) for r in m.itertuples()]
    m["v2_in_topld"] = [in_topld(r.chrom, r.pos2) for r in m.itertuples()]
    m["topld_pair_state"] = np.where(
        m.topld_dprime.notna(), "covered",
        np.where(m.v1_in_topld & m.v2_in_topld,
                 "both_variants_in_panel_pair_below_r2_0.2", "variant_absent_from_topld"))

    m = m.merge(panel, on=["signal_uid", "v1", "v2"], how="left")

    # ---- TOP-LD sign convention check against the phased panel
    both = m[(m.topld_corr_sign.notna()) & (m.panel_r.notna()) & (m.panel_r.abs() > 1e-9)].copy()
    topld_alleles = inf.set_index(["chrom", "pos38"])[["topld_ref", "topld_alt"]].to_dict("index") \
        if len(inf) else {}
    def sign_expected(r):
        # TOP-LD sign is assumed to be for ALT-ALT as listed in the TOP-LD info file.
        # Flip it once for each variant whose TOP-LD ALT is the GRCh38 REF.
        f = 1
        for pos, ref38, alt38 in ((r.pos1, r.ref1, r.alt1), (r.pos2, r.ref2, r.alt2)):
            a = topld_alleles.get((r.chrom, pos))
            if a is None:
                return np.nan
            if a["topld_alt"] == alt38 and a["topld_ref"] == ref38:
                continue
            if a["topld_alt"] == ref38 and a["topld_ref"] == alt38:
                f = -f
            else:
                return np.nan
        return f
    both["orient_mult"] = [sign_expected(r) for r in both.itertuples()]
    ok = both[both.orient_mult.notna()].copy()
    ok["topld_sign_num"] = np.where(ok.topld_corr_sign.str.strip() == "+", 1, -1)
    ok["topld_sign_grch38"] = ok.topld_sign_num * ok.orient_mult
    ok["panel_sign"] = np.sign(ok.panel_r)
    sign_conc = float((ok.topld_sign_grch38 == ok.panel_sign).mean()) if len(ok) else np.nan
    ok[["chrom", "v1", "v2", "topld_corr_sign", "orient_mult", "topld_sign_grch38",
        "panel_r", "panel_sign", "topld_r2", "panel_r2", "topld_dprime", "panel_Dprime"]].to_csv(
        TAB / "topld_sign_convention_check.tsv", sep="\t", index=False)
    log(f"TOP-LD sign concordance with phased panel: {sign_conc} on n={len(ok)}")

    # ---- arm plausibility
    m["dprime_best"] = m.topld_dprime.where(m.topld_dprime.notna(), m.panel_Dprime)
    m["dprime_source"] = np.where(m.topld_dprime.notna(), "topld",
                                  np.where(m.panel_Dprime.notna(), "panel_1kg_phase1_eur", "none"))
    m["arm_v1_freq"] = m.panel_f_alt1_ref2
    m["arm_v2_freq"] = m.panel_f_ref1_alt2
    m["arm_min_freq"] = m[["arm_v1_freq", "arm_v2_freq"]].min(axis=1)
    m["arm_near_unobserved"] = m.arm_min_freq < ARM_UNOBSERVED
    m["arm_exactly_zero"] = np.where(
        m.panel_n_alt1_ref2.isna(), np.nan,
        ((m.panel_n_alt1_ref2 == 0) | (m.panel_n_ref1_alt2 == 0)).astype(float))
    m["arm_near_unobserved"] = m.arm_near_unobserved.where(m.arm_min_freq.notna())
    m.to_csv(TAB / "pair_phase_substrate.tsv", sep="\t", index=False)

    ase = stage_f(pairs)

    # ---- GNMT
    p5only = pairs[~pairs.is_gnmt]
    gpos = {GNMT[0][1], GNMT[1][1]}
    gnmt = dict(
        pair=f"{GNMT[0][2]} chr6:{GNMT[0][1]} / {GNMT[1][2]} chr6:{GNMT[1][1]} (GRCh38)",
        in_p5_pair_table=bool(((p5only.chrom == "chr6") &
                               (p5only.pos1.isin(gpos) | p5only.pos2.isin(gpos))).any()),
        gene_gnmt_in_p5_pairs=int((p5only.gene.astype(str) == "GNMT").sum()),
        row=m[m.is_gnmt].to_dict("records"),
        topld_info=[inf[(inf.chrom == "chr6") & (inf.pos38 == p)].to_dict("records")
                    for _, p, _, _, _ in GNMT],
        liftover={rs: lv[(lv.chrom == c) & (lv.pos38 == p)].to_dict("records")
                  for c, p, rs, _, _ in GNMT},
    )
    (TAB / "gnmt_pair.json").write_text(json.dumps(gnmt, indent=1, default=str))

    # ---- F1.1 verdict: P5 pairs only, the GNMT row never enters a denominator
    m = m.copy()
    mp = m[~m.is_gnmt]
    n_all = len(mp)
    n_cov = int((mp.topld_pair_state == "covered").sum())
    n_panel = int(mp.panel_Dprime.notna().sum())
    def frac(mask, denom):
        return (int(mask.sum()), denom, float(mask.sum() / denom) if denom else np.nan)
    v_primary = frac((mp.dprime_best > 0.9).fillna(False), n_all)
    v_topld_only = frac((mp.topld_dprime > 0.9).fillna(False), n_cov)
    v_panel_only = frac((mp.panel_Dprime > 0.9).fillna(False), n_panel)
    v_topld_over_all = frac((mp.topld_dprime > 0.9).fillna(False), n_all)

    summary = dict(
        seed=SEED,
        n_pairs_table=n_all,
        n_pairs_scored_by_p5=int(mp.p5_scored.sum()),
        n_pairs_in_table_never_scored=int((~mp.p5_scored).sum()),
        n_distinct_variants=int(len(pd.unique(pd.concat([p5only.v1, p5only.v2])))),
        phase_retained_chroms=int(phase.explicitly_phased_hardcalls.sum()),
        phase_probe_chroms=len(phase),
        liftover=lv.allele_state.value_counts().to_dict(),
        liftover_crosswalk=lv.crosswalk_state.value_counts().to_dict(),
        topld_pair_state=mp.topld_pair_state.value_counts().to_dict(),
        panel_state=mp.panel_state.value_counts().to_dict(),
        topld_sign_concordance_with_phased_panel=sign_conc,
        topld_sign_check_n=int(len(ok)),
        F1_1_primary_all_pairs=dict(zip(("n_dprime_gt_0.9", "denominator", "fraction"), v_primary)),
        F1_1_topld_covered_only=dict(zip(("n_dprime_gt_0.9", "denominator", "fraction"), v_topld_only)),
        F1_1_panel_only=dict(zip(("n_dprime_gt_0.9", "denominator", "fraction"), v_panel_only)),
        F1_1_topld_numerator_over_all_pairs=dict(
            zip(("n_dprime_gt_0.9", "denominator", "fraction"), v_topld_over_all)),
        F1_1_met_primary=bool(v_primary[2] >= 0.40),
        arm_near_unobserved_n=int(mp.arm_near_unobserved.fillna(False).sum()),
        arm_near_unobserved_denominator=n_panel,
        arm_exactly_zero_n=int(mp.arm_exactly_zero.fillna(False).sum()),
        arm_min_freq_median=float(mp.arm_min_freq.median(skipna=True)),
        arm_min_freq_lt_0005=int((mp.arm_min_freq < 0.005).fillna(False).sum()),
        arm_near_unobserved_among_dprime_gt_09=int(
            (mp.arm_near_unobserved.fillna(False) & (mp.dprime_best > 0.9).fillna(False)).sum()),
        ase_pairs_both_sites_measured=ase.groupby("dataset").both_sites_measured.sum().to_dict(),
        ase_pairs_with_any_double_het_donor=ase[ase.n_donors_het_at_both > 0]
            .groupby("dataset").size().to_dict(),
        ase_pairs_le_500bp_both_sites_measured=ase[(ase.separation_bp <= 500)]
            .groupby("dataset").both_sites_measured.sum().to_dict(),
        ase_pairs_le_500bp_with_double_het=ase[(ase.separation_bp <= 500) &
                                               (ase.n_donors_het_at_both > 0)]
            .groupby("dataset").size().to_dict(),
    )
    (TAB / "f1_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    log(json.dumps(summary, indent=1, default=str))

    # ---- provenance
    inputs = [P5 / "haplotype_pairs.tsv", P5 / "haplotype_effects.tsv",
              P5 / "haplotype_prespec.json", P5 / "haplotype_summary.json",
              ASE / "allelic_donor_counts.tsv.gz", ASE / "gse244832_allelic_donor_counts.tsv.gz",
              ASE / "allelic_sites.tsv", ASE / "gse244832_allelic_sites.tsv",
              CHAIN, PROJ / "scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md",
              PROJ / "scripts/analysis/alphagenome_program/f1_phase.py"]
    inputs += [KG / f"{c}_eur.pgen" for c in CHROMS]
    inputs += [TOPLD / f"EUR_{c}_Dprime_extended_0.2_1000000.csv.gz" for c in CHROMS]
    inputs += [TOPLD / f"EUR_{c}_Dprime_extended_0.2_1000000_info.csv.gz" for c in CHROMS]
    mrows = []
    for p in inputs:
        if p.exists():
            mrows.append(dict(path=str(p), bytes=p.stat().st_size, sha256=sha256(p)))
        else:
            mrows.append(dict(path=str(p), bytes=-1, sha256="MISSING"))
    pd.DataFrame(mrows).to_csv(PROV / "MANIFEST.tsv", sep="\t", index=False)
    with open(PROV / "pip_freeze.txt", "w") as fh:
        fh.write(subprocess.run([sys.executable, "-m", "pip", "freeze"],
                                capture_output=True, text=True).stdout)
    (PROV / "plink2_version.txt").write_text(
        subprocess.run([PLINK2, "--version"], capture_output=True, text=True).stdout)
    log("done")


if __name__ == "__main__":
    main()
