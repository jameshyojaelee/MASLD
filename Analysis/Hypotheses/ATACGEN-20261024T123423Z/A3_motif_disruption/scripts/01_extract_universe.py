#!/usr/bin/env python3
"""A3 step 1: variant universe from the study-level SuSiE fits (hg19).

Universe = biallelic SNVs in study-by-locus windows that hold >= 1 gated credible set
(run_19783321). Per variant: max PIP over every row of susie_all_results.csv, any CS
membership, max gated-member PIP, assigned study, effect allele and corrected beta.
"""
import hashlib
import os
import sys

import numpy as np
import pandas as pd

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
FM = f"{BASE}/GWAS/finemapping"
OUT = os.environ.get("A3_OUT", f"{BASE}/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption")
MEMBERS = f"{FM}/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv"
SUSIE_ALL = os.environ.get("A3_SUSIE_ALL", f"{FM}/results/susie_all_results.csv")
REGISTRY = f"{FM}/config/gwas_registry.tsv"
R2 = f"{FM}/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
PANEL_AF = f"{FM}/data/ld_ref/panel_af"
PAZOKI = "2022_36402844_PDFF_EUR"  # SuSiE fits predate the 2026-10-10 label fix: negate beta

BASES = {"A", "C", "G", "T"}
PAIR_CODE = {("A", "C"): 0, ("A", "G"): 1, ("A", "T"): 2, ("C", "G"): 3, ("C", "T"): 4, ("G", "T"): 5}
CODE_PAIR = {v: k for k, v in PAIR_CODE.items()}


def sha256(path, block=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def pair_code(a1, a2):
    lo = np.where(a1 < a2, a1, a2)
    hi = np.where(a1 < a2, a2, a1)
    s = pd.Series(lo) + pd.Series(hi)
    m = {"AC": 0, "AG": 1, "AT": 2, "CG": 3, "CT": 4, "GT": 5}
    return s.map(m).to_numpy()


def make_key(chrom, pos, code):
    return (chrom.astype(np.int64) * 1_000_000_000 + pos.astype(np.int64)) * 8 + code.astype(np.int64)


def main():
    os.makedirs(f"{OUT}/work", exist_ok=True)
    os.makedirs(f"{OUT}/results", exist_ok=True)

    reg = pd.read_csv(REGISTRY, sep="\t", dtype=str)
    anc = dict(zip(reg.study_name, reg.ancestry))

    mem = pd.read_csv(MEMBERS, sep="\t", dtype=str)
    mem = mem[mem.pip_sum_gate_passed.str.lower() == "true"].copy()
    print(f"gated member rows: {len(mem)}; studies: {mem.study.nunique()}; "
          f"CS instances: {mem.credible_set_uid.nunique()}", flush=True)

    r2_studies = set(pd.read_csv(R2, usecols=["gwas_name"], dtype=str).gwas_name.unique())
    missing = sorted(set(mem.study) - r2_studies)
    print(f"r2 studies: {len(r2_studies)}; member studies not in r2: {missing}", flush=True)
    if missing:
        sys.exit("member study missing from COLOC r2")

    mem["susie_pip"] = mem.susie_pip.astype(float)
    snv = mem.allele1.isin(BASES) & mem.allele2.isin(BASES)
    mem = mem[snv & mem.chromosome.isin([str(i) for i in range(1, 23)])].copy()
    mem["chr"] = mem.chromosome.astype(int)
    mem["pos19"] = mem.position.astype(int)
    mem["code"] = pair_code(mem.allele1.to_numpy(), mem.allele2.to_numpy())
    mem["key"] = make_key(mem.chr.to_numpy(), mem.pos19.to_numpy(), mem.code.to_numpy())
    gated_windows = set(zip(mem.study, mem.locus))
    mem_rowkeys = set(zip(mem.study, mem.locus, mem.key))
    print(f"gated SNV member rows: {len(mem)}; unique SNVs: {mem.key.nunique()}; "
          f"gated windows: {len(gated_windows)}", flush=True)

    usecols = ["chromosome", "position", "allele1", "allele2", "beta", "locus", "PIP", "CS", "study"]
    dtypes = {"chromosome": str, "position": "Int64", "allele1": str, "allele2": str,
              "beta": float, "locus": str, "PIP": float, "CS": float, "study": str}
    agg_parts = []
    member_rows = []
    window_keys = []
    n_rows = 0
    for i, ch in enumerate(pd.read_csv(SUSIE_ALL, usecols=usecols, dtype=dtypes, chunksize=5_000_000)):
        n_rows += len(ch)
        ch = ch[ch.chromosome.isin([str(c) for c in range(1, 23)])]
        ch = ch[ch.allele1.isin(BASES) & ch.allele2.isin(BASES) & ch.position.notna()]
        chrom = ch.chromosome.astype(int).to_numpy()
        pos = ch.position.astype(np.int64).to_numpy()
        code = pair_code(ch.allele1.to_numpy(), ch.allele2.to_numpy())
        key = make_key(chrom, pos, code)
        ch = ch.assign(key=key)
        in_window = pd.MultiIndex.from_arrays([ch.study, ch.locus]).isin(list(gated_windows))
        window_keys.append(np.unique(key[in_window]))
        # rows that are gated members: keep effect allele and beta
        is_mem = pd.MultiIndex.from_arrays([ch.study, ch.locus, ch.key]).isin(list(mem_rowkeys))
        if is_mem.any():
            member_rows.append(ch.loc[is_mem, ["key", "study", "locus", "allele1", "allele2", "beta", "PIP", "CS"]])
        cs_any = (ch.CS.fillna(0) > 0).to_numpy()
        part = pd.DataFrame({"key": key, "pip": ch.PIP.to_numpy(), "cs_any": cs_any,
                             "pip_nonmiss": ch.PIP.notna().to_numpy(), "study": ch.study.to_numpy()})
        part = part.sort_values("pip", ascending=False, na_position="last")
        g = part.groupby("key", sort=False)
        agg = pd.DataFrame({"max_pip_all": g.pip.max(), "cs_any": g.cs_any.any(),
                            "n_nonmiss": g.pip_nonmiss.sum(), "study_of_max": g.study.first()})
        agg_parts.append(agg.reset_index())
        print(f"chunk {i}: rows so far {n_rows}", flush=True)

    wk = np.unique(np.concatenate(window_keys))
    print(f"total rows read: {n_rows}; unique SNV keys in gated windows: {len(wk)}", flush=True)
    agg = pd.concat(agg_parts, ignore_index=True)
    agg = agg[agg.key.isin(wk)]
    agg = agg.sort_values("max_pip_all", ascending=False, na_position="last")
    g = agg.groupby("key", sort=False)
    uni = pd.DataFrame({"max_pip_all": g.max_pip_all.max(), "cs_any": g.cs_any.any(),
                        "n_nonmiss": g.n_nonmiss.sum(), "study_of_max": g.study_of_max.first()}).reset_index()
    del agg, agg_parts

    mr = pd.concat(member_rows, ignore_index=True)
    mr = mr.merge(mem[["study", "locus", "key", "susie_pip", "tier", "trait", "ancestry"]],
                  on=["study", "locus", "key"], how="inner")
    print(f"member rows matched in susie_all_results: {len(mr)} of {len(mem)}", flush=True)
    mr["beta_corr"] = np.where(mr.study == PAZOKI, -mr.beta, mr.beta)
    mr["risk_allele"] = np.where(mr.beta_corr > 0, mr.allele1, mr.allele2)
    mr.loc[mr.beta_corr == 0, "risk_allele"] = None

    # orientation QC at rs738409 (hg19 22:44324727 C/G); G is the PDFF/ALT-raising allele
    rs = mr[mr.key == make_key(np.array([22]), np.array([44324727]), np.array([3]))[0]]
    rs = rs.drop_duplicates("study")[["study", "allele1", "allele2", "beta", "beta_corr", "risk_allele", "tier", "trait"]]
    rs.to_csv(f"{OUT}/results/orientation_check_rs738409.tsv", sep="\t", index=False)
    print(rs.to_string(), flush=True)
    bad = rs[(rs.tier == "1") & (rs.risk_allele != "G")]
    if len(bad):
        sys.exit(f"rs738409 orientation failure in direct-trait studies: {bad.study.tolist()}")

    mr = mr.sort_values("susie_pip", ascending=False)
    gm = mr.groupby("key", sort=False)
    eur = mr[mr.ancestry == "EUR"].groupby("key").susie_pip.max().rename("max_gated_pip_eur")
    risk_n = mr.dropna(subset=["risk_allele"]).groupby("key").risk_allele.nunique().rename("n_distinct_risk_alleles")
    memv = pd.DataFrame({
        "max_gated_pip": gm.susie_pip.max(),
        "n_gated_rows": gm.size(),
        "n_gated_studies": gm.study.nunique(),
        "assigned_study": gm.study.first(),
        "assigned_tier": gm.tier.first(),
        "assigned_trait": gm.trait.first(),
        "effect_allele": gm.allele1.first(),
        "other_allele": gm.allele2.first(),
        "beta_corr": gm.beta_corr.first(),
        "risk_allele": gm.risk_allele.first(),
    }).join(eur).join(risk_n).reset_index()

    uni = uni.merge(memv, on="key", how="left")
    uni["is_gated_member"] = uni.max_gated_pip.notna()
    uni["vclass"] = "other"
    uni.loc[uni.is_gated_member & (uni.max_gated_pip >= 0.1), "vclass"] = "high"
    uni.loc[uni.is_gated_member & (uni.max_gated_pip < 0.1), "vclass"] = "member_low"
    null = (~uni.is_gated_member) & (~uni.cs_any) & (uni.n_nonmiss > 0) & (uni.max_pip_all < 0.01)
    uni.loc[null, "vclass"] = "null"
    uni.loc[uni.vclass == "null", "assigned_study"] = uni.loc[uni.vclass == "null", "study_of_max"]
    uni["assigned_ancestry"] = uni.assigned_study.map(anc)

    k = uni.key.to_numpy()
    code = k % 8
    rest = k // 8
    uni["chr"] = rest // 1_000_000_000
    uni["pos19"] = rest % 1_000_000_000
    pairs = np.array([CODE_PAIR[c] for c in range(6)])
    uni["a_lo"] = pairs[code, 0]
    uni["a_hi"] = pairs[code, 1]

    # panel MAF in the assigned study's ancestry (sensitivity B only)
    uni["maf"] = np.nan
    popmap = {"EUR": "eur", "EAS": "eas", "AFR": "afr", "AMR": "amr", "SAS": "sas"}
    for a, pop in popmap.items():
        sel = uni.assigned_ancestry == a
        if not sel.any():
            continue
        need = uni.loc[sel, ["key", "chr"]]
        for c in sorted(need.chr.unique()):
            f = f"{PANEL_AF}/{pop}/chr{c}.af.tsv.gz"
            if not os.path.exists(f):
                continue
            p = pd.read_csv(f, sep="\t", usecols=["chromosome", "position", "bim_a1", "bim_a2", "af_a1"],
                            dtype={"chromosome": int, "position": np.int64, "bim_a1": str, "bim_a2": str, "af_a1": float})
            p = p[p.bim_a1.isin(BASES) & p.bim_a2.isin(BASES)]
            pk = make_key(p.chromosome.to_numpy(), p.position.to_numpy(), pair_code(p.bim_a1.to_numpy(), p.bim_a2.to_numpy()))
            maf = np.minimum(p.af_a1.to_numpy(), 1 - p.af_a1.to_numpy())
            s = pd.Series(maf, index=pk)
            s = s[~s.index.duplicated()]
            idx = need.index[need.chr == c]
            uni.loc[idx, "maf"] = s.reindex(uni.loc[idx, "key"].to_numpy()).to_numpy()
        print(f"MAF {a}: {uni.loc[sel, 'maf'].notna().mean():.3f} non-missing", flush=True)

    uni = uni.drop(columns=["study_of_max"])
    uni.to_csv(f"{OUT}/work/universe_hg19.tsv.gz", sep="\t", index=False, compression="gzip")
    summ = uni.groupby("vclass").agg(n=("key", "size")).reset_index()
    summ.to_csv(f"{OUT}/results/universe_class_counts_hg19.tsv", sep="\t", index=False)
    print(summ.to_string(), flush=True)

    with open(f"{OUT}/results/input_sha256.tsv", "a") as fh:
        for p in [MEMBERS, SUSIE_ALL, REGISTRY, R2]:
            fh.write(f"{p}\t{sha256(p)}\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
