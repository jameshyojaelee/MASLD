#!/usr/bin/env python3
"""DIR6 step 1 (PRESPEC_DIR6.md sections 3, 4, 6, 7): raw GWAS rows at the H26 instruments.

For the 286 H26 calls: the instrument (GWAS lead of the defining trait, trait-raising allele) from the
direction table; liver-fat (Pazoki 2022, corrected copy, GRCh37) and cirrhosis (Ghouse 2024, GRCh38)
rows, direct or by 1kGP EUR proxy at r2 >= 0.8; Pan-UKBB EUR ALT/AST/GGT; co-effect rows (GLGC TG,
Yengo BMI, Manning HOMA-IR). Control variants are added. Allele orientation asserts run first.

Usage: python -I 01_extract.py <project_root> <out_dir>
Output: instruments.tsv, rows.tsv (one row per id x trait), proxies.tsv, controls_raw.tsv
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])
MR = ROOT / "GWAS/MR_Data"
H26 = ROOT / "Analysis/Hypotheses/H26-20261007T210504Z/run_22092011/h26_gene_direction_table.tsv"
DIRT = ROOT / ("Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/"
               "coloc_direction/full-20260930T120450Z/direction_by_trait.tsv")
RELEASE = "candidate_r2_20260924T053244Z"
PDFF = MR / "GCST90267352_PDFF_Pazoki2022.tsv.gz"
GHOUSE = MR / "Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz"
PANUKBB = {"ALT": "30620-ALT", "AST": "30650-AST", "GGT": "30730-GGT"}
GLGC = MR / "GLGC/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz"
YENGO = MR / "Yengo_2018_BMI_GWAS.txt.gz"
MANNING = MR / "MAGIC/MAGIC_Manning_et_al_HOMAIR_MainEffect.txt.gz"
KGP = ROOT / "data/external/allelic_refs/kgp_highcov"
KING = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/"
            "executions/in-biopsy-allelic-tags-v2-20260929T134419Z/founders/exclude_related.txt")
GENELOC = ROOT / "tools/magma/ref/NCBI37.3.gene.loc"
WINDOW = 250_000
R2_MIN = 0.8
COMP = str.maketrans("ACGT", "TGCA")

# control variants (GRCh37 positions as in ENZFAT1 20_extract_loci.py and PRESPEC_ENZFAT1.md)
CONTROLS = [("ctrl_PNPLA3_rs738409", "22", 44324727, "rs738409"),
            ("ctrl_TM6SF2_rs58542926", "19", 19379549, "rs58542926"),
            ("ctrl_GCKR_rs1260326", "2", 27730940, "rs1260326"),
            ("ctrl_HSD17B13_rs72613567", "4", 88231393, "rs72613567"),
            ("ctrl_HSD17B13_rs6834314", "4", 88213808, "rs6834314")]


def log(*a):
    print(*a, flush=True)


def pair_eq(a1, a2, b1, b2):
    return ((a1 == b1) & (a2 == b2)) | ((a1 == b2) & (a2 == b1))


def palindromic(a, b):
    return {a, b} in ({"A", "T"}, {"C", "G"})


def scan(path, usecols, keyfn, wanted, sep="\t", **kw):
    hdr = pd.read_csv(path, sep=sep, nrows=0, **kw).columns
    names = [c.lstrip("#") for c in hdr]
    out = []
    for ch in pd.read_csv(path, sep=sep, names=names, header=0, usecols=usecols, dtype=str,
                          chunksize=2_000_000, **kw):
        ch = ch.assign(key=keyfn(ch))
        out.append(ch[ch["key"].isin(wanted)])
    return pd.concat(out, ignore_index=True)


# ----------------------------------------------------------------------------- instruments
def instruments():
    h = pd.read_csv(H26, sep="\t", dtype=str)
    h = h[h["lineA_status"] == "call"]
    t = pd.read_csv(DIRT, sep="\t", dtype=str, low_memory=False)
    t = t[t["release"] == RELEASE]
    cols = ["ensembl", "trait", "hit1", "hit1_grch38", "gwas_allele1", "gwas_allele2", "trait_raising_allele",
            "reported_study", "study_superpop", "direction_sign"]
    m = h.merge(t[cols], left_on=["ensembl", "lineA_defining_trait"], right_on=["ensembl", "trait"], how="left")
    assert len(m) == 286 and m["hit1"].notna().all(), "instrument join failed"
    assert (m["lineA_call"].value_counts().to_dict() == {"inhibit": 149, "activate": 137})
    sign = m["direction_sign"].astype(float)
    assert ((sign > 0) == (m["lineA_call"] == "inhibit")).all()
    m["chr"] = m["hit1"].str.split(":").str[0]
    m["pos37"] = m["hit1"].str.split(":").str[1].astype(int)
    g38 = m["hit1_grch38"].str.split(":")
    assert (g38.str[2] == "+").all()
    m["pos38"] = g38.str[1].astype(int)
    m["a_raise"] = m["trait_raising_allele"].str.upper()
    a1, a2 = m["gwas_allele1"].str.upper(), m["gwas_allele2"].str.upper()
    assert ((m["a_raise"] == a1) | (m["a_raise"] == a2)).all()
    m["a_other"] = np.where(m["a_raise"] == a1, a2, a1)
    m["vid"] = m["chr"] + ":" + m["pos37"].astype(str) + ":" + m["a_raise"] + ":" + m["a_other"]
    keep = ["ensembl", "symbol", "primary_evidence_class", "lineA_call", "lineA_defining_trait", "lineA_ggt_only",
            "published_coding_mechanism", "reported_study", "study_superpop", "chr", "pos37", "pos38", "a_raise",
            "a_other", "vid"]
    return m[keep].rename(columns={"lineA_call": "call", "lineA_defining_trait": "defining_trait"})


# ----------------------------------------------------------------------------- 1kGP LD
def eur_samples():
    p = pd.read_csv(KGP / "20130606_g1k_3202_samples_ped_population.txt", sep=r"\s+", dtype=str)
    f = p[(p["Superpopulation"] == "EUR") & (p["FatherID"] == "0") & (p["MotherID"] == "0")]["SampleID"]
    excl = set(KING.read_text().split())
    return [s for s in f if s not in excl]


def window_haps(chrom, pos, samples, work):
    """Biallelic records within +/- WINDOW: list of (pos, ref, alt, haplotype vector ALT=1, -1 missing)."""
    vcf = KGP / f"1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    sfile = work / "eur_samples.txt"
    if not sfile.exists():
        sfile.write_text("\n".join(samples) + "\n")
    reg = f"chr{chrom}:{max(1, pos - WINDOW)}-{pos + WINDOW}"
    cmd = ["bcftools", "query", "-r", reg, "-S", str(sfile), "--force-samples", "-i", "N_ALT=1",
           "-f", "%POS\t%REF\t%ALT\t[%GT,]\n", str(vcf)]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines()
    out = []
    for line in res:
        p, ref, alt, gts = line.split("\t")
        b = np.frombuffer(gts.encode(), dtype=np.uint8)
        if len(b) % 4:  # non-diploid or unexpected genotype strings: skip the record
            continue
        b = b.reshape(-1, 4)
        phased = b[:, 1] == ord("|")
        a = b[:, [0, 2]].astype(np.int16) - ord("0")
        ok = phased[:, None] & (a >= 0) & (a <= 1)
        h = np.where(ok, a, -1).astype(np.int8).reshape(-1)
        out.append((int(p), ref.upper(), alt.upper(), h))
    return out


def r_of(x, y):
    ok = (x >= 0) & (y >= 0)
    if ok.sum() < 50:
        return np.nan
    x, y = x[ok].astype(float), y[ok].astype(float)
    if x.std() == 0 or y.std() == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def proxies_for(inst, samples, work):
    """Instrument EUR frequency and candidate proxies (r2 >= R2_MIN) with the allele on the instrument
    allele's haplotype."""
    recs = window_haps(inst.chr, inst.pos38, samples, work)
    x, eaf = None, np.nan
    for pos, ref, alt, h in recs:
        if pos != inst.pos38:
            continue
        if (inst.a_raise, inst.a_other) == (alt, ref):
            x = h.copy()
        elif (inst.a_raise, inst.a_other) == (ref, alt):
            x = h.copy(); x[h >= 0] = 1 - h[h >= 0]
        if x is not None:
            eaf = float(x[x >= 0].mean())
            break
    if x is None:
        return eaf, "instrument allele pair not in 1kGP", []
    cands = []
    for pos, ref, alt, h in recs:
        if pos == inst.pos38:
            continue
        af = h[h >= 0].mean() if (h >= 0).any() else np.nan
        if not np.isfinite(af) or min(af, 1 - af) < 0.01:
            continue
        r = r_of(x, h)
        if not np.isfinite(r) or r * r < R2_MIN:
            continue
        al, ol = (alt, ref) if r > 0 else (ref, alt)
        cands.append(dict(pos38=pos, ref=ref, alt=alt, r=r, r2=r * r, dist=abs(pos - inst.pos38),
                          p_al=al, p_ol=ol, p_maf=min(af, 1 - af)))
    cands.sort(key=lambda d: (-round(d["r2"], 6), d["dist"]))
    return eaf, "ok", cands


# ----------------------------------------------------------------------------- main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    work = OUT / "work"; work.mkdir(exist_ok=True)
    inst = instruments()
    uniq = inst.drop_duplicates("vid")[["vid", "chr", "pos37", "pos38", "a_raise", "a_other"]].reset_index(drop=True)
    log(f"286 calls, {len(uniq)} unique instrument variant-allele rows, {inst[['chr', 'pos37']].drop_duplicates().shape[0]} unique positions")

    # ---- liver fat: instrument positions (GRCh37) and control rsIDs; also GPT/GOT1 windows later
    key37 = set(uniq["chr"] + ":" + uniq["pos37"].astype(str)) | {f"{c}:{p}" for _, c, p, _ in CONTROLS}
    pdff_cols = ["variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error"]
    f = scan(PDFF, pdff_cols, lambda d: d["chromosome"] + ":" + d["base_pair_location"], key37)
    f = f.rename(columns={"variant_id": "rsid", "effect_allele": "ea", "other_allele": "oa", "standard_error": "se"})
    f["ea"], f["oa"] = f["ea"].str.upper(), f["oa"].str.upper()
    log(f"liver-fat rows at instrument/control positions: {len(f)}")

    # ---- cirrhosis: instrument GRCh38 positions and control rsIDs
    key38 = set(uniq["chr"] + ":" + uniq["pos38"].astype(str))
    ctrl_rs = {r for *_, r in CONTROLS}
    g_cols = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "rsid"]
    g = scan(GHOUSE, g_cols, lambda d: np.where(d["rsid"].isin(ctrl_rs), d["rsid"], d["chromosome"] + ":" + d["base_pair_location"]),
             key38 | ctrl_rs)
    g = g.rename(columns={"effect_allele": "ea", "other_allele": "oa", "standard_error": "se"})
    g["ea"], g["oa"] = g["ea"].str.upper(), g["oa"].str.upper()
    log(f"cirrhosis rows at instrument positions/control rsIDs: {len(g)}")

    # ---- Pan-UKBB (GRCh37) at instruments and controls; GPT/GOT1 windows for control leads
    gl = pd.read_csv(GENELOC, sep=r"\s+", header=None, names=["id", "chr", "start", "end", "strand", "gene"], dtype=str)
    win = {}
    for gene in ["GPT", "GOT1"]:
        r = gl[gl["gene"] == gene].iloc[0]
        win[gene] = (r["chr"], int(r["start"]) - 100_000, int(r["end"]) + 100_000)
    winkeys = {f"{c}:{p}" for c, a, b in win.values() for p in range(a, b + 1)}
    pk = key37 | winkeys
    pan = {}
    for tr, code in PANUKBB.items():
        u = scan(MR / f"PanUKBB/biomarkers-{code}-both_sexes-irnt.tsv.bgz",
                 ["chr", "pos", "ref", "alt", "af_EUR", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "low_confidence_EUR"],
                 lambda d: d["chr"] + ":" + d["pos"], pk, compression="gzip")
        u = u[u["low_confidence_EUR"].str.lower() == "false"].rename(
            columns={"alt": "ea", "ref": "oa", "beta_EUR": "beta", "se_EUR": "se"})
        pan[tr] = u
        log(f"Pan-UKBB {tr} rows: {len(u)}")

    # GPT (ALT) and GOT1 (AST) control leads: smallest p, MAF >= 1%, present in the liver-fat file with the same pair
    fw = scan(PDFF, pdff_cols, lambda d: d["chromosome"] + ":" + d["base_pair_location"], winkeys)
    fw = fw.rename(columns={"variant_id": "rsid", "effect_allele": "ea", "other_allele": "oa", "standard_error": "se"})
    fw["ea"], fw["oa"] = fw["ea"].str.upper(), fw["oa"].str.upper()
    extra_ctrl = []
    for gene, tr in [("GPT", "ALT"), ("GOT1", "AST")]:
        c, a, b = win[gene]
        u = pan[tr]
        u = u[(u["chr"] == c) & u["pos"].astype(int).between(a, b)].dropna(subset=["af_EUR", "neglog10_pval_EUR"])
        af = u["af_EUR"].astype(float)
        u = u[(af >= 0.01) & (af <= 0.99)]
        mm = u.merge(fw[["key", "rsid", "ea", "oa"]], on="key", suffixes=("", "_f"))
        mm = mm[pair_eq(mm["ea"], mm["oa"], mm["ea_f"], mm["oa_f"])]
        mm = mm[~mm["key"].duplicated(keep=False)]
        lead = mm.loc[mm["neglog10_pval_EUR"].astype(float).idxmax()]
        log(f"{gene} control lead: {lead['key']} {lead['rsid']} {lead['oa']}/{lead['ea']} {tr} -log10p {lead['neglog10_pval_EUR']}")
        extra_ctrl.append((f"ctrl_{gene}_lead", c, int(lead["pos"]), lead["rsid"]))
    f = pd.concat([f, fw[fw["key"].isin({f"{c}:{p}" for _, c, p, _ in extra_ctrl})]], ignore_index=True)
    ctrl_all = CONTROLS + extra_ctrl
    ctrl_rs2 = {r for *_, r in extra_ctrl}
    g2 = scan(GHOUSE, g_cols, lambda d: d["rsid"], ctrl_rs2)
    g2 = g2.rename(columns={"effect_allele": "ea", "other_allele": "oa", "standard_error": "se"})
    g2["ea"], g2["oa"] = g2["ea"].str.upper(), g2["oa"].str.upper()
    g = pd.concat([g, g2], ignore_index=True)

    rows = []

    def add(id_, trait, source, ea, oa, beta, se, extra=None):
        d = dict(id=id_, trait=trait, source=source, ea=ea, oa=oa, beta=float(beta), se=float(se))
        d.update(extra or {})
        rows.append(d)

    # ---- controls: rows by GRCh37 position (liver fat, Pan-UKBB) or rsID (cirrhosis), alleles from the liver-fat row
    ctrl_rows = []
    for cid, c, p, rs in ctrl_all:
        fr = f[(f["key"] == f"{c}:{p}") & (f["rsid"] == rs)]
        if len(fr) != 1:
            log(f"control {cid}: {len(fr)} liver-fat rows at {c}:{p} {rs}"); continue
        fr = fr.iloc[0]
        add(cid, "pdff", "direct", fr["ea"], fr["oa"], fr["beta"], fr["se"])
        gr = g[(g["rsid"] == rs) & pair_eq(g["ea"], g["oa"], fr["ea"], fr["oa"])]
        if len(gr) == 1:
            add(cid, "cirr", "direct", gr.iloc[0]["ea"], gr.iloc[0]["oa"], gr.iloc[0]["beta"], gr.iloc[0]["se"])
        for tr in PANUKBB:
            u = pan[tr]
            ur = u[(u["key"] == f"{c}:{p}") & pair_eq(u["ea"], u["oa"], fr["ea"], fr["oa"])]
            if len(ur) == 1:
                add(cid, tr, "direct", ur.iloc[0]["ea"], ur.iloc[0]["oa"], ur.iloc[0]["beta"], ur.iloc[0]["se"])
        ctrl_rows.append(dict(id=cid, chr=c, pos37=p, rsid=rs))
    pd.DataFrame(ctrl_rows).to_csv(OUT / "controls_variants.tsv", sep="\t", index=False)

    # ---- allele orientation asserts (section 7.1) on raw rows
    R = pd.DataFrame(rows)

    def eff(cid, trait, allele):
        r = R[(R["id"] == cid) & (R["trait"] == trait)]
        assert len(r) == 1, f"assert row missing: {cid} {trait}"
        r = r.iloc[0]
        return r["beta"] if r["ea"] == allele else -r["beta"]

    checks = [("ctrl_PNPLA3_rs738409", "pdff", "G"), ("ctrl_PNPLA3_rs738409", "cirr", "G"),
              ("ctrl_TM6SF2_rs58542926", "pdff", "T"), ("ctrl_GCKR_rs1260326", "pdff", "T"),
              ("ctrl_PNPLA3_rs738409", "ALT", "G")]
    ck = []
    for cid, tr, a in checks:
        v = eff(cid, tr, a)
        ck.append(dict(variant=cid, trait=tr, allele=a, effect=v, pass_=v > 0))
    ck = pd.DataFrame(ck).rename(columns={"pass_": "pass"})
    ck.to_csv(OUT / "orientation_asserts.tsv", sep="\t", index=False)
    log(ck.to_string())
    if not ck["pass"].all():
        sys.exit("ALLELE ORIENTATION ASSERT FAILED; stopping")

    # ---- instruments: direct rows
    usable = {}  # vid -> {"pdff": ..., "cirr": ...} reason or row
    pan_af = {}
    for u in uniq.itertuples(index=False):
        k37, k38 = f"{u.chr}:{u.pos37}", f"{u.chr}:{u.pos38}"
        for tr in PANUKBB:
            p = pan[tr]
            pr = p[(p["key"] == k37) & pair_eq(p["ea"], p["oa"], u.a_raise, u.a_other)]
            if len(pr) == 1:
                add(u.vid, tr, "direct", pr.iloc[0]["ea"], pr.iloc[0]["oa"], pr.iloc[0]["beta"], pr.iloc[0]["se"])
                a = float(pr.iloc[0]["af_EUR"]) if pr.iloc[0]["af_EUR"] == pr.iloc[0]["af_EUR"] else np.nan
                pan_af[u.vid] = a if pr.iloc[0]["ea"] == u.a_raise else 1 - a
        fr = f[(f["key"] == k37) & pair_eq(f["ea"], f["oa"], u.a_raise, u.a_other)]
        gr = g[(g["key"] == k38) & pair_eq(g["ea"], g["oa"], u.a_raise, u.a_other)]
        usable[u.vid] = dict(pdff=fr.iloc[0] if len(fr) == 1 else None, cirr=gr.iloc[0] if len(gr) == 1 else None,
                             rsid=(fr.iloc[0]["rsid"] if len(fr) == 1 else (gr.iloc[0]["rsid"] if len(gr) == 1 else None)))

    # ---- 1kGP EUR frequency (palindrome rule) and proxies for every unique instrument
    samples = eur_samples()
    log(f"1kGP EUR founders minus KING-related: {len(samples)}")
    prox, info = {}, []
    for u in uniq.itertuples(index=False):
        eaf, status, cands = proxies_for(u, samples, work)
        prox[u.vid] = cands
        fr_af = eaf if np.isfinite(eaf) else pan_af.get(u.vid, np.nan)
        pal = palindromic(u.a_raise, u.a_other)
        pal_block = pal and (not np.isfinite(fr_af) or min(fr_af, 1 - fr_af) > 0.4)
        info.append(dict(vid=u.vid, kgp_status=status, eaf_raise_eur=eaf, eaf_raise_panukbb=pan_af.get(u.vid, np.nan),
                         palindromic=pal, palindrome_blocked=pal_block, n_proxy_candidates=len(cands)))
    info = pd.DataFrame(info).set_index("vid")
    log(f"instruments in 1kGP: {(info['kgp_status'] == 'ok').sum()}; palindrome-blocked: {info['palindrome_blocked'].sum()}")

    # ---- second pass for proxies where needed
    need = {}
    for vid, d in usable.items():
        blocked = info.loc[vid, "palindrome_blocked"]
        for tr in ["pdff", "cirr"]:
            if d[tr] is None or blocked:
                need[(vid, tr)] = [c for c in prox[vid] if not (palindromic(c["p_al"], c["p_ol"]) and c["p_maf"] > 0.4)]
    chrom_of = uniq.set_index("vid")["chr"].to_dict()
    pk38 = {f"{chrom_of[v]}:{c['pos38']}" for (v, _), cs in need.items() for c in cs}
    gp = scan(GHOUSE, g_cols, lambda d: d["chromosome"] + ":" + d["base_pair_location"], pk38) if pk38 else g.iloc[:0]
    gp = gp.rename(columns={"effect_allele": "ea", "other_allele": "oa", "standard_error": "se"})
    gp["ea"], gp["oa"] = gp["ea"].str.upper(), gp["oa"].str.upper()
    prs = set(gp["rsid"].dropna())
    fp = scan(PDFF, pdff_cols, lambda d: d["variant_id"], prs) if prs else f.iloc[:0]
    fp = fp.rename(columns={"variant_id": "rsid", "effect_allele": "ea", "other_allele": "oa", "standard_error": "se"})
    fp["ea"], fp["oa"] = fp["ea"].str.upper(), fp["oa"].str.upper()
    log(f"proxy pass: {len(need)} instrument x outcome gaps, {len(pk38)} candidate positions, "
        f"{len(gp)} Ghouse rows, {len(fp)} liver-fat rows by rsID")

    plog = []
    for u in uniq.itertuples(index=False):
        d = usable[u.vid]
        blocked = info.loc[u.vid, "palindrome_blocked"]
        for tr in ["pdff", "cirr"]:
            if d[tr] is not None and not blocked:
                r = d[tr]
                add(u.vid, tr, "direct", r["ea"], r["oa"], r["beta"], r["se"], dict(proxy_r2=1.0))
                continue
            chosen = None
            for c in need.get((u.vid, tr), []):
                k = f"{u.chr}:{c['pos38']}"
                grow = gp[(gp["key"] == k) & pair_eq(gp["ea"], gp["oa"], c["p_al"], c["p_ol"])]
                if len(grow) != 1:
                    continue
                if tr == "cirr":
                    chosen = (c, grow.iloc[0]); break
                rs = grow.iloc[0]["rsid"]
                frow = fp[(fp["rsid"] == rs) & pair_eq(fp["ea"], fp["oa"], c["p_al"], c["p_ol"])]
                if len(frow) == 1:
                    chosen = (c, frow.iloc[0]); break
            if chosen is None:
                reason = ("palindromic, EUR MAF > 0.4" if blocked else "variant not in outcome file") + \
                         ("; no usable proxy at r2 >= 0.8" if info.loc[u.vid, "kgp_status"] == "ok" else
                          f"; {info.loc[u.vid, 'kgp_status']}")
                plog.append(dict(vid=u.vid, trait=tr, status="untestable", reason=reason))
                continue
            c, r = chosen
            # outcome effect of the proxy allele on the instrument allele's haplotype
            b = float(r["beta"]) if r["ea"] == c["p_al"] else -float(r["beta"])
            add(u.vid, tr, "proxy", u.a_raise, u.a_other, b, r["se"],
                dict(proxy_r2=c["r2"], proxy_pos38=c["pos38"], proxy_allele=c["p_al"], proxy_rsid=r.get("rsid")))
            plog.append(dict(vid=u.vid, trait=tr, status="proxy", reason=f"r2 {c['r2']:.3f} at {u.chr}:{c['pos38']} (GRCh38)"))
    pd.DataFrame(plog).to_csv(OUT / "proxies.tsv", sep="\t", index=False)

    # ---- co-effects at the instrument (no proxy)
    k37u = set(uniq["chr"] + ":" + uniq["pos37"].astype(str))
    gl_ = scan(GLGC, ["rsID", "CHROM", "POS_b37", "REF", "ALT", "METAL_Effect", "METAL_StdErr"],
               lambda d: d["CHROM"].str.strip() + ":" + d["POS_b37"].str.strip(), k37u)
    y = scan(YENGO, ["CHR", "POS", "SNP", "Tested_Allele", "Other_Allele", "BETA", "SE"],
             lambda d: d["CHR"].str.strip() + ":" + d["POS"].str.strip(), k37u)
    rs_by_vid = {v: d["rsid"] for v, d in usable.items() if isinstance(d["rsid"], str)}
    mn = scan(MANNING, ["Snp", "effect_allele", "other_allele", "MainEffects", "MainSE"], lambda d: d["Snp"],
              set(rs_by_vid.values()), sep=r"\s+")
    for u in uniq.itertuples(index=False):
        k = f"{u.chr}:{u.pos37}"
        for tr, df, ea, oa, b, s in [("tg_glgc", gl_, "ALT", "REF", "METAL_Effect", "METAL_StdErr"),
                                     ("bmi_yengo", y, "Tested_Allele", "Other_Allele", "BETA", "SE")]:
            r = df[(df["key"] == k) & pair_eq(df[ea].str.upper(), df[oa].str.upper(), u.a_raise, u.a_other)]
            if len(r) == 1:
                add(u.vid, tr, "direct", r.iloc[0][ea].upper(), r.iloc[0][oa].upper(), r.iloc[0][b], r.iloc[0][s])
        rs = rs_by_vid.get(u.vid)
        if rs:
            r = mn[(mn["key"] == rs) & pair_eq(mn["effect_allele"].str.upper(), mn["other_allele"].str.upper(), u.a_raise, u.a_other)]
            if len(r) == 1:
                add(u.vid, "homair_manning", "direct", r.iloc[0]["effect_allele"].upper(), r.iloc[0]["other_allele"].upper(),
                    r.iloc[0]["MainEffects"], r.iloc[0]["MainSE"])

    R = pd.DataFrame(rows)
    dup = R.duplicated(["id", "trait"], keep=False)
    if dup.any():
        log(f"dropping {dup.sum()} duplicated id x trait rows"); R = R[~dup]
    R.to_csv(OUT / "rows.tsv", sep="\t", index=False)
    inst = inst.merge(info.reset_index(), on="vid", how="left")
    inst["rsid"] = inst["vid"].map(rs_by_vid)
    inst.to_csv(OUT / "instruments.tsv", sep="\t", index=False)
    log(R.groupby(["trait", "source"])["id"].nunique().to_string())
    log("STEP1 DONE")


if __name__ == "__main__":
    main()
