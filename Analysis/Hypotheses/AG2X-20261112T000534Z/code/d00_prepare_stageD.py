#!/usr/bin/env python3
"""Stage D preparation (no hosted call): locus selection by the prespecified rule, motif scan of REF and ALT,
measured data from files on disk, and the hosted-call plan.

Selection rule (PRESPEC_AG2X_stageD.md section 1):
  tier 1 loci of stageB_tier1_locus_table.tsv; a credible-set member that is sequence-supported in Stage B,
  noncoding, with max PIP >= 0.5. Excluded if a coding member (class coding_protein_altering) with PIP >= 0.5 is
  in a credible set where the candidate itself has PIP >= 0.5 ("qualifying set"). Every set containing the
  candidate is also checked and reported. Amendment (orchestrator, 2026-11-12, before any Stage D call): also
  excluded if a protein-altering member (class coding_protein_altering) with PIP >= 0.5 in any credible set of
  the same locus lies within 50 kb. At most 8 loci; SORT1 rs12740374 added as the positive control.
Noncoding: Stage B class 'noncoding'; for members with no class (lifted or TRAIL4-top members), no overlap
  with a GENCODE v44 CDS (cellranger-arc GRCh38-2024-A genes.gtf.gz), labelled 'noncoding_by_gtf'.
Motif scan: the ATACGEN A3 PWMs (JASPAR2024 CORE vertebrates, 879 motifs, A3 log2 PWMs re-referenced to the A3
  background) and A3 p < 5e-5 thresholds; score = max log2 score over windows that cover the variant and lie
  inside the 64-bp window, both strands. Gained/lost: max(s_ref, s_alt) >= threshold and |delta| >= 2 (A3
  disruption rule); strong if |delta| >= 3. TF expression in hepatocytes: gene passes edgeR filterByExpr in
  GSE296875 donor hepatocyte pseudobulk (ATACGEN A4 posthoc_tested_genes_hepatocyte.tsv); A3 itself uses no
  expression source.
Every output row carries candidate_not_adopted = True.
"""
from __future__ import annotations

import gzip
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

B = L.PKG / "stageB"
D = L.PKG / "stageD"
IN = D / "inputs"
REPO = L.REPO
GTF = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/genes/genes.gtf.gz"
A3 = REPO / "Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption"
A4_GENES = REPO / "Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation/work/posthoc_tested_genes_hepatocyte.tsv"
TRAIL4_SETS = REPO / "Analysis/Hypotheses/TRAIL4-20261111T025501Z/out_23190150/susie_credible_sets.tsv"
FM = REPO / "GWAS/finemapping/results"
CURRIN = REPO / "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1"
D1 = FM / "alphagenome_program/d1-lineage-allelic-20260914T194213Z/tables/lineage_site_summary.tsv"
MPRA = FM / "seqfunc/mpra_benchmark/mpra_substrate_truth.tsv"
COLOC = FM / "susie_coloc_r2_20260924T053244Z"
SORT1 = {"uid": "chr1:109274968:G:T", "locus_id": "DIG1_6", "label": "SORT1 (control)", "hg19": "1:109817590"}
MAX_LOCI = 8
HALF = 32  # 64-bp window: positions pos-31 .. pos+32 (1-based); the variant is at 0-based offset 31
CALLS = {"tracks_conditions": 3, "ism_ref_conditions": 2, "ism_ref_bp": 64, "ism_alt_conditions": 1, "ism_alt_bp": 32}
CAP_D = 6000


def flip(u):
    c, p, r, a = u.split(":")
    return f"{c}:{p}:{a}:{r}"


def cds_intervals(chroms):
    out = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            if f[2] == "CDS" and f[0] in chroms:
                out.setdefault(f[0], []).append((int(f[3]), int(f[4])))
    return out


def select():
    t1 = pd.read_csv(B / "tables/stageB_tier1_locus_table.tsv", sep="\t")
    t1_loci = set(t1.locus_id)
    vf = pd.read_csv(B / "tables/stageB_variant_fields.tsv.gz", sep="\t",
                     usecols=["uid", "roles", "loci", "sequence_supported", "supported_conditions",
                              "accessibility_sign_alt_minus_ref", "coding_pip05"])
    mv = pd.read_csv(B / "inputs/member_variants.tsv", sep="\t", dtype=str)
    mv["max_pip"] = mv.max_pip.astype(float)
    imap = dict(pd.read_csv(B / "inputs/indel_uid_map.tsv", sep="\t").values)
    mv["uid"] = mv.uid.map(lambda u: imap.get(u, u))
    cs = pd.read_csv(B / "inputs/credible_set_member_rows.tsv.gz", sep="\t", dtype=str)
    cs["pip"] = cs.pip.astype(float)
    cs["uid"] = cs.uid.map(lambda u: imap.get(u, u))
    cls = dict(zip(mv.uid, mv["class"].fillna("")))
    cs["class"] = cs.uid.map(cls).fillna("")

    sup = vf[(vf.sequence_supported == True) & vf.roles.str.contains(r"(^|;)member($|;)")]  # noqa: E712
    cand = sup.merge(mv[["uid", "k", "max_pip", "class", "Consequence", "primary_locus", "traits", "studies",
                         "n_sets", "member_sources"]], on="uid", how="left")
    cand = cand[cand.primary_locus.isin(t1_loci)].copy()
    chroms = set(cand.uid.str.split(":").str[0])
    cds = cds_intervals(chroms)

    def nc_status(r):
        if r["class"] == "noncoding":
            return "noncoding"
        if isinstance(r["class"], str) and r["class"] and r["class"] != "nan":
            return r["class"]
        c, p = r.uid.split(":")[:2]
        p = int(p)
        hit = any(s <= p <= e for s, e in cds.get(c, []))
        return "coding_by_gtf_CDS_overlap" if hit else "noncoding_by_gtf"
    cand["class"] = cand["class"].fillna("")
    cand["noncoding_status"] = cand.apply(nc_status, axis=1)
    all_sup_t1 = cand.copy()
    cand = cand[(cand.max_pip >= 0.5) & cand.noncoding_status.isin(["noncoding", "noncoding_by_gtf"])].copy()

    tr = pd.read_csv(TRAIL4_SETS, sep="\t", dtype=str)
    tr["set_id"] = "TRAIL4|" + tr.test + "|" + tr.cs
    tr_size = dict(zip(tr.set_id, tr["size"]))

    rows, set_rows = [], []
    for _, r in cand.iterrows():
        mine = cs[cs.uid == r.uid]
        best = mine.sort_values("pip", ascending=False).iloc[0]
        qual = set(mine.loc[mine.pip >= 0.5, "set_id"])
        allsets = set(mine.set_id)
        members = cs[cs.set_id.isin(allsets)]
        cod = members[(members["class"] == "coding_protein_altering") & (members.pip >= 0.5)]
        cod_q = cod[cod.set_id.isin(qual)]
        cod_any = cod[cod.set_id.isin(allsets)]
        locus_cod = cs[(cs.locus_id == r.primary_locus) & (cs["class"] == "coding_protein_altering") & (cs.pip >= 0.5)]
        pos = int(r.uid.split(":")[1])
        lc = (locus_cod.groupby("uid").pip.max().reset_index()
              .assign(dist_kb=lambda x: (x.uid.str.split(":").str[1].astype(int) - pos).abs() / 1000))
        trail_sets = sorted(s for s in allsets if s.startswith("TRAIL4|"))
        for s in sorted(allsets):
            m = members[members.set_id == s]
            set_rows.append({"uid": r.uid, "set_id": s, "candidate_pip": float(mine.loc[mine.set_id == s, "pip"].max()),
                             "qualifying_set": s in qual, "n_members_saved": int(m.uid.nunique()),
                             "trail4_set_size": tr_size.get(s, ""),
                             "coding_members_pip05": ";".join(sorted(set(cod.loc[cod.set_id == s, "uid"]))),
                             "candidate_not_adopted": True})
        rows.append({
            "locus_id": r.primary_locus, "uid": r.uid, "hg19_key": r.k, "max_pip": r.max_pip,
            "max_pip_set": best.set_id, "max_pip_trait": best.trait, "max_pip_study": best.study,
            "max_pip_ancestry": best.ancestry, "max_pip_method": best.method, "member_source": best.member_source,
            "n_sets": len(allsets), "n_qualifying_sets": len(qual), "noncoding_status": r.noncoding_status,
            "consequence": r.Consequence, "supported_conditions": r.supported_conditions,
            "stageB_sign": r.accessibility_sign_alt_minus_ref,
            "coding_pip05_in_qualifying_set": ";".join(sorted(set(cod_q.uid))),
            "coding_pip05_in_any_set_of_candidate": ";".join(sorted(set(cod_any.uid))),
            "coding_pip05_elsewhere_in_locus": ";".join(f"{u} (PIP {p:.3f}; {d:.1f} kb)" for u, p, d in lc.values),
            "trail4_sets_top_only": ";".join(f"{s} (size {tr_size.get(s, 'NA')})" for s in trail_sets),
            "excluded_coding_in_qualifying_set": bool(len(cod_q)),
            "excluded_protein_altering_within_50kb": bool((lc.dist_kb <= 50).any()),
            "excluded": bool(len(cod_q)) or bool((lc.dist_kb <= 50).any()), "candidate_not_adopted": True})
    sel = pd.DataFrame(rows)
    ex_coding = all_sup_t1[(all_sup_t1.max_pip >= 0.5) & ~all_sup_t1.noncoding_status.isin(["noncoding", "noncoding_by_gtf"])]
    kept = sel[~sel.excluded].copy()
    loci_order = (kept.groupby("locus_id").max_pip.max().sort_values(ascending=False).index.tolist())
    if len(loci_order) > MAX_LOCI:  # prespecified: keep the 8 loci with the highest max PIP
        loci_order = loci_order[:MAX_LOCI]
    kept = kept[kept.locus_id.isin(loci_order)]
    kept = kept.sort_values(["locus_id", "uid"])
    ctrl = pd.DataFrame([{"locus_id": SORT1["locus_id"], "uid": SORT1["uid"], "hg19_key": SORT1["hg19"] + ":G:T",
                          "role": "control_SORT1", "candidate_not_adopted": True}])
    kept["role"] = "selected"
    final = pd.concat([ctrl, kept], ignore_index=True)
    # label for panel folder
    t1lab = dict(zip(t1.locus_id, t1.genes.fillna("")))
    final["genes_stageB"] = final.locus_id.map(t1lab).fillna("SORT1, CELSR2, PSRC1, PSMA5, SARS1")
    bad = L.assert_ref(final.uid.tolist())
    notes = {"tier1_loci": len(t1_loci), "tier1_sequence_supported_members": int(len(all_sup_t1)),
             "maxpip05_noncoding_candidates": int(len(sel)), "excluded_by_coding_set_rule": sel.loc[sel.excluded_coding_in_qualifying_set, "uid"].tolist(),
             "excluded_by_protein_altering_50kb_rule": sel.loc[sel.excluded_protein_altering_within_50kb, "uid"].tolist(),
             "maxpip05_not_noncoding": ex_coding[["uid", "primary_locus", "max_pip", "noncoding_status"]].to_dict("records"),
             "loci_selected": loci_order, "variants_selected": kept.uid.tolist(), "ref_mismatch": bad,
             "all_tier1_supported": all_sup_t1[["uid", "primary_locus", "max_pip", "noncoding_status"]].to_dict("records")}
    return final, sel, pd.DataFrame(set_rows), notes


# ------------------------------------------------------------------ motifs
def load_motifs():
    pw = pd.read_csv(A3 / "work/motif_pwm_log2.tsv", sep="\t")
    meta = pd.read_csv(A3 / "work/motif_meta_thresholds.tsv", sep="\t")
    clus = pd.read_csv(A3 / "results/motif_clusters.tsv", sep="\t", usecols=["motif_id", "cluster_id", "cluster_label",
                                                                             "cluster_family"])
    meta = meta.merge(clus, on="motif_id", how="left")
    mats = {m: g.sort_values("pos")[list("ACGT")].to_numpy(float) for m, g in pw.groupby("motif_id", sort=False)}
    return meta, mats


def onehot(s):
    x = np.zeros((len(s), 4))
    for i, ch in enumerate(s):
        j = "ACGT".find(ch)
        if j >= 0:
            x[i, j] = 1
    return x


def window_scores(seq, pwm):
    """All windows inside seq, both strands: returns list of (start, end_exclusive, strand, score)."""
    X = onehot(seq)
    Lm = pwm.shape[0]
    rc = pwm[::-1, ::-1]
    out = []
    for s in range(0, len(seq) - Lm + 1):
        w = X[s:s + Lm]
        out.append((s, s + Lm, "+", float((w * pwm).sum())))
        out.append((s, s + Lm, "-", float((w * rc).sum())))
    return out


def motif_scan(final):
    import pyfaidx
    fa = pyfaidx.Fasta(str(L.FASTA), as_raw=True, sequence_always_upper=True)
    meta, mats = load_motifs()
    genes = pd.read_csv(A4_GENES, sep="\t", usecols=["symbol"]).symbol.str.upper()
    expressed = set(genes)
    var_rows, hit_rows = [], []
    for u in final.uid:
        c, p, r, a = L.parse_uid(u)
        s0 = p - HALF  # 0-based start; variant at 0-based offset HALF - 1 = 31
        ref = fa[c][s0:s0 + 2 * HALF]
        off = HALF - 1
        assert ref[off] == r, (u, ref[off], r)
        alt = ref[:off] + a + ref[off + 1:]
        for _, m in meta.iterrows():
            pwm = mats[m.motif_id]
            t = m["threshold_p5e-5"]
            best = {}
            for allele, seq in (("ref", ref), ("alt", alt)):
                ws = window_scores(seq, pwm)
                cov = [w for w in ws if w[0] <= off < w[1]]
                best[allele] = max(w[3] for w in cov) if cov else np.nan
                for w in ws:
                    if w[3] >= t:
                        hit_rows.append({"uid": u, "allele": allele, "motif_id": m.motif_id, "name": m["name"],
                                         "cluster_id": m.cluster_id, "cluster_family": m.cluster_family,
                                         "start_rel": w[0] - off, "end_rel": w[1] - 1 - off, "strand": w[2],
                                         "score": w[3], "threshold": t, "covers_variant": w[0] <= off < w[1]})
            if np.isnan(best["ref"]):
                continue
            d = best["alt"] - best["ref"]
            if max(best.values()) >= t and abs(d) >= 2:
                tfs = [x.upper() for x in str(m["name"]).split("::")]
                var_rows.append({"uid": u, "motif_id": m.motif_id, "name": m["name"], "family": m.family,
                                 "cluster_id": m.cluster_id, "cluster_label": m.cluster_label,
                                 "cluster_family": m.cluster_family, "s_ref": best["ref"], "s_alt": best["alt"],
                                 "delta_alt_minus_ref": d, "threshold": t, "direction": "gained" if d > 0 else "lost",
                                 "strong_abs_delta_ge3": abs(d) >= 3,
                                 "tf_expressed_hepatocyte": ";".join(f"{x}:{'yes' if x in expressed else 'no'}" for x in tfs),
                                 "candidate_not_adopted": True})
    var = pd.DataFrame(var_rows)
    hits = pd.DataFrame(hit_rows)
    # one row per (variant, motif cluster): the member with the largest |delta|
    if len(var):
        cl = (var.assign(ad=var.delta_alt_minus_ref.abs()).sort_values("ad", ascending=False)
              .drop_duplicates(["uid", "cluster_id"]).drop(columns="ad"))
        nmem = var.groupby(["uid", "cluster_id"]).size().rename("n_member_motifs_disrupted").reset_index()
        cl = cl.merge(nmem, on=["uid", "cluster_id"])
    else:
        cl = var
    seqs = []
    for u in final.uid:
        c, p, r, a = L.parse_uid(u)
        ref = fa[c][p - HALF:p + HALF]
        seqs.append({"uid": u, "start0": p - HALF, "end0": p + HALF, "variant_offset0": HALF - 1, "ref_seq": ref,
                     "alt_seq": ref[:HALF - 1] + a + ref[HALF:]})
    return var, cl, hits, pd.DataFrame(seqs)


# ------------------------------------------------------------------ measured data
def measured(final):
    uids = set(final.uid) | {flip(u) for u in final.uid}
    out = []
    lead = pd.read_csv(CURRIN / "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz",
                       sep="\t", usecols=["lead_variant_ID", "peak_ID", "beta", "p_val_nominal", "q_val", "EA", "NEA",
                                          "distance_between_leadVariant_peakCenter"])
    prox = pd.read_csv(CURRIN / "caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed.gz", sep="\t",
                       usecols=["proxy_ID", "lead_ID", "caPeak", "proxy_caQTL_beta"])
    peaks = pd.read_csv(CURRIN / "supplementalData1_liver_ATAC_peaks.bed.gz", sep="\t", header=None, comment="#",
                        usecols=[0, 1, 2], names=["chr", "start", "end"])
    d1 = pd.read_csv(D1, sep="\t")
    d1 = d1[d1.lineage == "Hepatocyte"]
    mp = pd.read_csv(MPRA, sep="\t", dtype=str)
    mp["uid"] = mp.chr + ":" + mp.pos_hg38 + ":" + mp.ref_hg38.fillna("") + ":" + mp.alt_hg38.fillna("")
    sp = pd.read_csv(COLOC / "susie_coloc_signal_pairs.tsv.gz", sep="\t", dtype={"hit1": str})
    sp = sp[sp["PP.H4"] >= 0.5].copy()
    sp["c19"] = sp.hit1.str.split(":").str[0]
    sp["p19"] = sp.hit1.str.split(":").str[1].astype(int)
    ab = pd.read_csv(COLOC / "susie_coloc_all_gwas.csv", usecols=["gwas_name", "gene", "chr", "PP.H4.abf", "top_snp"],
                     dtype={"top_snp": str})
    ab = ab[ab["PP.H4.abf"] >= 0.5].dropna(subset=["top_snp"]).copy()
    ab["c19"] = ab.top_snp.str.split(":").str[0]
    ab["p19"] = pd.to_numeric(ab.top_snp.str.split(":").str[1], errors="coerce")
    for _, r in final.iterrows():
        u = r.uid
        c, p = u.split(":")[:2]
        p = int(p)
        rec = {"uid": u, "locus_id": r.locus_id}
        # caQTL (raw lead beta and proxy beta = hg38 ALT effect; Stage B section 3)
        eff = []
        for _, x in lead[lead.lead_variant_ID.isin([u, flip(u)])].iterrows():
            s = 1 if x.lead_variant_ID == u else -1
            eff.append(f"lead {x.peak_ID} beta_alt {s * x.beta:.3f} (p {x.p_val_nominal:.2g}, q {x.q_val:.2g})")
        for _, x in prox[prox.proxy_ID.isin([u, flip(u)]) & (prox.proxy_ID != prox.lead_ID)].iterrows():
            s = 1 if x.proxy_ID == u else -1
            eff.append(f"proxy of {x.lead_ID} {x.caPeak} beta_alt {s * x.proxy_caQTL_beta:.3f}")
        rec["currin_caqtl"] = "; ".join(eff) if eff else "not measured (no caQTL lead or r2>0.8 proxy row)"
        rec["currin_caqtl_beta_alt_values"] = ";".join(
            [f"{(1 if x.lead_variant_ID == u else -1) * x.beta:.4f}" for _, x in lead[lead.lead_variant_ID.isin([u, flip(u)])].iterrows()] +
            [f"{(1 if x.proxy_ID == u else -1) * x.proxy_caQTL_beta:.4f}" for _, x in
             prox[prox.proxy_ID.isin([u, flip(u)]) & (prox.proxy_ID != prox.lead_ID)].iterrows()])
        pk = peaks[(peaks.chr == c) & (peaks.start < p) & (peaks.end >= p)]
        rec["in_currin_liver_atac_peak"] = bool(len(pk))
        # hepatocyte allelic imbalance (D1; POOLED_30_DONORS is the Stage B arm; per-cohort rows reported)
        ai = []
        for _, x in d1[d1.uid.isin([u, flip(u)])].iterrows():
            s = 1 if x.uid == u else -1
            ai.append(f"{x.scope}: log2(ALT/REF) {s * x.mean_log2_alt_over_ref:.3f} (SE {x.se_log2:.3f}, "
                      f"{int(x.n_het_donors)} het donors, z {s * x.mean_log2_alt_over_ref / x.se_log2:.2f})")
        pooled = d1[d1.uid.isin([u, flip(u)]) & (d1.scope == "POOLED_30_DONORS")]
        rec["hep_ai"] = "; ".join(ai) if ai else "not measured (no heterozygous-donor row)"
        if len(pooled):
            x = pooled.iloc[0]
            s = 1 if x.uid == u else -1
            rec["hep_ai_pooled_log2_alt_ref"] = s * x.mean_log2_alt_over_ref
            rec["hep_ai_pooled_z"] = s * x.mean_log2_alt_over_ref / x.se_log2
            rec["hep_ai_pooled_n_donors"] = int(x.n_het_donors)
        m = mp[mp.uid.isin([u, flip(u)])]
        rec["mpra"] = ("; ".join(f"{x.context}: barcode_mapped {x.mpra_barcode_mapped}, DAV {x.dav}, log2FC {x.mpra_log2fc}"
                                 for _, x in m.iterrows()) if len(m) else "not measured (not in the MPRA substrate table)")
        # COLOC r2 (A3 Q4 rule): PP.H4 >= 0.5 with GWAS lead within 250 kb (hg19)
        c19, p19 = r.hg19_key.split(":")[:2]
        p19 = int(p19)
        s1 = sp[(sp.c19 == c19) & ((sp.p19 - p19).abs() <= 250_000)]
        g1 = (s1.groupby("gene").agg(pp=("PP.H4", "max"), n=("gwas", "nunique"),
                                     eligible=("eligible", lambda x: bool(x.astype(str).str.upper().eq("TRUE").any())))
              .sort_values("pp", ascending=False))
        rec["coloc_r2_susie_genes"] = ("; ".join(f"{g} (max PP.H4 {x.pp:.2f}, {x.n} studies, eligible {x.eligible})"
                                                 for g, x in g1.iterrows()) or "none within 250 kb")
        a1 = ab[(ab.c19 == c19) & ((ab.p19 - p19).abs() <= 250_000)]
        g2 = a1.groupby("gene")["PP.H4.abf"].agg(["max", "count"]).sort_values("max", ascending=False)
        rec["coloc_r2_abf_genes"] = ("; ".join(f"{g} (max PP.H4.abf {x['max']:.2f}, {int(x['count'])} studies)"
                                               for g, x in g2.iterrows()) or "none within 250 kb")
        rec["candidate_not_adopted"] = True
        out.append(rec)
    return pd.DataFrame(out)


def calls_plan(final):
    n = len(final)
    per = (CALLS["tracks_conditions"] + CALLS["ism_ref_conditions"] * CALLS["ism_ref_bp"] * 3
           + CALLS["ism_alt_conditions"] * CALLS["ism_alt_bp"] * 3)
    rows = [{"item": "tracks (predict_variant; one call per variant x condition; REF and ALT in one response)",
             "per_variant": CALLS["tracks_conditions"], "total": n * CALLS["tracks_conditions"]},
            {"item": "ISM on REF, 64 bp x 3 substitutions, bulk liver (ATAC, DNASE) and hepatocyte (SC_ATAC)",
             "per_variant": CALLS["ism_ref_conditions"] * CALLS["ism_ref_bp"] * 3,
             "total": n * CALLS["ism_ref_conditions"] * CALLS["ism_ref_bp"] * 3},
            {"item": "ISM on ALT background, central 32 bp x 3, hepatocyte (SC_ATAC)",
             "per_variant": CALLS["ism_alt_conditions"] * CALLS["ism_alt_bp"] * 3,
             "total": n * CALLS["ism_alt_conditions"] * CALLS["ism_alt_bp"] * 3},
            {"item": "total", "per_variant": per, "total": n * per}]
    t = pd.DataFrame(rows)
    t["variants"] = n
    t["stageD_cap"] = CAP_D
    return t


def main():
    IN.mkdir(parents=True, exist_ok=True)
    (D / "motifs").mkdir(exist_ok=True)
    (D / "measured").mkdir(exist_ok=True)
    final, sel, sets, notes = select()
    final.to_csv(IN / "stageD_variants.tsv", sep="\t", index=False)
    sel.to_csv(IN / "selection_candidates.tsv", sep="\t", index=False)
    sets.to_csv(IN / "selection_sets_checked.tsv", sep="\t", index=False)
    json.dump(notes, open(IN / "selection_notes.json", "w"), indent=1, default=str)
    L.log(json.dumps({k: v for k, v in notes.items() if k != "all_tier1_supported"}, default=str))
    var, cl, hits, seqs = motif_scan(final)
    var.to_csv(D / "motifs/motif_changes_all_motifs.tsv", sep="\t", index=False)
    cl.to_csv(D / "motifs/motif_changes_by_cluster.tsv", sep="\t", index=False)
    hits.to_csv(D / "motifs/motif_hits_64bp.tsv.gz", sep="\t", index=False)
    seqs.to_csv(D / "motifs/sequences_64bp.tsv", sep="\t", index=False)
    me = measured(final)
    me.to_csv(D / "measured/measured_at_variants.tsv", sep="\t", index=False)
    plan = calls_plan(final)
    plan.to_csv(IN / "calls_plan.tsv", sep="\t", index=False)
    L.log(plan.to_string())
    L.log("no hosted call made")


if __name__ == "__main__":
    main()
