#!/usr/bin/env python3
"""Ancestry step 3: tag REF/ALT heterozygote frequency and tag-lead LD in each 1kGP superpopulation.

Inputs
  --ref-vcf      1kGP high-coverage phased hardcalls, founders only, biallelic SNV records at tag and
                 lead positions (anc1 <out_vcf>).
  --tagspan-vcf  every panel record (SNV, indel, SV, split multiallelic) whose REF span covers a tag
                 position (anc1 <out_tagspan_vcf>).
  --tags         tags_final.tsv (gene, GTEx lead_variant_id, tag chrom/pos/ref/alt, r_eur, r_eas, orientation).
  --ped          1kGP ped file.
  --unrelated    kgp_unrelated.keep from anc2 (founders minus KING >= 0.0884 relatives).
  --exclude-bed  MHC / IG / TR spans (0-based BED); tags inside are flagged (in_mhc_ig_tr_span), not dropped.
Outputs (in --out, all outcome-free population summaries)
  tag_2pq_by_superpop.tsv   one row per tag SNV. Per group (AFR, AMR, EAS, EUR, SAS and JPT; unrelated
                            founders): n_hap, alt_freq (tag record), p_ref, p_alt, other_alt_freq,
                            deleted_freq, conflict_freq, het_ref_alt = 2 * p_ref * p_alt (the column the
                            T2 het rule should use), obs_het_ref_alt (observed REF/ALT individuals).
                            multiallelic_in_1kgp and max_other_allele_freq (other base or deletion, max
                            over the five superpopulations).
  tag_ld_by_superpop.tsv    one row per tags_final row: signed r and r2 to the lead per superpopulation,
                            sign agreement with the tag orientation, and flags.
  tags_multiallelic_in_1kgp.tsv  tag rows whose other-base-or-deletion frequency is > 0.01 in any
                            superpopulation.
  tag_ld_summary.json       counts, and the re-derivation check of r_eur / r_eas against tags_final.
Rules
  * Haplotype state at the tag base: start REF; a carried record that puts the tag ALT base there -> ALT;
    another base -> OTHER; a deletion of the base -> DEL; two different non-REF states on one haplotype
    -> CONFLICT. Record effect at offset o = tag_pos - POS: base ALT[o] if o < len(ALT), else deleted.
    p_ref, p_alt, other, deleted and conflict frequencies sum to 1.
  * r is the Pearson correlation of phased haplotype alleles (ALT = 1), i.e. plink2 --r2-phased,
    never an EM estimate from unphased genotypes. Records with any unphased call are skipped (as T1).
  * flag_r2_lt_0.8_afr_amr_sas: r2 < 0.8 in AFR, AMR or SAS. low_r2_EUR_larger_ref /
    low_r2_EAS_larger_ref / low_r2_eur_eas_larger_ref: r2 < 0.8 in EUR / EAS / either, using this
    step's unrelated founders (larger than the T1 sample). r is undefined when the tag or lead is
    monomorphic in that superpopulation; that is flagged separately (undefined_<POP>), not as low r2.
  * Check: r_eur and r_eas are recomputed on the T1 sample definition (founders that are not
    parents) and compared with tags_final.tsv.
"""
import argparse
import bisect
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
EXTRA_POPS = ["JPT"]  # the three Japanese cohorts project to EAS and sit nearest JPT on PC1-10
FLAG_POPS = ["AFR", "AMR", "SAS"]
MIN_R2 = 0.8
OTHER_ALLELE_REPORT = 0.01
REF, ALT, OTHER, DEL, CONFLICT = 0, 1, 2, 3, 4


def query(vcf):
    """-> (samples, list of (chrom, pos, ref, alt, haps n_samples x 2 int8), skipped counts)."""
    samples = subprocess.run(["bcftools", "query", "-l", vcf], capture_output=True, text=True, check=True).stdout.split()
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", vcf]
    recs, skipped = [], {"unphased_or_missing": 0, "multiallelic_record": 0}
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc:
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if "," in f[3]:
                skipped["multiallelic_record"] += 1
                continue
            gts = f[4:]
            if any(len(g) != 3 or g[1] != "|" or g[0] not in "01" or g[2] not in "01" for g in gts):
                skipped["unphased_or_missing"] += 1
                continue
            h = np.array([[g[0], g[2]] for g in gts], dtype=np.int8)
            recs.append((f[0], int(f[1]), f[2], f[3], h))
    return samples, recs, skipped


def freq(h):
    return float(h.mean()) if h.size else np.nan


def corr(t, l):
    t, l = t.reshape(-1).astype(float), l.reshape(-1).astype(float)
    if t.std() == 0 or l.std() == 0:
        return np.nan
    return float(np.corrcoef(t, l)[0, 1])


def effect(rec, tag):
    """State a carried record puts on the tag base, or None (base unchanged); 'mismatch' if REF disagrees."""
    _, rpos, rref, ralt, _ = rec
    _, tpos, tref, talt = tag
    o = tpos - rpos
    if o < 0 or o >= len(rref):
        return None
    if rref[o] != tref:
        return "mismatch"
    if o < len(ralt):
        b = ralt[o]
        return None if b == tref else (ALT if b == talt else OTHER)
    return DEL


def main():
    ap = argparse.ArgumentParser()
    for k in ["ref-vcf", "tagspan-vcf", "tags", "ped", "unrelated", "exclude-bed", "out"]:
        ap.add_argument(f"--{k}", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    ped = pd.read_csv(a.ped, sep=r"\s+", dtype=str, keep_default_na=False)
    parents = set(ped["FatherID"]) | set(ped["MotherID"])
    unrel = set(pd.read_csv(a.unrelated, sep="\t", dtype=str)["#IID"])
    sp = ped.set_index("SampleID")["Superpopulation"]
    popl = ped.set_index("SampleID")["Population"]
    founders_not_parents = set(ped[(ped["FatherID"] == "0") & (ped["MotherID"] == "0")]["SampleID"]) - parents

    samples, recs, skipped = query(a.ref_vcf)
    haps = {(c, p, r, al): h for c, p, r, al, h in recs}
    idx = {s: i for i, s in enumerate(samples)}
    groups = {p: np.array([idx[s] for s in samples if s in unrel and sp[s] == p]) for p in SUPERPOPS}
    groups.update({p: np.array([idx[s] for s in samples if s in unrel and popl[s] == p]) for p in EXTRA_POPS})
    t1_groups = {p: np.array([idx[s] for s in samples if s in founders_not_parents and sp[s] == p]) for p in ["EUR", "EAS"]}

    span_samples, span_recs, span_skipped = query(a.tagspan_vcf)
    span_recs = list({(c, p, r, al): (c, p, r, al, h) for c, p, r, al, h in span_recs}.values())  # dedup
    sidx = {s: i for i, s in enumerate(span_samples)}
    span_groups = {p: np.array([sidx[s] for s in span_samples if s in unrel and (sp[s] == p or popl[s] == p)])
                   for p in SUPERPOPS + EXTRA_POPS}
    assert all(len(span_groups[p]) == len(groups[p]) for p in groups), "tag-span VCF sample set differs"

    tags = pd.read_csv(a.tags, sep="\t", dtype={"pos": int}, keep_default_na=False)
    lead = tags["lead_variant_id"].str.split("_", expand=True)
    tags["lead_key"] = list(zip(lead[0], lead[1].astype(int), lead[2], lead[3]))
    tags["tag_key"] = list(zip(tags["chrom"], tags["pos"], tags["ref"], tags["alt"]))

    # records covering each tag base
    by_chrom = {}
    for r in span_recs:
        by_chrom.setdefault(r[0], []).append(r)
    tag_pos_sorted = {c: sorted({k[1] for k in tags["tag_key"] if k[0] == c}) for c in tags["chrom"].unique()}
    covering = {}
    for c, rs in by_chrom.items():
        tp = tag_pos_sorted.get(c, [])
        for r in rs:
            lo = bisect.bisect_left(tp, r[1])
            hi = bisect.bisect_right(tp, r[1] + len(r[2]) - 1)
            for p in tp[lo:hi]:
                covering.setdefault((c, p), []).append(r)

    # REF/ALT heterozygote frequency per unique tag SNV
    rows, mismatches, no_tag_record = [], 0, 0
    n_span = len(span_samples)
    for key in dict.fromkeys(tags["tag_key"]):
        row = dict(zip(["chrom", "pos", "ref", "alt"], key))
        state = np.zeros((n_span, 2), np.int8)
        found_tag, n_rec = False, 0
        for r in covering.get((key[0], key[1]), []):
            e = effect(r, key)
            if e == "mismatch":
                mismatches += 1
                continue
            if e is None:
                continue
            n_rec += 1
            found_tag |= (r[1], r[2], r[3]) == key[1:]
            carr = r[4] == 1
            clash = carr & (state != REF) & (state != e)
            state = np.where(carr, e, state).astype(np.int8)
            state[clash] = CONFLICT
        no_tag_record += not found_tag
        row["in_panel"] = found_tag
        row["records_changing_tag_base"] = n_rec
        h_tag = haps.get(key)
        other_max = 0.0
        for p in SUPERPOPS + EXTRA_POPS:
            s = state[span_groups[p]] if found_tag else np.empty((0, 2), np.int8)  # no tag record: all NaN
            fr = {c: float((s == c).mean()) if s.size else np.nan for c in (REF, ALT, OTHER, DEL, CONFLICT)}
            obs = float((np.sort(s, axis=1) == [REF, ALT]).all(axis=1).mean()) if s.size else np.nan
            row[f"n_hap_{p}"] = s.size
            row[f"alt_freq_{p}"] = freq(h_tag[groups[p]]) if h_tag is not None else np.nan
            row[f"p_ref_{p}"], row[f"p_alt_{p}"] = fr[REF], fr[ALT]
            row[f"other_alt_freq_{p}"], row[f"deleted_freq_{p}"], row[f"conflict_freq_{p}"] = fr[OTHER], fr[DEL], fr[CONFLICT]
            row[f"het_ref_alt_{p}"] = 2 * fr[REF] * fr[ALT]
            row[f"obs_het_ref_alt_{p}"] = obs
            if p in SUPERPOPS and found_tag:
                other_max = max(other_max, fr[OTHER] + fr[DEL] + fr[CONFLICT])
        row["max_other_allele_freq"] = other_max
        row["multiallelic_in_1kgp"] = other_max > 0
        rows.append(row)
    tpq = pd.DataFrame(rows)
    bed = pd.read_csv(a.exclude_bed, sep="\t", header=None, dtype={0: str}).iloc[:, :3]
    in_span = np.zeros(len(tpq), bool)
    for _, b in bed.iterrows():
        in_span |= ((tpq["chrom"] == b[0]) & (tpq["pos"] - 1 >= b[1]) & (tpq["pos"] - 1 < b[2])).to_numpy()
    tpq["in_mhc_ig_tr_span"] = in_span
    tpq.to_csv(out / "tag_2pq_by_superpop.tsv", sep="\t", index=False, float_format="%.6g")

    # tag-lead r per superpopulation
    rows = []
    for _, t in tags.iterrows():
        h_t, h_l = haps.get(t["tag_key"]), haps.get(t["lead_key"])
        r = {c: t[c] for c in ["gene_id", "gene_name", "lead_variant_id", "chrom", "pos", "ref", "alt",
                               "orientation", "tag_is_lead", "r_eur", "r_eas"]}
        r["in_panel"] = h_t is not None and h_l is not None
        for p in SUPERPOPS:
            rr = corr(h_t[groups[p]], h_l[groups[p]]) if r["in_panel"] else np.nan
            r[f"r_{p}"], r[f"r2_{p}"] = rr, rr ** 2
            r[f"undefined_{p}"] = bool(np.isnan(rr))
            r[f"sign_flip_{p}"] = bool(not np.isnan(rr) and np.sign(rr) != int(t["orientation"]))
        for p, col in [("EUR", "r_eur"), ("EAS", "r_eas")]:
            r[f"r_{p}_t1_samples"] = corr(h_t[t1_groups[p]], h_l[t1_groups[p]]) if r["in_panel"] else np.nan
        r["flag_r2_lt_0.8_afr_amr_sas"] = bool(any(r[f"r2_{p}"] < MIN_R2 for p in FLAG_POPS if not np.isnan(r[f"r2_{p}"])))
        r["low_r2_pops"] = ",".join(p for p in FLAG_POPS if not np.isnan(r[f"r2_{p}"]) and r[f"r2_{p}"] < MIN_R2)
        r["undefined_pops"] = ",".join(p for p in FLAG_POPS if r[f"undefined_{p}"])
        r["sign_flip_pops"] = ",".join(p for p in SUPERPOPS if r[f"sign_flip_{p}"])
        r["low_r2_EUR_larger_ref"] = bool(r["r2_EUR"] < MIN_R2)
        r["low_r2_EAS_larger_ref"] = bool(r["r2_EAS"] < MIN_R2)
        r["low_r2_eur_eas_larger_ref"] = r["low_r2_EUR_larger_ref"] or r["low_r2_EAS_larger_ref"]
        r["r2_min_EUR_EAS"] = np.nanmin([r["r2_EUR"], r["r2_EAS"]]) if r["in_panel"] else np.nan
        rows.append(r)
    ld = pd.DataFrame(rows)
    ld.to_csv(out / "tag_ld_by_superpop.tsv", sep="\t", index=False, float_format="%.6g")

    # tags with a common other base or deletion at the tag position
    other_cols = [f"{k}_{p}" for p in SUPERPOPS for k in ("other_alt_freq", "deleted_freq", "p_ref", "alt_freq", "het_ref_alt")]
    multi = tpq[tpq["max_other_allele_freq"] > OTHER_ALLELE_REPORT]
    genes = tags.groupby("tag_key")["gene_name"].agg(lambda g: ",".join(sorted(set(g)))).to_dict()
    multi = multi.assign(gene_names=[genes[k] for k in zip(multi["chrom"], multi["pos"], multi["ref"], multi["alt"])])
    multi[["chrom", "pos", "ref", "alt", "gene_names", "max_other_allele_freq"] + other_cols].to_csv(
        out / "tags_multiallelic_in_1kgp.tsv", sep="\t", index=False, float_format="%.4g")

    def low_block(col, r2col):
        f = ld[ld[col]]
        per_gene = ld.groupby("gene_id")[col].all()
        return {"rows": int(len(f)), "genes": int(f["gene_id"].nunique()), "genes_all_tags_affected": int(per_gene.sum()),
                "r2_range": [round(float(f[r2col].min()), 4), round(float(f[r2col].max()), 4)] if len(f) else None}

    d_eur = (ld["r_EUR_t1_samples"] - ld["r_eur"].astype(float)).abs()
    d_eas = (ld["r_EAS_t1_samples"] - ld["r_eas"].astype(float)).abs()
    ok_all = (ld[[f"r2_{p}" for p in SUPERPOPS]] >= MIN_R2).all(axis=1).groupby(ld["gene_id"]).any()
    old_2pq = {p: 2 * tpq[f"alt_freq_{p}"] * (1 - tpq[f"alt_freq_{p}"]) for p in SUPERPOPS}
    info = {
        "tags_rows": int(len(tags)), "tag_snvs": int(len(tpq)), "genes": int(tags["gene_id"].nunique()),
        "panel_records_skipped": skipped, "tagspan_records_skipped": span_skipped,
        "tagspan_records_used": int(len(span_recs)), "tagspan_ref_base_mismatches": mismatches,
        "tag_snvs_without_own_record_in_tagspan": no_tag_record,
        "tag_rows_missing_tag_or_lead_in_panel": int((~ld["in_panel"]).sum()),
        "unrelated_founders_by_group": {p: int(len(g)) for p, g in groups.items()},
        "t1_founders_not_parents": {p: int(len(g)) for p, g in t1_groups.items()},
        "recheck_max_abs_diff_r_eur_vs_tags_final": float(d_eur.max()),
        "recheck_max_abs_diff_r_eas_vs_tags_final": float(d_eas.max()),
        "tag_rows_flag_r2_lt_0.8_afr_amr_sas": int(ld["flag_r2_lt_0.8_afr_amr_sas"].sum()),
        "tag_rows_r2_lt_0.8_by_pop": {p: int((ld[f"r2_{p}"] < MIN_R2).sum()) for p in SUPERPOPS},
        "low_r2_larger_ref": {"EUR": low_block("low_r2_EUR_larger_ref", "r2_EUR"),
                              "EAS": low_block("low_r2_EAS_larger_ref", "r2_EAS"),
                              "EUR_or_EAS": low_block("low_r2_eur_eas_larger_ref", "r2_min_EUR_EAS")},
        "tag_rows_undefined_by_pop": {p: int(ld[f"undefined_{p}"].sum()) for p in SUPERPOPS},
        "tag_rows_sign_flip_by_pop": {p: int(ld[f"sign_flip_{p}"].sum()) for p in SUPERPOPS},
        "genes_with_a_tag_r2_ge_0.8_in_all_five": int(ok_all.sum()),
        "het_column_for_T2": "het_ref_alt_<POP> = 2 * p_ref * p_alt (REF and ALT at the tag base; other bases and deletions excluded)",
        "tag_snvs_in_mhc_ig_tr_span": int(tpq["in_mhc_ig_tr_span"].sum()),
        "tag_snvs_in_mhc_ig_tr_span_list": [f"{c}:{q}:{r}:{al}" for c, q, r, al in
                                            tpq.loc[tpq["in_mhc_ig_tr_span"], ["chrom", "pos", "ref", "alt"]].itertuples(index=False)],
        "tag_snvs_multiallelic_in_1kgp_any": int(tpq["multiallelic_in_1kgp"].sum()),
        f"tag_snvs_other_allele_freq_gt_{OTHER_ALLELE_REPORT}_any_superpop": int(len(multi)),
        f"tag_snvs_other_allele_freq_gt_{OTHER_ALLELE_REPORT}_by_superpop": {
            p: int((tpq[f"other_alt_freq_{p}"] + tpq[f"deleted_freq_{p}"] + tpq[f"conflict_freq_{p}"] > OTHER_ALLELE_REPORT).sum())
            for p in SUPERPOPS},
        "tag_snvs_het_ref_alt_differs_from_alt_only_2pq_by_gt_0.01": {
            p: int((np.abs(tpq[f"het_ref_alt_{p}"] - old_2pq[p]) > 0.01).sum()) for p in SUPERPOPS},
        "median_het_ref_alt_by_group": {p: round(float(tpq[f"het_ref_alt_{p}"].median()), 4) for p in SUPERPOPS + EXTRA_POPS},
    }
    (out / "tag_ld_summary.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
