#!/usr/bin/env python3
"""E1: reporter-to-endogenous bridge table (IMPLEMENTATION_SPEC.md section 6, E1 only).

Joins the 4,359 reconstructed GSE281364 MPRA paired SNVs to every endogenous liver
label on hg38 chromosome, position and alleles, trying both allele orders, recording
for every matched row which orientation matched and what sign convention the label
uses.

ORIENT ONCE: `label_value_oriented_to_mpra_alt` already carries the sign correction
(`orientation_multiplier`). Downstream consumers use that column as-is and must never
apply the multiplier a second time.

Outputs (in --outdir):
  bridge_variants.tsv          one row per (MPRA element x endogenous label)
  bridge_counts.tsv            per label coverage summary
  orientation_proof.tsv        Currin lead / LD-proxy beta vs nominal beta
  position_only_matches.tsv    same position, non-matching alleles (diagnostic)
  scan_stats.tsv               per-chromosome nominal-file scan statistics
  predictions.tsv              E1.1-E1.3 verdicts
  MANIFEST.tsv, pip_freeze.txt, run_config.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

PROJ = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SEED = 20260914  # recorded per program rule 9; E1 itself has no stochastic step.
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
AMBIG = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

P_ELEMENTS = f"{PROJ}/Analysis/MASLD_Model_Benchmark/executions/gse281364-outcome-blind-splits-21068885/split/elements.tsv"
P_GROUPMAP = f"{PROJ}/Analysis/MASLD_Model_Benchmark/executions/gse281364-borzoi-native-fixture-21083008/fixture/source_group_map.tsv"
P_NOMINAL_DIR = f"{PROJ}/data/external/currin_2025_caqtl"
P_LEADS = f"{PROJ}/GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz"
P_PROXY = f"{PROJ}/GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed.gz"
P_ASE_1367 = f"{PROJ}/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/allelic_sites.tsv"
P_ASE_4832 = f"{PROJ}/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/gse244832_allelic_sites.tsv"
P_BROADAWAY = f"{PROJ}/GWAS/finemapping/results/seqfunc/broadaway_benchmark_truth.tsv"
P_SQTL = f"{PROJ}/data/external/gtex_v8_liver_sqtl/GTEx_Analysis_v8_sQTL/Liver.v8.sqtl_signifpairs.txt.gz"
P_QUERYSET = f"{PROJ}/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables/query_set.tsv"

NOMINAL_TMPL = (
    "caQTL_hg38_{c}_with-varBeta-and-MAF_FastQTL-WASP-normal-nominal-"
    "updatedConsensusPeaks-25ATAC-PCs-minMAC10-results_1mbFromPeakCenters.txt.gz"
)
CHROMS = [f"chr{i}" for i in range(1, 23)]

AWK_PROG = r"""
BEGIN { FS="\t"; OFS="\t"; nrow=0; nuniq=0; nbad=0 }
NR==FNR {
    if ($1=="M") mp[$2]=1;
    else if ($1=="L") lv[$2]=1;
    else if ($1=="P") pv[$2]=1;
    next
}
FNR==1 { next }
{
    nrow++
    n=split($1, v, ":")
    if (n!=4) { nbad++; next }
    if (!($1 in seen)) { seen[$1]=1; nuniq++ }
    k = v[1] ":" v[2]
    if (k in mp) print "M", $0
    vp = $1 "|" $2
    if (vp in lv) print "L", $0
    if (vp in pv) print "P", $0
}
END { print "S", CHROM, nrow, nuniq, nbad }
"""

BRIDGE_COLS = [
    "element_id", "contig", "variant_pos1", "mpra_ref", "mpra_alt", "strand_ambiguous",
    "outer_locus_sequence_group_id", "outer_fold", "block_1mb", "block_lr239",
    "label", "label_source_path", "source_variant_id", "source_ref", "source_alt",
    "match_orientation", "orientation_multiplier",
    "label_value_source_scale", "label_value_oriented_to_mpra_alt", "label_se", "label_p",
    "label_unit", "label_sign_convention_source", "label_sign_convention_oriented",
    "label_feature", "label_n_tests", "label_extra",
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] {msg}", flush=True)


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def as_int_str(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").astype("Int64").astype(str)


def annotate_orientation(df: pd.DataFrame, sref: str, salt: str,
                         mref: str = "mpra_ref", malt: str = "mpra_alt") -> pd.DataFrame:
    """Add match_orientation and orientation_multiplier.

    direct              source (REF,ALT) == MPRA (ref,alt)              multiplier +1
    swapped             source (REF,ALT) == MPRA (alt,ref)              multiplier -1
    complement          complemented source matches directly            multiplier +1 (flagged)
    complement_swapped  complemented source matches swapped             multiplier -1 (flagged)
    no_match            alleles incompatible (includes indel alleles)   multiplier 0

    Priority is direct > swapped > complement > complement_swapped, so a strand-ambiguous
    SNV (A/T, C/G) is always resolved as a same-strand match and flagged separately.
    """
    df = df.copy()
    if len(df) == 0:
        df["match_orientation"] = pd.Series(dtype=object)
        df["orientation_multiplier"] = pd.Series(dtype=float)
        return df
    sr, sa = df[sref].astype(str), df[salt].astype(str)
    mr, ma = df[mref].astype(str), df[malt].astype(str)
    csr, csa = sr.map(COMP), sa.map(COMP)
    direct = (sr == mr) & (sa == ma)
    swapped = (sr == ma) & (sa == mr) & ~direct
    comp = (csr == mr) & (csa == ma) & ~direct & ~swapped
    comp_sw = (csr == ma) & (csa == mr) & ~direct & ~swapped & ~comp
    df["match_orientation"] = np.select(
        [direct, swapped, comp, comp_sw],
        ["direct", "swapped", "complement", "complement_swapped"], default="no_match")
    df["orientation_multiplier"] = np.select(
        [direct, swapped, comp, comp_sw], [1.0, -1.0, 1.0, -1.0], default=0.0)
    return df


def make_rows(el: pd.DataFrame, matched: pd.DataFrame, *, label: str, source_path: str,
              value_col: str, se_col, p_col, unit: str, convention: str,
              feature_col, n_tests_col, extra_col, svid_col: str,
              sref: str, salt: str) -> pd.DataFrame:
    """Assemble bridge rows. `value_col` must already be expressed per the SOURCE alt allele."""
    if len(matched) == 0:
        return pd.DataFrame(columns=BRIDGE_COLS)
    m = matched.merge(el, on="element_id", how="left", validate="many_to_one")
    out = pd.DataFrame({
        "element_id": m.element_id, "contig": m.contig, "variant_pos1": m.variant_pos1,
        "mpra_ref": m.mpra_ref, "mpra_alt": m.mpra_alt, "strand_ambiguous": m.strand_ambiguous,
        "outer_locus_sequence_group_id": m.outer_locus_sequence_group_id,
        "outer_fold": m.outer_fold, "block_1mb": m.block_1mb, "block_lr239": m.block_lr239,
        "label": label, "label_source_path": source_path,
        "source_variant_id": m[svid_col], "source_ref": m[sref], "source_alt": m[salt],
        "match_orientation": m.match_orientation,
        "orientation_multiplier": m.orientation_multiplier,
        "label_value_source_scale": pd.to_numeric(m[value_col], errors="coerce"),
        "label_value_oriented_to_mpra_alt": (pd.to_numeric(m[value_col], errors="coerce")
                                             * m.orientation_multiplier),
        "label_se": pd.to_numeric(m[se_col], errors="coerce") if se_col else np.nan,
        "label_p": pd.to_numeric(m[p_col], errors="coerce") if p_col else np.nan,
        "label_unit": unit, "label_sign_convention_source": convention,
        "label_sign_convention_oriented": "positive = higher signal for the MPRA genomic_alt allele",
        "label_feature": m[feature_col].astype(str) if feature_col else "",
        "label_n_tests": m[n_tests_col] if n_tests_col else 1,
        "label_extra": m[extra_col].astype(str) if extra_col else "",
    })
    return out[BRIDGE_COLS]


def scan_one(job):
    chrom, path, lookup, awkfile, outdir = job
    out = os.path.join(outdir, "scan", f"{chrom}.tsv")
    cmd = "set -o pipefail; zcat {f} | mawk -v CHROM={c} -f {a} {l} - > {o}".format(
        f=shlex.quote(path), c=chrom, a=shlex.quote(awkfile),
        l=shlex.quote(lookup), o=shlex.quote(out))
    t0 = time.time()
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    return chrom, r.returncode, round(time.time() - t0, 1), r.stderr[-2000:]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skip-scan", action="store_true",
                    help="reuse an existing scan/ directory (restart aid)")
    ap.add_argument("--skip-manifest", action="store_true",
                    help="smoke-test only: skip the input sha256 manifest")
    args = ap.parse_args()
    outdir = args.outdir
    os.makedirs(os.path.join(outdir, "scan"), exist_ok=True)
    log(f"outdir={outdir} seed={SEED} threads={args.threads}")

    # ---------------------------------------------------------------- MPRA side
    el = pd.read_csv(P_ELEMENTS, sep="\t")
    assert len(el) == 4359, len(el)
    el = el.rename(columns={"genomic_ref": "mpra_ref", "genomic_alt": "mpra_alt"})
    n_indel = int(((el.mpra_ref.str.len() != 1) | (el.mpra_alt.str.len() != 1)).sum())
    el = el[(el.mpra_ref.str.len() == 1) & (el.mpra_alt.str.len() == 1)].copy()
    el["strand_ambiguous"] = [(r, a) in AMBIG for r, a in zip(el.mpra_ref, el.mpra_alt)]
    el["block_1mb"] = el.contig + "~" + (el.variant_pos1 // 1_000_000).astype(str)
    gm = pd.read_csv(P_GROUPMAP, sep="\t")
    el = el.merge(gm[["source_locus_group_id", "borzoi_long_range_group_id"]].rename(
        columns={"source_locus_group_id": "outer_locus_sequence_group_id",
                 "borzoi_long_range_group_id": "block_lr239"}),
        on="outer_locus_sequence_group_id", how="left", validate="many_to_one")
    assert el.block_lr239.notna().all()
    el["key"] = el.contig + ":" + el.variant_pos1.astype(str)
    el = el[["element_id", "contig", "variant_pos1", "mpra_ref", "mpra_alt", "strand_ambiguous",
             "outer_locus_sequence_group_id", "outer_fold", "block_1mb", "block_lr239", "key"]]
    el_nokey = el.drop(columns=["key"])
    log(f"MPRA SNVs {len(el)} (indels normalized out: {n_indel}); "
        f"1-Mb blocks {el.block_1mb.nunique()}; 239-blocks {el.block_lr239.nunique()}; "
        f"strand-ambiguous {int(el.strand_ambiguous.sum())}")

    # -------------------------------------------- Currin lead + proxy (in memory)
    leads = pd.read_csv(P_LEADS, sep="\t")
    lv = leads.lead_variant_ID.str.split(":", expand=True)
    leads["lc"], leads["lp"], leads["lref"], leads["lalt"] = lv[0], lv[1], lv[2], lv[3]
    leads["key"] = leads.lc + ":" + leads.lp
    prox = pd.read_csv(P_PROXY, sep="\t")
    pv = prox.proxy_ID.str.split(":", expand=True)
    prox["pc"], prox["pp"], prox["pref"], prox["palt"] = pv[0], pv[1], pv[2], pv[3]
    prox["key"] = prox.pc + ":" + prox.pp

    # ------------------------------------------------- streaming scan of nominal
    lookup = os.path.join(outdir, "scan", "lookup.tsv")
    awkfile = os.path.join(outdir, "scan", "scan.awk")
    if not args.skip_scan:
        with open(lookup, "w") as fh:
            for k in el.key:
                fh.write(f"M\t{k}\n")
            for vid, pk in zip(leads.lead_variant_ID, leads.peak_ID):
                fh.write(f"L\t{vid}|{pk}\n")
            for vid, pk in zip(prox.proxy_ID, prox.caPeak):
                fh.write(f"P\t{vid}|{pk}\n")
        with open(awkfile, "w") as fh:
            fh.write(AWK_PROG)
        jobs = []
        for c in CHROMS:
            p = os.path.join(P_NOMINAL_DIR, NOMINAL_TMPL.format(c=c))
            assert os.path.exists(p), p
            jobs.append((c, p, lookup, awkfile, outdir))
        jobs.sort(key=lambda j: -os.path.getsize(j[1]))
        log(f"scanning {len(jobs)} Currin nominal files with {args.threads} workers")
        with Pool(args.threads) as pool:
            for chrom, rc, secs, err in pool.imap_unordered(scan_one, jobs):
                log(f"  scan {chrom} rc={rc} {secs}s {err.strip()[:200]}")
                if rc != 0:
                    raise RuntimeError(f"scan failed for {chrom}: {err}")

    ncols = ["variant", "peak", "distance_from_peakCenter", "pvalue", "beta", "varbeta",
             "imputation_R2", "MAF"]
    parts, stats = [], []
    for c in CHROMS:
        d = pd.read_csv(os.path.join(outdir, "scan", f"{c}.tsv"), sep="\t", header=None,
                        names=["tag"] + ncols, dtype=str, keep_default_na=False)
        s = d[d.tag == "S"].iloc[0]
        stats.append({"chrom": c, "rows_scanned": int(s["peak"]),
                      "distinct_variants": int(s["distance_from_peakCenter"]),
                      "malformed_variant_ids": int(s["pvalue"])})
        parts.append(d[d.tag != "S"])
    scan = pd.concat(parts, ignore_index=True)
    del parts
    st = pd.DataFrame(stats)
    st.to_csv(os.path.join(outdir, "scan_stats.tsv"), sep="\t", index=False)
    log(f"nominal rows scanned {st.rows_scanned.sum():,}; distinct variants "
        f"{st.distinct_variants.sum():,}; matched/verify rows kept {len(scan):,}")
    for col in ["distance_from_peakCenter", "pvalue", "beta", "varbeta", "imputation_R2", "MAF"]:
        scan[col] = pd.to_numeric(scan[col], errors="coerce")
    sv = scan.variant.str.split(":", expand=True)
    scan["sc"], scan["sp"], scan["sref"], scan["salt"] = sv[0], sv[1], sv[2], sv[3]
    scan["key"] = scan.sc + ":" + scan.sp

    bridge_parts, counts = [], []

    def add_count(label, rows, side_total, side_matched, side_unit, notes=""):
        counts.append({
            "label": label,
            "mpra_variants_matched": int(rows.element_id.nunique()) if len(rows) else 0,
            "blocks_1mb": int(rows.block_1mb.nunique()) if len(rows) else 0,
            "blocks_lr239": int(rows.block_lr239.nunique()) if len(rows) else 0,
            "n_direct": int((rows.match_orientation == "direct").sum()) if len(rows) else 0,
            "n_swapped": int((rows.match_orientation == "swapped").sum()) if len(rows) else 0,
            "n_complement": int(rows.match_orientation.isin(
                ["complement", "complement_swapped"]).sum()) if len(rows) else 0,
            "n_strand_ambiguous_matched": int(rows.strand_ambiguous.sum()) if len(rows) else 0,
            "endogenous_side_unit": side_unit,
            "endogenous_side_total": side_total,
            "endogenous_side_matched": side_matched,
            "endogenous_side_share": (round(side_matched / side_total, 8) if side_total else np.nan),
            "notes": notes,
        })

    # -------------------------------------------------- labels 1/2: Currin nominal
    nom = scan[scan.tag == "M"][["variant", "peak", "distance_from_peakCenter", "pvalue",
                                 "beta", "varbeta", "imputation_R2", "MAF", "sref", "salt",
                                 "key"]].copy()
    nom = nom.merge(el, on="key", how="inner")
    nom = annotate_orientation(nom, "sref", "salt")
    pos_only = nom[nom.match_orientation == "no_match"]
    pos_only[["element_id", "mpra_ref", "mpra_alt", "variant", "sref", "salt", "peak",
              "pvalue", "beta"]].to_csv(
        os.path.join(outdir, "position_only_matches.tsv"), sep="\t", index=False)
    log(f"nominal: {len(nom):,} variant-peak rows at MPRA positions; "
        f"position-only (allele mismatch) {len(pos_only):,} on "
        f"{pos_only.element_id.nunique()} elements")
    nom = nom[nom.match_orientation != "no_match"].copy()
    nom["label_se"] = np.sqrt(nom.varbeta)
    nom["abs_dist"] = nom.distance_from_peakCenter.abs()
    nom = nom.merge(nom.groupby("element_id").size().rename("label_n_tests").reset_index(),
                    on="element_id", how="left")
    nom["label_extra"] = ("imputation_R2=" + nom.imputation_R2.round(4).astype(str)
                          + ";MAF=" + nom.MAF.round(5).astype(str)
                          + ";distance_from_peakCenter=" + as_int_str(nom.distance_from_peakCenter)
                          + ";min_pvalue_over_peaks=" + nom.groupby("element_id").pvalue
                          .transform("min").map(lambda x: f"{x:.3e}"))
    conv = "beta is the effect of the source ALT allele (FastQTL nominal output)"
    keep = ["element_id", "variant", "sref", "salt", "match_orientation", "orientation_multiplier",
            "beta", "label_se", "pvalue", "peak", "label_n_tests", "label_extra"]
    minp = nom.sort_values(["element_id", "pvalue", "abs_dist"]).groupby(
        "element_id", as_index=False).first()
    near = nom.sort_values(["element_id", "abs_dist", "pvalue"]).groupby(
        "element_id", as_index=False).first()
    for lab, sel in (("currin_nominal_caqtl_minp_peak", minp),
                     ("currin_nominal_caqtl_nearest_peak", near)):
        rows = make_rows(el_nokey, sel[keep], label=lab,
                         source_path=f"{P_NOMINAL_DIR}/caQTL_hg38_chr*_...1mbFromPeakCenters.txt.gz",
                         value_col="beta", se_col="label_se", p_col="pvalue",
                         unit="caQTL beta, normalised ATAC peak signal per ALT allele copy",
                         convention=conv, feature_col="peak", n_tests_col="label_n_tests",
                         extra_col="label_extra", svid_col="variant", sref="sref", salt="salt")
        bridge_parts.append(rows)
        add_count(lab, rows, int(st.distinct_variants.sum()), int(sel.variant.nunique()),
                  "distinct variant IDs tested across the 22 nominal files",
                  "one peak per MPRA variant; selection rule named in the label. "
                  "min-p over peaks is a selected statistic, not a test")

    # ------------------------------------------------ orientation proof vs leads
    proof_rows = []
    lk = leads.copy()
    lk.index = lk.lead_variant_ID + "|" + lk.peak_ID
    lk = lk[~lk.index.duplicated()]
    lv_rows = scan[scan.tag == "L"].copy()
    lv_rows["vp"] = lv_rows.variant + "|" + lv_rows.peak
    lv_rows = lv_rows.drop_duplicates("vp")
    j = lv_rows.join(lk[["beta", "EA", "NEA", "allele_flip", "q_val"]].add_prefix("lead_"), on="vp")
    j = j[j.lead_beta.notna()].copy()
    j["ea_equals_ref"] = j.lead_EA == j.sref
    j["same_sign"] = np.sign(j.beta) == np.sign(j.lead_beta)
    j["abs_rel_diff"] = (j.beta - j.lead_beta).abs() / j.beta.abs().clip(lower=1e-12)
    for grp, sub in j.groupby("ea_equals_ref"):
        proof_rows.append({
            "comparison": "currin_lead_file_beta_vs_nominal_file_beta",
            "stratum": "EA==REF (allele_flip=yes)" if grp else "EA==ALT (allele_flip=no)",
            "n_pairs": len(sub),
            "n_beta_equal_rel_1e-4": int((sub.abs_rel_diff < 1e-4).sum()),
            "n_sign_agree": int(sub.same_sign.sum()),
            "pearson_beta": round(float(sub.beta.corr(sub.lead_beta)), 6) if len(sub) > 2 else np.nan,
        })
    px_rows = scan[scan.tag == "P"].copy()
    px_rows["vp"] = px_rows.variant + "|" + px_rows.peak
    px_rows = px_rows.drop_duplicates("vp")
    pk = prox.copy()
    pk.index = pk.proxy_ID + "|" + pk.caPeak
    pk = pk[~pk.index.duplicated()]
    jp = px_rows.join(pk[["proxy_caQTL_beta"]], on="vp")
    jp = jp[jp.proxy_caQTL_beta.notna()].copy()
    jp["abs_rel_diff"] = (jp.beta - jp.proxy_caQTL_beta).abs() / jp.beta.abs().clip(lower=1e-12)
    proof_rows.append({
        "comparison": "currin_LD_proxy_file_beta_vs_nominal_file_beta", "stratum": "all",
        "n_pairs": len(jp),
        "n_beta_equal_rel_1e-4": int((jp.abs_rel_diff < 1e-4).sum()),
        "n_sign_agree": int((np.sign(jp.beta) == np.sign(jp.proxy_caQTL_beta)).sum()),
        "pearson_beta": round(float(jp.beta.corr(jp.proxy_caQTL_beta)), 6) if len(jp) > 2 else np.nan,
    })
    pd.DataFrame(proof_rows).to_csv(os.path.join(outdir, "orientation_proof.tsv"),
                                    sep="\t", index=False)
    ecols = ["variant", "peak", "sref", "salt", "beta", "lead_beta", "lead_EA", "lead_NEA",
             "lead_allele_flip"]
    pd.concat([j[j.ea_equals_ref].head(4)[ecols], j[~j.ea_equals_ref].head(4)[ecols]]).to_csv(
        os.path.join(outdir, "orientation_proof_examples.tsv"), sep="\t", index=False)
    log("orientation proof:\n" + pd.DataFrame(proof_rows).to_string(index=False))

    # ---------------------------------------------------- label 3: Currin leads
    ld = leads.merge(el, on="key", how="inner")
    ld = annotate_orientation(ld, "lref", "lalt")
    ld = ld[ld.match_orientation != "no_match"].copy()
    ld = ld.sort_values(["element_id", "q_val"]).groupby("element_id", as_index=False).first()
    ld["label_extra"] = ("EA=" + ld.EA + ";NEA=" + ld.NEA + ";allele_flip=" + ld.allele_flip
                         + ";EUR_MAF=" + ld.EUR_MAF.astype(str))
    rows = make_rows(el_nokey,
                     ld[["element_id", "lead_variant_ID", "lref", "lalt", "match_orientation",
                         "orientation_multiplier", "beta", "q_val", "peak_ID", "label_extra"]],
                     label="currin_lead_caqtl", source_path=P_LEADS, value_col="beta",
                     se_col=None, p_col="q_val",
                     unit="caQTL beta, normalised ATAC peak signal per ALT allele copy",
                     convention=("beta follows the nominal-file ALT convention and is NOT "
                                 "re-oriented to the EA column; see orientation_proof.tsv"),
                     feature_col="peak_ID", n_tests_col=None, extra_col="label_extra",
                     svid_col="lead_variant_ID", sref="lref", salt="lalt")
    bridge_parts.append(rows)
    add_count("currin_lead_caqtl", rows, int(leads.lead_variant_ID.nunique()),
              int(ld.lead_variant_ID.nunique()), "distinct q<0.05 peak-lead variants",
              "LD-selected marginal associations, one lead per peak")

    # --------------------------------------------------- label 4: Currin proxies
    px = prox.merge(el, on="key", how="inner")
    px = annotate_orientation(px, "pref", "palt")
    px = px[px.match_orientation != "no_match"].copy()
    px = px.sort_values(["element_id", "proxy_ID"]).groupby("element_id", as_index=False).first()
    px["label_extra"] = "lead_ID=" + px.lead_ID + ";proxy_rsID=" + px.proxy_rsID.astype(str)
    rows = make_rows(el_nokey,
                     px[["element_id", "proxy_ID", "pref", "palt", "match_orientation",
                         "orientation_multiplier", "proxy_caQTL_beta", "caPeak", "label_extra"]],
                     label="currin_ld_proxy_caqtl", source_path=P_PROXY,
                     value_col="proxy_caQTL_beta", se_col=None, p_col=None,
                     unit="caQTL beta, normalised ATAC peak signal per ALT allele copy",
                     convention="nominal-file ALT convention (see orientation_proof.tsv)",
                     feature_col="caPeak", n_tests_col=None, extra_col="label_extra",
                     svid_col="proxy_ID", sref="pref", salt="palt")
    bridge_parts.append(rows)
    add_count("currin_ld_proxy_caqtl", rows, int(prox.proxy_ID.nunique()),
              int(px.proxy_ID.nunique()), "distinct in-peak LD proxies (r2>=0.8 with a lead)",
              "proxies of q<0.05 leads; not independent tests")

    # ------------------------------------------- labels 5/6: allelic imbalance
    for lab, path in (("allelic_imbalance_gse281367", P_ASE_1367),
                      ("allelic_imbalance_gse244832", P_ASE_4832)):
        a = pd.read_csv(path, sep="\t")
        a["key"] = a.chrom + ":" + a.pos.astype(str)
        m = a.merge(el, on="key", how="inner")
        m = annotate_orientation(m, "ref", "alt")
        m = m[m.match_orientation != "no_match"].drop_duplicates("element_id").copy()
        n_sites_matched = int(m.uid.nunique()) if len(m) else 0
        if len(m):
            m["label_extra"] = ("n_het_donors=" + m.n_het_donors.astype(str)
                                + ";total_reads=" + m.total_reads.astype(str)
                                + ";sd_log2=" + m.sd_log2.round(4).astype(str))
            m = m[["element_id", "uid", "ref", "alt", "match_orientation",
                   "orientation_multiplier", "mean_log2_alt_over_ref", "se_log2", "label_extra"]]
        rows = make_rows(el_nokey, m, label=lab, source_path=path,
                         value_col="mean_log2_alt_over_ref", se_col="se_log2", p_col=None,
                         unit="mean log2(ALT reads / REF reads) over heterozygous donors",
                         convention="positive = more reads carrying the source ALT allele",
                         feature_col=None, n_tests_col=None, extra_col="label_extra",
                         svid_col="uid", sref="ref", salt="alt") if len(m) else \
            pd.DataFrame(columns=BRIDGE_COLS)
        bridge_parts.append(rows)
        add_count(lab, rows, int(len(a)), n_sites_matched,
                  "measured allelic sites in the deposit",
                  "read-based het calls, no WASP or N-masked mapping-bias control")

    # ------------------------------------------------- label 7: Broadaway eQTL
    bw = pd.read_csv(P_BROADAWAY, sep="\t")
    bw["chrom"] = "chr" + bw.chr.astype(str)
    bw["key"] = bw.chrom + ":" + as_int_str(bw.pos_hg38)
    bw["beta_alt"] = np.where(bw.effect_allele == bw.alt, bw.eqtl_beta,
                              np.where(bw.effect_allele == bw.ref, -bw.eqtl_beta, np.nan))
    log(f"broadaway: effect_allele==alt in {int((bw.effect_allele == bw.alt).sum())} of "
        f"{len(bw)} rows; ==ref in {int((bw.effect_allele == bw.ref).sum())}")
    m = bw.merge(el, on="key", how="inner", suffixes=("_src", ""))
    m = annotate_orientation(m, "ref", "alt")
    m = m[(m.match_orientation != "no_match") & m.beta_alt.notna()].copy()
    m = m.merge(m.groupby("element_id").gene.nunique().rename("label_n_tests").reset_index(),
                on="element_id", how="left")
    m["label_extra"] = ("study=" + m.study.astype(str) + ";trait=" + m.trait.astype(str)
                        + ";pp4_best=" + pd.to_numeric(m.pp4_best, errors="coerce")
                        .round(4).astype(str))
    m = m.sort_values(["element_id", "eqtl_p"]).groupby("element_id", as_index=False).first()
    rows = make_rows(el_nokey,
                     m[["element_id", "variant_id", "ref", "alt", "match_orientation",
                        "orientation_multiplier", "beta_alt", "eqtl_se", "eqtl_p", "gene",
                        "label_n_tests", "label_extra"]],
                     label="broadaway_liver_eqtl", source_path=P_BROADAWAY, value_col="beta_alt",
                     se_col="eqtl_se", p_col="eqtl_p",
                     unit="liver eQTL beta (Broadaway meta-analysis) per ALT allele copy",
                     convention="eqtl_beta re-expressed per ALT using the file's effect_allele column",
                     feature_col="gene", n_tests_col="label_n_tests", extra_col="label_extra",
                     svid_col="variant_id", sref="ref", salt="alt")
    bridge_parts.append(rows)
    add_count("broadaway_liver_eqtl", rows,
              int(bw.drop_duplicates(["chrom", "pos_hg38", "ref", "alt"]).shape[0]),
              int(m.drop_duplicates(["chrom", "pos_hg38"]).shape[0]),
              "distinct hg38 variant-allele rows in the benchmark truth table",
              "colocalised liver eQTL rows only, not a genome-wide eQTL scan")

    # ---- audit: is each source's (REF,ALT) reference-oriented? Compare every source
    # against two independent reference-oriented panels (GTEx v8 b38 ids, Atlas query set)
    # at shared positions. A source that is reference-oriented matches "direct"; one that
    # stores effect/other alleles instead shows "swapped".
    audit = []
    qs_a = pd.read_csv(P_QUERYSET, sep="\t")
    qs_a["pkey"] = qs_a.hg38_chrom + ":" + qs_a.hg38_position_1based.astype(str)
    qs_a = qs_a.drop_duplicates("pkey")[["pkey", "hg38_ref", "hg38_alt"]]
    sq_a = pd.read_csv(P_SQTL, sep="\t", usecols=["variant_id"])
    sq_s = sq_a.variant_id.str.split("_", expand=True)
    sq_a = pd.DataFrame({"pkey": sq_s[0] + ":" + sq_s[1], "gref": sq_s[2], "galt": sq_s[3]})
    sq_a = sq_a.drop_duplicates("pkey")
    panels = {"atlas_query_set": (qs_a, "hg38_ref", "hg38_alt"),
              "gtex_v8_liver_sqtl": (sq_a, "gref", "galt")}
    sources = {
        "broadaway_liver_eqtl": (pd.DataFrame({
            "pkey": "chr" + bw.chr.astype(str) + ":" + as_int_str(bw.pos_hg38),
            "sr": bw.ref, "sa": bw.alt}).drop_duplicates("pkey")),
        "currin_lead_caqtl": (pd.DataFrame({
            "pkey": leads.key, "sr": leads.lref, "sa": leads.lalt}).drop_duplicates("pkey")),
        "allelic_imbalance_gse281367": None, "mpra_elements": None,
    }
    a1 = pd.read_csv(P_ASE_1367, sep="\t")
    sources["allelic_imbalance_gse281367"] = pd.DataFrame({
        "pkey": a1.chrom + ":" + a1.pos.astype(str), "sr": a1.ref,
        "sa": a1.alt}).drop_duplicates("pkey")
    sources["mpra_elements"] = pd.DataFrame({
        "pkey": el.key, "sr": el.mpra_ref, "sa": el.mpra_alt}).drop_duplicates("pkey")
    for sname, sdf in sources.items():
        for pname, (pdf, pr, pa) in panels.items():
            mm = sdf.merge(pdf, on="pkey", how="inner")
            if len(mm) == 0:
                continue
            mm = annotate_orientation(mm, "sr", "sa", pr, pa)
            audit.append({"source": sname, "reference_oriented_panel": pname, "n_shared": len(mm),
                          "n_direct": int((mm.match_orientation == "direct").sum()),
                          "n_swapped": int((mm.match_orientation == "swapped").sum()),
                          "n_other": int((~mm.match_orientation.isin(["direct", "swapped"])).sum())})
    pd.DataFrame(audit).to_csv(os.path.join(outdir, "source_orientation_audit.tsv"),
                               sep="\t", index=False)
    log("source orientation audit:\n" + pd.DataFrame(audit).to_string(index=False))
    del qs_a, sq_a, sources

    # -------------------------------------------------- label 8: GTEx liver sQTL
    sq = pd.read_csv(P_SQTL, sep="\t")
    sv2 = sq.variant_id.str.split("_", expand=True)
    sq["sc2"], sq["sp2"], sq["sref"], sq["salt"] = sv2[0], sv2[1], sv2[2], sv2[3]
    sq["key"] = sq.sc2 + ":" + sq.sp2
    m = sq.merge(el, on="key", how="inner")
    m = annotate_orientation(m, "sref", "salt")
    m = m[m.match_orientation != "no_match"].copy()
    m = m.merge(m.groupby("element_id").phenotype_id.nunique()
                .rename("label_n_tests").reset_index(), on="element_id", how="left")
    m["label_extra"] = "maf=" + m.maf.astype(str) + ";tss_distance=" + m.tss_distance.astype(str)
    m = m.sort_values(["element_id", "pval_nominal"]).groupby("element_id", as_index=False).first()
    rows = make_rows(el_nokey,
                     m[["element_id", "variant_id", "sref", "salt", "match_orientation",
                        "orientation_multiplier", "slope", "slope_se", "pval_nominal",
                        "phenotype_id", "label_n_tests", "label_extra"]],
                     label="gtex_v8_liver_sqtl", source_path=P_SQTL, value_col="slope",
                     se_col="slope_se", p_col="pval_nominal",
                     unit="GTEx v8 sQTL slope (intron excision ratio) per ALT allele copy",
                     convention="GTEx v8 slope is the effect of the ALT allele named in variant_id",
                     feature_col="phenotype_id", n_tests_col="label_n_tests",
                     extra_col="label_extra", svid_col="variant_id", sref="sref", salt="salt")
    bridge_parts.append(rows)
    add_count("gtex_v8_liver_sqtl", rows, int(sq.variant_id.nunique()),
              int(m.variant_id.nunique()), "distinct variants in Liver.v8.sqtl_signifpairs",
              "significant pairs only; coverage is conditional on sQTL significance")

    # ---------------------------------------------- label 9: Atlas query set
    qs = pd.read_csv(P_QUERYSET, sep="\t")
    qs["key"] = qs.hg38_chrom + ":" + qs.hg38_position_1based.astype(str)
    m = qs.merge(el, on="key", how="inner")
    m = annotate_orientation(m, "hg38_ref", "hg38_alt")
    m = m[m.match_orientation != "no_match"].drop_duplicates("element_id").copy()
    m["served"] = 1.0
    m["label_extra"] = ""
    rows = make_rows(el_nokey,
                     m[["element_id", "variant_uid", "hg38_ref", "hg38_alt", "match_orientation",
                        "orientation_multiplier", "served", "label_extra"]],
                     label="atlas_query_set_membership", source_path=P_QUERYSET,
                     value_col="served", se_col=None, p_col=None,
                     unit="membership indicator (1 = present in the Track 0 query set)",
                     convention="membership only, unsigned; orientation recorded for identity",
                     feature_col=None, n_tests_col=None, extra_col="label_extra",
                     svid_col="variant_uid", sref="hg38_ref", salt="hg38_alt")
    rows["label_value_oriented_to_mpra_alt"] = 1.0  # unsigned: never flipped
    bridge_parts.append(rows)
    add_count("atlas_query_set_membership", rows, int(qs.variant_uid.nunique()),
              int(m.variant_uid.nunique()), "Track 0 Atlas query-set variants",
              "membership flag only; oriented value fixed to 1 because membership is unsigned")

    # -------------------------------------------------------------- write output
    bridge = pd.concat(bridge_parts, ignore_index=True)
    bridge = bridge.sort_values(["label", "contig", "variant_pos1"])
    bridge.to_csv(os.path.join(outdir, "bridge_variants.tsv"), sep="\t", index=False)
    cdf = pd.DataFrame(counts)
    ai = bridge[bridge.label.str.startswith("allelic_imbalance")]
    cdf.loc[len(cdf)] = {
        "label": "allelic_imbalance_union", "mpra_variants_matched": int(ai.element_id.nunique()),
        "blocks_1mb": int(ai.block_1mb.nunique()), "blocks_lr239": int(ai.block_lr239.nunique()),
        "n_direct": int((ai.match_orientation == "direct").sum()),
        "n_swapped": int((ai.match_orientation == "swapped").sum()),
        "n_complement": int(ai.match_orientation.isin(["complement", "complement_swapped"]).sum()),
        "n_strand_ambiguous_matched": int(ai.strand_ambiguous.sum()),
        "endogenous_side_unit": "union of GSE281367 (435) and GSE244832 (279) site rows",
        "endogenous_side_total": 714,
        "endogenous_side_matched": int(ai.source_variant_id.nunique()),
        "endogenous_side_share": round(ai.source_variant_id.nunique() / 714, 8),
        "notes": "the two deposits share 214 uids; the 714 denominator is site rows, not unique sites",
    }
    cdf.to_csv(os.path.join(outdir, "bridge_counts.tsv"), sep="\t", index=False)
    log("bridge_counts:\n" + cdf.to_string(index=False))

    def matched(lbl):
        v = cdf.loc[cdf.label == lbl, "mpra_variants_matched"]
        return int(v.iloc[0]) if len(v) else 0

    n_nom, n_ai, n_bw = matched("currin_nominal_caqtl_minp_peak"), matched(
        "allelic_imbalance_union"), matched("broadaway_liver_eqtl")
    preds = pd.DataFrame([
        {"prediction": "E1.1",
         "statement": ">=1,000 of the 4,359 MPRA variants carry a Currin nominal caQTL test",
         "predicted": ">=1000", "observed": n_nom, "met": bool(n_nom >= 1000)},
        {"prediction": "E1.2", "statement": "<=20 carry a measured allelic-imbalance value",
         "predicted": "<=20", "observed": n_ai, "met": bool(n_ai <= 20)},
        {"prediction": "E1.3", "statement": ">=100 carry a Broadaway liver eQTL row",
         "predicted": ">=100", "observed": n_bw, "met": bool(n_bw >= 100)},
    ])
    preds.to_csv(os.path.join(outdir, "predictions.tsv"), sep="\t", index=False)
    log("predictions:\n" + preds.to_string(index=False))

    # block-count table for the E2 feasibility note
    blocks = (bridge.groupby("label")
              .agg(n_variants=("element_id", "nunique"),
                   blocks_1mb=("block_1mb", "nunique"),
                   blocks_lr239=("block_lr239", "nunique"),
                   median_abs_value=("label_value_oriented_to_mpra_alt",
                                     lambda s: float(np.nanmedian(np.abs(s))) if len(s) else np.nan))
              .reset_index())
    blocks.to_csv(os.path.join(outdir, "e2_block_feasibility.tsv"), sep="\t", index=False)
    log("block feasibility:\n" + blocks.to_string(index=False))

    # manifest + environment
    man = []
    for p in ([P_ELEMENTS, P_GROUPMAP, P_LEADS, P_PROXY, P_ASE_1367, P_ASE_4832,
               P_BROADAWAY, P_SQTL, P_QUERYSET]
              + [os.path.join(P_NOMINAL_DIR, NOMINAL_TMPL.format(c=c)) for c in CHROMS]):
        man.append({"path": p, "bytes": os.path.getsize(p),
                    "mtime_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                               time.gmtime(os.path.getmtime(p))),
                    "sha256": "skipped" if args.skip_manifest else sha256(p)})
    pd.DataFrame(man).to_csv(os.path.join(outdir, "MANIFEST.tsv"), sep="\t", index=False)
    with open(os.path.join(outdir, "pip_freeze.txt"), "w") as fh:
        fh.write(subprocess.run([sys.executable, "-m", "pip", "freeze"],
                                capture_output=True, text=True).stdout)
    with open(os.path.join(outdir, "run_config.json"), "w") as fh:
        json.dump({"seed": SEED, "stochastic_steps": "none; E1 is deterministic counting",
                   "python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
                   "threads": args.threads, "mpra_indels_removed": n_indel,
                   "mpra_snvs": int(len(el)),
                   "nominal_rows_scanned": int(st.rows_scanned.sum()),
                   "nominal_distinct_variants": int(st.distinct_variants.sum())}, fh, indent=2)
    log("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
