#!/usr/bin/env python3
"""
69_borzoi_magnitude_nom.py  --  AXIS 4: Borzoi MAGNITUDE nomination + gene
attribution for eQTL-ABSENT fine-mapped liver signals (direction-INDEPENDENT).

ADDITIVE ONLY -- new file. It IMPORTS (never edits) 64b_borzoi_eqtl_logsed.py and
reuses its gene-model / windowing / 4-fold-ensemble / fwd+RC / GTEx-liver-matched
logSED machinery VERBATIM. It does NOT touch 46d/78/27a or any existing 60-65 script.
Output is a NOMINATION / mechanism-class table only -- it is NEVER wired into the
convergence atlas or a scored channel (that would be circular with COLOC + epigenomic
evidence).

WHY MAGNITUDE, NOT DIRECTION
  The eQTL-DIRECTION gate FAILED for BOTH Borzoi and AlphaGenome (sign-auROC ~0.53-
  0.56, permutation p>=0.08). So we make NO directional claim. But Borzoi |logSED|
  tracks |eQTL beta| (Spearman 0.390, p=3.6e-8 in the 64b benchmark) -- a legitimate
  MAGNITUDE nomination (the published Sniff use). And Borzoi does NOT beat TSS-distance
  at PICKING the eGene, so we NEVER let the model pick a gene alone: we score a
  proximity-defined candidate set and report the best-MAGNITUDE gene as a HYPOTHESIS to
  be paired downstream with PoPS / TSS-distance.

WHAT THIS DOES
  For each independent eQTL-absent signal lead (from 69a_borzoi_magnitude_prep.py):
    1. ENUMERATE candidate genes: every protein_coding / lncRNA gene whose TSS lies
       within Borzoi's input receptive field of the variant (|var - TSS| < 262,144 bp).
    2. SCORE each candidate's ensemble |logSED| MAGNITUDE with the 64b recipe --
       gene-TSS-centered 524,288 bp window, 4-fold ensemble (replicate-0..3),
       forward+reverse-complement averaged within fold, GTEx-liver tissue-matched
       tracks (7563/7564/7565), per-fold log(sum_exon ALT +1) - log(sum_exon REF +1),
       folds averaged.  IDENTICAL geometry to the reference distribution below, so the
       magnitudes are directly comparable.
    3. NOMINATE the best-magnitude candidate gene per locus + its 4-fold sign-agreement
       (0-4) + the magnitude PERCENTILE of |logSED| against the colocalizing / eQTL-
       PRESENT lead distribution (64b's borzoi_eqtl_logsed_scores.tsv), which calibrates
       what "large effect" means on the same scale.

  CAVEAT (reported, not hidden): best_gene is an argmax over K candidates, so its
  percentile carries a best-of-K upward bias vs the single-eGene reference. We ALSO
  emit the magnitude + percentile of the substrate's pre-specified (VEP) gene, an
  apples-to-apples single-gene comparison, and a full per-candidate long sidecar.

OUTPUTS (results/seqfunc/)
  borzoi_magnitude_nominations.tsv   -- one row per locus (headline)
  borzoi_magnitude_candidates.tsv    -- one row per (locus, candidate gene) [long]
"""
import argparse
import csv
import gzip
import importlib.util
import os
import sys
import time

import numpy as np

# ---------------------------------------------------------------------------
# Import 64b machinery VERBATIM (no edits). Module name starts with a digit, so
# load it by path with importlib.
# ---------------------------------------------------------------------------
_SRC = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "borzoi64b", os.path.join(_SRC, "64b_borzoi_eqtl_logsed.py"))
b64 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(b64)

# constants + helpers reused verbatim
SEQ_LEN = b64.SEQ_LEN
N_HUMAN_TRACKS = b64.N_HUMAN_TRACKS
N_OUT_BINS = b64.N_OUT_BINS
BIN_BP = b64.BIN_BP
OUT_HALF = b64.OUT_HALF
IN_HALF = b64.IN_HALF
GTEX_LIVER_TRACK_INDICES = b64.GTEX_LIVER_TRACK_INDICES
FOLDS = b64.FOLDS
one_hot = b64.one_hot
rc_one_hot = b64.rc_one_hot
fetch_window = b64.fetch_window
find_axes = b64.find_axes
load_liver_tracks = b64.load_liver_tracks
build_rc_track_perm = b64.build_rc_track_perm
_attr = b64._attr
strip_ver = b64.strip_ver
gene_tss = b64.gene_tss
exon_union_bins = b64.exon_union_bins
resolve_hg38_alleles = b64.resolve_hg38_alleles


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Genome-wide gene models (protein_coding + lncRNA) -- same record schema as
# 64b.parse_gtf_genes, but biotype-filtered genome-wide instead of name-filtered,
# so we can ENUMERATE candidates. Reuses 64b's _attr / strip_ver / gene_tss /
# exon_union_bins for everything downstream.
# ---------------------------------------------------------------------------
WANT_BIOTYPES = ("protein_coding", "lncRNA")


def parse_all_genes(gtf_path):
    genes = {}  # (sid, chrom) -> record
    opener = gzip.open if gtf_path.endswith(".gz") else open
    n_lines = 0
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line and line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat = f[2]
            if feat not in ("gene", "transcript", "exon"):
                continue
            chrom = f[0]
            # primary assembly only (Borzoi FASTA has chr1..22,X,Y; drop *_alt/_random/M)
            if "_" in chrom or chrom in ("chrM", "chrEBV"):
                continue
            attrs = f[8]
            bt = _attr(attrs, "gene_type") or _attr(attrs, "gene_biotype")
            if bt not in WANT_BIOTYPES:
                continue
            gid = _attr(attrs, "gene_id")
            sid = strip_ver(gid)
            gname = _attr(attrs, "gene_name")
            strand = f[6]
            start = int(f[3])
            end = int(f[4])
            key = (sid, chrom)
            rec = genes.get(key)
            if rec is None:
                rec = {"sid": sid, "chrom": chrom, "strand": strand,
                       "gene_start": None, "gene_end": None, "gene_name": gname,
                       "biotype": bt, "mane_tss": None, "exons": []}
                genes[key] = rec
            if feat == "gene":
                rec["gene_start"] = start
                rec["gene_end"] = end
                rec["strand"] = strand
                rec["gene_name"] = gname
                rec["biotype"] = bt
            elif feat == "transcript":
                if 'tag "MANE_Select"' in attrs:
                    rec["mane_tss"] = start if strand == "+" else end
            elif feat == "exon":
                rec["exons"].append((start, end))
    # keep only records with a gene line + >=1 exon; index TSS per contig
    by_chr = {}
    n_ok = 0
    for (sid, chrom), rec in genes.items():
        if rec["gene_start"] is None or not rec["exons"]:
            continue
        tss1, strand = gene_tss(rec)
        rec["tss1"] = tss1
        by_chr.setdefault(chrom, []).append(rec)
        n_ok += 1
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda r: r["tss1"])
    log(f"GTF: {n_ok} protein_coding/lncRNA gene models (primary contigs); "
        f"contigs={len(by_chr)}")
    return by_chr


def candidates_in_window(by_chr, chrom, v0):
    """protein_coding/lncRNA genes with TSS within IN_HALF (input receptive field)
    of the 0-based variant position, nearest-first."""
    arr = by_chr.get(chrom, [])
    lo = v0 - IN_HALF
    hi = v0 + IN_HALF
    import bisect
    tss = [r["tss1"] - 1 for r in arr]  # 0-based
    i0 = bisect.bisect_left(tss, lo)
    i1 = bisect.bisect_right(tss, hi)
    cand = [arr[i] for i in range(i0, i1) if abs((arr[i]["tss1"] - 1) - v0) < IN_HALF]
    cand.sort(key=lambda r: abs((r["tss1"] - 1) - v0))
    return cand


# ---------------------------------------------------------------------------
def load_reference(path):
    """|borzoi_logsed_liver| of the eQTL-PRESENT (colocalizing) lead distribution
    scored by 64b -- non-anchor, non-NaN. Sorted ascending for percentile lookup."""
    vals = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("is_anchor") == "TRUE":
                continue
            v = r.get("borzoi_logsed_liver", "")
            if v in ("", "nan", "NaN"):
                continue
            try:
                vals.append(abs(float(v)))
            except ValueError:
                continue
    vals.sort()
    return np.array(vals, dtype=np.float64)


def pct_of(ref_sorted, x):
    """percentile (0-100) of x within ref_sorted (fraction of ref <= x)."""
    if len(ref_sorted) == 0 or not np.isfinite(x):
        return float("nan")
    return 100.0 * float(np.searchsorted(ref_sorted, x, side="right")) / len(ref_sorted)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--liver-tracks", required=True)
    ap.add_argument("--leads", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--folds", nargs="+", default=list(FOLDS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-candidates", required=True)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # ---- leads -------------------------------------------------------------
    leads = []
    with open(args.leads) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            leads.append(r)
    log(f"loaded {len(leads)} eQTL-absent signal leads")

    # ---- reference (eQTL-present |logSED| distribution) --------------------
    ref = load_reference(args.reference)
    log(f"reference |logSED_liver| distribution: n={len(ref)}  "
        f"p50={np.percentile(ref,50):.4f} p80={np.percentile(ref,80):.4f} "
        f"p90={np.percentile(ref,90):.4f} max={ref.max():.4f}")
    ref_p80 = float(np.percentile(ref, 80))

    # ---- gene models -------------------------------------------------------
    by_chr = parse_all_genes(args.gtf)

    # ---- liver tracks (verbatim 64b manifest handling) ---------------------
    liver = load_liver_tracks(args.liver_tracks)
    liver_idx = np.array([t["track_index"] for t in liver], dtype=int)
    pos_in_subset = {t["track_index"]: i for i, t in enumerate(liver)}
    rna_tracks = [t for t in liver if t["assay"] == "rna_liver"]
    rna_all_sub = np.array([pos_in_subset[t["track_index"]] for t in rna_tracks], dtype=int)
    gtex_sub = np.array([pos_in_subset[ti] for ti in GTEX_LIVER_TRACK_INDICES
                         if ti in pos_in_subset], dtype=int)
    if len(gtex_sub) == 0:
        log("FATAL: no GTEx-liver tracks resolved")
        sys.exit(3)
    rc_perm, n_sw, n_up = build_rc_track_perm(liver)
    log(f"liver tracks: {len(liver)} ({len(rna_tracks)} rna_liver); "
        f"GTEx-primary indices={[liver[j]['track_index'] for j in gtex_sub]}; "
        f"RC perm {n_sw} swapped / {n_up} unpaired")

    # ---- models (verbatim 64b) ---------------------------------------------
    import torch
    import pysam
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"torch {torch.__version__} | device={dev} | cuda={torch.cuda.is_available()}")
    from borzoi_pytorch import Borzoi
    models = []
    for m in args.folds:
        log(f"loading Borzoi fold: {m}")
        models.append(Borzoi.from_pretrained(m).to(dev).eval())
    n_folds = len(models)
    fasta = pysam.FastaFile(args.fasta)
    n_targets = N_HUMAN_TRACKS

    def predict_liver(model, oh):
        x = torch.from_numpy(oh).unsqueeze(0).to(dev)
        with torch.no_grad():
            try:
                out = model(x)
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                with torch.autocast("cuda", dtype=torch.bfloat16):
                    out = model(x)
        if isinstance(out, (tuple, list)):
            out = out[0]
        out = out.float().cpu()[0]
        mode, a, b = find_axes((1,) + tuple(out.shape), n_targets)
        pred = out.numpy() if mode == "bins_first" else out.numpy().T
        return pred[:, liver_idx]

    def rc_to_fwd_frame(pred_rc):
        return pred_rc[::-1, :][:, rc_perm]

    def fwd_rc_cov(model, oh_fwd, oh_rc):
        pf = predict_liver(model, oh_fwd)
        pr = rc_to_fwd_frame(predict_liver(model, oh_rc))
        return 0.5 * (pf + pr)

    def score_gene(chrom, v0, hg38_ref, hg38_alt, rec):
        """64b logSED recipe for ONE gene-TSS-centered window. Returns dict or None."""
        tss1, strand = gene_tss(rec)
        tss0 = tss1 - 1
        dist = abs(v0 - tss0)
        if dist >= IN_HALF:
            return None
        window = fetch_window(fasta, chrom, tss0, SEQ_LEN)
        vi = v0 - tss0 + IN_HALF
        win_base = window[vi]
        mism = (win_base != hg38_ref)
        alt_window = window[:vi] + hg38_alt + window[vi + 1:]
        ebins = np.array(sorted(exon_union_bins(rec, tss0)), dtype=int)
        if len(ebins) == 0:
            return {"logsed": float("nan"), "fold_agree": 0, "n_exon_bins": 0,
                    "tss_dist": dist, "strand": strand, "note": "no_exon_bins",
                    "perfold": ""}
        ref_oh = one_hot(window)
        alt_oh = one_hot(alt_window)
        ref_rc = rc_one_hot(ref_oh)
        alt_rc = rc_one_hot(alt_oh)
        logsed_pf = np.empty((n_folds, len(liver_idx)), dtype=np.float64)
        for fi, model in enumerate(models):
            ref_cov = fwd_rc_cov(model, ref_oh, ref_rc)
            alt_cov = fwd_rc_cov(model, alt_oh, alt_rc)
            ref_es = np.clip(ref_cov[ebins, :].sum(axis=0), 0.0, None)
            alt_es = np.clip(alt_cov[ebins, :].sum(axis=0), 0.0, None)
            logsed_pf[fi] = np.log(alt_es + 1.0) - np.log(ref_es + 1.0)
        logsed_ens = logsed_pf.mean(axis=0)
        primary = float(np.mean(logsed_ens[gtex_sub]))
        perfold_liver = logsed_pf[:, gtex_sub].mean(axis=1)
        s = np.sign(primary)
        fold_agree = int(np.sum(np.sign(perfold_liver) == s)) if s != 0 else 0
        return {"logsed": primary, "fold_agree": fold_agree, "n_exon_bins": len(ebins),
                "tss_dist": dist, "strand": strand,
                "note": ("window_base_mismatch" if mism else ""),
                "perfold": ",".join(f"{v:.6f}" for v in perfold_liver)}

    # ---- score every candidate for every lead ------------------------------
    nom_cols = ["locus_id", "variant_id_hg19", "chr", "pos_hg38", "hg38_ref", "hg38_alt",
                "lead_pip", "best_gene", "best_gene_ensembl", "best_gene_biotype",
                "borzoi_abs_logsed", "borzoi_logsed_signed", "magnitude_percentile",
                "exceeds_p80", "fold_agreement", "best_gene_tss_dist",
                "candidate_genes_scored", "n_candidates_in_window",
                "substrate_gene", "substrate_gene_abs_logsed",
                "substrate_gene_percentile", "note"]
    cand_cols = ["locus_id", "variant_id_hg19", "chr", "pos_hg38", "gene", "ensembl",
                 "biotype", "borzoi_logsed_signed", "borzoi_abs_logsed",
                 "magnitude_percentile", "fold_agreement", "tss_dist",
                 "n_exon_bins", "is_best", "is_substrate_gene", "note"]
    nom_rows = []
    cand_rows = []
    n_p80 = 0
    t_start = time.time()

    for li, lead in enumerate(leads):
        chrom = lead["chr"]
        pos1 = int(lead["pos_hg38"])
        v0 = pos1 - 1
        hg38_ref = lead["ref_hg38"].upper()
        hg38_alt = lead["alt_hg38"].upper()
        sub_ens = strip_ver(lead.get("ensembl_substrate", "") or "")
        sub_sym = lead.get("gene_substrate", "") or ""

        cands = candidates_in_window(by_chr, chrom, v0)
        n_in_window = len(cands)

        scored = []
        for rec in cands:
            r = score_gene(chrom, v0, hg38_ref, hg38_alt, rec)
            if r is None:
                continue
            is_sub = (rec["sid"] == sub_ens) or (rec["gene_name"] == sub_sym)
            abslog = abs(r["logsed"]) if np.isfinite(r["logsed"]) else float("nan")
            pctile = pct_of(ref, abslog)
            scored.append({
                "rec": rec, "logsed": r["logsed"], "abslog": abslog,
                "fold_agree": r["fold_agree"], "tss_dist": r["tss_dist"],
                "n_exon_bins": r["n_exon_bins"], "note": r["note"],
                "pctile": pctile, "is_sub": is_sub,
            })

        valid = [s for s in scored if np.isfinite(s["abslog"])]
        n_scored = len(valid)

        # substrate (pre-specified) gene comparison
        sub_hit = next((s for s in scored if s["is_sub"] and np.isfinite(s["abslog"])), None)

        if not valid:
            nom_rows.append({c: "" for c in nom_cols} | {
                "locus_id": lead["locus_id"], "variant_id_hg19": lead["variant_id_hg19"],
                "chr": chrom, "pos_hg38": pos1, "hg38_ref": hg38_ref, "hg38_alt": hg38_alt,
                "lead_pip": lead["lead_pip"], "candidate_genes_scored": 0,
                "n_candidates_in_window": n_in_window, "substrate_gene": sub_sym,
                "note": "no_scorable_candidate"})
            for s in scored:
                cand_rows.append(_cand_row(lead, chrom, pos1, s, False))
            log(f"[{li+1}/{len(leads)}] {lead['locus_id']} {lead['variant_id_hg19']} "
                f"NO scorable candidate (in_window={n_in_window})")
            continue

        best = max(valid, key=lambda s: s["abslog"])
        exceeds = best["abslog"] > ref_p80
        if exceeds:
            n_p80 += 1

        nom_rows.append({
            "locus_id": lead["locus_id"], "variant_id_hg19": lead["variant_id_hg19"],
            "chr": chrom, "pos_hg38": pos1, "hg38_ref": hg38_ref, "hg38_alt": hg38_alt,
            "lead_pip": lead["lead_pip"],
            "best_gene": best["rec"]["gene_name"], "best_gene_ensembl": best["rec"]["sid"],
            "best_gene_biotype": best["rec"]["biotype"],
            "borzoi_abs_logsed": f"{best['abslog']:.6f}",
            "borzoi_logsed_signed": f"{best['logsed']:.6f}",
            "magnitude_percentile": f"{best['pctile']:.2f}",
            "exceeds_p80": str(bool(exceeds)).upper(),
            "fold_agreement": best["fold_agree"],
            "best_gene_tss_dist": best["tss_dist"],
            "candidate_genes_scored": n_scored,
            "n_candidates_in_window": n_in_window,
            "substrate_gene": sub_sym,
            "substrate_gene_abs_logsed": (f"{sub_hit['abslog']:.6f}" if sub_hit else ""),
            "substrate_gene_percentile": (f"{sub_hit['pctile']:.2f}" if sub_hit else ""),
            "note": ("substrate_gene_not_in_candidates" if sub_hit is None else ""),
        })
        for s in scored:
            cand_rows.append(_cand_row(lead, chrom, pos1, s, s is best))

        if li % 10 == 0 or exceeds:
            log(f"[{li+1}/{len(leads)}] {lead['locus_id']} pip={lead['lead_pip']} "
                f"-> {best['rec']['gene_name']} |logSED|={best['abslog']:.4f} "
                f"pct={best['pctile']:.0f} agree={best['fold_agree']}/4 "
                f"cand={n_scored}/{n_in_window} {'**>p80**' if exceeds else ''} "
                f"[{time.time()-t_start:.0f}s]")

    # ---- write -------------------------------------------------------------
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=nom_cols, delimiter="\t")
        w.writeheader()
        w.writerows(nom_rows)
    with open(args.out_candidates, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cand_cols, delimiter="\t")
        w.writeheader()
        w.writerows(cand_rows)

    n_loci = len(nom_rows)
    n_scored_loci = sum(1 for r in nom_rows if r.get("candidate_genes_scored"))
    log("================= SUMMARY =================")
    log(f"eQTL-absent signal leads          : {len(leads)}")
    log(f"loci with >=1 scorable candidate  : {n_scored_loci}")
    log(f"candidate gene-scorings (long)    : {len(cand_rows)}")
    log(f"reference (eQTL-present) n / p80   : {len(ref)} / {ref_p80:.4f}")
    log(f"loci exceeding 80th magnitude pct : {n_p80}  "
        f"({100.0*n_p80/max(n_scored_loci,1):.1f}% of scored) => FUNCTIONAL-IMPACT NOMINATIONS")
    log(f"wrote {args.out} ({n_loci} loci)")
    log(f"wrote {args.out_candidates} ({len(cand_rows)} candidate rows)")
    log("==========================================")
    sys.exit(0 if n_scored_loci > 0 else 2)


def _cand_row(lead, chrom, pos1, s, is_best):
    rec = s["rec"]
    return {
        "locus_id": lead["locus_id"], "variant_id_hg19": lead["variant_id_hg19"],
        "chr": chrom, "pos_hg38": pos1,
        "gene": rec["gene_name"], "ensembl": rec["sid"], "biotype": rec["biotype"],
        "borzoi_logsed_signed": (f"{s['logsed']:.6f}" if np.isfinite(s["logsed"]) else "nan"),
        "borzoi_abs_logsed": (f"{s['abslog']:.6f}" if np.isfinite(s["abslog"]) else "nan"),
        "magnitude_percentile": (f"{s['pctile']:.2f}" if np.isfinite(s["pctile"]) else "nan"),
        "fold_agreement": s["fold_agree"], "tss_dist": s["tss_dist"],
        "n_exon_bins": s["n_exon_bins"],
        "is_best": str(bool(is_best)).upper(),
        "is_substrate_gene": str(bool(s["is_sub"])).upper(),
        "note": s["note"],
    }


if __name__ == "__main__":
    main()
