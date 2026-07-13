#!/usr/bin/env python3
"""
79_score_panels.py — Phase-5 HARDENED modality-matched feature scoring for the
caQTL (and sQTL) DIRECTION benchmark leads that lack sequence-model features.

The POWER fix. src/78 built the leakage-tiered, LD-de-pseudoreplicated, orientation-
audited SIGN-label tables (direction_labels_{caqtl,sqtl}.tsv). Most caQTL/sQTL leads
were never scored by the eQTL-substrate feature runs (62/64b/65/68 scored the Broadaway
eQTL substrate). This script scores the caQTL/sQTL LEADS with MODALITY-MATCHED heads:

  caQTL (ACCESSIBILITY phenotype -> accessibility features ONLY, never gene-expression):
     - ChromBPNet HepG2 ATAC   (signed cbp_logfc; reuse src/68 vendored variant-scorer)
     - AlphaGenome ATAC + DNASE heads (signed deltas; reuse src/65 machinery)
     - Borzoi liver ATAC + DNASE tracks (signed logSED-style delta; reuse src/62)
  sQTL (leafcutter intron-excision phenotype -> SPLICE features; LEAKAGE_EXPLORATORY):
     - AlphaGenome SpliceJunction + SpliceSiteUsage (signed-extreme raw + quantile)
     - SpliceAI DS_AG/AL/DG/DL + max_ds (magnitude; reuse src/61 CLI pattern)

NON-NEGOTIABLE protocol (every item is a 5-lens review fix):
  (1) LEAKAGE TIERS carried from the labels: caQTL=clean_secondary CERTIFIABLE;
      sQTL=leakage_exploratory (GTEx in Borzoi+AG training) -> EXPLORATORY / magnitude
      only, NEVER certified. Tier column propagated into every feature row.
  (2) FEATURE-MODALITY MATCHING: caQTL scored with ACCESSIBILITY heads (NOT gene-expr);
      sQTL scored with SPLICE heads. No RNA/expression head is applied to caQTL.
  (3) LEADS only (src/78 already restricted to is_signal_lead / lead-per-peak /
      lead-intron-per-variant); this script consumes those tables verbatim.
  (4) ORIENTATION: every model feature is scored strictly hg38 REF->ALT, then re-expressed
      into the LABEL frame (the src/78 ref->alt). AG/Borzoi resolve alleles against the
      FASTA and may swap; we record hg38_ref/hg38_alt + an orient_mult in {+1,-1} and the
      merge step multiplies signed features by it so a strand/allele swap is never
      silently sign-wrong. Strand-ambiguous A/T & C/G already hard-dropped in src/78.
      SORT1 rs12740374 (G>T, alt creates C/EBP site) unit-checked -> accessibility UP.
  (5) NEVER a scored convergence channel. caQTL/sQTL = same genetic axis; these features
      are labels/corroboration only, APPLY-ONLY firewall enforced downstream, never summed.

ADDITIVE ONLY. New file (script 79). Reads direction_labels_*.tsv + reuses scorer
machinery from 62/65/68/61; writes direction_features_{caqtl,sqtl}.tsv (+ per-panel
intermediates). Does NOT touch 74-77 / 46d / 78 / 27a / 60-73.

Panels (each its own env; dispatch via --panel):
  chrombpnet-prep      write caQTL VARLIST + genemap for the vendored variant-scorer
  chrombpnet-assemble  fold-mean annotations -> signed accessibility feature (rnaseq env)
  alphagenome          caQTL ATAC/DNASE + sQTL splice via AlphaGenome API (alphagenome env)
  borzoi               caQTL liver ATAC/DNASE signed delta (borzoi env, GPU)
  spliceai             sQTL SpliceAI DS (spliceai_throwaway env)
  merge                assemble direction_features_{caqtl,sqtl}.tsv (rnaseq env)

Author: coding-upgrade agent (Phase 5).
"""
import argparse
import importlib.util
import os
import sys
import time

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------- paths
ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
SRC = os.path.join(ROOT, "GWAS/finemapping/src")
FASTA = "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
TRACKS = os.path.join(SF, "tracks/targets_human.txt")
LIVER_TRACKS = os.path.join(SF, "tracks/liver_tracks_selected.csv")

LABELS_CAQTL = os.path.join(SF, "direction_labels_caqtl.tsv")
LABELS_SQTL = os.path.join(SF, "direction_labels_sqtl.tsv")

# per-panel intermediates
CAQTL_CBP = os.path.join(SF, "direction_features_caqtl.chrombpnet.tsv")
CAQTL_AG = os.path.join(SF, "direction_features_caqtl.alphagenome.tsv")
CAQTL_BZ = os.path.join(SF, "direction_features_caqtl.borzoi.tsv")
SQTL_AG = os.path.join(SF, "direction_features_sqtl.alphagenome.tsv")
SQTL_SAI = os.path.join(SF, "direction_features_sqtl.spliceai.tsv")
OUT_CAQTL = os.path.join(SF, "direction_features_caqtl.tsv")
OUT_SQTL = os.path.join(SF, "direction_features_sqtl.tsv")

# ChromBPNet vendored scorer assets (from src/68)
CBP_DIR = os.path.join(SF, "chrombpnet_hepg2")
CBP_WORK = os.path.join(SF, "chrombpnet_caqtl")   # our new working dir (isolated from 68's)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
BASES = "ACGT"

# SORT1 rs12740374 hg38 chr1:109274968 G>T ; alt T creates C/EBP site -> UP accessibility
SORT1_CHR, SORT1_POS = "1", 109274968
# liver ACCESSIBILITY track indices in targets_human.txt (from liver_tracks_selected.csv)
# ATAC: 2023(Fibro Liver Adrenal) 2035(Hepatocyte) ; DNASE: 1302 1303 1367 1510 1623


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def norm_chr(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def canon(chrom, pos, ref, alt):
    return f"{norm_chr(chrom)}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def load_module(num_name):
    """Import a digit-prefixed sibling script as a module (e.g. '62_borzoi_score')."""
    path = os.path.join(SRC, num_name + ".py")
    spec = importlib.util.spec_from_file_location(num_name.replace("-", "_"), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_labels(which, independent_only=True):
    """Load a src/78 label table. By default restrict to INDEPENDENT LEADS
    (is_peak_lead==True) — protocol item (3): score independent leads, NOT the
    LD-pseudoreplicated variant-level members. caQTL: 498k variant-leads collapse
    to ~11.9k independent peak-leads (one top variant per accessibility peak).
    sQTL/eQTL carry is_peak_lead==True for all lead rows."""
    fp = LABELS_CAQTL if which == "caqtl" else LABELS_SQTL
    if not os.path.exists(fp):
        sys.exit(f"FATAL: label table missing (run src/78 first): {fp}")
    df = pd.read_csv(fp, sep="\t", dtype=str)
    n_all = len(df)
    if independent_only and "is_peak_lead" in df.columns:
        df = df[df["is_peak_lead"].astype(str).str.upper() == "TRUE"].copy()
        log(f"[{which}] independent peak-leads {len(df)}/{n_all} "
            f"(LD-pseudoreplicated members dropped)")
    df["pos_hg38"] = df["pos_hg38"].astype(np.int64)
    df["chr"] = df["chr"].map(norm_chr)
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()
    df["sign"] = df["sign"].astype(float).astype(int)
    df["canon"] = [canon(c, p, r, a) for c, p, r, a in
                   zip(df["chr"], df["pos_hg38"], df["ref"], df["alt"])]
    # SNV-only for base-resolution accessibility models (indels handled only by AG/SpliceAI)
    df["is_snv"] = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
    return df


def subsample(df, max_n, tag):
    """Deterministic subsample for GPU-heavy panels: is_peak_lead first, then p asc."""
    if max_n <= 0 or len(df) <= max_n:
        return df, False
    d = df.copy()
    d["_pl"] = (d.get("is_peak_lead", "True").astype(str).str.upper() == "TRUE").astype(int)
    d["_p"] = pd.to_numeric(d.get("pvalue"), errors="coerce").fillna(1.0)
    d = d.sort_values(["_pl", "_p"], ascending=[False, True], kind="mergesort").head(max_n)
    log(f"[{tag}] SUBSAMPLED {len(df)} -> {max_n} leads (is_peak_lead + p-rank)")
    return d.drop(columns=["_pl", "_p"]), True


def resolve_orient(fasta_base, ref, alt):
    """Return (hg38_ref, hg38_alt, orient_mult, note).
    orient_mult multiplies a fasta-ref->fasta-alt SIGNED feature into the LABEL
    ref->alt frame: +1 if same direction, -1 if the alleles are swapped on hg38."""
    ref, alt, b = ref.upper(), alt.upper(), fasta_base.upper()
    if b == ref:
        return ref, alt, 1, ""
    if b == alt:
        return alt, ref, -1, "hg19_ref_alt_swapped_on_hg38"
    if b == COMP.get(ref):
        return COMP[ref], COMP[alt], 1, "strand_flip"
    if b == COMP.get(alt):
        return COMP[alt], COMP[ref], -1, "strand_flip_swapped"
    return b, None, 0, f"UNRESOLVED_fasta={b}_alleles={ref}/{alt}"


# ============================================================ ChromBPNet prep / assemble
def panel_chrombpnet_prep(args):
    df = load_labels("caqtl")
    df = df[df["is_snv"]].copy()          # ChromBPNet variant-scorer scores SNVs
    df, subsamp = subsample(df, args.max_variants, "chrombpnet")
    os.makedirs(CBP_WORK, exist_ok=True)
    varlist = os.path.join(CBP_WORK, "variants_input.chrombpnet.tsv")
    genemap = os.path.join(CBP_WORK, "variant_genemap.tsv")
    with open(varlist, "w") as vf:
        for _, r in df.iterrows():
            vid = f"chr{r['chr']}_{r['pos_hg38']}_{r['ref']}_{r['alt']}"
            vf.write(f"chr{r['chr']}\t{r['pos_hg38']}\t{r['ref']}\t{r['alt']}\t{vid}\n")
    gm = pd.DataFrame({
        "variant_id": [f"chr{r['chr']}_{r['pos_hg38']}_{r['ref']}_{r['alt']}"
                       for _, r in df.iterrows()],
        "gene": df["feature_id"].values,                 # Currin peak id
        "max_pip": np.nan,
        "var_class": "caqtl_lead",
        "variant_id_hg19": np.nan,
    })
    gm.to_csv(genemap, sep="\t", index=False)
    log(f"[chrombpnet-prep] wrote {len(df)} caQTL SNV leads (subsampled={subsamp})")
    log(f"  VARLIST : {varlist}")
    log(f"  GENEMAP : {genemap}")


def panel_chrombpnet_assemble(args):
    """Reuse src/68_chrombpnet_accessibility.py to assemble, then re-key to canon.
    ChromBPNet scores the VARLIST allele1->allele2 == label ref->alt directly, so
    orient_mult is +1 (no fasta resolution/swap inside the vendored scorer)."""
    ann = args.annotations
    gm = os.path.join(CBP_WORK, "variant_genemap.tsv")
    tmp = os.path.join(CBP_WORK, "chrombpnet_accessibility_caqtl.tsv")
    m68 = load_module("68_chrombpnet_accessibility")
    sys.argv = ["68", "--annotations", ann, "--genemap", gm,
                "--n_folds", str(args.n_folds), "--out", tmp]
    m68.main()
    d = pd.read_csv(tmp, sep="\t", dtype=str)
    d["canon"] = [canon(c, p, a1, a2) for c, p, a1, a2 in
                  zip(d["chr"], d["pos_hg38"], d["allele1"], d["allele2"])]
    out = pd.DataFrame({
        "canon": d["canon"],
        "cbp_logfc_labelframe": pd.to_numeric(d["cbp_logfc"], errors="coerce"),  # already ref->alt
        "cbp_abs_logfc": pd.to_numeric(d["cbp_abs_logfc"], errors="coerce"),
        "cbp_jsd": pd.to_numeric(d["cbp_jsd"], errors="coerce"),
        "cbp_active_quantile": pd.to_numeric(d["cbp_active_quantile"], errors="coerce"),
        "cbp_in_hepg2_peak": d["in_hepg2_peak"],
        "cbp_disrupted_tf_motif": d["disrupted_tf_motif"],
    })
    out.to_csv(CAQTL_CBP, sep="\t", index=False)
    log(f"[chrombpnet-assemble] wrote {CAQTL_CBP} ({len(out)} variants)")
    _sort1_check(out, "cbp_logfc_labelframe", "ChromBPNet")


# ============================================================ AlphaGenome (accessibility + splice)
def panel_alphagenome(args):
    m65 = load_module("65_alphagenome_score")
    import pysam
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
    import grpc

    api_key, key_src = m65.load_api_key(args)
    if api_key is None:
        log("READY — awaiting AlphaGenome API key (env ALPHA_GENOME_API_KEY or "
            "~/.alphagenome_key). Key never logged."); sys.exit(3)
    log(f"API key loaded from {key_src} (len={len(api_key)}; value not shown)")

    seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS[args.seq_length]
    R = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    acc_keys = [k for k in ("ATAC", "DNASE") if k in R]
    spl_keys = [k for k in ("SPLICE_JUNCTIONS", "SPLICE_SITE_USAGE", "SPLICE_SITES") if k in R]
    fasta = pysam.FastaFile(FASTA)
    model = dna_client.create(api_key)

    def score_one(chrom, pos1, hg38_ref, hg38_alt, vid, scorers):
        variant = genome.Variant(chromosome=chrom, position=pos1,
                                 reference_bases=hg38_ref, alternate_bases=hg38_alt, name=vid)
        interval = variant.reference_interval.resize(seq_len)
        last = None
        for attempt in range(args.max_retries):
            try:
                s = model.score_variant(interval=interval, variant=variant,
                                        variant_scorers=scorers,
                                        organism=dna_client.Organism.HOMO_SAPIENS)
                return variant_scorers.tidy_scores([s])
            except grpc.RpcError as e:
                last = e
                if e.code() not in (grpc.StatusCode.UNAVAILABLE,
                                    grpc.StatusCode.RESOURCE_EXHAUSTED,
                                    grpc.StatusCode.DEADLINE_EXCEEDED,
                                    grpc.StatusCode.INTERNAL):
                    raise
                w = min(30.0, 2.0 ** attempt)
                log(f"  transient {e.code()} ({attempt+1}/{args.max_retries}); backoff {w:.0f}s")
                time.sleep(w)
        raise last

    def run_panel(which, scorer_keys, is_acc):
        df = load_labels(which)
        if args.max_variants > 0 and len(df) > args.max_variants:
            df, _ = subsample(df, args.max_variants, f"alphagenome-{which}")
        scorers = [R[k] for k in scorer_keys]
        log(f"[alphagenome-{which}] scoring {len(df)} leads with "
            f"{','.join(scorer_keys)} @ {args.seq_length}")
        rows = []
        n_ok = n_err = n_unres = 0
        for i, (_, r) in enumerate(df.iterrows()):
            chrom = "chr" + r["chr"]; pos1 = int(r["pos_hg38"]); v0 = pos1 - 1
            ref, alt = r["ref"], r["alt"]
            vid = r["variant"]
            base = dict(canon=r["canon"], variant=vid, chr=r["chr"], pos_hg38=pos1,
                        label_ref=ref, label_alt=alt, sign=int(r["sign"]),
                        leakage_tier=r["leakage_tier"], orient_mult=1, note="")
            is_indel = (len(ref) != 1 or len(alt) != 1)
            try:
                fbase = fasta.fetch(chrom, v0, v0 + 1).upper()
            except Exception as e:
                base["note"] = f"fasta_fail:{e}"; rows.append(base); n_err += 1; continue
            if not is_indel:
                hg38_ref, hg38_alt, om, note = resolve_orient(fbase, ref, alt)
                if hg38_alt is None:
                    base["note"] = note; rows.append(base); n_unres += 1; continue
                base["orient_mult"] = om; base["note"] = note
            else:
                hg38_ref, hg38_alt = ref, alt
                base["orient_mult"] = 1  # indels: label frame == provided; AG scores as given
            try:
                tdf = score_one(chrom, pos1, hg38_ref, hg38_alt, vid, scorers)
                ext = m65.extract_scores(tdf, "", "", dbg=False)
                if is_acc:
                    base["ag_atac_delta"] = ext["ag_accessibility_delta"]
                    base["ag_dnase_delta"] = ext["ag_dnase_delta"]
                else:
                    base["ag_splice_score"] = ext["ag_splice_score"]
                    base["ag_splice_quantile"] = ext["ag_splice_quantile"]
                    base["ag_n_splice_rows"] = ext["n_splice_rows"]
                n_ok += 1
            except Exception as e:
                base["note"] = (base["note"] + ";" if base["note"] else "") + \
                    f"score_error:{type(e).__name__}:{e}"; n_err += 1
            rows.append(base)
            if i % 200 == 0:
                log(f"  [{which} {i+1}/{len(df)}] {vid} ok={n_ok} err={n_err} unres={n_unres}")
            time.sleep(args.sleep)
        out = pd.DataFrame(rows)
        outfp = CAQTL_AG if is_acc else SQTL_AG
        out.to_csv(outfp, sep="\t", index=False)
        log(f"[alphagenome-{which}] wrote {outfp} ({len(out)}); ok={n_ok} err={n_err} "
            f"unresolved={n_unres}")
        if is_acc:
            _sort1_check(out, "ag_atac_delta", "AlphaGenome-ATAC", orient=True)
        return out

    if args.which in ("caqtl", "both"):
        run_panel("caqtl", acc_keys, is_acc=True)
    if args.which in ("sqtl", "both"):
        run_panel("sqtl", spl_keys, is_acc=False)


# ============================================================ Borzoi (liver ATAC/DNASE)
def panel_borzoi(args):
    m62 = load_module("62_borzoi_score")
    import torch, pysam
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"torch {torch.__version__} device={dev}")
    from borzoi_pytorch import Borzoi
    model = Borzoi.from_pretrained(args.model).to(dev).eval()
    fasta = pysam.FastaFile(FASTA)
    targets = m62.load_targets(TRACKS)
    n_targets = len(targets)
    liver = m62.load_liver_tracks(LIVER_TRACKS)
    atac_idx = [t["track_index"] for t in liver if t["category"] == "ATAC"]
    dnase_idx = [t["track_index"] for t in liver if t["category"] == "DNASE"]
    log(f"Borzoi liver accessibility tracks: {len(atac_idx)} ATAC + {len(dnase_idx)} DNASE")

    def predict(seq):
        oh = m62.one_hot(seq)
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
        mode, a, b = m62.find_axes((1,) + tuple(out.shape), n_targets)
        return out.numpy() if mode == "bins_first" else out.numpy().T   # (bins, tracks)

    df = load_labels("caqtl")
    df = df[df["is_snv"]].copy()
    df, subsamp = subsample(df, args.max_variants, "borzoi")
    SEQ_LEN = m62.SEQ_LEN
    rows = []
    n_ok = n_unres = 0
    for i, (_, r) in enumerate(df.iterrows()):
        chrom = "chr" + r["chr"]; pos1 = int(r["pos_hg38"]); pos0 = pos1 - 1
        ref, alt = r["ref"], r["alt"]
        base = dict(canon=r["canon"], variant=r["variant"], chr=r["chr"], pos_hg38=pos1,
                    label_ref=ref, label_alt=alt, sign=int(r["sign"]),
                    leakage_tier=r["leakage_tier"], orient_mult=1, note="")
        try:
            seq = m62.fetch_window(fasta, chrom, pos0, SEQ_LEN)
        except Exception as e:
            base["note"] = f"fetch_fail:{e}"; rows.append(base); continue
        center = SEQ_LEN // 2
        fbase = seq[center].upper()
        hg38_ref, hg38_alt, om, note = resolve_orient(fbase, ref, alt)
        base["orient_mult"] = om; base["note"] = note
        if hg38_alt is None:
            rows.append(base); n_unres += 1; continue
        alt_seq = seq[:center] + hg38_alt + seq[center + 1:]
        rp = predict(seq); ap = predict(alt_seq)
        n_bins = rp.shape[0]; cbin = n_bins // 2
        lo = max(0, cbin - args.win_bins); hi = min(n_bins, cbin + args.win_bins + 1)
        dwin = (ap[lo:hi] - rp[lo:hi]).sum(axis=0)   # (tracks,) fasta ref->alt
        base["borzoi_atac_delta"] = float(np.nanmean([dwin[k] for k in atac_idx]))
        base["borzoi_dnase_delta"] = float(np.nanmean([dwin[k] for k in dnase_idx]))
        n_ok += 1
        rows.append(base)
        if i % 50 == 0:
            log(f"  [borzoi {i+1}/{len(df)}] ok={n_ok} unres={n_unres}")
    out = pd.DataFrame(rows)
    out.to_csv(CAQTL_BZ, sep="\t", index=False)
    log(f"[borzoi] wrote {CAQTL_BZ} ({len(out)}); ok={n_ok} unresolved={n_unres} "
        f"subsampled={subsamp}")
    _sort1_check(out, "borzoi_atac_delta", "Borzoi-ATAC", orient=True)


# ============================================================ SpliceAI (sQTL)
def panel_spliceai(args):
    import subprocess, pysam
    spliceai_bin = os.path.join(ROOT, ".envs/spliceai_throwaway/bin/spliceai")
    if not os.path.exists(spliceai_bin):
        log(f"SpliceAI CLI absent ({spliceai_bin}) -> deferred"); return
    df = load_labels("sqtl")
    if args.max_variants > 0 and len(df) > args.max_variants:
        df, _ = subsample(df, args.max_variants, "spliceai-sqtl")
    fa = pysam.FastaFile(FASTA)
    recs = []; key_for = {}
    for _, r in df.iterrows():
        contig = "chr" + r["chr"]; pos = int(r["pos_hg38"])
        try:
            gref = fa.fetch(contig, pos - 1, pos).upper()
        except Exception:
            continue
        a1, a2 = r["ref"], r["alt"]
        if a1 == gref:
            alt = a2
        elif a2 == gref:
            alt = a1
        elif COMP.get(a1) == gref:
            alt = COMP.get(a2)
        elif COMP.get(a2) == gref:
            alt = COMP.get(a1)
        else:
            continue
        if not alt or alt == gref:
            continue
        key = (contig, pos, gref, alt)
        key_for[r["canon"]] = key
        recs.append(key)
    fa.close()
    recs = sorted(set(recs), key=lambda k: (len(k[0]), k[0], k[1]))
    if not recs:
        log("SpliceAI: no sQTL records to score"); return
    vcf_in = os.path.join(SF, "spliceai_sqtl_input.vcf")
    vcf_out = os.path.join(SF, "spliceai_sqtl_output.vcf")
    with open(vcf_in, "w") as vf:
        vf.write("##fileformat=VCFv4.2\n")
        for c in sorted({k[0] for k in recs}, key=lambda c: (len(c), c)):
            vf.write(f"##contig=<ID={c}>\n")
        vf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for (c, p, gr, al) in recs:
            vf.write(f"{c}\t{p}\t{c}_{p}_{gr}_{al}\t{gr}\t{al}\t.\t.\t.\n")
    log(f"SpliceAI: wrote {len(recs)} sQTL records; running CLI (CPU)")
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="", TF_CPP_MIN_LOG_LEVEL="3")
    proc = subprocess.run([spliceai_bin, "-I", vcf_in, "-O", vcf_out,
                           "-R", FASTA, "-A", "grch38"],
                          env=env, capture_output=True, text=True)
    if proc.returncode != 0 or not os.path.exists(vcf_out):
        log(f"SpliceAI FAILED rc={proc.returncode}\n{proc.stderr[-1200:]}"); return
    ds_by = {}
    with open(vcf_out) as vf:
        for line in vf:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, info = f[0], int(f[1]), f[3], f[4], f[7]
            if "SpliceAI=" not in info:
                continue
            best = None
            for field in info.split(";"):
                if not field.startswith("SpliceAI="):
                    continue
                for ann in field[len("SpliceAI="):].split(","):
                    p = ann.split("|")
                    if len(p) < 6:
                        continue
                    try:
                        ds = [float(p[2]), float(p[3]), float(p[4]), float(p[5])]
                    except ValueError:
                        continue
                    if best is None or max(ds) > max(best):
                        best = ds
            if best is not None:
                ds_by[(chrom, pos, ref, alt)] = best
    rows = []
    for _, r in df.iterrows():
        rec = dict(canon=r["canon"], variant=r["variant"], leakage_tier=r["leakage_tier"],
                   spliceai_DS_AG="", spliceai_DS_AL="", spliceai_DS_DG="",
                   spliceai_DS_DL="", spliceai_max_ds="")
        key = key_for.get(r["canon"])
        if key and key in ds_by:
            ag, al, dg, dl = ds_by[key]
            rec.update(spliceai_DS_AG=f"{ag:.4f}", spliceai_DS_AL=f"{al:.4f}",
                       spliceai_DS_DG=f"{dg:.4f}", spliceai_DS_DL=f"{dl:.4f}",
                       spliceai_max_ds=f"{max(ag, al, dg, dl):.4f}")
        rows.append(rec)
    pd.DataFrame(rows).to_csv(SQTL_SAI, sep="\t", index=False)
    log(f"SpliceAI: wrote {SQTL_SAI} ({len(rows)} rows; "
        f"{sum(1 for x in rows if x['spliceai_max_ds'])} scored)")


# ============================================================ merge
def panel_merge(args):
    # ---- caQTL: label base + accessibility features, re-expressed into label frame ----
    lab = load_labels("caqtl")
    base = lab[["canon", "variant", "chr", "pos_hg38", "ref", "alt", "sign", "effect",
                "maf", "tss_or_peak_dist", "feature_id", "pvalue", "is_peak_lead",
                "leakage_tier"]].copy()
    n_feat = {}
    if os.path.exists(CAQTL_CBP):
        c = pd.read_csv(CAQTL_CBP, sep="\t")
        base = base.merge(c, on="canon", how="left")
        n_feat["chrombpnet(cbp_logfc)"] = c["cbp_logfc_labelframe"].notna().sum()
    if os.path.exists(CAQTL_AG):
        a = pd.read_csv(CAQTL_AG, sep="\t")
        a["ag_atac_delta_labelframe"] = pd.to_numeric(a.get("ag_atac_delta"), errors="coerce") \
            * pd.to_numeric(a["orient_mult"], errors="coerce")
        a["ag_dnase_delta_labelframe"] = pd.to_numeric(a.get("ag_dnase_delta"), errors="coerce") \
            * pd.to_numeric(a["orient_mult"], errors="coerce")
        base = base.merge(a[["canon", "ag_atac_delta_labelframe",
                             "ag_dnase_delta_labelframe", "orient_mult"]]
                          .rename(columns={"orient_mult": "ag_orient_mult"}),
                          on="canon", how="left")
        n_feat["alphagenome(ag_atac)"] = a["ag_atac_delta_labelframe"].notna().sum()
    if os.path.exists(CAQTL_BZ):
        b = pd.read_csv(CAQTL_BZ, sep="\t")
        b["borzoi_atac_delta_labelframe"] = pd.to_numeric(b.get("borzoi_atac_delta"),
                                                          errors="coerce") \
            * pd.to_numeric(b["orient_mult"], errors="coerce")
        b["borzoi_dnase_delta_labelframe"] = pd.to_numeric(b.get("borzoi_dnase_delta"),
                                                           errors="coerce") \
            * pd.to_numeric(b["orient_mult"], errors="coerce")
        base = base.merge(b[["canon", "borzoi_atac_delta_labelframe",
                            "borzoi_dnase_delta_labelframe", "orient_mult"]]
                          .rename(columns={"orient_mult": "borzoi_orient_mult"}),
                          on="canon", how="left")
        n_feat["borzoi(borzoi_atac)"] = b["borzoi_atac_delta_labelframe"].notna().sum()
    base["modality"] = "accessibility"
    base.to_csv(OUT_CAQTL, sep="\t", index=False)
    log(f"WROTE {OUT_CAQTL} ({len(base)} caQTL leads)")
    for k, v in n_feat.items():
        log(f"   caQTL feature coverage {k}: {v}/{len(base)}")
    _report_sign_concord(base, "caqtl")

    # ---- sQTL: label base + splice features (magnitude; leakage_exploratory) ----
    labs = load_labels("sqtl")
    bs = labs[["canon", "variant", "chr", "pos_hg38", "ref", "alt", "sign", "effect",
               "maf", "tss_or_peak_dist", "feature_id", "pvalue", "leakage_tier"]].copy()
    sn = {}
    if os.path.exists(SQTL_AG):
        a = pd.read_csv(SQTL_AG, sep="\t")
        keep = [c for c in ("canon", "ag_splice_score", "ag_splice_quantile",
                            "ag_n_splice_rows") if c in a.columns]
        bs = bs.merge(a[keep], on="canon", how="left")
        sn["alphagenome(ag_splice)"] = pd.to_numeric(
            a.get("ag_splice_score"), errors="coerce").notna().sum()
    if os.path.exists(SQTL_SAI):
        s = pd.read_csv(SQTL_SAI, sep="\t")
        bs = bs.merge(s[["canon", "spliceai_max_ds", "spliceai_DS_AG", "spliceai_DS_AL",
                        "spliceai_DS_DG", "spliceai_DS_DL"]], on="canon", how="left")
        sn["spliceai(max_ds)"] = pd.to_numeric(
            s.get("spliceai_max_ds"), errors="coerce").notna().sum()
    bs["modality"] = "splice"
    bs["leakage_caveat"] = "GTEx in Borzoi+AG training -> EXPLORATORY, magnitude-only, never certified"
    bs.to_csv(OUT_SQTL, sep="\t", index=False)
    log(f"WROTE {OUT_SQTL} ({len(bs)} sQTL leads)")
    for k, v in sn.items():
        log(f"   sQTL feature coverage {k}: {v}/{len(bs)}")


def _sort1_check(df, col, name, orient=False):
    hit = df[df["canon"].astype(str).str.startswith(f"{SORT1_CHR}:{SORT1_POS}:")]
    if len(hit) == 0:
        log(f"   SORT1 unit-check [{name}]: not in this panel's lead set"); return
    r = hit.iloc[0]
    v = pd.to_numeric(pd.Series([r[col]]), errors="coerce").iloc[0]
    if orient and "orient_mult" in r:
        v = v * float(r["orient_mult"])
    ok = "PASS(up)" if (v == v and v > 0) else "check"
    log(f"   SORT1 unit-check [{name}] label-frame {col}={v} -> {ok} "
        f"(expect UP accessibility; alt T creates C/EBP site)")


def _report_sign_concord(base, which):
    """Descriptive only (NOT a scored channel): sign agreement of each accessibility
    feature with the measured caQTL sign, among leads that carry the feature."""
    for col in ("cbp_logfc_labelframe", "ag_atac_delta_labelframe",
                "borzoi_atac_delta_labelframe"):
        if col not in base.columns:
            continue
        v = pd.to_numeric(base[col], errors="coerce")
        m = v.notna() & base["sign"].notna()
        if m.sum() == 0:
            continue
        agree = (np.sign(v[m]) == np.sign(base["sign"][m])).mean()
        log(f"   [{which}] sign-concordance {col} vs measured caQTL sign: "
            f"{agree:.3f} (n={int(m.sum())})  [descriptive corroboration only]")


# ============================================================ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True,
                    choices=["chrombpnet-prep", "chrombpnet-assemble", "alphagenome",
                             "borzoi", "spliceai", "merge"])
    ap.add_argument("--which", default="both", choices=["caqtl", "sqtl", "both"],
                    help="alphagenome panel: which label set(s) to score")
    ap.add_argument("--max-variants", type=int, default=0,
                    help="cap leads (GPU panels); 0=all. Deterministic peak-lead+p subsample")
    ap.add_argument("--seq-length", default="SEQUENCE_LENGTH_500KB",
                    choices=["SEQUENCE_LENGTH_16KB", "SEQUENCE_LENGTH_100KB",
                             "SEQUENCE_LENGTH_500KB", "SEQUENCE_LENGTH_1MB"])
    ap.add_argument("--model", default="johahi/borzoi-replicate-0")
    ap.add_argument("--win-bins", type=int, default=16)
    ap.add_argument("--n-folds", type=int, default=5)
    ap.add_argument("--annotations", default=None,
                    help="chrombpnet-assemble: vendored fold-mean annotations tsv")
    ap.add_argument("--key-file", default=None)
    ap.add_argument("--sleep", type=float, default=0.15)
    ap.add_argument("--max-retries", type=int, default=5)
    args = ap.parse_args()
    os.makedirs(SF, exist_ok=True)

    if args.panel == "chrombpnet-prep":
        panel_chrombpnet_prep(args)
    elif args.panel == "chrombpnet-assemble":
        if not args.annotations:
            sys.exit("chrombpnet-assemble requires --annotations")
        panel_chrombpnet_assemble(args)
    elif args.panel == "alphagenome":
        panel_alphagenome(args)
    elif args.panel == "borzoi":
        panel_borzoi(args)
    elif args.panel == "spliceai":
        panel_spliceai(args)
    elif args.panel == "merge":
        panel_merge(args)


if __name__ == "__main__":
    main()
