#!/usr/bin/env python3
"""
64b_borzoi_eqtl_logsed.py  --  DEFINITIVE on-recipe Borzoi GENE-LEVEL logSED eQTL
scorer (Stage-2a seqfunc). Pre-registered on-recipe rerun of the FAILED single-fold
linear-SED preview (64_borzoi_eqtl_sed.py, sign-auROC 0.526).

This is the ONE decisive Borzoi direction gate. ADDITIVE ONLY -- new file; it does
NOT edit any existing script (62/64/46d/78/27a). It reuses 64_borzoi_eqtl_sed.py's
gene-model / windowing / polarity machinery verbatim and layers the four canonical
on-recipe fixes on top.

THE FOUR ON-RECIPE FIXES (all applied; the whole point of this rerun):
  1. logSED, not linear SED.  Per RNA track, per fold:
        logSED = log(sum_exon ALT_cov + 1) - log(sum_exon REF_cov + 1)
     (pseudocount 1; the Borzoi/Kelley canonical, bioRxiv 2023.08.30.555582).
     Computed per track, THEN aggregated over tracks.
  2. TISSUE-MATCHED aggregation.  The truth is measured liver microarray eQTL, so
     the matched signal is the GTEx-LIVER RNA track(s): indices 7563/7564/7565
     (identifier GTEX-*, description "RNA:liver", UNSTRANDED -- verified from
     tracks/liver_tracks_selected.csv). PRIMARY score = mean logSED over the GTEx
     liver tracks (borzoi_logsed_liver). SECONDARY = mean logSED over ALL 58
     rna_liver tracks (borzoi_logsed_allliver) -- reported so we can SEE whether the
     58 heterogeneous (mostly HepG2) tracks wash out the matched signal.
  3. 4-FOLD ENSEMBLE.  Average across johahi/borzoi-replicate-{0,1,2,3}.
     Convention (stated): per fold we average the FORWARD + REVERSE-COMPLEMENT
     predicted coverage, exon-sum, and compute that fold's logSED; the four fold
     logSED values are then AVERAGED (the baskerville per-model-score ensemble).
     (Track-mean and fold-mean are linear and commute; only the within-fold
     log(sum+1) is nonlinear, hence computed per fold.) This also yields the four
     per-fold signs for the free confidence flag below.
  4. FORWARD + REVERSE-COMPLEMENT.  Each sequence is scored in both orientations and
     averaged (baskerville borzoi_snp_gene convention). The RC prediction is mapped
     back to the forward genomic frame: bin axis reversed AND stranded (+/-) track
     pairs swapped (unstranded tracks -- incl. all three GTEx liver -- map to self,
     so the PRIMARY is unaffected by any strand-pairing detail). Shifts skipped.

BONUS CONFIDENCE SIGNAL: fold_sign_agreement (0-4) = number of the four per-fold
liver logSED values whose sign equals the ensemble (primary) sign. 4/4 => all folds
agree on direction. borzoi_logsed_per_fold carries the four fold values.

Everything else is IDENTICAL to 64_borzoi_eqtl_sed.py: gene-TSS-centered 524,288 bp
window, |variant - TSS| < 262,144 in-window filter, exon-union output-bin mapping,
STRICT hg38REF -> hg38ALT polarity (hg38_ref/hg38_alt recorded, NOT pre-oriented; the
benchmark orients to the eQTL effect allele), gencode v49 gene resolution, SNV-only
(indels skipped), SORT1 rs12740374 G>T positive-control anchor.

Borzoi geometry (fixed): input 524,288 bp; output 6,144 bins x 32 bp = central
196,608 bp; 7,611 human tracks.
"""

import argparse
import csv
import gzip
import os
import sys
import time

import numpy as np

# ---- Borzoi geometry (fixed by the trained model) --------------------------
SEQ_LEN = 524_288
N_HUMAN_TRACKS = 7_611
N_OUT_BINS = 6_144
BIN_BP = 32
OUT_SPAN = N_OUT_BINS * BIN_BP          # 196,608 bp
OUT_HALF = OUT_SPAN // 2                # 98,304 bp
IN_HALF = SEQ_LEN // 2                  # 262,144 bp
BASES = "ACGT"
BASE_IDX = {b: i for i, b in enumerate(BASES)}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

# GTEx-liver RNA tracks (tissue-matched to the liver-eQTL truth). VERIFIED from
# tracks/liver_tracks_selected.csv: three GTEX-* identifiers, assay rna_liver,
# description "RNA:liver", UNSTRANDED (no +/- suffix).
GTEX_LIVER_TRACK_INDICES = (7563, 7564, 7565)

FOLDS = ("johahi/borzoi-replicate-0", "johahi/borzoi-replicate-1",
         "johahi/borzoi-replicate-2", "johahi/borzoi-replicate-3")

# Positive-control anchor: SORT1 rs12740374 chr1:109274968 hg38 G>T -- minor T
# allele CREATES a C/EBP site and INCREASES hepatic SORT1 expression (Musunuru
# 2010). => logSED(G->T) on liver RNA tracks MUST be POSITIVE.
SORT1_ANCHOR = {
    "variant_id": "SORT1_rs12740374",
    "chr": "1",
    "pos_hg38": 109_274_968,
    "ref": "G",
    "alt": "T",
    "gene": "SORT1",
    "ensembl": "ENSG00000134243",
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---- reused verbatim from 64_borzoi_eqtl_sed.py ---------------------------
def one_hot(seq):
    """seq (str, len SEQ_LEN) -> float32 (4, SEQ_LEN), ACGT order; N -> zeros."""
    arr = np.zeros((4, len(seq)), dtype=np.float32)
    idx = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    for b, i in BASE_IDX.items():
        arr[i] = (idx == ord(b)).astype(np.float32)
    return arr


def rc_one_hot(oh):
    """Reverse-complement of a (4, L) ACGT one-hot: reverse length, swap
    A<->T (rows 0<->3) and C<->G (rows 1<->2). All-zero (N) columns stay zero."""
    return oh[[3, 2, 1, 0], ::-1].copy()


def fetch_window(fasta, chrom, pos0, seq_len):
    """seq_len bp so 0-based pos0 lands at exact center seq_len//2. N-pad edges."""
    half = seq_len // 2
    start = pos0 - half
    end = start + seq_len
    clen = fasta.get_reference_length(chrom)
    fstart = max(0, start)
    fend = min(clen, end)
    core = fasta.fetch(chrom, fstart, fend).upper()
    left_pad = fstart - start
    right_pad = end - fend
    seq = ("N" * left_pad) + core + ("N" * right_pad)
    assert len(seq) == seq_len, f"window len {len(seq)} != {seq_len}"
    return seq


def find_axes(out_shape, n_targets):
    a, b = out_shape[1], out_shape[2]
    if a == n_targets:
        return ("tracks_first", a, b)
    if b == n_targets:
        return ("bins_first", a, b)
    if a > b:
        return ("tracks_first", a, b)
    return ("bins_first", a, b)


# ---- track manifest --------------------------------------------------------
def load_liver_tracks(path):
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            ident = r["identifier"]
            if ident.endswith("+"):
                strand = "+"
            elif ident.endswith("-"):
                strand = "-"
            else:
                strand = "."
            out.append({
                "track_index": int(r["track_index"]),
                "identifier": ident,
                "assay": r["assay"],
                "category": r["category"],
                "description": r["description"],
                "strand": strand,
            })
    return out


def build_rc_track_perm(liver):
    """Permutation p over the liver subset so that, after reversing the RC
    prediction's bin axis, out[:, j] = rc_binrev[:, p[j]] restores the forward
    genomic frame. Stranded (+/-) identifier pairs are swapped; unstranded tracks
    (incl. all three GTEx liver) map to self."""
    n = len(liver)
    stem_strand_to_pos = {}
    for j, t in enumerate(liver):
        ident = t["identifier"]
        if ident.endswith("+") or ident.endswith("-"):
            stem_strand_to_pos[(ident[:-1], ident[-1])] = j
    perm = np.arange(n, dtype=int)
    n_swapped = 0
    n_unpaired = 0
    for j, t in enumerate(liver):
        ident = t["identifier"]
        if ident.endswith("+"):
            partner = stem_strand_to_pos.get((ident[:-1], "-"))
        elif ident.endswith("-"):
            partner = stem_strand_to_pos.get((ident[:-1], "+"))
        else:
            partner = None  # unstranded -> self
        if partner is not None:
            perm[j] = partner
            n_swapped += 1
        elif ident.endswith("+") or ident.endswith("-"):
            n_unpaired += 1  # stranded but partner not in subset -> self (rare)
    return perm, n_swapped, n_unpaired


# ---- GTF gene model (verbatim from 64) ------------------------------------
def _attr(attrs, key):
    tok = f'{key} "'
    i = attrs.find(tok)
    if i < 0:
        return None
    i += len(tok)
    j = attrs.find('"', i)
    if j < 0:
        return None
    return attrs[i:j]


def strip_ver(ens):
    return ens.split(".")[0] if ens else ens


def parse_gtf_genes(gtf_path, want_ens, want_sym):
    want_ens = set(want_ens)
    want_sym = set(want_sym)
    genes = {}
    sym2sid = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    n_lines = 0
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line and line[0] == "#":
                continue
            n_lines += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat = f[2]
            if feat not in ("gene", "transcript", "exon"):
                continue
            attrs = f[8]
            gid = _attr(attrs, "gene_id")
            sid = strip_ver(gid)
            gname = _attr(attrs, "gene_name")
            if sid not in want_ens and (gname not in want_sym):
                continue
            chrom = f[0]
            strand = f[6]
            start = int(f[3])
            end = int(f[4])
            key = (sid, chrom)
            rec = genes.get(key)
            if rec is None:
                rec = {
                    "sid": sid, "chrom": chrom, "strand": strand,
                    "gene_start": None, "gene_end": None,
                    "gene_name": gname, "mane_tss": None,
                    "exons": [],
                }
                genes[key] = rec
                if gname:
                    sym2sid.setdefault(gname, set()).add(sid)
            if feat == "gene":
                rec["gene_start"] = start
                rec["gene_end"] = end
                rec["strand"] = strand
                rec["gene_name"] = gname
            elif feat == "transcript":
                if 'tag "MANE_Select"' in attrs:
                    rec["mane_tss"] = start if strand == "+" else end
            elif feat == "exon":
                rec["exons"].append((start, end))
    log(f"GTF: parsed {n_lines:,} feature lines; matched "
        f"{len(genes)} (gene,chrom) records for "
        f"{len({k[0] for k in genes})} distinct gene ids")
    return genes, sym2sid


def resolve_gene(genes, sym2sid, ensembl, symbol, want_chrom):
    sid = strip_ver(ensembl)
    rec = genes.get((sid, want_chrom))
    if rec is not None:
        return rec, "ensembl_chrom"
    cands = [genes[k] for k in genes if k[0] == sid]
    if cands:
        cands.sort(key=lambda r: (0 if r["chrom"] == want_chrom else
                                  1 if r["chrom"].startswith("chr") and "_" not in r["chrom"]
                                  else 2))
        return cands[0], "ensembl_anychrom"
    for s in sym2sid.get(symbol, ()):
        rec = genes.get((s, want_chrom))
        if rec is not None:
            return rec, "symbol_chrom"
    return None, "unresolved"


def gene_tss(rec):
    strand = rec["strand"]
    if rec["mane_tss"] is not None:
        return rec["mane_tss"], strand
    if strand == "+":
        return rec["gene_start"], strand
    return rec["gene_end"], strand


def exon_union_bins(rec, tss0):
    out_start = tss0 - OUT_HALF
    bins = set()
    for (s1, e1) in rec["exons"]:
        es0 = s1 - 1
        ee0 = e1
        first = (es0 - out_start) // BIN_BP
        last = (ee0 - 1 - out_start) // BIN_BP
        if last < 0 or first >= N_OUT_BINS:
            continue
        first = max(0, first)
        last = min(N_OUT_BINS - 1, last)
        bins.update(range(first, last + 1))
    return bins


def resolve_hg38_alleles(fasta_base, ref, alt):
    ref = ref.upper(); alt = alt.upper(); b = fasta_base.upper()
    if b == ref:
        return ref, alt, ""
    if b == alt:
        return alt, ref, "hg19_ref_alt_swapped_on_hg38"
    if b == COMP.get(ref):
        return COMP[ref], COMP[alt], "strand_flip"
    if b == COMP.get(alt):
        return COMP[alt], COMP[ref], "strand_flip_swapped"
    return b, None, f"UNRESOLVED_fasta={b}_alleles={ref}/{alt}"


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--liver-tracks", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--substrate", default=None)
    ap.add_argument("--folds", nargs="+", default=list(FOLDS),
                    help="HF model ids of the 4 replicate folds to ensemble")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # ---- load leads --------------------------------------------------------
    leads = []
    with open(args.truth) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("is_signal_lead") == "TRUE":
                leads.append(row)
    log(f"loaded {len(leads)} independent lead signals (is_signal_lead==TRUE)")

    # ---- substrate cross-check (optional) ----------------------------------
    sub = {}
    if args.substrate and os.path.exists(args.substrate):
        with open(args.substrate) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                sub.setdefault(row["variant_id_hg19"], row)
        log(f"substrate cross-check table: {len(sub)} variants")

    # ---- GTF gene models ---------------------------------------------------
    want_ens = {strip_ver(l["ensembl"]) for l in leads} | {SORT1_ANCHOR["ensembl"]}
    want_sym = {l["gene"] for l in leads} | {SORT1_ANCHOR["gene"]}
    log(f"parsing GTF for {len(want_ens)} gene ids ...")
    genes, sym2sid = parse_gtf_genes(args.gtf, want_ens, want_sym)

    # ---- tracks ------------------------------------------------------------
    liver = load_liver_tracks(args.liver_tracks)
    liver_idx = np.array([t["track_index"] for t in liver], dtype=int)
    pos_in_subset = {t["track_index"]: i for i, t in enumerate(liver)}
    rna_tracks = [t for t in liver if t["assay"] == "rna_liver"]
    acc_tracks = [t for t in liver if t["assay"] == "accessibility"]
    acc_sub = np.array([pos_in_subset[t["track_index"]] for t in acc_tracks], dtype=int)
    rna_all_sub = np.array([pos_in_subset[t["track_index"]] for t in rna_tracks], dtype=int)

    # tissue-matched GTEx-liver subset positions (verify presence)
    gtex_sub = []
    for ti in GTEX_LIVER_TRACK_INDICES:
        if ti in pos_in_subset:
            gtex_sub.append(pos_in_subset[ti])
        else:
            log(f"WARNING: GTEx liver track_index {ti} not in manifest subset")
    gtex_sub = np.array(gtex_sub, dtype=int)
    if len(gtex_sub) == 0:
        log("FATAL: no GTEx-liver tracks resolved; cannot compute tissue-matched primary")
        sys.exit(3)
    gtex_descs = [liver[j]["description"] for j in gtex_sub]
    log(f"liver tracks: {len(liver)} ({len(rna_tracks)} rna_liver / "
        f"{len(acc_tracks)} accessibility)")
    log(f"PRIMARY tissue-matched = {len(gtex_sub)} GTEx-liver tracks "
        f"(indices {[liver[j]['track_index'] for j in gtex_sub]}; "
        f"unstranded={[liver[j]['strand']=='.' for j in gtex_sub]}); desc={gtex_descs[0]!r}")

    rc_perm, n_swapped, n_unpaired = build_rc_track_perm(liver)
    log(f"RC track permutation: {n_swapped} stranded tracks swap to partner, "
        f"{n_unpaired} stranded-but-unpaired -> self, "
        f"{len(liver)-n_swapped-n_unpaired} unstranded -> self "
        f"(GTEx primary all unstranded -> self)")

    # ---- models ------------------------------------------------------------
    import torch
    import pysam
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"torch {torch.__version__} | device={dev} | cuda={torch.cuda.is_available()}")
    if dev == "cuda":
        log(f"GPU: {torch.cuda.get_device_name(0)} | "
            f"{torch.cuda.get_device_properties(0).total_memory/1e9:.1f} GB")
    from borzoi_pytorch import Borzoi
    models = []
    for m in args.folds:
        log(f"loading Borzoi fold weights: {m}")
        models.append(Borzoi.from_pretrained(m).to(dev).eval())
    n_folds = len(models)
    log(f"loaded {n_folds}-fold ensemble.")

    fasta = pysam.FastaFile(args.fasta)
    n_targets = N_HUMAN_TRACKS

    def predict_liver(model, oh):
        """Run one Borzoi fold on a (4, L) one-hot; return (bins, n_liver) sliced
        to the liver subset, in the sequence's own (input) frame."""
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
        return pred[:, liver_idx]                     # (bins, n_liver)

    def rc_to_fwd_frame(pred_rc):
        """Map an RC-input prediction back to the forward genomic frame:
        reverse the bin axis and swap stranded track pairs."""
        return pred_rc[::-1, :][:, rc_perm]

    def fwd_rc_cov(model, oh_fwd, oh_rc):
        """Forward+RC-averaged coverage for one fold, forward genomic frame."""
        pf = predict_liver(model, oh_fwd)
        pr = rc_to_fwd_frame(predict_liver(model, oh_rc))
        return 0.5 * (pf + pr)

    # ---- score -------------------------------------------------------------
    rows = []
    work = list(leads)
    lead_vids = {l["variant_id"] for l in leads}
    if SORT1_ANCHOR["variant_id"] not in lead_vids:
        work.append({**SORT1_ANCHOR, "is_signal_lead": "ANCHOR"})

    n_scored = 0
    n_unresolved_gene = 0
    n_out_of_window = 0
    n_indel = 0
    n_fasta_mismatch = 0
    sub_agree = 0
    sub_checked = 0
    sort1_logsed = None

    cols = ["variant_id", "chr", "pos_hg38", "gene", "ensembl",
            "hg38_ref", "hg38_alt",
            "borzoi_logsed_liver", "borzoi_logsed_allliver", "borzoi_acc_delta",
            "fold_sign_agreement", "borzoi_logsed_per_fold",
            "gene_tss", "gene_strand", "gene_in_window", "variant_tss_dist",
            "n_exon_bins", "n_rna_tracks_used", "gene_resolution",
            "is_indel", "is_anchor", "note"]

    for k, lead in enumerate(work):
        vid = lead["variant_id"]
        chrom_raw = lead["chr"]
        chrom = chrom_raw if chrom_raw.startswith("chr") else "chr" + chrom_raw
        pos1 = int(lead["pos_hg38"])
        v0 = pos1 - 1
        ref = lead["ref"].upper()
        alt = lead["alt"].upper()
        gene = lead["gene"]
        ensembl = lead["ensembl"]
        is_anchor = lead.get("is_signal_lead") == "ANCHOR"
        is_indel = (len(ref) != 1 or len(alt) != 1)

        rec, how = resolve_gene(genes, sym2sid, ensembl, gene, chrom)
        base_row = {c: "" for c in cols}
        base_row.update({
            "variant_id": vid, "chr": chrom_raw, "pos_hg38": pos1,
            "gene": gene, "ensembl": ensembl,
            "gene_in_window": "FALSE", "gene_resolution": how,
            "is_indel": str(is_indel).upper(), "is_anchor": str(is_anchor).upper(),
            "n_rna_tracks_used": len(gtex_sub),
        })
        if rec is None or rec.get("gene_start") is None:
            base_row["note"] = "gene_unresolved_in_gtf"
            n_unresolved_gene += 1
            rows.append(base_row)
            continue

        tss1, strand = gene_tss(rec)
        tss0 = tss1 - 1
        dist = abs(v0 - tss0)
        base_row["gene_tss"] = tss1
        base_row["gene_strand"] = strand
        base_row["variant_tss_dist"] = dist

        if dist >= IN_HALF:
            base_row["note"] = "variant_outside_input_receptive_field"
            n_out_of_window += 1
            rows.append(base_row)
            continue
        base_row["gene_in_window"] = "TRUE"

        if is_indel:
            base_row["note"] = "indel_not_scored_as_snv"
            n_indel += 1
            rows.append(base_row)
            continue

        fasta_base = fasta.fetch(chrom, v0, v0 + 1).upper()
        hg38_ref, hg38_alt, note = resolve_hg38_alleles(fasta_base, ref, alt)
        base_row["hg38_ref"] = hg38_ref
        base_row["hg38_alt"] = hg38_alt if hg38_alt else ""
        if hg38_alt is None:
            base_row["note"] = note
            n_fasta_mismatch += 1
            rows.append(base_row)
            continue
        if note:
            base_row["note"] = note

        if vid in sub:
            s = sub[vid]
            sub_checked += 1
            if (s.get("ref_hg38", "").upper() == hg38_ref and
                    s.get("alt_hg38", "").upper() == hg38_alt):
                sub_agree += 1
            else:
                base_row["note"] = (base_row["note"] + ";" if base_row["note"] else "") + \
                    f"substrate_disagree(s_ref={s.get('ref_hg38')},s_alt={s.get('alt_hg38')})"

        # TSS-centered window; substitute at variant
        window = fetch_window(fasta, chrom, tss0, SEQ_LEN)
        vi = v0 - tss0 + IN_HALF
        win_base = window[vi]
        if win_base != hg38_ref:
            n_fasta_mismatch += 1
            base_row["note"] = (base_row["note"] + ";" if base_row["note"] else "") + \
                f"window_base={win_base}!=hg38_ref={hg38_ref}"
        alt_window = window[:vi] + hg38_alt + window[vi + 1:]

        # exon-overlapping output bins
        ebins = np.array(sorted(exon_union_bins(rec, tss0)), dtype=int)
        base_row["n_exon_bins"] = len(ebins)

        # one-hots (forward + reverse-complement) for ref and alt
        ref_oh = one_hot(window)
        alt_oh = one_hot(alt_window)
        ref_oh_rc = rc_one_hot(ref_oh)
        alt_oh_rc = rc_one_hot(alt_oh)

        t0 = time.time()
        if len(ebins) == 0:
            base_row["note"] = (base_row["note"] + ";" if base_row["note"] else "") + \
                "no_exon_bins_in_output_window"
            base_row["borzoi_logsed_liver"] = float("nan")
            base_row["borzoi_logsed_allliver"] = float("nan")
            base_row["borzoi_acc_delta"] = float("nan")
            base_row["fold_sign_agreement"] = 0
            base_row["borzoi_logsed_per_fold"] = ""
            rows.append(base_row)
            continue

        # per-fold logSED per liver track (forward+RC coverage averaged within fold)
        logsed_perfold = np.empty((n_folds, len(liver_idx)), dtype=np.float64)
        acc_delta_accum = 0.0
        vbin = (v0 - (tss0 - OUT_HALF)) // BIN_BP
        acc_valid = (0 <= vbin < N_OUT_BINS and len(acc_sub) > 0)
        for fi, model in enumerate(models):
            ref_cov = fwd_rc_cov(model, ref_oh, ref_oh_rc)   # (bins, n_liver)
            alt_cov = fwd_rc_cov(model, alt_oh, alt_oh_rc)
            ref_es = np.clip(ref_cov[ebins, :].sum(axis=0), 0.0, None)  # (n_liver,)
            alt_es = np.clip(alt_cov[ebins, :].sum(axis=0), 0.0, None)
            logsed_perfold[fi] = np.log(alt_es + 1.0) - np.log(ref_es + 1.0)
            if acc_valid:
                acc_delta_accum += float(np.mean((alt_cov[vbin] - ref_cov[vbin])[acc_sub]))

        # ensemble = mean over folds of per-fold logSED (per track), then aggregate
        logsed_ens = logsed_perfold.mean(axis=0)                     # (n_liver,)
        primary = float(np.mean(logsed_ens[gtex_sub]))              # tissue-matched
        allliver = float(np.mean(logsed_ens[rna_all_sub]))          # all-58 robustness
        perfold_liver = logsed_perfold[:, gtex_sub].mean(axis=1)    # (n_folds,)
        s_primary = np.sign(primary)
        fold_agree = int(np.sum(np.sign(perfold_liver) == s_primary)) if s_primary != 0 else 0
        acc_delta = (acc_delta_accum / n_folds) if acc_valid else float("nan")

        base_row["borzoi_logsed_liver"] = primary
        base_row["borzoi_logsed_allliver"] = allliver
        base_row["borzoi_acc_delta"] = acc_delta
        base_row["fold_sign_agreement"] = fold_agree
        base_row["borzoi_logsed_per_fold"] = ",".join(f"{v:.6f}" for v in perfold_liver)

        rows.append(base_row)
        n_scored += 1
        if is_anchor:
            sort1_logsed = primary
        if k % 25 == 0 or is_anchor:
            log(f"[{k+1}/{len(work)}] {vid} {gene} {hg38_ref}>{hg38_alt} "
                f"strand={strand} tss_dist={dist} exon_bins={len(ebins)} "
                f"logSED_liver={primary:+.5f} (all58={allliver:+.5f}) "
                f"foldagree={fold_agree}/{n_folds} acc={acc_delta:+.5f} "
                f"[{time.time()-t0:.1f}s]")

    # ---- write -------------------------------------------------------------
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"WROTE {args.out} ({len(rows)} rows)")

    # ---- summary -----------------------------------------------------------
    log("================= SUMMARY =================")
    log(f"ensemble folds             : {n_folds} ({', '.join(args.folds)})")
    log(f"lead pairs (excl anchor)   : {len(leads)}")
    log(f"scored (logSED computed)   : {n_scored}")
    log(f"gene unresolved in GTF     : {n_unresolved_gene}")
    log(f"variant out of input window: {n_out_of_window}")
    log(f"indel (unscored as SNV)    : {n_indel}")
    log(f"fasta/window base warnings : {n_fasta_mismatch}")
    if sub_checked:
        log(f"substrate cross-check      : {sub_agree}/{sub_checked} hg38 alleles agree")
    if sort1_logsed is not None:
        verdict = "POSITIVE (PASS)" if sort1_logsed > 0 else "NON-POSITIVE (FAIL)"
        log(f"SORT1 ANCHOR logSED(G->T)  : {sort1_logsed:+.6f}  => {verdict}")
    else:
        log("SORT1 ANCHOR               : NOT SCORED")
    log("==========================================")

    anchor_ok = (sort1_logsed is not None and sort1_logsed > 0)
    wrote_leads = (n_scored - (1 if anchor_ok else 0)) > 0
    sys.exit(0 if (wrote_leads and anchor_ok) else 2)


if __name__ == "__main__":
    main()
