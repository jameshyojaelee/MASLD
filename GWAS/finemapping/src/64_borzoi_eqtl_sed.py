#!/usr/bin/env python3
"""
64_borzoi_eqtl_sed.py  --  Borzoi GENE-LEVEL SED eQTL scorer (Stage-2a seqfunc).

Scores each independent Broadaway lead (variant x gene) signal with Borzoi by
computing the signed change in predicted liver RNA coverage summed over the
gene's exonic output bins (SED = Sum of Expression Differences), the correct
eQTL-scoring quantity -- NOT the variant-center delta the smoke test (62) used.

Reuses 62_borzoi_score.py machinery (one_hot / fetch_window / find_axes /
predict). ADDITIVE ONLY -- new file (script number 64); touches no existing
script and writes only under results/seqfunc/.

METHOD (per (variant, gene) lead pair):
  a. Gene exon coords + TSS + strand from gencode v49 GTF (match ensembl,
     version-stripped; fall back to gene symbol). TSS = MANE_Select transcript
     TSS if resolvable, else gene-boundary TSS (start if +, end if -).
     Exon set = union of all exons across all of the gene's transcripts.
  b. CENTER the 524,288 bp input window on the gene TSS (not the variant) so the
     gene body lands in the central 196,608 bp OUTPUT window. Require
     |variant_pos - TSS| < 262,144 (input half-width); else gene_in_window=FALSE
     and the pair is left unscored.
  c. REF = hg38 reference (FASTA). ALT = hg38 non-reference allele substituted at
     the variant. POLARITY: strictly hg38REF -> hg38ALT (genomic reference
     orientation); hg38_ref/hg38_alt recorded. The stats stage orients to the
     eQTL effect allele -- the scorer does a clean genomic ref->alt SED only.
  d. SED = sum over OUTPUT BINS overlapping the gene's exons of
     (alt_cov - ref_cov), averaged over the sense-strand rna_liver tracks.
     bin i covers [TSS - 98304 + i*32, +32). Positive SED => hg38ALT increases
     predicted gene liver expression. Secondary: accessibility delta at the
     VARIANT's output bin (hepatocyte/HepG2 ATAC+DNase).

Borzoi geometry (fixed): input 524,288 bp; output 6,144 bins x 32 bp = central
196,608 bp; 7,611 human tracks; model johahi/borzoi-replicate-0.
"""

import argparse
import csv
import gzip
import os
import sys
import time

import numpy as np

# ---- Borzoi geometry (fixed by the trained model) --------------------------
SEQ_LEN = 524_288           # model input receptive field
N_HUMAN_TRACKS = 7_611      # human head
N_OUT_BINS = 6_144          # output bins
BIN_BP = 32                 # bp per output bin
OUT_SPAN = N_OUT_BINS * BIN_BP          # 196,608 bp central output window
OUT_HALF = OUT_SPAN // 2                # 98,304 bp
IN_HALF = SEQ_LEN // 2                  # 262,144 bp
BASES = "ACGT"
BASE_IDX = {b: i for i, b in enumerate(BASES)}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

# Positive-control anchor (added as an extra scored row regardless of the leads).
#   SORT1 rs12740374 chr1:109274968 hg38 G>T -- minor T allele CREATES a C/EBP
#   site and INCREASES hepatic SORT1 expression (Musunuru et al. Nature 2010).
#   => Borzoi SED(G->T) on SORT1 RNA tracks MUST be POSITIVE.
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


# ---- reused from 62_borzoi_score.py ---------------------------------------
def one_hot(seq):
    """seq (str, len SEQ_LEN) -> float32 (4, SEQ_LEN), ACGT order; N -> zeros."""
    arr = np.zeros((4, len(seq)), dtype=np.float32)
    idx = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    for b, i in BASE_IDX.items():
        arr[i] = (idx == ord(b)).astype(np.float32)
    return arr


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
    """Decide which non-batch axis is tracks (==n_targets) vs bins."""
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
    """liver_tracks_selected.csv -> list of dicts (with strand suffix parsed)."""
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


# ---- GTF gene model --------------------------------------------------------
def _attr(attrs, key):
    """Extract value of `key "value";` from a GTF attribute string, or None."""
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
    """Single streaming pass. For every gene whose version-stripped gene_id is in
    want_ens (or whose gene_name is in want_sym), collect per (stripped_id, chrom):
      - gene record: strand, gene_start(1b), gene_end(1b)
      - MANE_Select transcript TSS (1-based genomic), if present
      - union of all exon intervals (1-based inclusive)
    Returns dict keyed by (stripped_id, chrom) and a symbol->set(stripped_id) map.
    """
    want_ens = set(want_ens)
    want_sym = set(want_sym)
    genes = {}   # (sid, chrom) -> record
    sym2sid = {}  # gene_name -> set(sid)

    opener = gzip.open if gtf_path.endswith(".gz") else open
    n_lines = 0
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line and line[0] == "#":
                continue
            n_lines += 1
            # cheap prefilter: only parse lines that mention a wanted ENSG or symbol
            # (attribute field carries gene_id/gene_name). Full parse below.
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
                tags = attrs  # tag "MANE_Select"
                if 'tag "MANE_Select"' in tags:
                    rec["mane_tss"] = start if strand == "+" else end
            elif feat == "exon":
                rec["exons"].append((start, end))
    log(f"GTF: parsed {n_lines:,} feature lines; matched "
        f"{len(genes)} (gene,chrom) records for "
        f"{len({k[0] for k in genes})} distinct gene ids")
    return genes, sym2sid


def resolve_gene(genes, sym2sid, ensembl, symbol, want_chrom):
    """Pick the gene record for this lead. Prefer ensembl match on the expected
    chromosome; else ensembl on any chrom; else symbol on expected chrom."""
    sid = strip_ver(ensembl)
    # ensembl on expected chrom
    rec = genes.get((sid, want_chrom))
    if rec is not None:
        return rec, "ensembl_chrom"
    # ensembl on any chrom
    cands = [genes[k] for k in genes if k[0] == sid]
    if cands:
        # prefer a primary chromosome (chr1..chrM) over patches/scaffolds
        cands.sort(key=lambda r: (0 if r["chrom"] == want_chrom else
                                  1 if r["chrom"].startswith("chr") and "_" not in r["chrom"]
                                  else 2))
        return cands[0], "ensembl_anychrom"
    # symbol fallback
    for s in sym2sid.get(symbol, ()):
        rec = genes.get((s, want_chrom))
        if rec is not None:
            return rec, "symbol_chrom"
    return None, "unresolved"


def gene_tss(rec):
    """Return (tss_1based, strand). MANE if present, else gene-boundary."""
    strand = rec["strand"]
    if rec["mane_tss"] is not None:
        return rec["mane_tss"], strand
    if strand == "+":
        return rec["gene_start"], strand
    return rec["gene_end"], strand


def exon_union_bins(rec, tss0):
    """Set of output-bin indices [0,6144) overlapping the union of the gene's
    exons. tss0 = 0-based TSS. bin i covers [tss0 - OUT_HALF + i*32, +32)."""
    out_start = tss0 - OUT_HALF          # 0-based genomic coord of bin 0 left edge
    bins = set()
    for (s1, e1) in rec["exons"]:
        es0 = s1 - 1                     # 1-based inclusive -> 0-based
        ee0 = e1                         # 0-based half-open end
        first = (es0 - out_start) // BIN_BP
        last = (ee0 - 1 - out_start) // BIN_BP
        if last < 0 or first >= N_OUT_BINS:
            continue
        first = max(0, first)
        last = min(N_OUT_BINS - 1, last)
        bins.update(range(first, last + 1))
    return bins


def resolve_hg38_alleles(fasta_base, ref, alt):
    """Return (hg38_ref, hg38_alt, note). fasta_base is ground truth reference."""
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
    ap.add_argument("--truth", required=True,
                    help="broadaway_benchmark_truth.tsv (uses is_signal_lead==TRUE)")
    ap.add_argument("--substrate", default=None,
                    help="variant_substrate_hg38.tsv (optional cross-check)")
    ap.add_argument("--model", default="johahi/borzoi-replicate-0")
    ap.add_argument("--out", required=True, help="output TSV path")
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
    # position of each liver track within the subset array
    pos_in_subset = {t["track_index"]: i for i, t in enumerate(liver)}
    rna_tracks = [t for t in liver if t["assay"] == "rna_liver"]
    acc_tracks = [t for t in liver if t["assay"] == "accessibility"]
    acc_sub = np.array([pos_in_subset[t["track_index"]] for t in acc_tracks], dtype=int)
    log(f"liver tracks: {len(liver)} ({len(rna_tracks)} rna_liver / "
        f"{len(acc_tracks)} accessibility); "
        f"rna strand split: +={sum(t['strand']=='+' for t in rna_tracks)} "
        f"-={sum(t['strand']=='-' for t in rna_tracks)} "
        f".={sum(t['strand']=='.' for t in rna_tracks)}")

    def sense_rna_subset(strand):
        """Indices (into liver subset) of sense-strand rna tracks for `strand`.
        sense = matching-strand tracks + unstranded (total-RNA) tracks."""
        keep = [t for t in rna_tracks if t["strand"] == strand or t["strand"] == "."]
        if not keep:  # degenerate: fall back to all rna tracks
            keep = rna_tracks
        return np.array([pos_in_subset[t["track_index"]] for t in keep], dtype=int)

    # ---- model -------------------------------------------------------------
    import torch
    import pysam
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"torch {torch.__version__} | device={dev} | cuda={torch.cuda.is_available()}")
    if dev == "cuda":
        log(f"GPU: {torch.cuda.get_device_name(0)} | "
            f"{torch.cuda.get_device_properties(0).total_memory/1e9:.1f} GB")
    from borzoi_pytorch import Borzoi
    log(f"loading Borzoi weights: {args.model}")
    model = Borzoi.from_pretrained(args.model).to(dev).eval()
    log("model loaded.")

    fasta = pysam.FastaFile(args.fasta)
    n_targets = N_HUMAN_TRACKS

    def predict_liver(seq_str):
        """Run Borzoi; return (n_bins, n_liver_tracks) sliced to liver tracks."""
        oh = one_hot(seq_str)
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
        pred = out.numpy() if mode == "bins_first" else out.numpy().T  # (bins,tracks)
        return pred[:, liver_idx]                                      # (bins, 65)

    # ---- score -------------------------------------------------------------
    rows = []
    work = list(leads)
    # append SORT1 anchor as an extra row (flagged) unless already a lead
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
    sort1_sed = None

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
        base_row = {
            "variant_id": vid, "chr": chrom_raw, "pos_hg38": pos1,
            "gene": gene, "ensembl": ensembl,
            "hg38_ref": "", "hg38_alt": "",
            "borzoi_sed_rna": "", "borzoi_acc_delta": "",
            "gene_tss": "", "gene_strand": "",
            "gene_in_window": "FALSE", "variant_tss_dist": "",
            "n_exon_bins": "", "n_rna_tracks_used": "",
            "borzoi_sed_rna_alltracks": "", "gene_resolution": how,
            "is_indel": str(is_indel).upper(), "is_anchor": str(is_anchor).upper(),
            "note": "",
        }
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

        # hg38 ref/alt from FASTA ground truth
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

        # substrate cross-check
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
        vi = v0 - tss0 + IN_HALF           # index of variant within the window
        win_base = window[vi]
        if win_base != hg38_ref:
            # boundary N-pad or mismatch; record but still attempt
            n_fasta_mismatch += 1
            base_row["note"] = (base_row["note"] + ";" if base_row["note"] else "") + \
                f"window_base={win_base}!=hg38_ref={hg38_ref}"
        alt_window = window[:vi] + hg38_alt + window[vi + 1:]

        t0 = time.time()
        ref_pred = predict_liver(window)       # (bins, 65)
        alt_pred = predict_liver(alt_window)
        dpred = alt_pred - ref_pred            # (bins, 65)

        # exon-overlapping output bins
        ebins = exon_union_bins(rec, tss0)
        ebins_arr = np.array(sorted(ebins), dtype=int)
        base_row["n_exon_bins"] = len(ebins_arr)

        rna_sense = sense_rna_subset(strand)
        rna_all = np.array([pos_in_subset[t["track_index"]] for t in rna_tracks], dtype=int)
        base_row["n_rna_tracks_used"] = len(rna_sense)

        if len(ebins_arr) == 0:
            base_row["note"] = (base_row["note"] + ";" if base_row["note"] else "") + \
                "no_exon_bins_in_output_window"
            sed_sense = float("nan")
            sed_all = float("nan")
        else:
            # per-track SED = sum over exon bins of (alt-ref); then mean over tracks
            per_track = dpred[ebins_arr, :].sum(axis=0)      # (65,)
            sed_sense = float(np.mean(per_track[rna_sense]))
            sed_all = float(np.mean(per_track[rna_all]))
        base_row["borzoi_sed_rna"] = sed_sense
        base_row["borzoi_sed_rna_alltracks"] = sed_all

        # accessibility delta at the variant's output bin
        vbin = (v0 - (tss0 - OUT_HALF)) // BIN_BP
        if 0 <= vbin < N_OUT_BINS and len(acc_sub) > 0:
            acc_delta = float(np.mean(dpred[vbin, acc_sub]))
        else:
            acc_delta = float("nan")
        base_row["borzoi_acc_delta"] = acc_delta

        rows.append(base_row)
        n_scored += 1
        if is_anchor:
            sort1_sed = sed_sense
        if k % 25 == 0 or is_anchor:
            log(f"[{k+1}/{len(work)}] {vid} {gene} {hg38_ref}>{hg38_alt} "
                f"strand={strand} tss_dist={dist} exon_bins={len(ebins_arr)} "
                f"SED_rna={sed_sense:+.5f} (all={sed_all:+.5f}) acc={acc_delta:+.5f} "
                f"[{time.time()-t0:.1f}s]")

    # ---- write -------------------------------------------------------------
    cols = ["variant_id", "chr", "pos_hg38", "gene", "ensembl",
            "hg38_ref", "hg38_alt", "borzoi_sed_rna", "borzoi_acc_delta",
            "gene_tss", "gene_strand", "gene_in_window", "variant_tss_dist",
            "n_exon_bins", "n_rna_tracks_used", "borzoi_sed_rna_alltracks",
            "gene_resolution", "is_indel", "is_anchor", "note"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"WROTE {args.out} ({len(rows)} rows)")

    # ---- summary -----------------------------------------------------------
    log("================= SUMMARY =================")
    log(f"lead pairs (excl anchor)   : {len(leads)}")
    log(f"scored (SED computed)      : {n_scored}")
    log(f"gene unresolved in GTF     : {n_unresolved_gene}")
    log(f"variant out of input window: {n_out_of_window}")
    log(f"indel (unscored as SNV)    : {n_indel}")
    log(f"fasta/window base warnings : {n_fasta_mismatch}")
    if sub_checked:
        log(f"substrate cross-check      : {sub_agree}/{sub_checked} hg38 alleles agree")
    if sort1_sed is not None:
        verdict = "POSITIVE (PASS)" if sort1_sed > 0 else "NON-POSITIVE (FAIL)"
        log(f"SORT1 ANCHOR SED(G->T)     : {sort1_sed:+.6f}  => {verdict}")
    else:
        log("SORT1 ANCHOR               : NOT SCORED")
    log("==========================================")

    anchor_ok = (sort1_sed is not None and sort1_sed > 0)
    wrote_leads = (n_scored - (1 if anchor_ok else 0)) > 0
    sys.exit(0 if (wrote_leads and anchor_ok) else 2)


if __name__ == "__main__":
    main()
