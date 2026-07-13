#!/usr/bin/env python3
"""
65_alphagenome_score.py  --  AlphaGenome variant scorer (Stage-2b seqfunc).

Independent SECOND sequence-to-function model alongside Borzoi (64_borzoi_eqtl_sed.py).
Calls the AlphaGenome API (Google DeepMind, Avsec et al. Nature 2026) to score each
Broadaway eQTL lead variant on its eGene with four heads, mirroring the Borzoi output
schema so the downstream concordance join (65b) is trivial:

  (i)   ag_gene_lfc         SIGNED gene-expression effect -- GeneMaskLFCScorer
        (RECOMMENDED_VARIANT_SCORERS['RNA_SEQ'], aggregation DIFF_LOG2_SUM =
        log2(sum ALT exon coverage) - log2(sum REF exon coverage)) averaged over
        LIVER RNA tracks. This is the AlphaGenome analog of Borzoi's logSED and the
        one compared to measured Broadaway eQTL sign.
  (ii)  ag_splice_score     splice-junction + splice-site-usage disruption
        (SpliceJunctionScorer + GeneMaskSplicingScorer) -- patches the SpliceAI
        HSD17B13 miss.
  (iii) ag_accessibility_delta  chromatin accessibility delta (CenterMaskScorer on
        ATAC; DNASE recorded separately) over liver/hepatocyte tracks.
  (iv)  ag_gene_quantile / ag_splice_quantile  quantile_score = calibration vs the
        common-variant (gnomAD MAF>0.01) background.

POLARITY CONTRACT (must match Borzoi/Broadaway exactly -- the sign-flip trap):
  Score strictly hg38 REF -> hg38 ALT (genomic reference orientation). RECORD the
  hg38_ref/hg38_alt used. Do NOT pre-orient to the eQTL effect allele -- the stats
  step (65b) orients. GeneMaskLFCScorer raw_score is ALT-REF in log2, so
  POSITIVE ag_gene_lfc => hg38ALT increases liver gene expression.

ADDITIVE ONLY -- new file (script 65). Touches no existing script; writes only under
results/seqfunc/. Reuses the Borzoi hg38 reference FASTA and the same benchmark truth.

RUN-PATH: AlphaGenome API (free non-commercial key). Compute is Google-side, so NO GPU
is needed locally -- only gRPC egress to gdmscience.googleapis.com:443. Key is read
from $ALPHA_GENOME_API_KEY or a keyfile (default ~/.alphagenome_key, chmod 600); it is
never logged and never written to the repo.
"""

import argparse
import csv
import os
import re
import sys
import time

# ----------------------------------------------------------------------------- constants
BASES = "ACGT"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

# Liver / hepatocyte ontology (UBERON:0002107 liver; CL:0000182 hepatocyte). Used to
# subset AlphaGenome tracks to the eQTL-relevant tissue, analogous to Borzoi's
# liver_tracks_selected.csv.
LIVER_ONTOLOGY = {"UBERON:0002107", "CL:0000182"}
LIVER_NAME_RE = re.compile(r"\bliver\b|hepatocyte", re.IGNORECASE)

# Positive-control + disease anchors (appended regardless of the leads). hg38 coords
# verified against the reference FASTA before use.
#   SORT1 rs12740374 chr1:109274968 G>T -- minor T CREATES a C/EBP site and INCREASES
#     hepatic SORT1 (Musunuru 2010) => ag_gene_lfc(G->T) MUST be POSITIVE.
#   PNPLA3 rs738409 chr22:43928847 C>G -- I148M, the canonical MASLD risk allele.
#   HSD17B13 rs72613567 chr4:87310240 T>TA -- splice-donor insertion (dbSNP
#     g.87310241dup; FASTA base@87310240=T, @87310241=A). SpliceAI underscored it;
#     the splice heads are tested here. AlphaGenome scores indels (Borzoi could not).
ANCHORS = [
    {"variant_id": "SORT1_rs12740374", "chr": "1", "pos_hg38": 109274968,
     "ref": "G", "alt": "T", "gene": "SORT1", "ensembl": "ENSG00000134243",
     "anchor_name": "SORT1_posctrl"},
    {"variant_id": "PNPLA3_rs738409", "chr": "22", "pos_hg38": 43928847,
     "ref": "C", "alt": "G", "gene": "PNPLA3", "ensembl": "ENSG00000100344",
     "anchor_name": "PNPLA3_I148M"},
    {"variant_id": "HSD17B13_rs72613567", "chr": "4", "pos_hg38": 87310240,
     "ref": "T", "alt": "TA", "gene": "HSD17B13", "ensembl": "ENSG00000170509",
     "anchor_name": "HSD17B13_splice"},
]

OUT_COLS = [
    "variant_id", "chr", "pos_hg38", "gene", "ensembl", "hg38_ref", "hg38_alt",
    "ag_gene_lfc", "ag_gene_quantile", "ag_splice_score", "ag_accessibility_delta",
    "ag_gene_lfc_alltracks", "ag_splice_quantile", "ag_dnase_delta",
    "n_liver_rna_tracks", "n_splice_rows", "gene_matched",
    "is_indel", "is_anchor", "anchor_name", "note",
]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def strip_ver(ens):
    return ens.split(".")[0] if ens else ens


# ------------------------------------------------------------------ hg38 allele resolve
def resolve_hg38_alleles(fasta_base, ref, alt):
    """SNV hg38 orientation, identical logic to 64_borzoi_eqtl_sed.py.
    Returns (hg38_ref, hg38_alt, note). fasta_base is the reference ground truth."""
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


# ------------------------------------------------------------------ key handling
def load_api_key(args):
    k = os.environ.get("ALPHA_GENOME_API_KEY")
    if k and k.strip():
        return k.strip(), "env:ALPHA_GENOME_API_KEY"
    kf = args.key_file or os.path.expanduser("~/.alphagenome_key")
    if os.path.exists(kf):
        with open(kf) as f:
            v = f.read().strip()
        if v:
            return v, f"keyfile:{kf}"
    return None, None


# ------------------------------------------------------------------ tidy-score helpers
def _liver_mask(df):
    """Boolean Series: track is liver/hepatocyte by ontology, gtex tissue, or name."""
    import pandas as pd  # noqa
    m = df["ontology_curie"].isin(LIVER_ONTOLOGY) if "ontology_curie" in df else False
    for col in ("gtex_tissue", "biosample_name", "track_name"):
        if col in df.columns:
            hit = df[col].astype(str).str.contains(LIVER_NAME_RE)
            m = hit if m is False else (m | hit)
    return m


def _gene_mask(df, ensembl, symbol):
    """Boolean Series: row corresponds to the target eGene (versioned-id or symbol)."""
    sid = strip_ver(ensembl)
    m = False
    if "gene_id" in df.columns:
        m = df["gene_id"].astype(str).map(strip_ver) == sid
    if "gene_name" in df.columns:
        hit = df["gene_name"].astype(str) == symbol
        m = hit if m is False else (m | hit)
    return m


def _signed_extreme(series):
    """Return the signed value with the largest absolute magnitude (NaN-safe)."""
    import numpy as np
    s = series.dropna()
    if len(s) == 0:
        return float("nan")
    return float(s.iloc[s.abs().values.argmax()])


def extract_scores(df, ensembl, symbol, dbg=False):
    """Flatten one variant's tidy_scores DataFrame into the output-row numeric fields."""
    import numpy as np
    out = dict(ag_gene_lfc=float("nan"), ag_gene_quantile=float("nan"),
               ag_splice_score=float("nan"), ag_splice_quantile=float("nan"),
               ag_accessibility_delta=float("nan"), ag_dnase_delta=float("nan"),
               ag_gene_lfc_alltracks=float("nan"),
               n_liver_rna_tracks=0, n_splice_rows=0, gene_matched="FALSE", note="")
    if df is None or len(df) == 0:
        out["note"] = "empty_tidy_scores"
        return out

    if dbg:
        log(f"  tidy cols: {list(df.columns)}")
        if "output_type" in df:
            log(f"  output_types: {sorted(df['output_type'].astype(str).unique())}")

    ot = df["output_type"].astype(str) if "output_type" in df else None
    gmask = _gene_mask(df, ensembl, symbol)
    lmask = _liver_mask(df)

    # ---- (i) gene-expression LFC (RNA_SEQ, GeneMaskLFCScorer) ------------------
    if ot is not None:
        rna = df[ot.str.contains("RNA_SEQ", case=False, na=False)]
        if len(rna):
            rna_gene = rna[_gene_mask(rna, ensembl, symbol)]
            out["gene_matched"] = "TRUE" if len(rna_gene) else "FALSE"
            # all-tissue robustness aggregate (mirrors borzoi_sed_rna_alltracks)
            src_all = rna_gene if len(rna_gene) else rna
            out["ag_gene_lfc_alltracks"] = float(np.nanmean(
                np.asarray(src_all["raw_score"], dtype=float)))
            rna_liver = rna_gene[_liver_mask(rna_gene)] if len(rna_gene) else \
                rna[_liver_mask(rna)]
            out["n_liver_rna_tracks"] = int(len(rna_liver))
            if len(rna_liver):
                out["ag_gene_lfc"] = float(np.nanmean(
                    np.asarray(rna_liver["raw_score"], dtype=float)))
                if "quantile_score" in rna_liver:
                    out["ag_gene_quantile"] = float(np.nanmean(
                        np.asarray(rna_liver["quantile_score"], dtype=float)))
            else:
                # no liver track for this gene -> fall back to all-tissue mean, flagged
                out["ag_gene_lfc"] = out["ag_gene_lfc_alltracks"]
                out["note"] = "liver_fallback_alltracks"

        # ---- (ii) splice (junctions + site usage) -----------------------------
        spl = df[ot.str.contains("SPLICE", case=False, na=False)]
        out["n_splice_rows"] = int(len(spl))
        if len(spl):
            spl_gene = spl[_gene_mask(spl, ensembl, symbol)]
            src = spl_gene if len(spl_gene) else spl          # junctions may be gene-agnostic
            out["ag_splice_score"] = _signed_extreme(src["raw_score"].astype(float))
            if "quantile_score" in src:
                out["ag_splice_quantile"] = _signed_extreme(src["quantile_score"].astype(float))

        # ---- (iii) accessibility (ATAC / DNASE, liver) ------------------------
        for otype, key in (("ATAC", "ag_accessibility_delta"), ("DNASE", "ag_dnase_delta")):
            sub = df[ot.str.fullmatch(otype, case=False, na=False)]
            if len(sub):
                sub_l = sub[_liver_mask(sub)]
                use = sub_l if len(sub_l) else sub
                out[key] = float(np.nanmean(np.asarray(use["raw_score"], dtype=float)))
    return out


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True,
                    help="broadaway_benchmark_truth.tsv (is_signal_lead==TRUE)")
    ap.add_argument("--fasta", required=True, help="hg38 genome.fa (SNV hg38 orient)")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--borzoi-scores", default=None,
                    help="borzoi_eqtl_scores.tsv (optional hg38-allele cross-check)")
    ap.add_argument("--key-file", default=None,
                    help="API key file (default ~/.alphagenome_key)")
    ap.add_argument("--seq-length", default="SEQUENCE_LENGTH_1MB",
                    choices=["SEQUENCE_LENGTH_16KB", "SEQUENCE_LENGTH_100KB",
                             "SEQUENCE_LENGTH_500KB", "SEQUENCE_LENGTH_1MB"])
    ap.add_argument("--limit", type=int, default=0, help="score only first N leads (smoke)")
    ap.add_argument("--anchors-only", action="store_true",
                    help="score only the 3 anchors (smoke)")
    ap.add_argument("--sleep", type=float, default=0.2, help="seconds between API calls")
    ap.add_argument("--max-retries", type=int, default=5)
    ap.add_argument("--dry-run", action="store_true",
                    help="validate client + inputs, do NOT call the API")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # ---- key --------------------------------------------------------------
    api_key, key_src = load_api_key(args)
    if api_key is None and not args.dry_run:
        log("=" * 62)
        log("READY -- awaiting user API key. No scoring performed.")
        log("Provide a free non-commercial AlphaGenome key one of two ways:")
        log("  export ALPHA_GENOME_API_KEY=<key>       # env var, OR")
        log("  printf %s <key> > ~/.alphagenome_key && chmod 600 ~/.alphagenome_key")
        log("Register at https://deepmind.google.com/science/alphagenome/ (accept ToS).")
        log("Then re-run this script / resubmit the sbatch. Key is never logged or committed.")
        log("=" * 62)
        sys.exit(3)
    if api_key is not None:
        log(f"API key loaded from {key_src} (len={len(api_key)}; value not shown)")

    # ---- imports (fail loudly if env wrong) -------------------------------
    import numpy as np       # noqa
    import pandas as pd
    import pysam
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
    import grpc

    seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS[args.seq_length]
    R = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    # four heads: signed gene LFC (RNA_SEQ) + splice (junctions + site usage) + accessibility
    scorer_keys = ["RNA_SEQ", "SPLICE_JUNCTIONS", "SPLICE_SITE_USAGE", "SPLICE_SITES",
                   "ATAC", "DNASE"]
    selected = [R[k] for k in scorer_keys if k in R]
    log(f"seq_length={args.seq_length} ({seq_len:,} bp); scorers: "
        + ", ".join(type(s).__name__ for s in selected))

    # ---- leads + anchors --------------------------------------------------
    leads = []
    if not args.anchors_only:
        with open(args.truth) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("is_signal_lead") == "TRUE":
                    leads.append({
                        "variant_id": row["variant_id"], "chr": row["chr"],
                        "pos_hg38": int(row["pos_hg38"]), "ref": row["ref"],
                        "alt": row["alt"], "gene": row["gene"],
                        "ensembl": row["ensembl"], "anchor_name": "",
                    })
        if args.limit and args.limit > 0:
            leads = leads[:args.limit]
    lead_vids = {l["variant_id"] for l in leads}
    work = list(leads)
    for a in ANCHORS:
        if a["variant_id"] not in lead_vids:
            work.append(dict(a))
    log(f"scoring {len(work)} variants ({len(leads)} leads + "
        f"{len(work)-len(leads)} anchors)")

    # ---- borzoi hg38-allele cross-check table -----------------------------
    bz = {}
    if args.borzoi_scores and os.path.exists(args.borzoi_scores):
        with open(args.borzoi_scores) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                bz[r["variant_id"] + "|" + strip_ver(r.get("ensembl", ""))] = r
        log(f"borzoi cross-check table: {len(bz)} rows")

    fasta = pysam.FastaFile(args.fasta)

    if args.dry_run:
        # Build every Variant/Interval to validate coords & alleles; no API call.
        n_ok = n_bad = 0
        for lead in work:
            chrom = lead["chr"] if str(lead["chr"]).startswith("chr") else "chr" + str(lead["chr"])
            pos1 = int(lead["pos_hg38"]); v0 = pos1 - 1
            ref, alt = lead["ref"].upper(), lead["alt"].upper()
            is_indel = (len(ref) != 1 or len(alt) != 1)
            try:
                fbase = fasta.fetch(chrom, v0, v0 + 1).upper()
                if not is_indel:
                    hg38_ref, hg38_alt, _ = resolve_hg38_alleles(fbase, ref, alt)
                else:
                    hg38_ref, hg38_alt = ref, alt
                v = genome.Variant(chromosome=chrom, position=pos1,
                                   reference_bases=hg38_ref,
                                   alternate_bases=(hg38_alt or "N"),
                                   name=lead["variant_id"])
                _ = v.reference_interval.resize(seq_len)
                n_ok += 1
            except Exception as e:  # noqa
                n_bad += 1
                log(f"  DRY-BAD {lead['variant_id']}: {type(e).__name__}: {e}")
        log(f"DRY-RUN: {n_ok} variants build OK, {n_bad} failed. Client + inputs validated. "
            f"No API call made.")
        sys.exit(0 if n_bad == 0 else 2)

    # ---- model ------------------------------------------------------------
    model = dna_client.create(api_key)
    log("dna_client created; scoring...")

    def score_one(chrom, pos1, hg38_ref, hg38_alt, vid):
        variant = genome.Variant(chromosome=chrom, position=pos1,
                                 reference_bases=hg38_ref, alternate_bases=hg38_alt,
                                 name=vid)
        interval = variant.reference_interval.resize(seq_len)
        last = None
        for attempt in range(args.max_retries):
            try:
                scores = model.score_variant(interval=interval, variant=variant,
                                             variant_scorers=selected,
                                             organism=dna_client.Organism.HOMO_SAPIENS)
                return variant_scorers.tidy_scores([scores])
            except grpc.RpcError as e:  # transient -> backoff & retry
                last = e
                code = e.code()
                transient = code in (grpc.StatusCode.UNAVAILABLE,
                                     grpc.StatusCode.RESOURCE_EXHAUSTED,
                                     grpc.StatusCode.DEADLINE_EXCEEDED,
                                     grpc.StatusCode.INTERNAL)
                if not transient:
                    raise
                wait = min(30.0, 2.0 ** attempt)
                log(f"  transient {code} (attempt {attempt+1}/{args.max_retries}); "
                    f"backoff {wait:.0f}s")
                time.sleep(wait)
        raise last

    rows = []
    n_scored = n_err = n_fasta_mismatch = 0
    sort1_lfc = hsd_splice_q = None
    dbg_done = False

    for k, lead in enumerate(work):
        vid = lead["variant_id"]
        chrom_raw = str(lead["chr"])
        chrom = chrom_raw if chrom_raw.startswith("chr") else "chr" + chrom_raw
        pos1 = int(lead["pos_hg38"]); v0 = pos1 - 1
        ref, alt = lead["ref"].upper(), lead["alt"].upper()
        gene, ensembl = lead["gene"], lead["ensembl"]
        is_anchor = bool(lead.get("anchor_name"))
        is_indel = (len(ref) != 1 or len(alt) != 1)

        base = {c: "" for c in OUT_COLS}
        base.update(variant_id=vid, chr=chrom_raw, pos_hg38=pos1, gene=gene,
                    ensembl=ensembl, is_indel=str(is_indel).upper(),
                    is_anchor=str(is_anchor).upper(),
                    anchor_name=lead.get("anchor_name", ""), gene_matched="FALSE",
                    n_liver_rna_tracks=0, n_splice_rows=0)

        fbase = fasta.fetch(chrom, v0, v0 + 1).upper()
        note = ""
        if not is_indel:
            hg38_ref, hg38_alt, note = resolve_hg38_alleles(fbase, ref, alt)
            if hg38_alt is None:
                base["hg38_ref"] = hg38_ref; base["note"] = note
                n_fasta_mismatch += 1; rows.append(base); continue
        else:
            hg38_ref, hg38_alt = ref, alt
            if fbase != ref[0].upper():
                note = f"indel_fasta_base={fbase}!=ref0={ref[0]}"
        base["hg38_ref"] = hg38_ref; base["hg38_alt"] = hg38_alt

        # borzoi hg38-allele cross-check
        bkey = vid + "|" + strip_ver(ensembl)
        if bkey in bz:
            br, ba_ = bz[bkey].get("hg38_ref", ""), bz[bkey].get("hg38_alt", "")
            if br and ba_ and (br.upper() != hg38_ref or ba_.upper() != hg38_alt):
                note = (note + ";" if note else "") + \
                    f"borzoi_allele_disagree(bz={br}/{ba_})"

        try:
            df = score_one(chrom, pos1, hg38_ref, hg38_alt, vid)
            ext = extract_scores(df, ensembl, gene, dbg=(not dbg_done))
            dbg_done = True
            for kk in ("ag_gene_lfc", "ag_gene_quantile", "ag_splice_score",
                       "ag_accessibility_delta", "ag_gene_lfc_alltracks",
                       "ag_splice_quantile", "ag_dnase_delta",
                       "n_liver_rna_tracks", "n_splice_rows", "gene_matched"):
                base[kk] = ext[kk]
            if ext["note"]:
                note = (note + ";" if note else "") + ext["note"]
            n_scored += 1
        except Exception as e:  # noqa
            note = (note + ";" if note else "") + f"score_error:{type(e).__name__}:{e}"
            n_err += 1

        base["note"] = note
        rows.append(base)

        if is_anchor and lead.get("anchor_name") == "SORT1_posctrl":
            sort1_lfc = base["ag_gene_lfc"]
        if is_anchor and lead.get("anchor_name") == "HSD17B13_splice":
            hsd_splice_q = base["ag_splice_quantile"]

        if is_anchor or k % 20 == 0:
            log(f"[{k+1}/{len(work)}] {vid} {gene} {hg38_ref}>{hg38_alt} "
                f"lfc={base['ag_gene_lfc']} q={base['ag_gene_quantile']} "
                f"splice={base['ag_splice_score']} atac={base['ag_accessibility_delta']} "
                f"{'['+base['note']+']' if base['note'] else ''}")
        time.sleep(args.sleep)

    # ---- write ------------------------------------------------------------
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"WROTE {args.out} ({len(rows)} rows)")

    # ---- summary ----------------------------------------------------------
    log("================= SUMMARY =================")
    log(f"variants total             : {len(work)} ({len(leads)} leads)")
    log(f"scored (API returned)      : {n_scored}")
    log(f"scoring errors             : {n_err}")
    log(f"fasta-mismatch (unscored)  : {n_fasta_mismatch}")

    def _isnum(x):
        try:
            return x != "" and x == x and float(x) == float(x)
        except Exception:
            return False

    if sort1_lfc is not None and _isnum(sort1_lfc):
        v = float(sort1_lfc)
        log(f"SORT1 ANCHOR ag_gene_lfc   : {v:+.5f}  => "
            f"{'POSITIVE (PASS)' if v > 0 else 'NON-POSITIVE (FAIL)'}")
    else:
        log(f"SORT1 ANCHOR ag_gene_lfc   : {sort1_lfc} (not numeric)")
    if hsd_splice_q is not None:
        log(f"HSD17B13 splice quantile   : {hsd_splice_q}  (nonzero => AG flags the splice variant)")
    log("==========================================")

    anchor_ok = (sort1_lfc is not None and _isnum(sort1_lfc) and float(sort1_lfc) > 0)
    sys.exit(0 if (n_scored > 0 and anchor_ok) else 2)


if __name__ == "__main__":
    main()
