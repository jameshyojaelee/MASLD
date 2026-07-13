#!/usr/bin/env bash
#SBATCH --job-name=cas13_offtarget_qc
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=cas13_offtarget_qc.%j.out
#SBATCH --error=cas13_offtarget_qc.%j.err

module load bowtie/1.3.1-GCC-11.3.0

echo "Bowtie path:"
which bowtie
bowtie --version

echo "Bowtie-build path:"
which bowtie-build

set -euo pipefail

usage() {
  cat <<'USAGE'
Cas13 guide transcriptome off-target QC using exhaustive Bowtie1 short-sequence alignment.

This script aligns target-like guide sequences to a transcriptome FASTA and summarizes
intended hits and off-target hits at 0/1/2/3 mismatches.

Recommended for your current library:
  - Use the same mouse transcriptome annotation version used for guide design.
  - Provide a transcriptome FASTA, not a genome FASTA.
  - Provide the matching GTF if possible for transcript -> gene mapping.
  - The script uses target_seq if present. Otherwise it uses reverse-complement(guide_seq),
    because your library's guide_seq is the spacer/complement and target_seq is transcript-like.

Required dependencies:
  - bowtie 1.x
  - bowtie-build
  - python3
Optional dependency:
  - Rscript + tidyverse for PNG plots

Basic run:
  bash cas13_transcriptome_offtarget_qc.sh \
    --library library.csv \
    --transcriptome-fasta Mus_musculus.transcripts.fa.gz \
    --gtf Mus_musculus.annotation.gtf.gz \
    --outdir cas13_offtarget_qc \
    --threads 8

SLURM run:
  sbatch cas13_transcriptome_offtarget_qc.sh \
    --library /path/to/library.csv \
    --transcriptome-fasta /path/to/transcripts.fa.gz \
    --gtf /path/to/annotation.gtf.gz \
    --outdir /path/to/cas13_offtarget_qc \
    --threads 8

Options:
  --library FILE              Guide library CSV. Required.
  --transcriptome-fasta FILE  Transcriptome FASTA/cDNA FASTA, optionally .gz. Required.
  --gtf FILE                  Matching GTF/GFF-style annotation, optionally .gz. Optional but recommended.
  --outdir DIR                Output directory. Default: cas13_transcriptome_offtarget_qc
  --threads N                 Bowtie threads. Default: 8
  --mismatches N              Max mismatches. Bowtie1 -v supports 0-3. Default: 3
  --query-mode MODE           auto|target_seq|guide_seq|guide_rc. Default: auto
                              auto uses target_seq if available, else reverse-complement(guide_seq).
  --both-strands              Align both query orientations. Default is forward-only with --norc,
                              which is appropriate when query sequence is target-like and the
                              transcriptome FASTA contains sense transcript sequences.
  --keep-sam                  Keep uncompressed SAM after parsing. Default: remove SAM.
  --skip-index                Reuse an existing Bowtie index in outdir/bowtie_index if present.
  --no-plots                  Do not try to make optional R plots.
  -h, --help                  Show help.

Major outputs:
  guide_metadata.tsv
  transcript_map.tsv
  offtarget_alignments.tsv.gz
  guide_offtarget_summary.tsv
  gene_offtarget_summary.tsv
  high_risk_guides.tsv
  offtarget_gene_pairs.tsv
  off_target_qc_summary.txt
  plots/*.png if R/tidyverse are available

Interpretation notes:
  - "Off-target" means an alignment to a transcript/gene other than the intended gene/transcript.
  - Intended is called by matching transcript IDs from tx_id_set, gene_id_mouse, or gene_symbol_mouse.
  - This is transcriptome-level sequence matching. It does not account for RNA expression level,
    RNA structure, accessibility, Cas13 ortholog-specific mismatch tolerance, RNA editing, or indels.
  - Bowtie1 -v mode is ungapped; this script detects substitution mismatches, not bulges/indels.
USAGE
}

LIBRARY=""
TRANSCRIPTOME_FASTA=""
GTF=""
OUTDIR="cas13_transcriptome_offtarget_qc"
THREADS=8
MISMATCHES=3
QUERY_MODE="auto"
BOTH_STRANDS=false
KEEP_SAM=false
SKIP_INDEX=false
RUN_PLOTS=true

while [[ $# -gt 0 ]]; do
  case "$1" in
    --library) LIBRARY="$2"; shift 2 ;;
    --transcriptome-fasta) TRANSCRIPTOME_FASTA="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --mismatches) MISMATCHES="$2"; shift 2 ;;
    --query-mode) QUERY_MODE="$2"; shift 2 ;;
    --both-strands) BOTH_STRANDS=true; shift 1 ;;
    --keep-sam) KEEP_SAM=true; shift 1 ;;
    --skip-index) SKIP_INDEX=true; shift 1 ;;
    --no-plots) RUN_PLOTS=false; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$LIBRARY" || -z "$TRANSCRIPTOME_FASTA" ]]; then
  echo "ERROR: --library and --transcriptome-fasta are required." >&2
  usage
  exit 1
fi

if [[ ! -f "$LIBRARY" ]]; then
  echo "ERROR: library file not found: $LIBRARY" >&2
  exit 1
fi

if [[ ! -f "$TRANSCRIPTOME_FASTA" ]]; then
  echo "ERROR: transcriptome FASTA not found: $TRANSCRIPTOME_FASTA" >&2
  exit 1
fi

if [[ -n "$GTF" && ! -f "$GTF" ]]; then
  echo "ERROR: GTF file not found: $GTF" >&2
  exit 1
fi

if [[ "$MISMATCHES" -lt 0 || "$MISMATCHES" -gt 3 ]]; then
  echo "ERROR: --mismatches must be 0, 1, 2, or 3 for Bowtie1 -v mode." >&2
  exit 1
fi

case "$QUERY_MODE" in
  auto|target_seq|guide_seq|guide_rc) ;;
  *) echo "ERROR: --query-mode must be auto, target_seq, guide_seq, or guide_rc." >&2; exit 1 ;;
esac

command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found." >&2; exit 1; }
command -v bowtie >/dev/null 2>&1 || { echo "ERROR: bowtie not found. Load a Bowtie1 module, e.g. module load bowtie/1.3.1" >&2; exit 1; }
command -v bowtie-build >/dev/null 2>&1 || { echo "ERROR: bowtie-build not found. Load a Bowtie1 module." >&2; exit 1; }

mkdir -p "$OUTDIR" "$OUTDIR/reference" "$OUTDIR/bowtie_index" "$OUTDIR/logs" "$OUTDIR/plots"

GUIDE_FASTA="$OUTDIR/guide_queries.fa"
GUIDE_META="$OUTDIR/guide_metadata.tsv"
TRANSCRIPT_MAP="$OUTDIR/transcript_map.tsv"
REF_FASTA="$OUTDIR/reference/transcriptome.fa"
INDEX_PREFIX="$OUTDIR/bowtie_index/transcriptome"
SAM="$OUTDIR/bowtie_alignments.sam"
ALIGN_TSV_GZ="$OUTDIR/offtarget_alignments.tsv.gz"
GUIDE_SUMMARY="$OUTDIR/guide_offtarget_summary.tsv"
GENE_SUMMARY="$OUTDIR/gene_offtarget_summary.tsv"
HIGH_RISK="$OUTDIR/high_risk_guides.tsv"
GENE_PAIRS="$OUTDIR/offtarget_gene_pairs.tsv"
SUMMARY_TXT="$OUTDIR/off_target_qc_summary.txt"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log "Preparing guide query FASTA from library: $LIBRARY"
python3 - "$LIBRARY" "$GUIDE_FASTA" "$GUIDE_META" "$QUERY_MODE" <<'PY'
import csv
import re
import sys
from collections import defaultdict

library_csv, guide_fasta, guide_meta, query_mode = sys.argv[1:5]

CANDIDATES = {
    "guide_id": ["guide_id", "sgRNA_id", "grna_id", "id"],
    "gene_id": ["gene_id_mouse", "gene_id", "target_gene_id"],
    "gene_symbol": ["gene_symbol_mouse", "gene_symbol", "target_gene", "gene"],
    "guide_seq": ["guide_seq", "spacer", "spacer_seq", "guide_sequence", "sgRNA_sequence"],
    "target_seq": ["target_seq", "target_sequence"],
    "tx_id_set": ["tx_id_set", "transcript_id_set", "transcripts", "tx_ids"],
}

def choose(header, key, required=True):
    for c in CANDIDATES[key]:
        if c in header:
            return c
    if required:
        raise SystemExit(f"Missing required column for {key}; tried {CANDIDATES[key]}")
    return None

def clean_seq(s):
    if s is None:
        return ""
    s = str(s).upper().replace("U", "T")
    s = re.sub(r"[^ACGTN]", "", s)
    return s

def rc(s):
    tab = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(tab)[::-1].upper()

def clean_id(s):
    s = "" if s is None else str(s).strip()
    return re.sub(r"\.\d+$", "", s)

def split_tx_ids(s):
    if s is None:
        return []
    parts = re.split(r"[|;,\s]+", str(s).strip())
    out = []
    for p in parts:
        p = p.strip()
        if not p or p.upper() in {"NA", "NAN", "NULL"}:
            continue
        out.append(p)
    return out

def safe_qname(s):
    s = str(s).strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.:+\-=|]", "_", s)
    return s

with open(library_csv, newline="") as f:
    reader = csv.DictReader(f)
    header = reader.fieldnames or []
    col_guide_id = choose(header, "guide_id")
    col_gene_id = choose(header, "gene_id", required=False)
    col_gene_symbol = choose(header, "gene_symbol")
    col_guide_seq = choose(header, "guide_seq")
    col_target_seq = choose(header, "target_seq", required=False)
    col_tx_id_set = choose(header, "tx_id_set", required=False)

    rows = list(reader)

seen_qnames = defaultdict(int)
n_written = 0
n_skipped = 0

with open(guide_fasta, "w") as fa, open(guide_meta, "w", newline="") as meta_f:
    fields = [
        "guide_id", "original_guide_id", "intended_gene_id", "intended_gene_id_clean",
        "intended_gene_symbol", "intended_gene_symbol_upper", "guide_seq", "target_seq",
        "query_seq", "query_seq_source", "query_length", "tx_id_set", "tx_ids_clean"
    ]
    writer = csv.DictWriter(meta_f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()

    for row in rows:
        original_id = row.get(col_guide_id, "")
        qname_base = safe_qname(original_id)
        seen_qnames[qname_base] += 1
        qname = qname_base if seen_qnames[qname_base] == 1 else f"{qname_base}__dup{seen_qnames[qname_base]}"

        guide_seq = clean_seq(row.get(col_guide_seq, ""))
        target_seq = clean_seq(row.get(col_target_seq, "")) if col_target_seq else ""

        if query_mode == "auto":
            if target_seq:
                query_seq = target_seq
                source = "target_seq"
            else:
                query_seq = rc(guide_seq)
                source = "guide_seq_reverse_complement"
        elif query_mode == "target_seq":
            query_seq = target_seq
            source = "target_seq"
        elif query_mode == "guide_seq":
            query_seq = guide_seq
            source = "guide_seq_as_given"
        elif query_mode == "guide_rc":
            query_seq = rc(guide_seq)
            source = "guide_seq_reverse_complement"
        else:
            raise SystemExit(f"Unknown query mode: {query_mode}")

        if not query_seq or re.search(r"[^ACGT]", query_seq):
            n_skipped += 1
            continue

        intended_gene_id = row.get(col_gene_id, "") if col_gene_id else ""
        intended_gene_symbol = row.get(col_gene_symbol, "")
        tx_id_set = row.get(col_tx_id_set, "") if col_tx_id_set else ""
        tx_ids_clean = sorted(set(clean_id(x) for x in split_tx_ids(tx_id_set)))

        fa.write(f">{qname}\n{query_seq}\n")
        writer.writerow({
            "guide_id": qname,
            "original_guide_id": original_id,
            "intended_gene_id": intended_gene_id,
            "intended_gene_id_clean": clean_id(intended_gene_id),
            "intended_gene_symbol": intended_gene_symbol,
            "intended_gene_symbol_upper": str(intended_gene_symbol).upper(),
            "guide_seq": guide_seq,
            "target_seq": target_seq,
            "query_seq": query_seq,
            "query_seq_source": source,
            "query_length": len(query_seq),
            "tx_id_set": tx_id_set,
            "tx_ids_clean": "|".join(tx_ids_clean),
        })
        n_written += 1

sys.stderr.write(f"Wrote {n_written} guide query sequences to {guide_fasta}\n")
if n_skipped:
    sys.stderr.write(f"Skipped {n_skipped} rows with missing/ambiguous query sequences.\n")
PY

log "Preparing transcriptome FASTA copy/symlink"
if [[ "$TRANSCRIPTOME_FASTA" == *.gz ]]; then
  gzip -dc "$TRANSCRIPTOME_FASTA" > "$REF_FASTA"
else
  # Use symlink to avoid copying a large FASTA when possible.
  ln -sf "$(readlink -f "$TRANSCRIPTOME_FASTA")" "$REF_FASTA"
fi

log "Parsing transcriptome FASTA/GTF transcript map"
python3 - "$REF_FASTA" "${GTF:-}" "$TRANSCRIPT_MAP" <<'PY'
import csv
import gzip
import re
import sys

ref_fasta, gtf_path, out_path = sys.argv[1:4]

def open_text(path):
    if not path:
        return None
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def clean_id(s):
    s = "" if s is None else str(s).strip()
    return re.sub(r"\.\d+$", "", s)

def parse_attrs(attr):
    d = {}
    # GTF style: key "value";
    for key, val in re.findall(r'([A-Za-z0-9_]+)\s+"([^"]*)"', attr):
        d[key] = val
    # GFF-ish: key=value;key=value
    if not d:
        for item in attr.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                d[k.strip()] = v.strip()
    return d

def parse_fasta_header(header):
    # header does not include leading >
    first_token = header.split()[0]
    pipe = first_token.split("|")
    transcript_id = pipe[0]
    gene_id = ""
    gene_symbol = ""
    gene_biotype = ""
    transcript_biotype = ""

    # GENCODE-like pipe headers: transcript|gene|...|transcript_name|gene_name|...
    if len(pipe) >= 2 and re.match(r"ENS[A-Z]*G", pipe[1]):
      gene_id = pipe[1]
    if len(pipe) >= 7:
      gene_symbol = pipe[6]
    elif len(pipe) >= 6:
      gene_symbol = pipe[5]

    # Ensembl-like key:value fields after first token.
    for key, val in re.findall(r"([A-Za-z_]+):([^\s]+)", header):
        kl = key.lower()
        if kl == "gene":
            gene_id = val
        elif kl in {"gene_symbol", "gene_name"}:
            gene_symbol = val
        elif kl in {"gene_biotype", "gene_type"}:
            gene_biotype = val
        elif kl in {"transcript_biotype", "transcript_type"}:
            transcript_biotype = val

    return {
        "ref_name": first_token,
        "transcript_id": transcript_id,
        "transcript_id_clean": clean_id(transcript_id),
        "fasta_gene_id": gene_id,
        "fasta_gene_id_clean": clean_id(gene_id),
        "fasta_gene_symbol": gene_symbol,
        "fasta_gene_symbol_upper": gene_symbol.upper(),
        "fasta_gene_biotype": gene_biotype,
        "fasta_transcript_biotype": transcript_biotype,
        "full_fasta_header": header,
    }

gtf_by_tx = {}
if gtf_path:
    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_attrs(parts[8])
            tx = attrs.get("transcript_id") or attrs.get("transcriptId") or attrs.get("ID", "")
            if not tx:
                continue
            tx_clean = clean_id(tx)
            if tx_clean in gtf_by_tx:
                continue
            gene_id = attrs.get("gene_id") or attrs.get("geneId") or attrs.get("gene", "")
            gene_symbol = attrs.get("gene_name") or attrs.get("gene_symbol") or attrs.get("gene") or ""
            gene_biotype = attrs.get("gene_type") or attrs.get("gene_biotype") or attrs.get("gene_bio_type") or ""
            transcript_biotype = attrs.get("transcript_type") or attrs.get("transcript_biotype") or ""
            gtf_by_tx[tx_clean] = {
                "gtf_transcript_id": tx,
                "gtf_gene_id": gene_id,
                "gtf_gene_id_clean": clean_id(gene_id),
                "gtf_gene_symbol": gene_symbol,
                "gtf_gene_symbol_upper": gene_symbol.upper(),
                "gtf_gene_biotype": gene_biotype,
                "gtf_transcript_biotype": transcript_biotype,
            }

records = []
with open(ref_fasta, "r") as fh:
    for line in fh:
        if line.startswith(">"):
            rec = parse_fasta_header(line[1:].rstrip("\n"))
            gtf = gtf_by_tx.get(rec["transcript_id_clean"], {})
            gene_id = gtf.get("gtf_gene_id") or rec.get("fasta_gene_id") or ""
            gene_symbol = gtf.get("gtf_gene_symbol") or rec.get("fasta_gene_symbol") or ""
            gene_biotype = gtf.get("gtf_gene_biotype") or rec.get("fasta_gene_biotype") or ""
            transcript_biotype = gtf.get("gtf_transcript_biotype") or rec.get("fasta_transcript_biotype") or ""
            rec.update(gtf)
            rec.update({
                "gene_id": gene_id,
                "gene_id_clean": clean_id(gene_id),
                "gene_symbol": gene_symbol,
                "gene_symbol_upper": gene_symbol.upper(),
                "gene_biotype": gene_biotype,
                "transcript_biotype": transcript_biotype,
            })
            records.append(rec)

fields = [
    "ref_name", "transcript_id", "transcript_id_clean", "gene_id", "gene_id_clean",
    "gene_symbol", "gene_symbol_upper", "gene_biotype", "transcript_biotype",
    "fasta_gene_id", "fasta_gene_symbol", "gtf_gene_id", "gtf_gene_symbol", "full_fasta_header"
]
with open(out_path, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    for rec in records:
        writer.writerow(rec)

sys.stderr.write(f"Wrote transcript map for {len(records)} FASTA records to {out_path}\n")
if gtf_path:
    sys.stderr.write(f"Parsed GTF transcript annotations for {len(gtf_by_tx)} transcript IDs\n")
PY

if [[ "$SKIP_INDEX" == false || ! -f "${INDEX_PREFIX}.1.ebwt" ]]; then
  log "Building Bowtie1 transcriptome index"
  bowtie-build -f "$REF_FASTA" "$INDEX_PREFIX" > "$OUTDIR/logs/bowtie_build.log" 2>&1
else
  log "Reusing existing Bowtie index: $INDEX_PREFIX"
fi

BOWTIE_STRAND_ARGS=(--norc)
if [[ "$BOTH_STRANDS" == true ]]; then
  BOWTIE_STRAND_ARGS=()
fi

log "Aligning guide queries to transcriptome with Bowtie1 -v $MISMATCHES -a"
log "Strand mode: $([[ "$BOTH_STRANDS" == true ]] && echo 'both orientations' || echo 'forward-only query alignment via --norc')"
bowtie \
  -f \
  -v "$MISMATCHES" \
  -a \
  --best \
  -S \
  -p "$THREADS" \
  "${BOWTIE_STRAND_ARGS[@]}" \
  "$INDEX_PREFIX" \
  "$GUIDE_FASTA" \
  "$SAM" \
  2> "$OUTDIR/logs/bowtie_align.log"

log "Parsing SAM and summarizing off-targets"
python3 - "$SAM" "$GUIDE_META" "$TRANSCRIPT_MAP" "$ALIGN_TSV_GZ" "$GUIDE_SUMMARY" "$GENE_SUMMARY" "$HIGH_RISK" "$GENE_PAIRS" "$SUMMARY_TXT" "$MISMATCHES" <<'PY'
import csv
import gzip
import math
import re
import sys
from collections import defaultdict, Counter

sam_path, guide_meta_path, transcript_map_path, align_out_gz, guide_summary_path, gene_summary_path, high_risk_path, gene_pairs_path, summary_txt_path, max_mismatches_s = sys.argv[1:11]
MAX_MISMATCHES = int(max_mismatches_s)
THRESHOLDS = list(range(MAX_MISMATCHES + 1))

def clean_id(s):
    s = "" if s is None else str(s).strip()
    return re.sub(r"\.\d+$", "", s)

def split_bar(s):
    if not s:
        return set()
    return set(x for x in str(s).split("|") if x)

def parse_nm(tags):
    nm = None
    xm = None
    for t in tags:
        if t.startswith("NM:i:"):
            try: nm = int(t.split(":")[-1])
            except ValueError: pass
        elif t.startswith("XM:i:"):
            try: xm = int(t.split(":")[-1])
            except ValueError: pass
    if nm is not None:
        return nm
    if xm is not None:
        return xm
    return None

def median(vals):
    vals = sorted(vals)
    n = len(vals)
    if n == 0:
        return "NA"
    mid = n // 2
    if n % 2:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2

def fmt_set(vals, max_items=40):
    vals = sorted(x for x in vals if x)
    if len(vals) > max_items:
        return "|".join(vals[:max_items]) + f"|...(+{len(vals)-max_items})"
    return "|".join(vals)

# Read guide metadata.
guides = {}
with open(guide_meta_path, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        row["tx_ids_clean_set"] = split_bar(row.get("tx_ids_clean", ""))
        guides[row["guide_id"]] = row

# Read transcript map.
transcripts = {}
with open(transcript_map_path, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        transcripts[row["ref_name"]] = row

# Per-guide accumulator.
def init_guide_acc():
    return {
        "n_alignments_total": 0,
        "n_intended_alignments_total": 0,
        "best_intended_mismatches": None,
        "best_any_mismatches": None,
        "hits_by_exact_mismatch": Counter(),
        "intended_hits_by_exact_mismatch": Counter(),
        "offtarget_hits_by_exact_mismatch": Counter(),
        "total_hits_le": {t: 0 for t in THRESHOLDS},
        "intended_hits_le": {t: 0 for t in THRESHOLDS},
        "offtarget_hits_le": {t: 0 for t in THRESHOLDS},
        "offtarget_tx_le": {t: set() for t in THRESHOLDS},
        "offtarget_gene_le": {t: set() for t in THRESHOLDS},
        "offtarget_gene_symbol_le": {t: set() for t in THRESHOLDS},
        "intended_tx_le": {t: set() for t in THRESHOLDS},
    }

guide_acc = {gid: init_guide_acc() for gid in guides}

# Guide/off-target-gene accumulator.
gene_pair_acc = {}

def pair_key(g, tx):
    intended_gene_id = g.get("intended_gene_id_clean", "")
    intended_symbol = g.get("intended_gene_symbol", "")
    target_gene_id = tx.get("gene_id_clean", "") or tx.get("gene_id", "") or "unknown_gene_id"
    target_symbol = tx.get("gene_symbol", "") or "unknown_symbol"
    return (g["guide_id"], intended_gene_id, intended_symbol, target_gene_id, target_symbol)

align_fields = [
    "guide_id", "original_guide_id", "intended_gene_id", "intended_gene_symbol",
    "query_seq", "query_seq_source", "query_length",
    "transcript_ref_name", "transcript_id", "transcript_id_clean",
    "target_gene_id", "target_gene_symbol", "target_gene_biotype", "target_transcript_biotype",
    "alignment_pos_1based", "alignment_strand", "cigar", "mismatches",
    "is_intended_transcript", "is_intended_gene", "is_off_target"
]

n_sam_records = 0
n_written_alignments = 0
n_unknown_guides = 0
n_unknown_transcripts = 0

with gzip.open(align_out_gz, "wt", newline="") as out_f:
    writer = csv.DictWriter(out_f, fieldnames=align_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()

    with open(sam_path, "r") as sam:
        for line in sam:
            if not line or line.startswith("@"):
                continue
            n_sam_records += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            qname, flag_s, rname, pos_s, mapq, cigar = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
            if rname == "*":
                continue
            g = guides.get(qname)
            if g is None:
                n_unknown_guides += 1
                continue
            tx = transcripts.get(rname)
            if tx is None:
                n_unknown_transcripts += 1
                tx = {
                    "ref_name": rname,
                    "transcript_id": rname.split("|")[0],
                    "transcript_id_clean": clean_id(rname.split("|")[0]),
                    "gene_id": "",
                    "gene_id_clean": "",
                    "gene_symbol": "",
                    "gene_symbol_upper": "",
                    "gene_biotype": "",
                    "transcript_biotype": "",
                }

            try:
                flag = int(flag_s)
            except ValueError:
                flag = 0
            strand = "-" if (flag & 16) else "+"
            mismatches = parse_nm(parts[11:])
            if mismatches is None:
                # Bowtie1 SAM should report XM/NM tags; if absent, mark as NA and skip threshold summaries.
                continue

            guide_tx_ids = g.get("tx_ids_clean_set", set())
            tx_clean = tx.get("transcript_id_clean", "")
            is_intended_tx = bool(tx_clean and tx_clean in guide_tx_ids)

            intended_gene_id = g.get("intended_gene_id_clean", "")
            intended_symbol_upper = g.get("intended_gene_symbol_upper", "")
            target_gene_id = tx.get("gene_id_clean", "")
            target_symbol_upper = tx.get("gene_symbol_upper", "")
            is_intended_gene = bool(
                is_intended_tx or
                (intended_gene_id and target_gene_id and intended_gene_id == target_gene_id) or
                (intended_symbol_upper and target_symbol_upper and intended_symbol_upper == target_symbol_upper)
            )
            is_off_target = not is_intended_gene

            acc = guide_acc[qname]
            acc["n_alignments_total"] += 1
            acc["hits_by_exact_mismatch"][mismatches] += 1
            if acc["best_any_mismatches"] is None or mismatches < acc["best_any_mismatches"]:
                acc["best_any_mismatches"] = mismatches

            if is_intended_gene:
                acc["n_intended_alignments_total"] += 1
                acc["intended_hits_by_exact_mismatch"][mismatches] += 1
                if acc["best_intended_mismatches"] is None or mismatches < acc["best_intended_mismatches"]:
                    acc["best_intended_mismatches"] = mismatches
            else:
                acc["offtarget_hits_by_exact_mismatch"][mismatches] += 1

            for t in THRESHOLDS:
                if mismatches <= t:
                    acc["total_hits_le"][t] += 1
                    if is_intended_gene:
                        acc["intended_hits_le"][t] += 1
                        acc["intended_tx_le"][t].add(tx_clean or tx.get("transcript_id", ""))
                    else:
                        acc["offtarget_hits_le"][t] += 1
                        acc["offtarget_tx_le"][t].add(tx_clean or tx.get("transcript_id", ""))
                        off_gene = target_gene_id or tx.get("gene_id", "") or tx_clean or rname
                        acc["offtarget_gene_le"][t].add(off_gene)
                        if tx.get("gene_symbol", ""):
                            acc["offtarget_gene_symbol_le"][t].add(tx.get("gene_symbol", ""))

            if is_off_target:
                pk = pair_key(g, tx)
                if pk not in gene_pair_acc:
                    gene_pair_acc[pk] = {
                        "n_alignment_records": 0,
                        "min_mismatches": mismatches,
                        "transcripts": set(),
                    }
                gene_pair_acc[pk]["n_alignment_records"] += 1
                gene_pair_acc[pk]["min_mismatches"] = min(gene_pair_acc[pk]["min_mismatches"], mismatches)
                gene_pair_acc[pk]["transcripts"].add(tx_clean or tx.get("transcript_id", ""))

            writer.writerow({
                "guide_id": qname,
                "original_guide_id": g.get("original_guide_id", qname),
                "intended_gene_id": g.get("intended_gene_id", ""),
                "intended_gene_symbol": g.get("intended_gene_symbol", ""),
                "query_seq": g.get("query_seq", ""),
                "query_seq_source": g.get("query_seq_source", ""),
                "query_length": g.get("query_length", ""),
                "transcript_ref_name": rname,
                "transcript_id": tx.get("transcript_id", ""),
                "transcript_id_clean": tx_clean,
                "target_gene_id": tx.get("gene_id", ""),
                "target_gene_symbol": tx.get("gene_symbol", ""),
                "target_gene_biotype": tx.get("gene_biotype", ""),
                "target_transcript_biotype": tx.get("transcript_biotype", ""),
                "alignment_pos_1based": pos_s,
                "alignment_strand": strand,
                "cigar": cigar,
                "mismatches": mismatches,
                "is_intended_transcript": str(is_intended_tx),
                "is_intended_gene": str(is_intended_gene),
                "is_off_target": str(is_off_target),
            })
            n_written_alignments += 1

# Write guide summary.
guide_summary_fields = [
    "guide_id", "original_guide_id", "intended_gene_id", "intended_gene_symbol",
    "query_seq", "query_seq_source", "query_length", "tx_id_set",
    "n_alignments_total", "n_intended_alignments_total", "best_any_mismatches",
    "best_intended_mismatches", "no_intended_hit_le_max_mismatches"
]
for t in THRESHOLDS:
    guide_summary_fields += [
        f"total_hits_le{t}mm", f"intended_hits_le{t}mm", f"offtarget_hits_le{t}mm",
        f"offtarget_transcripts_le{t}mm", f"offtarget_genes_le{t}mm",
    ]
for t in THRESHOLDS:
    guide_summary_fields += [f"all_hits_exact_{t}mm", f"intended_hits_exact_{t}mm", f"offtarget_hits_exact_{t}mm"]
guide_summary_fields += [
    "offtarget_gene_symbols_le1mm", "offtarget_gene_symbols_le2mm",
    "flag_exact_offtarget", "flag_offtarget_le1mm", "flag_offtarget_le2mm", "flag_no_intended_hit"
]

summary_rows = []
with open(guide_summary_path, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=guide_summary_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for gid, g in guides.items():
        acc = guide_acc[gid]
        no_intended = acc["n_intended_alignments_total"] == 0
        row = {
            "guide_id": gid,
            "original_guide_id": g.get("original_guide_id", gid),
            "intended_gene_id": g.get("intended_gene_id", ""),
            "intended_gene_symbol": g.get("intended_gene_symbol", ""),
            "query_seq": g.get("query_seq", ""),
            "query_seq_source": g.get("query_seq_source", ""),
            "query_length": g.get("query_length", ""),
            "tx_id_set": g.get("tx_id_set", ""),
            "n_alignments_total": acc["n_alignments_total"],
            "n_intended_alignments_total": acc["n_intended_alignments_total"],
            "best_any_mismatches": "NA" if acc["best_any_mismatches"] is None else acc["best_any_mismatches"],
            "best_intended_mismatches": "NA" if acc["best_intended_mismatches"] is None else acc["best_intended_mismatches"],
            "no_intended_hit_le_max_mismatches": str(no_intended),
        }
        for t in THRESHOLDS:
            row[f"total_hits_le{t}mm"] = acc["total_hits_le"][t]
            row[f"intended_hits_le{t}mm"] = acc["intended_hits_le"][t]
            row[f"offtarget_hits_le{t}mm"] = acc["offtarget_hits_le"][t]
            row[f"offtarget_transcripts_le{t}mm"] = len(acc["offtarget_tx_le"][t])
            row[f"offtarget_genes_le{t}mm"] = len(acc["offtarget_gene_le"][t])
        for t in THRESHOLDS:
            row[f"all_hits_exact_{t}mm"] = acc["hits_by_exact_mismatch"][t]
            row[f"intended_hits_exact_{t}mm"] = acc["intended_hits_by_exact_mismatch"][t]
            row[f"offtarget_hits_exact_{t}mm"] = acc["offtarget_hits_by_exact_mismatch"][t]
        row["offtarget_gene_symbols_le1mm"] = fmt_set(acc["offtarget_gene_symbol_le"].get(1, set())) if 1 in acc["offtarget_gene_symbol_le"] else ""
        row["offtarget_gene_symbols_le2mm"] = fmt_set(acc["offtarget_gene_symbol_le"].get(2, set())) if 2 in acc["offtarget_gene_symbol_le"] else ""
        row["flag_exact_offtarget"] = str(acc["offtarget_hits_le"].get(0, 0) > 0)
        row["flag_offtarget_le1mm"] = str(acc["offtarget_hits_le"].get(1, 0) > 0)
        row["flag_offtarget_le2mm"] = str(acc["offtarget_hits_le"].get(2, 0) > 0)
        row["flag_no_intended_hit"] = str(no_intended)
        writer.writerow(row)
        summary_rows.append(row)

# High-risk guide table.
with open(high_risk_path, "w", newline="") as out:
    fields = guide_summary_fields
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    def risk_key(r):
        return (
            r["flag_exact_offtarget"] == "True",
            r["flag_offtarget_le1mm"] == "True",
            r["flag_offtarget_le2mm"] == "True",
            r["flag_no_intended_hit"] == "True",
            int(r.get("offtarget_genes_le2mm", 0)),
            int(r.get("offtarget_hits_le2mm", 0)),
        )
    for row in sorted(
        [r for r in summary_rows if r["flag_exact_offtarget"] == "True" or r["flag_offtarget_le1mm"] == "True" or r["flag_offtarget_le2mm"] == "True" or r["flag_no_intended_hit"] == "True"],
        key=risk_key,
        reverse=True,
    ):
        writer.writerow(row)

# Gene-pair table.
pair_fields = [
    "guide_id", "intended_gene_id", "intended_gene_symbol", "offtarget_gene_id", "offtarget_gene_symbol",
    "min_mismatches", "n_alignment_records", "n_offtarget_transcripts", "offtarget_transcripts"
]
with open(gene_pairs_path, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=pair_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for (guide_id, intended_gene_id, intended_symbol, target_gene_id, target_symbol), acc in sorted(gene_pair_acc.items(), key=lambda kv: (kv[1]["min_mismatches"], -kv[1]["n_alignment_records"], kv[0])):
        writer.writerow({
            "guide_id": guide_id,
            "intended_gene_id": intended_gene_id,
            "intended_gene_symbol": intended_symbol,
            "offtarget_gene_id": target_gene_id,
            "offtarget_gene_symbol": target_symbol,
            "min_mismatches": acc["min_mismatches"],
            "n_alignment_records": acc["n_alignment_records"],
            "n_offtarget_transcripts": len(acc["transcripts"]),
            "offtarget_transcripts": fmt_set(acc["transcripts"], max_items=80),
        })

# Gene-level summary from guide summaries.
gene_groups = defaultdict(list)
for row in summary_rows:
    key = (row.get("intended_gene_id", ""), row.get("intended_gene_symbol", ""))
    gene_groups[key].append(row)

gene_fields = [
    "intended_gene_id", "intended_gene_symbol", "n_guides",
    "n_guides_no_intended_hit", "n_guides_with_exact_offtarget",
    "n_guides_with_offtarget_le1mm", "n_guides_with_offtarget_le2mm", f"n_guides_with_offtarget_le{MAX_MISMATCHES}mm",
    "median_offtarget_genes_le2mm", "max_offtarget_genes_le2mm",
    "median_offtarget_hits_le2mm", "max_offtarget_hits_le2mm",
    "guides_with_offtarget_le1mm", "guides_no_intended_hit"
]
with open(gene_summary_path, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=gene_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for (gid, sym), rows in sorted(gene_groups.items(), key=lambda kv: kv[0][1]):
        off_genes_le2 = [int(r.get("offtarget_genes_le2mm", 0)) for r in rows]
        off_hits_le2 = [int(r.get("offtarget_hits_le2mm", 0)) for r in rows]
        guides_le1 = [r["guide_id"] for r in rows if r.get("flag_offtarget_le1mm") == "True"]
        guides_no_intended = [r["guide_id"] for r in rows if r.get("flag_no_intended_hit") == "True"]
        writer.writerow({
            "intended_gene_id": gid,
            "intended_gene_symbol": sym,
            "n_guides": len(rows),
            "n_guides_no_intended_hit": sum(r.get("flag_no_intended_hit") == "True" for r in rows),
            "n_guides_with_exact_offtarget": sum(r.get("flag_exact_offtarget") == "True" for r in rows),
            "n_guides_with_offtarget_le1mm": sum(r.get("flag_offtarget_le1mm") == "True" for r in rows),
            "n_guides_with_offtarget_le2mm": sum(r.get("flag_offtarget_le2mm") == "True" for r in rows),
            f"n_guides_with_offtarget_le{MAX_MISMATCHES}mm": sum(int(r.get(f"offtarget_hits_le{MAX_MISMATCHES}mm", 0)) > 0 for r in rows),
            "median_offtarget_genes_le2mm": median(off_genes_le2),
            "max_offtarget_genes_le2mm": max(off_genes_le2) if off_genes_le2 else 0,
            "median_offtarget_hits_le2mm": median(off_hits_le2),
            "max_offtarget_hits_le2mm": max(off_hits_le2) if off_hits_le2 else 0,
            "guides_with_offtarget_le1mm": "|".join(guides_le1),
            "guides_no_intended_hit": "|".join(guides_no_intended),
        })

# Text summary.
n_guides = len(summary_rows)
summary = []
summary.append("Cas13 transcriptome off-target QC summary")
summary.append(f"Guides queried: {n_guides}")
summary.append(f"SAM alignment records parsed: {n_sam_records}")
summary.append(f"Alignment rows written: {n_written_alignments}")
summary.append(f"Unknown guide IDs in SAM: {n_unknown_guides}")
summary.append(f"Unknown transcript IDs in transcript map: {n_unknown_transcripts}")
summary.append("")
summary.append("Guide-level headline counts:")
summary.append(f"Guides with no intended transcript/gene hit <= {MAX_MISMATCHES} mismatches: {sum(r['flag_no_intended_hit'] == 'True' for r in summary_rows)}")
for t in THRESHOLDS:
    n_guides_off = sum(int(r.get(f"offtarget_hits_le{t}mm", 0)) > 0 for r in summary_rows)
    n_guides_off_gene = sum(int(r.get(f"offtarget_genes_le{t}mm", 0)) > 0 for r in summary_rows)
    summary.append(f"Guides with any off-target alignment <= {t} mismatches: {n_guides_off} ({n_guides_off / n_guides * 100:.2f}%)")
    summary.append(f"Guides with any off-target gene <= {t} mismatches: {n_guides_off_gene} ({n_guides_off_gene / n_guides * 100:.2f}%)")
summary.append("")
summary.append("Output files:")
summary.append(f"- {align_out_gz}")
summary.append(f"- {guide_summary_path}")
summary.append(f"- {gene_summary_path}")
summary.append(f"- {high_risk_path}")
summary.append(f"- {gene_pairs_path}")

with open(summary_txt_path, "w") as out:
    out.write("\n".join(summary) + "\n")

sys.stderr.write("\n".join(summary) + "\n")
PY

if [[ "$KEEP_SAM" == true ]]; then
  log "Keeping SAM file: $SAM"
else
  rm -f "$SAM"
  log "Removed temporary SAM file. Use --keep-sam if you want to retain it."
fi

if [[ "$RUN_PLOTS" == true ]]; then
  if command -v Rscript >/dev/null 2>&1; then
    log "Trying to make optional R plots"
    cat > "$OUTDIR/make_offtarget_plots.R" <<'RSCRIPT'
suppressPackageStartupMessages({
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    stop("tidyverse is not installed; skipping plots")
  }
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[[1]]
plot_dir <- file.path(outdir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

guide_summary <- readr::read_tsv(file.path(outdir, "guide_offtarget_summary.tsv"), show_col_types = FALSE)
gene_summary <- readr::read_tsv(file.path(outdir, "gene_offtarget_summary.tsv"), show_col_types = FALSE)

comma_number <- scales::label_comma(accuracy = 1)
qc_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "grey35"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      plot.margin = margin(10, 32, 10, 10)
    )
}

save_plot <- function(p, name, width = 9, height = 6) {
  ggsave(file.path(plot_dir, name), p, width = width, height = height, dpi = 300)
}

truthy <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")

threshold_cols <- names(guide_summary)[stringr::str_detect(names(guide_summary), "^offtarget_genes_le[0-9]+mm$")]
threshold_counts <- purrr::map_dfr(threshold_cols, function(col) {
  thr <- stringr::str_match(col, "le([0-9]+)mm")[,2]
  tibble(
    threshold = paste0("<=", thr, " mismatches"),
    threshold_num = as.integer(thr),
    n_guides = sum(guide_summary[[col]] > 0, na.rm = TRUE),
    pct_guides = mean(guide_summary[[col]] > 0, na.rm = TRUE) * 100
  )
}) %>% arrange(threshold_num) %>% mutate(threshold = factor(threshold, levels = threshold))

p1 <- ggplot(threshold_counts, aes(x = threshold, y = n_guides)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(comma_number(n_guides), "\n", sprintf("%.1f%%", pct_guides))), vjust = -0.25, size = 3.3) +
  scale_y_continuous(labels = comma_number, expand = expansion(mult = c(0, 0.22))) +
  labs(
    title = "Guides with transcriptome off-target genes",
    subtitle = "A guide is counted if it aligns to at least one non-intended gene at or below the mismatch threshold.",
    x = NULL,
    y = "Number of guides"
  ) +
  qc_theme()
save_plot(p1, "guides_with_offtarget_genes_by_mismatch.png", width = 9.5, height = 6)

if ("offtarget_genes_le2mm" %in% names(guide_summary)) {
  dist2 <- guide_summary %>%
    mutate(bucket = case_when(
      offtarget_genes_le2mm >= 10 ~ "10+",
      TRUE ~ as.character(offtarget_genes_le2mm)
    )) %>%
    count(bucket, name = "n") %>%
    mutate(bucket = factor(bucket, levels = c(as.character(0:9), "10+")))
  p2 <- ggplot(dist2, aes(x = bucket, y = n)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = comma_number(n)), vjust = -0.25, size = 3.2) +
    scale_y_continuous(labels = comma_number, expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "Off-target gene burden per guide at <=2 mismatches",
      x = "Distinct off-target genes per guide",
      y = "Number of guides"
    ) +
    qc_theme()
  save_plot(p2, "offtarget_gene_burden_le2mm.png", width = 9.5, height = 6)
}

status_counts <- guide_summary %>%
  mutate(status = if_else(truthy(flag_no_intended_hit), "No intended hit", "Intended hit found")) %>%
  count(status, name = "n") %>%
  mutate(status = factor(status, levels = c("Intended hit found", "No intended hit")))

p3 <- ggplot(status_counts, aes(x = status, y = n)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma_number(n)), vjust = -0.25, size = 3.4) +
  scale_y_continuous(labels = comma_number, expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Intended transcript/gene hit recovery",
    subtitle = "Intended hits are called using tx_id_set, gene_id_mouse, or gene_symbol_mouse.",
    x = NULL,
    y = "Number of guides"
  ) +
  qc_theme()
save_plot(p3, "intended_hit_recovery.png", width = 8, height = 6)

flag_counts <- tibble(
  flag = c("Exact off-target", "Off-target <=1 mm", "Off-target <=2 mm", "No intended hit"),
  n = c(
    sum(truthy(guide_summary$flag_exact_offtarget), na.rm = TRUE),
    sum(truthy(guide_summary$flag_offtarget_le1mm), na.rm = TRUE),
    sum(truthy(guide_summary$flag_offtarget_le2mm), na.rm = TRUE),
    sum(truthy(guide_summary$flag_no_intended_hit), na.rm = TRUE)
  )
) %>% mutate(flag = factor(flag, levels = rev(flag)))

p4 <- ggplot(flag_counts, aes(x = flag, y = n)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma_number(n)), hjust = -0.12, size = 3.4) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = comma_number, expand = expansion(mult = c(0, 0.22))) +
  labs(
    title = "High-level off-target QC flags",
    x = NULL,
    y = "Number of guides"
  ) +
  qc_theme()
save_plot(p4, "offtarget_qc_flag_counts.png", width = 9.5, height = 6)

pdf(file.path(outdir, "offtarget_qc_plots.pdf"), width = 9.5, height = 6)
print(p1)
if (exists("p2")) print(p2)
print(p3)
print(p4)
dev.off()
RSCRIPT

    if Rscript "$OUTDIR/make_offtarget_plots.R" "$OUTDIR" > "$OUTDIR/logs/offtarget_plots.log" 2>&1; then
      log "Optional plots written to $OUTDIR/plots and $OUTDIR/offtarget_qc_plots.pdf"
    else
      log "Optional R plotting failed; alignment summaries are still complete. See $OUTDIR/logs/offtarget_plots.log"
    fi
  else
    log "Rscript not found; skipping optional plots."
  fi
fi

log "Done. Main summary: $SUMMARY_TXT"
