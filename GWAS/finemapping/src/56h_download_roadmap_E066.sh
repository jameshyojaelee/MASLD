#!/bin/bash
# 56h_download_roadmap_E066.sh
# Download Roadmap Epigenomics E066 (Liver) chromHMM 15-state + H3K27ac + H3K4me3 peaks.
# NOTE: E066 (Liver) was NOT assayed for DNase by Roadmap (chromHMM is imputed/core-marks).
# We substitute H3K27ac broadPeak (active enhancers/promoters) and H3K4me3 narrowPeak (active TSS)
# as direct empirical regulatory annotations. Both are hg19 and need liftOver to hg38.
# chromHMM is already hg38-lifted by Roadmap.

set -eo pipefail

OUTDIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/roadmap_liver_E066
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Source: Roadmap Epigenomics Project (egg2.wustl.edu)
CHMM_URL="https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E066_15_coreMarks_hg38lift_segments.bed.gz"
H3K27AC_URL="https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E066-H3K27ac.broadPeak.gz"
H3K4ME3_URL="https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E066-H3K4me3.narrowPeak.gz"

CHMM_BED="$OUTDIR/E066_15_chromHMM_hg38.bed.gz"
H3K27AC_HG19="$OUTDIR/E066_H3K27ac_hg19.broadPeak.gz"
H3K4ME3_HG19="$OUTDIR/E066_H3K4me3_hg19.narrowPeak.gz"
H3K27AC_HG38="$OUTDIR/E066_H3K27ac_hg38.broadPeak.gz"
H3K4ME3_HG38="$OUTDIR/E066_H3K4me3_hg38.narrowPeak.gz"

download_or_skip() {
    local url=$1 out=$2 label=$3
    if [[ ! -s "$out" ]]; then
        echo "[$(date)] Downloading $label: $url"
        curl --retry 5 --retry-delay 30 --connect-timeout 60 -L -o "$out" "$url"
        # Validate by extracting first 3 lines into tmp file (avoids SIGPIPE issues with pipefail)
        local first_lines
        first_lines=$(zcat "$out" 2>/dev/null | { head -3; cat >/dev/null; } || true)
        if ! echo "$first_lines" | grep -q "^chr"; then
            echo "ERROR: $label download appears malformed (no 'chr' rows)" >&2
            echo "First lines were:" >&2
            echo "$first_lines" >&2
            rm -f "$out"
            exit 1
        fi
        echo "  $label validated; first row: $(echo "$first_lines" | head -1)"
    else
        echo "$label already downloaded: $out"
    fi
}

download_or_skip "$CHMM_URL"    "$CHMM_BED"     "chromHMM E066 (hg38)"
download_or_skip "$H3K27AC_URL" "$H3K27AC_HG19" "H3K27ac broadPeak (hg19)"
download_or_skip "$H3K4ME3_URL" "$H3K4ME3_HG19" "H3K4me3 narrowPeak (hg19)"

# Liftover hg19 -> hg38
CHAIN=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl/hg19ToHg38.over.chain
if [[ ! -s "$CHAIN" ]]; then
    echo "ERROR: chain file not found: $CHAIN" >&2
    exit 1
fi

if [[ ! -s "$H3K27AC_HG38" || ! -s "$H3K4ME3_HG38" ]]; then
    echo "[$(date)] Lifting H3K27ac and H3K4me3 hg19 -> hg38..."
    eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
    micromamba activate rnaseq
    Rscript --vanilla - <<'EOF'
suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
})
outdir <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/roadmap_liver_E066"
chain_f <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl/hg19ToHg38.over.chain"
chain <- import.chain(chain_f)

lift_peak <- function(hg19_f, hg38_f, label) {
    if (file.exists(hg38_f) && file.info(hg38_f)$size > 0) {
        cat(sprintf("%s already lifted: %s\n", label, hg38_f)); return(invisible())
    }
    peaks <- read.table(gzfile(hg19_f), sep="\t", header=FALSE, stringsAsFactors=FALSE)
    n_cols <- ncol(peaks)
    base_cols <- c("chrom","start","end","name","score","strand","signalValue","pValue","qValue","peak")
    colnames(peaks) <- base_cols[seq_len(n_cols)]
    cat(sprintf("%s hg19 peaks: %d (cols=%d)\n", label, nrow(peaks), n_cols))
    gr <- GRanges(seqnames=peaks$chrom,
                  ranges=IRanges(start=peaks$start+1L, end=peaks$end),
                  strand="*",
                  name=peaks$name, score=peaks$score,
                  signalValue=peaks$signalValue, pValue=peaks$pValue, qValue=peaks$qValue)
    gr_hg38 <- unlist(liftOver(gr, chain))
    cat(sprintf("%s hg38 peaks after liftOver: %d (loss %.2f%%)\n",
                label, length(gr_hg38), 100*(1-length(gr_hg38)/length(gr))))
    df <- data.frame(chrom=as.character(seqnames(gr_hg38)),
                     start=start(gr_hg38)-1L, end=end(gr_hg38),
                     name=gr_hg38$name, score=gr_hg38$score, strand="+",
                     signalValue=gr_hg38$signalValue,
                     pValue=gr_hg38$pValue, qValue=gr_hg38$qValue,
                     stringsAsFactors=FALSE)
    df <- df[df$chrom %in% paste0("chr", c(1:22,"X","Y")),]
    cat(sprintf("%s hg38 after chr filter: %d\n", label, nrow(df)))
    gz <- gzfile(hg38_f, "w")
    write.table(df, gz, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    close(gz)
    cat("Wrote:", hg38_f, "\n")
}

lift_peak(file.path(outdir, "E066_H3K27ac_hg19.broadPeak.gz"),
          file.path(outdir, "E066_H3K27ac_hg38.broadPeak.gz"),
          "H3K27ac")
lift_peak(file.path(outdir, "E066_H3K4me3_hg19.narrowPeak.gz"),
          file.path(outdir, "E066_H3K4me3_hg38.narrowPeak.gz"),
          "H3K4me3")
EOF
else
    echo "H3K27ac + H3K4me3 already lifted to hg38"
fi

# Summary
echo ""
echo "Roadmap E066 files in $OUTDIR:"
ls -lh "$OUTDIR"
echo "[$(date)] Done."
