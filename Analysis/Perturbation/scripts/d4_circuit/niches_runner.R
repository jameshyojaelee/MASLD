#!/usr/bin/env Rscript
# D4 NICHES runner (R). Subagent: circuit-runner-niches.
#
# Called via subprocess from a thin Python wrapper to keep the arm pipeline
# uniform. Reads hit-list CSV, runs NICHES on hepatocyte + receiver atlases,
# emits a JSON matching ModelRunOutput schema with CircuitPrediction entries.
#
# USAGE:
#   Rscript niches_runner.R --modality zero_shot --context all \
#                           --hits ../../data/hits/d4_circuit_ligands.csv \
#                           --out-dir ../../results/d4_circuit
#
# TODO[circuit-runner-niches] blocks below.

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
})

option_list <- list(
    make_option("--modality", type = "character", default = "zero_shot"),
    make_option("--context", type = "character", default = "all"),
    make_option(
        "--hits", type = "character",
        default = "../../data/hits/d4_circuit_ligands.csv"
    ),
    make_option(
        "--out-dir", type = "character",
        default = "../../results/d4_circuit"
    ),
    make_option("--dry-run", action = "store_true", default = FALSE),
    make_option("--max-hits", type = "integer", default = NA)
)
opt <- parse_args(OptionParser(option_list = option_list))

hits <- read.csv(opt$hits, stringsAsFactors = FALSE)
if (!is.na(opt$`max-hits`)) hits <- head(hits, opt$`max-hits`)
message(sprintf("[niches] loaded %d hits", nrow(hits)))

# TODO[circuit-runner-niches]: load hep + receiver Seurat objects and run NICHES.
# library(NICHES); niches_obj <- RunNICHES(seurat_obj, ...)
# For each (ligand, receiver_cell_type) pair, extract receiver response genes.

predictions <- vector("list", nrow(hits))
for (i in seq_len(nrow(hits))) {
    row <- hits[i, ]
    predictions[[i]] <- list(
        sender_gene = row$ligand,
        receiver_cell_type = as.character(row$receiver_cell_type),
        receptor_gene = if ("receptor" %in% names(row)) row$receptor else "—",
        receiver_response_genes = list("R_PLACEHOLDER"),
        response_magnitude = 0.0,
        propagation_method = "niches"
    )
}

envelope <- list(
    arm = "D4",
    model = "niches",
    modality = opt$modality,
    context = opt$context,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    checkpoint_hash = "niches_R",
    n_predictions = length(predictions),
    predictions = predictions,
    runtime_seconds = 0.0,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", NA),
    notes = "D4 circuit / NICHES (R)"
)

out_dir <- opt$`out-dir`
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(
    out_dir,
    sprintf("niches_%s_%s.json", opt$modality, opt$context)
)
write_json(envelope, out_path, auto_unbox = TRUE, pretty = TRUE)
message(sprintf("[niches] wrote %s", out_path))
