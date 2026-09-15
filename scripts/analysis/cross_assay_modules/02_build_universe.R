#!/usr/bin/env Rscript
# 02: the shared gene universe.
#
# This is the step that makes the family testable everywhere. A gene enters the
# universe only if bulk RNA, the single-nucleus atlas, both Visium sources and
# the liver proteome all measure it. Every later "untestable" cell is therefore
# a statement about a module's members in one assay, not about the design.
#
# The universe is built with no outcome in hand. Nothing below reads a stage,
# a NAS score or a diagnosis.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))

contract <- cam_contract()
out <- cam_dir("universe")
cam_assert(file.exists(file.path(cam_out_root(), "inputs", "READY")),
           "Run 01_freeze_inputs.R first")

# --- bulk RNA ---------------------------------------------------------------
# The tested set of the five-cohort substrate, in the symbol space every other
# assay uses. The ENSG-to-symbol pairing is taken from the substrate's own
# results file so it cannot drift from the matrix the modules are built on.
cam_say("reading bulk gene universe")
bulk <- unique(fread(cam_input("bulk_stage_results", contract),
                     select = c("gene_id_versioned", "gene_name")))
cam_assert(nrow(bulk) == 23370L,
           paste0("Expected 23370 bulk genes, got ", nrow(bulk)))
bulk_symbols <- unique(bulk[!is.na(gene_name) & gene_name != "", gene_name])
cam_say("bulk symbols: ", length(bulk_symbols))

# --- snRNA atlas ------------------------------------------------------------
# All-cells pseudobulk is the sum over lineages of the same run columns. A gene
# counts as measured only if it is non-zero on average in at least two datasets,
# so a single deep library cannot admit a gene the atlas does not really see.
cam_say("reading snRNA pseudobulk")
pb_dir <- cam_input("snrna_pseudobulk_dir", contract)
pb_files <- list.files(pb_dir, pattern = "^pseudobulk_.*\\.tsv\\.gz$", full.names = TRUE)
cam_assert(length(pb_files) >= 5L, "Too few pseudobulk lineage files")
sn_meta <- fread(cam_input("snrna_metadata", contract))
# Lineage files do not share a run roster: a run with no B cells has no B-cell
# column. Accumulating onto the first file's columns would silently reduce the
# atlas to the runs of its smallest lineage, so the target matrix is the UNION
# of run columns and each lineage adds into its own.
run_union <- sort(unique(unlist(lapply(pb_files, function(f) {
  names(data.table::fread(f, nrows = 0L))[-1L]
}))))
gene_ref <- data.table::fread(pb_files[[1L]], select = 1L)[[1L]]
all_cells <- matrix(0, nrow = length(gene_ref), ncol = length(run_union),
                    dimnames = list(gene_ref, run_union))
for (f in pb_files) {
  x <- fread(f)
  setnames(x, 1L, "gene")
  cam_assert(identical(x$gene, gene_ref),
             paste0("Pseudobulk gene order differs in ", basename(f)))
  m <- as.matrix(x[, -1L]); rownames(m) <- x$gene
  all_cells[, colnames(m)] <- all_cells[, colnames(m)] + m
}
cam_say("all-cells pseudobulk: ", nrow(all_cells), " genes x ", ncol(all_cells), " runs")
run_dataset <- sn_meta[match(colnames(all_cells), sample), dataset]
cam_assert(sum(is.na(run_dataset)) == 0L, "Some pseudobulk runs have no dataset label")
datasets <- unique(run_dataset)
detected_per_dataset <- vapply(datasets, function(ds) {
  cols <- which(run_dataset == ds)
  rowMeans(all_cells[, cols, drop = FALSE]) > 0
}, logical(nrow(all_cells)))
sn_detected <- rowSums(detected_per_dataset) >= contract$universe$snrna_min_datasets_detected
snrna_symbols <- rownames(all_cells)[sn_detected]
cam_say("snRNA symbols detected in >=2 datasets: ", length(snrna_symbols))
cam_write_rds(all_cells, file.path(out, "snrna_all_cells_pseudobulk.rds"))
cam_write_tsv(data.table(run = colnames(all_cells), dataset = run_dataset),
              file.path(out, "snrna_run_dataset.tsv"))

# --- Visium -----------------------------------------------------------------
# Both sources must measure the gene. One-source coverage would make the
# spatial column depend on which section happened to be sequenced deeper.
cam_say("reading Visium feature lists")
vis_g <- readLines(file.path(out, "measured_visium_guilliams.txt"))
vis_v <- readLines(file.path(out, "measured_visium_vu.txt"))
visium_symbols <- intersect(vis_g, vis_v)
cam_say("Visium symbols in both sources: ", length(visium_symbols))

# --- liver proteome ---------------------------------------------------------
# Detection, not imputation. A protein quantified in under 70 percent of the 58
# participants cannot order them, and filling the gap would invent the ordering
# the test is about.
cam_say("reading liver proteome")
prot <- fread(cam_input("proteome_liver", contract))
prot_vals <- as.matrix(prot[, 4:ncol(prot)])
detect_frac <- rowMeans(!is.na(prot_vals))
prot_genes <- vapply(strsplit(prot$Genes, ";", fixed = TRUE), function(g) {
  g <- trimws(g[nzchar(trimws(g))]); if (length(g)) g[1] else NA_character_
}, character(1))
keep_prot <- detect_frac >= contract$universe$proteome_min_detection_fraction &
  !is.na(prot_genes)
proteome_symbols <- unique(prot_genes[keep_prot])
cam_say("proteome symbols at >=70% detection: ", length(proteome_symbols))

# --- the universe -----------------------------------------------------------
universe <- Reduce(intersect, list(bulk_symbols, snrna_symbols,
                                   visium_symbols, proteome_symbols))
universe <- sort(unique(universe))
cam_say("UNIVERSE SIZE: ", length(universe))
cam_assert(length(universe) >= contract$universe$min_universe_size,
  paste0("Universe of ", length(universe), " is below the prespecified floor of ",
         contract$universe$min_universe_size, "; the analysis stops here"))

cam_write_tsv(data.table(gene_symbol = universe), file.path(out, "universe.tsv"))

# --- coverage of the universe in every assay, including the non-reference ones
cosmx <- readLines(file.path(out, "measured_cosmx_govaere.txt"))
# The plasma deposit carries UniProt accessions only. The released protein
# instrument already resolved accession to symbol against this same cohort, so
# its axis is the crosswalk rather than a second ad hoc parse.
plasma <- fread(cam_input("proteome_plasma", contract), select = "ProteinAccessions")
uni <- fread(cam_input("uniprot_symbol_axis", contract),
             select = c("protein_accession", "gene_symbol"))
acc2sym <- setNames(uni$gene_symbol, uni$protein_accession)
plasma_first_acc <- vapply(strsplit(plasma$ProteinAccessions, ";", fixed = TRUE),
                           function(a) { a <- trimws(a[nzchar(trimws(a))]); if (length(a)) a[1] else NA_character_ },
                           character(1))
plasma_genes <- unique(stats::na.omit(unname(acc2sym[plasma_first_acc])))
cam_say("plasma accessions resolved to symbols: ", length(plasma_genes),
        " of ", nrow(plasma))

coverage <- data.table(
  assay = c("bulk_rna_five_cohort", "snrna_atlas", "visium_guilliams", "visium_vu",
            "proteome_liver", "cosmx_govaere", "proteome_plasma"),
  role = c(rep("universe_reference", 2), rep("universe_reference", 2),
           "universe_reference", "coverage_only", "coverage_only"),
  n_features_measured = c(length(bulk_symbols), length(snrna_symbols),
                          length(vis_g), length(vis_v), length(proteome_symbols),
                          length(cosmx), length(plasma_genes)),
  n_universe_measured = c(
    length(intersect(universe, bulk_symbols)),
    length(intersect(universe, snrna_symbols)),
    length(intersect(universe, vis_g)),
    length(intersect(universe, vis_v)),
    length(intersect(universe, proteome_symbols)),
    length(intersect(universe, cosmx)),
    length(intersect(universe, plasma_genes))
  )
)
coverage[, fraction_universe_measured := n_universe_measured / length(universe)]
cam_assert_no_prohibited_columns(coverage, contract)
cam_write_tsv(coverage, file.path(out, "assay_universe_coverage.tsv"))

cam_write_json(list(
  universe_size = length(universe),
  bulk_symbols = length(bulk_symbols),
  snrna_symbols = length(snrna_symbols),
  visium_intersection = length(visium_symbols),
  proteome_symbols = length(proteome_symbols),
  cosmx_universe_fraction = coverage[assay == "cosmx_govaere", fraction_universe_measured],
  plasma_universe_fraction = coverage[assay == "proteome_plasma", fraction_universe_measured],
  floor = contract$universe$min_universe_size
), file.path(out, "universe_summary.json"))

writeLines("universe built", file.path(out, "READY"))
cam_say("02 complete")
