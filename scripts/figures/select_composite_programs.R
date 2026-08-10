#!/usr/bin/env Rscript
# select_composite_programs.R — HISTORICAL exploratory row-selection producer.
# The main Figure 4C row identities are now versioned at
# Analysis/Multimodal_Program_Projection/config/panel4c_fixed_rows.tsv and this
# script must not be rerun during production. It is retained only for provenance.
# The original data-driven selection replaced hand-curated candidate lists + top-6
# by |protein_logFC| rule, which admitted non-significant (GSTM1 padj 0.56, FASN 0.055) and
# direction-discordant (NID2 mRNA-/protein+) picks.
#
# RULE (per program): candidates = anchor-set members measured in >=15 PXD051911 patients AND
#   protein padj < 0.05 AND mRNA-protein direction-concordant; rank by |protein_t| (significance-aware
#   moderated effect; |logFC| tiebreak); take top 5. Ranking on |t| rather than |logFC| demotes
#   barely-significant large-|logFC| noise (e.g. CKB padj~0.05) in favour of robustly validated proteins.
#   GLOBAL de-duplication in priority order (up-block programs first) so each gene occupies one row.
# Anchors are MSigDB Reactome/Hallmark sets (reproducible) except the lipid-droplet program, which is a
# curated cell-structural set (lipid droplet is an organelle, not a metabolic pathway).
#
# Mito/OXPHOS is included in the ENRICHMENT table only, NOT as gene rows: it is the 2nd-strongest
# pathway-level signal (NES ~ -2.55) but a coordinated small-effect shift across ETC subunits with no
# individually significant, well-quantified, concordant protein representatives in this DIA-MS — forcing
# gene rows yields up/discordant picks. It is shown honestly as a pathway-level NES bar instead.
#
# Historical outputs: FIG4_DIR/data/{composite_program_genes,composite_program_enrichment}.csv.
# The production composite no longer reads either file.
suppressPackageStartupMessages({ library(data.table); library(msigdbr); library(fgsea) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
if (toupper(Sys.getenv("ALLOW_PANEL4C_RESELECTION", "FALSE")) != "TRUE") {
  stop(
    "Panel 4C selection is frozen. Production reads ",
    "Analysis/Multimodal_Program_Projection/config/panel4c_fixed_rows.tsv. ",
    "Set ALLOW_PANEL4C_RESELECTION=TRUE only for an explicitly versioned exploratory reselection."
  )
}
DATA_DIR <- file.path(BASE, "figures/main/fig4_validation/data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
N_PER <- 5L

# ── data: PXD051911 protein DE + mRNA concordance, and coverage from the abundance matrix ─────────────
con <- fread(file.path(BASE, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv"))
con <- unique(con[dataset == "PXD051911" & is.finite(protein_logFC)], by = "gene")
ab <- fread(file.path(BASE, "data/PXD051911/liver_protein_quant.txt"))
ab <- ab[Genes != "" & !is.na(Genes)]; ab[, gene := tstrsplit(Genes, ";", fixed = TRUE)[[1]]]
ab <- ab[!duplicated(gene)]
scols <- setdiff(names(ab), c("ProteinAccessions", "Genes", "ProteinDescriptions", "gene"))
abmat <- as.matrix(ab[, ..scols]); rownames(abmat) <- ab$gene; mode(abmat) <- "numeric"
covered <- rownames(abmat)[rowSums(!is.na(abmat)) >= 15]

R <- function(gs) unique(as.data.table(msigdbr(species = "Homo sapiens", collection = "C2",
                                               subcollection = "CP:REACTOME"))[gs_name == gs, gene_symbol])
H <- function(gs) unique(as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))[gs_name == gs, gene_symbol])

# ── program anchors. First 5 = gene-row programs (up-block then down-block); Mito = enrichment-only ────
anchors <- list(
  "Lipid droplet" = c("PLIN1","PLIN2","PLIN3","PLIN4","PLIN5","ABHD5","CIDEA","CIDEB","CIDEC","G0S2",
                       "HILPDA","ACLY","FASN","ACACA","SCD","DGAT1","DGAT2","LPIN1"),
  "ECM / BM"      = R("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"),
  # Immune infiltration is a cell-population signal, not one pathway -> curated liver MASH immune/
  # inflammation markers (broad Hallmark ALLOGRAFT_REJECTION/IFN-gamma sets contaminate with housekeeping
  # genes e.g. SEC24C/DCTN4 that rank high on |t| but are not immune biology).
  "Immune / inflammation" = c("HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","CD74",
                       "CD68","CD163","AIF1","TYROBP","FCER1G","LGALS3","CTSS","LAPTM5","MARCO","VSIG4",
                       "CAPG","SPP1","CASP1","PYCARD","S100A8","S100A9","PTPRC","LCP1","CORO1A","ITGAL",
                       "C1QA","C1QB","C1QC","C3"),
  # Detox = xenobiotic/drug metabolism: Phase I + Phase II only (BIOLOGICAL_OXIDATIONS also pulls in
  # methylation/one-carbon enzymes e.g. MAT1A/NNMT that belong with amino-acid metabolism).
  "Detox"         = unique(c(R("REACTOME_PHASE_I_FUNCTIONALIZATION_OF_COMPOUNDS"),
                             R("REACTOME_PHASE_II_CONJUGATION_OF_COMPOUNDS"))),
  "AA catabolism" = R("REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES"),
  "Mito / OXPHOS" = unique(c(R("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"),
                             R("REACTOME_COMPLEX_I_BIOGENESIS"),
                             R("REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE"))))
row_programs <- names(anchors)[1:5]

# ── selection (5 gene-row programs, global dedup) ─────────────────────────────────────────────────────
used <- character(0); picks <- list()
for (i in seq_along(row_programs)) {
  p <- row_programs[i]
  d <- con[gene %in% anchors[[p]] & gene %in% covered & !gene %in% used &
           is.finite(protein_padj) & protein_padj < 0.05]
  d[, conc := is.finite(bulk_logFC) & sign(bulk_logFC) == sign(protein_logFC)]
  d <- d[conc == TRUE][order(-abs(protein_t), -abs(protein_logFC))]
  sel <- head(d, N_PER); used <- c(used, sel$gene)
  # Display order: up-regulated proteins first (then |t|), so a direction-discordant member of a
  # down-program (GGT1 in Detox) rises to the block TOP, extending the contiguous up-block formed by
  # the up-programs (Lipid/ECM/Immune) directly above it. Gene SELECTION above is still |t|-ranked;
  # this only re-sorts the rows for the shared vertical axis. Mono-direction programs are unaffected.
  sel[, .updown := as.integer(protein_logFC < 0)]           # up=0 before down=1
  setorder(sel, .updown)                                     # radix sort is stable -> |t| order preserved within each group
  sel[, .updown := NULL]
  sel[, `:=`(program = p, prog_order = i, gene_order = seq_len(.N))]
  picks[[p]] <- sel[, .(program, prog_order, gene_order, gene,
                        protein_logFC = round(protein_logFC, 4), protein_t = round(protein_t, 3),
                        protein_padj = signif(protein_padj, 4), bulk_logFC = round(bulk_logFC, 4))]
  cat(sprintf("%-14s pool(measured,sig,conc,unused)=%2d  picked=%d: %s\n",
              p, nrow(d), nrow(sel), paste(sel$gene, collapse = ", ")))
}
genes_dt <- rbindlist(picks)
stopifnot(nrow(genes_dt) == length(row_programs) * N_PER, !anyDuplicated(genes_dt$gene))
fwrite(genes_dt, file.path(DATA_DIR, "composite_program_genes.csv"))

# ── pathway-level enrichment for all 6 programs (incl. Mito) — for the honest NES strip ────────────────
de <- fread(file.path(BASE, "Analysis/Proteomics/results/protein_differential_results_v3.csv"))
de <- de[dataset == "PXD051911" & is.finite(t) & gene != ""][order(-abs(t))][!duplicated(gene)]
rk <- sort(setNames(de$t, de$gene), decreasing = TRUE)
sets <- lapply(anchors, function(s) intersect(s, names(rk)))
fg <- fgsea(pathways = sets, stats = rk, minSize = 8, maxSize = 600, eps = 0)
enr <- fg[, .(program = pathway, NES = round(NES, 3), padj = signif(padj, 4), size,
              direction = ifelse(NES > 0, "up", "down"),
              representable = pathway %in% row_programs)]
enr <- enr[match(names(anchors), program)]           # preserve program order (Mito last)
fwrite(enr, file.path(DATA_DIR, "composite_program_enrichment.csv"))
cat("\n===== pathway-level NES (all 6 programs) =====\n"); print(enr)
cat("\n[select] wrote:", file.path(DATA_DIR, "composite_program_genes.csv"), "and composite_program_enrichment.csv\n")
