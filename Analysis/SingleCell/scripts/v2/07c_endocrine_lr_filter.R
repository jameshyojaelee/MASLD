#!/usr/bin/env Rscript
# ============================================================================
# 07c_endocrine_lr_filter.R
#
# Addresses critique #6: SAA1-class endocrine "LR pairs" dominate top hits.
# SAA1 reaches receptors via bloodstream, not paracrine signaling. LIANA's
# rank_aggregate cannot distinguish hepatocyte-secreted-to-blood ligands from
# locally-acting paracrine ligands.
#
# This script builds a per-LR-pair `endocrine_suspect` boolean from:
#   (1) OmniPath UniProt_topology annotations: ligands tagged `Secreted` AND
#       known to circulate in plasma.
#   (2) Manually-curated hepatokine + acute-phase + classical-blood-protein
#       list (overrides #1 in case OmniPath query fails).
#
# Inputs:
#   - Either: a merged per-donor LR table (default:
#       results_gpu_v2_phase05/ccc/stage_trajectory_v2/all_donor_lr_scores_v2.tsv.gz)
#     Reads only the unique (ligand_complex, receptor_complex) combinations.
#
# Output:
#   - results_gpu_v2_phase05/ccc/stage_trajectory_v3/endocrine_lr_flags.tsv
#     Columns: ligand_complex, receptor_complex, lr_pair, ligand,
#              ligand_is_secreted_omnipath, ligand_is_hepatokine_curated,
#              endocrine_suspect
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

V2_STAGE_DIR <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2"
)
V3_STAGE_DIR <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3"
)
dir.create(V3_STAGE_DIR, recursive = TRUE, showWarnings = FALSE)

MERGED_TSV <- file.path(V2_STAGE_DIR, "all_donor_lr_scores_v2.tsv.gz")
OUT_FLAGS <- file.path(V3_STAGE_DIR, "endocrine_lr_flags.tsv")

# ----------------------------------------------------------------------------
# (1) Manually-curated hepatokine + acute-phase + classical-blood-protein list
#     Sourced from canonical hepatokine reviews (e.g. Stefan & Häring 2013,
#     Meex & Watt 2017) + Kleiner-2005 acute phase markers + standard
#     liver-secreted plasma proteins (coagulation factors, complement,
#     transport proteins, lipoproteins, lipid-handling enzymes).
# ----------------------------------------------------------------------------
HEPATOKINE_CURATED <- c(
  # Acute phase proteins (positive APR)
  "SAA1", "SAA2", "SAA4", "CRP", "AGT", "HP", "FGA", "FGB", "FGG",
  "AHSG", "SERPINA1", "SERPINA3", "SERPING1", "ORM1", "ORM2", "F2",
  # Complement components + regulators (all hepatocyte-secreted)
  "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2",
  "C3", "C4A", "C4B", "C4BPA", "C4BPB", "C5", "C6", "C7",
  "C8A", "C8B", "C8G", "C9",
  "CFB", "CFD", "CFH", "CFI", "CFP", "MBL2",
  # Coagulation factors
  "F2", "F5", "F7", "F8", "F9", "F10", "F11", "F12", "F13A1", "F13B",
  "PROC", "PROS1", "PLG", "SERPINF2", "SERPINF1",
  # Lipoproteins
  "APOA1", "APOA2", "APOA4", "APOA5", "APOB", "APOC1", "APOC2", "APOC3",
  "APOC4", "APOD", "APOE", "APOF", "APOH", "APOL1", "APOM",
  "CETP", "LCAT", "LPL", "LIPC", "LIPG",
  # Transport proteins
  "ALB", "TF", "TTR", "RBP4", "GC", "SHBG", "TBG", "TBPA",
  "VTN", "CP", "LBP", "AFP", "A2M", "AMBP", "AOAH",
  # Inter-alpha-trypsin-inhibitor family (hepatocyte-secreted, plasma)
  "ITIH1", "ITIH2", "ITIH3", "ITIH4", "ITIH5",
  # Insulin-like growth factor / IGFBP family (hepatic-dominant secretion)
  "IGF1", "IGF2", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5",
  # Hepatokines / metabolism-modulating
  "ANGPTL3", "ANGPTL4", "ANGPTL6", "ANGPTL8", "FGF21", "FGF19",
  "HSDL2", "FETUB", "LECT2", "SELENOP", "RETN", "ADIPOQ",
  # Liver enzymes that leak / are secreted
  "PON1", "PON2", "PON3",
  # Hepcidin / iron
  "HAMP", "TFRC",
  # Kallikrein / kininogen
  "KNG1", "KLKB1",
  # MASP family
  "MASP1", "MASP2",
  # Bilirubin / heme transport
  "HPX", "HRG"
)
HEPATOKINE_CURATED <- sort(unique(HEPATOKINE_CURATED))
cat(sprintf("[curated] %d hepatokine/acute-phase/blood-protein genes\n",
            length(HEPATOKINE_CURATED)))

# ----------------------------------------------------------------------------
# (2) Try OmniPath UniProt_topology annotations (Secreted topology label).
#     Fails-soft to a NULL with a clear warning so the curated list still runs.
# ----------------------------------------------------------------------------
fetch_omnipath_secreted <- function() {
  if (!requireNamespace("OmnipathR", quietly = TRUE)) {
    cat("[omnipath] OmnipathR not installed; skipping\n")
    return(NULL)
  }
  ann <- tryCatch(
    {
      OmnipathR::import_omnipath_annotations(
        resources = "UniProt_topology"
      )
    },
    error = function(e) {
      cat(sprintf("[omnipath] query failed: %s\n", e$message))
      NULL
    }
  )
  if (is.null(ann)) return(NULL)
  ann <- as.data.table(ann)
  if (!all(c("genesymbol", "label", "value") %in% names(ann))) {
    cat("[omnipath] unexpected schema; columns:",
        paste(names(ann), collapse = ", "), "\n")
    return(NULL)
  }
  secreted_genes <- unique(ann[label == "topology" & value == "Secreted",
                               genesymbol])
  cat(sprintf("[omnipath] %d genes with Secreted topology\n",
              length(secreted_genes)))
  secreted_genes
}

secreted_omni <- fetch_omnipath_secreted()
if (is.null(secreted_omni)) {
  cat("[omnipath] using curated hepatokine list as sole signal\n")
  secreted_omni <- character(0)
}

# ----------------------------------------------------------------------------
# (3) Build per-LR-pair flags
# ----------------------------------------------------------------------------
if (!file.exists(MERGED_TSV)) {
  stop("MERGED_TSV not found at ", MERGED_TSV,
       " -- run 07_chain_346_349_v2 (or 07b_chain_v3_hardened) first to produce it")
}

cat(sprintf("[input] reading %s\n", MERGED_TSV))
lr_long <- fread(MERGED_TSV,
                 select = c("ligand_complex", "receptor_complex"))
lr_unique <- unique(lr_long)
cat(sprintf("[input] %d unique (ligand_complex, receptor_complex)\n",
            nrow(lr_unique)))

# Decompose complex ligands into all member genes -- if ANY component is
# endocrine, flag the whole complex.
split_complex <- function(s) {
  if (is.na(s) || s == "") return(character(0))
  unlist(strsplit(as.character(s), "_", fixed = TRUE))
}

flag_one <- function(lc) {
  members <- split_complex(lc)
  is_secreted <- any(members %in% secreted_omni)
  is_hepatokine <- any(members %in% HEPATOKINE_CURATED)
  list(
    ligand = paste(members, collapse = "+"),
    ligand_is_secreted_omnipath = is_secreted,
    ligand_is_hepatokine_curated = is_hepatokine
  )
}

lr_unique[, c("ligand", "ligand_is_secreted_omnipath",
              "ligand_is_hepatokine_curated") :=
            rbindlist(lapply(ligand_complex, flag_one))]
lr_unique[, endocrine_suspect := ligand_is_secreted_omnipath |
                                 ligand_is_hepatokine_curated]
lr_unique[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]

cat("\n[flags] breakdown:\n")
cat(sprintf("  total unique LR pairs:               %d\n", nrow(lr_unique)))
cat(sprintf("  ligand in OmniPath Secreted:         %d\n",
            sum(lr_unique$ligand_is_secreted_omnipath, na.rm = TRUE)))
cat(sprintf("  ligand in curated hepatokine list:   %d\n",
            sum(lr_unique$ligand_is_hepatokine_curated, na.rm = TRUE)))
cat(sprintf("  endocrine_suspect (either):          %d  (%.1f%%)\n",
            sum(lr_unique$endocrine_suspect, na.rm = TRUE),
            100 * mean(lr_unique$endocrine_suspect, na.rm = TRUE)))

# Sanity log: top hepatokine ligands present in our LR list
top_hepatokine <- lr_unique[ligand_is_hepatokine_curated == TRUE,
                            .N, by = ligand][order(-N)]
cat("\n[sanity] top 10 hepatokine ligands appearing in LR pairs:\n")
print(head(top_hepatokine, 10))

setcolorder(lr_unique,
            c("lr_pair", "ligand_complex", "receptor_complex", "ligand",
              "ligand_is_secreted_omnipath", "ligand_is_hepatokine_curated",
              "endocrine_suspect"))
fwrite(lr_unique, OUT_FLAGS, sep = "\t")
cat(sprintf("\n[output] %s (%d rows)\n", OUT_FLAGS, nrow(lr_unique)))
cat("[done]\n")
