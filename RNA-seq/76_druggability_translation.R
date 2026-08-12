#!/usr/bin/env Rscript
# ============================================================================
# 76_druggability_translation.R
# Tier-1 atlas layer: druggability / tractability / clinical-translation bundle.
#
# Produces gene-level columns (keyed on human_symbol) that reframe the drug
# section from "0/10 pass our strict criteria" into a transparent, externally
# anchored druggability ladder: clinical phase x target development level x
# genetic support x population constraint.
#
# Sources (all public, no COLOC, no LD):
#   - gnomAD v4.1 constraint    -> gnomad_loeuf, gnomad_mis_z, gnomad_loeuf_constrained
#   - Pharos / TCRD (GraphQL)    -> pharos_tdl (Tclin/Tchem/Tbio/Tdark), pharos_fam, pharos_novelty
#   - Curated MASH/MASLD pipeline (ClinicalTrials.gov-anchored gene<->drug map)
#                                -> clintrial_max_phase, clintrial_status,
#                                   clintrial_active, clintrial_n_drugs,
#                                   clintrial_drugs, clintrial_mechanism
#
# The derived `genetic_support` flag + Minikel-2024 framing is computed in
# 27a_assemble_evidence_atlas.R (post-assembly block), where the COLOC/TWAS
# columns already exist. This script writes only source-derived annotations.
#
# Output: RNA-seq/results/multi_evidence/druggability_atlas_columns.tsv
# Merged by 27a's post-assembly block (drop-then-merge on human_symbol).
# ============================================================================

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
EXT  <- file.path(BASE, "data/external")
OUT  <- file.path(BASE, "RNA-seq/results/multi_evidence/druggability_atlas_columns.tsv")

cat("== 76_druggability_translation ==\n")

# ---------------------------------------------------------------------------
# 1. gnomAD v4.1 constraint (LOEUF = lof.oe_ci.upper; missense Z) -- safety axis
#    One row per gene via the MANE Select transcript.
# ---------------------------------------------------------------------------
gn_path <- file.path(EXT, "gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv")
gn <- fread(gn_path, select = c("gene", "mane_select",
                                "lof.oe_ci.upper", "mis.z_score"))
# fread auto-detects "true"/"false" as logical; coerce robustly either way.
gn <- gn[as.logical(mane_select) %in% TRUE]
gn <- gn[!duplicated(gene)]
gn <- gn[, .(human_symbol = gene,
             gnomad_loeuf = suppressWarnings(as.numeric(`lof.oe_ci.upper`)),
             gnomad_mis_z = suppressWarnings(as.numeric(`mis.z_score`)))]
# LOEUF < 0.6 = LoF-intolerant -> inhibiting the target carries higher a-priori
# safety risk (gnomAD constraint working group convention).
gn[, gnomad_loeuf_constrained := !is.na(gnomad_loeuf) & gnomad_loeuf < 0.6]
cat(sprintf("  gnomAD: %d MANE genes (LOEUF non-NA: %d)\n",
            nrow(gn), sum(!is.na(gn$gnomad_loeuf))))

# ---------------------------------------------------------------------------
# 2. Pharos / TCRD Target Development Level (Tclin > Tchem > Tbio > Tdark)
# ---------------------------------------------------------------------------
ph_path <- file.path(EXT, "druggability/pharos_tdl.tsv")
if (file.exists(ph_path)) {
  ph <- fread(ph_path)
  ph <- ph[!duplicated(pharos_sym) & pharos_sym != ""]
  ph <- ph[, .(human_symbol = pharos_sym,
               pharos_tdl, pharos_fam,
               pharos_novelty = suppressWarnings(as.numeric(pharos_novelty)))]
  cat(sprintf("  Pharos: %d targets; TDL = %s\n", nrow(ph),
              paste(sprintf("%s:%d", names(table(ph$pharos_tdl)),
                            as.integer(table(ph$pharos_tdl))), collapse = " ")))
} else {
  cat("  ! Pharos TDL file missing; skipping (run fetch_pharos_tdl.py)\n")
  ph <- data.table(human_symbol = character(), pharos_tdl = character(),
                   pharos_fam = character(), pharos_novelty = numeric())
}

# ---------------------------------------------------------------------------
# 3. Curated MASH/MASLD clinical pipeline (gene <-> drug, ClinicalTrials-anchored)
#    Target<->drug mapping is stable knowledge; auto-mapping CT.gov interventions
#    to gene targets is unreliable, so a curated map is the correct approach.
# ---------------------------------------------------------------------------
ct_path <- file.path(EXT, "druggability/mash_clinical_pipeline.tsv")
ct <- fread(ct_path, fill = TRUE)  # some rows omit the trailing `note` column
ct_agg <- ct[, .(
  clintrial_max_phase = max(max_phase, na.rm = TRUE),
  clintrial_status    = fifelse(any(status == "approved"), "approved",
                         fifelse(any(status == "active"),   "active", "discontinued")),
  clintrial_n_drugs   = uniqueN(drug),
  clintrial_drugs     = paste(unique(drug), collapse = ";"),
  clintrial_mechanism = paste(unique(mechanism), collapse = ";")
), by = .(human_symbol = gene_symbol)]
ct_agg[, clintrial_active := clintrial_status %in% c("approved", "active")]
cat(sprintf("  ClinicalTrials map: %d target genes (approved=%d active=%d disc=%d)\n",
            nrow(ct_agg),
            sum(ct_agg$clintrial_status == "approved"),
            sum(ct_agg$clintrial_status == "active"),
            sum(ct_agg$clintrial_status == "discontinued")))

# ---------------------------------------------------------------------------
# 4. Merge (union of gnomAD + Pharos; clinical map left-joined)
# ---------------------------------------------------------------------------
out <- merge(gn, ph, by = "human_symbol", all = TRUE)
out <- merge(out, ct_agg, by = "human_symbol", all.x = TRUE)
out[is.na(clintrial_active), clintrial_active := FALSE]
out[is.na(clintrial_max_phase), clintrial_max_phase := NA_real_]

# ---------------------------------------------------------------------------
# 4b. Systematic drug-development-status classification (Wave-4 reframe).
#     Genome-wide drug_dev_status (OpenTargets + DGIdb + Pharos + CT.gov +
#     curated MASH pipeline). MASLD-specific status comes from curated + CT.gov.
#     OUTER-joined so every classified gene is represented in the output.
#     See docs/archive/relocated_2026-07-08/manuscript/working/drug_dev_status_reframe_spec.md.
# ---------------------------------------------------------------------------
dds_path <- file.path(EXT, "drug_targets/drug_target_classification.tsv")
if (file.exists(dds_path)) {
  dds <- fread(dds_path, select = c("symbol", "drug_dev_status",
                                    "max_phase_any", "max_phase_masld",
                                    "has_masld_trial", "masld_preclinical_evidence"))
  dds <- dds[symbol != "" & !duplicated(symbol)]
  setnames(dds, "symbol", "human_symbol")
  out <- merge(out, dds, by = "human_symbol", all = TRUE)  # OUTER: keep all classified genes
  cat(sprintf("  drug_dev_status: %d classified genes joined; status = %s\n",
              nrow(dds),
              paste(sprintf("%s:%d", names(table(dds$drug_dev_status)),
                            as.integer(table(dds$drug_dev_status))), collapse = " ")))
} else {
  cat("  ! drug_target_classification.tsv missing; skipping drug_dev_status block\n")
}

out[is.na(clintrial_active), clintrial_active := FALSE]
setorder(out, human_symbol)

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUT, sep = "\t")
cat(sprintf("WROTE %s : %d genes x %d cols\n", OUT, nrow(out), ncol(out)))

# Sanity anchors
for (g in c("THRB", "PNPLA3", "GLP1R", "FGF21", "NR1H4")) {
  r <- out[human_symbol == g]
  if (nrow(r))
    cat(sprintf("  [%s] TDL=%s LOEUF=%.2f phase=%s status=%s drugs=%s\n",
                g, r$pharos_tdl[1], r$gnomad_loeuf[1],
                as.character(r$clintrial_max_phase[1]), r$clintrial_status[1],
                substr(r$clintrial_drugs[1], 1, 40)))
}

# drug_dev_status sanity anchor (Wave-4 reframe)
if ("drug_dev_status" %in% names(out)) {
  for (g in c("THRB", "RORA", "HKDC1", "NR1H4", "GLP1R")) {
    r <- out[human_symbol == g]
    if (nrow(r))
      cat(sprintf("  [%s] drug_dev_status=%s max_phase_any=%s max_phase_masld=%s\n",
                  g, r$drug_dev_status[1],
                  as.character(r$max_phase_any[1]),
                  as.character(r$max_phase_masld[1])))
  }
}
