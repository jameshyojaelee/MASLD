#!/usr/bin/env Rscript
# 15: one row per module per assay, with two evidence states.
#
# This is the only place direction enters and the only place a state is
# assigned. The rule has these properties, each of which blocks a way of
# overclaiming:
#
#   Direction is applied ONCE. Assay scripts save the raw standardized slope;
#   effect_oriented = effect_raw * direction. v1 multiplied by the direction at
#   scoring and again here, which for the 101 decreasing modules turned
#   agreement into the wrong tail.
#
#   Two questions, two states. association_state reads a two-sided BH q and
#   the sign of the oriented effect. specificity_state reads a BH-corrected
#   competitive q. Neither is folded into the other and nothing counts them.
#
#   A negative is a bound, not an absence of significance. tested_negative
#   requires the BH-corrected equivalence p at the prespecified SESOI. Every
#   other non-significant cell prints its one-sided 95 percent upper bound.
#
#   Coverage and applicability are decided before any statistic. An assay
#   without an endpoint of the module's family is not_applicable.
#
#   Caps are stated. A proteome cell cannot replicate unless it survives the
#   ten intensity descriptors; an ATAC column with too few wells carrying both
#   classes cannot replicate; a Vu spatial cell is always source_dependent.
#
# The self-tests run first, and the script refuses to write a row if the rule
# does not behave as documented.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))

contract <- cam_contract()
cam_state_self_test()
cam_say("state-rule self-tests passed")
out <- cam_dir("panel")
A <- file.path(cam_out_root(), "assays")
N <- file.path(cam_out_root(), "nulls")

labels <- fread(file.path(cam_out_root(), "discovery", "module_discovery_labels.tsv"))
registry <- fread(file.path(cam_out_root(), "modules", "module_registry.tsv"))
module_ids <- labels$module_id
n_family <- length(module_ids)
direction <- setNames(labels$direction, labels$module_id)[module_ids]
family_of <- setNames(labels$endpoint_family, labels$module_id)[module_ids]
fam_map <- contract$endpoint_families
proxy_assays <- names(contract$endpoint_families$proxy_endpoints)
sesoi <- contract$negatives$sesoi_standardized_slope
sesoi2 <- contract$negatives$sesoi_sensitivity

read_if <- function(f) if (file.exists(file.path(A, f))) fread(file.path(A, f)) else NULL
nulls <- if (file.exists(file.path(N, "competitive_nulls_per_assay.tsv")))
  fread(file.path(N, "competitive_nulls_per_assay.tsv")) else NULL
cal <- if (file.exists(file.path(N, "competitive_calibration.json")))
  jsonlite::fromJSON(file.path(N, "competitive_calibration.json")) else list(calibrated = NA)
spec_descriptive_only <- !isTRUE(cal$calibrated)
atac_summary <- jsonlite::fromJSON(file.path(A, "atac_summary.json"))
atac_ok <- atac_summary$wells_with_both_classes >= contract$assays$atac_gse296875$min_wells_with_both_classes

assay_spec <- list(
  bulk_gse268273 = list(file = "external_bulk_results.tsv", question = "stage association"),
  bulk_gse276114 = list(file = "external_bulk_results.tsv", question = "stage association"),
  proteome_liver_pxd051911 = list(file = "proteome_results.tsv", question = "histology association"),
  snrna_all_cells = list(file = "snrna_results.tsv", question = "stage association"),
  h3k27ac_gse267145 = list(file = "h3k27ac_results.tsv", question = "histology association"),
  atac_gse296875 = list(file = "atac_results.tsv", question = "fibrosis association")
)

na_row <- function() list(testable = NA, n_units = NA_integer_, df = NA_real_, n_members = NA_integer_,
  n_measured = NA_integer_, fraction_measured = NA_real_, beta = NA_real_, se = NA_real_,
  p = NA_real_, endpoint = NA_character_, endpoint_note = NA_character_, unit = NA_character_,
  cx_p = NA_real_, n_amb = NA_integer_, lowo_min = NA_real_, lowo_max = NA_real_)

rows <- list()
for (nm in names(assay_spec)) {
  spec <- assay_spec[[nm]]
  d <- read_if(spec$file)
  if (is.null(d)) { cam_say("missing ", spec$file, "; ", nm, " omitted"); next }
  d <- d[assay == nm]
  cells <- lapply(seq_len(n_family), function(i) {
    mid <- module_ids[i]; fam <- family_of[[i]]
    ep <- fam_map[[fam]][[nm]]
    if (is.null(ep)) {
      # No endpoint of this family here, but the module IS measured (or not);
      # coverage is a property of the assay, not of the question asked of it.
      r0 <- d[module_id == mid][1]
      cov <- na_row()
      cov$testable <- isTRUE(r0$testable); cov$n_members <- r0$n_members
      cov$n_measured <- r0$n_measured; cov$fraction_measured <- r0$fraction_measured
      cov$unit <- r0$unit
      return(c(list(applicable = FALSE), cov))
    }
    r <- d[module_id == mid & endpoint == ep]
    cam_assert(nrow(r) == 1L, paste0("expected one row for ", nm, " ", mid, " ", ep, ", got ", nrow(r)))
    list(applicable = TRUE, testable = isTRUE(r$testable), n_units = r$n_units, df = r$df,
         n_members = r$n_members, n_measured = r$n_measured, fraction_measured = r$fraction_measured,
         beta = r$beta, se = r$se,
         # ATAC: the association p consumed here is the within-well permutation p.
         p = if (nm == "atac_gse296875") r$within_well_p else r$p_two_sided,
         endpoint = ep, endpoint_note = r$endpoint_note, unit = r$unit,
         cx_p = if ("complexity_adjusted_p" %in% names(r)) r$complexity_adjusted_p else NA_real_,
         n_amb = if ("n_members_from_multigene_groups" %in% names(r)) r$n_members_from_multigene_groups else NA_integer_,
         lowo_min = if ("lowo_min_beta" %in% names(r)) r$lowo_min_beta else NA_real_,
         lowo_max = if ("lowo_max_beta" %in% names(r)) r$lowo_max_beta else NA_real_)
  })
  g <- function(k) unlist(lapply(cells, function(x) { v <- x[[k]]; if (is.null(v)) NA else v }))
  applicable <- as.logical(g("applicable")); testable <- as.logical(g("testable"))
  testable[is.na(testable)] <- FALSE
  usable <- applicable & testable
  beta <- as.numeric(g("beta")); se <- as.numeric(g("se")); df <- as.numeric(g("df"))
  p_assoc <- as.numeric(g("p"))
  eff_or <- beta * direction
  q_assoc <- ml_complete_bh(fifelse(usable, p_assoc, NA_real_), n_family)
  eq_p <- cam_equivalence_p(eff_or, se, df, sesoi)
  eq_q <- ml_complete_bh(fifelse(usable, eq_p, NA_real_), n_family)
  eq_q2 <- ml_complete_bh(fifelse(usable, cam_equivalence_p(eff_or, se, df, sesoi2), NA_real_), n_family)
  ub <- cam_upper_bound(eff_or, se, df)
  cx_p <- as.numeric(g("cx_p"))
  cx_q <- ml_complete_bh(fifelse(usable, cx_p, NA_real_), n_family)
  complexity_robust <- if (nm == "proteome_liver_pxd051911") is.finite(cx_q) & cx_q < 0.05 else NA
  cap <- rep(NA_character_, n_family); cap_reason <- rep(NA_character_, n_family)
  if (nm == "proteome_liver_pxd051911") {
    cap[usable & !complexity_robust] <- "indeterminate"
    cap_reason[usable & !complexity_robust] <- "association not retained at BH q<0.05 under ten intensity-complexity descriptors"
  }
  if (nm == "atac_gse296875" && !atac_ok) {
    cap[usable] <- "indeterminate"
    cap_reason[usable] <- paste0("only ", atac_summary$wells_with_both_classes, " wells carry both fibrosis classes")
  }
  nn <- if (!is.null(nulls)) nulls[assay == nm][match(module_ids, module_id)] else NULL
  comp_p <- if (is.null(nn)) rep(NA_real_, n_family) else nn$competitive_connected_p
  comp_q <- if (is.null(nn)) rep(NA_real_, n_family) else nn$competitive_connected_q
  cm_q <- if (is.null(nn)) rep(NA_real_, n_family) else nn$coherence_matched_q
  cm_gap <- if (is.null(nn)) rep(NA_real_, n_family) else nn$coherence_gap
  a_state <- vapply(seq_len(n_family), function(i) cam_association_state(
    applicable[i], testable[i], q_assoc[i], eff_or[i], eq_q[i], cap[i]), character(1))
  s_state <- vapply(seq_len(n_family), function(i) cam_specificity_state(usable[i], comp_q[i]), character(1))
  rows[[length(rows) + 1L]] <- data.table(
    module_id = module_ids, assay = nm, evidence_axis = "association", question = spec$question,
    unit = g("unit"), n_units = as.integer(g("n_units")),
    endpoint = g("endpoint"), endpoint_family = unname(family_of), endpoint_note = g("endpoint_note"),
    applicable = applicable, endpoint_is_proxy = nm %in% proxy_assays, testable = testable,
    n_members = as.integer(g("n_members")), n_measured = as.integer(g("n_measured")),
    fraction_measured = as.numeric(g("fraction_measured")),
    n_members_from_multigene_groups = as.integer(g("n_amb")),
    effect_raw = beta, effect_oriented = eff_or, effect_se = se, df = df,
    association_p = p_assoc, association_q = q_assoc,
    equivalence_p = eq_p, equivalence_q = eq_q, equivalence_q_sesoi_0p3 = eq_q2,
    upper_bound_95 = ub,
    competitive_p = comp_p, competitive_q = comp_q, coherence_matched_q = cm_q, coherence_gap = cm_gap,
    complexity_adjusted_q = cx_q, complexity_robust = complexity_robust, cap_reason = cap_reason,
    lowo_min_beta = as.numeric(g("lowo_min")), lowo_max_beta = as.numeric(g("lowo_max")),
    association_state = a_state, specificity_state = s_state,
    specificity_descriptive_only = spec_descriptive_only,
    coherence_state = NA_character_, display_state = a_state,
    composition_adjustment = "none")
}

# Spatial: a different question (coherence), its own state column, and the Vu
# cap applied to every testable cell rather than only to significant ones.
sp <- read_if("spatial_results.tsv")
if (!is.null(sp)) {
  for (tag in unique(sp$assay)) {
    s <- sp[assay == tag][match(module_ids, module_id)]
    testable <- isTRUE_vec(s$testable)
    q <- ml_complete_bh(fifelse(testable, s$empirical_p, NA_real_), n_family)
    c_state <- vapply(seq_len(n_family), function(i) cam_coherence_state(
      testable[i], q[i], source_dependent = (tag == "visium_vu")), character(1))
    rows[[length(rows) + 1L]] <- data.table(
      module_id = module_ids, assay = tag, evidence_axis = "spatial_coherence",
      question = "spatial coherence vs matched genes",
      unit = s$unit, n_units = as.integer(s$n_units),
      endpoint = "residual_morans_I", endpoint_family = NA_character_,
      endpoint_note = "residual Moran's I per section, collapsed to donor; not a stage endpoint",
      applicable = FALSE, endpoint_is_proxy = FALSE, testable = testable,
      n_members = s$n_members, n_measured = s$n_measured, fraction_measured = s$fraction_measured,
      n_members_from_multigene_groups = NA_integer_,
      effect_raw = s$z_vs_null, effect_oriented = NA_real_, effect_se = NA_real_, df = NA_real_,
      association_p = s$empirical_p, association_q = q,
      equivalence_p = NA_real_, equivalence_q = NA_real_, equivalence_q_sesoi_0p3 = NA_real_,
      upper_bound_95 = NA_real_,
      competitive_p = NA_real_, competitive_q = NA_real_, coherence_matched_q = NA_real_, coherence_gap = NA_real_,
      complexity_adjusted_q = NA_real_, complexity_robust = NA,
      cap_reason = fifelse(tag == "visium_vu", "arrays cannot be joined to donors", NA_character_),
      lowo_min_beta = NA_real_, lowo_max_beta = NA_real_,
      association_state = "not_applicable", specificity_state = "not_assessed",
      specificity_descriptive_only = NA, coherence_state = c_state, display_state = c_state,
      composition_adjustment = "cell_abundance_factors")
  }
}

# Coverage-only assays: no statistic, and the state says why.
for (f in c("cosmx_coverage.tsv", "plasma_coverage.tsv")) {
  cv <- read_if(f); if (is.null(cv)) next
  cv <- cv[match(module_ids, module_id)]
  st <- fifelse(isTRUE_vec(cv$testable), "indeterminate", "untestable")
  rows[[length(rows) + 1L]] <- data.table(
    module_id = module_ids, assay = cv$assay[1], evidence_axis = "coverage_only",
    question = "coverage only", unit = NA_character_, n_units = NA_integer_,
    endpoint = NA_character_, endpoint_family = NA_character_, endpoint_note = "no statistic computed",
    applicable = FALSE, endpoint_is_proxy = FALSE, testable = isTRUE_vec(cv$testable),
    n_members = cv$n_members, n_measured = cv$n_measured, fraction_measured = cv$fraction_measured,
    n_members_from_multigene_groups = NA_integer_,
    effect_raw = NA_real_, effect_oriented = NA_real_, effect_se = NA_real_, df = NA_real_,
    association_p = NA_real_, association_q = NA_real_,
    equivalence_p = NA_real_, equivalence_q = NA_real_, equivalence_q_sesoi_0p3 = NA_real_,
    upper_bound_95 = NA_real_,
    competitive_p = NA_real_, competitive_q = NA_real_, coherence_matched_q = NA_real_, coherence_gap = NA_real_,
    complexity_adjusted_q = NA_real_, complexity_robust = NA, cap_reason = NA_character_,
    lowo_min_beta = NA_real_, lowo_max_beta = NA_real_,
    association_state = st, specificity_state = "not_assessed",
    specificity_descriptive_only = NA, coherence_state = NA_character_, display_state = st,
    composition_adjustment = "none")
}

panel <- rbindlist(rows, fill = TRUE)

# Overlap with the 117 prespecified programs: an annotation, not evidence that
# the family is a new coordinate system (set size and fragmentation alone can
# lower Jaccard).
hs <- fread(cam_input("hotspot_membership", contract))
hs <- hs[!is.na(mapped_symbol) & mapped_symbol != ""]
programs <- split(hs$mapped_symbol, hs$program_uid)
memb <- fread(file.path(cam_out_root(), "modules", "module_membership.tsv"))
mods_list <- split(memb$gene_symbol, memb$module_id)
jac <- t(vapply(mods_list, function(M) {
  v <- vapply(programs, function(P) length(intersect(M, P)) / length(union(M, P)), numeric(1))
  c(max(v), which.max(v))
}, numeric(2)))
annot <- data.table(module_id = names(mods_list), max_jaccard_vs_117_programs = jac[, 1],
                    nearest_prespecified_program = names(programs)[jac[, 2]])
cam_write_tsv(annot, file.path(out, "module_program_overlap.tsv"))

panel <- merge(panel, labels[, .(module_id, direction, endpoint_family_label = endpoint_family,
                                 discovery_label, discovery_supported, discovery_effect)],
               by = "module_id", sort = FALSE)
panel <- merge(panel, registry[, .(module_id, n_genes, min_loco_heldout_score_spearman,
                                   score_stability_flag)], by = "module_id", sort = FALSE)
panel <- merge(panel, annot, by = "module_id", sort = FALSE)
cam_assert_no_prohibited_columns(panel, contract)
cam_assert(!any(panel$display_state == "supported"), "old state vocabulary leaked into the panel")
cam_assert(!any(is.na(panel$fraction_measured)), "a cell lost its coverage value")
cam_write_tsv(panel, file.path(out, "module_by_assay_states.tsv"))

# Paired transfer is its own family and its own column block.
pt <- read_if("paired_transfer_results.tsv")
if (!is.null(pt)) {
  pt <- merge(pt, labels[, .(module_id, discovery_label, discovery_supported, endpoint_family)],
              by = "module_id", sort = FALSE)
  cam_assert_no_prohibited_columns(pt, contract)
  cam_write_tsv(pt, file.path(out, "paired_transfer_panel.tsv"))
}

# Descriptive counts per assay. No family-level test: the modules share
# participants, a graph and their null draws, so a binomial count is not a
# test and is not printed.
assoc <- panel[evidence_axis == "association"]
counts <- assoc[, .(
  n_applicable = sum(applicable), n_testable = sum(applicable & testable),
  n_replicates = sum(association_state == "replicates"),
  n_discordant = sum(association_state == "discordant"),
  n_tested_negative = sum(association_state == "tested_negative"),
  n_indeterminate = sum(association_state == "indeterminate"),
  n_untestable = sum(association_state == "untestable"),
  n_not_applicable = sum(association_state == "not_applicable"),
  n_specificity_exceeds_bh = sum(specificity_state == "exceeds_coexpressed_sets"),
  n_competitive_p05_uncorrected_descriptive = sum(competitive_p < 0.05, na.rm = TRUE),
  expected_at_0p05_if_independent_descriptive = 0.05 * sum(is.finite(competitive_p)),
  n_replicates_and_exceeds = sum(association_state == "replicates" &
                                 specificity_state == "exceeds_coexpressed_sets"),
  n_meeting_coverage_rule_all_modules = sum(testable),
  endpoint_is_proxy = endpoint_is_proxy[1]), by = assay]
cam_assert_no_prohibited_columns(counts, contract)
cam_write_tsv(counts, file.path(out, "association_counts_by_assay.tsv"))
state_tab <- dcast(panel[, .N, by = .(assay, display_state)], assay ~ display_state, value.var = "N", fill = 0)
cam_write_tsv(state_tab, file.path(out, "state_counts_by_assay.tsv"))
cam_say("association counts by assay:"); print(counts)

# The six prespecified example modules for Figure 5K: top-3 |discovery effect|
# among discovery-supported modules in each family. Rule fixed 2026-09-12.
ex <- labels[discovery_supported == TRUE][order(endpoint_family, -discovery_effect)]
ex <- ex[, head(.SD, 3L), by = endpoint_family][, .(endpoint_family, module_id, discovery_effect, direction)]
cam_write_tsv(ex, file.path(out, "prespecified_example_modules.tsv"))

cam_write_json(list(
  state_self_test = TRUE, direction_applied_once_in = "15_integrate_states.R",
  n_modules = n_family, n_assay_columns = length(unique(panel$assay)), n_cells = nrow(panel),
  sesoi = sesoi, sesoi_sensitivity = sesoi2,
  competitive_calibrated = cal$calibrated, specificity_descriptive_only = spec_descriptive_only,
  atac_wells_with_both_classes = atac_summary$wells_with_both_classes, atac_column_capped = !atac_ok,
  display_states = as.list(table(panel$display_state)),
  association_counts = counts,
  n_tested_negative_sesoi_0p3 = assoc[applicable & testable & association_q >= 0.05 &
                                        equivalence_q_sesoi_0p3 < 0.05, .N],
  modules_replicating_in_any_assay = length(unique(assoc[association_state == "replicates", module_id])),
  modules_replicating_in_any_non_proxy_assay = length(unique(assoc[association_state == "replicates" & !endpoint_is_proxy, module_id])),
  proteome_modules_with_multigene_member = assoc[assay == "proteome_liver_pxd051911" & n_members_from_multigene_groups > 0, .N],
  proteome_replicates_with_multigene_member = assoc[assay == "proteome_liver_pxd051911" & association_state == "replicates" & n_members_from_multigene_groups > 0, .N],
  modules_discordant_in_any_assay = length(unique(assoc[association_state == "discordant", module_id])),
  discovery_supported_modules = sum(labels$discovery_supported),
  score_unstable_modules = sum(registry$score_stability_flag == "score_unstable"),
  max_jaccard_vs_117_programs = max(annot$max_jaccard_vs_117_programs),
  median_jaccard_vs_117_programs = stats::median(annot$max_jaccard_vs_117_programs),
  prespecified_example_modules = ex$module_id
), file.path(out, "panel_summary.json"))
writeLines("panel integrated", file.path(out, "READY"))
cam_say("15 complete")
