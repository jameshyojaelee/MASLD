#!/usr/bin/env Rscript
# 13: CosMx and the other coverage-only assays.
#
# A 968-gene targeted panel cannot test a whole-module hypothesis, and the
# Resource's own accounting says so: it adequately observes 12 of the 117
# frozen programs. Reporting it as a blank column would hide the reason. This
# step produces the coverage number and nothing else, so the panel can show why
# the cell is empty rather than leaving the reader to guess.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))

contract <- cam_contract()
out <- cam_dir("assays")
membership <- fread(file.path(cam_out_root(), "modules", "module_membership.tsv"))
modules <- split(membership$gene_symbol, membership$module_id)

cosmx <- readLines(file.path(cam_out_root(), "universe", "measured_cosmx_govaere.txt"))
rows <- rbindlist(lapply(names(modules), function(mid) {
  tb <- cam_testability(modules[[mid]], cosmx, contract)
  data.table(assay = "cosmx_govaere", module_id = mid, role = "coverage_only",
             n_members = tb$n_members, n_measured = tb$n_measured,
             fraction_measured = tb$fraction_measured, testable = tb$testable,
             reason = "targeted 968-gene panel; whole-module scoring not supported")
}))
cam_assert_no_prohibited_columns(rows, contract)
cam_write_tsv(rows, file.path(out, "cosmx_coverage.tsv"))
cam_write_json(list(
  panel_size = length(cosmx),
  modules_meeting_coverage_rule = sum(rows$testable),
  median_fraction_measured = stats::median(rows$fraction_measured),
  max_fraction_measured = max(rows$fraction_measured)
), file.path(out, "cosmx_summary.json"))
writeLines("cosmx coverage done", file.path(out, "READY_cosmx"))
cam_say("13 complete: ", sum(rows$testable), " of ", nrow(rows),
        " modules meet the coverage rule on the CosMx panel")
