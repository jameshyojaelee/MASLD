library(data.table)
setDTthreads(1)
base <- Sys.getenv("MASLD_PROJECT_ROOT")
out <- Sys.getenv("DIAGNOSIS_REFIT_ROOT")
meta_path <- Sys.getenv("QUALIFIED_DONOR_METADATA")
old <- file.path(base, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot")
dir.create(file.path(out, "hotspot"), recursive=TRUE)
EXPECTED_PROGRAMS <- 117L
direction_of <- function(x) fifelse(is.na(x), "untestable", fifelse(x>0,"positive",fifelse(x<0,"negative","zero")))
# Evaluate only the unchanged model functions, never the producer's data/output code.
exprs <- parse(file.path(base,"Analysis/SingleCell/scripts/hotspot_modules/512_hotspot_v2_donor_refit.R"))
wanted <- c("fit_fixed_ols","fit_program_family","fit_context_effects","classify_source_stability")
for (e in exprs) if (is.call(e) && identical(e[[1]],as.name("<-")) && as.character(e[[2]]) %in% wanted) eval(e)
meta <- fread(meta_path)
if ("sample" %in% names(meta)) setnames(meta,"sample","donor")
stopifnot(!anyDuplicated(meta$donor))
meta[, stage_ordinal := match(disease_stage_coarse,c("Healthy","Steatosis","Steatohepatitis","Cirrhosis"))-1L]
eligible <- meta[tolower(as.character(exclude_stage_analysis)) %in% c("false","f","0") & disease_stage_coarse %in% c("Healthy","Steatosis","Steatohepatitis")]
scores <- function(name) {
 x <- fread(file.path(old,name)); x[, c("dataset","disease_stage_coarse","stage_ordinal","exclude_stage_analysis") := NULL]
 merge(x,eligible[,.(donor,dataset,disease_stage_coarse,stage_ordinal,exclude_stage_analysis)],by="donor")
}
p <- scores("donor_program_scores_primary.tsv"); eq <- scores("donor_program_scores_equal_run.tsv")
reg <- fread(file.path(old,"program_registry_v2.tsv")); programs <- unique(reg[,.(cell_type,module)])
stopifnot(nrow(programs)==117L,uniqueN(p$donor)==62L,!any(p$donor %in% c("SRR14602057","SRR14602058")))
a <- fit_program_family(p,programs,"stage_ordinal",TRUE,TRUE)
b <- fit_program_family(eq,programs,"stage_ordinal",TRUE,TRUE)
context <- fit_context_effects(p,programs,sort(unique(p$dataset)))
stability <- classify_source_stability(a,context)
for (i in seq_len(nrow(reg))) {
 j <- which(a$cell_type==reg$cell_type[i] & a$module==reg$module[i]); k <- which(b$cell_type==reg$cell_type[i] & b$module==reg$module[i])
 for (nm in names(a)) if (paste0("primary_",nm) %in% names(reg)) set(reg,i,paste0("primary_",nm),a[[nm]][j])
 for (nm in names(b)) if (paste0("equal_run_",nm) %in% names(reg)) set(reg,i,paste0("equal_run_",nm),b[[nm]][k])
 for (nm in setdiff(names(stability),c("cell_type","module"))) if (nm %in% names(reg)) set(reg,i,nm,stability[[nm]][j])
}
reg[, primary_selected := primary_estimable & is.finite(primary_qvalue) & primary_qvalue<.05 & is.finite(stability_median) & stability_median>=.5]
reg[, hc3_supported := primary_estimable & is.finite(primary_hc3_qvalue) & primary_hc3_qvalue<.05]
reg[, direction_agree := primary_estimable & equal_run_estimable & is.finite(primary_beta) & is.finite(equal_run_beta) & primary_beta!=0 & equal_run_beta!=0 & sign(primary_beta)==sign(equal_run_beta)]
reg[, robust_display := primary_selected & direction_agree & hc3_supported]
reg[, scoring_direction_agree := direction_agree]
reg[, selected_unstable := primary_estimable & is.finite(primary_qvalue) & primary_qvalue<.05 & (!is.finite(stability_median) | stability_median<.5)]
reg[, selected_hc3_fragile := primary_selected & scoring_direction_agree & !hc3_supported]
reg[, external_test_eligible := robust_display & cell_type=="hepatocytes"]
# Historical native flag meant ordinary nonsignificance, not adequate negative evidence.
reg[, source_native_tested_negative := tested_negative]
reg[, tested_nonsignificant := primary_estimable & is.finite(primary_qvalue) & primary_qvalue>=.05]
reg[, tested_negative := FALSE]
reg[, registry_state := fcase(!primary_estimable,"untestable",selected_unstable,"selected_unstable",primary_selected & !scoring_direction_agree,"selected_direction_discordant",selected_hc3_fragile,"selected_hc3_fragile",robust_display,"robust_display",tested_nonsignificant,"tested_nonsignificant",default="testable_unclassified")]
for (nm in c("program_membership_v2.tsv","READY")) file.copy(file.path(old,nm),file.path(out,"hotspot",nm))
for (nm in c("program_registry_v2","donor_program_scores_primary","donor_program_scores_equal_run","primary_effects","equal_effects","context_effects","source_stability")) {
 x <- switch(nm,program_registry_v2=reg,donor_program_scores_primary=p,donor_program_scores_equal_run=eq,primary_effects=a,equal_effects=b,context_effects=context,source_stability=stability)
 fwrite(x,file.path(out,"hotspot",paste0(nm,".tsv")),sep="\t")
}
fwrite(meta,file.path(out,"qualified_donor_metadata.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
