# H49 step C: score early-MASLD gene sets in steatogenic-drug PHH contrasts.
#
# Usage:
#   Rscript --vanilla h49_score.R <per_compound_dir> <gpl570_annot_gz> <zhao_supp5_xlsx>
#       <zhao_supp9_xlsx> <lincs_phh_tsv_gz> <dilirank_xlsx> <out_dir>
# Env: H49_N_NULL (default 10000; dry runs use fewer).
# Design and thresholds: Analysis/Hypotheses/H49-20261009T015504Z/PRESPEC_H49.md

suppressPackageStartupMessages({
  library(Matrix)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) stop("7 arguments required")
pc_dir <- args[1]; annot_gz <- args[2]; s5_xlsx <- args[3]; s9_xlsx <- args[4]
lincs_tsv <- args[5]; dili_xlsx <- args[6]; out_dir <- args[7]
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "h49_lib.R"))
if (dir.exists(out_dir) && length(list.files(out_dir))) stop("refusing to write into non-empty ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
N_NULL <- as.integer(Sys.getenv("H49_N_NULL", "10000"))
SEED <- 20260923L
wt <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE)
log_lines <- character(0)
say <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }

# ---------------------------------------------------------------- drug lists (prespec section 4)
STEAT <- c("amiodarone", "tamoxifen", "perhexiline", "tetracycline", "nefazodone",
           "carbamazepine", "enalapril", "isoniazid", "ketoconazole", "valproicacid")
STEAT_SEC <- c("acetaminophen", "ibuprofen", "naproxen", "aspirin", "dexamethasone")
NEG <- c("omeprazole", "caffeine", "chloramphenicol", "chlorpheniramine", "hydroxyzine", "theophylline")
PC_DRUGS <- c("amiodarone", "tetracycline")
PC_GENES <- c(SREBF1 = "6720", THRSP = "7069", ACLY = "47", FASN = "2194", SCD = "6319",
              PLIN4 = "729359", PLIN2 = "123")
ASSAY <- data.frame(drug = c("omeprazole", "rifampicin"), gene = c("CYP1A1", "CYP3A4"),
                    entrez = c("1543", "1576"), stringsAsFactors = FALSE)
PRIMARY <- "High_24h"
L_STEAT <- c("amiodarone", "tamoxifen", "carbamazepine", "enalapril", "atorvastatin",
             "raloxifene", "prednisolone", "propofol")
L_NEG_CLUSTER0 <- c("omeprazole", "loperamide")
L_NEG_EXCLUDE <- c("hydrocortisone", "cortisone-acetate")

# ---------------------------------------------------------------- array annotation
ann <- read.delim(gzfile(annot_gz), skip = 27, colClasses = "character", quote = "", comment.char = "")
ann <- ann[!grepl("^!", ann$ID), c("ID", "Gene.symbol", "Gene.ID")]
ann <- ann[ann$Gene.ID != "" & !grepl("///", ann$Gene.ID), ]
say("single-gene probe sets: ", nrow(ann), "; genes: ", length(unique(ann$Gene.ID)))

# ---------------------------------------------------------------- per-compound limma results
files <- sort(list.files(pc_dir, pattern = "\\.rds$", full.names = TRUE))
if (!length(files)) stop("no per-compound results")
res <- lapply(files, readRDS)
names(res) <- vapply(res, `[[`, "", "key")
say("compounds loaded: ", length(res))
probes <- rownames(res[[1]]$t)
for (r in res) stopifnot(identical(rownames(r$t), probes))
num <- Reduce(`+`, lapply(res, function(r) r$ctrl_mean * r$n_ctrl))
den <- sum(vapply(res, `[[`, 0, "n_ctrl"))
pooled_ctrl <- num / den
ann <- ann[ann$ID %in% probes, ]
ann$ctrl <- pooled_ctrl[ann$ID]
ann <- ann[order(ann$Gene.ID, -ann$ctrl, ann$ID), ]
sel <- ann[!duplicated(ann$Gene.ID), ]
genes <- sel$Gene.ID
sym <- setNames(sel$Gene.symbol, genes)
gene_ctrl <- setNames(sel$ctrl, genes)
wt(sel[, c("Gene.ID", "Gene.symbol", "ID", "ctrl")], "probe_selection.tsv")
say("array genes (one probe set each): ", length(genes))

cond_id <- character(0); cond_cmp <- character(0); cond_con <- character(0)
Zl <- list(); Ll <- list(); Ql <- list()
veh_id <- character(0); Zv <- list()
for (k in names(res)) {
  r <- res[[k]]
  t <- r$t[sel$ID, , drop = FALSE]; df <- r$df[sel$ID]; lf <- r$logFC[sel$ID, , drop = FALSE]
  for (j in seq_len(ncol(t))) {
    cond_id <- c(cond_id, paste0(k, "|", colnames(t)[j])); cond_cmp <- c(cond_cmp, k); cond_con <- c(cond_con, colnames(t)[j])
    Zl[[length(Zl) + 1]] <- t_to_z(t[, j], df)
    Ll[[length(Ll) + 1]] <- lf[, j]
    Ql[[length(Ql) + 1]] <- p.adjust(t_to_p(t[, j], df), "BH")
  }
  if (!is.null(r$veh)) {
    veh_id <- c(veh_id, paste0(k, "|veh24"))
    Zv[[length(Zv) + 1]] <- t_to_z(r$veh$t[sel$ID], r$veh$df[sel$ID])
  }
}
Z <- do.call(cbind, Zl); L <- do.call(cbind, Ll); Q <- do.call(cbind, Ql)
dimnames(Z) <- dimnames(L) <- dimnames(Q) <- list(genes, cond_id)
ZV <- do.call(cbind, Zv); dimnames(ZV) <- list(genes, veh_id)
rm(Zl, Ll, Ql, Zv, res); invisible(gc())
say("conditions: ", ncol(Z), "; vehicle contrasts: ", ncol(ZV))
role <- ifelse(cond_cmp %in% STEAT, "steatogenic", ifelse(cond_cmp %in% STEAT_SEC, "steatogenic_secondary",
         ifelse(cond_cmp %in% NEG, "non_steatogenic_control", "other")))

# ---------------------------------------------------------------- human gene sets (prespec section 2)
s5 <- as.data.frame(readxl::read_excel(s5_xlsx))
names(s5)[1:4] <- c("ensg", "symbol", "desc", "entrez")
s5$entrez <- as.character(s5$entrez)
num5 <- function(x) suppressWarnings(as.numeric(x))
zhao_entrez <- unique(s5$entrez)
U <- genes[genes %in% zhao_entrez]
qS <- num5(s5$qval_Histo_Steatosis); cS <- num5(s5$coef_Histo_Steatosis)
qF <- num5(s5$qval_Histo_Fibrosis); cF <- num5(s5$coef_Histo_Fibrosis)
in_U <- function(e) U[U %in% unique(e)]
sets <- list(
  early_up = in_U(s5$entrez[!is.na(qS) & qS < 0.05 & cS > 0]),
  early_down = in_U(s5$entrez[!is.na(qS) & qS < 0.05 & cS < 0]))
s9 <- as.data.frame(readxl::read_excel(s9_xlsx))
sym2ent <- setNames(s5$entrez, s5$symbol)
sets$S9_up <- in_U(sym2ent[na.omit(s9[["Steatosis-specific upregulated genes"]])])
sets$S9_down <- in_U(sym2ent[na.omit(s9[["Steatosis-specific downregulated genes"]])])
sets$F5_up <- in_U(s5$entrez[!is.na(qF) & qF < 0.05 & cF > 0])
sets$F5_down <- in_U(s5$entrez[!is.na(qF) & qF < 0.05 & cF < 0])
U_expr <- U[gene_ctrl[U] >= median(gene_ctrl[U])]
sets$expr_up <- sets$early_up[sets$early_up %in% U_expr]
sets$expr_down <- sets$early_down[sets$early_down %in% U_expr]
say("U: ", length(U), "; U_expr: ", length(U_expr), "; set sizes: ",
    paste(names(sets), lengths(sets), sep = "=", collapse = ", "))
gs <- do.call(rbind, lapply(names(sets), function(s) data.frame(set = s, entrez = sets[[s]], symbol = sym[sets[[s]]])))
wt(gs, "gene_sets_used.tsv")

# ---------------------------------------------------------------- L1000 PHH (prespec sections 3-4)
lz <- read.delim(gzfile(lincs_tsv), colClasses = "character", check.names = FALSE)
lm_cols <- grep("^g[0-9]+$", names(lz), value = TRUE)
lm_genes <- sub("^g", "", lm_cols)
UL <- lm_genes[lm_genes %in% zhao_entrez]
# landmark sets come from the Zhao list directly (not via the array universe U)
lsets <- list(early_up = UL[UL %in% s5$entrez[!is.na(qS) & qS < 0.05 & cS > 0]],
              early_down = UL[UL %in% s5$entrez[!is.na(qS) & qS < 0.05 & cS < 0]])
say("L1000 UL: ", length(UL), "; up ", length(lsets$early_up), "; down ", length(lsets$early_down),
    "; PC genes on landmarks: ", sum(PC_GENES %in% lm_genes))

# ---------------------------------------------------------------- random draws (fixed order, prespec section 6)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(SEED)
draw_order <- c("early_up", "early_down", "S9_up", "S9_down", "F5_up", "F5_down", "expr_up", "expr_down")
# Null background (addendum 1): the set's universe minus that family's up and down
# genes, so a shift of the opposite-direction set cannot move the null.
bg <- list()
for (s in draw_order) {
  univ <- if (s %in% c("expr_up", "expr_down")) U_expr else U
  fam <- sub("_(up|down)$", "", s)
  bg[[s]] <- univ[!univ %in% c(sets[[paste0(fam, "_up")]], sets[[paste0(fam, "_down")]])]
}
bg_pc <- genes[!genes %in% PC_GENES]
bg_l <- UL[!UL %in% c(lsets$early_up, lsets$early_down)]
M <- list()
for (s in draw_order) M[[s]] <- draw_sets(length(bg[[s]]), length(sets[[s]]), N_NULL)
M_pc <- draw_sets(length(bg_pc), length(PC_GENES), N_NULL)
M_l <- list(early_up = draw_sets(length(bg_l), length(lsets$early_up), N_NULL),
            early_down = draw_sets(length(bg_l), length(lsets$early_down), N_NULL))
say("random draws: ", N_NULL, " per set, seed ", SEED, "; background sizes: ",
    paste(names(bg), lengths(bg), sep = "=", collapse = ", "), ", pc=", length(bg_pc), ", L1000=", length(bg_l))

score_set <- function(s, mat) {
  obs <- set_mean(mat, sets[[s]])
  null <- as.matrix(M[[s]] %*% mat[bg[[s]], , drop = FALSE])
  colnames(null) <- names(obs)
  list(S = obs, null_med = col_median(null), null_sd = apply(null, 2, sd),
       p = perm_p_two(obs, null))
}
frac_set <- function(s, ind) {
  obs <- colMeans(ind[sets[[s]], , drop = FALSE])
  null <- as.matrix(M[[s]] %*% ind[bg[[s]], , drop = FALSE])
  colnames(null) <- names(obs)
  list(F = obs, null_med = col_median(null), p = perm_p_greater(obs, null))
}

# ---------------------------------------------------------------- all conditions, primary sets
sc <- list(up = score_set("early_up", Z), down = score_set("early_down", Z))
conc <- function(dir, s) if (dir == "up") s$S - s$null_med else s$null_med - s$S
allc <- data.frame(condition = cond_id, compound = cond_cmp, contrast = cond_con, role = role,
                   S_up = sc$up$S, null_med_up = sc$up$null_med, p_up = sc$up$p, conc_up = conc("up", sc$up),
                   S_down = sc$down$S, null_med_down = sc$down$null_med, p_down = sc$down$p, conc_down = conc("down", sc$down))
wt(allc, "all_conditions_primary_sets.tsv")

# ---------------------------------------------------------------- primary family
prim <- paste0(STEAT, "|", PRIMARY)
missing_prim <- setdiff(prim, cond_id)
if (length(prim) - length(missing_prim) < 8) stop("stop rule: fewer than 8 primary drugs processed: missing ", paste(missing_prim, collapse = ","))
prim <- intersect(prim, cond_id)
negc <- intersect(paste0(NEG, "|", PRIMARY), cond_id)
IU <- (L > 0 & Q < 0.05) * 1; ID <- (L < 0 & Q < 0.05) * 1
fr <- list(up = frac_set("early_up", IU[, c(prim, negc), drop = FALSE]),
           down = frac_set("early_down", ID[, c(prim, negc), drop = FALSE]))

# vehicle versus vehicle
vs <- list(up = score_set("early_up", ZV), down = score_set("early_down", ZV))
veh <- data.frame(contrast = veh_id, S_up = vs$up$S, p_up = vs$up$p, conc_up = conc("up", vs$up),
                  S_down = vs$down$S, p_down = vs$down$p, conc_down = conc("down", vs$down))
wt(veh, "vehicle_vs_vehicle.tsv")
veh_sum <- data.frame(direction = c("up", "down"), n_contrasts = nrow(veh),
                      frac_p_lt_0.05 = c(mean(veh$p_up < 0.05), mean(veh$p_down < 0.05)),
                      q95_abs_dev = c(quantile(abs(vs$up$S - vs$up$null_med), 0.95), quantile(abs(vs$down$S - vs$down$null_med), 0.95)))
veh_sum$pass <- veh_sum$frac_p_lt_0.05 <= 0.10
wt(veh_sum, "vehicle_summary.tsv")
say("vehicle FPR up ", round(veh_sum$frac_p_lt_0.05[1], 3), ", down ", round(veh_sum$frac_p_lt_0.05[2], 3))

mk_rows <- function(conds, dir) {
  s <- sc[[dir]]; f <- fr[[dir]]
  data.frame(drug = sub("\\|.*$", "", conds), condition = conds, direction = dir,
             n_genes = length(sets[[paste0("early_", dir)]]),
             S = s$S[conds], null_median = s$null_med[conds], null_sd = s$null_sd[conds],
             conc_S = conc(dir, s)[conds], p_two = s$p[conds],
             F = f$F[conds], F_null_median = f$null_med[conds], p_F = f$p[conds], row.names = NULL)
}
pt_ <- rbind(mk_rows(prim, "up"), mk_rows(prim, "down"))
stopifnot(!anyNA(pt_[, c("S", "null_median", "p_two", "F", "p_F")]))
pt_$concordant <- pt_$conc_S > 0
pt_$q_BH <- p.adjust(pt_$p_two, "BH")
pt_$q_F_BH <- p.adjust(pt_$p_F, "BH")
nt <- rbind(mk_rows(negc, "up"), mk_rows(negc, "down"))
wt(nt, "negative_control_drugs.tsv")
neg_max <- tapply(nt$conc_S, nt$direction, max)
pt_$neg_ctrl_max_conc_S <- neg_max[pt_$direction]
pt_$exceeds_all_neg <- pt_$conc_S > pt_$neg_ctrl_max_conc_S
pt_$vehicle_rule_needed <- !veh_sum$pass[match(pt_$direction, veh_sum$direction)]
pt_$vehicle_q95 <- veh_sum$q95_abs_dev[match(pt_$direction, veh_sum$direction)]
pt_$passes_vehicle_rule <- ifelse(pt_$vehicle_rule_needed, abs(pt_$S - pt_$null_median) > pt_$vehicle_q95, TRUE)
pt_$significant_concordant <- pt_$concordant & pt_$q_BH < 0.05 & pt_$passes_vehicle_rule
pt_$significant_anticoncordant <- !pt_$concordant & pt_$q_BH < 0.05 & pt_$passes_vehicle_rule
wt(pt_, "primary_tests.tsv")

# class test: steatogenic vs non-steatogenic, exact permutation
cls <- do.call(rbind, lapply(c("up", "down"), function(dir) {
  x <- c(pt_$conc_S[pt_$direction == dir], nt$conc_S[nt$direction == dir])
  a <- c(rep(TRUE, sum(pt_$direction == dir)), rep(FALSE, sum(nt$direction == dir)))
  e <- exact_perm_diff(x, a)
  data.frame(direction = dir, n_steatogenic = sum(a), n_controls = sum(!a),
             mean_conc_S_steatogenic = mean(x[a]), mean_conc_S_controls = mean(x[!a]),
             diff = e$obs, p_exact = e$p, n_relabelings = e$n_perm)
}))
cls$q_BH <- p.adjust(cls$p_exact, "BH")
wt(cls, "class_test.tsv")

# ---------------------------------------------------------------- positive control and assay check
pc_conds <- intersect(paste0(PC_DRUGS, "|", PRIMARY), cond_id)
pc_rows <- match(PC_GENES, genes)
pc_obs <- colMeans(Z[pc_rows, pc_conds, drop = FALSE])
pc_null <- as.matrix(M_pc %*% Z[bg_pc, pc_conds, drop = FALSE])
pc <- data.frame(condition = pc_conds, n_genes = length(pc_rows), S = pc_obs,
                 null_median = col_median(pc_null), p_greater = perm_p_greater(pc_obs, pc_null))
pc$q_BH <- p.adjust(pc$p_greater, "BH")
pc_pass <- any(pc$q_BH < 0.05)
wt(pc, "positive_control.tsv")
pcg <- do.call(rbind, lapply(pc_conds, function(cn) data.frame(condition = cn, gene = names(PC_GENES), entrez = PC_GENES,
          logFC = L[pc_rows, cn], z = Z[pc_rows, cn], q_within = Q[pc_rows, cn])))
wt(pcg, "positive_control_genes.tsv")
ac <- do.call(rbind, lapply(seq_len(nrow(ASSAY)), function(i) {
  cn <- paste0(cmp_key(ASSAY$drug[i]), "|", PRIMARY); g <- ASSAY$entrez[i]
  if (!(cn %in% cond_id) || !(g %in% genes)) return(data.frame(condition = cn, gene = ASSAY$gene[i], logFC = NA, q_within = NA, pass = FALSE))
  data.frame(condition = cn, gene = ASSAY$gene[i], logFC = L[g, cn], q_within = Q[g, cn], pass = L[g, cn] > 0 & Q[g, cn] < 0.05)
}))
assay_pass <- all(ac$pass)
wt(ac, "assay_check.tsv")
say("positive control pass: ", pc_pass, "; assay check pass: ", assay_pass)

# ---------------------------------------------------------------- secondary families (descriptive)
sec_family <- function(su, sd_, label) {
  out <- do.call(rbind, lapply(list(c("up", su), c("down", sd_)), function(v) {
    s <- score_set(v[2], Z[, prim, drop = FALSE])
    data.frame(family = label, drug = sub("\\|.*$", "", prim), direction = v[1], n_genes = length(sets[[v[2]]]),
               S = s$S, null_median = s$null_med, conc_S = if (v[1] == "up") s$S - s$null_med else s$null_med - s$S,
               p_two = s$p, row.names = NULL)
  }))
  out$q_BH <- p.adjust(out$p_two, "BH")
  out
}
sec <- rbind(sec_family("S9_up", "S9_down", "steatosis_specific_S9"),
             sec_family("F5_up", "F5_down", "fibrosis_F5"),
             sec_family("expr_up", "expr_down", "expressed_universe"))
wt(sec, "secondary_families.tsv")

# ---------------------------------------------------------------- descriptive: background ranks at High 24 h
hi <- allc[allc$contrast == PRIMARY, ]
hi$rank_conc_up <- rank(-hi$conc_up); hi$rank_conc_down <- rank(-hi$conc_down)
hi$n_compounds <- nrow(hi)
wt(hi[order(hi$role, hi$compound), ], "background_ranks_high24h.tsv")

# ---------------------------------------------------------------- L1000 PHH (secondary)
lz$key <- tolower(lz$pert_iname)
lz$cc <- num5(lz$distil_cc_q75)
pick <- function(nm) {
  s <- lz[lz$key == nm & abs(num5(lz$pert_dose) - 10) < 1e-6 & num5(lz$pert_time) == 24 & lz$pert_time_unit == "h", ]
  if (!nrow(s)) return(NULL)
  s[order(-s$cc, s$sig_id), ][1, ]
}
dili <- as.data.frame(readxl::read_excel(dili_xlsx, skip = 1))
salts <- "[[:space:]]+(hydrochloride|hcl|sulfate|sodium|potassium|maleate|citrate|mesylate|tartrate|bitartrate|acetate|phosphate|fumarate|succinate|besylate|bromide|hydrobromide|dihydrochloride|monohydrate|calcium|magnesium|nitrate|lactate|gluconate|decanoate|palmitate|glycinate|disodium|tromethamine|meglumine|hyclate|estolate)$"
normn <- function(s) { s <- tolower(trimws(gsub("-", " ", s))); for (i in 1:3) s <- sub(salts, "", s); s }
vno <- unique(normn(dili$CompoundName[tolower(dili[["vDILI-Concern"]]) == "vno-dili-concern"]))
phh_names <- unique(lz$key)
l_neg_names <- sort(unique(c(L_NEG_CLUSTER0, phh_names[normn(phh_names) %in% vno])))
l_neg_names <- setdiff(l_neg_names, c(L_NEG_EXCLUDE, L_STEAT))
l_st <- do.call(rbind, lapply(L_STEAT, pick)); l_ng <- do.call(rbind, lapply(l_neg_names, pick))
say("L1000 steatogenic signatures: ", nrow(l_st), "; control signatures: ", nrow(l_ng))
lmat <- function(df) {
  m <- t(matrix(num5(as.matrix(df[, lm_cols])), nrow = nrow(df)))
  dimnames(m) <- list(lm_genes, df$key)
  m
}
LZ <- cbind(lmat(l_st), lmat(l_ng))
lsc <- lapply(c(up = "early_up", down = "early_down"), function(s) {
  obs <- colMeans(LZ[lsets[[s]], , drop = FALSE])
  null <- as.matrix(M_l[[s]] %*% LZ[bg_l, , drop = FALSE])
  list(S = obs, null_med = col_median(null), p = perm_p_two(obs, null))
})
lt <- do.call(rbind, lapply(c("up", "down"), function(dir) {
  s <- lsc[[dir]]
  data.frame(drug = colnames(LZ), role = ifelse(colnames(LZ) %in% L_STEAT, "steatogenic", "non_steatogenic_control"),
             sig_id = c(l_st$sig_id, l_ng$sig_id), direction = dir, n_genes = length(lsets[[paste0("early_", dir)]]),
             S = s$S, null_median = s$null_med, conc_S = if (dir == "up") s$S - s$null_med else s$null_med - s$S,
             p_two = s$p, row.names = NULL)
}))
lt$q_BH <- NA_real_
st_rows <- lt$role == "steatogenic"
lt$q_BH[st_rows] <- p.adjust(lt$p_two[st_rows], "BH")
lt$frac_controls_ge <- vapply(seq_len(nrow(lt)), function(i) {
  cc <- lt$conc_S[lt$role != "steatogenic" & lt$direction == lt$direction[i]]; mean(cc >= lt$conc_S[i]) }, 0)
wt(lt, "lincs_phh_tests.tsv")

# ---------------------------------------------------------------- decision (prespec section 8)
any_conc <- any(pt_$significant_concordant)
spec <- pt_$significant_concordant & pt_$exceeds_all_neg
answer <- if (any_conc && any(spec)) {
  "no: early genes are induced by at least one steatogenic drug, more than by every non-steatogenic control (shared early response)"
} else if (any_conc) {
  "no for the literal claim: early genes are induced by at least one steatogenic drug, but not more than by non-steatogenic controls (not specific)"
} else if (pc_pass && assay_pass) {
  "yes: no steatogenic drug induces early-MASLD genes in the concordant direction at BH q < 0.05; controls passed"
} else {
  "untestable: no concordant induction, but the positive control or assay check failed"
}
say("answer: ", answer)
summ <- list(
  answer = answer, n_null = N_NULL, seed = SEED, primary_condition = PRIMARY,
  universe_U = length(U), universe_U_expr = length(U_expr), set_sizes = as.list(lengths(sets)),
  null_background_sizes = c(as.list(lengths(bg)), list(positive_control = length(bg_pc), L1000 = length(bg_l))),
  n_conditions = ncol(Z), n_vehicle_contrasts = ncol(ZV),
  primary_family_size = nrow(pt_), primary_drugs = sub("\\|.*$", "", prim), missing_primary = missing_prim,
  n_sig_concordant = sum(pt_$significant_concordant), n_sig_anticoncordant = sum(pt_$significant_anticoncordant),
  n_sig_concordant_specific = sum(spec), positive_control_pass = pc_pass, assay_check_pass = assay_pass,
  vehicle_pass = setNames(as.list(veh_sum$pass), veh_sum$direction),
  lincs_family_size = sum(st_rows), lincs_n_controls = nrow(l_ng),
  lincs_pc_genes_on_landmarks = sum(PC_GENES %in% lm_genes))
writeLines(toJSON(summ, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(out_dir, "summary.json"))
writeLines(log_lines, file.path(out_dir, "score_log.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_score.txt"))
