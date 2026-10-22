# H49 dry run: build synthetic per-compound limma results and a synthetic L1000
# PHH table with planted effects, to check h49_score.R end to end.
# Usage: Rscript --vanilla dryrun_synthetic.R <gpl570_annot_gz> <zhao_supp5_xlsx> <landmark_order_tsv> <out_dir>
# Planted: amiodarone High_24h early-up t + 5 (concordant); valproic acid High_24h
# early-down t + 3 (anti-concordant); tetracycline High_24h lipogenic genes t + 6;
# omeprazole CYP1A1 and rifampicin CYP3A4 t = 15; L1000 tamoxifen early-up z + 2.
args <- commandArgs(trailingOnly = TRUE)
annot_gz <- args[1]; s5_xlsx <- args[2]; lm_tsv <- args[3]; out <- args[4]
if (dir.exists(out)) stop("refusing to overwrite ", out)
dir.create(file.path(out, "per_compound"), recursive = TRUE)
set.seed(7)
ann <- read.delim(gzfile(annot_gz), skip = 27, colClasses = "character", quote = "", comment.char = "")
ann <- ann[!grepl("^!", ann$ID), c("ID", "Gene.ID")]
probes <- ann$ID
s5 <- as.data.frame(readxl::read_excel(s5_xlsx)); s5$entrez <- as.character(s5[[4]])
q <- as.numeric(s5$qval_Histo_Steatosis); cf <- as.numeric(s5$coef_Histo_Steatosis)
up <- s5$entrez[q < 0.05 & cf > 0]; dn <- s5$entrez[q < 0.05 & cf < 0]
pr_of <- function(e) which(ann$Gene.ID %in% e)
pc_genes <- c("6720", "7069", "47", "2194", "6319", "729359", "123")
cmps <- c("amiodarone", "tamoxifen", "perhexiline", "tetracycline", "nefazodone", "carbamazepine",
          "enalapril", "isoniazid", "ketoconazole", "valproic acid", "acetaminophen", "ibuprofen",
          "naproxen", "aspirin", "dexamethasone", "omeprazole", "caffeine", "chloramphenicol",
          "chlorpheniramine", "hydroxyzine", "theophylline", "rifampicin", "alpidem", "ethionine")
cons <- c("Low_8h", "Middle_8h", "High_8h", "Low_24h", "Middle_24h", "High_24h")
ctrl_mean <- setNames(rnorm(length(probes), 7, 2), probes)
for (cn in cmps) {
  key <- gsub("[^a-z0-9]", "", tolower(cn))
  t <- matrix(rt(length(probes) * length(cons), 14), ncol = length(cons), dimnames = list(probes, cons))
  if (key == "amiodarone") t[pr_of(up), "High_24h"] <- t[pr_of(up), "High_24h"] + 5
  if (key == "valproicacid") t[pr_of(dn), "High_24h"] <- t[pr_of(dn), "High_24h"] + 3
  if (key == "tetracycline") t[pr_of(pc_genes), "High_24h"] <- t[pr_of(pc_genes), "High_24h"] + 6
  if (key == "omeprazole") t[pr_of("1543"), "High_24h"] <- 15
  if (key == "rifampicin") t[pr_of("1576"), "High_24h"] <- 15
  df <- setNames(rep(14, length(probes)), probes)
  vt <- setNames(rt(length(probes), 14), probes)
  res <- list(compound = cn, key = key, exp_id = "0", n_arrays = 14, groups = NULL, doses = NULL,
              t = t, logFC = t * 0.2, df = df, ctrl_mean = ctrl_mean + rnorm(length(probes), 0, 0.1), n_ctrl = 6,
              veh = list(t = vt, logFC = vt * 0.2, df = df, barcodes = c("a", "b")))
  saveRDS(res, file.path(out, "per_compound", paste0(key, ".rds")))
}
lm <- read.delim(lm_tsv, colClasses = "character")$pr_gene_id
names_l <- c("amiodarone", "tamoxifen", "carbamazepine", "enalapril", "atorvastatin", "raloxifene",
             "prednisolone", "propofol", "omeprazole", "loperamide", "caffeine", "lidocaine", "naloxone",
             "hydrocortisone", "digoxin")
zl <- matrix(rnorm(length(names_l) * length(lm)), nrow = length(names_l), dimnames = list(NULL, paste0("g", lm)))
zl[names_l == "tamoxifen", paste0("g", lm[lm %in% up])] <- zl[names_l == "tamoxifen", paste0("g", lm[lm %in% up])] + 2
meta <- data.frame(sig_id = paste0("SYN_", seq_along(names_l)), pert_id = "BRD", pert_iname = names_l, pert_type = "trt_cp",
                   cell_id = "PHH", pert_dose = "10.0", pert_dose_unit = "µM", pert_time = "24", pert_time_unit = "h",
                   distil_cc_q75 = "0.3", distil_nsample = "3", check.names = FALSE)
write.table(cbind(meta, zl), gzfile(file.path(out, "lincs_phh_synthetic.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
cat("synthetic inputs written to", out, "\n")
