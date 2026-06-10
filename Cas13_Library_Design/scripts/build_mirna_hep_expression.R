suppressMessages(library(data.table))

## ============ key normalization helpers (all validated) ============
# strip a terminal hairpin-copy index only when a digit (the family number) remains in the stem
.strip_copy <- function(x) {
  stem <- sub("-?[0-9]+$", "", x)
  ifelse((stem != x) & grepl("[0-9]", stem), stem, x)
}
# TCGA hairpin id  ->  family key  (e.g. hsa-mir-194-1 -> mir194 ; hsa-mir-122 -> mir122)
fk_from_hairpin <- function(x) {
  x <- tolower(x); x <- sub("^hsa-", "", x)
  stem <- sub("-[0-9]+$", "", x)
  drop <- (stem != x) & grepl("[0-9]", stem)
  out <- ifelse(drop, stem, x)
  gsub("-", "", out)
}
# conserved mirbase_family (miR-122, let-7a) -> family key
fk_from_family <- function(f) {
  f <- tolower(f); f <- sub("^mir-", "mir", f); f <- sub("^let-", "let", f)
  gsub("-", "", f)
}
# human gene symbol (MIRLET7A1, MIR194-1, MIR92A1) -> family key. Strip dash-copy BEFORE removing dashes.
fk_from_symbol <- function(s) {
  s <- tolower(s)
  s <- sub("^mirlet", "let", s)        # MIRLET7A1 -> let7a1
  s <- sub("^mir", "mir", s)
  s <- sub("-[0-9]+$", "", s)          # MIR194-1 -> mir194  (dash-delimited copy)
  s <- gsub("-", "", s)
  .strip_copy(s)                       # let7a1->let7a, mir103a1->mir103a, mir92a1->mir92a, mir122 stays
}

## ============ 1. Candidate set: 325 conserved mouse miRNA genes ============
master <- fread("data/external/orthologs/master_ortholog_table.tsv.gz")
cand <- master[mouse_biotype == "miRNA" & (tier_H_mirbase == 1 | tier_H_mirgenedb == 1),
               .(mouse_ensembl, mouse_symbol, human_symbol_master = human_symbol,
                 mirbase_family_m = mirbase_family, mirbase_arm_m = mirbase_arm)]
cand <- unique(cand, by = "mouse_ensembl")

## ============ 2. L2 / L2b layers ============
L2 <- fread("data/external/orthologs/layers/L2_mirbase.tsv")
L2 <- unique(L2[tier_H_mirbase == 1, .(mouse_ensembl, human_symbol, mirbase_family,
                                       mirbase_arm, mirbase_human_mature_name)], by = "mouse_ensembl")
L2b <- fread("data/external/orthologs/layers/L2b_mirgenedb.tsv")
L2b <- unique(L2b[tier_H_mirgenedb == 1, .(mouse_ensembl, human_symbol_b = human_symbol,
                                           mirgenedb_family, mirgenedb_arm)], by = "mouse_ensembl")

ann <- merge(cand, L2,  by = "mouse_ensembl", all.x = TRUE)
ann <- merge(ann, L2b, by = "mouse_ensembl", all.x = TRUE)

ann[, gene_symbol_human := fifelse(!is.na(human_symbol) & human_symbol != "", human_symbol,
                            fifelse(!is.na(human_symbol_b) & human_symbol_b != "", human_symbol_b, human_symbol_master))]
ann[, mirbase_family := fifelse(!is.na(mirbase_family) & mirbase_family != "", mirbase_family,
                         fifelse(!is.na(mirbase_family_m) & mirbase_family_m != "", mirbase_family_m, mirgenedb_family))]
ann[, mirbase_arm := fifelse(!is.na(mirbase_arm) & mirbase_arm != "", mirbase_arm,
                      fifelse(!is.na(mirbase_arm_m) & mirbase_arm_m != "", mirbase_arm_m, mirgenedb_arm))]

## ============ 3. Liver reference (TCGA-LIHC normal) aggregated to family ============
liver <- fread("data/external/mirna_liver_atlas/tcga_lihc_normal/tcga_lihc_normal_hairpin_rpm.csv")
liver[, fkey := fk_from_hairpin(hairpin_id)]
liver_fam <- liver[, .(liver_rpm_sum = sum(mean_rpm)), by = fkey]

# Try three keys in order; first that hits the liver panel wins.
ann[, k_family := fk_from_family(mirbase_family)]
ann[, k_symbol := fk_from_symbol(gene_symbol_human)]
liver_lut <- setNames(liver_fam$liver_rpm_sum, liver_fam$fkey)
ann[, rpm_fam := liver_lut[k_family]]
ann[, rpm_sym := liver_lut[k_symbol]]
ann[, liver_rpm_sum := fifelse(!is.na(rpm_fam), rpm_fam, rpm_sym)]
ann[, liver_key_used := fifelse(!is.na(rpm_fam), k_family,
                         fifelse(!is.na(rpm_sym), k_symbol, NA_character_))]

## ============ 4. Host-gene inheritance (GENCODE v49 intragenic same-strand) ============
hostmap <- fread("data/external/mirna_liver_atlas/mirna_host_gene_map.csv")
hep <- fread("Cas13_Library_Design/data/hep_specificity.csv")
hostmap <- merge(hostmap, hep[, .(gene_symbol, host_hep_substrate = hep_substrate)],
                 by.x = "host_name", by.y = "gene_symbol", all.x = TRUE)
ann <- merge(ann, hostmap[, .(mir_name, host_gene = host_name, host_hep_substrate)],
             by.x = "gene_symbol_human", by.y = "mir_name", all.x = TRUE)

## ============ 5. Flag: liver_atlas > host_gene > curated_literature > none ============
LIVER_RPM_THRESHOLD <- 10
ann[, liver_expr_value := round(liver_rpm_sum, 3)]
ann[, liver_expr_unit  := fifelse(!is.na(liver_rpm_sum), "RPM_mean_TCGA_LIHC_normal", NA_character_)]

curated_raw <- c("miR-122","miR-192","miR-194","miR-148a","miR-21","miR-22","miR-26a","miR-26b",
  "miR-30a","miR-30c","miR-30d","miR-30e","miR-99a","miR-99b","miR-27a","miR-27b","miR-378a",
  "miR-451a","miR-126","let-7a","let-7b","let-7c","let-7d","let-7e","let-7f","let-7g","let-7i",
  "miR-101","miR-29a","miR-29b","miR-29c","miR-103a","miR-107","miR-143","miR-145","miR-199a",
  "miR-214","miR-181a","miR-16","miR-23a","miR-23b","miR-24","miR-125b","miR-100","miR-10a",
  "miR-10b","miR-92a","miR-25","miR-93","miR-15a","miR-15b","miR-191","miR-320a")
curated_keys <- unique(c(vapply(curated_raw, fk_from_family, character(1))))
ann[, is_curated_liver := (k_family %in% curated_keys) | (k_symbol %in% curated_keys)]

ann[, host_hep_pass := !is.na(host_hep_substrate) & host_hep_substrate %in% c("high", "ambient_suspect")]
ann[, liver_pass    := !is.na(liver_rpm_sum) & liver_rpm_sum >= LIVER_RPM_THRESHOLD]

ann[, mirna_hep_expressed := FALSE]
ann[, mirna_hep_source    := "none"]
ann[is_curated_liver == TRUE, `:=`(mirna_hep_expressed = TRUE, mirna_hep_source = "curated_literature")]
ann[host_hep_pass    == TRUE, `:=`(mirna_hep_expressed = TRUE, mirna_hep_source = "host_gene")]
ann[liver_pass       == TRUE, `:=`(mirna_hep_expressed = TRUE, mirna_hep_source = "liver_atlas")]

## ============ 6. Final table ============
final <- ann[, .(mouse_ensembl, gene_symbol_human, mirbase_family, mirbase_arm,
                 liver_expr_value, liver_expr_unit, host_gene, host_hep_substrate,
                 mirna_hep_expressed, mirna_hep_source)]
final <- unique(final, by = "mouse_ensembl")
setorder(final, -mirna_hep_expressed, -liver_expr_value, na.last = TRUE)
stopifnot(nrow(final) == nrow(cand))
fwrite(final, "Cas13_Library_Design/data/mirna_hep_expression.csv")

cat("=== SUMMARY ===\n")
cat("total conserved mouse miRNA genes:", nrow(final), "\n\n")
cat("mirna_hep_source distribution:\n"); print(table(final$mirna_hep_source))
cat("\nmirna_hep_expressed:\n"); print(table(final$mirna_hep_expressed))
cat("\nliver_atlas: any RPM mapped:", sum(!is.na(final$liver_expr_value)),
    " | pass (RPM>=10):", sum(final$liver_expr_value >= 10, na.rm = TRUE), "\n")
cat("families still unmapped to liver:", sum(is.na(final$liver_expr_value)), "\n")
cat("\n--- miR-122 row ---\n"); print(final[gene_symbol_human == "MIR122"])
cat("\n--- still-unmapped (no liver RPM) families ---\n")
print(sort(unique(final[is.na(liver_expr_value), mirbase_family])))
cat("\n--- head(8) ---\n"); print(head(final, 8))
