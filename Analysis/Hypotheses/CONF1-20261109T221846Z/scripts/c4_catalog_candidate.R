# CONF1 C4 (PRESPEC_C3_C4.md): candidate Catalog fields. Reads the adopted table; writes a new candidate table.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
if (is.na(OUT) || dir.exists(OUT)) stop("output directory missing or already exists")
dir.create(OUT)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CAT <- file.path(ROOT, "RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/catalog_annotated_gene_table.tsv.gz")
LIAB <- file.path(ROOT, "RNA-seq/results/evidence_passports/candidates/gene-catalog-r8-anchor-disease-sidecar-20261025T131731Z/catalog_r8_anchor_disease_liability.tsv")
wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")
wr(data.table(path = c(CAT, LIAB), sha256 = vapply(c(CAT, LIAB), function(f) sub(" .*", "", system2("sha256sum", shQuote(f), stdout = TRUE)), "")), "inputs_sha256.tsv")
x <- fread(CAT, select = c("ensembl_id", "symbol", "primary_evidence_class", "genetics_call_state", "genetics_testability_state",
                           "r8_rna_call_state", "rna_canonical_call", "rna_log2fc", "rna_bh_q", "coding_cs_max_pip",
                           "coding_cs_top_variant", "coding_cs_top_study", "coding_cs_top_consequence", "coding_cs_top_trait_tier"))
stopifnot(nrow(x) == 28233, !anyDuplicated(x$ensembl_id))
li <- fread(LIAB, select = c("ensembl_id", "sidecar_class", "anchor_trait"))
stopifnot(!anyDuplicated(li$ensembl_id))
x[, regulatory := genetics_call_state == "supported"]
x[, coding := !is.na(coding_cs_max_pip) & coding_cs_max_pip >= 0.5]
x[, state_call := r8_rna_call_state == "supported"]
x[, genetic_mechanism_candidate := fifelse(regulatory & coding, "both", fifelse(regulatory, "regulatory_colocalization",
                                   fifelse(coding, "coding_fine_mapped", "none")))]
x[, candidate_class := primary_evidence_class]
x[primary_evidence_class == "concordant", candidate_class := "genetic_and_state"]
x[coding & primary_evidence_class == "unresolved", candidate_class := "genetically_anchored"]
x[coding & primary_evidence_class == "established_state_associated", candidate_class := "genetic_and_state"]
x <- merge(x, li, by = "ensembl_id", all.x = TRUE)
x[, disease_liability_candidate := fifelse(!is.na(sidecar_class), sidecar_class,
  fifelse(coding & coding_cs_top_trait_tier == "tier1_direct_masld_pdff", "coding_masld_or_liver_fat_gwas",
  fifelse(coding, "coding_liver_enzyme_gwas_only", "not_assessed")))]
x[, class_changed := candidate_class != primary_evidence_class & primary_evidence_class != "concordant"]
fwrite(x, file.path(OUT, "catalog_candidate_fields.tsv.gz"), sep = "\t", na = "NA")
wr(x[, .N, by = .(adopted = primary_evidence_class, candidate = candidate_class)][order(adopted, candidate)], "class_transition.tsv")
wr(x[class_changed == TRUE | coding == TRUE, .(symbol, adopted = primary_evidence_class, candidate = candidate_class, genetic_mechanism_candidate,
     disease_liability_candidate, coding_cs_max_pip, coding_cs_top_variant, coding_cs_top_consequence, coding_cs_top_trait_tier,
     rna_log2fc, rna_bh_q)][order(candidate, symbol)], "reclassified_genes.tsv")
KNOWN <- c("PNPLA3", "TM6SF2", "MTARC1", "GCKR", "APOE", "SERPINA1", "HSD17B13", "MBOAT7", "GPAM", "MTTP", "CIDEB", "PSD3", "TRIB1",
           "APOB", "PPP1R3B", "GLP1R", "THRB")
wr(x[symbol %in% KNOWN, .(symbol, adopted = primary_evidence_class, candidate = candidate_class, genetic_mechanism_candidate,
     disease_liability_candidate, coding_cs_max_pip, rna_log2fc, rna_bh_q)][match(KNOWN, symbol, nomatch = 0)], "known_gene_check.tsv")
wr(x[candidate_class %in% c("genetically_anchored", "genetic_and_state"), .N, by = .(candidate_class, genetic_mechanism_candidate,
     disease_liability_candidate)][order(candidate_class, -N)], "liability_by_class.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
