suppressPackageStartupMessages({library(data.table);library(readxl)})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V <- file.path(BASE,"data/external/vacca_2024")

# MOESM6 NES
p6 <- as.data.table(read_excel(file.path(V,"42255_2024_1043_MOESM6_ESM.xlsx"),sheet="NES"))
cat("=== MOESM6 NES dims:",dim(p6),"===\n")
cat("col1 name:",names(p6)[1],"\n")
cat("ALL column names:\n"); print(names(p6))
cat("\nN pathways (rows):",nrow(p6),"\n")
cat("\nfirst 20 pathway names (col1):\n"); print(head(p6[[1]],20))
cat("\nALL pathway names matching metabolic keywords:\n")
pw <- p6[[1]]
mk <- grepl("LIPID|GLUCOS|INSULIN|PPAR|FATTY|GLYCO|CHOLESTEROL|METABOL|BILE|STEROID|FAT|GLUCONEO|PYRUVATE|CITRATE|OXIDATIV|BETA.?ALAN|ADIPO|PROPANOATE|BUTANOATE|RETINOL|ARACHIDONIC|LINOLEIC|TRYPTOPHAN|TYROSINE",toupper(pw))
print(pw[mk])
cat("\nN metabolic-keyword pathways:",sum(mk),"\n")

# MOESM9 DHPS
m9 <- as.data.table(read_excel(file.path(V,"42255_2024_1043_MOESM9_ESM.xlsx"),col_names=FALSE,skip=2))
cat("\n=== MOESM9 DHPS dims:",dim(m9),"===\n")
cat("first 6 rows, first 10 cols:\n")
print(m9[1:6,1:min(10,ncol(m9))])
