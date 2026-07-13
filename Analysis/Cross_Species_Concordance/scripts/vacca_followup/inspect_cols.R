suppressPackageStartupMessages({library(readxl);library(data.table)})
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
s4 <- as.data.table(read_excel(file.path(BASE,"data/external/vacca_2024/42255_2024_1043_MOESM4_ESM.xlsx"), sheet="Table S4"))
cat("S4 ncol", ncol(s4), "nrow", nrow(s4), "\n")
cat("--- all col names ---\n"); print(names(s4))
cat("\n--- human progression cols ---\n")
print(grep("UCAM|EPoS|Severe|Mild|Moderate", names(s4), value=TRUE))
