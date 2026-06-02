if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="http://cran.us.r-project.org")

BiocManager::install(c("slingshot", "tradeSeq", "viper", "progeny", "signatureSearch", "ExperimentHub"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools", repos="http://cran.us.r-project.org")

options(timeout=9999999)
# MR packages removed 2026-04-22 (MR ditched from paper; see archive/mr_ditched_2026-04-22/)
# devtools::install_github("MRCIEU/TwoSampleMR")
# devtools::install_github("WSpiller/MVMR")
devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")

install.packages("coloc", repos="http://cran.us.r-project.org")

print("Installation script parsed successfully")
