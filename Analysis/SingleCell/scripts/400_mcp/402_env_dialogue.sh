#!/bin/bash
# 402_env_dialogue.sh — install DIALOGUE (Jerby-Arnon 2022) + PMA + MOFA2 into celltype_bio env.
# celltype_bio already has R 4.4.3, lme4, lmerTest, Seurat, decoupleR, edgeR, limma.
# Adds: DIALOGUE (github: livnatje/DIALOGUE), PMA (Witten sparse PCA for DIALOGUE's PMD step),
# MOFA2 (benchmark against DIALOGUE), matrixStats (already present), plyr, psych, unikn, mixtools.
set -euo pipefail

ENV_NAME="${ENV_NAME:-celltype_bio}"
echo "[402] Installing DIALOGUE + PMA + MOFA2 into: ${ENV_NAME}"

micromamba run -n "${ENV_NAME}" R --no-save <<'EOF'
options(Ncpus = 4, repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c("PMA", "plyr", "psych", "unikn", "mixtools", "statmod", "remotes")
missing <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("MOFA2", quietly = TRUE)) BiocManager::install("MOFA2", update = FALSE, ask = FALSE)

if (!requireNamespace("DIALOGUE", quietly = TRUE)) {
  remotes::install_github("livnatje/DIALOGUE", upgrade = "never", quiet = FALSE)
}

cat("\n---- versions ----\n")
for (p in c("DIALOGUE", "PMA", "lme4", "lmerTest", "MOFA2", "matrixStats", "plyr", "psych", "Seurat", "edgeR", "limma", "decoupleR")) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("%-14s %s\n", p, as.character(packageVersion(p))))
  } else {
    cat(sprintf("%-14s NOT INSTALLED\n", p))
  }
}

library(DIALOGUE)
cat("\nDIALOGUE loaded OK.\n")
EOF

echo "[402] DONE."
