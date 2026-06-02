#!/usr/bin/env Rscript
# Install InstaPrism and BayesPrism packages
# Run this first to set up the environment

# Force conda library path to prevent user library contamination
conda_lib <- Sys.getenv("R_LIBS_SITE")
if (nzchar(conda_lib) && dir.exists(conda_lib)) {
  .libPaths(conda_lib)
}

# Check if the packages are already installed
packages_needed <- c("InstaPrism", "BayesPrism")
packages_installed <- sapply(packages_needed, requireNamespace, quietly = TRUE)

if (!all(packages_installed)) {
  # Install devtools if needed
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  
  # Install InstaPrism (faster alternative to BayesPrism)
  if (!packages_installed["InstaPrism"]) {
    cat("Installing InstaPrism...\n")
    tryCatch({
      devtools::install_github("humengying0907/InstaPrism", upgrade = "never")
      cat("InstaPrism installed successfully!\n")
    }, error = function(e) {
      cat("Failed to install InstaPrism:", conditionMessage(e), "\n")
    })
  }
  
  # Install BayesPrism (original, slower)
  if (!packages_installed["BayesPrism"]) {
    cat("Installing BayesPrism...\n")
    tryCatch({
      devtools::install_github("Danko-Lab/BayesPrism", upgrade = "never")
      cat("BayesPrism installed successfully!\n")
    }, error = function(e) {
      cat("Failed to install BayesPrism:", conditionMessage(e), "\n")
    })
  }
}

# Verify installation
cat("\n=== Verification ===\n")
for (pkg in packages_needed) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, ": installed\n")
  } else {
    cat(pkg, ": NOT FOUND\n")
  }
}
