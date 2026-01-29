#!/usr/bin/env Rscript
# setup_r_environment.R
# Installs all required R packages for the combined simulation

cat("Setting up R environment for combined simulation...\n\n")

# Create a local library directory if it doesn't exist
lib_path <- "~/R_libs/pop_uncertainty"
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)

cat("Installing packages to:", lib_path, "\n\n")

# Set library path
.libPaths(c(lib_path, .libPaths()))

# List of required packages
required_packages <- c(
  "ggplot2",
  "tidyr",
  "patchwork",
  "dplyr",
  "knitr",
  "kableExtra",
  "MASS",
  "nimble",
  "coda",
  "extraDistr",
  "igraph",
  "RColorBrewer",
  "nimbleNoBounds",
  "pracma",
  "R.utils"
)

# Function to install package if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, lib = lib_path, repos = "https://cloud.r-project.org/")
    cat(pkg, "installed successfully.\n\n")
  } else {
    cat(pkg, "already installed.\n")
  }
}

# Install all packages
cat("Checking and installing required packages...\n")
cat("This may take 10-20 minutes for first-time setup.\n\n")

for (pkg in required_packages) {
  tryCatch({
    install_if_missing(pkg)
  }, error = function(e) {
    cat("ERROR installing", pkg, ":", conditionMessage(e), "\n")
  })
}

# Verify installation
cat("\n========================================\n")
cat("Verifying package installation...\n")
cat("========================================\n\n")

all_installed <- TRUE
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "\n")
  } else {
    cat("✗", pkg, "- FAILED TO INSTALL\n")
    all_installed <- FALSE
  }
}

cat("\n========================================\n")
if (all_installed) {
  cat("SUCCESS: All packages installed!\n")
  cat("========================================\n\n")
  cat("To use these packages, add this to your R scripts or .Rprofile:\n")
  cat(".libPaths(c('", lib_path, "', .libPaths()))\n", sep = "")
  cat("\nOr set environment variable in your SLURM script:\n")
  cat("export R_LIBS_USER=", lib_path, "\n", sep = "")
} else {
  cat("WARNING: Some packages failed to install.\n")
  cat("========================================\n")
  cat("You may need to:\n")
  cat("1. Check if dependencies are available on your system\n")
  cat("2. Load additional modules (e.g., GSL, GDAL)\n")
  cat("3. Install system-level dependencies\n")
}

cat("\nLibrary path set to:", lib_path, "\n")
cat("Session info:\n")
print(sessionInfo())
