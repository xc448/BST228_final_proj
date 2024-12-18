#!/usr/bin/env Rscript

# Set up logging
log_file <- "./output/decon_control.log"
error_file <- "./output/decon_control_errors.log"
log_con <- file(log_file, open = "wt")
error_con <- file(error_file, open = "wt")

sink(log_con, append = TRUE)
sink(error_con, type = "message", append = TRUE)
on.exit({
  sink()
  sink(type = "message")
  close(log_con)
  close(error_con)
})

cat("Starting control deconvolution script...\n")

# Set working directory
setwd("~/Documents/R_working/bst228")  # Adjust as needed
if (!dir.exists("./output")) dir.create("./output")

# Load control data
cat("Loading input data for control samples...\n")
if (!file.exists("./data/Exp_control.rds") || !file.exists("./data/ref_control.rds")) {
  stop("Required input files for control samples not found in './data'. Ensure 'Exp_control.rds' and 'ref_control.rds' are present.")
}
Exp_control <- readRDS("./data/Exp_control.rds")
ref_control <- readRDS("./data/ref_control.rds")

# Check data structure
cat("Inspecting data structure for control samples...\n")
cat("Exp_control dimensions: ", dim(Exp_control), "\n")
cat("ref_control dimensions: ", dim(ref_control), "\n")

# Validate input dimensions
cat("Validating input dimensions...\n")
if (ncol(ref_control) != nrow(Exp_control)) {
  stop("Dimension mismatch: Number of rows in 'ref_control' must match the number of columns in 'Exp_control'.")
}

# Source Gibbs sampler
cat("Sourcing Gibbs sampler function...\n")
if (!file.exists("./scripts/gibbs.R")) {
  stop("The Gibbs sampler script './scripts/gibbs.R' not found. Ensure the file exists and is accessible.")
}
source("./scripts/gibbs.R")

# Run Gibbs sampler
cat("Running Gibbs sampler for control samples...\n")
result <- tryCatch({
  sample.Z.theta_multi(
    X = Exp_control, 
    phi = ref_control, 
    alpha = 10e-8,    # Dirichlet prior
    n_iter = 5000, 
    warmup = 2500, 
    thin = 2, 
    n_cores = 10, 
    chains = 5, 
    compute.elbo = F, 
    save_Z = F
  )
}, error = function(e) {
  message("Error occurred during Gibbs sampling for control samples: ", e$message)
  stop(e)
})

# Save results
cat("Saving results for control samples...\n")
saveRDS(result, file = "./output/bulk_control_output.rds")
cat("Control deconvolution completed successfully!\n")