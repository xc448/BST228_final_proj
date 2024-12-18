#!/usr/bin/env Rscript

# Set up logging
log_file <- "./output/decon_case.log"
error_file <- "./output/decon_case_errors.log"
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

cat("Starting deconvolution script...\n")

# Set working directory
setwd("~/Documents/R_working/bst228")  # Adjust as needed
if (!dir.exists("./output")) dir.create("./output")

# Load data
cat("Loading input data...\n")
Exp_case <- readRDS("./data/Exp_case.rds")
ref_case <- readRDS("./data/ref_case.rds")

cat("Validating input dimensions...\n")
if (ncol(ref_case) != nrow(Exp_case)) {
  stop("Dimension mismatch: ref_case rows must match Exp_case columns.")
}

# Source Gibbs sampler
cat("Sourcing Gibbs sampler function...\n")
source("./scripts/gibbs.R")

# Run Gibbs sampler
cat("Running Gibbs sampler...\n")
result <- tryCatch({
  sample.Z.theta_multi(
    X = Exp_case, 
    phi = ref_case, 
    alpha = 10e-8, 
    n_iter = 5000, 
    warmup = 2500, 
    thin = 2, 
    chains = 5, 
    n_cores = 10, 
    compute.elbo = F, 
    save_Z = F
  )
}, error = function(e) {
  message("Error occurred during Gibbs sampling: ", e$message)
  stop(e)
})

# Save results
cat("Saving results...\n")
saveRDS(result, file = "./output/bulk_case_output.rds")
cat("Deconvolution completed successfully!\n")