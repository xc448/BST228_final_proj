# here we define two functions to calculate MCMC diagnostics for valid parameters 
# and to plot MCMC diagnostics for a specific sample and cell type.

library(coda)

mcmc_diagnostics_metrics <- function(theta_n, theta_n_samples, 
                                            sample_names = NULL, cell_type_names = NULL, 
                                            threshold = 1e-8) {
  # Args:
  #   theta_n: A 3D array (N x K x chains) where proportions are stored.
  #   theta_n_samples: A 4D array (iter, N, K, chains) with MCMC samples.
  #   sample_names: A vector of sample names (length N).
  #   cell_type_names: A vector of cell type names (length K).
  #   threshold: Values below this threshold are treated as zero.
  #
  # Returns:
  #   A list of diagnostics (Rhat, ESS) for non-zero parameters from invalid expression data.
  
  # Step 0: Assign real names
  if (is.null(sample_names)) sample_names <- dimnames(theta_n)[[1]]
  if (is.null(cell_type_names)) cell_type_names <- dimnames(theta_n)[[2]]
  if (is.null(sample_names) || is.null(cell_type_names)) stop("Sample and cell type names required.")
  
  # Step 1: Identify non-zero proportions
  nonzero_indices <- which(apply(theta_n, c(1, 2), mean) > threshold, arr.ind = TRUE)
  
  if (nrow(nonzero_indices) == 0) {
    message("No non-zero proportions found.")
    return(NULL)
  }
  
  # Step 2: Run diagnostics for non-zero proportions
  diagnostics <- list()
  
  for (i in 1:nrow(nonzero_indices)) {
    n_idx <- nonzero_indices[i, 1]
    k_idx <- nonzero_indices[i, 2]
    sample_name <- sample_names[n_idx]
    cell_type_name <- cell_type_names[k_idx]
    
    # Extract samples across all chains for the parameter
    samples_list <- list()
    for (ch in 1:dim(theta_n_samples)[4]) {
      chain_samples <- theta_n_samples[, n_idx, k_idx, ch]
      chain_samples <- chain_samples[is.finite(chain_samples)]  # Clean invalid values
      if (length(chain_samples) > 0) {
        samples_list[[ch]] <- mcmc(chain_samples)  # Convert to MCMC object
      }
    }
    
    # Check if valid chains exist
    if (length(samples_list) < 2) {
      warning(sprintf("Not enough valid chains for Sample '%s', Cell Type '%s'.", sample_name, cell_type_name))
      next
    }
    
    # Combine chains into an mcmc.list
    mcmc_list <- mcmc.list(samples_list)
    
    # Calculate diagnostics safely
    rhat <- tryCatch({
      rhat_res <- gelman.diag(mcmc_list, autoburnin = FALSE)$psrf
      rhat_res[1, "Point est."]  # Extract the point estimate
    }, error = function(e) NA)
    
    ess <- tryCatch(effectiveSize(mcmc_list), error = function(e) NA)
    
    # Store diagnostics
    diag_key <- paste0(sample_name, "_", cell_type_name)
    diagnostics[[diag_key]] <- list(
      rhat = rhat,
      effective_size = ess
    )
  }
  
  return(diagnostics)
}


plot_mcmc_diagnostics <- function(res, sample_name, cell_type_name, sample_names, cell_type_names) {
  # Args:
  #   res: The result list containing 'theta_n_samples' (MCMC samples array).
  #   sample_name: The name of the sample to plot.
  #   cell_type_name: The name of the cell type to plot.
  #   sample_names: Vector of sample names (length N).
  #   cell_type_names: Vector of cell type names (length K).
  #
  # Generates trace plots and ACF plots for the specified sample and cell type.
  
  # Find indices for sample and cell type
  n_idx <- which(sample_names == sample_name)
  k_idx <- which(cell_type_names == cell_type_name)
  
  # Extract MCMC samples for the selected parameter
  num_kept <- dim(res$theta_n_samples)[1]  # Number of kept samples
  chains <- dim(res$theta_n_samples)[4]    # Number of chains
  param_samples <- array(NA, dim = c(num_kept, chains))
  
  for (ch in 1:chains) {
    param_samples[, ch] <- res$theta_n_samples[, n_idx, k_idx, ch]
  }
  
  # Convert to mcmc.list
  mcmc_chains <- lapply(1:chains, function(ch) mcmc(param_samples[, ch]))
  mcmc_list <- mcmc.list(mcmc_chains)
  
  # Open a new device to ensure the plot is displayed
  dev.new()  # Works interactively in RStudio or R terminal
  
  # Generate Trace Plot
  par(mfrow = c(1, 2))  # Set layout: 1 row, 2 columns
  plot(mcmc_list, main = sprintf("Trace Plot: %s - %s", sample_name, cell_type_name),density=F)
  
  # Generate ACF Plot for the first chain
  acf(as.numeric(mcmc_list[[1]]), main = sprintf("ACF Plot: %s - %s (Chain 1)", sample_name, cell_type_name))
  
  par(mfrow = c(1, 1))  # Reset layout
}
