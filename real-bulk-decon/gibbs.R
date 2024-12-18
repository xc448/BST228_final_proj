rdirichlet <- function(alpha) {
  x <- rgamma(length(alpha), alpha)
  return(x / sum(x))
}

library(parallel)
library(coda)  # For MCMC diagnostics (if not installed, install.packages("coda"))

# Single-sample, single-chain sampler
single_sample_sampler <- function(X_n, phi, alpha, n_iter, warmup = 0, thin = 1, compute.elbo = FALSE, save_Z = FALSE) {
  G <- length(X_n)
  K <- nrow(phi)
  
  # Compute how many samples will actually be kept
  num_kept <- floor((n_iter - warmup) / thin)
  
  # Initial values: start at 1/K
  theta_n_i <- rep(1 / K, K)
  Z_n_i <- matrix(NA, nrow = G, ncol = K)
  
  if (num_kept > 0) {
    theta_samples <- matrix(NA, nrow = num_kept, ncol = K)
    Z_samples <- if (save_Z) array(NA, dim = c(num_kept, G, K)) else NULL
  } else {
    theta_samples <- matrix(NA, nrow = 0, ncol = K)
    Z_samples <- if (save_Z) array(NA, dim = c(0, G, K)) else NULL
  }
  
  theta_n_sum <- rep(0, K)
  theta_n2_sum <- rep(0, K)
  Z_n_sum <- matrix(0, nrow = G, ncol = K)
  multinom_coef_sum <- 0
  
  max_size <- 1e6
  keep_idx <- 0
  
  for (iter in 1:n_iter) {
    # Sample Z
    prob.mat <- phi * theta_n_i
    for (g in 1:G) {
      col_prob <- prob.mat[, g]
      psum <- sum(col_prob)
      if (psum <= 0 || is.na(psum)) {
        col_prob <- rep(1/K, K)
      } else {
        col_prob <- col_prob / psum
      }
      
      size_for_multinom <- X_n[g]
      if (size_for_multinom > max_size) {
        scale_factor <- max_size / size_for_multinom
        prob.tmp <- col_prob * scale_factor
        ptmp_sum <- sum(prob.tmp)
        if (ptmp_sum <= 0 || is.na(ptmp_sum)) {
          prob.tmp <- rep(1/K, K)
        } else {
          prob.tmp <- prob.tmp / ptmp_sum
        }
        Z_n_i[g, ] <- as.numeric(rmultinom(n = 1, size = max_size, prob = prob.tmp))
      } else {
        Z_n_i[g, ] <- as.numeric(rmultinom(n = 1, size = size_for_multinom, prob = col_prob))
      }
    }
    
    # Sample theta
    Z_nk_i <- colSums(Z_n_i)
    theta_n_i <- rdirichlet(Z_nk_i + alpha)
    
    if (iter > warmup && ((iter - warmup) %% thin == 0)) {
      keep_idx <- keep_idx + 1
      theta_samples[keep_idx, ] <- theta_n_i
      if (save_Z) Z_samples[keep_idx, , ] <- Z_n_i
      
      # Accumulate sums
      theta_n_sum <- theta_n_sum + theta_n_i
      theta_n2_sum <- theta_n2_sum + theta_n_i^2
      Z_n_sum <- Z_n_sum + Z_n_i
      
      if (compute.elbo) {
        mcoef_num <- sum(lfactorial(Z_nk_i))
        mcoef_den <- sum(lfactorial(Z_n_i))
        
        if (any(is.na(mcoef_num), is.na(mcoef_den)) ||
            any(is.infinite(mcoef_num), is.infinite(mcoef_den))) {
          multinom_coef_i <- 0
        } else {
          multinom_coef_i <- mcoef_num - mcoef_den
        }
        
        if (is.nan(multinom_coef_i) || is.infinite(multinom_coef_i)) {
          multinom_coef_i <- 0
        }
        
        multinom_coef_sum <- multinom_coef_sum + multinom_coef_i
      }
    }
  }
  
  if (num_kept > 0) {
    theta_n_mean <- theta_n_sum / num_kept
    var_theta <- (theta_n2_sum / num_kept) - (theta_n_mean^2)
    theta_cv_val <- sqrt(var_theta) / theta_n_mean
    gibbs_constant <- multinom_coef_sum / num_kept
    Z_n_mean <- Z_n_sum / num_kept
  } else {
    theta_n_mean <- rep(NA, K)
    theta_cv_val <- rep(NA, K)
    gibbs_constant <- NA
    Z_n_mean <- matrix(NA, nrow = G, ncol = K)
  }
  
  list(
    theta_n = theta_n_mean,
    theta_samples = theta_samples,
    Z_samples = Z_samples,
    Z_n_mean = Z_n_mean,
    theta_cv = theta_cv_val,
    gibbs_constant = gibbs_constant
  )
}

# Function to loop over all bulk samples and multiple chains
sample.Z.theta_multi <- function(X, phi, alpha, n_iter, warmup = 0, thin = 1,
                                 n_cores = 2, chains = 1, compute.elbo = FALSE, save_Z = FALSE) {
  G <- nrow(X)
  N <- ncol(X)
  K <- nrow(phi)
  
  if (nrow(phi) != K || ncol(phi) != G) stop("Dimension mismatch: phi must have K rows and G columns.")
  if (nrow(X) != G) stop("Dimension mismatch: X must have G rows.")
  
  phi <- t(t(phi) / colSums(phi))
  
  num_kept <- floor((n_iter - warmup) / thin)
  
  if (num_kept > 0) {
    theta_n_samples <- array(NA, dim = c(num_kept, N, K, chains),
                             dimnames = list(NULL, colnames(X), rownames(phi), NULL))
  } else {
    theta_n_samples <- array(NA, dim = c(0, N, K, chains))
  }
  
  theta_n <- array(NA, dim = c(N, K, chains),
                   dimnames = list(colnames(X), rownames(phi), NULL))
  theta_cv <- array(NA, dim = c(N, chains), dimnames = list(colnames(X), NULL))
  gibbs_constants <- array(NA, dim = c(N, chains), dimnames = list(colnames(X), NULL))
  
  if (save_Z) {
    if (num_kept > 0) {
      Z_n_samples <- array(NA, dim = c(num_kept, G, K, N, chains),
                           dimnames = list(NULL, rownames(X), rownames(phi), colnames(X), NULL))
    } else {
      Z_n_samples <- array(NA, dim = c(0, G, K, N, chains))
    }
  } else {
    Z_n_samples <- NULL
  }
  
  Z_n_mean_all <- array(NA, dim = c(G, K, N, chains),
                        dimnames = list(rownames(X), rownames(phi), colnames(X), NULL))
  
  results <- mclapply(seq_len(N), function(n_idx) {
    chain_results <- lapply(seq_len(chains), function(chain_idx) {
      message(sprintf("Processing sample %d/%d, chain %d/%d", n_idx, N, chain_idx, chains))
      single_sample_sampler(X[, n_idx], phi, alpha, n_iter, warmup, thin, compute.elbo, save_Z)
    })
    chain_results
  }, mc.cores = n_cores)
  
  # Combine results
  for (n_idx in seq_len(N)) {
    for (chain_idx in seq_len(chains)) {
      res <- results[[n_idx]][[chain_idx]]
      if (!is.null(res)) {
        theta_n[n_idx, , chain_idx] <- res$theta_n
        theta_cv[n_idx, chain_idx] <- mean(res$theta_cv, na.rm = TRUE)
        gibbs_constants[n_idx, chain_idx] <- res$gibbs_constant
        Z_n_mean_all[, , n_idx, chain_idx] <- res$Z_n_mean
        
        if (num_kept > 0) {
          theta_n_samples[, n_idx, , chain_idx] <- res$theta_samples
          if (save_Z) {
            Z_n_samples[, , , n_idx, chain_idx] <- res$Z_samples
          }
        }
      }
    }
  }
  
  output <- list(
    theta_n = theta_n,
    theta_n_samples = theta_n_samples,
    theta_cv = theta_cv,
    gibbs_constants = gibbs_constants,
    Z_n = Z_n_mean_all
  )
  
  if (save_Z) {
    output$Z_n_samples <- Z_n_samples
  }
  
  return(output)
}