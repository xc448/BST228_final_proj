# here define a few functions to help with data manipulation
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(reshape2)

### Function to process sample proportions ###

process_sample_proportions <- function(df) {
  # Ensure sample names are accessible as a column
  df <- df %>%
    rownames_to_column(var = "sample")
  
  # Reshape the data into long format
  df_long <- df %>%
    pivot_longer(-sample, names_to = "cell_ref", values_to = "proportion") %>%
    mutate(
      # Extract cell type from column names (e.g., Astrocytes from "Astrocytes_ADxxxx")
      cell_type = str_extract(cell_ref, "^[^_]+"),
      # Extract reference ID (e.g., ADxxxx from "Astrocytes_ADxxxx")
      reference = str_extract(cell_ref, "AD[0-9]+")
    )
  
  return(df_long)  # Return the reshaped data without applying filters
}


### Function to dynamically group cell types by prefix ###
aggregate_by_prefix <- function(Z_n, patterns) {
  aggregated <- list()
  for (pattern in patterns) {
    # Identify cell types matching the current pattern
    idx <- grep(paste0("^", pattern), dimnames(Z_n)[[2]])  # Match cell types starting with pattern
    if (length(idx) > 0) {
      # Sum across the matching cell types
      aggregated[[pattern]] <- apply(Z_n[, idx, , drop = FALSE], c(1, 3), sum)
    }
  }
  return(aggregated)
}


### Function to reshape a list of matrices into a tidy data frame ###

reshape_expression_list <- function(expression_list) {
  # Args:
  #   expression_list: A named list of matrices, where each matrix contains gene expression values.
  #
  # Returns:
  #   A tidy data frame with columns: Gene, Sample, CellType, and Expression.
  
  # Combine list of matrices into a single data frame
  reshaped_df <- do.call(rbind, lapply(names(expression_list), function(cell_type) {
    data.frame(
      Gene = rep(rownames(expression_list[[cell_type]]), ncol(expression_list[[cell_type]])),
      Sample = rep(colnames(expression_list[[cell_type]]), each = nrow(expression_list[[cell_type]])),
      CellType = cell_type,
      Expression = as.vector(expression_list[[cell_type]])
    )
  }))
  
  return(reshaped_df)
}

