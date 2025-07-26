```R
library(glmnet)
library(mclust)
library(dplyr)
library(tidyr)

methylation_df <- read.csv(file = 'Model/final_m5C_data.csv', row.names = 1)
expression_df <- read.csv(file = 'Model/expr_data.csv', row.names = 1)
normal_samples <- c('nor_1', 'nor_2', 'nor_3', 'nor_4')
cancer_samples <- c('can_1', 'can_2', 'can_3', 'can_4')

output_file <- "Model/predictive_genes_output.txt"

write_output <- function(file, text) {
  write(text, file = file, append = TRUE)
}

identify_predictive_genes <- function(methylation_df, expression_df) {
  predictive_genes <- c()
  gene_names <- rownames(methylation_df)
  file.remove(output_file)
  
  for (gene in gene_names) {
    gene_methylation <- as.numeric(methylation_df[gene, ])
    expression_levels <- as.numeric(expression_df[gene, ])
    if (length(gene_methylation) == length(expression_levels)) {
      output_text <- paste("Gene:", gene, "\n")
      output_text <- paste(output_text, "Methylation:\n")
      output_text <- paste(output_text, paste(sprintf("%.3f", gene_methylation), collapse = " "), "\n")
      output_text <- paste(output_text, "Expression:\n")
      output_text <- paste(output_text, paste(sprintf("%.3f", expression_levels), collapse = " "), "\n")
      
      if (!all(is.na(expression_levels)) && !all(expression_levels == 0)) {
        
        if (sd(gene_methylation) != 0 && sd(expression_levels) != 0) {
          correlation <- cor(gene_methylation, expression_levels, use = "complete.obs")
          output_text <- paste(output_text, sprintf("Correlation: %.3f\n", correlation))
         
          if (!is.na(correlation) && correlation > 0.6) {
            predictive_genes <- c(predictive_genes, gene)
            output_text <- paste(output_text, "Predictive Gene: Added\n")
          }
        } else {
          output_text <- paste(output_text, paste("Skipping gene due to zero standard deviation:", gene, "\n"))
        }
      }
      
      write_output(output_file, output_text)
      
    } else {
      write_output(output_file, paste("Skipping gene due to dimension mismatch:", gene, "\n"))
    }
  }
  return(predictive_genes)
}

identify_methylation_states <- function(methylation_df, predictive_genes) {
  methylation_states <- list()
  
  for (gene in predictive_genes) {
    gene_methylation <- as.numeric(methylation_df[gene, ])
    
    print(paste("Processing gene:", gene))
    print("Methylation data:")
    print(gene_methylation)

    gene_methylation <- as.vector(gene_methylation)
    
    if (length(gene_methylation) < 2) {
      next
    }
    
    tryCatch({
      fit <- Mclust(gene_methylation, G = 1:3)
      methylation_states[[gene]] <- fit
    }, error = function(e) {
      print(paste("Error in Mclust for gene:", gene))
      print(e)
    })
  }
  
  return(methylation_states)
}


formaldehyde_weight_function <- function(state_indices, gene, gene_state, cancer_samples) {
  cancer_state_methylation <- gene_state$data[state_indices]
  state_mean_methylation <- mean(cancer_state_methylation, na.rm = TRUE)
  state_complexity <- length(unique(gene_state$classification))
  weight <- max(1, log10(state_mean_methylation + 1)) * state_complexity
  return(weight)
}

calculate_dm_values_weighted <- function(methylation_df, normal_samples, cancer_samples, genes, methylation_states, weight_function) {
  dm_values <- list()
  
  for (gene in genes) {
    if (!(gene %in% rownames(methylation_df))) {
      warning(paste("Gene", gene, "not found in methylation_df"))
      dm_values[[gene]] <- NA
      next
    }
    normal_means <- rowMeans(methylation_df[gene, normal_samples, drop = FALSE], na.rm = TRUE)
    
    if (length(normal_means) == 0) {
      warning(paste("No normal samples for gene", gene))
      dm_values[[gene]] <- NA
      next
    }
    
    if (!(gene %in% names(methylation_states))) {
      warning(paste("Methylation states for gene", gene, "not found"))
      dm_values[[gene]] <- NA
      next
    }
    
    states <- unique(methylation_states[[gene]]$classification)

    weighted_dm_value <- 0
    total_weight <- 0
    
    for (state in states) {
      state_indices <- which(methylation_states[[gene]]$classification == state)
      
      if (length(state_indices) == 0) {
        next
      }

      cancer_state_samples <- cancer_samples[state_indices]

      valid_samples <- cancer_state_samples[cancer_state_samples %in% colnames(methylation_df)]
      
      if (length(valid_samples) == 0) {
        next
      }
      cancer_state_samples_state <- methylation_df[gene, valid_samples, drop = FALSE]
      cancer_state_samples_state_str <- capture.output(print(cancer_state_samples_state))
      cat("State:", state, "\n")
      cat("State indices:", state_indices, "\n")
      cat("Cancer state samples:", valid_samples, "\n")
      cat("Cancer state samples data:\n")
      cat(cancer_state_samples_state_str, sep = "\n")
      cancer_state_max <- apply(cancer_state_samples_state, 2, function(x) max(x, na.rm = TRUE))
      
      if (length(cancer_state_max) == 0 || all(is.na(cancer_state_max))) {
        next
      }
      state_dm_value <- mean(cancer_state_max, na.rm = TRUE) - mean(normal_means, na.rm = TRUE)
      state_weight <- weight_function(state_indices, gene, methylation_states[[gene]], cancer_samples)
      
      weighted_dm_value <- weighted_dm_value + (state_dm_value * state_weight)
      total_weight <- total_weight + state_weight
    }
    if (total_weight > 0) {
      dm_values[[gene]] <- weighted_dm_value / total_weight
    } else {
      dm_values[[gene]] <- NA
    }
  }
  
  return(dm_values)
}

predictive_genes <- identify_predictive_genes(methylation_df, expression_df)
methylation_states <- identify_methylation_states(methylation_df, predictive_genes)
dm_values <- calculate_dm_values_weighted(methylation_df, normal_samples, cancer_samples, predictive_genes, methylation_states, formaldehyde_weight_function)

head(predictive_genes)
head(methylation_states)
head(dm_values)
```

