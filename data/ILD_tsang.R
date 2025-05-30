# Load required libraries
library(ape)          # Phylogenetic trees & distance matrices
library(phangorn)     # UPGMA clustering & RF distance
library(ggplot2)      # Visualization
library(Biostrings)   # FASTA file parsing
library(seqinr)       # Alternative sequence handling
library(tidyverse)    # Data manipulation

ILD_analysis <- function(infile, output_directory = ".") {
  
  # Split input filenames
  files <- unlist(strsplit(infile, ","))
  fas1 <- files[1]
  fas2 <- files[2]
  
  # Read sequences
  if (!file.exists(fas1) | !file.exists(fas2)) {
    stop("Error: One or both input FASTA files do not exist.")
  }

  sequences1 <- read.dna(fas1, format = "fasta")
  sequences2 <- read.dna(fas2, format = "fasta")

  # Ensure valid sequences
  if (nrow(sequences1) == 0 | nrow(sequences2) == 0) {
    stop("Error: One or both input files contain no valid sequences.")
  }
  
  # Number of taxa
  taxon_number <- nrow(sequences1)
  
  # Define splitting points
  s1 <- round(taxon_number * 1 / 3)
  s2 <- round(taxon_number * 2 / 4)
  s3 <- taxon_number  # Total sequences
  
  cat("Your sequence is processing into 3 subsets:\n")
  cat(s1, " separated (1/3 of total sequences)\n")
  cat(s2, " separated (1/2 of total sequences)\n")
  cat(s3, " total sequences\n\n")
  
  split_sizes <- c(s1, s2, s3)
  results <- data.frame(
    Matching_Trees = numeric(),
    Probability_Match = numeric(),
    Phylo_Info_Content = numeric(),
    Num_Sequences = numeric()
  )
  
  # Process each subset
  for (split_size in split_sizes) {
    n_splits <- as.integer(split_size)
    
    # Skip iteration if split size is invalid
    if (n_splits <= 1) next

    # Subset sequences
    splits1 <- sequences1[1:n_splits, , drop = FALSE]
    splits2 <- sequences2[1:n_splits, , drop = FALSE]
    
    # Compute pairwise distances
    distances1 <- dist.dna(splits1, model = "raw")
    distances2 <- dist.dna(splits2, model = "raw")
    
    # Build trees using UPGMA
    tree1 <- upgma(distances1)
    tree2 <- upgma(distances2)
    
    # Ensure trees are valid
    if (!inherits(tree1, "phylo") | !inherits(tree2, "phylo")) {
      warning("Warning: One or both UPGMA trees were not created properly. Skipping subset.")
      next
    }
    
    # Compute Robinson-Foulds distance using phangorn (Fix for `@` error)
    rf_distance <- tryCatch({
      RF.dist(tree1, tree2)  # Corrected function
    }, error = function(e) {
      warning("Robinson-Foulds distance calculation failed: ", e$message)
      return(NA)
    })
    
    # Ensure nonzero values for probability match
    probability_match <- ifelse(is.na(rf_distance) | rf_distance == 0, 0.0001, rf_distance / n_splits)
    
    # Fix log(0) issue by setting min limit for probability match
    phylo_info <- ifelse(probability_match > 0, -log2(probability_match), 0)
    
    # Append results
    results <- rbind(results, data.frame(
      Matching_Trees = rf_distance,
      Probability_Match = probability_match,
      Phylo_Info_Content = phylo_info,
      Num_Sequences = n_splits
    ))
  }
  
  # Print results
  cat("\n_________________________________________________________________\n")
  cat("Matching Trees\tP(Match in random tree)\tPhylogenetic Info (bits)\tNumber of Sequences\n")
  print(results, row.names = FALSE)
  
  # Calculate and print averages
  avg_matching_trees <- mean(results$Matching_Trees, na.rm = TRUE)
  avg_probability_match <- mean(results$Probability_Match, na.rm = TRUE)
  avg_phylo_info <- mean(results$Phylo_Info_Content, na.rm = TRUE)
  
  cat("\n_________________________________________________________________\n")
  cat("Average Results:\n")
  cat("Matching trees:", round(avg_matching_trees, 2), "\n")
  cat("P(Match in random tree):", round(avg_probability_match, 4), "\n")
  cat("Phylogenetic information content:", round(avg_phylo_info, 2), "bits\n")
  cat("Total number of sequences:", taxon_number, "\n")
  cat("_________________________________________________________________\n")
  
  # Save results to file
  outfile <- paste0(sub("\\..*$", "", fas1), "_", sub("\\..*$", "", fas2), "_results.txt")
  filename <- file.path(output_directory, outfile)
  write.table(results, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Results saved to:", filename, "\n")
  
  # Decision based on probability match threshold
  if (!is.na(avg_probability_match) & avg_probability_match <= 0.05) {
    cat(fas1, " can be merged with ", fas2, "\n")
  } else {
    cat(fas1, " cannot be merged with ", fas2, "\n")
  }
}
