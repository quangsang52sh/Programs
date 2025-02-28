library(ape)          # For tree construction and distance matrix
library(phangorn)     # For UPGMA clustering
library(phylotaR)     # For Robinson-Foulds distance calculation
library(ggplot2)      # For visualization
library(Biostrings)   # For FASTA file parsing
library(seqinr)       # Alternative sequence handling
library(tidyverse)        # Data manipulation

ILD_analysis <- function(infile, output_directory = ".") {
  
  # Split input filenames
  files <- unlist(strsplit(infile, ","))
  fas1 <- files[1]
  fas2 <- files[2]
  
  # Read sequences
  sequences1 <- read.dna(fas1, format = "fasta")
  sequences2 <- read.dna(fas2, format = "fasta")
  
  # Number of taxa
  taxon_number <- nrow(sequences1)
  
  # Define splitting points
  s1 <- round(taxon_number * 1 / 3, 0)
  s2 <- round(taxon_number * 2 / 4, 0)
  s3 <- round(taxon_number * 4 / 4, 0)
  
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
    
    # Subset sequences
    splits1 <- sequences1[1:n_splits, ]
    splits2 <- sequences2[1:n_splits, ]
    
    # Compute pairwise distances
    distances1 <- dist.dna(splits1, model = "raw")
    distances2 <- dist.dna(splits2, model = "raw")
    
    # Build trees using UPGMA
    tree1 <- upgma(distances1)
    tree2 <- upgma(distances2)
    
    # Compute Robinson-Foulds distance
    rf_distance <- calcDstRF(tree1, tree2)
    
    # Compute probability and phylogenetic information content
    probability_match <- rf_distance / n_splits
    phylo_info <- -log2(probability_match)
    
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
  avg_matching_trees <- mean(results$Matching_Trees)
  avg_probability_match <- mean(results$Probability_Match)
  avg_phylo_info <- mean(results$Phylo_Info_Content)
  avg_seqs_info <- mean(results$Num_Sequences)
  
  cat("\n_________________________________________________________________\n")
  cat("Average Results:\n")
  cat("Matching trees:", round(avg_matching_trees, 2), "\n")
  cat("P(Match in random tree):", round(avg_probability_match, 4), "\n")
  cat("Phylogenetic information content:", round(avg_phylo_info, 2), "bits\n")
  cat("Total number of sequences:", taxon_number, "\n")
  cat("_________________________________________________________________\n")
  
  # Save results to file
  outfile <- paste0(sub("\\..*$", "", fas1), "_", sub("\\..*$", "", fas2), "_results.fas")
  filename <- file.path(output_directory, outfile)
  write.table(results, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Results saved to:", filename, "\n")
  
  # Decision based on probability match threshold
  if (avg_probability_match <= 0.05) {
    cat(fas1, " can be merged with ", fas2, "\n")
  } else {
    cat(fas1, " cannot be merged with ", fas2, "\n")
  }
}
