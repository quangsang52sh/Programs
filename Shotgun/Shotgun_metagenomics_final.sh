#!/bin/bash

# Variables
input_fastq="shrimp_gut.fastq"
output_trimmed="shrimp_gut.trimmed.fastq"
adapters_file="adapters.fa"
log_file="trimming_log.txt"
fastqc_output="fastqc_output"
metaphlan_output="shrimp_gut_metaphlan_output"
blastn_db="shrimp_gut_blastn_db"
copepod_output="shrimp_gut_copepod_output"
amrplusplus_output="shrimp_gut_amrplusplus_output"
lefse_output="shrimp_gut_lefse_output"
rda_input="rda_input.csv"
rda_output="rda_output"
heatmap_input="heatmap_input.csv"
heatmap_output="heatmap_output.png"

# Step 1: Remove adapters with Trimmomatic
echo "Step 1: Removing adapters with Trimmomatic..."
trimmomatic SE -phred33 "$input_fastq" "$output_trimmed" \
    ILLUMINACLIP:"$adapters_file":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> "$log_file"

# Step 2: Quality control
echo "Step 2: Running FastQC for quality control..."
fastqc "$output_trimmed" -o "$fastqc_output"

# Step 3: Analyze metagenomics with Metaphlan 3.0
echo "Step 3: Analyzing metagenomics with Metaphlan 3.0..."
metaphlan "$output_trimmed" --input_type fastq \
    --output_file "$metaphlan_output"

# Step 4: Prepare reference database and use blastn
echo "Step 4: Preparing reference database and using blastn..."
makeblastdb -in shrimp_reference_sequences.fasta -dbtype nucl -out "$blastn_db"
blastn -query shrimp_query_sequence.fasta -db "$blastn_db" \
    -out shrimp_blastn_results.txt -outfmt "6 qseqid sseqid pident length evalue"

# Step 5: Analyze diversity of gut shrimp microbiota with micom
echo "Step 5: Analyzing diversity of gut shrimp microbiota with micom..."
micom community --table $copepod_output/otu_table.txt \
    --metadata $copepod_output/metadata.txt \
    --out_dir $copepod_output/micom_output

# Step 6: Analyze antibiotic resistance genes with AMRPlusPlus
echo "Step 6: Analyzing antibiotic resistance genes with AMRPlusPlus..."
amrplusplus -i "$output_trimmed" -o "$amrplusplus_output"

# Step 7: Perform LEfSe analysis
echo "Step 7: Performing LEfSe analysis..."
run_lefse.py "$copepod_output/otu_table.txt" "$lefse_output"

# Step 8: Perform Redundancy Analysis (RDA)
echo "Step 8: Performing Redundancy Analysis (RDA)..."
Rscript -e "
library(vegan)
data <- read.csv('$rda_input', row.names = 1)
env <- data.frame(data[, c('Env1', 'Env2')])
matrix <- data.frame(data[, -c(1, 2)])
rda_result <- rda(matrix ~ ., data = env)
pdf('$rda_output/rda_plot.pdf')
plot(rda_result, main = 'RDA Analysis')
dev.off()
write.csv(summary(rda_result), file = '$rda_output/rda_summary.csv')
"

# Step 9: Generate heatmap
echo "Step 9: Generating heatmap..."
python3 -c "
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv('$heatmap_input', index_col=0)
plt.figure(figsize=(10, 8))
sns.heatmap(data, cmap='viridis', annot=False)
plt.title('Heatmap of Microbial Abundance')
plt.savefig('$heatmap_output')
"

# Step 10: Generate a final summary report
echo "Step 10: Generating a final summary report..."
echo "Date: $(date)" > shrimp_gut_final_summary_report.txt
echo "Input FASTQ: $input_fastq" >> shrimp_gut_final_summary_report.txt
echo "Trimmed Output: $output_trimmed" >> shrimp_gut_final_summary_report.txt
echo "Adapters File: $adapters_file" >> shrimp_gut_final_summary_report.txt
echo "Metaphlan Output: $metaphlan_output" >> shrimp_gut_final_summary_report.txt
echo "Blastn Database: $blastn_db" >> shrimp_gut_final_summary_report.txt
echo "Shrimp Gut Micom Output: $copepod_output/micom_output" >> shrimp_gut_final_summary_report.txt
echo "AMRPlusPlus Output: $amrplusplus_output" >> shrimp_gut_final_summary_report.txt
echo "LEfSe Output: $lefse_output" >> shrimp_gut_final_summary_report.txt
echo "RDA Output: $rda_output" >> shrimp_gut_final_summary_report.txt
echo "Heatmap Output: $heatmap_output" >> shrimp_gut_final_summary_report.txt
