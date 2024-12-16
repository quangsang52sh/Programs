#!/bin/bash

# Variables
INPUT_FASTQ="shrimp_gut.fastq"
TRIMMED_FASTQ="shrimp_gut.trimmed.fastq"
ADAPTERS_FILE="adapters.fa"
LOG_FILE="trimming_log.txt"
FASTQC_DIR="fastqc_output"
METAPHLAN_OUT="shrimp_gut_metaphlan_output.txt"
BLASTN_DB="shrimp_gut_blastn_db"
BLASTN_RESULTS="shrimp_blastn_results.txt"
MICOM_DIR="shrimp_gut_micom_output"
AMRPLUSPLUS_DIR="shrimp_gut_amrplusplus_output"
LEFSE_OUT="shrimp_gut_lefse_output"
RDA_INPUT="rda_input.csv"
RDA_OUTPUT_DIR="rda_output"
HEATMAP_INPUT="heatmap_input.csv"
HEATMAP_OUTPUT="heatmap_output.png"
FINAL_REPORT="shrimp_gut_final_summary_report.txt"

# Step 1: Remove adapters with Trimmomatic
echo "Step 1: Trimming adapters from FASTQ file..."
trimmomatic SE -phred33 "$INPUT_FASTQ" "$TRIMMED_FASTQ" \
    ILLUMINACLIP:"$ADAPTERS_FILE":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> "$LOG_FILE"

# Step 2: Quality control using FastQC
echo "Step 2: Running FastQC for quality control..."
mkdir -p "$FASTQC_DIR"
fastqc "$TRIMMED_FASTQ" -o "$FASTQC_DIR"

# Step 3: Metagenomics analysis using Metaphlan
echo "Step 3: Analyzing metagenomics data with Metaphlan..."
metaphlan "$TRIMMED_FASTQ" --input_type fastq --output_file "$METAPHLAN_OUT"

# Step 4: BLAST analysis
echo "Step 4: Preparing reference database and running BLASTn..."
makeblastdb -in shrimp_reference_sequences.fasta -dbtype nucl -out "$BLASTN_DB"
blastn -query shrimp_query_sequence.fasta -db "$BLASTN_DB" \
    -out "$BLASTN_RESULTS" -outfmt "6 qseqid sseqid pident length evalue"

# Step 5: Microbiota diversity analysis with Micom
echo "Step 5: Analyzing microbiota diversity using Micom..."
mkdir -p "$MICOM_DIR"
micom community --table "$MICOM_DIR/otu_table.txt" \
    --metadata "$MICOM_DIR/metadata.txt" \
    --out_dir "$MICOM_DIR"

# Step 6: Antibiotic resistance analysis with AMRPlusPlus
echo "Step 6: Identifying antibiotic resistance genes with AMRPlusPlus..."
mkdir -p "$AMRPLUSPLUS_DIR"
amrplusplus -i "$TRIMMED_FASTQ" -o "$AMRPLUSPLUS_DIR"

# Step 7: LEfSe analysis
echo "Step 7: Running LEfSe analysis..."
run_lefse.py "$MICOM_DIR/otu_table.txt" "$LEFSE_OUT"

# Step 8: Redundancy Analysis (RDA)
echo "Step 8: Performing Redundancy Analysis (RDA)..."
mkdir -p "$RDA_OUTPUT_DIR"
Rscript -e "
library(vegan)
data <- read.csv('$RDA_INPUT', row.names = 1)
env <- data.frame(data[, c('Env1', 'Env2')])
matrix <- data.frame(data[, -c(1, 2)])
rda_result <- rda(matrix ~ ., data = env)
pdf('$RDA_OUTPUT_DIR/rda_plot.pdf')
plot(rda_result, main = 'RDA Analysis')
dev.off()
write.csv(summary(rda_result), file = '$RDA_OUTPUT_DIR/rda_summary.csv')
"

# Step 9: Heatmap generation
echo "Step 9: Creating a heatmap of microbial abundance..."
python3 -c "
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv('$HEATMAP_INPUT', index_col=0)
plt.figure(figsize=(10, 8))
sns.heatmap(data, cmap='viridis', annot=False)
plt.title('Heatmap of Microbial Abundance')
plt.savefig('$HEATMAP_OUTPUT')
"

# Step 10: Generate final summary report
echo "Step 10: Compiling final summary report..."
{
    echo "Date: $(date)"
    echo "Input FASTQ: $INPUT_FASTQ"
    echo "Trimmed FASTQ: $TRIMMED_FASTQ"
    echo "Adapters File: $ADAPTERS_FILE"
    echo "Metaphlan Output: $METAPHLAN_OUT"
    echo "BLASTn Results: $BLASTN_RESULTS"
    echo "Micom Output Directory: $MICOM_DIR"
    echo "AMRPlusPlus Output Directory: $AMRPLUSPLUS_DIR"
    echo "LEfSe Output: $LEFSE_OUT"
    echo "RDA Output Directory: $RDA_OUTPUT_DIR"
    echo "Heatmap Output: $HEATMAP_OUTPUT"
} > "$FINAL_REPORT"

echo "Pipeline completed successfully. Summary report saved to $FINAL_REPORT."
