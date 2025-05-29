#!/bin/bash

# Memory-optimized shell script to run PICS calculation for chunks of index_data_file

# Paths to required files and directories
#LD_BLOCK_FILE="/mnt/project/analyses_KJ/analysis_codes/coloc_pics_base/ld_blocks_with_ids.bed" # LD block definitions
#LD_DIR="/mnt/project/analyses_KJ/ld_matrices" # Directory containing LD matrix files
#INDEX_DATA_FILE="/mnt/project/analyses_KJ/analysis_codes/coloc_pics_base/hyprcoloc_lead_pairs_b38.tsv" # Input file with index data
INPUT_DIR="/mnt/project/analyses_KJ/analysis_scripts/coloc_pics_base/gwas_input"
OUTPUT_DIR="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/coloc_gwas_eqtls/pics_gwas_vs_eqtls/gwas_pics"

# Determine the total number of lines in the index data file
#TOTAL_LINES=$(wc -l < "$INDEX_DATA_FILE")
START_TIME=$(date +%Y-%m-%d\ %H:%M:%S)

echo "Starting PICS calculation at: $START_TIME"

# Loop through each chunk and run the PICS calculation
#for CHUNK_FILE in "$INPUT_DIR"/*.tsv; do
CHUNK_FILE="/mnt/project/analyses_KJ/analysis_scripts/coloc_pics_base/gwas_input/chunk_1.tsv"
    # Extract chunk identifier for naming output files
#    CHUNK_ID=$(basename "$CHUNK_FILE" _top_hits_susie_results_overlap_topmed_sites.tsv)
    CHUNK_ID=$(basename "$CHUNK_FILE" .tsv)
    CHUNK_NAME=$(basename "$CHUNK_FILE" )

    echo "Processing chunk: $CHUNK_ID"

    # Run Swiss Army Knife tool with PICS calculation
      dx run \
      --instance-type=mem3_ssd1_v2_x8 \
      --priority="normal" \
      --name "run_pics_${CHUNK_ID}" \
      --tag "pics_calc_gwas" \
      --brief \
      -y \
      swiss-army-knife \
      -iimage_file="project-Gz0KVk8JPx8q9yBzpgXX6y8K:/project-resources/docker_images/tidyverse_docker_w_susie_coloc.tar.gz" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/analysis_scripts/coloc_pics_base/pics_calc/run_pics.R" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/analysis_scripts/coloc_pics_base/pics_calc/ld_blocks_with_ids.bed" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/analysis_scripts/coloc_pics_base/gwas_input/${CHUNK_NAME}" \
      --destination="${OUTPUT_DIR}" \
      -icmd="Rscript run_pics.R \
        --ld_block_file ld_blocks_with_ids.bed \
        --index_data_file  ${CHUNK_NAME} \
        --output gwas_pics_results_${CHUNK_ID}.txt"

    if [ $? -ne 0 ]; then
        echo "Error occurred while processing chunk: $CHUNK_ID"
        exit 1
    fi

    # Optional: Clean up chunk file after processing (uncomment if needed)
    # rm "$CHUNK_FILE"
#done

END_TIME=$(date +%Y-%m-%d\ %H:%M:%S)
echo "PICS calculation completed at: $END_TIME"

# Calculate elapsed time
START_SECONDS=$(date -d "$START_TIME" +%s)
END_SECONDS=$(date -d "$END_TIME" +%s)
ELAPSED_SECONDS=$((END_SECONDS - START_SECONDS))
ELAPSED_MINUTES=$((ELAPSED_SECONDS / 60))

echo "Total computation time: $ELAPSED_MINUTES minutes ($ELAPSED_SECONDS seconds)"
echo "Results saved in: $OUTPUT_DIR"
