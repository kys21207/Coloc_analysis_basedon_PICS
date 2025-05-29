#!/bin/bash

# Paths to required files and directories
GWAS_INPUT_DIR="/mnt/project/analyses_KJ/coloc_gwas_eqtls/pics_gwas_vs_eqtls/gwas_pics"
#EQTLS_INPUT_DIR="/mnt/project/analyses_KJ/coloc_gwas_eqtls/pics_gwas_vs_eqtls/eqtls_pics"
EQTLS_INPUT_DIR="/mnt/project/project-resources/ebi_eqtls_catalog/CS_sets_by_susie"
OUTPUT_DIR="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/coloc_gwas_eqtls/pics_gwas_vs_eqtls/coloc_gwas_eqtls_pics"

START_TIME=$(date +%Y-%m-%d\ %H:%M:%S)

echo "Starting PICS Coloc calculation at: $START_TIME"

# Loop through each chunk and run the PICS calculation
for EQTL_FILE in "$EQTLS_INPUT_DIR"/*.tsv.gz; do
#EQTL_FILE="/mnt/project/analyses_KJ/coloc_gwas_eqtls/pics_gwas_vs_eqtls/eqtls_pics/eqtl_pics_results_B_intermediate.txt"
    # Extract chunk identifier for naming output files
    EQTL_ID=$(basename "$EQTL_FILE" .tsv.gz)
    EQTL_NAME=$(basename "$EQTL_FILE" )
#    EQTL_ID=${EQTL_NAME#eqtl_pics_results_}  # Remove prefix
#    EQTL_ID=${EQTL_ID%.txt}                   # Remove suffix

    echo "Processing eQTL: $EQTL_ID"

    # Run Swiss Army Knife tool with PICS calculation
      dx run \
      --instance-type=mem3_ssd1_v2_x4 \
      --priority="normal" \
      --name "run_coloc_with_${EQTL_ID}" \
      --tag "coloc_calc_with_pics" \
      --brief \
      -y \
      swiss-army-knife \
      -iimage_file="project-Gz0KVk8JPx8q9yBzpgXX6y8K:project-resources/docker_images/tidyverse_docker_w_susie_coloc.tar.gz" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/analysis_scripts/coloc_pics_base/coloc_calc/main_coloc_pics_run.R" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/coloc_gwas_eqtls/pics_gwas_vs_eqtls/gwas_pics/gwas_pics_results_combined.txt" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:project-resources/ebi_eqtls_catalog/CS_sets_by_susie/${EQTL_NAME}" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:project-resources/ebi_eqtls_catalog/CS_sets_by_susie/eQTL-Catalogue-resources.txt" \
      --destination="${OUTPUT_DIR}" \
      -icmd="Rscript main_coloc_pics_run.R \
        --gwas_data_file gwas_pics_results_combined.txt \
        --eqtl_data_file ${EQTL_NAME} \
        --eqtl_info_file eQTL-Catalogue-resources.txt" 

    if [ $? -ne 0 ]; then
        echo "Error occurred while processing eQTL: $EQTL_ID"
        exit 1
    fi

done

END_TIME=$(date +%Y-%m-%d\ %H:%M:%S)
echo "PICS Coloc calculation completed at: $END_TIME"

# Calculate elapsed time
START_SECONDS=$(date -d "$START_TIME" +%s)
END_SECONDS=$(date -d "$END_TIME" +%s)
ELAPSED_SECONDS=$((END_SECONDS - START_SECONDS))
ELAPSED_MINUTES=$((ELAPSED_SECONDS / 60))

echo "Total computation time: $ELAPSED_MINUTES minutes ($ELAPSED_SECONDS seconds)"
echo "Results saved in: $OUTPUT_DIR"
