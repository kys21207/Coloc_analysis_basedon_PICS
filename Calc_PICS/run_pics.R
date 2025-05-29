
suppressPackageStartupMessages({
library(optparse)
library(tidyverse)
})

# Load the PICS function
source("/mnt/project/analyses_KJ/analysis_scripts/coloc_pics_base/pics_calc/pics_calc_EUR_custom.R")

# Path to files and directories
#ld_block_file <- "./ld_blocks_with_ids.bed"  # Your LD block definitions
ld_dir <- "/mnt/project/analyses_KJ/ld_matrices"          # Directory containing *.ld.gz and *.vars files
#index_data_file <- "./hyprcoloc_lead_pair_for_coloc.tsv"
#index_data_file <- "./test_example100.tsv"
#output_dir <- "/opt/notebooks/pics"

# Parse command line arguments
option_list <- list(
  make_option("--ld_block_file", type="character", default=NULL, 
              help="Path to the ld block info file", metavar="CHARACTER"),
  make_option("--index_data_file", type="character", 
              default=NULL,
              help="Path to eqtl index file", metavar="CHARACTER"),
  make_option("--output", type="character", default=NULL,
              help="Output file name (if not specified, will be generated from phenotype name)", metavar="CHARACTER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

index_data_path     <- opt$index_data_file
ld_block_path       <- opt$ld_block_file
output_file         <- opt$output

#index_data_path <-"/opt/notebooks/test_data.txt"
#ld_block_path <- "/opt/notebooks/pics_calc/ld_blocks_with_ids.bed"
# Run PICS calculation with memory optimization parameters
pics_results <- pics_calc(
  index_data_file = index_data_path,
  ld_block_file = ld_block_path,
  ld_dir = ld_dir,
  nCORE_pics = 1           # Single-threaded to prevent memory issues
       # Switch to row-by-row processing if matrices exceed 1GB
)

# Save the results to a file
# Use fwrite for memory efficiency with large result sets
data.table::fwrite(pics_results, output_file, sep = "\t")
cat(paste("Results saved to:", output_file, "\n"))

