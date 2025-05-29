suppressPackageStartupMessages({
library(optparse)
library(tidyverse)
library(data.table)
library(stringr)
})

# Source the helper functions and parallel colocalization function
source("/mnt/project/analyses_KJ/analysis_scripts/coloc_pics_base/coloc_calc/utils_pics_coloc_calc.R")

# Parse command line arguments
option_list <- list(
  make_option("--gwas_data_file", type="character", default=NULL, 
              help="Path to the gwas data file (pics)", metavar="CHARACTER"),
  make_option("--eqtl_data_file", type="character", 
              default=NULL,
              help="Path to eqtl data file (pics)", metavar="CHARACTER"),
  make_option("--eqtl_info_file", type="character", 
              default=NULL,
              help="Path to eqtl info file", metavar="CHARACTER")
#  make_option("--output", type="character", default=NULL,
#              help="Output file name (if not specified, will be generated from phenotype name)", metavar="CHARACTER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

gwas_data_path     <- opt$gwas_data_file
eqtl_data_path     <- opt$eqtl_data_file
eqtl_info_path     <- opt$eqtl_info_file

#gwas_data_path<-"/mnt/project/analyses_KJ/coloc_gwas_eqtls/pics_gwas_vs_eqtls/gwas_pics/gwas_pics_results_combined.txt"
#eqtl_data_path<-"/mnt/project/project-resources/ebi_eqtls_catalog/CS_sets_by_susie/QTD000606.credible_sets.tsv.gz"
#eqtl_data_path<-"/mnt/project/analyses_KJ/coloc_gwas_eqtls/pics_gwas_vs_eqtls/eqtls_pics/eqtl_pics_results_B_intermediate.txt"
#eqtl_info_path<-"/mnt/project/project-resources/ebi_eqtls_catalog/CS_sets_by_susie/eQTL-Catalogue-resources.txt"

eqtl_info <- read.table(eqtl_info_path, header=T, sep="\t")

if(grepl("credible", basename(eqtl_data_path))){
  eqtl_prefix <- sub("\\.credible_sets\\.tsv\\.gz$", "", basename(eqtl_data_path))
  eqtl_info <- subset(eqtl_info, dataset_id == eqtl_prefix)
  eqtl_id <- eqtl_info$sample_group[1]
  eqtl_data <- convert_qtd_to_pics(eqtl_data_path)
}else if(grepl("pics", basename(eqtl_data_path))){
  eqtl_data <- fread(eqtl_data_path)
  eqtl_id <- sub(".*results_(.+)\\.txt$", "\\1", basename(eqtl_data_path))
  eqtl_data$gene_id <- str_extract(eqtl_data$indication, "ENSG\\d+")
}


# Read the data files
gwas_data <- fread(gwas_data_path)
gwas_id <- sub("\\_results\\_combined\\.txt$", "", basename(gwas_data_path))

# Find overlaps 
overlap_results <- find_overlaps_fast(gwas_data, eqtl_data)

# Run the improved parallel analysis that returns a data frame of successful results
coloc_results_df <- run_coloc_analysis_parallel(overlap_results, gwas_data, eqtl_data)

# Save the results to a file
if (nrow(coloc_results_df) > 0) {
  # Add timestamp to the output file
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_file <- paste0("coloc_pics_",gwas_id,"_vs_",eqtl_id,".txt")
  
  # Save the data frame to a file using write.table
  write.table(coloc_results_df, output_file, 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  cat(sprintf("Saved %d successful colocalization results to %s\n", nrow(coloc_results_df), output_file))
  
  # Sort by PP_H4 (posterior probability) to find the strongest colocalizations
#  coloc_results_df <- coloc_results_df[order(coloc_results_df$PP_H4, decreasing = TRUE), ]
  
  # Display the top results
#  cat("Top colocalization results:\n")
#  print(head(coloc_results_df, 10))
  
  # Additional analysis on successful results
#  significant_coloc <- coloc_results_df[coloc_results_df$PP_H4 > 0.5, ]
#  if (nrow(significant_coloc) > 0) {
#    cat(sprintf("Found %d significant colocalizations with PP.H4 > 0.5\n", nrow(significant_coloc)))
#  }
} else {
  cat("No colocalization results to save.\n")
}
