# Function to match files and identify unmatched records
match_files <- function(file1_path, file2_path, output_dir, file_type = "csv") {
  
  # Read files based on file type
  if (file_type == "csv") {
    df1 <- read.csv(file1_path, stringsAsFactors = FALSE)
    df2 <- read.csv(file2_path, stringsAsFactors = FALSE)
  } else if (file_type == "tsv") {
    df1 <- read.delim(file1_path, stringsAsFactors = FALSE)
    df2 <- read.delim(file2_path, stringsAsFactors = FALSE)
  } else if (file_type == "rds") {
    df1 <- readRDS(file1_path)
    df2 <- readRDS(file2_path)
  } else {
    stop("Unsupported file type. Use 'csv', 'tsv', or 'rds'.")
  }
  
  # Columns to match on
  match_columns <- c("gene_id", "gwas_indication", "gwas_index_snp", 
                      "eqtl_index_snp")
  
  # Check if all match columns exist in both datasets
  missing_cols_df1 <- match_columns[!match_columns %in% names(df1)]
  missing_cols_df2 <- match_columns[!match_columns %in% names(df2)]
  
  if (length(missing_cols_df1) > 0 || length(missing_cols_df2) > 0) {
    stop("Some matching columns are missing from the datasets.")
  }
  
  # Add source identifiers
  df1$source <- "pics"
  df2$source <- "susie"
  
  # Create a combined key for matching
  df1$match_key <- apply(df1[, match_columns], 1, paste, collapse = "_")
  df2$match_key <- apply(df2[, match_columns], 1, paste, collapse = "_")
  
  # Find matched and unmatched records
  matched_keys <- intersect(df1$match_key, df2$match_key)
  only_in_df1_keys <- setdiff(df1$match_key, df2$match_key)
  only_in_df2_keys <- setdiff(df2$match_key, df1$match_key)
  
  # Extract the records
  matched_df1 <- df1[df1$match_key %in% matched_keys, ]
  matched_df2 <- df2[df2$match_key %in% matched_keys, ]
  only_in_df1 <- df1[df1$match_key %in% only_in_df1_keys, ]
  only_in_df2 <- df2[df2$match_key %in% only_in_df2_keys, ]
  
  # Remove the temporary match_key column
  matched_df1$match_key <- NULL
  matched_df2$match_key <- NULL
  only_in_df1$match_key <- NULL
  only_in_df2$match_key <- NULL
  
  # Create a matched dataset (using full_join to keep all columns)
  matched_all <- full_join(matched_df1, matched_df2, by = match_columns, suffix = c("_file1", "_file2"))
  
  # Save the results
  if (file_type == "csv") {
    write.csv(matched_all, file.path(output_dir, "matched_records.csv"), row.names = FALSE)
    write.csv(only_in_df1, file.path(output_dir, "only_in_file1.csv"), row.names = FALSE)
    write.csv(only_in_df2, file.path(output_dir, "only_in_file2.csv"), row.names = FALSE)
  } else if (file_type == "tsv") {
    write.table(matched_all, file.path(output_dir, "matched_records.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(only_in_df1, file.path(output_dir, "only_in_file1.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(only_in_df2, file.path(output_dir, "only_in_file2.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  } else if (file_type == "rds") {
    saveRDS(matched_all, file.path(output_dir, "matched_records.rds"))
    saveRDS(only_in_df1, file.path(output_dir, "only_in_file1.rds"))
    saveRDS(only_in_df2, file.path(output_dir, "only_in_file2.rds"))
  }
  
  # Return a list with all datasets
  return(list(
    matched = matched_all,
    only_in_file1 = only_in_df1,
    only_in_file2 = only_in_df2
  ))
}

# Example usage
file1_path <- "/opt/notebooks/coloc_calc/coloc_results_20250528_140910.txt"   # Replace with actual file path
file2_path <- "/opt/notebooks/coloc_calc/coloc_results_20250528_141428.txt"  # Replace with actual file path
output_dir <- "/opt/notebooks"       # Replace with desired output directory

# Call the function
results <- match_files(file1_path, file2_path, output_dir, file_type = "tsv")

# Print a summary of the results
cat("Summary of Results:\n")
cat("- Matched records:", nrow(results$matched), "rows\n")
cat("- Only in file 1:", nrow(results$only_in_file1), "rows\n")
cat("- Only in file 2:", nrow(results$only_in_file2), "rows\n")

jpeg("myplot.jpg", width=800, height=600)
plot(dat1$PP_H4_file1,dat1$PP_H4_file2)
dev.off()
