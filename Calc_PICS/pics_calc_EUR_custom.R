#==========================================================================================================================
#' @title Cross-platform memory-optimized PICS calculation using existing LD matrices
#' 
#' @description The PICS algorithm calculates the most likely causal SNPs given the observed association signal
#'   at a locus, optimized for memory efficiency with large LD matrices.
#'   For an associated locus, enter the most highly-associated SNP 
#'   (referred to as the index SNP) and the strength of association. Using pre-computed LD matrices
#'   for EUR population, the algorithm identifies the SNPs that are most likely to be the 
#'   causal variants responsible for the association (PICS_Probability).
#'
#'   See \strong{Genetic and Epigenetic Fine-Mapping of Causal Variants in Autoimmune Disease} 
#'   by Kyle Kai-How Farh,et al. Nature 518, 337–343 (19 February 2015)
#' 
#' @param index_data_file Path to file containing index SNPs data with columns: snpID, pval, indication
#'   - snpID: SNP identifier in format chr:pos:ref:alt (e.g., 1:123456:A:C)
#'   - pval: p-value of association for the SNP
#'   - indication: disease/trait being studied
#' @param ld_block_file Path to the LD block definition file (BED format)
#' @param ld_dir Directory containing LD matrices (*.ld.gz) and variant info (*.vars)
#' @param nCORE_pics Define the number of cores for multiple tasks
#'   
#' @return A data frame containing the result of the PICS calculation
#' 
#' @author Kijoung Song (Modified by: ksong0924)
#' @date 2025-05-16
#' 
#' @export
                        
pics_calc_old <- function(index_data_file, ld_block_file, ld_dir, nCORE_pics = 1,
                      chunk_size = 1000, max_matrix_gb = 2) {
  # Load required packages
  if(!require(dplyr)) install.packages("dplyr")
  if(!require(data.table)) install.packages("data.table")
  
  library(dplyr)
  library(data.table)
  
  # Cross-platform memory management function
  clean_memory <- function() {
    # Just call garbage collector - works on all platforms
    gc(verbose = FALSE)
    message("Memory cleaned")
  }
  
  # Start with clean memory
  clean_memory()
  
  # Process index SNPs
  message("Processing index SNPs...")
  
  # Read index data with data.table for memory efficiency
  index.data <- fread(index_data_file, showProgress = FALSE)
  
  # Verify that required columns exist
  required_cols <- c("snpID", "pval", "indication")
  missing_cols <- required_cols[!required_cols %in% names(index.data)]
  
  if(length(missing_cols) > 0) {
    stop("Missing required columns in index data: ", paste(missing_cols, collapse=", "))
  }
  
  # Process SNPs in chunks if the dataset is large
  message(paste("Processing", nrow(index.data), "index SNPs"))
  
  # Parse SNP IDs to extract chromosome and position more efficiently
  # Using data.table syntax for faster performance
  index.data[, `:=`(
    raw_chrom = sapply(strsplit(as.character(snpID), ":"), function(x) x[1]),
    pos = as.numeric(sapply(strsplit(as.character(snpID), ":"), function(x) x[2])),
    ref = sapply(strsplit(as.character(snpID), ":"), function(x) if(length(x) > 2) x[3] else NA),
    alt = sapply(strsplit(as.character(snpID), ":"), function(x) if(length(x) > 3) x[4] else NA)
  )]
  
  index.data[, `:=`(
    chrom = raw_chrom,
    chr_for_block = paste0("chr", raw_chrom),
    ancestry = "EUR",
    rsid1 = snpID
  )]
  
  # Clean memory after initial processing
  clean_memory()
  
  message("Reading LD block file...")
  ld_blocks <- fread(ld_block_file, header = FALSE, showProgress = FALSE)
  setnames(ld_blocks, c("chrom", "start", "end", "block_id"))
  
  # Create index for faster lookups
  setkey(ld_blocks, chrom)
  
  #==========================================================================
  # Find which LD block contains each index SNP - simplified approach
  message("Determining LD blocks for index SNPs...")
  index.data[, block_id := NA_character_]
  
  # Process in batches to avoid memory issues with large datasets
  batch_size <- min(chunk_size, nrow(index.data))
  num_batches <- ceiling(nrow(index.data) / batch_size)
  
  for(batch_idx in 1:num_batches) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, nrow(index.data))
    
    # Process this batch
    for(i in start_idx:end_idx) {
      chr <- index.data$chr_for_block[i]
      pos <- index.data$pos[i]
      
      # Get potential blocks for this chromosome
      chr_blocks <- ld_blocks[chrom == chr]
      if(nrow(chr_blocks) > 0) {
        # Find block containing this position
        block <- chr_blocks[start <= pos & end >= pos]
        if(nrow(block) > 0) {
          index.data$block_id[i] <- block$block_id[1]
        }
      }
    }
    
    # Clean memory after each large batch
    if(batch_idx %% 10 == 0) {
      message(paste("Processed", end_idx, "of", nrow(index.data), "SNPs"))
      clean_memory()
    }
  }
  
  # Remove SNPs that don't fall within an LD block
  index.data <- index.data[!is.na(block_id)]
  
  if(nrow(index.data) == 0) {
    stop("None of the index SNPs fall within defined LD blocks.")
  }
  
  # Exclude chromosomes "X" and "XY"
  index.data <- index.data[!(chrom %in% c("X", "Y", "XY"))]
  clean_memory()
  
  #==========================================================================
  # Optimization: Group index SNPs by LD block
  message("Grouping index SNPs by LD block for efficient processing...")
  index.data <- as.data.frame(index.data)  # Convert back for compatibility
  block_groups <- split(index.data, index.data$block_id)
  
  # Ensure each block_group has all required columns
  required_group_cols <- c("snpID", "pval", "indication", "chrom", "pos")
  
  # Free original data to save memory
  rm(index.data)
  clean_memory()
  
  # Process each block group sequentially to avoid memory issues
  all_results <- list()
  block_counter <- 0
  total_blocks <- length(block_groups)
  
  for(block_name in names(block_groups)) {
    block_counter <- block_counter + 1
    block_group <- block_groups[[block_name]]
    block_id <- unique(block_group$block_id)
    message(paste("Processing LD block", block_counter, "of", total_blocks, ":", block_id))
    
    # Verify required columns exist in this block_group
    missing_cols <- required_group_cols[!required_group_cols %in% names(block_group)]
    if(length(missing_cols) > 0) {
      message(paste("Missing columns in block", block_id, ":", paste(missing_cols, collapse=", ")))
      message("Skipping this block due to missing columns")
      next
    }
    
   on.exit({
     if (exists("block_group", inherits = FALSE)) {
       rm(block_group)
     }
     clean_memory()
   }, add = TRUE)
    
    tryCatch({
      # Construct filenames for LD matrix and variants
      ld_file <- file.path(ld_dir, paste0("ld_matrix_block_", block_id, ".ld.gz"))
      vars_file <- file.path(ld_dir, paste0("ld_matrix_block_", block_id, ".vars"))
      
      # Skip if files don't exist
      if(!file.exists(ld_file) || !file.exists(vars_file)) {
        message("Missing LD files for block ", block_id)
        next
      }
      
      # Read variant information for this block (no header)
      message("Reading variant information...")
      vars_raw <- readLines(vars_file)
      
      if(length(vars_raw) == 0) {
        message("Empty variants file for block ", block_id)
        next
      }
      
      # Check the potential matrix size
      matrix_size_gb <- length(vars_raw)^2 * 8 / (1024^3)  # Size in GB
      if(matrix_size_gb > max_matrix_gb) {
        message(paste("Warning: LD matrix for block", block_id, "estimated size:", 
                     round(matrix_size_gb, 2), "GB. Using chunk processing."))
      }
      
      # Process variants efficiently
      variant_parts <- strsplit(vars_raw, ":")
      vars_chrom <- sapply(variant_parts, function(x) x[1])
      vars_position <- as.numeric(sapply(variant_parts, function(x) x[2]))
      
      # Create a lookup table for positions to row indices
      pos_to_index <- data.frame(
        position = vars_position,
        row_index = 1:length(vars_position),
        stringsAsFactors = FALSE
      )
      
      block_results <- NULL
      
      # If the matrix is too large, process index SNPs one by one
      if(matrix_size_gb > max_matrix_gb) {
        message("Processing large LD matrix by index SNP...")
        
        # Process one index SNP at a time
        for(i in 1:nrow(block_group)) {
          index_pos <- block_group$pos[i]
          index_rows <- which(vars_position == index_pos)
          
          if(length(index_rows) == 0) {
            message(paste("Position not found in vars file:", index_pos))
            next
          }
          
          # If multiple variants at the same position, use the first one
          if(length(index_rows) > 1) {
            message(paste("Multiple variants at position", index_pos, "- using the first one"))
            index_row <- index_rows[1]
          } else {
            index_row <- index_rows
          }
          
          # Read only the specific row from the LD matrix
          cmd <- paste0("zcat ", shQuote(ld_file), " | sed -n '", index_row, "p'")
          ld_values <- NULL
          
          tryCatch({
            # Try direct reading with system command - faster on Unix systems
            ld_line <- system(cmd, intern = TRUE)
            if(length(ld_line) > 0) {
              ld_values <- as.numeric(strsplit(ld_line, ",")[[1]])
            }
          }, error = function(e) {
            message("System command failed, falling back to full matrix read")
          })
          
          # Fallback if system command fails
          if(is.null(ld_values)) {
            # Read full matrix but only keep the row we need
            message("Reading full LD matrix for row extraction...")
            ld_matrix <- as.matrix(fread(ld_file))
            ld_values <- ld_matrix[index_row, ]
            rm(ld_matrix)
            clean_memory()
          }
          
          # Replace NaN values with 0
          ld_values[is.na(ld_values)] <- 0
          
          # Create data frame with LD information - only keep SNPs with r2 > 0.5 to save memory
          keep_indices <- which(ld_values > 0.5)
          
          if(length(keep_indices) > 0) {
            snp_ld <- data.frame(
              chrom1 = block_group$chrom[i],
              pos1 = index_pos,
              chrom2 = vars_chrom[keep_indices],
              pos2 = vars_position[keep_indices],
              r2 = ld_values[keep_indices],
              snpID = vars_raw[keep_indices],  # Changed from variant_id to snpID for consistency
              index_snp = block_group$snpID[i], # Changed from index_snp_id to index_snp for final output
              pval = block_group$pval[i],
              indication = block_group$indication[i],
              ancestry = "EUR",
              stringsAsFactors = FALSE
            )
            
            # Calculate PICS scores directly for this SNP
            if(nrow(snp_ld) > 0) {
              snp_ld$pval <- ifelse(snp_ld$pval == 0, 1e-250, snp_ld$pval)
              snp_ld$Mean <- snp_ld$r2 * (-log10(snp_ld$pval))
              snp_ld$SD <- sqrt(1 - snp_ld$r2^3.2) * sqrt(-log10(snp_ld$pval))/2
              snp_ld$SD <- ifelse(snp_ld$SD < 0.057, 0.057, snp_ld$SD)
              snp_ld$prob <- dnorm(-log10(snp_ld$pval), snp_ld$Mean, snp_ld$SD)
              
              # Calculate probability sum and PICS score
              prob_sum <- sum(snp_ld$prob)
              snp_ld$pics <- snp_ld$prob / prob_sum
              
              # Add to results
              block_results <- rbind(block_results, snp_ld)
            }
          }
          
          # Clean memory after each SNP to prevent accumulation
          rm(ld_values, keep_indices)
          if(i %% 5 == 0) clean_memory()
        }
        
      } else {
        # Standard processing for smaller matrices
        message("Reading full LD matrix...")
        ld_matrix <- as.matrix(fread(ld_file))
        
        # Replace any NaN values with 0 in the entire matrix
        if(any(is.na(ld_matrix))) {
          ld_matrix[is.na(ld_matrix)] <- 0
        }
        
        # Process each index SNP within this block
        for(i in 1:nrow(block_group)) {
          index_pos <- block_group$pos[i]
          index_row <- which(vars_position == index_pos)
          
          if(length(index_row) == 0) {
            message(paste("Position not found in vars file:", index_pos))
            next
          }
          
          # If multiple variants at the same position, use the first one
          if(length(index_row) > 1) {
            message(paste("Multiple variants at position", index_pos, "- using the first one"))
            index_row <- index_row[1]
          }
          
          # Extract LD values for this index SNP
          ld_values <- ld_matrix[index_row, ]
          
          # Only keep SNPs with r2 > 0.5 to save memory
          keep_indices <- which(ld_values > 0.5)
          
          if(length(keep_indices) > 0) {
            # Create data frame with LD information - only for SNPs with r2 > 0.5
            snp_ld <- data.frame(
              chrom1 = block_group$chrom[i],
              pos1 = index_pos,
              chrom2 = vars_chrom[keep_indices],
              pos2 = vars_position[keep_indices],
              r2 = ld_values[keep_indices],
              snpID = vars_raw[keep_indices],  # Changed from variant_id to snpID for consistency
              index_snp = block_group$snpID[i], # Changed from index_snp_id to index_snp for final output
              pval = block_group$pval[i],
              indication = block_group$indication[i],
              ancestry = "EUR",
              stringsAsFactors = FALSE
            )
            
            # Calculate PICS directly
            if(nrow(snp_ld) > 0) {
              snp_ld$pval <- ifelse(snp_ld$pval == 0, 1e-250, snp_ld$pval)
              snp_ld$Mean <- snp_ld$r2 * (-log10(snp_ld$pval))
              snp_ld$SD <- sqrt(1 - snp_ld$r2^3.2) * sqrt(-log10(snp_ld$pval))/2
              snp_ld$SD <- ifelse(snp_ld$SD < 0.057, 0.057, snp_ld$SD)
              snp_ld$prob <- dnorm(-log10(snp_ld$pval), snp_ld$Mean, snp_ld$SD)
              
              # Calculate probability sum and PICS score
              prob_sum <- sum(snp_ld$prob)
              snp_ld$pics <- snp_ld$prob / prob_sum
              
              # Add to results
              block_results <- rbind(block_results, snp_ld)
            }
          }
        }
        
        # Free memory
        rm(ld_matrix)
        clean_memory()
      }
      
      # If we have results for this block, add them to our list
      if(!is.null(block_results) && nrow(block_results) > 0) {
        # Select only the required columns for the final output
        # Make sure to use column names that actually exist in block_results
        result_cols <- c("snpID", "pics", "ancestry", "indication", "index_snp", "r2")
        
        # Check which columns exist in our results
        available_cols <- result_cols[result_cols %in% names(block_results)]
        
        if(length(available_cols) < length(result_cols)) {
          missing <- result_cols[!result_cols %in% available_cols]
          message(paste("Warning: Missing columns in block", block_id, "results:", 
                        paste(missing, collapse=", ")))
        }
        
        # Only select columns that actually exist
        if(length(available_cols) > 0) {
          all_results[[length(all_results) + 1]] <- block_results[, available_cols, drop=FALSE]
        } else {
          message("No valid columns found in block results, skipping")
        }
      }
      
      # Clean up
      rm(block_results, vars_raw, variant_parts, vars_chrom, vars_position)
      clean_memory()
      
    }, error = function(e) {
      message(paste("Error processing block", block_id, ":", e$message))
    })
    
    # Remove this block group to save memory
    block_groups[[block_name]] <- NULL
    
    # Clean memory every few blocks
    if(block_counter %% 5 == 0 || block_counter == total_blocks) {
      clean_memory()
    }
  }
  
  # Clean up remaining objects
  rm(block_groups)
  clean_memory()
  
  # Check if we have any results at all
  if(length(all_results) == 0) {
    stop("No SNPs in LD were found. Please check your input data and LD files.")
  }
  
  # Combine all results efficiently
  message("Combining results from all blocks...")
  
  # Process in smaller batches to avoid memory issues
  batch_size <- 10
  final_results <- NULL
  
  for(i in seq(1, length(all_results), by = batch_size)) {
    end_idx <- min(i + batch_size - 1, length(all_results))
    
    # Only combine if we have results in this batch
    if(length(all_results[i:end_idx]) > 0 && !all(sapply(all_results[i:end_idx], is.null))) {
      # Filter out NULL entries
      valid_results <- all_results[i:end_idx][!sapply(all_results[i:end_idx], is.null)]
      
      if(length(valid_results) > 0) {
        # Check if all data frames have the same columns
        all_cols <- unique(unlist(lapply(valid_results, colnames)))
        
        # Create a batch result with all columns
        batch_results <- NULL
        
        # Process each result to ensure column consistency
        for(j in 1:length(valid_results)) {
          result_j <- valid_results[[j]]
          
          # Make sure all columns exist
          for(col in all_cols) {
            if(!(col %in% names(result_j))) {
              result_j[[col]] <- NA
            }
          }
          
          # Add to batch results
          if(is.null(batch_results)) {
            batch_results <- result_j
          } else {
            batch_results <- rbind(batch_results, result_j)
          }
        }
        
        # Add to final results
        if(is.null(final_results)) {
          final_results <- batch_results
        } else {
          # Ensure columns match between final_results and batch_results
          for(col in setdiff(names(batch_results), names(final_results))) {
            final_results[[col]] <- NA
          }
          for(col in setdiff(names(final_results), names(batch_results))) {
            batch_results[[col]] <- NA
          }
          
          final_results <- rbind(final_results, batch_results)
        }
      }
    }
    
    # Clean up batch
    all_results[i:end_idx] <- NULL
    clean_memory()
  }
  
  # Final sort and clean
  if(!is.null(final_results) && nrow(final_results) > 0) {
    # Sort by indication and decreasing PICS score
    # Only sort if these columns exist
    if(all(c("index_snp", "indication", "pics") %in% names(final_results))) {
      final_results <- final_results[order(final_results$index_snp, final_results$indication, -final_results$pics), ]
    }
    
    # Make sure the requested columns are present (with NAs if necessary)
    expected_cols <- c("snpID", "pics", "ancestry", "indication", "index_snp", "r2")
    for(col in expected_cols) {
      if(!(col %in% names(final_results))) {
        final_results[[col]] <- NA
        message(paste("Warning: Creating empty column", col, "in final results"))
      }
    }
    
    # Select only the expected columns in the requested order
    final_results <- final_results[, expected_cols]
  }
  
  message("PICS calculation completed.")
  return(final_results)
}
                                  
                                  
#=============================================================================================================================  
#' int_coloc_pics_lite 
#'   Test for colocalization of two PICS sets
#' 
#' @return only H3 & H4 posteriors
#' @param data1,data2  PICS sets from read.pics or download.pics
#' @param pics1,pics2  column name to pull PICS prob from. Default = "PICS_probability"
#' @param rsid1,rsid2  column name to pull rsid from.      Default = "Linked_SNP"
#' @param rounded   Decimal points to round posteriors to
#' @param priorc1   Prior probability for colocalization with siganl for data1  Default = 1e-4
#' @param priorc2   Prior probability for colocalization with siganl for data2 Default = 1e-4
#' @param priorc12  Prior probability for colocalization of both signals.   Default = 1e-5
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}  


int_coloc_pics_lite <- function(data1,
                                data2,
                                pics1    = "PICS_probability", # column header for poster probabilities in data1
                                pics2    = "PICS_probability", # column header for poster probabilities in data2
                                rsid1    = "Linked_SNP",       # column header for snps in LD in data1
                                rsid2    = "Linked_SNP",       # column header for snps in LD in data2
                                rounded  = 6,
                                priorc1  = 1e-4, 
                                priorc2  = 1e-4, 
                                priorc12 = 1e-5
) {
  stopifnot(exists("data1") & exists("data2"))
  if(is.logical(data1)){
    if(is.na(data1)){
      gtx_warn("int_coloc_pics_lite: data1 is NA, skipping coloc.")
      return(list(results = NA, nvariants = NA))
    }
  }
  if(is.logical(data2)){
    if(is.na(data2)){
      gtx_warn("int_coloc_pics_lite: data2 is NA, skipping coloc.")
      return(list(results = NA, nvariants = NA))
    }
  }
  pics <- int_harmonize_pics(data1, 
                             data2, 
                             opts <- data.frame(rsid1 = rsid1, rsid2 = rsid2, pics1 = pics1, pics2 = pics2, stringsAsFactors = FALSE))
  
  nv <- dim(pics)[1]
  res <- data.frame(prior = norm(c(priorc1*priorc2*nv*(nv - 1), priorc12*nv)),
                    bf    = c((sum(pics[[1]])*sum(pics[[2]]) - sum(pics[[1]]*pics[[2]]))/(nv*(nv - 1)), 
                              sum(pics[[1]]*pics[[2]])/nv))
  res$bf <- res$bf/res$bf[1]
  res$posterior <- norm(res$prior*res$bf)
  if (is.finite(rounded)) {
    res$posterior = round(res$posterior, rounded)
  }
  return(res$posterior[2])
}

#' int_pico_lite 
#'     Simplify int_coloc_pics_lite function so that we can use map2
int_pico_lite <- function(x,y) int_coloc_pics_lite1(x,y,pics1="pics",pics2="pics",rsid1="pos",rsid2="pos")
                     
  
#' @title int_harmonize_pics
#' @description  int_harmonize_pics
#' @param data1 TBD
#' @param data2 TBD
#' @param opts TBD
int_harmonize_pics <- function(data1,
                               data2, 
                               opts = data.frame(pics1 = "PICS_probability",
                                                 pics2 = "PICS_probability",
                                                 rsid1 = "Linked_SNP",
                                                 rsid2 = "Linked_SNP",
                                                 stringsAsFactors = FALSE)
){
  ids <- unique(c(data1[[opts$rsid1]], data2[[opts$rsid2]]))
  tmp <- as.data.frame(matrix(data = NA, nrow = length(ids), ncol = 2))
  pp1 <- if (opts$pics1==opts$pics2) paste(opts$pics1, ".1", sep = "") else opts$pics1
  pp2 <- if (opts$pics1==opts$pics2) paste(opts$pics2, ".2", sep = "") else opts$pics2
  colnames(tmp) <- c(pp1, pp2)
  for(n in 1:length(ids)){
    tmp[[pp1]][n] <- if(!is.na(match(ids[n], data1[[opts$rsid1]]))) data1[which(data1[[opts$rsid1]]==ids[n]),][[opts$pics1]][1] else 0
    tmp[[pp2]][n] <- if(!is.na(match(ids[n], data2[[opts$rsid2]]))) data2[which(data2[[opts$rsid2]]==ids[n]),][[opts$pics2]][1] else 0 
  }
  res <- as.data.frame(cbind(
    norm1(tmp[[pp1]], log = FALSE),
    norm1(tmp[[pp2]], log = FALSE)
  ))
  colnames(res) <- c(pp1, pp2)
  rownames(res) <- ids
  return(res)
}



#' Combine Shape Objects
#'
#'   Combine shape objects into one shape object. It works analogous to rbind.
#'
#' @author  Martijn Tennekes
#'

#' @title int_sbind
#' @description int_sbind
#' @param x TBD
#' @param y TBD
#' @param fill TBD
int_sbind = function(x, y, fill=NA) {
  int_sbind.fill = function(d, cols){ 
    for(c in cols)
      d[[c]] = fill
    d
  }
  
  x = int_sbind.fill(x, setdiff(names(y),names(x)))
  y = int_sbind.fill(y, setdiff(names(x),names(y)))
  
  rbind(x, y)
}

#' An Alternative To The Internal Do.Call
#'  
#'   The do.call can be somewhat slow, especially when working with large objects. 
#'   This function is based upon the suggestions from Hadley Wickham on the R mailing list, 
#'   see here. Also thanks to Tommy at StackOverflow for suggesting how to handle double 
#' . and triple colon operators, ::, further enhancing the function.
#'  
#' @param what   either a function or a non-empty character string naming the function to be called.
#' @param args   a list of arguments to the function call. The names attribute of args gives the argument names.
#' @param quote  a logical value indicating whether to quote the arguments.
#' @param envir  an environment within which to evaluate the call. This will be most useful if what is a character string and the arguments are symbols or quoted expressions.
#'
#' @author  Max Gorden 
#' 
#

int_fastDoCall <- function(what, args, quote = FALSE, envir = parent.frame()){
  if (quote)
    args <- lapply(args, enquote)
  
  if (is.null(names(args))){
    argn <- args
    args <- list()
  }else{
    # Add all the named arguments
    argn <- lapply(names(args)[names(args) != ""], as.name)
    names(argn) <- names(args)[names(args) != ""]
    # Add the unnamed arguments
    argn <- c(argn, args[names(args) == ""])
    args <- args[names(args) != ""]
  }
  
  if (class(what) == "character"){
    if(is.character(what)){
      fn <- strsplit(what, "[:]{2,3}")[[1]]
      what <- if(length(fn)==1) {
        get(fn[[1]], envir=envir, mode="function")
      } else {
        get(fn[[2]], envir=asNamespace(fn[[1]]), mode="function")
      }
    }
    call <- as.call(c(list(what), argn))
  }else if (class(what) == "function"){
    f_name <- deparse(substitute(what))
    call <- as.call(c(list(as.name(f_name)), argn))
    args[[f_name]] <- what
  }else if (class(what) == "name"){
    call <- as.call(c(list(what, argn)))
  }
  
  eval(call,
       envir = args,
       enclos = envir)
}
  
#' norm1
#' @author: unknown
norm <- function(x, log = FALSE) {
 if (all(is.na(x))) return(x)
  if (log) {
    x <- x - max(x, na.rm = TRUE)
    x <- exp(x)
  } else {

  stopifnot(all(x >= 0, na.rm = TRUE))
  x <- x / max(x, na.rm = TRUE)
  return(x / sum(x, na.rm = TRUE))
  }
}

    #==========================================================================================================================
#' @title PICS calculation using existing LD matrices for EUR ancestry 
#' 
#' @description The PICS algorithm calculates the most likely causal SNPs given the observed association signal
#'   at a locus. For an associated locus, enter the most highly-associated SNP 
#'   (referred to as the index SNP) and the strength of association. Using pre-computed LD matrices
#'   for EUR population, the algorithm identifies the SNPs that are most likely to be the 
#'   causal variants responsible for the association (PICS_Probability).
#'
#'   See \strong{Genetic and Epigenetic Fine-Mapping of Causal Variants in Autoimmune Disease} 
#'   by Kyle Kai-How Farh,et al. Nature 518, 337–343 (19 February 2015)
#' 
#' @param index_data_file Path to file containing index SNPs data with columns: snpID, pval, indication
#'   - snpID: SNP identifier in format chr:pos:ref:alt (e.g., 1:123456:A:C)
#'   - pval: p-value of association for the SNP
#'   - indication: disease/trait being studied
#' @param ld_block_file Path to the LD block definition file (BED format)
#' @param ld_dir Directory containing LD matrices (*.ld.gz) and variant info (*.vars)
#' @param nCORE_pics Define the number of cores for multiple tasks
#'   
#' @return A data frame containing the result of the PICS calculation
#' 
#' @author Kijoung Song (Modified by: ksong0924)
#' @date 2025-05-16
#' 
#' @export
                        
pics_calc <- function(index_data_file, ld_block_file, ld_dir, nCORE_pics = 1, 
                     max_matrix_size_gb = 2, chunk_size = 1000) {
  # Load required packages
  if(!require(dplyr)) install.packages("dplyr")
  if(!require(data.table)) install.packages("data.table")
  if(!require(foreach)) install.packages("foreach")
  if(!require(doParallel)) install.packages("doParallel")
  
  library(dplyr)
  library(data.table)
  library(foreach)
  library(doParallel)
  
  # Create a data frame to store missing positions
  missing_positions <- data.frame(
    snpID = character(),
    indication = character(),
    block_id = character(),
    timestamp = character(),
    stringsAsFactors = FALSE
  )
  
  # Process index SNPs
  message("Processing index SNPs...")
  
  # Parse SNP IDs to extract chromosome and position
  index.data <- fread(index_data_file)
  index.data <- index.data %>%
    mutate(
      raw_chrom = sapply(strsplit(as.character(snpID), ":"), function(x) x[1]),
      pos = as.numeric(sapply(strsplit(as.character(snpID), ":"), function(x) x[2])),
      ref = sapply(strsplit(as.character(snpID), ":"), function(x) if(length(x) > 2) x[3] else NA),
      alt = sapply(strsplit(as.character(snpID), ":"), function(x) if(length(x) > 3) x[4] else NA),
      chrom = raw_chrom,
      chr_for_block = paste0("chr", raw_chrom)
    )
  
  # Add constant EUR ancestry for consistency
  index.data$ancestry <- "EUR"
  
  # Add rsid1 column using the original snpID
  index.data$rsid1 <- index.data$snpID
  
  message("Reading LD block file...")
  ld_blocks <- fread(ld_block_file, header = FALSE)
  colnames(ld_blocks) <- c("chrom", "start", "end", "block_id")
  
  # Find which LD block contains each index SNP
  message("Determining LD blocks for index SNPs...")
  index.data$block_id <- NA
  
  for(i in 1:nrow(index.data)) {
    chr <- index.data$chr_for_block[i]
    pos <- index.data$pos[i]
    
    block <- ld_blocks[ld_blocks$chrom == chr & 
                      ld_blocks$start <= pos & 
                      ld_blocks$end >= pos, ]
    
    if(nrow(block) > 0) {
      index.data$block_id[i] <- block$block_id[1]
    }
  }
  
  # Remove SNPs that don't fall within an LD block
  index.data <- index.data[!is.na(index.data$block_id), ]
  
  if(nrow(index.data) == 0) {
    stop("None of the index SNPs fall within defined LD blocks.")
  }
  
  # Exclude chromosomes "X" and "XY"
  index.data <- index.data %>% filter(!(chrom %in% c("X", "Y", "XY")))
  
  # Group index SNPs by LD block
  message("Grouping index SNPs by LD block for efficient processing...")
  block_groups <- split(index.data, index.data$block_id)
  
  all_results <- list()
  
  for (block_name in names(block_groups)) {
    block_group <- block_groups[[block_name]]
    block_id <- unique(block_group$block_id)
    message(paste("Processing LD block:", block_id))
    
    tryCatch({
      # Construct filenames for LD matrix and variants
      ld_file <- file.path(ld_dir, paste0("ld_matrix_block_", block_id, ".ld.gz"))
      vars_file <- file.path(ld_dir, paste0("ld_matrix_block_", block_id, ".vars"))
      
      # Skip if files don't exist
      if(!file.exists(ld_file) || !file.exists(vars_file)) {
        message("Missing LD files for block ", block_id)
        next
      }
      
      # Check file size before loading
      file_size_gb <- file.size(ld_file) / (1024^3)
      message(paste("LD matrix file size:", round(file_size_gb, 2), "GB"))
      
      # Read variant information for this block
      vars_raw <- readLines(vars_file)
      
      if(length(vars_raw) == 0) {
        message("Empty variants file for block ", block_id)
        next
      }
      
      # Parse variant information into a data frame
      vars <- data.frame(
        variant_id = vars_raw,
        stringsAsFactors = FALSE
      )
      
      # Split variant ID into components
      variant_parts <- strsplit(vars_raw, ":")
      vars$chrom <- sapply(variant_parts, function(x) x[1])
      vars$position <- as.numeric(sapply(variant_parts, function(x) x[2]))
      
      # Process ref/alt information
      component_count <- sapply(variant_parts, length)
      if(any(component_count < 4)) {
        message("Warning: Some variants in block ", block_id, " don't have ref:alt information")
        
        vars$ref <- NA
        vars$alt <- NA
        
        vars$ref[component_count >= 3] <- sapply(variant_parts[component_count >= 3], function(x) x[3])
        vars$alt[component_count >= 4] <- sapply(variant_parts[component_count >= 4], function(x) x[4])
        
        missing_ref_alt <- which(component_count < 4)
        for(i in missing_ref_alt) {
          if(!is.na(vars$ref[i]) && !is.na(vars$alt[i])) {
            vars$variant_id[i] <- paste0(vars$chrom[i], ":", vars$position[i], ":", vars$ref[i], ":", vars$alt[i])
          }
        }
      } else {
        vars$ref <- sapply(variant_parts, function(x) x[3])
        vars$alt <- sapply(variant_parts, function(x) x[4])
      }
      
      # Create data frame to store block results
      block_results <- NULL
      
      # MEMORY MANAGEMENT: Handle large matrices differently
      if (file_size_gb > max_matrix_size_gb) {
        message("Large LD matrix detected, using chunked processing...")
        
        # Get positions of all SNPs in this block group
        block_positions <- block_group$pos
        
        # Pre-filter variants file to find positions of interest
        relevant_indices <- integer(0)
        for(pos in block_positions) {
          indices <- which(vars$position == pos)
          if(length(indices) > 0) {
            relevant_indices <- c(relevant_indices, indices[1]) # Take first if multiple
          } else {
            # Record missing position
            missing_idx <- which(block_group$pos == pos)
            if(length(missing_idx) > 0) {
              current_snpID <- block_group$snpID[missing_idx[1]]
              current_indication <- block_group$indication[missing_idx[1]]
              
              message(paste("Position not found in vars file. snpID:", current_snpID, 
                            "| indication:", current_indication))
              
              missing_positions <- rbind(missing_positions, data.frame(
                snpID = current_snpID,
                indication = current_indication,
                block_id = block_id,
                timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        if(length(relevant_indices) == 0) {
          message("No relevant positions found in this block, skipping...")
          next
        }
        
        # Process each index SNP in chunks
        for(idx in relevant_indices) {
          # Get the corresponding index SNP info
          curr_pos <- vars$position[idx]
          block_idx <- which(block_group$pos == curr_pos)
          if(length(block_idx) == 0) next
          
          # Process this specific index SNP
          message(paste("Processing index SNP at position", curr_pos))
          
          # Read only the needed row from the LD matrix
          # Use system command with awk for efficient row extraction
          tmp_file <- tempfile()
          cmd <- paste0("zcat ", ld_file, " | awk 'NR==", idx, "{print}' > ", tmp_file)
          system(cmd)
          
          if(file.size(tmp_file) == 0) {
            message("Failed to extract row from LD matrix")
            next
          }
          
          # Read the extracted row
          ld_values <- scan(tmp_file, quiet = TRUE)
          unlink(tmp_file) # Clean up temp file
          
          # Create a data frame with LD information for this index SNP
          snp_ld <- data.frame(
            chrom1 = block_group$chrom[block_idx],
            pos1 = curr_pos,
            chrom2 = vars$chrom,
            pos2 = vars$position,
            r2 = ld_values,
            variant_id = vars$variant_id,
            block_id = block_id,
            index_snp_id = block_group$snpID[block_idx],
            stringsAsFactors = FALSE
          )
          
          # Join with index SNP information
          snp_ld <- snp_ld %>% 
            inner_join(block_group[block_idx, ], by = c("chrom1" = "chrom", "pos1" = "pos"))
          
          if(is.null(block_results)) {
            block_results <- snp_ld
          } else {
            block_results <- bind_rows(block_results, snp_ld)
          }
        }
      } else {
        # Regular processing for smaller matrices
        message("LD matrix size is within memory limits, using standard processing...")
        
        # Read LD matrix for this block
        ld_matrix <- as.matrix(fread(ld_file))
        
        # Ensure matrix dimensions match
        if(nrow(ld_matrix) != length(vars_raw)) {
          message(paste("Mismatch between LD matrix dimensions and variant count for block", block_id))
          message(paste("LD matrix rows:", nrow(ld_matrix), "Variants:", length(vars_raw)))
          next
        }
        
        # Process each index SNP within this block
        for(i in 1:nrow(block_group)) {
          index_pos <- block_group$pos[i]
          index_row <- which(vars$position == index_pos)
          
          if(length(index_row) == 0) {
            # Record missing position
            current_snpID <- block_group$snpID[i]
            current_indication <- block_group$indication[i]
            
            message(paste("Position not found in vars file. snpID:", current_snpID, 
                          "| indication:", current_indication))
            
            missing_positions <- rbind(missing_positions, data.frame(
              snpID = current_snpID,
              indication = current_indication,
              block_id = block_id,
              timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              stringsAsFactors = FALSE
            ))
            
            next
          }
          
          # If multiple variants at the same position, use the first one
          if(length(index_row) > 1) {
            message(paste("Multiple variants at position", index_pos, "- using the first one"))
            index_row <- index_row[1]
          }
          
          # Extract LD values for this index SNP
          ld_values <- ld_matrix[index_row, ]
          
          # Create data frame with LD information
          snp_ld <- data.frame(
            chrom1 = block_group$chrom[i],
            pos1 = index_pos,
            chrom2 = vars$chrom,
            pos2 = vars$position,
            r2 = ld_values,
            variant_id = vars$variant_id,
            block_id = block_id,
            index_snp_id = block_group$snpID[i],
            stringsAsFactors = FALSE
          )
          
          # Join with index SNP information
          snp_ld <- snp_ld %>% 
            inner_join(block_group[i, ], by = c("chrom1" = "chrom", "pos1" = "pos"))
          
          if(is.null(block_results)) {
            block_results <- snp_ld
          } else {
            block_results <- bind_rows(block_results, snp_ld)
          }
        }
      }
      
      # Check if we have any results
      if(is.null(block_results) || nrow(block_results) == 0) {
        message("No SNPs in LD found for block ", block_id)
        next
      }
      
      # Filter SNPs in LD with r2 > 0.5 (Key step for PICS calculation)
      block_results <- block_results %>% 
        filter(r2 > 0.5) %>%
        mutate(pval = ifelse(pval == 0, 1e-250, pval)) %>%
        mutate(Mean = r2 * (-log10(pval)),
               SD = sqrt(1 - r2^3.2) * sqrt(-log10(pval))/2, 
               SD = ifelse(SD < 0.057, 0.057, SD),            
               prob = dnorm(-log10(pval), Mean, SD))
      
      # Calculate probability sum for normalization
      prob.sum <- block_results %>% 
        group_by(chrom1, pos1, pval) %>% 
        summarise(prob_sum = sum(prob), .groups = "drop")
      
      # Calculate PICS score (normalized probability)
      block_results <- block_results %>% 
        inner_join(prob.sum, by = c("chrom1", "pos1", "pval")) %>%
        mutate(pics = prob/prob_sum)
      
      # Add this block's results to our list
      all_results[[length(all_results) + 1]] <- block_results
      
      # Clean up to free memory
      gc()
      
    }, error = function(e) {
      message(paste("Error processing block", block_id, ":", e$message))
    })
  }
  
  # Write missing positions to a log file if any were found
  if(nrow(missing_positions) > 0) {
    log_filename <- paste0("missing_positions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
    message(paste("Writing", nrow(missing_positions), "missing positions to", log_filename))
    write.table(missing_positions, file = log_filename, sep = "\t", 
                row.names = FALSE, quote = FALSE)
  }
  
  # Check if we have any results at all
  if(length(all_results) == 0) {
    stop("No SNPs in LD were found. Please check your input data and LD files.")
  }
  
  # Combine all results
  all_results_df <- bind_rows(all_results)
  
  if(nrow(all_results_df) == 0) {
    stop("No SNPs passed the r2 > 0.5 filter. Please check your LD matrices.")
  }
  
  # Create final PICS result with the requested format
  message("Creating final PICS results...")
  
  # Get a unique list of index SNP groups
  unique_index_snps <- unique(all_results_df$index_snp_id)
  
  # Process each index SNP separately and combine results
  final_results <- data.frame()
  
  for(curr_index_snp in unique_index_snps) {
    # Filter results for the current index SNP only
    snp_results <- all_results_df %>%
      filter(index_snp_id == curr_index_snp) %>%
      # Select and rename columns for final output
      select(
        snpID = variant_id,
        pics,
        ancestry,
        indication,
        index_snp = index_snp_id,
        r2
      ) %>%
      # Remove duplicates
      distinct() %>%
      # Sort by indication and descending PICS score
      arrange(indication, desc(pics))
    
    # Add to final results
    final_results <- bind_rows(final_results, snp_results)
  }
  
  message("PICS calculation completed.")
  return(final_results)
}
