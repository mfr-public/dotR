# File: R/dnadot.R
# This script contains the R wrapper functions for the dotR package, which
# provides methods for estimating census population size (Nc) from SNP data.

# ============================================================================
# --- Package Setup and Imports ---
# ============================================================================

#' @useDynLib dotR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import stringr
#' @import stats
#' @import graphics
#' @import vcfR
#' @import doParallel
#' @importFrom foreach %dopar% %do% foreach registerDoSEQ
#' @importFrom methods new
NULL

# A simple helper function to provide a default value if an object is NULL.
# This is used for setting default grid search parameters.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ============================================================================
# --- Data Preparation and Utility Functions ---
# ============================================================================

#' Convert a data frame to a custom dotR genind object
#'
#' This function takes a standard data frame of genotypes and converts it into
#' a special 'genind' object optimized for use with the dotR C++ backend.
#' The key difference is that the genotype matrix (`@tab`) stores integer allele
#' indices (0-based) instead of allele counts, which is required by the C++ functions.
#'
#' @param df A data frame where rows are individuals and columns are loci.
#' @param sep The character separating alleles within a genotype (e.g., "/").
#' @param ploidy The ploidy of the organism (default is 2).
#' @param type The type of marker, passed to adegenet::df2genind (default "codom").
#' @return A genind object with a custom integer-based `@tab` slot.
#' @export
df2genind_dotR <- function(df, sep = "", ploidy = 2, type = "codom") {
  # Ensure the adegenet package is available, as it's essential for this conversion.
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop("Package 'adegenet' is required. Please install it.", call. = FALSE)
  }
  
  # First, create a standard genind object to easily extract allele information.
  genind_std <- adegenet::df2genind(
    X = df,
    sep = sep,
    ploidy = ploidy,
    type = type
  )

  # Get basic dimensions.
  n_ind <- adegenet::nInd(genind_std)
  n_loc <- adegenet::nLoc(genind_std)
  
  # Initialize the new numeric matrix that will hold 0-based allele indices.
  # Each locus will have 'ploidy' number of columns.
  gen_numeric_matrix <- matrix(as.integer(NA), nrow = n_ind, ncol = n_loc * ploidy)
  rownames(gen_numeric_matrix) <- adegenet::indNames(genind_std)
  
  # Loop through each locus to perform the conversion.
  for (j in 1:n_loc) {
    locus_genotypes_char <- as.character(df[[j]])
    locus_unique_alleles <- adegenet::alleles(genind_std)[[j]]
    
    # Proceed only if alleles are defined for the locus.
    if (length(locus_unique_alleles) > 0) {
      # Loop through each individual.
      for (i in 1:n_ind) {
        genotype_string <- locus_genotypes_char[i]
        
        # Handle missing data.
        if (is.na(genotype_string) || genotype_string == "") {
          alleles_in_genotype <- rep(NA_integer_, ploidy)
        } else {
          # Split the genotype string (e.g., "0101") into individual alleles ("01", "01").
          alleles_in_genotype_char <- strsplit(genotype_string, sep)[[1]]
          # Match the character alleles to the list of unique alleles to get their
          # 1-based index, then subtract 1 to make it 0-based for C++.
          alleles_in_genotype <- match(alleles_in_genotype_char, locus_unique_alleles) - 1L
        }
        
        # Place the 0-based indices into the correct columns of our new matrix.
        gen_numeric_matrix[i, ((j-1)*ploidy + 1):(j*ploidy)] <- alleles_in_genotype
      }
    }
  }

  # Replace the original @tab slot with our new integer-based matrix.
  genind_std@tab <- gen_numeric_matrix
  
  return(genind_std)
}


#' Read a VCF file and convert it to a dotR genind object
#'
#' A convenience wrapper that reads a VCF file, extracts genotype information,
#' and uses `df2genind_dotR` to create the custom genind object required by the package.
#'
#' @param vcf_file Path to the VCF file.
#' @param ploidy The ploidy of the organism (default is 2).
#' @return A genind object ready for use with `dnadot_snp`.
#' @export
read_vcf_to_genind_dotR <- function(vcf_file, ploidy = 2) {
  if (!requireNamespace("vcfR", quietly = TRUE)) stop("Package 'vcfR' is required.")
  if (!file.exists(vcf_file)) stop("VCF file not found: ", vcf_file)
  
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
  if (nrow(vcf@fix) == 0) stop("VCF file has no variant data.")
  
  gt <- vcfR::extract.gt(vcf, element = "GT", as.character = TRUE)
  if (is.null(gt)) stop("No GT field in VCF.")
  
  gt_t <- t(gt) # Transpose to have individuals as rows, loci as columns.
  
  # Create robust locus names, as VCF IDs can be missing or duplicated.
  loc_names <- vcf@fix[, "ID"]
  na_ids <- is.na(loc_names) | loc_names == "."
  if (any(na_ids) || any(duplicated(loc_names[!na_ids]))) {
      loc_names <- make.unique(paste0(vcf@fix[, "CHROM"], "_", vcf@fix[, "POS"]), sep = "_")
  }
  colnames(gt_t) <- loc_names
  
  # Create a data frame to hold the genotypes in the format "Allele1Allele2".
  df <- as.data.frame(matrix(NA_character_, nrow = nrow(gt_t), ncol = ncol(gt_t), dimnames = dimnames(gt_t)))
  ref <- vcf@fix[, "REF"]; alt <- vcf@fix[, "ALT"]
  
  # Loop through each locus (column).
  for (j in 1:ncol(gt_t)) {
      # Create a mapping from numeric VCF codes (0, 1, 2...) to actual allele strings (A, T, C...).
      alleles_map <- c(ref[j], strsplit(alt[j], ",")[[1]])
      # Loop through each individual (row).
      for (i in 1:nrow(gt_t)) {
          # Skip missing data ('.')
          if (is.na(gt_t[i, j]) || grepl("\\.", gt_t[i, j])) next
          # Convert numeric genotype (e.g., "0/1") to character genotype (e.g., "AT").
          indices <- as.integer(strsplit(gt_t[i, j], "[/|]")[[1]]) + 1
          df[i, j] <- paste(sort(alleles_map[indices]), collapse = "")
      }
  }
  
  # Use the main conversion function to create the final object.
  df2genind_dotR(df, sep = "", ploidy = ploidy, type = "codom")
}

# ============================================================================
# --- Main User-Facing Function ---
# ============================================================================

#' Estimate Census Size (Nc) from SNP data
#'
#' This is the main entry point for the dotR package. It takes a genind object
#' and applies a specified Minimum Distance Estimation (MDE) method to estimate Nc.
#'
#' @param genind_object A genind object, preferably created with `df2genind_dotR` or `read_vcf_to_genind_dotR`.
#' @param method The estimation method to use. Options are:
#'   "wasserstein", "sherwin", "nll" (for MDE methods), "bhd", "betabin" (fixed Fis),
#'   "betabin_dynamic" (dynamic Fis), and "mle".
#' @param ... Additional arguments passed to the specific internal method function.
#' @return A list containing the final Nc estimate (`N_est`) and detailed results from the analysis.
#' @export
dnadot_snp <- function(genind_object, method = "wasserstein", ...) {
    # --- Input Validation ---
    if (!inherits(genind_object, "genind")) stop("Input must be a 'genind' object.")
    if (adegenet::nInd(genind_object) == 0) {
        warning("Input genind_object has 0 individuals. Returning NA.")
        return(list(N_est = NA_real_, details = data.frame()))
    }

    # Package all arguments into a list for easy passing to internal functions.
    args <- list(genind_object = genind_object, ...)
    
    # --- Method Dispatch ---
    # Route the call to the correct internal function based on the 'method' argument.
    if (method == "betabin") {
        return(do.call(.dnadot_betabin_internal, args))
    } else if (method == "betabin_dynamic") {
        return(do.call(.dnadot_betabin_dynamic_internal, args))
    } else if (method == "bhd") {
        return(do.call(.dnadot_bhd_internal, args))
    } else if (method == "mle") {
        return(do.call(.dnadot_mle_internal, args))
    } else if (method %in% c("wasserstein", "sherwin", "nll")) {
        # Map user-friendly names to the internal C++ discrepancy types.
        discrepancy_map <- c("wasserstein" = "wasserstein", 
                               "sherwin" = "absolute_diff", 
                               "nll" = "neg_log_likelihood")
        args$discrepancy_type <- discrepancy_map[method]
        
        # Handle specific arguments for the "sherwin" method's regularization.
        if ("sherwin_reg_strength" %in% names(args)) args$reg_strength <- args$sherwin_reg_strength
        if ("sherwin_reg_N_threshold" %in% names(args)) args$reg_N_threshold <- args$sherwin_reg_N_threshold
        
        return(do.call(.dnadot_mde_internal, args))
    } else {
        stop("Invalid method. Must be 'wasserstein', 'betabin', 'betabin_dynamic', 'bhd', 'sherwin', 'nll', or 'mle'.")
    }
}

# ============================================================================
# --- Internal Helper and Setup Functions ---
# ============================================================================

#' Manually Calculate Fis from a dotR genind object
#'
#' This function bypasses the problematic `adegenet::summary.genind` S4 method
#' by calculating Hobs, Hexp, and Fis directly from the integer-based genotype
#' matrix in a dotR custom genind object. This is a more robust solution that
#' avoids environment-specific method dispatch failures.
#'
#' @param dotR_genind The custom genind object used by dotR.
#' @return A numeric vector of Fis values, one for each locus.
.calculate_fis_manually <- function(dotR_genind) {
    n_ind <- adegenet::nInd(dotR_genind)
    n_loc <- adegenet::nLoc(dotR_genind)
    ploidy <- adegenet::ploidy(dotR_genind)[1]
    
    # Check for diploidy, as the Hobs calculation assumes this.
    if (ploidy != 2) {
        stop("Manual Fis calculation currently only supports diploid data.", call. = FALSE)
    }
    
    fis_per_locus <- numeric(n_loc)
    
    for (l in 1:n_loc) {
        # Extract the integer allele data for the current locus.
        start_col <- (l - 1) * ploidy + 1
        end_col <- l * ploidy
        locus_data <- dotR_genind@tab[, start_col:end_col, drop = FALSE]
        
        # --- 1. Calculate Observed Heterozygosity (Hobs) ---
        # Hobs is the proportion of typed individuals that are heterozygous.
        n_typed_ind <- 0
        n_het <- 0
        for (i in 1:n_ind) {
            alleles <- locus_data[i, ]
            if (any(is.na(alleles))) next # Skip individuals with missing data at this locus.
            
            n_typed_ind <- n_typed_ind + 1
            # An individual is heterozygous if its two allele indices are different.
            if (alleles[1] != alleles[2]) {
                n_het <- n_het + 1
            }
        }
        hobs <- if (n_typed_ind > 0) n_het / n_typed_ind else 0
        
        # --- 2. Calculate Expected Heterozygosity (Hexp) ---
        # Hexp is calculated from allele frequencies (1 - sum of squared freqs).
        all_alleles_vector <- as.vector(locus_data)
        all_alleles_vector <- all_alleles_vector[!is.na(all_alleles_vector)]
        
        # If no alleles were typed for this locus, Fis is undefined.
        if (length(all_alleles_vector) == 0) {
            fis_per_locus[l] <- NA
            next
        }
        
        # Calculate allele frequencies from the counts.
        allele_counts <- table(all_alleles_vector)
        allele_freqs <- allele_counts / sum(allele_counts)
        
        hexp <- 1 - sum(allele_freqs^2)
        
        # --- 3. Calculate Fis ---
        # Fis = 1 - (Hobs / Hexp). Avoid division by zero.
        fis_per_locus[l] <- if (hexp > 0) 1 - (hobs / hexp) else 0
    }
    
    return(fis_per_locus)
}


#' Set up the parallel backend for analysis.
#'
#' This function checks if parallel processing is requested and possible,
#' and registers the appropriate backend with the `foreach` package.
#'
#' @param n_cores Number of cores to use.
#' @param n_items The number of items to iterate over (e.g., number of loci).
#' @return A boolean indicating whether the analysis will run in parallel.
.setup_parallel <- function(n_cores, n_items) {
    run_parallel <- (n_cores > 1 && n_items > 1)
    if (run_parallel) {
        if (!requireNamespace("doParallel", quietly = TRUE) || !requireNamespace("foreach", quietly = TRUE)) {
            warning("doParallel or foreach not found. Running sequentially.")
            run_parallel <- FALSE
        } else {
            doParallel::registerDoParallel(cores = n_cores)
        }
    }
    # If not running in parallel, register the sequential backend.
    if (!run_parallel) foreach::registerDoSEQ()
    return(run_parallel)
}

#' Set up the grid search parameters.
#'
#' This function defines the grid of hypothesized Nc and allele frequency (p)
#' values that will be searched to find the best estimate. It also calculates
#' the size of the jackknife subsamples.
#'
#' @param n_individuals Total number of individuals in the sample.
#' @param jackknife_proportion The proportion of individuals to use in each subsample.
#' @param num_N_hypothesized The number of Nc values to test.
#' @param N_try_min The minimum Nc value to test.
#' @param N_try_max The maximum Nc value to test.
#' @param num_p_hypothesized The number of p values to test.
#' @param p_hypothesized_range A vector with the min and max p to test.
#' @return A list containing the vectors of Nc and p values to try, and the subsample size.
.setup_grid_params <- function(n_individuals, jackknife_proportion, num_N_hypothesized, N_try_min, N_try_max, num_p_hypothesized, p_hypothesized_range) {
    # Calculate the number of individuals in each jackknife subsample.
    n_subsample_indiv <- floor(jackknife_proportion * n_individuals)
    
    # The smallest possible Nc is the size of the subsample.
    min_N_val <- max(1, n_subsample_indiv)
    if (jackknife_proportion == 1) min_N_val <- max(min_N_val, n_individuals)
    
    # Set sensible defaults for the Nc search range if not provided by the user.
    N_try_min_eff <- N_try_min %||% max(min_N_val, floor(n_individuals * 0.5), 2)
    N_try_max_eff <- N_try_max %||% floor(n_individuals * 3)
    
    # Ensure the max is greater than the min.
    if (N_try_max_eff <= N_try_min_eff) {
        N_try_max_eff <- N_try_min_eff + num_N_hypothesized * 2
        if (!is.null(N_try_min) || !is.null(N_try_max)) {
             warning("N_try_max was less than or equal to N_try_min. Adjusting N_try_max.")
        }
    }
    
    # Create the sequences of hypothesized values.
    N_try_values <- round(seq(N_try_min_eff, N_try_max_eff, length.out = num_N_hypothesized))
    p_try_values <- seq(p_hypothesized_range[1], p_hypothesized_range[2], length.out = num_p_hypothesized)
    
    list(N_try = N_try_values, p_try = p_try_values, n_subsample = n_subsample_indiv)
}

#' Process and summarize results from all loci.
#'
#' This function takes the raw list of results from the `foreach` loop,
#' filters out errors, selects the most informative loci/alleles, and
#' calculates the final, single estimate of Nc.
#'
#' @param all_loci_results A list containing the results from each locus-allele combination.
#' @param method_name The name of the method used, for reporting.
#' @return A list containing the final Nc estimate and detailed data frames.
.process_results <- function(all_loci_results, method_name) {
    # Remove any runs that resulted in an error.
    valid_results <- Filter(function(x) is.list(x) && !inherits(x, "error"), all_loci_results)
    
    if (length(valid_results) == 0) {
        warning("No valid results were obtained from any locus.", call. = FALSE)
        return(list(N_est = NA_real_, selected_loci_details = data.frame(), all_runs_summary = data.frame(), method = method_name))
    }
    
    # Combine the list of lists into a single, clean data frame.
    results_df <- do.call(rbind, lapply(valid_results, as.data.frame))

    # Filter for runs that produced a finite Nc estimate and a valid CV.
    # The Coefficient of Variation (CV) of subsample counts is a key metric for
    # selecting informative loci. A low CV suggests the allele count is stable
    # across subsamples, making it a reliable signal.
    results_df_filtered <- results_df[is.finite(results_df$cv_subsample_counts) & is.finite(results_df$N_est), ]
    
    N_est_final <- NA_real_
    selected_df <- data.frame()
    
    if (nrow(results_df_filtered) > 0) {
      # Select the top 10% of runs with the lowest CV.
      num_to_select <- max(1, floor(0.10 * nrow(results_df_filtered)))
      selected_df <- head(results_df_filtered[order(results_df_filtered$cv_subsample_counts), ], num_to_select)
      
      if (nrow(selected_df) > 0) {
        # The final estimate is the mean of the Nc estimates from these top runs.
        N_est_final <- mean(selected_df$N_est, na.rm = TRUE)
      }
    }
    
    return(list(N_est = N_est_final, selected_loci_details = selected_df, all_runs_summary = results_df, method = method_name))
}

# ============================================================================
# --- Internal Method-Specific Handlers ---
# These functions contain the core logic for each estimation method.
# ============================================================================

#' Internal handler for standard MDE methods (Wasserstein, Sherwin, NLL).
#' @private
.dnadot_mde_internal <- function(genind_object, discrepancy_type, jackknife_proportion = 0.8,
                                 num_N_hypothesized = 11, N_try_min = NULL, N_try_max = NULL,
                                 num_p_hypothesized = 20, p_hypothesized_range = c(0.01, 0.99), 
                                 reg_strength = 0.0, reg_N_threshold = Inf, n_cores = 1, ...) {
    
    method_name <- paste0("MDE (", discrepancy_type, ")")
    message("Running ", method_name, " method...")

    # --- Setup ---
    n_individuals <- adegenet::nInd(genind_object)
    ploidy_val <- unique(adegenet::ploidy(genind_object))
    if (length(ploidy_val) > 1) stop("Mixed ploidy not currently supported.")
    n_loci <- adegenet::nLoc(genind_object)
    
    # Extract necessary data from the genind object to pass to the loop.
    tab_matrix <- genind_object@tab
    loc_names_vec <- adegenet::locNames(genind_object)
    alleles_list <- adegenet::alleles(genind_object)

    # Define the grid search space.
    grid_params <- .setup_grid_params(n_individuals, jackknife_proportion, num_N_hypothesized, N_try_min, N_try_max, num_p_hypothesized, p_hypothesized_range)

    # Configure parallel processing.
    run_parallel <- .setup_parallel(n_cores, n_loci)
    if (run_parallel) on.exit(foreach::registerDoSEQ(), add = TRUE)

    # --- Main Loop ---
    # Use '%dopar%' for parallel or '%do%' for sequential execution.
    `%op%` <- if (run_parallel) foreach::`%dopar%` else foreach::`%do%`
    all_loci_results <- foreach::foreach(locus_idx = 1:n_loci, .combine = 'c', .errorhandling = 'pass') %op% {
        locus_alleles <- alleles_list[[locus_idx]]
        # Skip monomorphic loci.
        if (length(locus_alleles) < 2) return(NULL) 
        
        locus_results_list <- list()
        # Extract the genotype data for just the current locus.
        start_col <- (locus_idx - 1) * ploidy_val + 1
        end_col <- locus_idx * ploidy_val
        genotypes_locus <- tab_matrix[, start_col:end_col, drop = FALSE]
        
        # The analysis is run for each allele at the locus.
        for (allele_idx_zero_based in 0:(length(locus_alleles) - 1)) {
            # Ensure the matrix is of integer type for C++.
            if(!is.integer(genotypes_locus)) storage.mode(genotypes_locus) <- "integer"
            
            # Call the core C++ function to perform the MDE grid search.
            cpp_result <- dotR_dnadot_mde_cpp(genotypes_locus, n_individuals, grid_params$n_subsample,
                                              grid_params$N_try, grid_params$p_try, allele_idx_zero_based,
                                              discrepancy_type, reg_strength, reg_N_threshold,
                                              debug = FALSE)
            
            # --- Post-processing for this allele ---
            counts <- cpp_result$subsample_target_allele_counts
            mean_counts <- mean(counts, na.rm = TRUE)
            cv_val <- if (length(counts) > 1 && mean_counts > 0) sd(counts, na.rm = TRUE) / mean_counts else NA_real_
            
            allele_char <- locus_alleles[allele_idx_zero_based + 1]
            
            # Store results in a structured list.
            result_for_allele <- list(locus_name = loc_names_vec[locus_idx],
                                      target_allele = allele_char, 
                                      N_est = cpp_result$N_est, 
                                      p_est = cpp_result$p_est,
                                      min_discrepancy = cpp_result$min_discrepancy, 
                                      cv_subsample_counts = cv_val)
            
            locus_results_list[[length(locus_results_list) + 1]] <- result_for_allele
        }
        locus_results_list
    }

    # --- Finalization ---
    return(.process_results(all_loci_results, method_name))
}


#' Internal handler for the Buffered Histogram Discrepancy (BHD) method.
#' @private
.dnadot_bhd_internal <- function(genind_object, jackknife_proportion = 0.8, bhd_buffer_width = 1,
                                 num_N_hypothesized = 11, N_try_min = NULL, N_try_max = NULL,
                                 num_p_hypothesized = 20, p_hypothesized_range = c(0.01, 0.99), n_cores = 1, ...) {
    
    method_name <- "BHD"
    message("Running BHD method...")
    
    # --- Setup (similar to MDE) ---
    n_individuals <- adegenet::nInd(genind_object)
    ploidy_val <- unique(adegenet::ploidy(genind_object))
    if (length(ploidy_val) > 1) stop("Mixed ploidy not currently supported.")
    n_loci <- adegenet::nLoc(genind_object)
    tab_matrix <- genind_object@tab
    loc_names_vec <- adegenet::locNames(genind_object)
    alleles_list <- adegenet::alleles(genind_object)
    grid_params <- .setup_grid_params(n_individuals, jackknife_proportion, num_N_hypothesized, N_try_min, N_try_max, num_p_hypothesized, p_hypothesized_range)
    run_parallel <- .setup_parallel(n_cores, n_loci)
    if (run_parallel) on.exit(foreach::registerDoSEQ(), add = TRUE)

    # --- Main Loop ---
    `%op%` <- if (run_parallel) foreach::`%dopar%` else foreach::`%do%`
    all_loci_results <- foreach::foreach(locus_idx = 1:n_loci, .combine = 'c', .errorhandling = 'pass') %op% {
        locus_alleles <- alleles_list[[locus_idx]]
        if (length(locus_alleles) < 2) return(NULL)
        
        locus_results_list <- list()
        start_col <- (locus_idx - 1) * ploidy_val + 1
        end_col <- locus_idx * ploidy_val
        genotypes_locus <- tab_matrix[, start_col:end_col, drop = FALSE]
        
        for (allele_idx_zero_based in 0:(length(locus_alleles) - 1)) {
            if(!is.integer(genotypes_locus)) storage.mode(genotypes_locus) <- "integer"
            
            # Call the specific C++ function for the BHD method.
            cpp_result <- dotR_dnadot_bhd_cpp(genotypes_locus, n_individuals, grid_params$n_subsample,
                                              grid_params$N_try, grid_params$p_try, allele_idx_zero_based, as.integer(bhd_buffer_width))
            
            counts <- cpp_result$subsample_target_allele_counts
            mean_counts <- mean(counts, na.rm = TRUE)
            cv_val <- if (length(counts) > 1 && mean_counts > 0) sd(counts, na.rm = TRUE) / mean_counts else NA_real_
            
            allele_char <- locus_alleles[allele_idx_zero_based + 1]
            
            result_for_allele <- list(locus_name = loc_names_vec[locus_idx],
                                      target_allele = allele_char, 
                                      N_est = cpp_result$N_est, 
                                      p_est = cpp_result$p_est,
                                      min_discrepancy = cpp_result$min_discrepancy, 
                                      cv_subsample_counts = cv_val)
            locus_results_list[[length(locus_results_list) + 1]] <- result_for_allele
        }
        locus_results_list
    }

    return(.process_results(all_loci_results, method_name))
}

#' Internal handler for the Beta-binomial method with a fixed, global Fis.
#' @private
.dnadot_betabin_internal <- function(genind_object, jackknife_proportion = 0.8,
                                     num_N_hypothesized = 11, N_try_min = NULL, N_try_max = NULL,
                                     num_p_hypothesized = 20, p_hypothesized_range = c(0.01, 0.99), 
                                     n_cores = 1, ...) {
    
    method_name <- "MDE (Beta-binomial, fixed Fis)"
    message("Running ", method_name, " method...")

    # --- Fis Calculation (Robust Manual Method) ---
    # This now calls our new helper function instead of relying on summary().
    fis_per_locus <- .calculate_fis_manually(genind_object)
    
    # Calculate the global average Fis from the per-locus values.
    valid_fis <- fis_per_locus[is.finite(fis_per_locus)]
    global_fis <- if (length(valid_fis) > 0) mean(valid_fis, na.rm = TRUE) else 0.0
    global_fis <- max(0, global_fis) # Fis cannot be negative in the model.
    message(paste("Using Global Fis estimate:", round(global_fis, 4)))

    # --- Setup (similar to MDE) ---
    n_individuals <- adegenet::nInd(genind_object)
    ploidy_val <- unique(adegenet::ploidy(genind_object))
    if (length(ploidy_val) > 1) stop("Mixed ploidy not currently supported.")
    n_loci <- adegenet::nLoc(genind_object)
    tab_matrix <- genind_object@tab
    loc_names_vec <- adegenet::locNames(genind_object)
    alleles_list <- adegenet::alleles(genind_object)
    grid_params <- .setup_grid_params(n_individuals, jackknife_proportion, num_N_hypothesized, N_try_min, N_try_max, num_p_hypothesized, p_hypothesized_range)
    run_parallel <- .setup_parallel(n_cores, n_loci)
    if (run_parallel) on.exit(foreach::registerDoSEQ(), add = TRUE)

    # --- Main Loop ---
    `%op%` <- if (run_parallel) foreach::`%dopar%` else foreach::`%do%`
    all_loci_results <- foreach::foreach(locus_idx = 1:n_loci, .combine = 'c', .errorhandling = 'pass') %op% {
        locus_alleles <- alleles_list[[locus_idx]]
        if (length(locus_alleles) < 2) return(NULL) 
        
        locus_results_list <- list()
        start_col <- (locus_idx - 1) * ploidy_val + 1
        end_col <- locus_idx * ploidy_val
        genotypes_locus <- tab_matrix[, start_col:end_col, drop = FALSE]
        
        for (allele_idx_zero_based in 0:(length(locus_alleles) - 1)) {
            if(!is.integer(genotypes_locus)) storage.mode(genotypes_locus) <- "integer"
            
            # Call the C++ function for beta-binomial, passing the fixed global Fis.
            cpp_result <- dotR_dnadot_betabin_cpp(genotypes_locus, n_individuals, grid_params$n_subsample,
                                                  grid_params$N_try, grid_params$p_try, allele_idx_zero_based,
                                                  global_fis)
            
            counts <- cpp_result$subsample_target_allele_counts
            mean_counts <- mean(counts, na.rm = TRUE)
            cv_val <- if (length(counts) > 1 && mean_counts > 0) sd(counts, na.rm = TRUE) / mean_counts else NA_real_
            
            allele_char <- locus_alleles[allele_idx_zero_based + 1]
            
            result_for_allele <- list(locus_name = loc_names_vec[locus_idx],
                                      target_allele = allele_char, N_est = cpp_result$N_est, p_est = cpp_result$p_est,
                                      min_discrepancy = cpp_result$min_discrepancy, cv_subsample_counts = cv_val)
            locus_results_list[[length(locus_results_list) + 1]] <- result_for_allele
        }
        locus_results_list
    }

    return(.process_results(all_loci_results, method_name))
}


#' Internal handler for the Beta-binomial method with a dynamic Fis.
#' This is the most complex method, involving a two-step calibration process.
#' @private
.dnadot_betabin_dynamic_internal <- function(genind_object, jackknife_proportion = 0.8,
                                             num_N_hypothesized = 11, N_try_min = NULL, N_try_max = NULL,
                                             num_p_hypothesized = 20, p_hypothesized_range = c(0.01, 0.99), 
                                             n_cores = 1, ...) {
    
    method_name <- "MDE (Beta-binomial, dynamic Fis)"
    message("Running ", method_name, " method...")

    # --- STEP 1: CALIBRATE THE MODEL ---
    message("Step 1: Calibrating model...")
    
    # --- Fis Calculation (Robust Manual Method) ---
    fis_per_locus <- .calculate_fis_manually(genind_object)
    
    # Calculate the observed global Fis from the data.
    valid_fis <- fis_per_locus[is.finite(fis_per_locus)]
    global_fis_obs <- if (length(valid_fis) > 0) mean(valid_fis, na.rm = TRUE) else 0.0
    
    if (global_fis_obs <= 0) {
        stop("Dynamic beta-binomial method requires a positive global Fis estimate from the data. The observed Fis was <= 0.", call. = FALSE)
    }
    message(paste("  - Observed Global Fis:", round(global_fis_obs, 4)))

    # Get a quick, initial guess for Nc using a faster method (Wasserstein).
    message("  - Getting initial Nc guess using 'wasserstein' method...")
    initial_guess_args <- list(genind_object = genind_object,
                               discrepancy_type = "wasserstein", 
                               jackknife_proportion = jackknife_proportion, num_N_hypothesized = num_N_hypothesized,
                               N_try_min = N_try_min, N_try_max = N_try_max, num_p_hypothesized = num_p_hypothesized,
                               p_hypothesized_range = p_hypothesized_range, n_cores = n_cores)
    initial_res <- do.call(.dnadot_mde_internal, initial_guess_args)
    Nc_initial <- initial_res$N_est

    if (is.na(Nc_initial) || !is.finite(Nc_initial)) {
        stop("Could not obtain a valid initial Nc estimate from the 'wasserstein' method to calibrate the model.", call. = FALSE)
    }
    message(paste("  - Initial Nc guess:", round(Nc_initial)))

    # Using Wright's Island Model formula (Fis = 1 / (4Nm + 1)), solve for m.
    migration_rate_m <- (1 - global_fis_obs) / (4 * Nc_initial * global_fis_obs)
    if (is.na(migration_rate_m) || !is.finite(migration_rate_m) || migration_rate_m <= 0) {
        stop("Could not calculate a valid migration rate 'm' for calibration. Check observed Fis and initial Nc values.", call. = FALSE)
    }
    message(paste("  - Calibrated migration rate (m):", format(migration_rate_m, scientific = TRUE, digits = 4)))
    
    
    # --- STEP 2: RUN THE DYNAMIC GRID SEARCH ---
    message("Step 2: Running dynamic grid search...")
    
    # --- Setup (similar to other methods) ---
    n_individuals <- adegenet::nInd(genind_object)
    ploidy_val <- unique(adegenet::ploidy(genind_object))
    if (length(ploidy_val) > 1) stop("Mixed ploidy not currently supported.")
    n_loci <- adegenet::nLoc(genind_object)
    tab_matrix <- genind_object@tab
    loc_names_vec <- adegenet::locNames(genind_object)
    alleles_list <- adegenet::alleles(genind_object)
    grid_params <- .setup_grid_params(n_individuals, jackknife_proportion, num_N_hypothesized, N_try_min, N_try_max, num_p_hypothesized, p_hypothesized_range)
    run_parallel <- .setup_parallel(n_cores, n_loci)
    if (run_parallel) on.exit(foreach::registerDoSEQ(), add = TRUE)

    # --- Main Loop ---
    `%op%` <- if (run_parallel) foreach::`%dopar%` else foreach::`%do%`
    all_loci_results <- foreach::foreach(locus_idx = 1:n_loci, .combine = 'c', .errorhandling = 'pass') %op% {
        locus_alleles <- alleles_list[[locus_idx]]
        if (length(locus_alleles) < 2) return(NULL) 
        
        locus_results_list <- list()
        start_col <- (locus_idx - 1) * ploidy_val + 1
        end_col <- locus_idx * ploidy_val
        genotypes_locus <- tab_matrix[, start_col:end_col, drop = FALSE]
        
        for (allele_idx_zero_based in 0:(length(locus_alleles) - 1)) {
            if(!is.integer(genotypes_locus)) storage.mode(genotypes_locus) <- "integer"
            
            # Call the C++ function for the dynamic beta-binomial method.
            cpp_result <- dotR_dnadot_betabin_dynamic_cpp(genotypes_locus, n_individuals, grid_params$n_subsample,
                                                          grid_params$N_try, grid_params$p_try, allele_idx_zero_based,
                                                          migration_rate_m)
            
            counts <- cpp_result$subsample_target_allele_counts
            mean_counts <- mean(counts, na.rm = TRUE)
            cv_val <- if (length(counts) > 1 && mean_counts > 0) sd(counts, na.rm = TRUE) / mean_counts else NA_real_
            
            allele_char <- locus_alleles[allele_idx_zero_based + 1]
            
            result_for_allele <- list(locus_name = loc_names_vec[locus_idx],
                                      target_allele = allele_char, N_est = cpp_result$N_est, p_est = cpp_result$p_est,
                                      min_discrepancy = cpp_result$min_discrepancy, cv_subsample_counts = cv_val)
            locus_results_list[[length(locus_results_list) + 1]] <- result_for_allele
        }
        locus_results_list
    }

    return(.process_results(all_loci_results, method_name))
}


#' Internal handler for the Maximum Likelihood Estimation (MLE) method.
#' NOTE: This method is flagged as experimental and potentially flawed.
#' @private
.dnadot_mle_internal <- function(genind_object, jackknife_proportion = 0.8,
                                 mle_target_alleles = NULL, num_N_for_aic = 30,
                                 N_hypothesized_range_factor = 3, n_cores = 1, ...) {
    
    message("Running MLE/AIC method...")
    warning("CRITICAL WARNING: The MLE implementation in dotR is considered experimental and may be statistically flawed.", call. = FALSE)

    # --- Setup ---
    n_individuals <- adegenet::nInd(genind_object)
    ploidy_val <- unique(adegenet::ploidy(genind_object))
    if (length(ploidy_val) > 1) stop("Mixed ploidy not currently supported.")
    n_loci <- adegenet::nLoc(genind_object)
    
    # Auto-select the most common alleles as targets if none are provided.
    if (is.null(mle_target_alleles)) {
        all_alleles_flat <- unlist(adegenet::alleles(genind_object))
        allele_counts <- table(all_alleles_flat)
        mle_target_alleles <- names(sort(allele_counts, decreasing = TRUE))[1:min(length(allele_counts), 5)]
        message("Auto-selecting target alleles for MLE: ", paste(mle_target_alleles, collapse=", "))
    }

    if (length(mle_target_alleles) == 0) {
        warning("No valid MLE target alleles remain, returning empty results.")
        return(list(method_used = "mle_aic", mle_aic_results = list()))
    }

    run_parallel <- .setup_parallel(n_cores, length(mle_target_alleles))
    if (run_parallel) on.exit(foreach::registerDoSEQ(), add = TRUE)

    # Define the likelihood function to be optimized.
    likelihood_func_mle <- function(N_param, target_allele_char, all_loci_genotypes_tab, 
                                    alleles_list, current_ploidy, current_n_loci, 
                                    current_n_individuals, current_jackknife_proportion) {
      
      min_N_param <- max(1, floor(current_jackknife_proportion * current_n_individuals))
      if (N_param < min_N_param || !is.finite(N_param)) return(Inf)
      
      total_log_likelihood_sum <- 0
      n_subsample_mle <- floor(current_jackknife_proportion * current_n_individuals)
      if (n_subsample_mle < 1) return(Inf)
      
      # This likelihood calculation sums across all loci for a single target allele.
      for (locus_idx_mle in 1:current_n_loci) {
          locus_alleles_here <- alleles_list[[locus_idx_mle]]
          target_allele_idx_mle <- match(target_allele_char, locus_alleles_here) - 1L
          
          if (is.na(target_allele_idx_mle)) next # Skip if allele not at this locus.
          
          start_col <- (locus_idx_mle - 1) * current_ploidy + 1
          end_col <- locus_idx_mle * current_ploidy
          genotypes_current_locus <- all_loci_genotypes_tab[, start_col:end_col, drop = FALSE]

          # Estimate allele frequency from the full sample.
          K_target_in_sample_locus <- sum(genotypes_current_locus == target_allele_idx_mle, na.rm = TRUE)
          total_alleles_in_sample_locus <- sum(!is.na(genotypes_current_locus))
          pop_allele_freq <- if (total_alleles_in_sample_locus > 0) K_target_in_sample_locus / total_alleles_in_sample_locus else 0
          
          # Hypothesize the total number of target alleles in the population.
          K_genes_hypothesized_N <- round(pop_allele_freq * N_param * current_ploidy)
          
          # Sum the log-likelihood across all jackknife subsamples.
          for (start_idx_mle in 1:(current_n_individuals - n_subsample_mle + 1)) {
            subsample_genotypes_locus <- genotypes_current_locus[start_idx_mle:(start_idx_mle + n_subsample_mle - 1), , drop = FALSE]
            k_observed <- sum(subsample_genotypes_locus == target_allele_idx_mle, na.rm = TRUE)
            n_drawn <- sum(!is.na(subsample_genotypes_locus))
            if (n_drawn == 0) next
            
            # Calculate probability using the hypergeometric distribution.
            prob <- dotR_hypergeo_cpp(k_observed, n_drawn, K_genes_hypothesized_N, as.integer(N_param * current_ploidy))
            total_log_likelihood_sum <- total_log_likelihood_sum + (if (prob > 1e-12) log(prob) else -Inf)
          }
      }
      # Return the negative log-likelihood for minimization.
      return(if(is.finite(total_log_likelihood_sum)) -total_log_likelihood_sum else Inf)
    }

    # A wrapper function to perform the full analysis for a single target allele.
    mle_allele_analysis_func <- function(target_allele_char) {
        n_subsample_indiv <- floor(n_individuals * jackknife_proportion)
        lower_bound <- max(1, n_subsample_indiv)
        if (jackknife_proportion == 1) lower_bound <- max(lower_bound, n_individuals)
        upper_bound <- floor(n_individuals * N_hypothesized_range_factor * 2)
        if (upper_bound <= lower_bound) upper_bound <- lower_bound + 100
        
        # Find the MLE for Nc using a numerical optimizer.
        initial_N <- lower_bound + (upper_bound - lower_bound) / 3
        mle_N_val <- NA_real_
        
        optim_res <- try(stats::optim(par=initial_N, fn=likelihood_func_mle,
                                      target_allele_char=target_allele_char, all_loci_genotypes_tab=genind_object@tab,
                                      alleles_list = adegenet::alleles(genind_object),
                                      current_ploidy=ploidy_val, current_n_loci=n_loci, current_n_individuals=n_individuals,
                                      current_jackknife_proportion=jackknife_proportion, method="Brent", lower=lower_bound, upper=upper_bound), silent=TRUE)
        
        if (!inherits(optim_res, "try-error") && optim_res$convergence == 0) mle_N_val <- optim_res$par

        # --- AIC Model Averaging ---
        # Calculate AIC scores for a grid of Nc values to get a model-averaged estimate.
        N_aic_grid <- round(seq(lower_bound, upper_bound, length.out = num_N_for_aic))
        log_lik_vals <- sapply(N_aic_grid, function(N) -likelihood_func_mle(N, target_allele_char, genind_object@tab, adegenet::alleles(genind_object), ploidy_val, n_loci, n_individuals, jackknife_proportion))
        
        aic_vals <- 2 * 1 - 2 * log_lik_vals # k=1 parameter (Nc)
        
        averaged_N <- NA_real_; aic_weights <- rep(0, length(N_aic_grid))
        finite_aic <- is.finite(aic_vals)
        if(any(finite_aic)){
            min_aic <- min(aic_vals[finite_aic], na.rm=TRUE)
            delta_aic <- aic_vals[finite_aic] - min_aic
            exp_delta <- exp(-0.5 * delta_aic)
            sum_exp_delta <- sum(exp_delta, na.rm=TRUE)
            if(is.finite(sum_exp_delta) && sum_exp_delta > 0){
                aic_weights[finite_aic] <- exp_delta / sum_exp_delta
                averaged_N <- sum(N_aic_grid[finite_aic] * aic_weights[finite_aic], na.rm=TRUE)
            }
        }
        
        list(mle_N=mle_N_val, averaged_N=averaged_N, hypothesized_N_for_aic=N_aic_grid, 
             log_likelihoods_for_aic=log_lik_vals, aic_weights=aic_weights)
    }

    # Run the analysis for each target allele in parallel.
    `%op%` <- if (run_parallel) foreach::`%dopar%` else foreach::`%do%`
    mle_results_list <- foreach::foreach(target_allele = mle_target_alleles, .combine = 'list', .multicombine = TRUE, .errorhandling = 'pass') %op% {
        mle_allele_analysis_func(target_allele)
    }

    if (length(mle_results_list) == 0) {
        return(list(method_used = "mle_aic", mle_aic_results = list()))
    }
    
    try({
       names(mle_results_list) <- mle_target_alleles
    }, silent = TRUE)
    
    valid_mle_results <- Filter(function(x) !inherits(x, "error") && is.list(x), mle_results_list)

    return(list(method_used = "mle_aic", mle_aic_results = valid_mle_results))
}