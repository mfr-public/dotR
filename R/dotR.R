# File: R/dotR.R
# This script contains the R wrapper functions for the dotR package.

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
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ============================================================================
# --- Data Preparation and Utility Functions (Ported from Beta) ---
# ============================================================================

#' Convert a data frame to a custom dotR genind object
#' (Ported from dnadot-beta.R)
#' @export
df2genind_dotR <- function(df, sep = "", ploidy = 2, type = "codom") {
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop("Package 'adegenet' is required. Please install it.", call. = FALSE)
  }
  
  genind_std <- adegenet::df2genind(X = df, sep = sep, ploidy = ploidy, type = type)

  n_ind <- adegenet::nInd(genind_std); n_loc <- adegenet::nLoc(genind_std)
  gen_numeric_matrix <- matrix(as.integer(NA), nrow = n_ind, ncol = n_loc * ploidy)
  rownames(gen_numeric_matrix) <- adegenet::indNames(genind_std)
  
  for (j in 1:n_loc) {
    locus_genotypes_char <- as.character(df[[j]])
    locus_unique_alleles <- adegenet::alleles(genind_std)[[j]]
    
    if (length(locus_unique_alleles) > 0) {
      for (i in 1:n_ind) {
        genotype_string <- locus_genotypes_char[i]
        if (is.na(genotype_string) || genotype_string == "") {
          alleles_in_genotype <- rep(NA_integer_, ploidy)
        } else {
          alleles_in_genotype_char <- strsplit(genotype_string, sep)[[1]]
          # Convert to 0-based index for C++.
          alleles_in_genotype <- match(alleles_in_genotype_char, locus_unique_alleles) - 1L
        }
        gen_numeric_matrix[i, ((j-1)*ploidy + 1):(j*ploidy)] <- alleles_in_genotype
      }
    }
  }
  genind_std@tab <- gen_numeric_matrix
  
  # Ensure population information is present, as required by dotR methods.
  if (is.null(adegenet::pop(genind_std))) {
      adegenet::pop(genind_std) <- as.factor(rep("Pop1", n_ind))
  }
  
  return(genind_std)
}

#' Read a VCF file and convert it to a dotR genind object
#' (Ported from dnadot-beta.R)
#' @export
read_vcf_to_genind_dotR <- function(vcf_file, ploidy = 2) {
  if (!requireNamespace("vcfR", quietly = TRUE)) stop("Package 'vcfR' is required.")
  if (!file.exists(vcf_file)) stop("VCF file not found: ", vcf_file)
  
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
  if (nrow(vcf@fix) == 0) stop("VCF file has no variant data.")
  
  gt <- vcfR::extract.gt(vcf, element = "GT", as.character = TRUE)
  if (is.null(gt)) stop("No GT field in VCF.")
  
  gt_t <- t(gt)
  loc_names <- vcf@fix[, "ID"]
  na_ids <- is.na(loc_names) | loc_names == "."
  if (any(na_ids) || any(duplicated(loc_names[!na_ids]))) {
      loc_names <- make.unique(paste0(vcf@fix[, "CHROM"], "_", vcf@fix[, "POS"]), sep = "_")
  }
  colnames(gt_t) <- loc_names
  
  df <- as.data.frame(matrix(NA_character_, nrow = nrow(gt_t), ncol = ncol(gt_t), dimnames = dimnames(gt_t)))
  ref <- vcf@fix[, "REF"]; alt <- vcf@fix[, "ALT"]
  
  for (j in 1:ncol(gt_t)) {
      alleles_map <- c(ref[j], strsplit(alt[j], ",")[[1]])
      for (i in 1:nrow(gt_t)) {
          if (is.na(gt_t[i, j]) || grepl("\\.", gt_t[i, j])) next
          indices <- as.integer(strsplit(gt_t[i, j], "[/|]")[[1]]) + 1
          df[i, j] <- paste(sort(alleles_map[indices]), collapse = "")
      }
  }
  df2genind_dotR(df, sep = "", ploidy = ploidy, type = "codom")
}


# ============================================================================
# --- Main User-Facing Function  ---
# ============================================================================

#' Estimate Census Size (Nc) from SNP data
#'
#' @param method The estimation method to use. Options are:
#'   "bhg_3d" (Beta-Hypergeometric 3D Regularized - Recommended for structured populations),
#'   "wasserstein", "sherwin", "nll", "bhd", "betabin", "betabin_dynamic", and "mle".
#' @export
dnadot_snp <- function(genind_object, method = "wasserstein", ...) {
    # --- Input Validation ---
    if (!inherits(genind_object, "genind")) stop("Input must be a 'genind' object.")
    if (adegenet::nInd(genind_object) == 0) {
        warning("Input genind_object has 0 individuals. Returning NA.")
        return(list(N_est = NA_real_, details = data.frame()))
    }

    args <- list(genind_object = genind_object, ...)
    
    # --- Method Dispatch (Merged Logic) ---
    if (method == "bhg_3d") { # NEW: BHG 3D Dispatch
        return(do.call(.dnadot_bhg_3d_internal, args))
    } else if (method == "betabin") {
        return(do.call(.dnadot_betabin_internal, args))
    } else if (method == "betabin_dynamic") {
        return(do.call(.dnadot_betabin_dynamic_internal, args))
    } else if (method == "bhd") {
        return(do.call(.dnadot_bhd_internal, args))
    } else if (method == "mle") {
        return(do.call(.dnadot_mle_internal, args))
    } else if (method %in% c("wasserstein", "sherwin", "nll")) {
        # Logic ported from dnadot-beta.R
        discrepancy_map <- c("wasserstein" = "wasserstein", 
                               "sherwin" = "absolute_diff", 
                               "nll" = "neg_log_likelihood")
        args$discrepancy_type <- discrepancy_map[method]
        
        if ("sherwin_reg_strength" %in% names(args)) args$reg_strength <- args$sherwin_reg_strength
        if ("sherwin_reg_N_threshold" %in% names(args)) args$reg_N_threshold <- args$sherwin_reg_N_threshold
        
        return(do.call(.dnadot_mde_internal, args))
    } else {
        stop("Invalid method. Must be 'bhg_3d', 'wasserstein', 'betabin', 'betabin_dynamic', 'bhd', 'sherwin', 'nll', or 'mle'.")
    }
}

# ============================================================================
# --- Internal Helper and Setup Functions (Ported and Merged) ---
# ============================================================================

#' Manually Calculate Fis from a dotR genind object
#' (Ported from dnadot-beta.R)
.calculate_fis_manually <- function(dotR_genind) {
    n_ind <- adegenet::nInd(dotR_genind)
    n_loc <- adegenet::nLoc(dotR_genind)
    ploidy <- adegenet::ploidy(dotR_genind)[1]
    
    if (ploidy != 2) {
        stop("Manual Fis calculation currently only supports diploid data.", call. = FALSE)
    }
    
    fis_per_locus <- numeric(n_loc)
    
    for (l in 1:n_loc) {
        start_col <- (l - 1) * ploidy + 1
        end_col <- l * ploidy
        locus_data <- dotR_genind@tab[, start_col:end_col, drop = FALSE]
        
        n_typed_ind <- 0; n_het <- 0
        for (i in 1:n_ind) {
            alleles <- locus_data[i, ]
            if (any(is.na(alleles))) next
            n_typed_ind <- n_typed_ind + 1
            if (alleles[1] != alleles[2]) n_het <- n_het + 1
        }
        hobs <- if (n_typed_ind > 0) n_het / n_typed_ind else 0
        
        all_alleles_vector <- as.vector(locus_data)
        all_alleles_vector <- all_alleles_vector[!is.na(all_alleles_vector)]
        
        if (length(all_alleles_vector) == 0) {
            fis_per_locus[l] <- NA; next
        }
        
        allele_counts <- table(all_alleles_vector)
        allele_freqs <- allele_counts / sum(allele_counts)
        hexp <- 1 - sum(allele_freqs^2)
        
        fis_per_locus[l] <- if (hexp > 0) 1 - (hobs / hexp) else 0
    }
    return(fis_per_locus)
}

#' Set up the parallel backend for analysis.
#' (Ported from dnadot-beta.R, slightly improved)
.setup_parallel <- function(n_cores, n_items) {
    run_parallel <- (n_cores > 1 && n_items > 1)
    if (run_parallel) {
        if (!requireNamespace("doParallel", quietly = TRUE) || !requireNamespace("foreach", quietly = TRUE)) {
            warning("doParallel or foreach not found. Running sequentially.")
            run_parallel <- FALSE
        } else {
            # Check available cores and use the minimum of requested and available
            available_cores <- parallel::detectCores()
            cores_to_use <- min(n_cores, available_cores)
            if (cores_to_use < n_cores) {
                message(sprintf("Note: Requested %d cores, using available %d cores.", n_cores, cores_to_use))
            }
            doParallel::registerDoParallel(cores = cores_to_use)
        }
    }
    if (!run_parallel) foreach::registerDoSEQ()
    return(run_parallel)
}

#' Set up the grid search parameters (2D for existing methods).
#' (Ported from dnadot-beta.R)
.setup_grid_params <- function(n_individuals, jackknife_proportion, num_N_hypothesized, N_try_min, N_try_max, num_p_hypothesized, p_hypothesized_range) {
    n_subsample_indiv <- floor(jackknife_proportion * n_individuals)
    min_N_val <- max(1, n_subsample_indiv)
    if (jackknife_proportion == 1) min_N_val <- max(min_N_val, n_individuals)
    
    N_try_min_eff <- N_try_min %||% max(min_N_val, floor(n_individuals * 0.5), 2)
    # Increased default max N slightly for better coverage
    N_try_max_eff <- N_try_max %||% floor(n_individuals * 5) 
    
    if (N_try_max_eff <= N_try_min_eff) {
        N_try_max_eff <- N_try_min_eff + num_N_hypothesized * 2
        if (!is.null(N_try_min) || !is.null(N_try_max)) {
             warning("N_try_max was less than or equal to N_try_min. Adjusting N_try_max.")
        }
    }
    
    N_try_values <- unique(round(seq(N_try_min_eff, N_try_max_eff, length.out = num_N_hypothesized)))
    p_try_values <- seq(p_hypothesized_range[1], p_hypothesized_range[2], length.out = num_p_hypothesized)
    
    list(N_try = N_try_values, p_try = p_try_values, n_subsample = n_subsample_indiv)
}

#' Set up the 3D grid search parameters (N, p, Theta).
#' (From BHG_dnadot_unfinished.R, relies on .setup_grid_params)
.setup_grid_params_3d <- function(n_individuals, jackknife_proportion, 
                                  num_N_hypothesized, N_try_min, N_try_max, 
                                  num_p_hypothesized, p_hypothesized_range,
                                  num_theta_hypothesized, theta_hypothesized_range) {
    
    # Use the existing 2D setup function to get N and p grids
    grid_2d <- .setup_grid_params(n_individuals, jackknife_proportion, num_N_hypothesized, N_try_min, N_try_max, num_p_hypothesized, p_hypothesized_range)
    
    # Define the theta grid
    # Ensure theta is strictly between 0 and 1 for numerical stability.
    min_theta <- max(1e-6, theta_hypothesized_range[1])
    max_theta <- min(1.0 - 1e-6, theta_hypothesized_range[2])
    
    if (min_theta >= max_theta) {
        stop("Invalid theta hypothesized range. Min must be less than Max, and both must be strictly between 0 and 1.")
    }
    
    theta_try_values <- seq(min_theta, max_theta, length.out = num_theta_hypothesized)
    
    list(N_try = grid_2d$N_try, 
         p_try = grid_2d$p_try, 
         theta_try = theta_try_values,
         n_subsample = grid_2d$n_subsample)
}


#' Process and summarize results from all loci.
#' (Ported from dnadot-beta.R)
.process_results <- function(all_loci_results, method_name) {
    valid_results <- Filter(function(x) is.list(x) && !inherits(x, "error"), all_loci_results)
    
    if (length(valid_results) == 0) {
        warning("No valid results were obtained from any locus.", call. = FALSE)
        return(list(N_est = NA_real_, selected_loci_details = data.frame(), all_runs_summary = data.frame(), method = method_name))
    }
    
    results_df <- do.call(rbind, lapply(valid_results, as.data.frame))
    results_df_filtered <- results_df[is.finite(results_df$cv_subsample_counts) & is.finite(results_df$N_est), ]
    
    N_est_final <- NA_real_
    selected_df <- data.frame()
    
    if (nrow(results_df_filtered) > 0) {
      # Select the top 10% of runs with the lowest CV.
      num_to_select <- max(1, floor(0.10 * nrow(results_df_filtered)))
      selected_df <- head(results_df_filtered[order(results_df_filtered$cv_subsample_counts), ], num_to_select)
      
      if (nrow(selected_df) > 0) {
        N_est_final <- mean(selected_df$N_est, na.rm = TRUE)
      }
    }
    
    return(list(N_est = N_est_final, selected_loci_details = selected_df, all_runs_summary = results_df, method = method_name))
}

# ============================================================================
# --- Internal Method-Specific Handlers (Merged) ---
# ============================================================================

# ----------------------------------------------------------------------------
# --- NEW BHG 3D Method (From BHG_dnadot_unfinished.R) ---
# ----------------------------------------------------------------------------

#' Internal handler for the 3D Beta-Hypergeometric (BHG) method (NEW).
#' @private
.dnadot_bhg_3d_internal <- function(genind_object, jackknife_proportion = 0.8,
                                    num_N_hypothesized = 15, N_try_min = NULL, N_try_max = NULL,
                                    num_p_hypothesized = 25, p_hypothesized_range = c(0.01, 0.99),
                                    num_theta_hypothesized = 20, 
                                    theta_hypothesized_range = c(0.001, 0.5),
                                    use_regularization = TRUE,
                                    n_cores = 1, ...) {
    
    method_name <- "MDE (BHG 3D Regularized)"
    message("Running ", method_name, " method...")

    # --- STEP 1: Calculate Empirical Theta (Fis) for Regularization ---
    
    # Initialize defaults (Infinite variance means zero penalty)
    theta_global_hat <- 0.0
    theta_global_hat_model <- 1e-6 # Default fallback if estimate is negative
    theta_sigma_sq <- Inf 

    if (use_regularization) {
        fis_per_locus <- .calculate_fis_manually(genind_object)
        valid_fis <- fis_per_locus[is.finite(fis_per_locus)]
        
        if (length(valid_fis) > 1) {
            # Calculate the mean (theta_hat) on the raw data
            theta_global_hat <- mean(valid_fis, na.rm = TRUE)
            
            # Constrain the model input. BHG requires theta > 0.
            theta_global_hat_model <- max(1e-6, theta_global_hat)

            # Calculate the uncertainty (Squared Standard Error of the Mean)
            # Sigma^2 = SE^2 = (SD / sqrt(n_loci))^2
            sd_fis <- sd(valid_fis, na.rm = TRUE)
            
            if (is.na(sd_fis) || sd_fis == 0) {
                # If SD is 0, we have high confidence. Set a very small variance (1e-12).
                theta_sigma_sq <- 1e-12 
            } else {
                theta_sigma_sq <- (sd_fis^2) / length(valid_fis)
            }
            message(sprintf("Empirical Theta (Fis) estimate: %.4f (Model Anchor: %.4f)", theta_global_hat, theta_global_hat_model))
            message(sprintf("Theta Uncertainty (Sigma^2): %.6g", theta_sigma_sq))

        } else {
            message("Note: Insufficient data (<2 valid loci) to calculate robust empirical Fis (theta). Running without strong regularization.")
        }
    } else {
        message("Regularization disabled by user request.")
    }

    # --- STEP 2: Setup Grid and Parameters ---
    n_individuals <- adegenet::nInd(genind_object)
    ploidy_val <- unique(adegenet::ploidy(genind_object))
    if (length(ploidy_val) > 1) stop("Mixed ploidy not currently supported.")
    n_loci <- adegenet::nLoc(genind_object)
    
    tab_matrix <- genind_object@tab
    loc_names_vec <- adegenet::locNames(genind_object)
    alleles_list <- adegenet::alleles(genind_object)

    # Define the 3D grid search space.
    grid_params <- .setup_grid_params_3d(n_individuals, jackknife_proportion, 
                                         num_N_hypothesized, N_try_min, N_try_max, 
                                         num_p_hypothesized, p_hypothesized_range,
                                         num_theta_hypothesized, theta_hypothesized_range)

    # Configure parallel processing.
    run_parallel <- .setup_parallel(n_cores, n_loci)
    if (run_parallel) on.exit(foreach::registerDoSEQ(), add = TRUE)

    # --- STEP 3: Main Loop (3D MDE) ---
    `%op%` <- if (run_parallel) foreach::`%dopar%` else foreach::`%do%`
    all_loci_results <- foreach::foreach(locus_idx = 1:n_loci, .combine = 'c', .errorhandling = 'pass') %op% {
        locus_alleles <- alleles_list[[locus_idx]]
        if (length(locus_alleles) < 2) return(NULL) 
        
        locus_results_list <- list()
        start_col <- (locus_idx - 1) * ploidy_val + 1
        end_col <- locus_idx * ploidy_val
        
        # Basic check for index issues
        if (start_col > ncol(tab_matrix) || end_col > ncol(tab_matrix)) return(NULL)

        genotypes_locus <- tab_matrix[, start_col:end_col, drop = FALSE]
        
        for (allele_idx_zero_based in 0:(length(locus_alleles) - 1)) {
            if(!is.integer(genotypes_locus)) storage.mode(genotypes_locus) <- "integer"
            
            # Call the new C++ function for the 3D BHG MDE
            cpp_result <- dotR_dnadot_bhg_3d_cpp(
                genotypes_locus, n_individuals, grid_params$n_subsample,
                grid_params$N_try, grid_params$p_try, grid_params$theta_try,
                allele_idx_zero_based,
                theta_global_hat_model, theta_sigma_sq # Pass regularization parameters
            )
            
            # --- Post-processing ---
            counts <- cpp_result$subsample_target_allele_counts
            mean_counts <- mean(counts, na.rm = TRUE)
            cv_val <- if (length(counts) > 1 && mean_counts > 0) sd(counts, na.rm = TRUE) / mean_counts else NA_real_
            
            allele_char <- locus_alleles[allele_idx_zero_based + 1]
            
            result_for_allele <- list(locus_name = loc_names_vec[locus_idx],
                                      target_allele = allele_char, 
                                      N_est = cpp_result$N_est, 
                                      p_est = cpp_result$p_est,
                                      theta_est = cpp_result$theta_est, # NEW
                                      min_discrepancy = cpp_result$min_discrepancy, 
                                      cv_subsample_counts = cv_val)
            
            locus_results_list[[length(locus_results_list) + 1]] <- result_for_allele
        }
        locus_results_list
    }

    # --- STEP 4: Finalization ---
    return(.process_results(all_loci_results, method_name))
}

# ----------------------------------------------------------------------------
# --- Existing Methods (Ported from dnadot-beta.R) ---
# ----------------------------------------------------------------------------

#' Internal handler for standard MDE methods (Wasserstein, Sherwin, NLL).
#' (Ported from dnadot-beta.R)
#' @private
.dnadot_mde_internal <- function(genind_object, discrepancy_type, jackknife_proportion = 0.8,
                                 num_N_hypothesized = 11, N_try_min = NULL, N_try_max = NULL,
                                 num_p_hypothesized = 20, p_hypothesized_range = c(0.01, 0.99), 
                                 reg_strength = 0.0, reg_N_threshold = Inf, n_cores = 1, ...) {
    
    method_name <- paste0("MDE (", discrepancy_type, ")")
    
    # Add a flag to prevent recursive messaging when called by betabin_dynamic
    if (!("internal_call" %in% names(list(...))) || !list(...)$internal_call) {
      message("Running ", method_name, " method...")
    }

    # --- Setup ---
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
        
        if (start_col > ncol(tab_matrix) || end_col > ncol(tab_matrix)) return(NULL)

        genotypes_locus <- tab_matrix[, start_col:end_col, drop = FALSE]
        
        for (allele_idx_zero_based in 0:(length(locus_alleles) - 1)) {
            if(!is.integer(genotypes_locus)) storage.mode(genotypes_locus) <- "integer"
            
            cpp_result <- dotR_dnadot_mde_cpp(genotypes_locus, n_individuals, grid_params$n_subsample,
                                              grid_params$N_try, grid_params$p_try, allele_idx_zero_based,
                                              discrepancy_type, reg_strength, reg_N_threshold,
                                              debug = FALSE)
            
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

# Still to commit
# ... [Insert implementations of .dnadot_bhd_internal, .dnadot_betabin_internal, .dnadot_betabin_dynamic_internal, .dnadot_mle_internal here] ...
