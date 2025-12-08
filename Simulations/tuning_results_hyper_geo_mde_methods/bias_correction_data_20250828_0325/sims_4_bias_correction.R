# R SCRIPT FOR BIAS CORRECTION CALIBRATION
# This script runs a targeted set of simulations to generate data for
# modeling the systematic bias of the wasserstein estimator at low sampling fractions.

# --- 0. Setup and Package Loading ---
message("--- dotR Bias Correction Simulation Script ---")

# Ensure all necessary packages are installed and loaded
packages_needed <- c("ggplot2", "tidyr", "dplyr", "adegenet", "foreach", "doParallel")
lapply(packages_needed, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
})

# Load the dotR package
if (!requireNamespace("dotRdebug", quietly = TRUE)) {
    stop("dotRdebug package not found. Please install it first.")
}
library(dotRdebug)

# --- 1. Simulation Configuration ---

# A. SIMULATION PARAMETERS: Focus on a wide range of Nc values and low sampling proportions.
param_sweep <- expand.grid(
  true_Nc = c(50, 100, 250, 500, 1000, 2000), # Expanded true census sizes
  n_sample_prop = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),     # Focus on low sampling proportions
  n_loci = c(500),                      # A reasonable number of loci
  AFS = c("uniform", "beta_L_shaped"),  # Test both allele frequency spectrums
  stringsAsFactors = FALSE
)

# B. OPTIMAL ESTIMATOR PARAMETERS (based on previous tuning)
OPTIMAL_METHOD <- "wasserstein"
OPTIMAL_NUM_N_HYPOTHESIZED <- 51
OPTIMAL_N_TRY_MAX_FACTOR <- 3

# C. GLOBAL settings
ploidy_sim <- 2
num_simulation_replicates <- 50 # A good number for stable modeling
default_n_cores_dnadot <- 1     # Cores for the C++ part

# D. SETUP PARALLEL BACKEND
n_cores_loop <- parallel::detectCores() - 1
message(paste("Registering", n_cores_loop, "cores for the main simulation loop."))
doParallel::registerDoParallel(cores = n_cores_loop)


# --- 2. Genotype Generation Function (Unchanged) ---
generate_true_population_genotypes <- function(true_Nc, n_loci, ploidy, AFS_shape) {
  if (ploidy != 2) stop("This genotype generation is designed for ploidy = 2.")
  if (AFS_shape == "beta_L_shaped") {
      locus_target_freqs <- rbeta(n_loci, 1, 5) 
  } else {
      locus_target_freqs <- runif(n_loci, 0.05, 0.95)
  }
  true_pop_df <- data.frame(matrix(NA_character_, nrow = true_Nc, ncol = n_loci))
  colnames(true_pop_df) <- paste0("Locus", 1:n_loci)
  for (locus in 1:n_loci) {
    alleles <- c("A", "T")
    p <- locus_target_freqs[locus]; q <- 1 - p
    probs <- c(p^2, 2*p*q, q^2)
    genos <- c("AA", "AT", "TT")
    true_pop_df[, locus] <- sample(genos, size = true_Nc, replace = TRUE, prob = probs)
  }
  return(true_pop_df)
}

# --- 3. Main Simulation Loop ---

results_list <- foreach(
    scenario_idx = 1:nrow(param_sweep),
    .combine = 'rbind',
    .packages = c("adegenet", "dotRdebug")
) %dopar% {
    
    params <- param_sweep[scenario_idx, ]
    replicate_results_df <- data.frame()

    for (rep_id in 1:num_simulation_replicates) {
        set.seed(scenario_idx * 1000 + rep_id)
        
        true_pop <- generate_true_population_genotypes(
          true_Nc = params$true_Nc, n_loci = params$n_loci, 
          ploidy = ploidy_sim, AFS_shape = params$AFS
        )
        n_sample <- floor(params$true_Nc * params$n_sample_prop)
        
        if (n_sample < 10) next
        
        sampled_df <- true_pop[sample(1:params$true_Nc, size = n_sample), , drop = FALSE]
        
        sim_genind <- try(
          df2genind_dotR(df = sampled_df, sep = "", ploidy = ploidy_sim),
          silent = TRUE
        )
        if (inherits(sim_genind, "try-error")) next

        args <- list(
            genind_object = sim_genind, 
            method = OPTIMAL_METHOD, 
            n_cores = default_n_cores_dnadot,
            num_N_hypothesized = OPTIMAL_NUM_N_HYPOTHESIZED,
            N_try_max = floor(n_sample * OPTIMAL_N_TRY_MAX_FACTOR)
        )
        
        res <- try(do.call(dnadot_snp, args), silent = TRUE)
        
        N_estimate <- if (!inherits(res, "try-error")) res$N_est else NA
        
        replicate_results_df <- rbind(replicate_results_df, data.frame(
            true_Nc = params$true_Nc,
            sample_prop = params$n_sample_prop,
            AFS = params$AFS,
            N_est = N_estimate
        ))
    }
    replicate_results_df
}

# Stop the parallel backend
doParallel::stopImplicitCluster()

if (is.null(results_list) || nrow(results_list) == 0) {
    stop("No results were generated from the simulation.")
}

# --- 4. Saving Results ---
output_dir <- paste0("bias_correction_data_", format(Sys.time(), "%Y%m%d_%H%M"))
if (!dir.exists(output_dir)) dir.create(output_dir)
output_file <- file.path(output_dir, "bias_calibration_data.csv")
write.csv(results_list, output_file, row.names = FALSE)

message(paste("\n\n--- Simulation Complete. ---"))
message(paste("Calibration data saved to:", output_file))

# --- 5. Quick Visualization of the Bias ---
# This plot helps visualize the relationship we need to model.
bias_plot <- ggplot(results_list, aes(x = N_est, y = true_Nc)) +
    geom_point(alpha = 0.3, aes(color = factor(sample_prop))) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
    geom_smooth(method = "loess", se = FALSE, aes(color = factor(sample_prop)), linewidth = 1.2) +
    facet_wrap(~AFS) +
    labs(title = "Systematic Bias: True Nc vs. Estimated Nc",
         subtitle = "The red dashed line is perfect 1:1 accuracy. The curved lines show the systematic underestimation.",
         x = "Estimated Nc (N_est)",
         y = "True Nc",
         color = "Sample Prop.") +
    theme_bw()

ggsave(file.path(output_dir, "bias_relationship_plot.png"), bias_plot, width = 12, height = 7)

message("A plot visualizing the bias has been saved in the results directory.")