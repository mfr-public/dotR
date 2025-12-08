# R SCRIPT FOR SIMULATING, TUNING, AND COMPARING dotR METHODS
# This script expands on the original to perform a structured parameter sweep
# for both the simulation conditions and the estimation method parameters.

# --- 0. Setup and Package Loading ---
message("--- dotR Simulation & Tuning Script ---")

# Ensure all necessary packages are installed and loaded
packages_needed <- c("ggplot2", "tidyr", "dplyr", "scales", "adegenet", "foreach", "doParallel")
lapply(packages_needed, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
})

# Load the dotR package itself
# If developing locally, use devtools::load_all() beforehand.
if (!requireNamespace("dotRdebug", quietly = TRUE)) {
    stop("dotRdebug package not found. Please install it before running simulations.")
}
library(dotRdebug)

# Helper function
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- 1. Simulation Configuration ---

# A. SIMULATION PARAMETERS: Defines the "true" populations we will simulate.
param_sweep <- expand.grid(
  true_Nc = c(100, 500),                # True census sizes
  n_sample_prop = c(0.3, 0.7),          # Proportion of Nc to sample
  n_loci = c(250),                      # Number of loci
  AFS = c("uniform", "beta_L_shaped"),  # Allele Frequency Spectrum
  stringsAsFactors = FALSE
)

# B. TUNING PARAMETERS: Defines the dnadot_snp parameters we want to test.
tuning_grid <- expand.grid(
    num_N_hypothesized = c(15, 51), # Test a coarse vs. a fine grid for Nc search
    N_try_max_factor = c(3, 6),     # Test a narrow vs. a wide search range for Nc (multiplier of n_sample)
    stringsAsFactors = FALSE
)

# C. METHODS to compare
methods_to_run <- c("wasserstein", "bhd", "sherwin", "nll") # Excluding MLE for now

# D. GLOBAL settings
ploidy_sim <- 2
num_simulation_replicates <- 25 # Increase for more robust results (e.g., 100+)
default_n_cores_dnadot <- 1     # Cores for the C++ part (keep at 1)

# E. SETUP PARALLEL BACKEND for the main simulation loop
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

# --- 3. Main Simulation Loop (Now iterates over tuning grid as well) ---

# We create a combined grid of all scenarios to test
full_simulation_grid <- tidyr::crossing(param_sweep, tuning_grid)

# Use foreach for efficient parallel execution
results_list <- foreach(
    scenario_idx = 1:nrow(full_simulation_grid),
    .combine = 'rbind',
    .packages = c("adegenet", "dotRdebug")
) %dopar% {
    
    # Extract current simulation and tuning parameters
    params <- full_simulation_grid[scenario_idx, ]
    
    replicate_results_df <- data.frame()

    for (rep_id in 1:num_simulation_replicates) {
        # A. Generate Data
        set.seed(scenario_idx * 1000 + rep_id) # Ensure reproducibility
        
        true_pop <- generate_true_population_genotypes(
          true_Nc = params$true_Nc, n_loci = params$n_loci, 
          ploidy = ploidy_sim, AFS_shape = params$AFS
        )
        n_sample <- floor(params$true_Nc * params$n_sample_prop)
        
        if (n_sample < 10) next # Skip if sample size is too low
        
        sampled_df <- true_pop[sample(1:params$true_Nc, size = n_sample), , drop = FALSE]
        
        sim_genind <- try(
          df2genind_dotR(df = sampled_df, sep = "", ploidy = ploidy_sim),
          silent = TRUE
        )
        if (inherits(sim_genind, "try-error")) next

        # B. Run all methods for the current tuning setting
        for (method_name in methods_to_run) {
            
            # Define arguments, now including the tuning parameters
            args <- list(
                genind_object = sim_genind, 
                method = method_name, 
                n_cores = default_n_cores_dnadot,
                num_N_hypothesized = params$num_N_hypothesized,
                N_try_max = floor(n_sample * params$N_try_max_factor)
            )
            
            res <- try(do.call(dnadot_snp, args), silent = TRUE)
            
            N_estimate <- if (!inherits(res, "try-error")) res$N_est else NA
            
            # C. Store results for this single run
            replicate_results_df <- rbind(replicate_results_df, data.frame(
                true_Nc = params$true_Nc,
                sample_prop = params$n_sample_prop,
                n_sample = n_sample,
                AFS = params$AFS,
                replicate = rep_id,
                method = method_name,
                num_N_hypothesized = params$num_N_hypothesized,
                N_try_max_factor = params$N_try_max_factor,
                N_est = N_estimate
            ))
        }
    }
    replicate_results_df # Return all results for this scenario
}

# Stop the parallel backend
doParallel::stopImplicitCluster()

if (is.null(results_list) || nrow(results_list) == 0) {
    stop("No results were generated from the simulation.")
}

# Calculate relative bias
results_list$bias <- (results_list$N_est - results_list$true_Nc) / results_list$true_Nc

# --- 4. Analysis and Visualization ---
message("\n\n--- Simulation Complete. Analyzing and plotting results. ---")

# Calculate summary statistics, now including tuning parameters
summary_stats <- results_list %>%
    group_by(true_Nc, sample_prop, AFS, method, num_N_hypothesized, N_try_max_factor) %>%
    summarise(
        mean_bias = mean(bias, na.rm = TRUE),
        median_bias = median(bias, na.rm = TRUE),
        rmse = sqrt(mean((N_est - first(true_Nc))^2, na.rm = TRUE)),
        cv = sd(N_est, na.rm = TRUE) / mean(N_est, na.rm = TRUE),
        .groups = 'drop'
    )

# --- 5. Saving Results ---
output_dir <- paste0("tuning_results_", format(Sys.time(), "%Y%m%d_%H%M"))
if (!dir.exists(output_dir)) dir.create(output_dir)
write.csv(results_list, file.path(output_dir, "tuning_raw_results.csv"), row.names = FALSE)
write.csv(summary_stats, file.path(output_dir, "tuning_summary_stats.csv"), row.names = FALSE)
message(paste("\nResults saved to", output_dir))

# --- 6. Plotting ---

# Prepare labels for faceting
summary_stats$Nc_Label <- factor(paste("True Nc:", summary_stats$true_Nc))
summary_stats$AFS_Label <- factor(paste("AFS:", summary_stats$AFS))
summary_stats$Prop_Label <- factor(paste("Sample Prop:", summary_stats$sample_prop))

# PLOT 1: How does grid density (num_N_hypothesized) affect bias?
plot_tuning_grid_density <- ggplot(summary_stats, 
       aes(x = factor(N_try_max_factor), y = median_bias, 
           color = factor(num_N_hypothesized), group = factor(num_N_hypothesized))) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    facet_grid(method ~ Nc_Label + AFS_Label + Prop_Label) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Effect of Tuning Parameters on Median Bias",
         subtitle = "Comparing coarse (15) vs. fine (51) Nc search grids across different search ranges",
         x = "Nc Search Range (Multiplier of Sample Size)",
         y = "Median Relative Bias",
         color = "Grid Density\n(num_N_hypothesized)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plot_tuning_bias.png"), plot_tuning_grid_density, width = 16, height = 10)


# PLOT 2: How does grid density affect precision (CV)?
plot_tuning_cv <- ggplot(summary_stats, 
       aes(x = factor(N_try_max_factor), y = cv, 
           color = factor(num_N_hypothesized), group = factor(num_N_hypothesized))) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    facet_grid(method ~ Nc_Label + AFS_Label + Prop_Label) +
    labs(title = "Effect of Tuning Parameters on Precision (CV)",
         subtitle = "Lower is better. Comparing coarse vs. fine Nc search grids.",
         x = "Nc Search Range (Multiplier of Sample Size)",
         y = "Coefficient of Variation (CV)",
         color = "Grid Density\n(num_N_hypothesized)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plot_tuning_cv.png"), plot_tuning_cv, width = 16, height = 10)


message("All done. Check the '", output_dir, "' directory for results and plots.")
