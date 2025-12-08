# R SCRIPT FOR SIMULATING, TUNING, AND COMPARING ALL dotR METHODS
# This script performs a structured parameter sweep for both the simulation
# conditions (e.g., Nc, Fis) and the estimation method parameters (e.g., grid density)
# to find optimal settings and compare method performance.

# --- 0. Setup and Package Loading ---
message("--- dotR Simulation, Tuning & Comparison Script ---")

# Ensure all necessary packages are installed and loaded
packages_needed <- c("ggplot2", "tidyr", "dplyr", "scales", "adegenet", "foreach", "doParallel")
lapply(packages_needed, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
})

# Load the dotR package itself
# If developing locally, use devtools::load_all() beforehand.
if (!requireNamespace("dotR", quietly = TRUE)) {
    stop("dotR package not found. Please install it before running simulations.")
}
library(dotR)

# --- 1. Simulation Configuration ---

# A. SIMULATION SCENARIOS: Defines the "true" populations we will simulate.
param_sweep <- expand.grid(
  true_Nc = c(150, 600),             # True census sizes to test
  target_fis = c(0.001, 0.05),         # Levels of inbreeding (0.0 = panmixia)
  n_sample_prop = c(0.1,0.2,0.3,0.5,0.7),           # Proportion of Nc to sample
  n_loci = c(1000,5000),                  # Number of loci
  stringsAsFactors = FALSE
)

# B. TUNING PARAMETERS: Defines the dnadot_snp parameters we want to test.
tuning_grid <- expand.grid(
    num_N_hypothesized = c(15, 51), # Test a coarse vs. a fine grid for Nc search
    N_try_max_factor = c(3, 6),     # Test a narrow vs. a wide search range for Nc (multiplier of n_sample)
    stringsAsFactors = FALSE
)

# C. METHODS to compare
methods_to_run <- c("wasserstein", "bhd", "betabin", "betabin_dynamic")

# D. GLOBAL settings
ploidy_sim <- 2
num_simulation_replicates <- 20 # Increase for more robust results (e.g., 50-100)
default_n_cores_dnadot <- 1     # Cores for the C++ part (best to keep at 1)

# E. SETUP PARALLEL BACKEND for the main simulation loop
n_cores_loop <- parallel::detectCores() - 1
if (is.na(n_cores_loop) || n_cores_loop < 1) n_cores_loop <- 1
message(paste("Registering", n_cores_loop, "cores for the main simulation loop."))
doParallel::registerDoParallel(cores = n_cores_loop)


# --- 2. Structured Genotype Generation Function ---
generate_structured_population <- function(true_Nc, n_loci, ploidy, target_fis) {
    if (ploidy != 2) stop("This simulation is designed for diploid organisms.")
    n_family_groups <- 25
    global_p_freqs <- runif(n_loci, 0.1, 0.9)
    if (target_fis > 0) {
        alpha_vals <- global_p_freqs * (1 - target_fis) / target_fis
        beta_vals <- (1 - global_p_freqs) * (1 - target_fis) / target_fis
    }
    family_p_freqs <- matrix(NA, nrow = n_family_groups, ncol = n_loci)
    for (l in 1:n_loci) {
        if (target_fis > 0) {
            family_p_freqs[, l] <- rbeta(n_family_groups, alpha_vals[l], beta_vals[l])
        } else {
            family_p_freqs[, l] <- rep(global_p_freqs[l], n_family_groups)
        }
    }
    Nc_per_group <- floor(true_Nc / n_family_groups)
    full_pop_df_list <- list()
    for (fam_group in 1:n_family_groups) {
        group_df <- data.frame(matrix(NA_character_, nrow = Nc_per_group, ncol = n_loci))
        colnames(group_df) <- paste0("Locus", 1:n_loci)
        for (l in 1:n_loci) {
            p <- family_p_freqs[fam_group, l]
            probs <- c(p^2, 2*p*(1-p), (1-p)^2)
            genos <- c("AA", "AT", "TT")
            group_df[, l] <- sample(genos, Nc_per_group, replace = TRUE, prob = probs)
        }
        full_pop_df_list[[fam_group]] <- group_df
    }
    full_pop_df <- do.call(rbind, full_pop_df_list)
    return(full_pop_df)
}

# --- 3. Main Simulation Loop ---
message("\n--- Starting Main Simulation & Tuning Loop ---")

# Create a combined grid of all scenarios to test
full_simulation_grid <- tidyr::crossing(param_sweep, tuning_grid)

results_list <- foreach(
    scenario_idx = 1:nrow(full_simulation_grid),
    .combine = 'rbind',
    .packages = c("adegenet", "dotR")
) %dopar% {
    
    params <- full_simulation_grid[scenario_idx, ]
    replicate_results_df <- data.frame()

    for (rep_id in 1:num_simulation_replicates) {
        set.seed(scenario_idx * 1000 + rep_id)
        true_pop <- generate_structured_population(
          true_Nc = params$true_Nc, n_loci = params$n_loci, 
          ploidy = ploidy_sim, target_fis = params$target_fis
        )
        n_sample <- floor(params$true_Nc * params$n_sample_prop)
        if (n_sample < 20) next
        
        sampled_df <- true_pop[sample(1:nrow(true_pop), size = n_sample), , drop = FALSE]
        sim_genind <- try(df2genind_dotR(df = sampled_df, sep = "", ploidy = ploidy_sim), silent = TRUE)
        if (inherits(sim_genind, "try-error")) next

        for (method_name in methods_to_run) {
            # Define arguments, now including the tuning parameters
            args <- list(
                genind_object = sim_genind, 
                method = method_name, 
                n_cores = default_n_cores_dnadot,
                num_N_hypothesized = params$num_N_hypothesized,
                N_try_max = floor(n_sample * params$N_try_max_factor) # Set max search range based on sample size
            )
            
            res <- try(do.call(dnadot_snp, args), silent = TRUE)
            N_estimate <- if (!inherits(res, "try-error")) res$N_est else NA
            
            replicate_results_df <- rbind(replicate_results_df, data.frame(
                true_Nc = params$true_Nc,
                target_fis = params$target_fis,
                n_sample = n_sample,
                replicate = rep_id,
                method = method_name,
                num_N_hypothesized = params$num_N_hypothesized,
                N_try_max_factor = params$N_try_max_factor,
                N_est = N_estimate
            ))
        }
    }
    replicate_results_df
}

doParallel::stopImplicitCluster()

if (is.null(results_list) || nrow(results_list) == 0) {
    stop("No results were generated from the simulation.")
}

results_list$bias <- (results_list$N_est - results_list$true_Nc) / results_list$true_Nc

# --- 4. Analysis and Visualization ---
message("\n--- Simulation Complete. Analyzing and plotting results. ---")

summary_stats <- results_list %>%
    group_by(true_Nc, target_fis, method, num_N_hypothesized, N_try_max_factor) %>%
    summarise(
        mean_bias = mean(bias, na.rm = TRUE),
        median_bias = median(bias, na.rm = TRUE),
        rmse = sqrt(mean((N_est - first(true_Nc))^2, na.rm = TRUE)),
        cv = sd(N_est, na.rm = TRUE) / mean(N_est, na.rm = TRUE),
        .groups = 'drop'
    )

# --- 5. Saving Results ---
output_dir <- paste0("tuning_and_comparison_results_", format(Sys.time(), "%Y%m%d_%H%M"))
if (!dir.exists(output_dir)) dir.create(output_dir)
write.csv(results_list, file.path(output_dir, "full_raw_results.csv"), row.names = FALSE)
write.csv(summary_stats, file.path(output_dir, "summary_stats.csv"), row.names = FALSE)
message(paste("\nResults saved to", output_dir))

# --- 6. Plotting ---

# Prepare labels for faceting
summary_stats$Nc_Label <- factor(paste("True Nc:", summary_stats$true_Nc))
summary_stats$Fis_Label <- factor(paste("Target Fis:", summary_stats$target_fis))

# PLOT 1: How do tuning parameters affect median bias?
plot_tuning_bias <- ggplot(summary_stats, 
       aes(x = factor(N_try_max_factor), y = median_bias, 
           color = factor(num_N_hypothesized), group = factor(num_N_hypothesized))) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    facet_grid(method ~ Nc_Label + Fis_Label, scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Effect of Tuning Parameters on Median Bias",
         subtitle = "Comparing coarse (e.g., 15) vs. fine (e.g., 51) Nc search grids across different search ranges",
         x = "Nc Search Range (Multiplier of Sample Size)",
         y = "Median Relative Bias",
         color = "Grid Density\n(num_N_hypothesized)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plot_tuning_bias.png"), plot_tuning_bias, width = 14, height = 10)


# PLOT 2: How do tuning parameters affect precision (CV)?
plot_tuning_cv <- ggplot(summary_stats, 
       aes(x = factor(N_try_max_factor), y = cv, 
           color = factor(num_N_hypothesized), group = factor(num_N_hypothesized))) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    facet_grid(method ~ Nc_Label + Fis_Label, scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Effect of Tuning Parameters on Precision (CV)",
         subtitle = "Lower is better. Comparing coarse vs. fine Nc search grids.",
         x = "Nc Search Range (Multiplier of Sample Size)",
         y = "Coefficient of Variation (CV)",
         color = "Grid Density\n(num_N_hypothesized)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plot_tuning_cv.png"), plot_tuning_cv, width = 14, height = 10)


message("All done. Check the '", output_dir, "' directory for results and plots.")
