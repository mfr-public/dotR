# R SCRIPT TO ANALYZE AND PLOT REPAIRED SIMULATION RESULTS
# This script reads the salvaged simulation data, calculates summary statistics,
# and generates plots to evaluate tuning parameters and compare method performance.

message("--- dotR Simulation Analysis & Plotting Script ---")

# --- 0. Setup ---
# Load necessary libraries for data manipulation and plotting
packages_needed <- c("ggplot2", "tidyr", "dplyr", "scales")
lapply(packages_needed, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

# --- 1. Load Repaired Data ---
# Define the path to your repaired results file.
repaired_results_path <- "repaired_full_raw_results.csv"

if (!file.exists(repaired_results_path)) {
  stop("Repaired results file not found: ", repaired_results_path,
       "\nPlease make sure the file is in your working directory or provide the full path.")
}

message("Loading repaired data from: ", repaired_results_path)
results_list <- read.csv(repaired_results_path)
message("Loaded ", nrow(results_list), " rows of repaired data.")

# --- 2. Calculate Summary Statistics ---
message("Calculating summary statistics...")

# Calculate bias relative to the true Nc
results_list$bias <- (results_list$N_est - results_list$true_Nc) / results_list$true_Nc

# Group by all simulation and tuning parameters to get detailed summaries.
summary_stats <- results_list %>%
  group_by(true_Nc, target_fis, n_sample_prop, n_loci, method, num_N_hypothesized, N_try_max_factor) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    median_bias = median(bias, na.rm = TRUE),
    rmse = sqrt(mean((N_est - first(true_Nc))^2, na.rm = TRUE)),
    cv = sd(N_est, na.rm = TRUE) / mean(N_est, na.rm = TRUE),
    .groups = 'drop'
  )

# --- 3. Saving Results ---
output_dir <- paste0("final_analysis_plots_", format(Sys.time(), "%Y%m%d_%H%M"))
if (!dir.exists(output_dir)) dir.create(output_dir)
write.csv(summary_stats, file.path(output_dir, "final_summary_stats.csv"), row.names = FALSE)
message(paste("\nAnalysis complete. Results and plots will be saved to:", output_dir))

# --- 4. Plotting ---
message("Generating plots...")

# Prepare clear labels for faceting in the plots
summary_stats$Nc_Label <- factor(paste("True Nc:", summary_stats$true_Nc))
summary_stats$Fis_Label <- factor(paste("Target Fis:", summary_stats$target_fis))
summary_stats$Prop_Label <- factor(paste("Sample Prop:", summary_stats$n_sample_prop))
summary_stats$Loci_Label <- factor(paste("Loci:", summary_stats$n_loci))

# PLOT 1: How do tuning parameters affect median bias?
# **IMPROVED**: Using color for sampling proportion and linetype for grid density for a clearer plot.
plot_tuning_bias <- ggplot(summary_stats,
                           aes(x = factor(N_try_max_factor), y = median_bias,
                               color = factor(n_sample_prop),
                               linetype = factor(num_N_hypothesized),
                               group = interaction(n_sample_prop, num_N_hypothesized))) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  facet_grid(method + Loci_Label ~ Nc_Label + Fis_Label, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Effect of Tuning Parameters on Median Bias",
       subtitle = "Lines show different sampling proportions; Linetype shows grid density (15 vs 51)",
       x = "Nc Search Range (Multiplier of Sample Size)",
       y = "Median Relative Bias",
       color = "Sample Prop.",
       linetype = "Grid Density") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0))

ggsave(file.path(output_dir, "plot_tuning_bias.png"), plot_tuning_bias, width = 14, height = 12, dpi=300)

# PLOT 2: How do tuning parameters affect precision (CV)?
# **IMPROVED**: Using color for sampling proportion and linetype for grid density.
plot_tuning_cv <- ggplot(summary_stats,
                         aes(x = factor(N_try_max_factor), y = cv,
                             color = factor(n_sample_prop),
                             linetype = factor(num_N_hypothesized),
                             group = interaction(n_sample_prop, num_N_hypothesized))) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  facet_grid(method + Loci_Label ~ Nc_Label + Fis_Label, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Effect of Tuning Parameters on Precision (CV)",
       subtitle = "Lower is better. Lines show sampling proportions; Linetype shows grid density.",
       x = "Nc Search Range (Multiplier of Sample Size)",
       y = "Coefficient of Variation (CV)",
       color = "Sample Prop.",
       linetype = "Grid Density") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0))

ggsave(file.path(output_dir, "plot_tuning_cv.png"), plot_tuning_cv, width = 14, height = 12, dpi=300)

# To make the next plots clearer, we filter for the best grid density we identified (51)
summary_stats_filtered <- summary_stats %>%
  filter(num_N_hypothesized == 51)

# PLOT 3: Effect of Number of Loci on Overall Performance (RMSE)
plot_loci_effect <- ggplot(summary_stats_filtered,
                           aes(x = factor(n_loci), y = rmse, color = method,
                               linetype = factor(N_try_max_factor))) +
  geom_line(aes(group = interaction(method, N_try_max_factor)), linewidth = 1.1) + # Group explicitly
  geom_point(size = 3) +
  facet_wrap(vars(Nc_Label, Fis_Label, Prop_Label), scales = "free_y", ncol = 4) +
  labs(
    title = "Effect of Number of Loci on Overall Performance (RMSE)",
    subtitle = "Lower RMSE is better. Comparing 1000 vs. 5000 loci across scenarios.",
    x = "Number of Loci",
    y = "Root Mean Square Error (RMSE)",
    color = "Method",
    linetype = "Nc Search Range\n(x Sample Size)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plot_loci_effect_rmse.png"), plot_loci_effect, width = 16, height = 12, dpi=300)

# PLOT 4: Effect of Sampling Proportion on Accuracy (Median Bias)
plot_sampling_effect <- ggplot(summary_stats_filtered,
                               aes(x = factor(n_sample_prop), y = median_bias, color = method,
                                   linetype = factor(N_try_max_factor))) +
  geom_line(aes(group = interaction(method, N_try_max_factor)), linewidth = 1.1) + # Group explicitly
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  facet_wrap(vars(Nc_Label, Fis_Label, Loci_Label), scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Effect of Sampling Proportion on Accuracy (Median Bias)",
    subtitle = "How performance changes as more of the population is sampled.",
    x = "Proportion of True Nc Sampled",
    y = "Median Relative Bias",
    color = "Method",
    linetype = "Nc Search Range\n(x Sample Size)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plot_sampling_effect_bias.png"), plot_sampling_effect, width = 16, height = 9, dpi=300)

message("\nAll done. Check the '", output_dir, "' directory for your final plots and summary statistics.")

