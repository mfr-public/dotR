# R SCRIPT FOR TESTING ALL dotR METHODS, INCLUDING BETA-BINOMIAL
# This version simulates a population with "soft" or "cryptic" structure
# (e.g., family groups, non-random mating) to generate realistic overdispersion (Fis > 0).

# --- 0. Setup and Package Loading ---
message("--- Comprehensive Test Script for dnadot ---")

# For development, ensure the local package source is loaded:
# devtools::load_all()

if (!requireNamespace("adegenet", quietly = TRUE)) {
  install.packages("adegenet")
}
library(adegenet)

# Ensure your package is loaded. Update "dotR" to your package name if different.
library(dotR) 

if (!exists("df2genind_dotR") || !exists("dnadot_snp")) {
  stop("dotR functions not found. Ensure the package is loaded.")
}

# --- 1. Define a Scenario with Cryptic Population Structure ---
message("\n--- Defining a structured population test scenario ---")
true_Nc_total <- 150       # Total census size of the population
n_family_groups <- 20      # Number of conceptual family lineages to simulate.
# A higher number models more subtle, continuous structure.
target_fis <- 0.03         # The desired level of overdispersion (Fis).
# A small positive value (0.01-0.05) is typical for cryptic structure.
n_loci_test <- 500         # Increased loci for a stable Fis estimate
n_sample_test <- 75        # Number of individuals to sample from the pooled population
jackknife_prop_test <- 0.8
ploidy_test <- 2

# --- 2. Generate Structured Genetic Data ---
# We simulate multiple "family groups" and then pool them. This mimics a single
# population with non-random mating, creating a Wahlund effect and thus Fis > 0.
message("Generating cryptically structured test data to produce Fis > 0...")
set.seed(123)

# 2a. Define the mean allele frequencies for the entire population
global_p_freqs <- runif(n_loci_test, 0.1, 0.9)

# 2b. Generate allele frequencies for each family group from a Beta distribution.
# The parameters of the Beta distribution are chosen to create the desired `target_fis`.
# This creates subtle variations in p around the global mean for each lineage.
alpha_vals <- global_p_freqs * (1 - target_fis) / target_fis
beta_vals <- (1 - global_p_freqs) * (1 - target_fis) / target_fis

family_p_freqs <- matrix(NA, nrow = n_family_groups, ncol = n_loci_test)
for (l in 1:n_loci_test) {
  # Draw a specific allele frequency for each family group for this locus
  family_p_freqs[, l] <- rbeta(n_family_groups, alpha_vals[l], beta_vals[l])
}

# 2c. Create individuals within each family group (assuming HWE within each lineage)
Nc_per_group <- floor(true_Nc_total / n_family_groups)
full_pop_df_list <- list()

for (fam_group in 1:n_family_groups) {
  group_df <- data.frame(matrix(NA_character_, nrow = Nc_per_group, ncol = n_loci_test))
  colnames(group_df) <- paste0("Locus", 1:n_loci_test)
  
  for (l in 1:n_loci_test) {
    p <- family_p_freqs[fam_group, l]
    probs <- c(p^2, 2*p*(1-p), (1-p)^2) # HWE probabilities for this lineage
    genos <- c("AA", "AT", "TT")
    group_df[, l] <- sample(genos, Nc_per_group, replace = TRUE, prob = probs)
  }
  full_pop_df_list[[fam_group]] <- group_df
}

# 2d. Combine all family groups into one large, intermixed population pool
full_pop_df <- do.call(rbind, full_pop_df_list)

# 2e. Sample from the pooled population
sampled_df_test <- full_pop_df[sample(1:nrow(full_pop_df), n_sample_test), ]

# --- 3. Create genind objects ---
message("\nCreating genind objects...")

# --- FIX: Create a STANDARD genind object FIRST for the diagnostic check ---
# The summary() function requires the standard adegenet @tab format, which our
# custom df2genind_dotR function bypasses for C++ compatibility.
message("  - Creating standard genind for diagnostics...")
genind_standard_for_diag <- adegenet::df2genind(sampled_df_test, sep = "", ploidy = ploidy_test)

# Now, create the custom-formatted genind object required by dnadot's C++ code.
message("  - Creating custom dnadot genind for analysis...")
genind_test <- df2genind_dotR(sampled_df_test, sep = "", ploidy = ploidy_test)
message("Test data and genind objects created successfully.")


# --- 4. Diagnostic Check ---
message("\n--- Diagnostic Check of Simulated Data ---")
# Run the summary on the STANDARD genind object.
gen_summary <- summary(genind_standard_for_diag)
# Calculate Fis per locus using the standard formula: 1 - (Hobs / Hexp)
fis_per_locus <- 1 - (gen_summary$Hobs / gen_summary$Hexp)
# Filter out non-finite values (e.g., from monomorphic loci where Hexp is 0)
observed_fis_values <- fis_per_locus[is.finite(fis_per_locus)]
observed_global_fis <- mean(observed_fis_values, na.rm = TRUE)

cat("Target Fis for simulation:", target_fis, "\n")
cat("Observed Global Fis from data:", round(observed_global_fis, 4), "\n")

if (observed_global_fis > 0.01) {
  message("SUCCESS: Cryptic population structure successfully simulated (Fis > 0).")
} else {
  warning("WARNING: Simulated data shows very low Fis. Beta-binomial tests may not be informative.")
}

# --- 5. Run Full Analysis on All Methods ---
methods_to_test <- c("wasserstein", "bhd", "betabin", "betabin_dynamic")
results_summary <- list()

for (method_name in methods_to_test) {
  message(paste0("\n--- Testing Method: ", method_name, " ---"))
  
  args <- list(
    genind_object = genind_test,
    method = method_name,
    jackknife_proportion = jackknife_prop_test,
    n_cores = 1 # Set to >1 for parallel execution
  )
  
  start_time <- Sys.time()
  res_test <- try(do.call(dnadot_snp, args), silent = FALSE)
  end_time <- Sys.time()
  
  if (!inherits(res_test, "try-error")) {
    message(paste0("Method ran successfully! (Time: ", round(difftime(end_time, start_time, units = "secs"), 2), "s)"))
    
    N_estimate <- res_test$N_est
    message("Selected Loci Details (Top 5):")
    print(head(res_test$selected_loci_details, 5))
    
    cat(paste0(method_name, " Estimated Nc: ", round(N_estimate, 2), "\n"))
    results_summary[[method_name]] <- N_estimate
    
  } else {
    message(paste0("\nMethod ", method_name, " failed with an error."))
    results_summary[[method_name]] <- NA
  }
}

message("\n--- Full test script finished ---")
message("Summary of Estimates (True Nc = ", true_Nc_total, "; Target Fis = ", target_fis, "):")
print(data.frame(Method = names(results_summary), N_est = unlist(results_summary)))
