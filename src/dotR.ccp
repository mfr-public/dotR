// File: src/dotR.cpp
// Complete C++ backend for dotR package.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string>
#include <utility> // Required for std::pair

#ifndef R_POSINF
#define R_POSINF R_PosInf
#endif

using namespace Rcpp;

// ============================================================================
// --- Global Constants for Optimization ---
// ============================================================================

// log(1e-16) approx -36.84. Used for Adaptive Windowing.
const double LOG_TOLERANCE = std::log(1e-16); 

// ============================================================================
// --- Probability Distribution Functions (Merged) ---
// ============================================================================

//' Calculate Hypergeometric Probability (Stable Version)
//' @keywords internal
// [[Rcpp::export]]
double dotR_hypergeo_cpp(int k_drawn, int n_draws, int K_successes_in_pop, int N_total_items_in_pop) {
    if (k_drawn < 0 || n_draws < 0 || K_successes_in_pop < 0 || N_total_items_in_pop < 0) return 0.0;
    if (k_drawn > K_successes_in_pop || k_drawn > n_draws) return 0.0;
    if (n_draws - k_drawn > N_total_items_in_pop - K_successes_in_pop) return 0.0;
    if (n_draws > N_total_items_in_pop) return 0.0;
    if (n_draws == 0) return (k_drawn == 0) ? 1.0 : 0.0;

    double log_prob = R::lchoose(static_cast<double>(K_successes_in_pop), static_cast<double>(k_drawn)) +
                      R::lchoose(static_cast<double>(N_total_items_in_pop - K_successes_in_pop), static_cast<double>(n_draws - k_drawn)) -
                      R::lchoose(static_cast<double>(N_total_items_in_pop), static_cast<double>(n_draws));

    return std::isfinite(log_prob) ? std::exp(log_prob) : 0.0;
}

//' Calculate Beta-binomial Probability
//' @keywords internal
// [[Rcpp::export]]
double dotR_betabin_cpp(int k, int n, double p, double theta) {
    if (k < 0 || n < 0 || p <= 0 || p >= 1 || theta < 0 || theta >= 1) {
        return 0.0;
    }
    if (k > n) return 0.0;
    
    // Optimization for low theta (Binomial approximation)
    if (theta < 1e-9) {
        return R::dbinom(k, n, p, false);
    }

    // Re-parameterization
    double alpha = p * ((1.0 - theta) / theta);
    double beta = (1.0 - p) * ((1.0 - theta) / theta);

    if (alpha <= 0 || beta <= 0) return 0.0;

    // Log-space calculation
    double log_prob = R::lchoose(n, k) + R::lbeta(k + alpha, n - k + beta) - R::lbeta(alpha, beta);
    
    return std::isfinite(log_prob) ? std::exp(log_prob) : 0.0;
}

// ============================================================================
// --- Beta-Hypergeometric PMF (Optimized Implementation) ---
// ============================================================================

// --- BHG Optimization Helpers ---

struct BHGConstants {
    int k, n, N;
    double alpha, beta;
    double log_beta_alpha_beta; 
    double log_binom_N_n;       
};

// Helper function to calculate L_K (Log-Joint Probability)
double calculate_LK(int K, const BHGConstants& consts) {
    // Check bounds
    if (K < consts.k || K > consts.N - (consts.n - consts.k)) {
        return -R_POSINF;
    }

    // Level 2: Log-Hypergeometric
    double log_hypergeo = R::lchoose(K, consts.k) + 
                          R::lchoose(consts.N - K, consts.n - consts.k) - 
                          consts.log_binom_N_n;

    // Level 1: Log-Beta-Binomial
    double log_betabin = R::lchoose(consts.N, K) + 
                         R::lbeta(K + consts.alpha, consts.N - K + consts.beta) - 
                         consts.log_beta_alpha_beta;

    return log_hypergeo + log_betabin;
}

// Function to find the mode (K_mode) and the maximum log-probability (M) using Hill Climbing
std::pair<int, double> find_mode(const BHGConstants& consts) {
    // Initialize near expected value E[K] ≈ Np
    double p = consts.alpha / (consts.alpha + consts.beta);
    int K_start = static_cast<int>(std::round(consts.N * p));
    
    int K_min = consts.k;
    int K_max = consts.N - (consts.n - consts.k);
    K_start = std::max(K_min, std::min(K_max, K_start));

    int K_mode = K_start;
    double M = calculate_LK(K_mode, consts);

    // Hill climbing: Search upwards
    int K = K_mode + 1;
    while (K <= K_max) {
        double LK_next = calculate_LK(K, consts);
        if (LK_next > M) {
            M = LK_next;
            K_mode = K;
            K++;
        } else {
            break; 
        }
    }

    // Hill climbing: Search downwards
    K = K_start - 1;
    while (K >= K_min) {
        double LK_prev = calculate_LK(K, consts);
        if (LK_prev > M) {
            M = LK_prev;
            K_mode = K;
            K--;
        } else {
            break;
        }
    }
    
    return {K_mode, M};
}

// --- END BHG Optimization Helpers ---


//' Calculate Beta-Hypergeometric Probability (Stabilized & Optimized Implementation)
//' @keywords internal
// [[Rcpp::export]]
double dotR_bhg_pmf_cpp(int k, int n, int N, double p, double theta) {

    // --- 1. Input Validation ---
    if (k < 0 || n <= 0 || N < n || p <= 0 || p >= 1 || theta < 0 || theta >= 1) {
        return 0.0;
    }
    if (k > n) return 0.0;

    // --- 2. Optimization for Low Theta (Panmixia) ---
    if (theta < 1e-9) {
        int K_hypothesized = static_cast<int>(std::round(p * N));
        return dotR_hypergeo_cpp(k, n, K_hypothesized, N); 
    }

    // --- 3. Calculate alpha and beta (Re-parameterization) ---
    double factor = (1.0 - theta) / theta;
    double alpha = p * factor;
    double beta = (1.0 - p) * factor;

    if (alpha <= 0 || beta <= 0) return 0.0;

    // --- 4. Initialize Constants Structure ---
    BHGConstants consts;
    consts.k = k;
    consts.n = n;
    consts.N = N;
    consts.alpha = alpha;
    consts.beta = beta;
    consts.log_beta_alpha_beta = R::lbeta(alpha, beta);
    consts.log_binom_N_n = R::lchoose(N, n);

    int K_min = k;
    int K_max = N - (n - k);

    // --- 5. Optimization Step 1: Find the Mode (Hill Climbing) ---
    std::pair<int, double> mode_info = find_mode(consts);
    int K_mode = mode_info.first;
    double M = mode_info.second; // Max log-probability

    if (!std::isfinite(M)) {
        return 0.0;
    }

    // --- 6. Optimization Step 2: Adaptive Window Summation (Integrated LSE) ---
    
    double sum_exp = 1.0; // Initialize LSE sum with the mode term

    // Search Downwards
    int K = K_mode - 1;
    while (K >= K_min) {
        double L_K = calculate_LK(K, consts);
        if (L_K - M < LOG_TOLERANCE) {
            break;
        }
        sum_exp += std::exp(L_K - M);
        K--;
    }

    // Search Upwards
    K = K_mode + 1;
    while (K <= K_max) {
        double L_K = calculate_LK(K, consts);
        if (L_K - M < LOG_TOLERANCE) {
            break;
        }
        sum_exp += std::exp(L_K - M);
        K++;
    }

    // --- 7. Final Calculation (LSE result) ---
    double log_pmf = M + std::log(sum_exp);

    return std::isfinite(log_pmf) ? std::exp(log_pmf) : 0.0;
}


// ============================================================================
// --- MDE Core Logic Functions (Unified) ---
// ============================================================================


//' NEW IMPLEMENTATION: 3D Beta-Hypergeometric MDE with Informed L2 Regularization
//' @keywords internal
// [[Rcpp::export]]
List dotR_dnadot_bhg_3d_cpp(
    IntegerMatrix genotypes_locus,
    int n_individuals,
    int n_subsample_individuals,
    NumericVector N_try_values,
    NumericVector p_try_values,
    NumericVector theta_try_values,
    int target_allele_idx,
    double theta_global_hat, // The empirical anchor point (theta_hat)
    double theta_sigma_sq    // The uncertainty (sigma^2)
) {
    // --- 1. Setup and Initialization ---
    int ploidy = (n_individuals > 0 && genotypes_locus.nrow() > 0) ? genotypes_locus.ncol() : 2;
    int n_gene_copies_per_subsample = n_subsample_individuals * ploidy;

    // Basic input validation
    if (n_individuals == 0 || n_subsample_individuals == 0) {
        return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL, _["theta_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    int num_jackknife_subsamples = n_individuals - n_subsample_individuals + 1;
    if (num_jackknife_subsamples <= 0) {
         return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL, _["theta_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    // --- 2. Calculate Observed Distribution (O) ---
    std::vector<double> observed_counts_O_bins(n_gene_copies_per_subsample + 1, 0.0);
    NumericVector subsample_target_allele_counts(num_jackknife_subsamples);

    for (int ss_idx = 0; ss_idx < num_jackknife_subsamples; ++ss_idx) {
        int count = 0;
        // Calculate counts in the subsample, handling potential missing data (NA)
        for (int i = 0; i < n_subsample_individuals; ++i) {
            for (int j = 0; j < ploidy; ++j) {
                int allele = genotypes_locus(ss_idx + i, j);
                if (!IntegerVector::is_na(allele)) {
                    if (allele == target_allele_idx) {
                        count++;
                    }
                }
            }
        }
        
        subsample_target_allele_counts[ss_idx] = count;
        if (count >= 0 && count <= n_gene_copies_per_subsample) {
            observed_counts_O_bins[count]++;
        }
    }

    // --- 3. 3D Grid Search ---
    double min_penalized_discrepancy = R_POSINF;
    double N_best = NA_REAL, p_best = NA_REAL, theta_best = NA_REAL;

    // Pre-calculate the regularization denominator (2 * sigma^2)
    // If theta_sigma_sq is Inf (no regularization), the denominator is Inf, making the penalty 0.
    double regularization_denom = (std::isfinite(theta_sigma_sq) && theta_sigma_sq > 0) ? (2.0 * theta_sigma_sq) : R_POSINF;

    for (int n_idx = 0; n_idx < N_try_values.size(); ++n_idx) {
        int N_try_ind = static_cast<int>(N_try_values[n_idx]);
        if (N_try_ind < n_subsample_individuals) continue;
        int N_try_alleles = N_try_ind * ploidy;

        for (int p_idx = 0; p_idx < p_try_values.size(); ++p_idx) {
            double p_try = p_try_values[p_idx];

            for (int t_idx = 0; t_idx < theta_try_values.size(); ++t_idx) {
                double theta_try = theta_try_values[t_idx];

                // --- 3a. Calculate Discrepancy (Wasserstein Distance) ---
                double current_discrepancy = 0.0;
                double cdf_obs = 0.0;
                double cdf_exp = 0.0;

                for (int k = 0; k <= n_gene_copies_per_subsample; ++k) {
                    double p_obs = observed_counts_O_bins[k] / num_jackknife_subsamples;
                    
                    // Calculate Expected probability (E) using the optimized BHG PMF
                    double p_exp = dotR_bhg_pmf_cpp(k, n_gene_copies_per_subsample, N_try_alleles, p_try, theta_try);

                    // Update CDFs and calculate Wasserstein distance
                    cdf_obs += p_obs;
                    cdf_exp += p_exp;
                    current_discrepancy += std::abs(cdf_obs - cdf_exp);
                }

                // --- 3b. Calculate Penalty (Informed L2) ---
                double penalty = 0.0;
                if (std::isfinite(regularization_denom)) {
                    // Penalty = (theta_try - theta_hat)^2 / (2 * sigma^2)
                    double diff = theta_try - theta_global_hat;
                    penalty = (diff * diff) / regularization_denom;
                }

                // --- 3c. Combine and Update ---
                double penalized_discrepancy = current_discrepancy + penalty;

                if (std::isfinite(penalized_discrepancy) && penalized_discrepancy < min_penalized_discrepancy) {
                    min_penalized_discrepancy = penalized_discrepancy;
                    N_best = N_try_values[n_idx];
                    p_best = p_try;
                    theta_best = theta_try;
                }
            }
        }
    }

    // --- 4. Return Results ---
    return List::create(
        _["N_est"] = N_best, 
        _["p_est"] = p_best,
        _["theta_est"] = theta_best,
        _["min_discrepancy"] = min_penalized_discrepancy, 
        _["subsample_target_allele_counts"] = subsample_target_allele_counts
    );
}


// (The following functions still need to be ported from dnadot-beta.cpp or removed from 0.4 as not needed - decision to be based on sim results)

//' Minimum Distance Estimation (MDE) Core Logic (C++ Core Logic)
//' @keywords internal
// [[Rcpp::export]]
List dotR_dnadot_mde_cpp(
    IntegerMatrix genotypes_locus,
    int n_individuals,
    // ... (Parameters identical to dnadot-beta.cpp)
) {
    // (Implementation identical to dnadot-beta.cpp)
    // ...
}


//' Buffered Histogram Discrepancy (BHD) Method (C++ Core Logic)
//' @keywords internal
// [[Rcpp::export]]
List dotR_dnadot_bhd_cpp(
    // ... (Parameters identical to dnadot-beta.cpp)
) {
    // (Implementation identical to dnadot-beta.cpp)
    // ...
}

//' MDE using Beta-binomial distribution (Fixed Fis)
//' @keywords internal
// [[Rcpp::export]]
List dotR_dnadot_betabin_cpp(
    // ... (Parameters identical to dnadot-beta.cpp)
) {
    // (Implementation identical to dnadot-beta.cpp)
    // ...
}

//' MDE using Beta-binomial distribution (Dynamic Fis)
//' @keywords internal
// [[Rcpp::export]]
List dotR_dnadot_betabin_dynamic_cpp(
    // ... (Parameters identical to dnadot-beta.cpp)
) {
    // (Implementation identical to dnadot-beta.cpp)
    // ...
}
