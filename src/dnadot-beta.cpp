// File: src/dnadot.cpp

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string>

#ifndef R_POSINF
#define R_POSINF R_PosInf
#endif

using namespace Rcpp;

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
    
    // If theta is very small, the distribution is approximately binomial.
    if (theta < 1e-9) {
        return R::dbinom(k, n, p, false);
    }

    double alpha = p * ((1.0 - theta) / theta);
    double beta = (1.0 - p) * ((1.0 - theta) / theta);

    if (alpha <= 0 || beta <= 0) return 0.0;

    double log_prob = R::lchoose(n, k) + R::lbeta(k + alpha, n - k + beta) - R::lbeta(alpha, beta);
    
    return std::isfinite(log_prob) ? std::exp(log_prob) : 0.0;
}


//' Minimum Distance Estimation (MDE) Core Logic (C++ Core Logic)
//' @keywords internal
// [[Rcpp::export]]
List dotR_dnadot_mde_cpp(
    IntegerMatrix genotypes_locus,
    int n_individuals,
    int n_subsample_individuals,
    NumericVector N_try_values,
    NumericVector p_try_values,
    int target_allele_idx,
    String discrepancy_type, 
    double reg_strength,
    double reg_N_threshold,
    bool debug = false
) {
    int ploidy = (n_individuals > 0 && genotypes_locus.nrow() > 0) ? genotypes_locus.ncol() : 2;
    int n_gene_copies_per_subsample = n_subsample_individuals * ploidy;

    if (n_individuals == 0 || n_subsample_individuals == 0) {
        return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    int num_jackknife_subsamples = n_individuals - n_subsample_individuals + 1;
    if (num_jackknife_subsamples <= 0) {
        return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    std::vector<double> observed_counts_O_bins(n_gene_copies_per_subsample + 1, 0.0);
    NumericVector subsample_target_allele_counts(num_jackknife_subsamples);

    for (int ss_idx = 0; ss_idx < num_jackknife_subsamples; ++ss_idx) {
        int count = 0;
        for (int i = 0; i < n_subsample_individuals; ++i) {
            for (int j = 0; j < ploidy; ++j) {
                if (!IntegerVector::is_na(genotypes_locus(ss_idx + i, j)) && genotypes_locus(ss_idx + i, j) == target_allele_idx) {
                    count++;
                }
            }
        }
        subsample_target_allele_counts[ss_idx] = count;
        if (count >= 0 && count <= n_gene_copies_per_subsample) {
            observed_counts_O_bins[count]++;
        }
    }

    double min_total_discrepancy = R_POSINF;
    double N_best = NA_REAL, p_best = NA_REAL;

    for (int n_idx = 0; n_idx < N_try_values.size(); ++n_idx) {
        int N_try_ind = static_cast<int>(N_try_values[n_idx]);
        if (N_try_ind < n_subsample_individuals) continue;
        int N_try_alleles = N_try_ind * ploidy;

        for (int p_idx = 0; p_idx < p_try_values.size(); ++p_idx) {
            double p_try = p_try_values[p_idx];
            int K_alleles = static_cast<int>(std::round(p_try * N_try_alleles));
            double current_discrepancy = 0.0;

            if (discrepancy_type == "neg_log_likelihood") {
                double log_lik = 0.0;
                bool possible = true;
                for (int k = 0; k <= n_gene_copies_per_subsample; ++k) {
                    if (observed_counts_O_bins[k] > 0) {
                        double p_exp = dotR_hypergeo_cpp(k, n_gene_copies_per_subsample, K_alleles, N_try_alleles);
                        if (p_exp > 1e-12) {
                            log_lik += observed_counts_O_bins[k] * std::log(p_exp);
                        } else {
                            possible = false;
                            break;
                        }
                    }
                }
                current_discrepancy = possible ? -log_lik : R_POSINF;

            } else if (discrepancy_type == "absolute_diff" || discrepancy_type == "wasserstein") {
                double cdf_obs = 0.0;
                double cdf_exp = 0.0;
                for (int k = 0; k <= n_gene_copies_per_subsample; ++k) {
                    double p_obs = observed_counts_O_bins[k] / num_jackknife_subsamples;
                    double p_exp = dotR_hypergeo_cpp(k, n_gene_copies_per_subsample, K_alleles, N_try_alleles);
                    if (discrepancy_type == "absolute_diff") {
                        current_discrepancy += std::abs(p_obs - p_exp);
                    } else { 
                        cdf_obs += p_obs;
                        cdf_exp += p_exp;
                        current_discrepancy += std::abs(cdf_obs - cdf_exp);
                    }
                }
                if (reg_strength > 0.0 && N_try_ind > reg_N_threshold && std::isfinite(reg_N_threshold)) {
                    if (N_try_ind > 0 && reg_N_threshold > 0) {
                        current_discrepancy += reg_strength * std::log(static_cast<double>(N_try_ind) / reg_N_threshold);
                    }
                }
            } else {
                stop("Invalid discrepancy_type provided.");
            }

            if (std::isfinite(current_discrepancy) && current_discrepancy < min_total_discrepancy) {
                min_total_discrepancy = current_discrepancy;
                N_best = N_try_values[n_idx];
                p_best = p_try;
            }
        }
    }

    return List::create(_["N_est"] = N_best, _["p_est"] = p_best,
                        _["min_discrepancy"] = min_total_discrepancy, 
                        _["subsample_target_allele_counts"] = subsample_target_allele_counts);
}


//' Buffered Histogram Discrepancy (BHD) Method (C++ Core Logic)
//' @keywords internal
// [[Rcpp::export]]
List dotR_dnadot_bhd_cpp(
    IntegerMatrix genotypes_locus,
    int n_individuals,
    int n_subsample_individuals,
    NumericVector N_try_values,
    NumericVector p_try_values,
    int target_allele_idx,
    int bhd_buffer_width
) {
    int ploidy = (n_individuals > 0 && genotypes_locus.nrow() > 0) ? genotypes_locus.ncol() : 2;
    int n_gene_copies_per_subsample = n_subsample_individuals * ploidy;

    if (n_individuals == 0 || n_subsample_individuals == 0) {
        return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    int num_jackknife_subsamples = n_individuals - n_subsample_individuals + 1;
    if (num_jackknife_subsamples <= 0) {
        return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    NumericVector subsample_target_allele_counts(num_jackknife_subsamples);
    for (int ss_idx = 0; ss_idx < num_jackknife_subsamples; ++ss_idx) {
        int count = 0;
        for (int i = 0; i < n_subsample_individuals; ++i) {
            for (int j = 0; j < ploidy; ++j) {
                if (!IntegerVector::is_na(genotypes_locus(ss_idx + i, j)) && genotypes_locus(ss_idx + i, j) == target_allele_idx) {
                    count++;
                }
            }
        }
        subsample_target_allele_counts[ss_idx] = count;
    }

    double min_total_discrepancy = R_POSINF;
    double N_best = NA_REAL, p_best = NA_REAL;

    for (int n_idx = 0; n_idx < N_try_values.size(); ++n_idx) {
        int N_try_ind = static_cast<int>(N_try_values[n_idx]);
        if (N_try_ind < n_subsample_individuals) continue;

        for (int p_idx = 0; p_idx < p_try_values.size(); ++p_idx) {
            double p_try = p_try_values[p_idx];

            std::vector<double> E_bins(n_gene_copies_per_subsample + 1, 0.0);
            for (int k = 0; k <= n_gene_copies_per_subsample; ++k) {
                E_bins[k] = dotR_hypergeo_cpp(k, n_gene_copies_per_subsample, 
                                              static_cast<int>(std::round(p_try * N_try_ind * ploidy)), 
                                              N_try_ind * ploidy);
            }

            double total_discrepancy = 0.0;
            for (int ss_idx = 0; ss_idx < num_jackknife_subsamples; ++ss_idx) {
                int k_obs = static_cast<int>(subsample_target_allele_counts[ss_idx]);
                
                std::vector<double> O_buffered(n_gene_copies_per_subsample + 1, 0.0);
                double normalizer = 0.0;
                for (int w = -bhd_buffer_width; w <= bhd_buffer_width; ++w) {
                    int bin = k_obs + w;
                    if (bin >= 0 && bin <= n_gene_copies_per_subsample) {
                        double weight = bhd_buffer_width - std::abs(w) + 1;
                        O_buffered[bin] = weight;
                        normalizer += weight;
                    }
                }
                if (normalizer > 0) {
                    for (size_t i = 0; i < O_buffered.size(); ++i) O_buffered[i] /= normalizer;
                }

                double subsample_d = 0.0;
                for (int k = 0; k <= n_gene_copies_per_subsample; ++k) {
                    subsample_d += std::abs(O_buffered[k] - E_bins[k]);
                }
                total_discrepancy += subsample_d;
            }

            if (std::isfinite(total_discrepancy) && total_discrepancy < min_total_discrepancy) {
                min_total_discrepancy = total_discrepancy;
                N_best = N_try_values[n_idx];
                p_best = p_try;
            }
        }
    }

    return List::create(_["N_est"] = N_best, _["p_est"] = p_best,
                        _["min_discrepancy"] = min_total_discrepancy, 
                        _["subsample_target_allele_counts"] = subsample_target_allele_counts);
}

//' MDE using Beta-binomial distribution (Fixed Fis)
//' @keywords internal
// [[Rcpp::export]]
List dotR_dnadot_betabin_cpp(
    IntegerMatrix genotypes_locus,
    int n_individuals,
    int n_subsample_individuals,
    NumericVector N_try_values,
    NumericVector p_try_values,
    int target_allele_idx,
    double global_fis
) {
    int ploidy = (n_individuals > 0 && genotypes_locus.nrow() > 0) ? genotypes_locus.ncol() : 2;
    int n_gene_copies_per_subsample = n_subsample_individuals * ploidy;

    if (n_individuals == 0 || n_subsample_individuals == 0) {
        return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    int num_jackknife_subsamples = n_individuals - n_subsample_individuals + 1;
    if (num_jackknife_subsamples <= 0) {
        return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    std::vector<double> observed_counts_O_bins(n_gene_copies_per_subsample + 1, 0.0);
    NumericVector subsample_target_allele_counts(num_jackknife_subsamples);

    for (int ss_idx = 0; ss_idx < num_jackknife_subsamples; ++ss_idx) {
        int count = 0;
        for (int i = 0; i < n_subsample_individuals; ++i) {
            for (int j = 0; j < ploidy; ++j) {
                if (!IntegerVector::is_na(genotypes_locus(ss_idx + i, j)) && genotypes_locus(ss_idx + i, j) == target_allele_idx) {
                    count++;
                }
            }
        }
        subsample_target_allele_counts[ss_idx] = count;
        if (count >= 0 && count <= n_gene_copies_per_subsample) {
            observed_counts_O_bins[count]++;
        }
    }

    double min_total_discrepancy = R_POSINF;
    double N_best = NA_REAL, p_best = NA_REAL;
    
    for (int n_idx = 0; n_idx < N_try_values.size(); ++n_idx) {
        int N_try_ind = static_cast<int>(N_try_values[n_idx]);
        if (N_try_ind < n_subsample_individuals) continue;

        for (int p_idx = 0; p_idx < p_try_values.size(); ++p_idx) {
            double p_try = p_try_values[p_idx];
            double current_discrepancy = 0.0;
            
            double cdf_obs = 0.0;
            double cdf_exp = 0.0;
            for (int k = 0; k <= n_gene_copies_per_subsample; ++k) {
                double p_obs = observed_counts_O_bins[k] / num_jackknife_subsamples;
                double p_exp = dotR_betabin_cpp(k, n_gene_copies_per_subsample, p_try, global_fis);
                
                cdf_obs += p_obs;
                cdf_exp += p_exp;
                current_discrepancy += std::abs(cdf_obs - cdf_exp);
            }

            if (std::isfinite(current_discrepancy) && current_discrepancy < min_total_discrepancy) {
                min_total_discrepancy = current_discrepancy;
                N_best = N_try_values[n_idx];
                p_best = p_try;
            }
        }
    }

    return List::create(_["N_est"] = N_best, _["p_est"] = p_best,
                        _["min_discrepancy"] = min_total_discrepancy, 
                        _["subsample_target_allele_counts"] = subsample_target_allele_counts);
}

//' --- NEW FUNCTION ---
//' MDE using Beta-binomial distribution (Dynamic Fis)
//' @keywords internal
// [[Rcpp::export]]
List dotR_dnadot_betabin_dynamic_cpp(
    IntegerMatrix genotypes_locus,
    int n_individuals,
    int n_subsample_individuals,
    NumericVector N_try_values,
    NumericVector p_try_values,
    int target_allele_idx,
    double migration_rate_m
) {
    int ploidy = (n_individuals > 0 && genotypes_locus.nrow() > 0) ? genotypes_locus.ncol() : 2;
    int n_gene_copies_per_subsample = n_subsample_individuals * ploidy;

    if (n_individuals == 0 || n_subsample_individuals == 0) {
        return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    int num_jackknife_subsamples = n_individuals - n_subsample_individuals + 1;
    if (num_jackknife_subsamples <= 0) {
        return List::create(_["N_est"] = NA_REAL, _["p_est"] = NA_REAL,
                            _["min_discrepancy"] = NA_REAL, _["subsample_target_allele_counts"] = NumericVector(0));
    }

    std::vector<double> observed_counts_O_bins(n_gene_copies_per_subsample + 1, 0.0);
    NumericVector subsample_target_allele_counts(num_jackknife_subsamples);

    for (int ss_idx = 0; ss_idx < num_jackknife_subsamples; ++ss_idx) {
        int count = 0;
        for (int i = 0; i < n_subsample_individuals; ++i) {
            for (int j = 0; j < ploidy; ++j) {
                if (!IntegerVector::is_na(genotypes_locus(ss_idx + i, j)) && genotypes_locus(ss_idx + i, j) == target_allele_idx) {
                    count++;
                }
            }
        }
        subsample_target_allele_counts[ss_idx] = count;
        if (count >= 0 && count <= n_gene_copies_per_subsample) {
            observed_counts_O_bins[count]++;
        }
    }

    double min_total_discrepancy = R_POSINF;
    double N_best = NA_REAL, p_best = NA_REAL;
    
    for (int n_idx = 0; n_idx < N_try_values.size(); ++n_idx) {
        int N_try_ind = static_cast<int>(N_try_values[n_idx]);
        if (N_try_ind < n_subsample_individuals) continue;

        // --- DYNAMIC PART: Calculate hypothesized Fis for this N_try ---
        double hypothesized_fis = 1.0 / (4.0 * N_try_ind * migration_rate_m + 1.0);
        // Ensure Fis is within a valid range [0, 1)
        hypothesized_fis = std::max(0.0, std::min(hypothesized_fis, 1.0 - 1e-9));

        for (int p_idx = 0; p_idx < p_try_values.size(); ++p_idx) {
            double p_try = p_try_values[p_idx];
            double current_discrepancy = 0.0;
            
            double cdf_obs = 0.0;
            double cdf_exp = 0.0;
            for (int k = 0; k <= n_gene_copies_per_subsample; ++k) {
                double p_obs = observed_counts_O_bins[k] / num_jackknife_subsamples;
                // Use Beta-binomial with the DYNAMICALLY calculated Fis
                double p_exp = dotR_betabin_cpp(k, n_gene_copies_per_subsample, p_try, hypothesized_fis);
                
                cdf_obs += p_obs;
                cdf_exp += p_exp;
                current_discrepancy += std::abs(cdf_obs - cdf_exp);
            }

            if (std::isfinite(current_discrepancy) && current_discrepancy < min_total_discrepancy) {
                min_total_discrepancy = current_discrepancy;
                N_best = N_try_values[n_idx];
                p_best = p_try;
            }
        }
    }

    return List::create(_["N_est"] = N_best, _["p_est"] = p_best,
                        _["min_discrepancy"] = min_total_discrepancy, 
                        _["subsample_target_allele_counts"] = subsample_target_allele_counts);
}
