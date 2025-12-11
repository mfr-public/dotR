# A Unified Framework for Estimating Census Size in Structured Populations

**3D Regularised MDE using the Beta-Hypergeometric Distribution**

**Author:** Mark Richardson
**Date:** 2025-10-31

## 1. Introduction and Motivation

The original `dnadot` method provides a Minimum Distance Estimation (MDE) framework for estimating census population size ($N_c$). It relies on the Hypergeometric distribution, which assumes a single, panmictic population. When population structure exists (e.g., Wahlund effect, cryptic relatedness), the observed allele counts exhibit overdispersion relative to the Hypergeometric expectation, leading to bias.

My previous attempts to address this used the Beta-Binomial model. However, the Beta-Binomial model inherently assumes sampling *with* replacement (or an infinite population size). When sampling *without* replacement from a finite population—as is necessary for estimating $N_c$—this assumption is violated. Furthermore, prior implementations using calibrated Beta-Binomial models introduced instability due to reliance on idealised population genetic models (e.g., Wright's Island Model, assuming $N_c \approx N_e$) and fragile two-step calibration processes.

This document outlines a unified framework that simultaneously accounts for the variance due to finite sampling without replacement and the variance due to population structure. The core of this approach is the **Beta-Hypergeometric (BHG)** distribution.

## 2. The Foundational Model: The Beta-Hypergeometric (BHG) Distribution

The BHG distribution is the mathematically correct model for sampling without replacement from a finite population where the underlying allele frequencies are themselves variable due to structure.

### 2.1. Conceptual Framework

The BHG distribution arises as a hierarchical mixture model, capturing two levels of stochasticity:

1. **Structural Stochasticity:** We model the variation in allele frequencies across the population using a Beta distribution. This reflects the outcome of genetic drift, migration, and mating patterns.

2. **Sampling Stochasticity:** When taking a finite sample from a finite population, we must account for the fact that removing individuals changes the remaining allele frequencies (sampling without replacement), modeled by the Hypergeometric distribution.

### 2.2. Parameterisation and the Beta Distribution

We aim to jointly estimate three core parameters. For clarity, we define $N_c$ as the census population size (number of individuals) and $N$ as the total number of gene copies (e.g., $N=2N_c$ for diploids).

The three parameters are:

* $N$: Total gene copies in the population.

* $p$: The mean allele frequency across the entire population.

* $\theta$: The overdispersion parameter (analogous to $F_{IS}$), quantifying the deviation from panmixia.

The Beta distribution, $Beta(\alpha, \beta)$, is uniquely suited to modeling probabilities bounded between 0 and 1. Following standard population genetics practice (e.g., the Balding-Nichols model), we re-parameterise the shape parameters $\alpha$ and $\beta$ using $p$ and $\theta$.

The parameters are linked by the following equations:

$$
p = \frac{\alpha}{\alpha+\beta}
$$

$$
\theta = \frac{1}{\alpha+\beta+1}
$$

Solving for $\alpha$ and $\beta$ yields:

$$
\alpha = p\left(\frac{1-\theta}{\theta}\right)
$$

$$
\beta = (1-p)\left(\frac{1-\theta}{\theta}\right)
$$

These parameters define the shape of the allele frequency distribution within the population being modeled.

### 2.3. The BHG Probability Mass Function (PMF)

The BHG PMF describes the probability of observing exactly $k$ copies of a target allele in a sample of size $n$ gene copies, given $N, \alpha,$ and $\beta$.

#### 2.3.1. Level 1: Population Realisation (Beta-Binomial)

The total count of the target allele in the entire population, $K$, is a latent variable. It follows a Beta-Binomial distribution, reflecting the realisation of the population structure:

$$
P(K|N, \alpha, \beta) = \binom{N}{K} \frac{B(K+\alpha, N-K+\beta)}{B(\alpha, \beta)}
$$

*(Where B is the Beta function).*

#### 2.3.2. Level 2: Sampling (Hypergeometric)

The observed count in the sample, $k$, conditional on the realised population count $K$, follows a Hypergeometric distribution:

$$
P(k|n, N, K) = \frac{\binom{K}{k}\binom{N-K}{n-k}}{\binom{N}{n}}
$$

#### 2.3.3. The Marginal BHG PMF

The full BHG PMF is obtained by marginalising (summing) over all possible values of the latent variable $K$. The summation limits for $K$ are constrained such that $K$ must be at least $k$, and at most $N-(n-k)$.

$$
P(k|n, N, \alpha, \beta) = \sum_{K=k}^{N-n+k} P(k|n, N, K) \cdot P(K|N, \alpha, \beta)
$$

This equation links the observed data $(k,n)$ to the parameters we wish to estimate $(N, p, \theta)$.

### 2.4. Computational Implementation and Optimization

The BHG PMF is computationally intensive because the marginalisation requires summing over the latent variable $K$, where each term involves evaluating computationally expensive Beta functions and combinatorial terms. These terms can span many orders of magnitude, leading to numerical instability if computed directly. To ensure feasibility and stability, the implementation employs several critical optimisations:

1. **Log-Space Calculations:** All intermediate probabilities are calculated in log-space to prevent arithmetic overflow and underflow.

2. **Log-Sum-Exp (LSE) Trick:** The summation over $K$ is performed using the LSE trick. We identify the maximum log-probability ($M$) and calculate the final log-probability as $M + \log(\sum \exp(L_K - M))$. This stabilizes the summation of terms of vastly different magnitudes.

3. **Adaptive Windowing:** Rather than summing over the entire range of $K$, the implementation uses a Hill Climbing algorithm to efficiently locate the mode of the distribution. It then sums only the numerically significant terms within a window around this mode (e.g., terms within $10^{-16}$ of the mode's probability). This drastically reduces computation time without loss of precision.

4. **Low Theta Optimization:** When $\theta$ approaches zero (near panmixia), the BHG distribution converges to the Hypergeometric distribution. The implementation detects this regime and switches to the faster Hypergeometric calculation.

## 3. The 3D Minimum Distance Estimation (MDE) Framework

We embed the BHG PMF into the MDE framework by expanding the search space to three dimensions.

### 3.1. The MDE Procedure

1. **Observe (O):** Generate the empirical distribution of allele counts (O) by analysing many jackknife subsamples of the input data.

2. **Hypothesize:** Create a 3D grid of parameter triplets: $(N_{try}, p_{try}, \theta_{try})$.

3. **Expect (E):** For each triplet, calculate the corresponding $\alpha$ and $\beta$, then calculate the Expected distribution (E) using the optimised BHG PMF.

4. **Compare:** Measure the distance between O and E using a robust discrepancy metric.

5. **Estimate:** The triplet that minimises the penalised distance (see Section 4) is chosen as the best estimate.

### 3.2. Discrepancy Metrics

We consider robust options such as the Wasserstein Distance ($D_W$) or the Buffered Histogram Discrepancy ($D_{BHD}$).

#### 3.2.1. Wasserstein Distance ($D_W$)

The Wasserstein distance (Earth Mover's Distance) compares the Cumulative Distribution Functions (CDFs).

$$
D_W = \sum_{k=0}^{n} |CDF_O(k)-CDF_E(k)|
$$

#### 3.2.2. Buffered Histogram Discrepancy ($D_{BHD}$)

The BHD attempts to account for the inherent noise in the observed data by "smoothing" the distribution O.

1. **Buffering:** For each jackknife observation $k_{obs}$, a triangular kernel is created. With a buffer width $w$:

   $$
   weight(k, k_{obs}, w) = \max(0, w - |k-k_{obs}| + 1)
   $$

2. **Smoothed O:** The final buffered observed distribution, $O_{B}$, is the normalised average of these kernels.

3. **Discrepancy:**

   $$
   D_{BHD} = \sum_{k=0}^{n} |O_B(k)-E(k)|
   $$

## 4. Regularisation: Addressing Parameter Identifiability

A fundamental challenge in this 3D framework is **weak parameter identifiability**. The variance observed in the jackknife samples is a convolution of Sampling Variance (increases as $N$ decreases) and Structural Variance (increases as $\theta$ increases). The model may struggle to distinguish between a large, highly structured population (High $N$, High $\theta$) and a smaller, less structured population (Low $N$, Low $\theta$), as both can produce similar observed variance.

We address this using **Penalized MDE** (Regularisation), leveraging genome-wide information to stabilise the estimation of $\theta$, thereby improving how we can identify $N$.

### 4.1. Empirical Anchoring

The key insight is that the overdispersion parameter $\theta$ can be estimated robustly using genome-wide data, independently of $N$.

1. **Global Estimate:** Calculate a robust, genome-wide empirical estimate of the overdispersion, denoted as $\hat{\theta}_{global}$. This is implemented as the mean $F_{IS}$ across all loci (e.g., using an approach analogous to Weir and Cockerham's estimator).

2. **Uncertainty Estimate:** Crucially, we quantify the uncertainty of this estimate, $\sigma^2_{\theta}$. This is implemented as the squared standard error of the mean $F_{IS}$ (SEM$^2$).

### 4.2. Regularisation Strategy: Informed L2 (Bayesian Soft Constraint)

We modify the MDE objective function by adding a penalty term that discourages the hypothesised $\theta_{try}$ from straying too far from the empirical evidence $\hat{\theta}_{global}$.

We utilise an "Informed L2" penalty, derived from a Bayesian perspective. We conceptualise this as imposing a Gaussian prior on $\theta$, centered at the empirical estimate and scaled by its uncertainty:

$$
P(\theta) \propto \exp\left(-\frac{(\theta - \hat{\theta}_{global})^2}{2\sigma^2_{\theta}}\right)
$$

The MDE objective function incorporates the penalty (which corresponds to the negative log-prior) as follows:

$$
D_{Penalized} = D_{Metric} + \frac{(\theta_{try} - \hat{\theta}_{global})^2}{2\sigma^2_{\theta}}
$$

(Where $D_{Metric}$ is either $D_W$ or $D_{BHD}$).

*Note on Scaling:* While this approach commendably removes the need for an arbitrary tuning parameter ($\lambda$), it combines an MDE metric with a term derived from a log-prior. MDE metrics do not always scale identically to log-likelihoods. However, the inclusion of the data-driven uncertainty term $\sigma^2_{\theta}$ adaptively handles the weighting, making the approach robust.

### 4.3. Advantages of the Informed Approach

This approach is superior to standard regularisation because it is adaptive and data-driven:

* **High Confidence Data:** If $\sigma^2_{\theta}$ is small (many loci, consistent $F_{IS}$), the penalty is strong. This tightly constrains $\theta_{try}$ near $\hat{\theta}_{global}$, forcing the model to explain the remaining variance via $N$, thus stabilising the estimate.

* **Low Confidence Data:** If $\sigma^2_{\theta}$ is large, the penalty is weak, allowing the MDE framework flexibility.

## 5. Conceptual Coherence and Assumptions

### 5.1. Advantages and Coherence

* **Unified Modeling:** The BHG distribution correctly and simultaneously models the underlying biological (structure) and sampling (finite population, without replacement) processes.

* **Joint Estimation:** It eliminates the circularity, reliance on idealized models (like Wright's Island Model), and the problematic assumption that $N_c \approx N_e$ inherent in the two-step calibration. It correctly models sampling variance based on $N_c$ (or $N$).

* **Empirically Informed Stabilization:** The Informed L2 regularisation uses the richness of the genomic data to adaptively stabilise the estimation process.

### 5.2. Population Genetic and Biological Assumptions

* **The Beta Distribution Assumption:** The framework assumes that the distribution of allele frequencies across the structured population is well-approximated by a Beta distribution. While flexible, this may not capture all forms of complex population structure (e.g., recent admixture).

* **Interpretation of** $\theta$**:** We assume that the empirically estimated $\hat{\theta}_{global}$ (e.g., $F_{IS}$) is a suitable measure of the overdispersion parameter driving the BHG process. $F_{IS}$ conflates various processes (Wahlund effect, true inbreeding), which the model simplifies into a single parameter.

* **Uniform Dispersion:** The model assumes that the level of overdispersion ($\theta$) is relatively uniform across the loci used. Standard SNP filtering practices are essential.

* **Closed Population:** The method assumes the population is closed during the sampling period.

* **Random Sampling and Equal Detectability:** It assumes that individuals are sampled randomly with respect to their genotype and that all individuals have an equal probability of being sampled.
