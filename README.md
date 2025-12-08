# dotR: A Unified Framework for Estimating Census Size

**dotR** provides a method for estimating census population size ($N_c$) in structured populations. It implements a **3D Regularised Minimum Distance Estimation (MDE)** framework using the **Beta-Hypergeometric (BHG)** distribution.

## Overview

Original methods for estimating $N_c$ (like `dnadot`) rely on the Hypergeometric distribution, which assumes a panmictic population. When population structure exists, these methods suffer from bias due to overdispersion. 

**dotR** addresses this by simultaneously accounting for:
1.  **Finite Sampling without Replacement:** Necessary for estimating $N_c$.
2.  **Population Structure:** Modeled via the Beta distribution to handle overdispersion ($\theta$ or $F_{IS}$).

## Key Features

* **Beta-Hypergeometric (BHG) Model:** Correctly models sampling variance in finite, structured populations without relying on infinite population assumptions (like Beta-Binomial).
* **3D MDE Framework:** Jointly estimates Total Gene Copies ($N$), Mean Allele Frequency ($p$), and Overdispersion ($\theta$).
* **Informed L2 Regularisation:** Uses genome-wide empirical $F_{IS}$ estimates to stabilise parameters, resolving identifiability issues between structural variance and sampling variance.
* **Robust Metrics:** Includes Wasserstein Distance and Buffered Histogram Discrepancy for robust estimation.

## Methodology

For a deep dive into the mathematical framework, probability mass functions, and the 3D regularization algorithms used in this package, please see the full technical documentation:

[**📄 Read the Full Methodological Framework (METHODS.md)**](./METHODS.md)

## Installation

```r
# install.packages("devtools")
# devtools::install_github("mfr-public/dotR")
