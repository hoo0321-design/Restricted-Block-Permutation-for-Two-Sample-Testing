## Restricted-Block-Permutation-for-Two-Sample-Testing

# Restricted Block Permutation Tests  
### Quantile-Block + Complementary Block-Pair Scheme for Fast & Valid Two-Sample Testing

This repository provides clean and reproducible R implementations of the **restricted permutation framework** introduced in our paper:

> **“Restricted Block Permutation for Two-Sample Testing” (2025)**  
> A structured permutation scheme using kernel-based quantile blocks and complementary block-pair swaps, achieving exact validity with significantly reduced variance.

## 📁 Repository Structure/mean_diff/
mean_diff.R # 1D mean-difference test
MMD.R # Multivariate MMD² test (Gaussian kernel)
README.md


### 1. `mean_diff.R`
A minimal 1D example showing:

- Quantile-block construction on pooled values  
- Complementary block-pair permutation generator  
- Generalized permutation p-value using σₘ ∘ σ₀⁻¹  
- Type I & power comparison with classical full-relabel permutation

### 2. `MMD.R`
A multivariate extension using:

- Gaussian RBF kernel  
- Unbiased MMD² (U-statistic)  
- Kernel-mean score to form quantile blocks  
- Complementary block-pair restricted permutations  

Includes an end-to-end demo using 10-dimensional Gaussian data.

---

## 🔍 Key Idea (Short Summary)

### 1. Quantile Blocks  
Scores are computed on the pooled samples:

- In 1D: raw values  
- In multivariate: kernel-mean score  

Then scores are cut into *B quantiles* → block labels 1..B.

### 2. Complementary Block-Pair Swapping  
Blocks are paired as:

\[
(1 \leftrightarrow B), (2 \leftrightarrow B-1), \dots
\]

Permutations are generated only via **cross-swaps between paired blocks**.

### 3. Generalized Permutation p-value  
Use i.i.d. restricted permutations:

\[
\sigma_0, \dots, \sigma_M \sim q_X
\]

Compare statistics under σₘ ∘ σ₀⁻¹ with +1 correction:

\[
p = \frac{1 + \#\{m : T_m \ge T_{\text{obs}}\}}{1 + M}.
\]

Provides **exact finite-sample validity** without doubling.

---

## ▶️ How to Run

```r
source("mean_diff/mean_diff_quantile_block.R")
source("mmd/mmd_quantile_block_complement.R")
