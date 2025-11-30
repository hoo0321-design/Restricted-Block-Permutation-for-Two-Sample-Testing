# ============================================================================
# Multivariate Two-Sample Test with Unbiased MMD^2
# Quantile-Block Restricted Permutations (Complementary Block Pairing) vs Classical
#
# This script implements a multivariate two-sample test based on the unbiased
# squared Maximum Mean Discrepancy (MMD^2) with a Gaussian RBF kernel:
#
#   T(X, Y) = MMD^2_unbiased(X, Y; k_gamma)
#
# We compare:
#   (i)  Classical full-relabeling permutation test
#   (ii) Restricted permutation test based on:
#        - Kernel-mean score based quantile blocks on the pooled data
#        - Complementary block-pair swaps (1 <-> B, 2 <-> B-1, ...)
#
# The restricted scheme:
#   - Pools all observations Z = (X, Y) and forms blocks via quantiles of
#     the kernel mean score s_i = (1/N) sum_j k(z_i, z_j).
#   - Defines complementary block pairs (1 <-> B, 2 <-> B-1, ...) over B blocks.
#   - Generates random permutations by repeatedly swapping indices between
#     complementary blocks, with disjoint swaps.
#   - Follows the Ramdas-style "generalized permutation test" framework:
#       σ_0, ..., σ_M ~ i.i.d. q_X
#       Compare T(σ_m ◦ σ_0^{-1} · Z) to T(Z) with +1 correction:
#         p = (1 + #extreme) / (1 + M).
#
# This file is meant as a simple multivariate example (d = 10 Gaussian)
# illustrating the same restricted permutation idea used in the paper.
# ============================================================================


# ===================== Utilities ============================================

library(MASS)  # for mvrnorm

# Ensure sample-by-feature matrix format
as_sample_matrix <- function(x) {
  if (is.vector(x)) {
    matrix(x, ncol = 1L)
  } else {
    as.matrix(x)
  }
}

# Row sampler that works safely with matrices
sample_rows <- function(M, n) {
  M <- as_sample_matrix(M)
  M[sample.int(nrow(M), n, replace = TRUE), , drop = FALSE]
}

# Pairwise squared Euclidean distances
pairwise_sq_dists <- function(X, Y = NULL) {
  X <- as.matrix(X)
  if (is.null(Y)) {
    Y <- X
  } else {
    Y <- as.matrix(Y)
  }
  Xn <- rowSums(X * X)
  Yn <- rowSums(Y * Y)
  outer(Xn, Yn, "+") - 2 * (X %*% t(Y))
}

# RBF (Gaussian) kernel: k(x,y) = exp(-gamma * ||x - y||^2)
rbf_kernel_matrix <- function(X, Y = NULL, gamma = 1.0) {
  D2 <- pairwise_sq_dists(X, Y)
  exp(-gamma * D2)
}

# Median heuristic for gamma using pooled data
# gamma = 1 / (2 * median(distance^2))
median_heuristic_gamma <- function(Z) {
  D2 <- pairwise_sq_dists(Z)
  d2 <- D2[upper.tri(D2, diag = FALSE)]
  d2 <- d2[d2 > 0]
  med <- if (length(d2)) stats::median(d2) else 1.0
  if (!is.finite(med) || med <= 0) med <- 1.0
  1.0 / (2.0 * med)
}

# Unbiased MMD^2 (U-statistic estimator)
mmd2_unbiased <- function(X, Y, gamma) {
  X <- as_sample_matrix(X)
  Y <- as_sample_matrix(Y)
  n1 <- nrow(X)
  n2 <- nrow(Y)
  if (n1 < 2 || n2 < 2) {
    stop("Need at least 2 samples per group for unbiased MMD^2")
  }
  
  Kxx <- rbf_kernel_matrix(X, X, gamma)
  Kyy <- rbf_kernel_matrix(Y, Y, gamma)
  Kxy <- rbf_kernel_matrix(X, Y, gamma)
  
  sum_Kxx_off <- sum(Kxx) - sum(diag(Kxx))
  sum_Kyy_off <- sum(Kyy) - sum(diag(Kyy))
  
  term_xx <- sum_Kxx_off / (n1 * (n1 - 1))
  term_yy <- sum_Kyy_off / (n2 * (n2 - 1))
  term_xy <- (2.0 / (n1 * n2)) * sum(Kxy)
  
  term_xx + term_yy - term_xy
}


# ===================== Kernel-Mean Quantile Blocks ==========================

# make_quantile_blocks_kernel_mean:
#   - pooled: matrix of pooled samples Z (N x d)
#   - num_blocks: desired number of blocks B
#   - gamma: RBF kernel parameter (if NULL, uses median heuristic)
#
# We:
#   1) Compute the Gram matrix K(Z, Z).
#   2) Define a kernel-mean score s_i = mean_j K_ij.
#   3) Cut the scores into B quantiles and assign block labels 1..B.
make_quantile_blocks_kernel_mean <- function(pooled,
                                             num_blocks = 3,
                                             gamma      = NULL) {
  Z <- as_sample_matrix(pooled)
  if (is.null(gamma)) {
    gamma <- median_heuristic_gamma(Z)
  }
  K_all <- rbf_kernel_matrix(Z, Z, gamma)
  score <- rowMeans(K_all)
  qs <- stats::quantile(
    score,
    probs = seq(0, 1, length.out = num_blocks + 1),
    type  = 7
  )
  findInterval(score, qs, rightmost.closed = TRUE, all.inside = TRUE)
}


# ===================== Complementary Block-Pair σ Generator ==================

# draw_sigma_from_qX_complement:
#   - block: integer vector of block labels in {1, ..., B}
#   - k_swaps:
#       * if k_swaps < 1, interpreted as a fraction of N (e.g., 0.25 -> ~0.25N swaps)
#       * if k_swaps >= 1, interpreted as an absolute number of swaps
#   - max_tries: safety cap to avoid infinite loops
#
# Scheme:
#   - Define complementary block pairs (1 <-> B), (2 <-> B-1), ...
#   - At each step, pick a pair (u, v), then randomly choose one unused index
#     from block u and one from block v, and swap their positions in σ.
#   - Each index can be swapped at most once (disjoint swaps).
#
# This induces a restricted permutation distribution q_X based on cross-swaps
# between complementary blocks.
draw_sigma_from_qX_complement <- function(block,
                                          k_swaps   = 0.2,
                                          max_tries = 5000L) {
  N <- length(block)
  B <- max(block)
  if (B < 2L) {
    # No pairing possible with a single block
    return(seq_len(N))
  }
  
  # Convert k_swaps to an integer swap count
  ks <- if (k_swaps < 1) {
    max(1L, round(k_swaps * N))
  } else {
    as.integer(k_swaps)
  }
  
  # Complementary pairs: (1 <-> B), (2 <-> B-1), ...
  pair_of <- function(b) B + 1L - b
  pairs <- lapply(seq_len(floor(B / 2)), function(b) c(b, pair_of(b)))
  
  # Indices in each block
  idx_by_block <- lapply(seq_len(B), function(b) which(block == b))
  
  # Start from identity permutation
  sigma <- seq_len(N)
  used  <- rep(FALSE, N)
  swaps <- 0L
  tries <- 0L
  
  while (swaps < ks && tries < max_tries) {
    tries <- tries + 1L
    
    # Pick a complementary block pair (u, v)
    pr <- sample(pairs, 1L)[[1]]
    u  <- pr[1]
    v  <- pr[2]
    
    # Candidate indices that have not been swapped yet
    cand_u <- idx_by_block[[u]][!used[idx_by_block[[u]]]]
    cand_v <- idx_by_block[[v]][!used[idx_by_block[[v]]]]
    if (!length(cand_u) || !length(cand_v)) {
      next
    }
    
    # Swap one index from each block
    i <- sample(cand_u, 1L)
    j <- sample(cand_v, 1L)
    if (i != j) {
      tmp        <- sigma[i]
      sigma[i]   <- sigma[j]
      sigma[j]   <- tmp
      used[c(i, j)] <- TRUE
      swaps <- swaps + 1L
    }
  }
  
  sigma
}


# ===================== Restricted p-value (Kernel-MMD) ======================

# p_value_quantile_blocks_kernel_raw:
#   - X, Y: matrices of samples (n1 x d, n2 x d)
#   - num_blocks: number of quantile blocks B
#   - k_swaps: controls number (or fraction) of cross-swaps
#   - M: number of restricted permutations (excluding σ0)
#
# Follows generalized permutation testing (Ramdas-style):
#   1) Compute observed T_obs = MMD^2_unbiased(X, Y; gamma).
#   2) Construct blocks on pooled Z via kernel-mean quantiles.
#   3) Draw σ0, ..., σM i.i.d. from q_X (complement-pair generator).
#   4) Use composition σm ◦ σ0^{-1} for the reference distribution and
#      compute T_m = T(σm ◦ σ0^{-1} · Z). Compare to T_obs.
#   5) p = (1 + # {m : T_m >= T_obs} ) / (1 + M).
p_value_quantile_blocks_kernel_raw <- function(X, Y,
                                               num_blocks = 3,
                                               k_swaps    = 0.2,
                                               M          = 2000) {
  X <- as_sample_matrix(X)
  Y <- as_sample_matrix(Y)
  n1 <- nrow(X)
  n2 <- nrow(Y)
  N  <- n1 + n2
  
  pooled <- rbind(X, Y)
  
  gamma <- median_heuristic_gamma(pooled)
  T_obs <- mmd2_unbiased(X, Y, gamma)
  
  block <- make_quantile_blocks_kernel_mean(
    pooled,
    num_blocks = num_blocks,
    gamma      = gamma
  )
  
  # σ0, ..., σM ~ i.i.d. from q_X
  sigmas <- replicate(
    M + 1L,
    draw_sigma_from_qX_complement(block, k_swaps = k_swaps),
    simplify = FALSE
  )
  sigma0_inv <- order(sigmas[[1L]])  # σ0^{-1}
  
  ge <- 0L
  for (m in 2:(M + 1L)) {
    comp <- sigmas[[m]][sigma0_inv]  # σm ◦ σ0^{-1}
    Zp   <- pooled[comp, , drop = FALSE]
    Xp   <- Zp[1:n1, , drop = FALSE]
    Yp   <- Zp[(n1 + 1):N, , drop = FALSE]
    Tm   <- mmd2_unbiased(Xp, Yp, gamma)
    if (Tm >= T_obs) {
      ge <- ge + 1L
    }
  }
  
  (1 + ge) / (1 + M)
}


# ===================== Type I / Power Wrappers (Restricted) =================

# Empirical Type I error (H0) for restricted kernel-MMD test
type1_mmd_kernel_blocks_raw <- function(pool,
                                        n          = 100,
                                        alpha      = 0.05,
                                        num_blocks = 3,
                                        k_swaps    = 0.2,
                                        M          = 2000,
                                        n_sim      = 200,
                                        seed       = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rej <- 0L
  for (s in seq_len(n_sim)) {
    X <- sample_rows(pool, n)
    Y <- sample_rows(pool, n)
    p <- p_value_quantile_blocks_kernel_raw(X, Y, num_blocks, k_swaps, M)
    if (p < alpha) {
      rej <- rej + 1L
    }
  }
  rej / n_sim
}

# Empirical power (H1) for restricted kernel-MMD test
power_mmd_kernel_blocks_raw <- function(g0, g1,
                                        n          = 100,
                                        alpha      = 0.05,
                                        num_blocks = 3,
                                        k_swaps    = 0.2,
                                        M          = 2000,
                                        n_sim      = 200,
                                        seed       = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rej <- 0L
  for (s in seq_len(n_sim)) {
    X <- sample_rows(g0, n)
    Y <- sample_rows(g1, n)
    p <- p_value_quantile_blocks_kernel_raw(X, Y, num_blocks, k_swaps, M)
    if (p < alpha) {
      rej <- rej + 1L
    }
  }
  rej / n_sim
}


# ===================== Classical Full-Relabeling MMD Test ===================

# mmd_power_classical_rows:
#   - group1, group2: matrices of samples from two distributions
#   - n: per-group sample size
#   - n_perm: number of full label permutations
#   - n_simulations: number of Monte Carlo runs
#
# For each run:
#   1) Draw X, Y of size n from group1, group2.
#   2) Compute observed MMD^2_unbiased(X, Y).
#   3) Generate n_perm full relabeling permutations and recompute MMD^2.
#   4) Compute classical permutation p-value with +1 correction.
mmd_power_classical_rows <- function(group1, group2,
                                     n             = 50,
                                     alpha         = 0.05,
                                     n_perm        = 2000,
                                     n_simulations = 300,
                                     seed          = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rej <- 0L
  for (sim in seq_len(n_simulations)) {
    X <- sample_rows(group1, n)
    Y <- sample_rows(group2, n)
    pooled <- rbind(X, Y)
    gamma  <- median_heuristic_gamma(pooled)
    obs    <- mmd2_unbiased(X, Y, gamma)
    N      <- 2L * n
    
    perm_stats <- replicate(n_perm, {
      idx <- sample.int(N)
      Zp  <- pooled[idx, , drop = FALSE]
      Xp  <- Zp[1:n, , drop = FALSE]
      Yp  <- Zp[(n + 1):N, , drop = FALSE]
      mmd2_unbiased(Xp, Yp, gamma)
    })
    
    pval <- (1 + sum(perm_stats >= obs)) / (1 + length(perm_stats))
    if (pval < alpha) {
      rej <- rej + 1L
    }
  }
  rej / n_simulations
}


# ===================== d = 10 Example: Type I & Power =======================

set.seed(1007)

# --- Data Generation (d = 10 Gaussian) ---
d     <- 10
Sigma <- diag(d)

# Null (Type I): both from N(0, I_d)
pool10 <- MASS::mvrnorm(40000, mu = rep(0, d), Sigma = Sigma)

# Alternative (Power): mean shift on the first coordinate
G0_10 <- MASS::mvrnorm(40000, mu = rep(0, d),           Sigma = Sigma)
G1_10 <- MASS::mvrnorm(40000, mu = c(0.4, rep(0, d - 1)), Sigma = Sigma)  # effect size can be adjusted

# --- Common settings ---
alpha   <- 0.05
M_perm  <- 100      # number of restricted permutations
n_sim   <- 100


# --- 1) Type I error comparison (larger n for stability) ---

n       <- 256      # per-group sample size
num_blk <- 5
k_swaps <- 0.2
n_perm_classic <- 100

type1_classic_10 <- mmd_power_classical_rows(
  pool10, pool10,
  n             = n,
  alpha         = alpha,
  n_perm        = n_perm_classic,
  n_simulations = n_sim,
  seed          = 1
)

type1_block_10 <- type1_mmd_kernel_blocks_raw(
  pool    = pool10,
  n       = n,
  alpha   = alpha,
  num_blocks = num_blk,
  k_swaps    = k_swaps,
  M          = M_perm,
  n_sim      = n_sim,
  seed       = 2
)

cat(sprintf("\n[d=10] Type I — classic: %.4f, block: %.4f (alpha = %.2f)\n",
            type1_classic_10, type1_block_10, alpha))


# --- 2) Power comparison (smaller n) ---

n       <- 128
num_blk <- 4
k_swaps <- 0.2
n_perm_classic <- 100

power_classic_10 <- mmd_power_classical_rows(
  G0_10, G1_10,
  n             = n,
  alpha         = alpha,
  n_perm        = n_perm_classic,
  n_simulations = n_sim,
  seed          = 3
)

power_block_10 <- power_mmd_kernel_blocks_raw(
  g0       = G0_10,
  g1       = G1_10,
  n        = n,
  alpha    = alpha,
  num_blocks = num_blk,
  k_swaps    = k_swaps,
  M          = M_perm,
  n_sim      = n_sim,
  seed       = 4
)

cat(sprintf("[d=10] Power  — classic: %.4f, block: %.4f\n",
            power_classic_10, power_block_10))

