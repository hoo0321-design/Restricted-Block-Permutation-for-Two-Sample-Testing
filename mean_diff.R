# ============================================================================
# 1D Mean Difference Test
# Quantile-Block Restricted Permutations (Complementary Block Pairing) vs Classical
#
# This script implements a simple 1D two-sample mean-difference test:
#   T(X, Y) = | mean(X) - mean(Y) |
#
# We compare:
#   (i)  Classical full-relabeling permutation test
#   (ii) Restricted permutation test based on
#        quantile-defined blocks and complementary block-pair swaps.
#
# The restricted scheme:
#   - Pools all observations and forms quantile blocks on the 1D values.
#   - Defines complementary block pairs (1 <-> B, 2 <-> B-1, ...) over B blocks.
#   - Generates random permutations by repeatedly swapping indices between
#     complementary blocks, with disjoint swaps.
#   - Follows the Ramdas et al. "generalized permutation test" framework:
#       σ_0, ..., σ_M ~ i.i.d. q_X
#       Compare T(σ_m ◦ σ_0^{-1} · Z) to T(Z) with +1 correction.
#
# This file is intended as a minimal 1D illustration of the block-restricted
# permutation idea used more generally for multivariate statistics (e.g., MMD^2).
# ============================================================================


# ===================== Utilities ============================================

# Sample n values (with replacement) from a numeric vector
sample_vals <- function(v, n) {
  v[sample.int(length(v), n, replace = TRUE)]
}

# Absolute mean-difference statistic: |mean(X) - mean(Y)|
mean_diff_1d <- function(x, y) {
  abs(mean(x) - mean(y))
}

# Quantile-based blocks on 1D pooled values.
# Returns an integer vector of block labels in {1, ..., num_blocks}.
make_quantile_blocks_1d <- function(z, num_blocks = 3) {
  z <- as.numeric(z)
  qs <- stats::quantile(
    z,
    probs = seq(0, 1, length.out = num_blocks + 1),
    type  = 7
  )
  findInterval(z, qs, rightmost.closed = TRUE, all.inside = TRUE)
}


# ===================== Complementary Block-Pair σ Generator ==================

# draw_sigma_complement_pairs:
#   - block: integer vector of block labels in {1, ..., B}
#   - k_swaps:
#       * if k_swaps < 1, interpreted as a fraction of N (e.g., 0.25 -> ~0.25N swaps)
#       * if k_swaps >= 1, interpreted as an absolute number of swaps
#   - max_tries: safety cap to avoid infinite loops
#
# The scheme:
#   - Define complementary pairs (1 <-> B), (2 <-> B-1), ...
#   - At each step, randomly pick one complementary pair (u, v),
#     then randomly choose one unused index from block u and one from block v,
#     and swap their positions in σ.
#   - Each index can be swapped at most once (disjoint swaps).
#
# This induces a restricted permutation distribution q_X supported on permutations
# reachable via cross-swaps between complementary blocks.
draw_sigma_complement_pairs <- function(block,
                                        k_swaps   = 0.2,
                                        max_tries = 8000L) {
  N <- length(block)
  B <- max(block)
  if (B < 2L) {
    # No pairing possible when there is only a single block
    return(seq_len(N))
  }
  
  # Convert k_swaps to an integer count
  ks <- if (k_swaps < 1) {
    max(1L, round(k_swaps * N))
  } else {
    as.integer(k_swaps)
  }
  
  # Complementary block pairs: (1 <-> B), (2 <-> B-1), ...
  pair_of <- function(b) B + 1L - b
  pairs <- lapply(seq_len(floor(B / 2)), function(b) c(b, pair_of(b)))
  
  # Indices belonging to each block
  idx_by_block <- lapply(seq_len(B), function(b) which(block == b))
  
  # Start from the identity permutation
  sigma <- seq_len(N)
  used  <- rep(FALSE, N)  # track indices already swapped
  swaps <- 0L
  tries <- 0L
  
  while (swaps < ks && tries < max_tries) {
    tries <- tries + 1L
    
    # 1) Randomly choose one complementary pair (u, v)
    pr <- sample(pairs, 1L)[[1]]
    u  <- pr[1]
    v  <- pr[2]
    
    # 2) Candidate indices from each block that have not yet been swapped
    cand_u <- idx_by_block[[u]][!used[idx_by_block[[u]]]]
    cand_v <- idx_by_block[[v]][!used[idx_by_block[[v]]]]
    if (!length(cand_u) || !length(cand_v)) {
      next
    }
    
    # 3) Sample one index from each block, and swap them
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


# ===================== Restricted p-value (Quantile Blocks) ==================

# pval_qblock_mean1d:
#   - X, Y: numeric vectors
#   - num_blocks: number of quantile blocks
#   - k_swaps: controls the number (or fraction) of cross-swaps
#   - M: number of restricted permutations (excluding σ0)
#
# Implementation follows Ramdas-style generalized permutation testing:
#   1) Draw σ0, ..., σM i.i.d. from q_X (our restricted permutation law).
#   2) Fix σ0 and compare T(σm ◦ σ0^{-1} · Z) with T(Z) for m = 1..M.
#   3) Use (1 + #extreme)/(1 + M) as the p-value (with +1 correction).
pval_qblock_mean1d <- function(X, Y,
                               num_blocks = 3,
                               k_swaps    = 0.2,
                               M          = 2000) {
  X <- as.numeric(X)
  Y <- as.numeric(Y)
  n1 <- length(X)
  N  <- n1 + length(Y)
  pooled <- c(X, Y)
  
  # 1) Observed statistic
  T_obs <- mean_diff_1d(X, Y)
  
  # 2) Label-free blocks on the pooled values
  block <- make_quantile_blocks_1d(pooled, num_blocks = num_blocks)
  
  # 3) σ0, σ1, ..., σM ~ i.i.d. from q_X (complement-pair generator)
  sigmas <- replicate(
    M + 1L,
    draw_sigma_complement_pairs(block, k_swaps = k_swaps),
    simplify = FALSE
  )
  sigma0_inv <- order(sigmas[[1L]])  # σ0^{-1}
  
  # 4) Compare T(σm ◦ σ0^{-1} · Z) with T_obs
  ge <- 0L
  for (m in 2:(M + 1L)) {
    comp <- sigmas[[m]][sigma0_inv]  # composition σm ◦ σ0^{-1}
    Zp   <- pooled[comp]
    Xp   <- Zp[1:n1]
    Yp   <- Zp[(n1 + 1):N]
    if (mean_diff_1d(Xp, Yp) >= T_obs) {
      ge <- ge + 1L
    }
  }
  
  (1 + ge) / (1 + M)
}


# ===================== Type I / Power Wrappers (Restricted) ==================

# Empirical Type I error under H0: X, Y ~ pool
type1_qblock_mean1d <- function(pool,
                                n         = 64,
                                alpha     = 0.05,
                                num_blocks = 3,
                                k_swaps    = 0.2,
                                M          = 2000,
                                n_sim      = 300,
                                seed       = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rej <- 0L
  for (s in seq_len(n_sim)) {
    X <- sample_vals(pool, n)
    Y <- sample_vals(pool, n)
    p <- pval_qblock_mean1d(X, Y, num_blocks, k_swaps, M)
    if (p < alpha) {
      rej <- rej + 1L
    }
  }
  rej / n_sim
}

# Empirical power under H1: X ~ g0, Y ~ g1
power_qblock_mean1d <- function(g0, g1,
                                n         = 64,
                                alpha     = 0.05,
                                num_blocks = 3,
                                k_swaps    = 0.2,
                                M          = 2000,
                                n_sim      = 300,
                                seed       = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rej <- 0L
  for (s in seq_len(n_sim)) {
    X <- sample_vals(g0, n)
    Y <- sample_vals(g1, n)
    p <- pval_qblock_mean1d(X, Y, num_blocks, k_swaps, M)
    if (p < alpha) {
      rej <- rej + 1L
    }
  }
  rej / n_sim
}


# ===================== Classical Full-Relabeling p-value ====================

# pval_classic_mean1d:
#   Standard Monte Carlo permutation test with full relabeling.
#   Uses the same statistic T(X, Y) = |mean(X) - mean(Y)| and the same
#   +1 correction: (1 + #extreme)/(1 + n_perm).
pval_classic_mean1d <- function(X, Y, n_perm = 2000) {
  X <- as.numeric(X)
  Y <- as.numeric(Y)
  n1 <- length(X)
  N  <- n1 + length(Y)
  pooled <- c(X, Y)
  
  T_obs <- mean_diff_1d(X, Y)
  
  perm_stats <- replicate(n_perm, {
    idx <- sample.int(N)
    Zp  <- pooled[idx]
    Xp  <- Zp[1:n1]
    Yp  <- Zp[(n1 + 1):N]
    mean_diff_1d(Xp, Yp)
  })
  
  (1 + sum(perm_stats >= T_obs)) / (1 + n_perm)
}

# Empirical Type I error for the classical permutation test
type1_classic_mean1d <- function(pool,
                                 n       = 64,
                                 alpha   = 0.05,
                                 n_perm  = 2000,
                                 n_sim   = 300,
                                 seed    = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rej <- 0L
  for (s in seq_len(n_sim)) {
    X <- sample_vals(pool, n)
    Y <- sample_vals(pool, n)
    p <- pval_classic_mean1d(X, Y, n_perm)
    if (p < alpha) {
      rej <- rej + 1L
    }
  }
  rej / n_sim
}

# Empirical power for the classical permutation test
power_classic_mean1d <- function(g0, g1,
                                 n       = 64,
                                 alpha   = 0.05,
                                 n_perm  = 2000,
                                 n_sim   = 300,
                                 seed    = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rej <- 0L
  for (s in seq_len(n_sim)) {
    X <- sample_vals(g0, n)
    Y <- sample_vals(g1, n)
    p <- pval_classic_mean1d(X, Y, n_perm)
    if (p < alpha) {
      rej <- rej + 1L
    }
  }
  rej / n_sim
}


# ===================== Example: Side-by-Side Comparison =====================

# This example compares Type I error and power of:
#   - Quantile-block restricted test with complementary block-pair permutations
#   - Classical full-relabeling permutation test
#
# Setting:
#   H0: X, Y ~ N(0, 1)
#   H1: X ~ N(0, 1), Y ~ N(0.4, 1) (mean shift)

set.seed(2024)

# Null and alternative pools
pool <- rnorm(10000, mean = 0, sd = 1)  # for H0
g0   <- rnorm(10000, mean = 0, sd = 1)
g1   <- rnorm(10000, mean = 0.4, sd = 1)

# Experiment configuration
n              <- 64
alpha          <- 0.05
num_blk        <- 3
k_sw           <- 0.2   # ~20% of indices cross-swapped between complementary blocks
M              <- 100   # number of restricted permutations
n_sim          <- 100   # simulation runs for Type I / power
n_perm_classic <- 100   # permutations for classical test

# Estimated Type I error (H0)
t1_q  <- type1_qblock_mean1d(pool, n, alpha,
                             num_blocks = num_blk,
                             k_swaps    = k_sw,
                             M          = M,
                             n_sim      = n_sim,
                             seed       = 1)

t1_cl <- type1_classic_mean1d(pool, n, alpha,
                              n_perm = n_perm_classic,
                              n_sim  = n_sim,
                              seed   = 1)

# Estimated power (H1)
pw_q  <- power_qblock_mean1d(g0, g1, n, alpha,
                             num_blocks = num_blk,
                             k_swaps    = k_sw,
                             M          = M,
                             n_sim      = n_sim,
                             seed       = 2)

pw_cl <- power_classic_mean1d(g0, g1, n, alpha,
                              n_perm = n_perm_classic,
                              n_sim  = n_sim,
                              seed   = 2)

cat(sprintf("[Quantile-Block (complement-pair)] Type I ≈ %.4f, Power ≈ %.4f\n",
            t1_q, pw_q))
cat(sprintf("[Classical full relabeling]       Type I ≈ %.4f, Power ≈ %.4f\n",
            t1_cl, pw_cl))

