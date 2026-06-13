# Diagnostic script to check if x variance shrinks during MCMC
# This helps identify whether the issue is in MALA (x sampling) or Gibbs (rho/kappa)

# Load package
library(spatiocompo)
library(Matrix)

set.seed(42)

# Small grid for fast testing
n_side <- 5
n <- n_side^2
d <- 2
D <- 3

# True parameters (realistic alpha)
kappa_true <- 0.5
rho_true <- matrix(c(1.0, 0.3, 0.3, 1.0), 2, 2)
alpha_true <- 5

# Create grid and precision
G <- create_Q(c(n_side, n_side), alpha = 1)
Q_marg <- Q_rhoxQ(NULL, kappa_true, G)$Q
rhoxQ_true <- Q_rhoxQ(rho_true, kappa_true, G)$rhoxQ

# Sample true x from GMRF prior
# x ~ N(0, (rhoxQ)^{-1})
# Use Cholesky: if L'L = rhoxQ, then x = L^{-1} z where z ~ N(0,I)
R_Q <- chol(rhoxQ_true)
z_rand <- rnorm(n * d)
x_true <- backsolve(R_Q, z_rand)
x_true_mat <- matrix(x_true, nrow = n, ncol = d)

cat("True x variance per component:\n")
cat("  Var(x[,1]) =", var(x_true_mat[,1]), "\n")
cat("  Var(x[,2]) =", var(x_true_mat[,2]), "\n")
cat("  Expected (approx 1/diag(Q) * rho_diag):", 1/mean(Matrix::diag(Q_marg)) * rho_true[1,1], "\n\n")

# Transform to compositional space and generate observations
z_true <- invalr(x_true_mat)
y <- matrix(0, n, D)
for (i in 1:n) {
  # Sample from Dirichlet(alpha * z_true[i,])
  y[i,] <- rdirichlet(1, alpha_true * z_true[i,])
}

# Setup for testing
A <- Matrix::Diagonal(n * d)
priors <- create_priors(G, d = 2)
priors$logstep.MALA <- log(0.1)  # Start with smaller step

# Test 1: Run MALA steps with FIXED kappa and rho (with step size adaptation)
cat("=== Test 1: MALA with fixed kappa, rho ===\n")
x_current <- matrix(rnorm(n * d, sd = 0.5), n, d)
alpha_current <- 10
logstep_current <- priors$logstep.MALA

n_iter <- 100
var_history <- matrix(0, n_iter, 2)
accept_count <- 0

for (i in 1:n_iter) {
  rhoxQ <- Q_rhoxQ(rho_true, kappa_true, G)$rhoxQ

  mala_result <- mala_step(
    x = as.vector(x_current),
    alpha = alpha_current,
    y = y,
    w = rep(1, n),
    Q = rhoxQ,
    A = A,
    priors = priors,
    d = d,
    logstep0 = logstep_current,  # Use adapted step
    iter = i
  )

  new <- mala_result$sample
  alpha_current <- new[length(new)]
  x_current <- matrix(new[1:(n*d)], n, d)
  accept_count <- accept_count + mala_result$count
  logstep_current <- mala_result$logstep  # Update step size

  var_history[i,1] <- var(x_current[,1])
  var_history[i,2] <- var(x_current[,2])
}

cat("  MALA acceptance rate:", accept_count/n_iter, "\n")
cat("  Initial var(x):", var_history[1,1], var_history[1,2], "\n")
cat("  Final var(x):  ", var_history[n_iter,1], var_history[n_iter,2], "\n")
cat("  Var change:    ", var_history[n_iter,1]/var_history[1,1], var_history[n_iter,2]/var_history[1,2], "\n")
cat("  Final log step:", logstep_current, "\n\n")

# Test 2: Check Fisher Information eigenvalues
cat("=== Test 2: Fisher Information structure ===\n")
x_test <- matrix(rnorm(n * d, sd = 0.5), n, d)
fi_result <- fisher_info(
  x = as.vector(x_test),
  alpha = 10,
  A = A,
  y = y,
  b = priors$alpha.a,
  c = priors$alpha.b,
  Q = Q_rhoxQ(rho_true, kappa_true, G)$rhoxQ,
  step = 0.1,
  d = d,
  w = rep(1, n)
)

RFI <- fi_result$RFI
FI_reconstructed <- t(RFI) %*% RFI

cat("  FI size:", nrow(RFI), "x", ncol(RFI), "\n")
cat("  RFI diag range:", range(diag(RFI)), "\n")
cat("  FI eigenvalue range:", range(eigen(FI_reconstructed, only.values=TRUE)$values), "\n")

# Check proposal variance
proposal_var <- solve(FI_reconstructed)
cat("  Proposal variance diag range:", range(diag(proposal_var)), "\n\n")

# Test 3: Check x'Qx symmetry
cat("=== Test 3: x'Qx computation ===\n")
x_mat <- matrix(rnorm(n * d), n, d)
xQx <- t(x_mat) %*% as.matrix(Q_marg) %*% x_mat
cat("  x'Qx:\n")
print(xQx)
cat("  Symmetric?", isSymmetric(xQx), "\n")
cat("  Off-diagonal difference:", abs(xQx[1,2] - xQx[2,1]), "\n\n")

# Test 4: Compare numerical and analytical gradient
cat("=== Test 4: Gradient check ===\n")
x_vec <- rnorm(n * d, sd = 0.5)
alpha_test <- 8
rhoxQ <- Q_rhoxQ(rho_true, kappa_true, G)$rhoxQ

result <- neg_log_lik(x_vec, alpha_test, A, y,
                      b = priors$alpha.a, c = priors$alpha.b,
                      Q = rhoxQ, d = d, w = rep(1, n))

# Numerical gradient for first few components
eps <- 1e-5
num_grad <- numeric(5)
for (j in 1:5) {
  x_plus <- x_vec; x_plus[j] <- x_plus[j] + eps
  x_minus <- x_vec; x_minus[j] <- x_minus[j] - eps
  L_plus <- neg_log_lik(x_plus, alpha_test, A, y, priors$alpha.a, priors$alpha.b, rhoxQ, d, rep(1,n))$L
  L_minus <- neg_log_lik(x_minus, alpha_test, A, y, priors$alpha.a, priors$alpha.b, rhoxQ, d, rep(1,n))$L
  num_grad[j] <- (L_plus - L_minus) / (2 * eps)
}

cat("  Analytical gradient (first 5):", round(result$dL[1:5], 4), "\n")
cat("  Numerical gradient (first 5): ", round(num_grad, 4), "\n")
cat("  Max difference:", max(abs(result$dL[1:5] - num_grad)), "\n\n")

# Test 5: Run short MCMC with full sampling
cat("=== Test 5: Short MCMC run ===\n")
theta0 <- list(
  alpha = 10,
  kappa = 0.3,
  rho = diag(2),
  x = matrix(0, n, d)
)

fit <- mcmc_sampling(theta0, priors, y, iter = 200, verbose = FALSE)

cat("  MALA acceptance:", fit$MALA$acc, "\n")
cat("  Gibbs acceptance:", fit$Gibbs$acc, "\n")
cat("  Alpha: true=", alpha_true, ", mean=", mean(fit$alpha[100:200]), "\n")
cat("  Kappa: true=", kappa_true, ", mean=", mean(fit$kappa[100:200]), "\n")

# Check rho recovery
rho_samples <- array(fit$rho[,100:200], dim=c(2,2,101))
rho_mean <- apply(rho_samples, c(1,2), mean)
cat("  Rho[1,1]: true=", rho_true[1,1], ", mean=", rho_mean[1,1], "\n")
cat("  Rho[2,2]: true=", rho_true[2,2], ", mean=", rho_mean[2,2], "\n")

# Check x variance over MCMC
x_var_start <- var(fit$x[,1])
x_var_end <- var(fit$x[,200])
cat("  Var(x) at iter 1:  ", x_var_start, "\n")
cat("  Var(x) at iter 200:", x_var_end, "\n")
cat("  Ratio:", x_var_end/x_var_start, "\n")

# Test 6: Start from TRUE values
cat("\n=== Test 6: MCMC starting from TRUE values ===\n")
theta0_true <- list(
  alpha = alpha_true,
  kappa = kappa_true,
  rho = rho_true,
  x = x_true_mat
)

fit_true <- mcmc_sampling(theta0_true, priors, y, iter = 500, verbose = FALSE)

cat("  MALA acceptance:", fit_true$MALA$acc, "\n")
cat("  Gibbs acceptance:", fit_true$Gibbs$acc, "\n")
cat("  Alpha: true=", alpha_true, ", mean=", mean(fit_true$alpha[250:500]), "\n")
cat("  Kappa: true=", kappa_true, ", mean=", mean(fit_true$kappa[250:500]), "\n")

rho_samples_true <- array(fit_true$rho[,250:500], dim=c(2,2,251))
rho_mean_true <- apply(rho_samples_true, c(1,2), mean)
cat("  Rho[1,1]: true=", rho_true[1,1], ", mean=", rho_mean_true[1,1], "\n")
cat("  Rho[2,2]: true=", rho_true[2,2], ", mean=", rho_mean_true[2,2], "\n")

cat("  Var(x) true:    ", var(x_true_mat[,1]), var(x_true_mat[,2]), "\n")
cat("  Var(x) at iter 1:  ", var(fit_true$x[,1]), "\n")
cat("  Var(x) at iter 500:", var(fit_true$x[,500]), "\n")

# Trace of key parameters
cat("\n  Alpha trace (every 50):", fit_true$alpha[seq(1, 500, by=50)], "\n")
cat("  Kappa trace (every 50):", fit_true$kappa[seq(1, 500, by=50)], "\n")

# Test 7: Check the priors
cat("\n=== Test 7: Prior check ===\n")
cat("  Alpha prior: Gamma(", priors$alpha.a, ",", priors$alpha.b, ")\n")
cat("    Prior mean:", priors$alpha.a / priors$alpha.b, "\n")
cat("  Kappa prior: Gamma(", priors$kappa.range, ",", priors$kappa.sigma, ")\n")
cat("    Prior mean:", priors$kappa.range / priors$kappa.sigma, "\n")
cat("  Rho prior: IW(", priors$rho.df, ", Sigma)\n")
cat("    Sigma =", priors$rho.Sigma, "\n")

# Test 8: Check x'Qx for true x and final x
cat("\n=== Test 8: x'Qx analysis ===\n")
Q_marg_final <- Q_rhoxQ(NULL, mean(fit_true$kappa[250:500]), G)$Q

# True x'Qx
xQx_true <- t(x_true_mat) %*% as.matrix(Q_marg_final) %*% x_true_mat
cat("  True x'Qx:\n")
print(round(xQx_true, 2))

# Final x (properly reshaped)
x_final <- matrix(fit_true$x[,500], n, d)
xQx_final <- t(x_final) %*% as.matrix(Q_marg_final) %*% x_final
cat("  Final x'Qx:\n")
print(round(xQx_final, 2))

cat("  True x variance per component:", var(x_true_mat[,1]), var(x_true_mat[,2]), "\n")
cat("  Final x variance per component:", var(x_final[,1]), var(x_final[,2]), "\n")

# Test 9: Longer run for convergence with variance tracking
cat("\n=== Test 9: Longer MCMC (2000 iter) with variance tracking ===\n")
fit_long <- mcmc_sampling(theta0_true, priors, y, iter = 2000, verbose = FALSE)

cat("  MALA acceptance:", fit_long$MALA$acc, "\n")
cat("  Gibbs acceptance:", fit_long$Gibbs$acc, "\n")
cat("  Alpha: true=", alpha_true, ", mean=", mean(fit_long$alpha[1000:2000]), "\n")
cat("  Kappa: true=", kappa_true, ", mean=", mean(fit_long$kappa[1000:2000]), "\n")

rho_long <- array(fit_long$rho[,1000:2000], dim=c(2,2,1001))
rho_mean_long <- apply(rho_long, c(1,2), mean)
cat("  Rho:\n")
print(round(rho_mean_long, 3))
cat("  True Rho:\n")
print(rho_true)

# Check off-diagonals
cat("  Rho[1,2]: true=", rho_true[1,2], ", mean=", rho_mean_long[1,2], "\n")

# Track x variance over time
var_x1_trace <- numeric(10)
var_x2_trace <- numeric(10)
checkpoints <- seq(200, 2000, by = 200)
for (i in seq_along(checkpoints)) {
  x_i <- matrix(fit_long$x[, checkpoints[i]], n, d)
  var_x1_trace[i] <- var(x_i[,1])
  var_x2_trace[i] <- var(x_i[,2])
}
cat("  Var(x[,1]) at checkpoints:", round(var_x1_trace, 4), "\n")
cat("  Var(x[,2]) at checkpoints:", round(var_x2_trace, 4), "\n")
cat("  (Checkpoints every 200 iter from 200 to 2000)\n")

# Test 10: Check if the issue is in rinvwishart
cat("\n=== Test 10: rinvwishart check ===\n")
# Expected value of IW(df, S) is S / (df - p - 1)
test_df <- 30
test_S <- matrix(c(25, 5, 5, 25), 2, 2)
iw_samples <- replicate(1000, rinvwishart(test_df, test_S), simplify = FALSE)
iw_mean <- Reduce(`+`, iw_samples) / length(iw_samples)
expected_mean <- test_S / (test_df - 2 - 1)
cat("  IW(30, [[25,5],[5,25]]) samples mean:\n")
print(round(iw_mean, 3))
cat("  Expected mean (S/(df-p-1)):\n")
print(round(expected_mean, 3))

# Test 11: MALA only (fix alpha, kappa, rho at true values)
cat("\n=== Test 11: MALA-only test (fixed alpha, kappa, rho) ===\n")
x_curr <- x_true_mat
alpha_fixed <- alpha_true
rhoxQ_fixed <- Q_rhoxQ(rho_true, kappa_true, G)$rhoxQ

var_trace <- matrix(0, 200, 2)
xQx_trace <- matrix(0, 200, 4)  # [1,1], [2,2], [1,2], [2,1]

for (i in 1:200) {
  # Only update x, keeping alpha fixed
  mala_result <- mala_step(
    x = as.vector(x_curr),
    alpha = alpha_fixed,
    y = y,
    w = rep(1, n),
    Q = rhoxQ_fixed,
    A = A,
    priors = priors,
    d = d,
    logstep0 = priors$logstep.MALA,
    iter = i
  )

  new <- mala_result$sample
  # Take x but keep alpha fixed
  x_curr <- matrix(new[1:(n*d)], n, d)

  var_trace[i,] <- c(var(x_curr[,1]), var(x_curr[,2]))
  xQx <- t(x_curr) %*% as.matrix(Q_marg) %*% x_curr
  xQx_trace[i,] <- c(xQx[1,1], xQx[2,2], xQx[1,2], xQx[2,1])
}

cat("  Var(x) at start:", var(x_true_mat[,1]), var(x_true_mat[,2]), "\n")
cat("  Var(x) at iter 50:", var_trace[50,], "\n")
cat("  Var(x) at iter 200:", var_trace[200,], "\n")
cat("  x'Qx[2,2] trace (every 25):", xQx_trace[seq(1,200,25), 2], "\n")

# Test 12: Check gradient direction
cat("\n=== Test 12: Gradient direction at x=0 ===\n")
x_zero <- rep(0, n*d)
result_zero <- neg_log_lik(x_zero, alpha_true, A, y, priors$alpha.a, priors$alpha.b,
                           rhoxQ_fixed, d, rep(1,n))

# Gradient of log-posterior (not negative)
grad_log_post <- -result_zero$dL[1:(n*d)]

cat("  |grad| at x=0:", sqrt(sum(grad_log_post^2)), "\n")
cat("  Mean grad for component 1:", mean(grad_log_post[1:n]), "\n")
cat("  Mean grad for component 2:", mean(grad_log_post[(n+1):(2*n)]), "\n")
cat("  (Positive = should move away from 0, Negative = should move toward 0)\n")

# Check gradient at true x
result_true <- neg_log_lik(as.vector(x_true_mat), alpha_true, A, y, priors$alpha.a, priors$alpha.b,
                            rhoxQ_fixed, d, rep(1,n))
grad_log_post_true <- -result_true$dL[1:(n*d)]
cat("  |grad| at x_true:", sqrt(sum(grad_log_post_true^2)), "\n")
cat("  Mean grad at x_true for comp 1:", mean(grad_log_post_true[1:n]), "\n")
cat("  Mean grad at x_true for comp 2:", mean(grad_log_post_true[(n+1):(2*n)]), "\n")

# Test 13: Prior-only sampling (verify GMRF sampling works)
cat("\n=== Test 13: Prior-only sampling ===\n")
# Sample from GMRF prior: x ~ N(0, (rho^{-1} ⊗ Q)^{-1})
Q_marg_true <- Q_rhoxQ(NULL, kappa_true, G)$Q
rhoxQ_true_test <- Q_rhoxQ(rho_true, kappa_true, G)$rhoxQ

# Compute expected variance from prior
# Var(x[,k]) at location i = rho[k,k] * (Q^{-1})[i,i]
Q_inv_diag <- diag(solve(as.matrix(Q_marg_true)))
expected_var_1 <- rho_true[1,1] * mean(Q_inv_diag)
expected_var_2 <- rho_true[2,2] * mean(Q_inv_diag)
cat("  Expected var from prior:", expected_var_1, expected_var_2, "\n")

# Sample from prior
n_samples <- 200
prior_var_1 <- numeric(n_samples)
prior_var_2 <- numeric(n_samples)
R_Q_test <- chol(rhoxQ_true_test)
for (i in 1:n_samples) {
  z <- rnorm(n * d)
  x_sample <- backsolve(R_Q_test, z)
  x_mat <- matrix(x_sample, n, d)
  prior_var_1[i] <- var(x_mat[,1])
  prior_var_2[i] <- var(x_mat[,2])
}
cat("  Sampled var from prior: mean =", mean(prior_var_1), mean(prior_var_2),
    ", sd =", sd(prior_var_1), sd(prior_var_2), "\n")

# Compare with true x
cat("  True x variance:", var(x_true_mat[,1]), var(x_true_mat[,2]), "\n")

# Test 14: Check if likelihood dominates prior
cat("\n=== Test 14: Likelihood vs Prior strength ===\n")
# Compute log-likelihood and log-prior at x_true
ll_true <- sum(lgamma(alpha_true) - rowSums(lgamma(alpha_true * z_true)) +
               rowSums((alpha_true * z_true - 1) * log(y)))
lp_true <- -0.5 * as.numeric(t(as.vector(x_true_mat)) %*% rhoxQ_true_test %*% as.vector(x_true_mat))
cat("  Log-likelihood at x_true:", ll_true, "\n")
cat("  Log-prior at x_true:", lp_true, "\n")
cat("  Ratio |prior/likelihood|:", abs(lp_true/ll_true), "\n")

# Compute at x = 0
z_zero <- invalr(matrix(0, n, d))
ll_zero <- sum(lgamma(alpha_true) - rowSums(lgamma(alpha_true * z_zero)) +
               rowSums((alpha_true * z_zero - 1) * log(y)))
lp_zero <- 0  # x'Qx = 0 when x = 0
cat("  Log-likelihood at x=0:", ll_zero, "\n")
cat("  Log-posterior ratio (true vs 0):", (ll_true + lp_true) - (ll_zero + lp_zero), "\n")

# Test 15: Larger grid (10x10)
cat("\n=== Test 15: Larger grid (10x10) ===\n")
set.seed(123)
n_side_large <- 10
n_large <- n_side_large^2

G_large <- create_Q(c(n_side_large, n_side_large), alpha = 1)
Q_marg_large <- Q_rhoxQ(NULL, kappa_true, G_large)$Q
rhoxQ_large <- Q_rhoxQ(rho_true, kappa_true, G_large)$rhoxQ

# Sample x from GMRF prior
R_Q_large <- chol(rhoxQ_large)
z_rand_large <- rnorm(n_large * d)
x_true_large <- backsolve(R_Q_large, z_rand_large)
x_true_mat_large <- matrix(x_true_large, nrow = n_large, ncol = d)

# Transform to compositional space and generate observations
z_true_large <- invalr(x_true_mat_large)
y_large <- matrix(0, n_large, D)
for (i in 1:n_large) {
  y_large[i,] <- rdirichlet(1, alpha_true * z_true_large[i,])
}

# Setup
A_large <- Matrix::Diagonal(n_large * d)
priors_large <- create_priors(G_large, d = 2)

theta0_large <- list(
  alpha = alpha_true,
  kappa = kappa_true,
  rho = rho_true,
  x = x_true_mat_large
)

fit_large <- mcmc_sampling(theta0_large, priors_large, y_large, iter = 1000, verbose = FALSE)

cat("  MALA acceptance:", fit_large$MALA$acc, "\n")
cat("  Gibbs acceptance:", fit_large$Gibbs$acc, "\n")
cat("  Alpha: true=", alpha_true, ", mean=", mean(fit_large$alpha[500:1000]), "\n")
cat("  Kappa: true=", kappa_true, ", mean=", mean(fit_large$kappa[500:1000]), "\n")

rho_large <- array(fit_large$rho[,500:1000], dim=c(2,2,501))
rho_mean_large <- apply(rho_large, c(1,2), mean)
cat("  Rho:\n")
print(round(rho_mean_large, 3))
cat("  True Rho:\n")
print(rho_true)

# x variance at end
x_final_large <- matrix(fit_large$x[,1000], n_large, d)
cat("  True x variance:", var(x_true_mat_large[,1]), var(x_true_mat_large[,2]), "\n")
cat("  Final x variance:", var(x_final_large[,1]), var(x_final_large[,2]), "\n")

# Test 16: Larger alpha (more informative likelihood)
cat("\n=== Test 16: Larger alpha=50 (10x10 grid) ===\n")
alpha_high <- 50
set.seed(456)

# Generate new data with higher alpha
y_high <- matrix(0, n_large, D)
for (i in 1:n_large) {
  y_high[i,] <- rdirichlet(1, alpha_high * z_true_large[i,])
}

theta0_high <- list(
  alpha = alpha_high,
  kappa = kappa_true,
  rho = rho_true,
  x = x_true_mat_large
)

fit_high <- mcmc_sampling(theta0_high, priors_large, y_high, iter = 1000, verbose = FALSE)

cat("  MALA acceptance:", fit_high$MALA$acc, "\n")
cat("  Gibbs acceptance:", fit_high$Gibbs$acc, "\n")
cat("  Alpha: true=", alpha_high, ", mean=", mean(fit_high$alpha[500:1000]), "\n")
cat("  Kappa: true=", kappa_true, ", mean=", mean(fit_high$kappa[500:1000]), "\n")

rho_high <- array(fit_high$rho[,500:1000], dim=c(2,2,501))
rho_mean_high <- apply(rho_high, c(1,2), mean)
cat("  Rho estimated:\n")
print(round(rho_mean_high, 3))
cat("  True Rho:\n")
print(rho_true)

x_final_high <- matrix(fit_high$x[,1000], n_large, d)
cat("  True x variance:", var(x_true_mat_large[,1]), var(x_true_mat_large[,2]), "\n")
cat("  Final x variance:", var(x_final_high[,1]), var(x_final_high[,2]), "\n")

# Test 17: Use smaller true rho (matching what the algorithm estimates)
cat("\n=== Test 17: Smaller true rho (closer to estimated) ===\n")
rho_small <- matrix(c(0.5, 0.05, 0.05, 0.5), 2, 2)
set.seed(789)

# Generate data with smaller rho
rhoxQ_small <- Q_rhoxQ(rho_small, kappa_true, G_large)$rhoxQ
R_Q_small <- chol(rhoxQ_small)
z_rand_small <- rnorm(n_large * d)
x_true_small <- backsolve(R_Q_small, z_rand_small)
x_true_mat_small <- matrix(x_true_small, nrow = n_large, ncol = d)
z_true_small <- invalr(x_true_mat_small)

y_small <- matrix(0, n_large, D)
for (i in 1:n_large) {
  y_small[i,] <- rdirichlet(1, alpha_true * z_true_small[i,])
}

theta0_small <- list(
  alpha = alpha_true,
  kappa = kappa_true,
  rho = rho_small,
  x = x_true_mat_small
)

fit_small <- mcmc_sampling(theta0_small, priors_large, y_small, iter = 1000, verbose = FALSE)

cat("  True rho:\n")
print(rho_small)
cat("  Alpha: true=", alpha_true, ", mean=", mean(fit_small$alpha[500:1000]), "\n")
cat("  Kappa: true=", kappa_true, ", mean=", mean(fit_small$kappa[500:1000]), "\n")

rho_fit_small <- array(fit_small$rho[,500:1000], dim=c(2,2,501))
rho_mean_small <- apply(rho_fit_small, c(1,2), mean)
cat("  Estimated rho:\n")
print(round(rho_mean_small, 3))

cat("  True x variance:", var(x_true_mat_small[,1]), var(x_true_mat_small[,2]), "\n")
x_final_small <- matrix(fit_small$x[,1000], n_large, d)
cat("  Final x variance:", var(x_final_small[,1]), var(x_final_small[,2]), "\n")

# Compute x'Qx for true and final x
Q_marg_small <- Q_rhoxQ(NULL, mean(fit_small$kappa[500:1000]), G_large)$Q
xQx_true_small <- t(x_true_mat_small) %*% as.matrix(Q_marg_small) %*% x_true_mat_small
xQx_final_small <- t(x_final_small) %*% as.matrix(Q_marg_small) %*% x_final_small
cat("  True x'Qx diag:", diag(xQx_true_small), "\n")
cat("  Final x'Qx diag:", diag(xQx_final_small), "\n")

# Expected rho from IW posterior
expected_rho <- (priors_large$rho.Sigma + xQx_final_small) / (priors_large$rho.df + n_large - d - 1)
cat("  Expected rho from x'Qx:\n")
print(round(expected_rho, 3))
