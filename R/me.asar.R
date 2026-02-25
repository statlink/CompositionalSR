################################################################################
# Marginal Effects for a-SAR Model
################################################################################
me.asar <- function(be, rho, mu, x, coords, k, cov_theta= NULL) {
  
  # Extract model components
  beta <- t( be )[, -1, drop = FALSE]        # d x p matrix

  n <- dim(mu)[1]
  D <- dim(mu)[2]
  d <- D - 1
  x <- as.matrix(x)
  p <- dim(x)[2]
  W <- CompositionalSR::contiguity(coords, k)
  
  # Compute S(rho)^{-1}
  In <- diag(n)
  S_inv <- solve(In - rho * W)
  # Compute transformed covariates
  x_tilde <- S_inv %*% x
  
  # Initialize arrays for marginal effects
  # ME[i, ell, k] = marginal effect at location i, component ell, covariate k
  ME_direct <- array(0, dim = c(n, D, p))
  ME_indirect <- array(0, dim = c(n, D, p))
  ME_total <- array(0, dim = c(n, D, p))
  
  # Compute all three types of marginal effects
  ME_direct <- .compute_direct_effects(mu, beta, S_inv, d, D, n, p)
  ME_indirect <- .compute_indirect_effects(mu, beta, S_inv, d, D, n, p)
  ME_total <- .compute_total_effects(mu, beta, S_inv, d, D, n, p)
  
  # Compute average marginal effects
  AME_direct <- t( apply(ME_direct, c(2, 3), mean) ) # D x p matrix
  AME_indirect <- t( apply(ME_indirect, c(2, 3), mean) )
  AME_total <- t( apply(ME_total, c(2, 3), mean) )
  
  # Set row and column names
  comp_names <- paste0("Y", 1:D)
  cov_names <- colnames(x)
  
  colnames(AME_direct) <- colnames(AME_indirect) <-  colnames(AME_total) <- comp_names
  rownames(AME_direct) <- rownames(AME_indirect) <- rownames(AME_total) <- cov_names
  
  # Prepare result list
  result <- list( me.dir = ME_direct, me.indir = ME_indirect, me.total = ME_total,
                  ame.dir = AME_direct, ame.indir = AME_indirect, ame.total = AME_total )
  
  # Compute standard errors for AMEs if covariance matrix provided
  if ( !is.null(cov_theta) ) {
    # Precompute S(rho)^{-1} * W * S(rho)^{-1} for efficiency
    S_inv_W_S_inv <- S_inv %*% W %*% S_inv
    # Compute standard errors for each type
    se_AME_direct <- .compute_ame_se( mu = mu, beta = beta,  rho = rho, S_inv = S_inv, 
      S_inv_W_S_inv = S_inv_W_S_inv, W = W, x = x, x_tilde = x_tilde,  cov_theta = cov_theta, 
      type = "direct", d = d, D = D, n = n, p = p )
    
    se_AME_indirect <- .compute_ame_se( mu = mu, beta = beta, rho = rho, S_inv = S_inv, 
      S_inv_W_S_inv = S_inv_W_S_inv, W = W, x = x, x_tilde = x_tilde, cov_theta= cov_theta, 
      type = "indirect", d = d, D = D, n = n, p = p )
    
    se_AME_total <- .compute_ame_se( mu = mu, beta = beta, rho = rho, S_inv = S_inv, 
      S_inv_W_S_inv = S_inv_W_S_inv, W = W, x = x, x_tilde = x_tilde, cov_theta = cov_theta, 
      type = "total", d = d, D = D, n = n, p = p )
    
    # Add to result
    result$se.amedir <- t( se_AME_direct )
    result$se.ameindir <- t( se_AME_indirect )
    result$se.ametotal <- t( se_AME_total )
  }
 
  result
}


################################################################################
# Helper Functions for Computing Marginal Effects
################################################################################

## Compute Direct Effects
.compute_direct_effects <- function(mu, beta, S_inv, d, D, n, p) {
  ME <- array(0, dim = c(n, D, p))
  # Extract diagonal of S(rho)^{-1}
  S_inv_diag <- diag(S_inv)
  
  for (i in 1:n) {
    for (k in 1:p) {
      # Compute partial derivative with respect to transformed covariate
      dmu_dtilde <- .compute_dmu_dtilde(mu[i, ], beta, d, D, k)
      # Direct effect = dmu/dtilde * [S(rho)^{-1}]_{ii}
      ME[i, , k] <- dmu_dtilde * S_inv_diag[i]
    }
  }
  
  ME
}

## Compute Indirect Effects
.compute_indirect_effects <- function(mu, beta, S_inv, d, D, n, p) {
  ME <- array(0, dim = c(n, D, p))
  
  for (i in 1:n) {
    # Sum of off-diagonal elements in row i
    S_inv_offdiag_sum <- sum(S_inv[i, ]) - S_inv[i, i]
    for (k in 1:p) {
      # Compute partial derivative with respect to transformed covariate
      dmu_dtilde <- .compute_dmu_dtilde(mu[i, ], beta, d, D, k)    
      # Indirect effect = dmu/dtilde * sum_{j != i} [S(rho)^{-1}]_{ij}
      ME[i, , k] <- dmu_dtilde * S_inv_offdiag_sum
    }
  }
  
  ME
}

## Compute Total Effects
.compute_total_effects <- function(mu, beta, S_inv, d, D, n, p) {
  ME <- array(0, dim = c(n, D, p))
  # Row sums of S(rho)^{-1}
  S_inv_rowsum <- rowSums(S_inv)
  
  for (i in 1:n) {
    for (k in 1:p) {
      # Compute partial derivative with respect to transformed covariate
      dmu_dtilde <- .compute_dmu_dtilde(mu[i, ], beta, d, D, k)
      # Total effect = dmu/dtilde * sum_j [S(rho)^{-1}]_{ij}
      ME[i, , k] <- dmu_dtilde * S_inv_rowsum[i]
    }
  }
  
  ME
}

## Compute derivative of mu with respect to transformed covariate
## This is the same formula as in the standard a-regression
.compute_dmu_dtilde <- function( mu_i, beta, d, D, k ) {
  dmu <- numeric(D)
  # For component 1 (reference)
  dmu[1] <- -mu_i[1] * sum(beta[, k] * mu_i[-1])
  # For components 2, ..., D
  for (ell in 2:D) {
    dmu[ell] <- mu_i[ell] * (beta[ell - 1, k] - sum(beta[, k] * mu_i[-1]))
  }
  
  dmu
}

################################################################################
# Standard Errors for AME via Delta Method
################################################################################

## Compute Standard Errors for Average Marginal Effects Only
.compute_ame_se <- function( mu, beta, rho, S_inv, S_inv_W_S_inv, W, x, x_tilde,
                             cov_theta, type, d, D, n, p ) {
  
  # Parameter vector: theta = (vec(beta), rho)
  n_params <- d * p + 1
  # Check covariance matrix dimension
  if (nrow(cov_theta) != n_params || ncol(cov_theta) != n_params) {
    stop( paste("Covariance matrix must be", n_params, "x", n_params) )
  }
  # Initialize array for standard errors of AME
  se_AME <- matrix(0, nrow = D, ncol = p)
  
  # Compute standard errors for each component and covariate
  for (ell in 1:D) {
    for (k in 1:p) {
      # Compute average Jacobian across all observations
      J_avg <- .compute_jacobian_ame( ell = ell, k = k, mu = mu, beta = beta, rho = rho, 
         S_inv = S_inv, S_inv_W_S_inv = S_inv_W_S_inv, W = W, x = x, x_tilde = x_tilde, 
         type = type, d = d, D = D, n = n, p = p )
      
      # Variance via delta method: J_avg * Cov * J_avg^T
      variance <- as.numeric( J_avg %*% cov_theta%*% t(J_avg) ) 
      se_AME[ell, k] <- sqrt( max(0, variance) )
    }
  }
  
  rownames(se_AME) <- paste0("Component_", 1:D)
  colnames(se_AME) <- colnames(x)
  se_AME
}

## Compute Average Jacobian for AME
## Jacobian of AME_{ell,k} = (1/n) * sum_i ME_{i,ell,k}
.compute_jacobian_ame <- function( ell, k, mu, beta, rho, S_inv, S_inv_W_S_inv,
                                    W, x, x_tilde, type, d, D, n, p ) {
  
  n_params <- d * p + 1
  J_sum <- matrix(0, nrow = 1, ncol = n_params)
  # Sum Jacobians across all observations
  for (i in 1:n) {
    J_i <- .compute_jacobian_me_single( i = i, ell = ell, k = k, mu = mu, beta = beta,
      rho = rho, S_inv = S_inv, S_inv_W_S_inv = S_inv_W_S_inv, W = W, x = x,
      x_tilde = x_tilde, type = type, d = d, D = D, n = n, p = p )
    J_sum <- J_sum + J_i
  }
  
  # Average Jacobian
  J_sum / n
}

## Compute Jacobian for a single observation's marginal effect
.compute_jacobian_me_single <- function( i, ell, k, mu, beta, rho, S_inv, 
                                         S_inv_W_S_inv, W, x, x_tilde,
                                         type, d, D, n, p ) {
  
  n_params <- d * p + 1
  J <- matrix(0, nrow = 1, ncol = n_params)
  
  # Get the multiplier based on type
  if (type == "direct") {
    multiplier <- S_inv[i, i]
    d_multiplier_drho <- S_inv_W_S_inv[i, i]
  } else if (type == "indirect") {
    multiplier <- sum(S_inv[i, ]) - S_inv[i, i]
    d_multiplier_drho <- sum(S_inv_W_S_inv[i, ]) - S_inv_W_S_inv[i, i]
  } else { # total
    multiplier <- sum(S_inv[i, ])
    d_multiplier_drho <- sum(S_inv_W_S_inv[i, ])
  }
  
  # Derivatives with respect to beta_{rs}
  for (r in 1:d) {
    for (s in 1:p) {
      param_idx <- (r - 1) * p + s
      # Compute d^2(mu_{iell}) / d(tilde_x_k) / d(beta_{rs})
      d2mu_dtilde_dbeta <- .compute_d2mu_dtilde_dbeta_single(
        mu_i = mu[i, ], beta = beta, k = k, r = r, s = s,
        x_tilde_is = x_tilde[i, s], d = d, D = D )
      
      J[1, param_idx] <- d2mu_dtilde_dbeta[ell] * multiplier
    }
  }
  
  # Derivative with respect to rho
  J[1, n_params] <- .compute_dME_drho_single(
     i = i, ell = ell, k = k, mu = mu, beta = beta, S_inv = S_inv, W = W, 
     x_tilde = x_tilde, multiplier = multiplier, d_multiplier_drho = d_multiplier_drho,
     d = d, D = D, n = n, p = p )
  
  J
}

## Compute second derivative d^2(mu_{iell}) / d(tilde_x_k) / d(beta_{rs})
.compute_d2mu_dtilde_dbeta_single <- function( mu_i, beta, k, r, s, x_tilde_is, d, D ) {
  
  # First compute d(mu_{ip}) / d(beta_{rs}) for all p
  dmu_dbeta <- numeric(D)
  # For reference component (p = 1)
  dmu_dbeta[1] <- -mu_i[1] * mu_i[r + 1] * x_tilde_is
  # For component p = r + 1
  dmu_dbeta[r + 1] <- mu_i[r + 1] * (1 - mu_i[r + 1]) * x_tilde_is
  # For other components
  for (p in 2:D) {
    if (p != r + 1) {
      dmu_dbeta[p] <- -mu_i[p] * mu_i[r + 1] * x_tilde_is
    }
  }
  
  # Now compute second derivative
  d2mu <- numeric(D)
  delta_rk <- as.numeric(r == k)
  delta_sk <- as.numeric(s == k)
  # For component ell = 1 (reference)
  term1 <- -dmu_dbeta[1] * sum(beta[, k] * mu_i[-1])
  term2 <- -mu_i[1] * sum(beta[, k] * dmu_dbeta[-1])
  term3 <- -mu_i[1] * delta_rk * mu_i[r + 1]
  d2mu[1] <- term1 + term2 + term3
  # For components ell = 2, ..., D
  for (ell in 2:D) {
    delta_ell_r <- as.numeric((ell - 1) == r)
    term1 <- dmu_dbeta[ell] * (beta[ell - 1, k] - sum(beta[, k] * mu_i[-1]))
    term2 <- mu_i[ell] * delta_ell_r * delta_sk
    term3 <- -mu_i[ell] * delta_rk * mu_i[r + 1]
    term4 <- -mu_i[ell] * sum(beta[, k] * dmu_dbeta[-1])
    d2mu[ell] <- term1 + term2 + term3 + term4
  }
  
  d2mu
}

## Compute derivative of marginal effect with respect to rho for single observation
.compute_dME_drho_single <- function( i, ell, k, mu, beta, S_inv, W, x_tilde,
                                      multiplier, d_multiplier_drho, d, D, n, p ) {
  
  # Compute d(mu_{iell}) / d(tilde_x_k)
  dmu_dtilde <- .compute_dmu_dtilde(mu[i, ], beta, d, D, k)
  # Compute d/d(rho) [d(mu_{iell}) / d(tilde_x_k)]
  d2mu_dtilde_drho <- .compute_d2mu_dtilde_drho_single(
    i = i, mu_i = mu[i, ], beta = beta, k = k, S_inv = S_inv, 
    W = W, x_tilde = x_tilde, d = d, D = D, n = n, p = p )
  # Combine terms
  term1 <- d2mu_dtilde_drho[ell] * multiplier
  term2 <- dmu_dtilde[ell] * d_multiplier_drho
  term1 + term2
}

## Compute second derivative d^2(mu_{iell}) / d(tilde_x_k) / d(rho)
.compute_d2mu_dtilde_drho_single <- function( i, mu_i, beta, k, S_inv, W, x_tilde, 
                                              d, D, n, p ) {
  
  # Compute d(tilde_x_i) / d(rho) = S(rho)^{-1} * W * tilde_x_i
  W_tilde_x_i <- as.vector((S_inv %*% W %*% x_tilde)[i, ])
  # Compute d(mu_{ip}) / d(rho) for all components p
  dmu_drho <- numeric(D)
  
  for (comp in 1:D) {
    if (comp == 1) {
      # Reference component
      sum_term <- 0
      for (j in 1:d)  sum_term <- sum_term + sum(beta[j, ] * W_tilde_x_i) * mu_i[j + 1]
      dmu_drho[1] <- -mu_i[1] * sum_term
    } else {
      # Other components
      direct_term <- sum(beta[comp - 1, ] * W_tilde_x_i)
      sum_term <- 0
      for (j in 1:d)  sum_term <- sum_term + sum(beta[j, ] * W_tilde_x_i) * mu_i[j + 1]
      dmu_drho[comp] <- mu_i[comp] * (direct_term - sum_term)
    }
  }
  
  # Now compute d^2(mu_{iell}) / d(tilde_x_k) / d(rho)
  d2mu <- numeric(D)
  # For component ell = 1 (reference)
  term1 <- -dmu_drho[1] * sum(beta[, k] * mu_i[-1])
  term2 <- -mu_i[1] * sum(beta[, k] * dmu_drho[-1])
  d2mu[1] <- term1 + term2
  
  # For components ell = 2, ..., D
  for (ell in 2:D) {
    term1 <- dmu_drho[ell] * (beta[ell - 1, k] - sum(beta[, k] * mu_i[-1]))
    term2 <- -mu_i[ell] * sum(beta[, k] * dmu_drho[-1])
    d2mu[ell] <- term1 + term2
  }
  
  d2mu
}
