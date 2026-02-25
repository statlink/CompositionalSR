alfareg.nr <- function(y, x, alpha = 1, beta_init = NULL, max_iter = 100,
                       tol = 1e-6, line_search = TRUE) {

  runtime <- proc.time()
  x <- model.matrix(y ~., data.frame(x))
  n <- dim(x)[1]  ;  D <- dim(y)[2]
  d <- D - 1
  p <- dim(x)[2]
  total_params <- d * p

  # Initialize parameters
  if ( is.null(beta_init) ) {
    beta_vec <- rep(0, total_params)
  } else  beta_vec <- beta_init

  # Precompute constants
  H <- Compositional::helm(D)  # d x D
  tH <- t(H)
  ya <- Compositional::alfa(y, alpha)$aff  # n x d
  D_over_alpha <- D / alpha
  # Precompute index lists
  beta_idx <- lapply(1:d, function(k) ((k-1)*p + 1):(k*p))

  # Main Newton-Raphson loop
  for ( iter in 1:max_iter ) {
    # ========== Compute fitted values and objective ==========
    beta_list <- matrix(beta_vec, ncol = d) ## lapply(beta_idx, function(idx)  beta_vec[idx])
    #mu <- inverse_additive_logit(X, beta_list)  # n x D
    mu <- cbind(1, exp(x %*% beta_list) )
    mu / Rfast::rowsums(mu)
    ma <- Compositional::alfa(mu, alpha)$aff    # n x d
    R_alpha <- ya - ma                # n x d
    obj_val <- sum( R_alpha^2 )
    # ========== Precompute shared quantities ==========
    mu_alpha <- mu^alpha              # n x D
    mu_alpha_m1 <- mu^(alpha - 1)     # n x D
    T_sum <- Rfast::rowsums(mu_alpha)        # n x 1
    T_sum_inv <- 1 / T_sum            # n x 1
    T_sum_inv2 <- T_sum_inv^2         # n x 1
    # Jacobian of power transformation (n x D)
    J_diag <- alpha * mu_alpha_m1 * T_sum_inv * (1 - mu_alpha * T_sum_inv)
    # Residuals projected through Helmert
    # R_alpha: (n x d), H: (d x D) => R_H: (n x D)
    R_H <- D_over_alpha * (R_alpha %*% H)
    # ========== Compute ALL J_mu matrices ==========
    J_mu_all <- array(0, dim = c(n, D, d))

    mu_1 <- mu[, 1]
    for (k in 1:d) {
      mu_k <- mu[, k + 1]
      J_mu_all[, 1, k] <-  -mu_1 * mu_k
      J_mu_all[, k + 1, k] <- mu_k * (1 - mu_k)

      if (D > 2) {
        for (pp in 2:D) {
          if (pp != (k + 1)) {
            J_mu_all[, pp, k] <-  -mu[, pp] * mu_k
          }
        }
      }
    }
    # ========== Compute dM_alpha for all components ==========
    dM_all <- array(0, dim = c(n, d, d))
    alpha_T_inv2 <-  -alpha * T_sum_inv2

    for (k in 1:d) {
      # Diagonal contribution (n x D)
      J_u_J_mu <- J_diag * J_mu_all[, , k]
      # Off-diagonal contribution
      sum_mu_alpha_m1_J_mu <- Rfast::rowsums(mu_alpha_m1 * J_mu_all[, , k])

      for (ell in 1:D) {
        off_contrib <- alpha_T_inv2 * mu_alpha[, ell] * sum_mu_alpha_m1_J_mu
        off_contrib <- off_contrib - alpha_T_inv2 * mu_alpha[, ell] * mu_alpha_m1[, ell] * J_mu_all[, ell, k]
        J_u_J_mu[, ell] <- J_u_J_mu[, ell] + off_contrib
      }
      dM_all[, , k] <- D_over_alpha * (J_u_J_mu %*% tH)
    }
    # ========== Compute gradient ==========
    gradient <- numeric(total_params)
    R_H_mu_alpha <- R_H * mu_alpha
    sum_R_H_mu_alpha <- Rfast::rowsums(R_H_mu_alpha)

    for ( k in 1:d ) {
      # Diagonal contribution
      w_k_diag <- Rfast::rowsums(R_H * J_diag * J_mu_all[, , k])
      # Off-diagonal contribution
      mu_alpha_m1_J_mu <- mu_alpha_m1 * J_mu_all[, , k]
      sum_mu_alpha_m1_J_mu <- Rfast::rowsums(mu_alpha_m1_J_mu)
      full_prod <- sum_R_H_mu_alpha * sum_mu_alpha_m1_J_mu
      diag_prod <- Rfast::rowsums(R_H_mu_alpha * mu_alpha_m1 * J_mu_all[, , k])
      w_k_offdiag <- -alpha * T_sum_inv2 * (full_prod - diag_prod)
      w_k <- w_k_diag + w_k_offdiag
      gradient[ beta_idx[[ k ]] ] <-  -crossprod(x, w_k)
    }
    grad_norm <- sqrt( sum(gradient^2) )
    # Check convergence
    if ( grad_norm < tol ) {
      runtime <- proc.time() - runtime
      return( list(runtime = runtime, iters = iter, objective = obj_val, be = beta_list,
                   est = mu, covb = solve(Hess) ) )
    }

    # ========== Compute Hessian ==========
    W_matrix <- matrix(0, nrow = n, ncol = d * d)
    idx <- 1
    for ( k in 1:d ) {
      for ( kp in 1:d ) {
        W_matrix[, idx] <- Rfast::rowsums(dM_all[, , k] * dM_all[, , kp])
        idx <- idx + 1
      }
    }

    Hess <- matrix(0, nrow = total_params, ncol = total_params)
    idx <- 1
    for ( k in 1:d ) {
      for ( kp in 1:d ) {
        x_weighted <- x * W_matrix[, idx]
        Hess_block <- crossprod(x_weighted, x)
        Hess[ beta_idx[[ k ]], beta_idx[[ kp ]] ] <- Hess_block
        idx <- idx + 1
      }
    }
    # Regularization
    Hess <- Hess + diag(1e-8, total_params)
    # ========== Solve Newton system ==========
    direction <- tryCatch({
      solve(Hess, -gradient)
    }, error = function(e) {
      -gradient / max(grad_norm, 1e-8)
    })
    # ========== Line search ==========
    if ( line_search ) {
      lambda <- 1.0
      c1 <- 1e-4
      rho <- 0.5
      max_bt <- 20

      grad_dot_dir <- sum(gradient * direction)

      for ( bt_iter in 1:max_bt ) {
        beta_new <- beta_vec + lambda * direction
        beta_list_new <- matrix(beta_new, ncol = d) ## beta_list_new <- lapply(beta_idx, function(idx) beta_new[idx])
        ##mu_new <- inverse_additive_logit(X, beta_list_new)
        mu_new <- cbind(1, exp(x %*% beta_list_new) )
        mu_new <- mu_new / Rfast::rowsums(mu_new)
        ma_new <- Compositional::alfa(mu_new, alpha)$aff
        obj_new <- sum( (ya - ma_new)^2 )

        if ( obj_new <= obj_val + c1 * lambda * grad_dot_dir ) {
          beta_vec <- beta_new
          break
        }
        lambda <- rho * lambda
      }
    } else  beta_vec <- beta_vec + direction
  }

  runtime <- proc.time() - runtime
  beta_list <- lapply(beta_idx, function(idx) beta_new[idx]) # beta_list <- lapply(beta_idx, function(idx) beta_vec[idx])

  list(runtime = runtime, iters = max_iter, objective = obj_val, be = beta_list, est = mu, covb = solve(Hess) )
}












