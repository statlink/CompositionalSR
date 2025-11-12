asar.grads <- function(y, x, a, rho, be, coords, k) {

  W <- CompositionalSR::contiguity(coords, k)
  x <- model.matrix( y~., data = as.data.frame(x) )
  n <- dim(x)[1]  ;  p <- dim(x)[2]
  D <- dim(y)[2]  ;  d <- D - 1
  # Compute S(rho) = I - rho*W and its inverse
  I_n <- diag(n)
  S_rho <- I_n - rho * W
  S_rho_inv <- solve(S_rho)
  # Compute transformed covariates
  x_tilde <- S_rho_inv %*% x
  H <- Compositional::helm(D)
  # Column order: rho, then component1 (all p covariates), component2, ..., component_d
  grad_matrix <- matrix(0, nrow = n, ncol = 1 + p * d)
  
  mu <- cbind(1, exp(x_tilde %*% be) )
  mu <- mu / Rfast::rowsums(mu)
  mu_alpha <- Compositional::alfa(mu, a)$aff   
  y_alpha <- Compositional::alfa(y, a)$aff
  # Residuals in alpha-space
  r_alpha <- y_alpha - mu_alpha
  
  # Compute gradients for each observation
  for ( i in 1:n ) {
    x_tilde_i <- x_tilde[i, ]
    mu_i <- mu[i, ]
    # Compute Jacobian of power transformation
    T_i <- sum(mu_i^a)
    J_u <- matrix(0, D, D)
    for ( ell in 1:D ) {
      for ( p_idx in 1:D ) {
        J_u[ell, p_idx] <- (a * mu_i[p_idx]^(a-1) / T_i) * 
          ( ifelse(ell == p_idx, 1, 0) - mu_i[ell]^a / T_i )
      }
    }
    # Store gradients in matrix (starting with rho)
    col_idx <- 1  # Start with rho
    # Compute gradient with respect to rho (first column)
    dx_tilde_drho <- S_rho_inv %*% W %*% x_tilde
    
    dmu_drho <- numeric(D)
    for ( comp in 1:D ) {
      if (comp == 1 ) {
        dmu_dcomp_drho <- 0
        for ( j in 1:d ) {
          dmu_dcomp_drho <- dmu_dcomp_drho - mu_i[1] * mu_i[j+1] * sum(be[, j] * dx_tilde_drho[i, ])
        }
        dmu_drho[comp] <- dmu_dcomp_drho
      } else {
        k <- comp - 1
        term1 <- mu_i[comp] * sum(be[, k] * dx_tilde_drho[i, ])
        term2 <-  -mu_i[comp] * sum( sapply(1:d, function(j) {
          mu_i[j+1] * sum( be[, j] * dx_tilde_drho[i, ]) } ) )
        dmu_drho[comp] <- term1 + term2
      }
    }
    
    dmu_alpha_drho <- numeric(d)
    for ( j in 1:d ) {
      sum_val <- 0
      for ( ell in 1:D ) {
        for ( p_idx in 1:D ) {
          sum_val <- sum_val + H[j, ell] * J_u[ell, p_idx] * dmu_drho[p_idx]
        }
      }
      dmu_alpha_drho[j] <- D/a * sum_val
    }
    
    grad_matrix[i, col_idx] <- -2 * sum(r_alpha[i, ] * dmu_alpha_drho)
    col_idx <- col_idx + 1
    # Compute gradients with respect to beta coefficients
    # Order: component1 (all p covariates), component2, ..., component_d
    for ( k in 1:d ) {  # For each component
      for ( s in 1:p ) {  # For each covariate
        # Compute derivative of mu with respect to beta_sk
        dmu_dbeta <- numeric(D)
        for ( comp in 1:D ) {
          if ( comp == 1 ) {
            dmu_dbeta[comp] <-  -mu_i[1] * mu_i[k+1] * x_tilde_i[s]
          } else if ( comp == k+1 ) {
            dmu_dbeta[comp] <- mu_i[comp] * (1 - mu_i[comp]) * x_tilde_i[s]
          } else {
            dmu_dbeta[comp] <- -mu_i[comp] * mu_i[k+1] * x_tilde_i[s]
          }
        }
        # Compute derivative of mu_alpha with respect to beta_sk
        dmu_alpha_dbeta <- numeric(d)
        for ( j in 1:d ) {
          sum_val <- 0
          for ( ell in 1:D ) {
            for ( p_idx in 1:D ) {
              sum_val <- sum_val + H[j, ell] * J_u[ell, p_idx] * dmu_dbeta[p_idx]
            }
          }
          dmu_alpha_dbeta[j] <- D/a * sum_val
        }
        
        grad_matrix[i, col_idx] <- -2 * sum(r_alpha[i, ] * dmu_alpha_dbeta)
        col_idx <- col_idx + 1
      }
    }
  }

  a2 <- colnames(x)
  if ( is.null(a2) ) {
    p <- dim(x)[2] - 1
     a2 <- c("constant", paste("X", 1:p, sep = "") )
  } else a2[1] <- "constant"
  a1 <- paste("Y", 2:D, sep = "")
  colnames(grad_matrix) <- c("rho", as.vector( t( outer(a1, a2, paste, sep = ":") ) ) )
  
  grad_matrix
}
