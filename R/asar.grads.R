asar.grads <- function(y, x, a, rho, be, coords, k) {

  n <- dim(y)[1]  ;  D <- dim(y)[2]
  d <- D - 1
  x <- model.matrix(y~., data = as.data.frame(x) )
  p <- dim(x)[2]
  total_params <- d * p + 1
  W <- CompositionalSR::contiguity(coords, k)

  B <- t(be)
  I_n <- diag(n)
  S_inv <- solve(I_n - rho * W)
  x_tilde <- S_inv %*% x
  zz <- cbind( 1, exp( x_tilde %*% be) )
  ta <- Rfast::rowsums(zz)
  mu <- round(zz / ta, 9)
  ya <- Compositional::alfa(y, a)$aff
  mua <- Compositional::alfa(mu, a)$aff
  R_a <- ya - mua
  H <- Compositional::helm(D)
  grad_matrix <- matrix(0, nrow = n, ncol = total_params)
  dx_tilde_drho <- S_inv %*% W %*% x_tilde

  # --- Compute gradient for each observation ---
  for ( i in 1:n ) {
    grad_params <- numeric(total_params)
    param_idx <- 1
    # --- Gradient w.r.t. each beta_k ---
    for ( k in 1:d ) {
      for ( s in 1:p ) {
        # Gradient w.r.t. beta_ks
        grad_val <- 0
        for ( j in 1:d ) {
          dmu_a_dbeta <- .compute_dmu_a_dbeta_single(mu[i,], x_tilde[i,], H, a, j, k, s)
          grad_val <- grad_val + R_a[i, j] * dmu_a_dbeta
        }
        grad_params[param_idx] <- grad_val
        param_idx <- param_idx + 1
      }
    }
    # --- Gradient w.r.t. rho ---
    grad_rho <- 0
    for ( j in 1:d ) {
      dmu_a_drho <- .compute_dmu_a_drho_single(mu[i,], x_tilde[i,], dx_tilde_drho[i,], B, H, a, j)
      grad_rho <- grad_rho + R_a[i, j] * dmu_a_drho
    }
    grad_params[param_idx] <- grad_rho
    grad_matrix[i, ] <- grad_params
  }

  grad_matrix <- grad_matrix[, c(param_idx, 1:c(d * p) ) ]

  a2 <- colnames(x)
  if ( is.null(a2) ) {
    p <- dim(x)[2] - 1
    a2 <- c("constant", paste("X", 1:p, sep = "") )
  } else a2  <- c("constant", colnames(x)[-1] )
  a1 <- paste("Y", 2:D, sep = "")
  nam <- as.vector( t( outer(a1, a2, paste, sep = ":") ) )
  colnames(grad_matrix) <- c("rho", nam)
  grad_matrix
}


# Compute derivative of mu_a[j] w.r.t. be for a single observation
.compute_dmu_a_dbeta_single <- function( mu_i, x_tilde_i, H, a, j, k, s ) {

  D <- length(mu_i)
  T_i <- sum(mu_i^a)
  result <- 0

  for (ell in 1:D) {
    for (p in 1:D) {
      # Derivative of power transformation: du_ell / dmu_p
      if (ell == p) {
        du_dmu <- (a * mu_i[p]^(a-1) / T_i) * (1 - mu_i[ell]^a / T_i)
      } else  du_dmu <- -(a * mu_i[ell]^a * mu_i[p]^(a-1)) / T_i^2
      # Derivative of multinomial logit: dmu_p / dbeta_ks
      if (p == 1) {
        # Reference category
        dmu_dbeta <- -mu_i[1] * mu_i[k+1] * x_tilde_i[s]
      } else if (p == k + 1) {
        # Same category as k
        dmu_dbeta <- mu_i[p] * (1 - mu_i[p]) * x_tilde_i[s]
      } else {
        # Other categories
        dmu_dbeta <- -mu_i[p] * mu_i[k+1] * x_tilde_i[s]
      }

      result <- result + H[j, ell] * du_dmu * dmu_dbeta
    }
  }

  (D / a) * result
}


# Compute derivative of mu_a[j] w.r.t. rho for a single observation
.compute_dmu_a_drho_single <- function( mu_i, x_tilde_i, dx_tilde_drho_i, B, H, a, j ) {

  D <- length(mu_i)
  d <- D - 1
  T_i <- sum(mu_i^a)
  result <- 0

  for ( ell in 1:D ) {
    for ( p in 1:D ) {
      # Derivative of power transformation
      if (ell == p) {
        du_dmu <- (a * mu_i[p]^(a-1) / T_i) * (1 - mu_i[ell]^a / T_i)
      } else  du_dmu <- -(a * mu_i[ell]^a * mu_i[p]^(a-1)) / T_i^2
      # Derivative of mu_p w.r.t. rho (through x_tilde)
      dmu_drho <- 0
      if (p == 1) {
        # Reference category
        for (k in 1:d)  dmu_drho <- dmu_drho - mu_i[1] * mu_i[k+1] * sum(B[k,] * dx_tilde_drho_i)
      } else {
        # Other categories
        k <- p - 1
        # Direct effect
        dmu_drho <- mu_i[p] * sum(B[k,] * dx_tilde_drho_i)
        # Indirect effect (through denominator)
        for (m in 1:d)  dmu_drho <- dmu_drho - mu_i[p] * mu_i[m+1] * sum(B[m,] * dx_tilde_drho_i)
      }
      result <- result + H[j, ell] * du_dmu * dmu_drho
    }
  }

  (D / a) * result
}


