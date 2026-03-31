ar.grads <- function(y, x, a, be) {

  n <- dim(y)[1]  ;  D <- dim(y)[2]
  d <- D - 1
  x <- model.matrix( y~., data = as.data.frame(x) )
  p <- ncol(x)

  mu <- cbind( 1, exp(x %*% be) )
  mu <- mu/Rfast::rowsums(mu)

  H <- Compositional::helm(D)
  # Transform observations
  y_a <- Compositional::alfa(y, a)$aff
  mu_a <- Compositional::alfa(mu, a)$aff
  r_a <- y_a - mu_a

  # Initialize storage
  # Each row will be vec(gradient_i), where gradient_i is p x d
  gradient <- matrix(0, n, p * d)
  weights_per_component <- matrix(0, d, n)
  # Compute gradient for each observation
  for (i in 1:n) {
    # Pre-compute Jacobian of power transformation for observation i
    Ju_i <- .compute_Ju(mu[i, ], a)
    # Store gradient as p x d matrix temporarily
    grad_i <- matrix(0, p, d)
    # For each component k (corresponding to b_k)
    for ( k in 1:d ) {
      # Compute Jacobian of multinomial logit
      Jmu_i_k <- .compute_Jmu(mu[i, ], k)
      # Compute weight for this observation and component
      # Formula from Section B.2.2
      weight_ik <- 0
      for ( m in 1:d ) {
        for ( ell in 1:D ) {
          for ( pp in 1:D ) {
            # Combine all terms: r_a x H x Ju x Jmu
            weight_ik <- weight_ik + r_a[i, m] * (D / a) * H[m, ell] * Ju_i[ell, pp] * Jmu_i_k[pp]
          }
        }
      }
      # Store weight
      weights_per_component[k, i] <- weight_ik
      # Gradient for component k: weight x covariate vector
      # grad_i[, k] is a p-vector (column k)
      grad_i[, k] <- weight_ik * x[i, ]
    }
    gradient[i, ] <- as.vector(grad_i)
  }

  a2 <- colnames(x)
  if ( is.null(a2) ) {
    p <- dim(x)[2] - 1
    a2 <- c("constant", paste("X", 1:p, sep = "") )
  } else a2[1] <- "constant"
  a1 <- paste("Y", 2:D, sep = "")
  colnames(gradient) <- as.vector( t( outer(a1, a2, paste, sep = ":") ) )

  gradient
}

.compute_Ju <- function(mu_i, a) {
  D <- length(mu_i)
  if (abs(a) < 1e-9) {
    # For a -> 0, use log-ratio derivative
    Ju <- matrix(0, D, D)
    for ( ell in 1:D ) {
      for ( p in 1:D ) {
        if ( ell == p ) {
          Ju[ell, p] <- 1 / mu_i[p] - 1 / D
        } else  Ju[ell, p] <-  -1 / D
      }
    }
    return(Ju)
  }
  # Power transformation Jacobian
  mu_a <- mu_i^a
  Ti <- sum(mu_a)
  Ju <- matrix(0, D, D)
  for ( ell in 1:D ) {
    for ( p in 1:D ) {
      if ( ell == p ) {
        # Diagonal: Section B.2.3
        Ju[ell, p] <- (a * mu_i[p]^(a - 1) / Ti) * (1 - mu_a[ell] / Ti)
      } else {
        # Off-diagonal
        Ju[ell, p] <- -(a * mu_a[ell] * mu_i[p]^(a - 1)) / Ti^2
      }
    }
  }
  return(Ju)
}

# Compute Jacobian of multinomial logit: 
# Returns D-dimensional vector (one value for each component p)
# k is 1-indexed (k=1 corresponds to b_1)
.compute_Jmu <- function(mu_i, k) {
  D <- length( mu_i )
  mu_k <- mu_i[k + 1]  # The (k+1)-th component
  Jmu <- numeric(D)
  for ( p in 1:D ) {
    if ( p == 1 ) {
    # Reference category
      Jmu[p] <-  -mu_i[1] * mu_k
    } else if (p == k + 1) {
      # Component k+1
      Jmu[p] <- mu_k * (1 - mu_k)
    } else {
      # Other components
      Jmu[p] <- -mu_i[p] * mu_k
    }
  }
  Jmu
}

