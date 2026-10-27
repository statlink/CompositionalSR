alfa.betest <- function(y, x, a, index = "all", R = 499) {
  
  ya <- Compositional::alfa(y, a)$aff
  
  if ( identical(class(index), "character" ) ) {
    pvalue <- .alfa.all(ya, x, a, R)
  } else {
    pvalue <- .alfa.index(ya, x, a, index = index, R)
  }
  pvalue 
}


.alfa.all <- function(ya, x, a, R = 499) {

  reg <- function(para, ya, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector(ya - ma)
  }

  ## Analytical Jacobian of reg()
  jac <- function(para, ya, ax, a, ha, d, D) {
    n <- dim(ax)[1]   ;   p <- dim(ax)[2]
    beta_idx  <- lapply(1:d, function(k) ((k - 1L) * p + 1L):(k * p))
    resid_idx <- lapply(1:d, function(j) ((j - 1L) * n + 1L):(j * n))

    be <- matrix(para, ncol = d)
    zz <- cbind(1, exp(ax %*% be))
    za <- zz / rowSums(zz)    # n x D
    zh <- za %*% ha           # n x d

    J <- matrix(0, nrow = n * d, ncol = d * p)
    Da <- D / a

    for (j in 1:d) {
      for (k in 1:d) {
        ## scalar weight per observation
        w <- Da * za[, k + 1L] * (ha[k + 1L, j] - zh[, j])
        ## block: n x p  (row-wise scaling of ax by w)
        J[ resid_idx[[j]], beta_idx[[k]] ] <- -w * ax
      }
    }
    J
  }

  x <- model.matrix(ya ~., data.frame(x) )
  ax <- a * x
  
  n <- dim(y)[1]  ;  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  ha  <- t( Compositional::helm(D) )
  ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
  suppressWarnings({
    mod <- minpack.lm::nls.lm(par = ini, fn = reg, jac = jac,
                              ya = ya, ax = ax, a = a, ha = ha, d = d, D = D)
    })
  m1 <- mod$deviance

  m1b <- numeric(R)  
  for ( i in 1:R ) {
    ind <- rangen::Sample.int(n, n)   
    yb <- ya[ind, ]   
    suppressWarnings({
    mod <- minpack.lm::nls.lm(par = ini, fn = reg, jac = jac,
                              ya = yb, ax = ax, a = a, ha = ha, d = d, D = D)
    })
    m1b[i] <- mod$deviance
  }
  ( sum(m1b < m1) + 1) / (R + 1 )
}


.alfa.index <- function(ya, x, a, index, R = 499) {
 
  reg <- function(para, ya, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector(ya - ma)
  }

  ## Analytical Jacobian of reg()
  jac <- function(para, ya, ax, a, ha, d, D) {
    n <- dim(ax)[1]   ;   p <- dim(ax)[2]
    beta_idx  <- lapply(1:d, function(k) ((k - 1L) * p + 1L):(k * p))
    resid_idx <- lapply(1:d, function(j) ((j - 1L) * n + 1L):(j * n))

    be <- matrix(para, ncol = d)
    zz <- cbind(1, exp(ax %*% be))
    za <- zz / rowSums(zz)    # n x D
    zh <- za %*% ha           # n x d

    J <- matrix(0, nrow = n * d, ncol = d * p)
    Da <- D / a

    for (j in 1:d) {
      for (k in 1:d) {
        ## scalar weight per observation
        w <- Da * za[, k + 1L] * (ha[k + 1L, j] - zh[, j])
        ## block: n x p  (row-wise scaling of ax by w)
        J[ resid_idx[[j]], beta_idx[[k]] ] <- -w * ax
      }
    }
    J
  }

  x <- model.matrix(ya ~., data.frame(x) )
  ax <- a * x
  
  n <- dim(y)[1]  ;  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  ha  <- t( Compositional::helm(D) )
  ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
  suppressWarnings({
    mod <- minpack.lm::nls.lm(par = ini, fn = reg, jac = jac,
                              ya = ya, ax = ax, a = a, ha = ha, d = d, D = D)
    })
  m1 <- mod$deviance
  
  index <- index + 1
  X <- x[, -index]
  aX <- ax[, -index]
  ini <- as.vector( solve(crossprod(X), crossprod(X, ya) ) )
  suppressWarnings({
    mod <- minpack.lm::nls.lm(par = ini, fn = reg, jac = jac,
                              ya = ya, ax = aX, a = a, ha = ha, d = d, D = D)
    })
  m0 <- mod$deviance
  
  tobs <- m0 - m1 
  
  m1b <- numeric(R)
  Z <- x
  aZ <- ax
  for ( i in 1:R ) {
    ind <- rangen::Sample.int(n, n)   
    Z[, index] <- x[ind, index]   
    aZ[, index] <- ax[ind, index] 
    ini <- as.vector( solve(crossprod(Z), crossprod(Z, ya) ) )
    suppressWarnings({
    mod <- minpack.lm::nls.lm(par = ini, fn = reg, jac = jac,
                              ya = ya, ax = aZ, a = a, ha = ha, d = d, D = D)
    })
    m1b[i] <- mod$deviance
  }
  
  ( sum( m0 - m1b > tobs) + 1 ) / (R + 1) 
}
  

 