################################
#### alfa-regression
#### Tsagris Michail 11/2015
#### mtsagris@yahoo.gr
#### References: Tsagris Michail (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean Journal of Statistics, 6(2): 47-57
################################
areg <- function(y, x, a, covb = FALSE, xnew = NULL, yb = NULL) {

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

  runtime <- proc.time()

  if ( is.null(yb) ) {
    ya <- Compositional::alfa(y, a)$aff
  } else  ya <- yb
  x <- model.matrix(ya ~., data.frame(x) )
  ax <- a * x

  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex

  if ( a <= 1e-5 ) {
    mod <- Compositional::comp.reg(y, x[, -1], yb = yb)
    be <- mod$be

  } else {
    ha  <- t( Compositional::helm(D) )
    ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
    suppressWarnings({
      mod <- minpack.lm::nls.lm(par = ini, fn = reg, jac = jac,
                                ya = ya, ax = ax, a = a, ha = ha, d = d, D = D)
    })
    be <- matrix(mod$par, ncol = d)

  }  ## end if (a == 0)

  runtime <- proc.time() - runtime

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- as.matrix(xnew)
    if ( dim(xnew)[1] == 1 )  xnew <- t(xnew)
    xnew <- model.matrix(~., data.frame(xnew) )
    est <- cbind( 1, exp(xnew %*% be) )
    est <- est/Rfast::rowsums(est)
  }

  if ( is.null( colnames(x) ) ) {
    p <- dim(x)[2] - 1
    rownames(be) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(be) <- c("constant", colnames(x)[-1] )
  }
  colnames(be) <- paste("Y", 2:D, sep = "")

  if ( covb ) {
    res <- optim( as.vector(be), .regar, ya = ya, ax = ax, a = a, ha = ha, d = d,
                  D = D, hessian = TRUE, control = list(maxit = 1000) )
    covbe <- solve(res$hessian)
    a2 <- colnames(x)
    a1 <- colnames(be)
    nam <- as.vector( t( outer(a1, a2, paste, sep = ":") ) )
    colnames(covbe) <- rownames(covbe) <- nam
  } else covbe <- NULL

  list(runtime = runtime, be = be, covbe = covbe, dev = mod$deviance, est = est)
}


.regar <- function(para, ya, ax, a, ha, d, D) {
  be <- matrix(para, ncol = d)
  zz <- cbind( 1, exp(ax %*% be) )
  ta <- rowSums(zz)
  za <- zz / ta
  ma <- ( D / a * za - 1/a ) %*% ha
  sum( (ya - ma)^2 )
}
