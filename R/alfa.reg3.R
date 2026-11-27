alfa.reg3 <- function(y, x, a = c(-1, 1), xnew = NULL) {

  reg <- function(para, ya, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector(ya - ma)
  }

  reg2 <- function(a, y, ha, x, xxinv, d, D) {
    ya <- Compositional::alfa(y, a)$aff
    ax <- a * x
    ini <- as.vector( xxinv %*% crossprod(x, ya) )
    suppressWarnings({
      mod <- minpack.lm::nls.lm(par = ini, fn = reg, ya = ya, ax = ax, a = a, ha = ha, d = d, D = D)
    })
    be <- matrix(mod$par, ncol = d)
    est <- cbind( 1, exp(x %*% be) )
    est <- est/Rfast::rowsums(est)
    sum(y * log(y / est), na.rm = TRUE)
  }

  x <- model.matrix(y ~., data.frame(x) )
  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  ha <- t( Compositional::helm(D) )
  if ( min(y) == 0 )  a <- c(0.0001, 1)
  ha <- t( Compositional::helm(D) )
  xxinv <- solve( crossprod(x) )

  runtime <- proc.time()
  a <- optimize(reg2, a, y = y, ha = ha, x = x, xxinv = xxinv, d = d, D = D)$minimum
  runtime <- proc.time() - runtime
  res <- CompositionalSR::areg(y = y, x = x[, -1], a = a, xnew = xnew)
  res$runtime <- res$runtime + runtime
  res$alfa <- a
  res
}



