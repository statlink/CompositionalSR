alfa.reg2 <- function(y, x, a, xnew = NULL, ncores = 1) {

  reg <- function(para, ya, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind(1, exp(ax %*% be))
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- (D / a * za - 1/a) %*% ha
    as.vector(ya - ma)
  }

  runtime <- proc.time()
  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  p <- dim(x)[2] - 1
  x <- model.matrix(y ~., data.frame(x))
  ha <- t( Compositional::helm(D) )

  if ( min(y ) == 0)  a <- a[a > 0]
  la <- length(a)
  xxinv <- solve(crossprod(x) )

  if ( is.null( colnames(x) ) ) {
    p <- dim(x)[2] - 1
    namx <- c("constant", paste("X", 1:p, sep = "") )
  } else  namx <- c("constant", colnames(x)[-1] )
  if ( is.null( colnames(y) ) ) {
    namy <- paste("Y", 2:D, sep = "")
  } else  namy <- colnames(y)[-1]

  # Function to fit a single alpha value
  fit_alpha <- function(i) {
    if ( abs(a[i]) <= 1e-5 ) {
      mod <- Compositional::comp.reg(y, x[, -1], yb = NULL)
      be <- mod$be
    } else {
      ya <- Compositional::alfa(y, a[i])$aff
      ax <- a[i] * x
      ini <- as.vector( xxinv %*% crossprod(x, ya) )
      suppressWarnings({
        mod <- minpack.lm::nls.lm( par = ini, fn = reg, ya = ya, ax = ax, a = a[i], ha = ha,
                                   d = d, D = D, control = minpack.lm::nls.lm.control(maxiter = 1000)  )
      })
      be <- matrix(mod$par, ncol = d)
    }

    est <- NULL
    if ( !is.null(xnew) ) {
      if ( dim(xnew)[1] == 1 )  xnew <- t(xnew)
      xnew <- model.matrix(~., data.frame(xnew) )
      est <- cbind( 1, exp(xnew %*% be) )
      est <- est / Rfast::rowsums(est)
    }
    rownames(be) <- namx
    colnames(be) <- namy
    list(be = be, est = est)
  }

  # Parallel or serial execution
  if ( ncores > 1 ) {
    cl <- parallel::makeCluster(ncores)
    # Export necessary objects and functions
    parallel::clusterExport( cl, varlist = c( "y", "x", "xxinv", "a", "xnew_prepared", "ha", "d",
                             "D", "namx", "namy", "reg", "fit_alpha" ), envir = environment() )
    # Load required packages on each worker
    parallel::clusterEvalQ(cl, {
      library(Compositional)
      library(minpack.lm)
      library(Rfast)
    })
    # Run parallel computation
    res <- parallel::parSapply(cl, 1:la, fit_alpha, simplify = FALSE, USE.NAMES = FALSE)
    # Stop cluster
    parallel::stopCluster(cl)

  } else  res <- sapply(1:la, fit_alpha, simplify = FALSE, USE.NAMES = FALSE)

  runtime <- proc.time() - runtime
  res$runtime <- runtime

  res
}
