alfa.reg2 <- function(y, x, a, xnew = NULL, parallel = FALSE, ncores = 1) {
  
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
  
  if (is.null(colnames(x))) {
    namx <- c("constant", paste("X", 1:p, sep = ""))
  } else {
    namx <- c("constant", colnames(x)[-1])
  }
  
  if (min(y) == 0)  a <- a[a > 0]
  la <- length(a)
  
  # Prepare xnew if provided
  xnew_prepared <- NULL
  if (!is.null(xnew)) {
    xnew <- as.matrix(xnew)
    if (dim(xnew)[1] == 1)  xnew <- t(xnew)
    xnew_prepared <- model.matrix(~., data.frame(xnew))
  }
  
  # Function to fit a single alpha value
  fit_alpha <- function(i) {
    if (a[i] <= 1e-4) {
      mod <- Compositional::comp.reg(y, x[, -1], yb = NULL)
      be <- mod$be
    } else {
      ya <- Compositional::alfa(y, a[i])$aff
      ax <- a[i] * x
      ini <- as.vector( solve(crossprod(x), crossprod(x, ya)) )
      suppressWarnings({
        mod <- minpack.lm::nls.lm(
          par = ini, fn = reg, ya = ya, ax = ax, a = a[i], 
          ha = ha, d = d, D = D,
          control = minpack.lm::nls.lm.control(maxiter = 10000)
        )
      })
      be <- matrix(mod$par, ncol = d)
    }
    
    est <- NULL
    if (!is.null(xnew_prepared)) {
      est <- cbind(1, exp(xnew_prepared %*% be))
      est <- est / Rfast::rowsums(est)
    }
    rownames(be) <- namx
    list(be = be, est = est)
  }
  
  # Parallel or serial execution
  if ( ncores > 1 ) {
    cl <- parallel::makeCluster(ncores)
    # Export necessary objects and functions
    parallel::clusterExport(cl, 
              varlist = c("y", "x", "a", "xnew_prepared", "ha", "d", "D", "namx", "reg", 
              "fit_alpha"), envir = environment()) 
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