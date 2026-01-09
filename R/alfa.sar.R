alfa.sar <- function(y, x, a, coords, k = 10, covb = FALSE, xnew = NULL, coordsnew, yb = NULL) {

  rega <- function(para, ya, W, ax, a, In, ha, d, D) {
    rho <- para[1]
    be <- matrix(para[-1], ncol = d)
    if ( abs(rho) > 1e-4 ) {
      sw <- In - rho * W
      zz <- cbind( 1, exp( solve(sw, ax %*% be) ) )
    } else  zz <- cbind( 1, exp( ax %*% be ) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector(ya - ma)
  }

  regarho <- function(para, sw, ya, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp( sw %*% (ax %*% be) ) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector(ya - ma)
  }

  runtime <- proc.time()

  W <- CompositionalSR::contiguity(coords, k)
  if ( is.null(yb) ) {
    ya <- Compositional::alfa(y, a)$aff
  } else  ya <- yb
  x <- model.matrix( ya ~., data.frame(x) )

  ax <- a * x

  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  n <- dim(x)[1]
  In <- diag(n)

  if ( a <= 1e-5 ) {
    mod <- Compositional::comp.reg(y, x[, -1], yb = yb)
    be <- mod$be

  } else {
    ha <- t( Compositional::helm(D) )
    p <- dim(x)[2]
    ini <- as.vector( CompositionalSR::alfa.reg(y, x[, -1], a, yb = ya)$be )
    r <- seq(-0.9, 0.9, by = 0.2)
    bes <- matrix(nrow = length(r), ncol = d * p)
    dev <- numeric( length(r) )
    suppressWarnings({
      for ( j in 1:length(r) ) {
        ## sw <- spatialreg::invIrW(x = W, rho = r[j], method = "solve", feasible = NULL)
        sw <- solve( In - r[j] * W )
        mod1 <- minpack.lm::nls.lm( par = ini, fn = regarho, ya = ya, sw = sw, ax = ax, a = a, ha = ha, d = d, D = D,
                                    control = minpack.lm::nls.lm.control(maxiter = 10000)  )
        bes[j, ] <- mod1$par
        dev[j] <- mod1$deviance
      }
    })
    ind <- which.min(dev)
    suppressWarnings({
      if ( r[ind] > 0 ) {
        bou <- c( r[ind] - 0.2, min(1, r[ind] + 0.2) )
      } else bou <- c( max(-1, r[ind] - 0.2), r[ind] + 0.2 )
      mod <- minpack.lm::nls.lm( par = c(r[ind], bes[ind, ] ), fn = rega, ya = ya, W = W, ax = ax, a = a, In = In, ha = ha, d = d, D = D,
                                 control = minpack.lm::nls.lm.control(maxiter = 10000), lower = c(bou[1], rep(-Inf, p * d) ),
                                 upper = c(bou[2], rep(Inf, p * d) ) )
    })
    rho <- mod$par[1]
    be <- matrix(mod$par[-1], ncol = d)
  }  ## end if (a == 0)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    Xaug <- rbind(x, xnew)
    Inm <- diag( dim(Xaug)[1] )
    coordsaug <- rbind(coords, coordsnew)
    Waug <- contiguity(coordsaug, k = k)
    ## sw <- spatialreg::invIrW(x = Waug, rho = rho, method = "solve", feasible = NULL)
    zz <- cbind( 1, exp( solve(Inm - rho * Waug, Xaug %*% be) ) )
    ta <- rowSums(zz)
    est <- zz / ta
    est <- round(est[-c(1:n), ], 9)
  } else {
    ## sw <- spatialreg::invIrW(x = W, rho = rho, method = "solve", feasible = NULL)
    zz <- cbind( 1, exp( solve(In - rho * W, x %*% be) ) )
    ta <- rowSums(zz)
    est <- zz / ta
    est <- round(est, 9)
  }

  runtime <- proc.time() - runtime

  p <- dim(x)[2] - 1
  if ( is.null( colnames(x) ) ) {
    rownames(be) <- c("constant", paste("X", 1:p, sep = "") )
  } else rownames(be)  <- c("constant", colnames(x)[-1] )
  colnames(be) <- paste("Y", 2:D, sep = "")

  if ( covb ) {
    res <- optim( c(rho, as.vector(be)), .regsar, ya = ya, W = W, ax = ax, a = a,
                  In = In, ha = ha, d = d, D = D, hessian = TRUE, method = "BFGS",
                  control = list(maxit = 1000) )
    covbe <- solve(res$hessian)
    diag(covbe) <- abs( diag(covbe) )
    a2 <- colnames(x)
    if ( is.null(a2) ) {
      p <- dim(x)[2] - 1
      a2 <- c("constant", paste("X", 1:p, sep = "") )
    } else a2[1] <- "constant"
    a1 <- paste("Y", 2:D, sep = "")
    nam <- as.vector( t( outer(a1, a2, paste, sep = ":") ) )
    colnames(covbe) <- rownames(covbe) <- c("rho", nam)
  } else covbe <- NULL

  list(runtime = runtime, rho = rho, be = be, covbe = covbe, dev = mod$deviance, est = est)
}


.regsar <- function(para, ya, W, ax, a, In, ha, d, D) {
  rho <- para[1]
  be <- matrix(para[-1], ncol = d)
  if ( abs(rho) > 1e-4 ) {
    sw <- In - rho * W
    zz <- cbind( 1, exp( solve(sw, ax %*% be) ) )
  } else  zz <- cbind( 1, exp( ax %*% be ) )
  ta <- rowSums(zz)
  za <- zz / ta
  ma <- ( D / a * za - 1/a ) %*% ha
  sum( (ya - ma)^2 )
}

