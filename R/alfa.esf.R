alfa.esf <- function(y, x, a, coords, model = "exp", covb = FALSE, xnew = NULL, coordsnew, yb = NULL) {

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

  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  if ( is.null(yb) ) {
    ya <- Compositional::alfa(y, a)$aff
  } else  ya <- yb
  x <- model.matrix( ya ~., data.frame(x) )

  esf <- spmoran::meigen_f(coords, model = model)
  X <- esf$sf

  if ( a <= 1e-5 ) {
    mod <- Compositional::comp.reg(y, X[, -1], yb = yb)
    be <- mod$be

  } else {
    ax <- a * x
    ha <- t( Compositional::helm(D) )
    ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
    suppressWarnings({
      mod <- minpack.lm::nls.lm(par = ini, fn = reg, jac = jac,
                                ya = ya, ax = ax, a = a, ha = ha, d = d, D = D)
    })

    be <- matrix(mod$par, ncol = d)
    est <- cbind(1, exp(x %*% be) )
    est <- est / Rfast::rowsums(est)
    kl <- y * log(y / est)
    kl[ is.infinite(kl) ] <- NA
    kld <- sum(kl, na.rm = TRUE)
    sela <- .areggompesf(y, x[, -1], X, a, ya, kld)[-1, 1]

    aX <- cbind(a * x, a * X[, sela] )
    ini <- as.vector( solve(crossprod(aX), crossprod(aX, ya) ) )
    suppressWarnings({
      mod <- minpack.lm::nls.lm(par = ini, fn = reg, jac =jac, ya = ya, ax = aX, a = a, ha = ha, d = d, D = D)
    })
    be <- matrix(mod$par, ncol = d)

  }  ## end if (a == 0)

  runtime <- proc.time() - runtime

  est <- NULL
  if ( !is.null(xnew) ) {
    esfnew <- NULL
    if ( length(sela) > 0 )  esfnew <- spmoran::meigen0(meig = esf, coords0 = coordsnew)$sf[, sela]
    xnew <- model.matrix(~., data.frame(xnew) )
    xnew <- cbind(xnew, esfnew)
    est <- cbind( 1, exp(xnew %*% be) )
    est <- est/Rfast::rowsums(est)
  }

  p <- dim(x)[2] - 1
  colnames(be) <- paste("Y", 2:D, sep = "")

  if ( length(sela) > 0 ) {
    if ( is.null( colnames(x) ) ) {
      rownames(be) <- c("constant", paste("X", 1:p, sep = ""), paste("SEF", sela, sep = "") )
    } else  rownames(be)  <- c("constant", colnames(x)[-1], paste("SEF", sela, sep = "") )
    gama <- be[-c( 1:(p + 1) ), , drop = FALSE]
    be <- be[1:(p + 1), ]
    X.esf <- X[, sela, drop = FALSE]
    colnames(X.esf) <- paste("ESF", sela, sep = "")
  } else {
    if ( is.null( colnames(x) ) ) {
      rownames(be) <- c("constant", paste("X", 1:p, sep = "") )
    } else  rownames(be) <- c("constant", colnames(x)[-1] )
    gama <- NULL
    X.esf <- NULL
  }

  if ( covb ) {
    res <- optim( as.vector( rbind(be, gama) ), .regar, ya = ya, ax = aX, a = a, ha = ha, d = d,
                  D = D, hessian = TRUE, method = "BFGS", control = list(maxit = 1000) )
    covbe <- solve(res$hessian)
    a2 <- c( colnames(x), colnames(X.esf) )
    a1 <- colnames(be)
    nam <- as.vector( t( outer(a1, a2, paste, sep = ":") ) )
    colnames(covbe) <- rownames(covbe) <- nam
  } else covbe <- NULL

  list(runtime = runtime, be = be, gama = gama, covbe = covbe, ESF = sela, X.esf = X.esf,
       dev = mod$deviance, est = est)
}




.areggompesf <- function(y, x, X, a, ya, kld) {

   KLD <- kld
   p <- dim(X)[2]
   klc <- numeric(p)
   ya <- Compositional::alfa(y, a)$aff
   ind <- 1:p
   for ( j in ind ) {
     klc[j] <- mean( Rfast::eachcol.apply(ya, X[, j])^2 )
   }
   sel <- which.max(klc)
   z <- cbind(x, X[, sel])
   mod <- CompositionalSR::areg(y, z, a, xnew = z)
   est <- mod$est
   kl <- y * log( y / est )
   kl[ is.infinite(kl) ] <- NA
   KLD <- c( KLD, sum(kl, na.rm = TRUE) )

   i <- 2
   while ( KLD[i] < KLD[i - 1] & length(sel) < p ) {
     i <- i + 1
     res <- ya - Compositional::alfa(est, a)$aff
     klc <- rep(-Inf, p)
     for ( j in ind[-sel] ) {
       klc[j] <- mean( Rfast::eachcol.apply(res, X[, j])^2 )
     }
     sel <- c(sel, which.max(klc) )
     z <- cbind(x, X[, sel])
     mod <- CompositionalSR::areg(y, z, a, xnew = z)
     est <- mod$est
     kl <- y * log( y / est )
     kl[ is.infinite(kl) ] <- NA
     KLD <- c( KLD, sum(kl, na.rm = TRUE) )
   }
   res <- cbind(c(0, sel), KLD)
   res <- res[-dim(res)[1], , drop = FALSE]
   colnames(res) <- c("Vars", "KLD")
   res
}














