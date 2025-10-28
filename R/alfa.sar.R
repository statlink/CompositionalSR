alfa.sar <- function(y, x, a, coords, k = 10, xnew = NULL, coordsnew, yb = NULL) {

  reg <- function(para, ystar, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector( ystar - x %*% be )
  }

  regalr <- function(para, ystar, x, d) {
    be <- matrix(para, ncol = d)
    ma <- x %*% be
    as.vector( ystar - x %*% be )
  }

  runtime <- proc.time()
  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  n <- dim(y)[1]

  W <- CompositionalSR::contiguity(coords, k)
  if ( is.null(yb) ) {
    ya <- Compositional::alfa(y, a)$aff
  } else  ya <- yb
  x <- model.matrix( ya ~., data.frame(x) )
  p <- dim(x)[2]

  # Log-determinant (efficient!)
  values <- eigen(W, symmetric = FALSE, only.values = TRUE)$values
  ini <- rnorm(d * p)

  rega <- function(rho, ya, Wya, ax, a, ha, d, D, n, ini) {
    ystar <- ya - rho * Wya
    suppressWarnings({
      mod <- minpack.lm::nls.lm( par = ini, fn = reg, ystar = ystar, ax = ax, a = a, ha = ha, d = d, D = D,
                                 control = minpack.lm::nls.lm.control(maxiter = 5000) )
    })
    0.5 * n * log(mod$deviance / n) - d * sum( log( Mod(1 - rho * values) ) )
  }

  rega0 <- function(rho, ya, Wya, x, d, n, ini) {
    ystar <- ya - rho * Wya
    suppressWarnings({
      mod <- minpack.lm::nls.lm( par = ini, fn = regalr, ystar = ystar, x = x, d = d,
                                 control = minpack.lm::nls.lm.control(maxiter = 5000) )
    })
    0.5 * n * log(mod$deviance / n) - d * sum( log( Mod(1 - rho * values) ) )
  }

  if ( a <= 1e-5 ) {
    if ( is.null(yb) ) {
      ya <- Compositional::alr(y)
    } else  ya <- yb
      Wya <- W %*% ya
      rho <- optimize(rega0, c(0, 1), ya = ya, Wya = Wya, x = x, d = d, n = n, ini = ini)$minimum
      ystar <- ya - rho * Wya
      suppressWarnings({
        mod <- minpack.lm::nls.lm( par = ini, fn = regalr, ystar = ystar, x = x, d = d,
                                   control = minpack.lm::nls.lm.control(maxiter = 5000) )
      })
    be <- matrix(mod$par, ncol = d)

  } else {
    ax <- a * x
    ha <- t( Compositional::helm(D) )
    Wya <- W %*% ya
    rho <- optimize( rega, c(0, 1), ya = ya, Wya = Wya, ax = ax, a = a, ha = ha, d = d,
                     D = D, n = n, ini = ini )$minimum
    ystar <- ya - rho * Wya
    suppressWarnings({
      mod <- minpack.lm::nls.lm( par = ini, fn = reg, ystar = ystar, ax = ax, a = a, ha = ha, d = d,
                                 D = D, control = minpack.lm::nls.lm.control(maxiter = 5000) )
    })
    be <- matrix(mod$par, ncol = d)
  }  ## end if (a == 0)

  runtime <- proc.time() - runtime

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )

    cordsnew <- pi * coordsnew / 180  ## from degrees to rads
    a1 <- sin(cordsnew[, 1])
    coordnew <- cbind( cos(cordsnew[, 1]), a1 * cos(cordsnew[, 2]), a1 * sin(cordsnew[, 2]) )
    cords <- pi * coords / 180  ## from degrees to rads
    a1 <- sin(cords[, 1])
    coord <- cbind( cos(cords[, 1]), a1 * cos(cords[, 2]), a1 * sin(cords[, 2]) )
    Wnew <- Rfast::dista(coordnew, coord, square = TRUE)

    b <- Rfast::rowOrder(Wnew)
    b[b > k + 1] <- 0
    b[b > 0] <- 1
    Wnew <- 1 / Wnew
    Wnew[ is.infinite(Wnew) ] <- 0
    Wnew <- b * Wnew
    b <- NULL
    Wnew <- Wnew / Rfast::rowsums(Wnew)
    mu <- cbind( 1, exp(xnew %*% be) )
    mu <- mu/Rfast::rowsums(mu)
    est <- ( rho * Wnew %*% y + mu ) / (1 + rho)
  }

  if ( is.null( colnames(x) ) ) {
    p <- dim(x)[2] - 1
    rownames(be) <- c("constant", paste("X", 1:p, sep = "") )
  } else rownames(be)  <- c("constant", colnames(x)[-1] )

  list(runtime = runtime, be = be, rho = rho, est = est)
}



