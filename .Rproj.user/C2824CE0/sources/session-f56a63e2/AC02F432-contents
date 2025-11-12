rob.alfareg <- function(y, x, a, loss = "welsh", xnew = NULL, yb = NULL) {

  reg <- function(para, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector(ma)
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
    ha <- t( Compositional::helm(D) )
    ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
    suppressWarnings({
    mod <- gslnls::gsl_nls( start = ini, y = as.vector(ya), fn = reg, ax = ax, a = a, ha = ha, d = d, D = D,
                            control = gslnls::gsl_nls_control(maxiter = 10000), algorithm = "lm", loss = loss )
    })
    be <- matrix(coef(mod), ncol = d)

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
  } else rownames(be)  <- c("constant", colnames(x)[-1] )
  colnames(be) <- paste("Y", 2:D, sep = "")

  list(runtime = runtime, be = be, est = est)
}


