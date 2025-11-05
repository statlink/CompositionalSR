################################
#### alfa-regression
#### Tsagris Michail 11/2015
#### mtsagris@yahoo.gr
#### References: Tsagris Michail (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean Journal of Statistics, 6(2): 47-57
################################
alfa.reg <- function(y, x, a, xnew = NULL, yb = NULL) {

  reg <- function(para, ya, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector(ya - ma)
  }

  runtime <- proc.time()

  if ( is.null(yb) ) {
    ya <- Compositional::alfa(y, a)$aff
  } else  ya <- yb
  x <- model.matrix(ya ~., data.frame(x) )
  ax <- a * x

  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex

  if ( a == 0 ) {
    mod <- Compositional::comp.reg(y, x[, -1], yb = yb)
    be <- mod$be
    if ( !is.null(seb) )  seb <- mod$seb
    runtime <- mod$runtime

  } else {
    ha <- t( Compositional::helm(D) )
    ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
    suppressWarnings({
      mod <- minpack.lm::nls.lm( par = ini, fn = reg, ya = ya, ax = ax, a = a, ha = ha, d = d, D = D,
                                 control = minpack.lm::nls.lm.control(maxiter = 10000) )
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
  } else rownames(be)  <- c("constant", colnames(x)[-1] )
  colnames(be) <- paste("Y", 2:D, sep = "")  

  list(runtime = runtime, be = be, dev = mod$deviance, est = est)
}



