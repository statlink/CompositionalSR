################################
#### Tuning the alfa in alfa-regression via K-fold cross validation
#### Tibshirani and Tibshirani method
#### Tsagris Michail 11/2015
#### mtsagris@yahoo.gr
#### References: Tsagris Michail (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean Journal of Statistics, 6(2): 47-57
################################
cv.alfareg <- function( y, x, a = seq(0.1, 1, by = 0.1), nfolds = 10, folds = NULL, nc = 1 ) {
  ## y is the compositional data (dependent variable)
  ## x is the independent variables
  ## a is a range of values of alpha
  ## nfolds is the number of folds for the K-fold cross validation
  ## nc is how many cores you want to use, default value is 1
  if ( min(y) == 0 )  a <- a[a > 0]
  la <- length(a)
  n <- dim(y)[1]
  ina <- 1:n
  x <- model.matrix( y~., data = as.data.frame(x) )
  x <- x[, -1, drop = FALSE]
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds, stratified = FALSE)
  nfolds <- length(folds)

  ## cross-validated KLD for the j-th alpha value, one value per fold
  cv_alpha <- function(j) {
    ytr <- Compositional::alfa(y, a[j])$aff
    kula <- numeric(nfolds)
    for ( i in 1:nfolds ) {
      xtest <- x[ folds[[ i ]], , drop = FALSE]
      ytest <- y[ folds[[ i ]], , drop = FALSE]
      xtrain <- x[ -folds[[ i ]], , drop = FALSE]
      ytrain <- y[ -folds[[ i ]], , drop = FALSE]
      yb <- ytr[ -folds[[ i ]], , drop = FALSE]
      yest <- CompositionalSR::areg(ytrain, xtrain, a[j], xnew = xtest, yb = yb)$est
      kl <- ytest * log(ytest / yest)
      kl[ is.infinite(kl) ] <- NA
      kula[i] <- sum(kl, na.rm = TRUE) / dim(yest)[1]
    }
    kula
  }

  apa <- proc.time()
  if ( nc <= 1 ) {
    kula <- sapply(1:la, cv_alpha)
  } else {
    cl <- parallel::makePSOCKcluster(nc)
    on.exit( parallel::stopCluster(cl), add = TRUE )
    parallel::clusterExport(cl, varlist = c("y", "x", "a", "folds", "nfolds"),
                            envir = environment() )
    kula <- parallel::parSapply(cl, 1:la, cv_alpha)
  }
  kula <- matrix(kula, nrow = nfolds, ncol = la)   ## nfolds x la, also when nfolds = 1

  kl <- Rfast::colmeans(kula)
  opt <- c( min(kl), a[ which.min(kl) ] )
  names(opt) <- c( "KLD", "alpha")
  apa <- proc.time() - apa

  list(runtime = apa, perf = kl, opt = opt)
}
