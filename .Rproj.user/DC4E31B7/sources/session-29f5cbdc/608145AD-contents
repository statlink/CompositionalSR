cv.alfapcreg <- function( y, x, a = seq(0.1, 1, by = 0.1), k = dim(x)[2] - 2, nfolds = 10, folds = NULL, seed = NULL ) {
  if ( min(y) == 0 )  a <- a[ a > 0 ]
  la <- length(a)
  n <- dim(y)[1]
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  apa <- proc.time()
  kula <- matrix(nrow = nfolds, ncol = k)
  akula <- matrix(nrow = la, ncol = k)
  rownames(akula) <- paste("alpha=", a, sep = "")
  colnames(akula) <- paste("PC", 1:k, sep = "")

  for ( j in 1:la ) {
    ytr <- Compositional::alfa(y, a[j])$aff
    xtr <- Compositional::alfa(x, a[j])$aff

    for ( i in 1:nfolds ) {
      ytrain <- y[-folds[[ i ]], ]
      yb <- ytr[ -folds[[ i ]], ]
      pca <- prcomp(xtr[ -folds[[ i ]], ], center = FALSE, scale. = FALSE)
      ytest <- y[ folds[[ i ]], ]
      for ( vim in 1:k ) {
        xtrain <- pca$x[, 1:vim]
        xtest <- xtr[ folds[[ i ]], , drop = FALSE] %*% pca$rotation[, 1:vim]
        yest <- CompositionalSR::areg(ytrain, xtrain, a[j], xnew = xtest, yb = yb)$est
        kl <- ytest * log(ytest / yest)
        kl[ is.infinite(kl) ] <- NA
        kula[i, vim] <- sum(kl, na.rm = TRUE) / dim(yest)[1]
      }
    }
    akula[j, ] <- Rfast::colmeans(kula)
  }

  opt <- which(akula == min(akula), arr.ind = TRUE)
  apa <- proc.time() - apa

  list(runtime = apa, perf = akula, kl = min(akula), opt_a = a[opt[, 1]], opt_k = opt[, 2] )
}
