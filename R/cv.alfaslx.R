cv.alfaslx <- function(y, x, a = seq(0.1, 1, by = 0.1), coords, k = 2:15, nfolds = 10, folds = NULL) {
  if ( min(y) == 0 )  a <- a[a>0]
  la <- length(a)
  lk <- length(k)
  n <- dim(y)[1]
  ina <- 1:n
  x <- as.matrix(x)

    apa <- proc.time()
    if ( is.null(folds) )  folds <- CompositionalSR::spat.folds(coords, nfolds = nfolds)

    nfolds <- length(folds)
    kul <- matrix(nrow = nfolds, ncol = lk)
    kula <- matrix(nrow = la, ncol = lk)
    rownames(kula) <- paste("alpha=", a, sep = "")
    colnames(kula) <- paste("k=", k, sep = "")

    for ( i in 1:la ) {
      ytr <- Compositional::alfa(y, a[i])$aff
      for ( m in 1:nfolds ) {
        xtrain <- x[folds[[ m ]][[ 1 ]], ]
        ytrain <- y[ folds[[ m ]][[ 1 ]], ]
        xtest <- x[ folds[[ m ]][[ 2 ]], ]
        ytest <- y[ folds[[ m ]][[ 2 ]], ]
        coordstrain <- coords[folds[[ m ]][[ 1 ]], ]
        coordstest <- coords[folds[[ m ]][[ 2 ]], , drop = FALSE]
        for ( j in 1:lk ) {
          mod <- CompositionalSR::alfa.slx( ytrain, xtrain, a[i], coords = coordstrain, k = k[j], xnew = xtest,
                                            coordsnew = coordstest, yb = ytr[folds[[ m ]][[ 1 ]], ] )
          yest <- mod$est
          kl <- ytest * log(ytest / yest)
          kl[ is.infinite(kl) ] <- NA
          kul[m, j] <- sum(kl, na.rm = TRUE) / dim(ytest)[1]
        }
      }
      kula[i, ] <- Rfast::colmeans(kul)
    }

    apa <- proc.time() - apa
    best <- which( kula == min(kula), arr.ind = TRUE )
    opt <- c( min(kula), a[ best[1] ], k[ best[2] ] )
    names(opt) <- c( "KLD", "alpha", "k")

  list(runtime = apa, perf = kula, opt = opt)
}



