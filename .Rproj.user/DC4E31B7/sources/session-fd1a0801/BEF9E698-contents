gwar.pred <- function(y, x, a, coords, h, xnew, coordsnew) {

  regwei <- function(para, ya, ax, a, wei, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector( wei * (ya - ma ) )
  }

  runtime <- proc.time()

  ##colnames(coords) <- c("latitude", "longitude")
  cords <- pi * coords / 180  ## from degrees to rads
  a1 <- sin(cords[, 1])
  coord <- rbind( cos(cords[, 1]), a1 * cos(cords[, 2]), a1 * sin(cords[, 2]) )

  if ( !is.matrix(coordsnew) )  coordsnew <- matrix(coordsnew, ncol = 2)
  cordsnew <- pi * coordsnew / 180  ## from degrees to rads
  a1 <- sin(cordsnew[, 1])
  coordnew <- rbind( cos(cordsnew[, 1]), a1 * cos(cordsnew[, 2]), a1 * sin(cordsnew[, 2]) )

  xnew <- as.matrix(xnew)
  if ( dim(xnew)[1] == 1 )  xnew <- t(xnew)
  nu <- dim(xnew)[1]  ## sample size

  la <- length(a)  ;  lh <- length(h)
  D <- dim(y)[2]  ;  d <- D - 1
  ha <- t( Compositional::helm(D) )
  x <- model.matrix( y ~., data.frame(x) )
  xxinv <- solve( crossprod(x) )

  est <- list()
  names <- paste("alpha=", a, sep = "")
  est <- sapply(names, function(x) NULL)
  for ( i in 1:la )  {
    est[[ i ]] <- list()
    for (j in 1:lh)  est[[ i ]][[ j ]] <- matrix(nrow = nu, ncol = D)
    names(est[[ i ]]) <- paste("h=", h, sep = "")
  }
  h <- h^2

  for ( i in 1:nu ) {
    di <- Rfast::eachcol.apply( coord, coordnew[, i] ) - 1
    for ( k in 1:la ) {
      ya <- Compositional::alfa(y, a[k])$aff
      ax <- a[k] * x
      ini <- as.vector( crossprod( xxinv, crossprod(x, ya) ) )
      for ( m in 1:lh ) {
        wei <- sqrt( exp(di / h[m]) )
        ind <- which(wei > 1e-6)
        suppressWarnings( {
        mod <- minpack.lm::nls.lm( par = ini, fn = regwei, ya = ya[ind, ], ax = ax[ind, ], a = a, wei = wei[ind],
                                   ha = ha, d = d, D = D, control = minpack.lm::nls.lm.control(maxiter = 10000) )
        } )
        be <- matrix(mod$par, ncol = d)
        est2 <- cbind( 1, exp(x[i, ] %*% be) )
        est[[ k ]][[ m ]][i, ] <- est2 / Rfast::rowsums(est2)
      }
    }
  }

  runtime <- proc.time() - runtime

  list(runtime = runtime, est = est)
}

