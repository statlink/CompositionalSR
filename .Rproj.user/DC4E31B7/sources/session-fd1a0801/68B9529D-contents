gwar <- function(y, x, a, coords, h, yb = NULL, nc = 1) {

    regwei <- function(para, ya, ax, a, wei, ha, d, D) {
      be <- matrix(para, ncol = d)
      zz <- cbind( 1, exp(ax %*% be) )
      ta <- rowSums(zz)
      za <- zz / ta
      ma <- ( D / a * za - 1/a ) %*% ha
      as.vector( wei * (ya - ma ) )
    }

    ## fits the local model at location i and returns c(be, est)
    locfit <- function(i, coord, ya, x, ax, a, h, ha, d, D, ini, regwei) {
      di <- Rfast::eachcol.apply( coord, coord[, i] ) - 1
      wei <- sqrt( exp(di / h) )
      ind <- which(wei > 1e-6)
      suppressWarnings( {
        mod <- minpack.lm::nls.lm( par = ini, fn = regwei, ya = ya[ind, ], ax = ax[ind, ], a = a,
                                   wei = wei[ind], ha = ha, d = d, D = D )
      } )
      be2 <- matrix(mod$par, ncol = d)
      est2 <- cbind( 1, exp(x[i, ] %*% be2) )
      c( as.vector(mod$par), est2 / sum(est2) )
    }
    ## detach the helpers from gwar's environment so that only what is passed
    ## as arguments is shipped to the workers (not the whole local frame)
    environment(regwei) <- globalenv()
    environment(locfit) <- globalenv()

    runtime <- proc.time()

    colnames(coords) <- c("latitude", "longitude")
    cords <- pi * coords / 180  ## from degrees to rads
    a1 <- sin(cords[, 1])
    coord <- rbind( cos(cords[, 1]), a1 * cos(cords[, 2]), a1 * sin(cords[, 2]) )

    if ( is.null(yb) ) {
      ya <- Compositional::alfa(y, a)$aff
    } else  ya <- yb
    x <- model.matrix( ya ~., data.frame(x) )
    ax <- a * x

    D <- dim(y)[2]  ;  p <- dim(x)[2]
    d <- D - 1  ## dimensionality of the simplex
    n <- dim(y)[1]  ## sample size
    ha <- t( Compositional::helm(D) )
    ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
    h <- h^2

    if ( nc <= 1 ) {
      ep <- lapply( 1:n, locfit, coord = coord, ya = ya, x = x, ax = ax, a = a, h = h,
                    ha = ha, d = d, D = D, ini = ini, regwei = regwei )
    } else {
      cl <- parallel::makePSOCKcluster(nc)
      on.exit( parallel::stopCluster(cl), add = TRUE )
      ep <- parallel::parLapply( cl, 1:n, locfit, coord = coord, ya = ya, x = x, ax = ax, a = a,
                                 h = h, ha = ha, d = d, D = D, ini = ini, regwei = regwei )
    }
    ep <- do.call(rbind, ep)
    be <- ep[, 1:(p * d), drop = FALSE]
    est <- ep[, -c( 1:(p * d) ), drop = FALSE]

    runtime <- proc.time() - runtime

    list(runtime = runtime, be = be, est = est)
  }
