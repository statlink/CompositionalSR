alfa.slx2 <- function(y, x, a, coords, k = 2:15, xnew, coordsnew, yb = NULL) {

  reg <- function(para, ya, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector(ya - ma)
  }

  runtime <- proc.time()

  D <- dim(y)[2]  ;  d <- D - 1
  if ( is.null(yb) ) {
    ya <- Compositional::alfa(y, a)$aff
  } else  ya <- yb
  ha <- t( Compositional::helm(D) )
  x <- model.matrix( ya ~., data.frame(x) )

  be <- gama <- est <- list()

  cords <- pi * coords / 180  ## from degrees to rads
  a1 <- sin(cords[, 1])
  coord <- cbind( cos(cords[, 1]), a1 * cos(cords[, 2]), a1 * sin(cords[, 2]) )
  W <- Rfast::Dist(coord, square = TRUE)
  B <- Rfast::rowOrder(W)
  W <- 1/W

  cordsnew <- pi * coordsnew / 180  ## from degrees to rads
  a1 <- sin(cordsnew[, 1])
  coordnew <- cbind( cos(cordsnew[, 1]), a1 * cos(cordsnew[, 2]), a1 * sin(cordsnew[, 2]) )
  cords <- pi * coords / 180  ## from degrees to rads
  a1 <- sin(cords[, 1])
  coord <- cbind( cos(cords[, 1]), a1 * cos(cords[, 2]), a1 * sin(cords[, 2]) )
  Wnew <- Rfast::dista(coordnew, coord, square = TRUE)
  Bnew <- Rfast::rowOrder(Wnew)
  Wnew <- 1 / Wnew
  Wnew[ is.infinite(Wnew) ] <- 0

  xnew <- model.matrix(~., data.frame(xnew) )

  for ( i in k )  {
    b <- B  
    b[b > i + 1] <- 0
    b[b > 0] <- 1
    w <- b * W
    b <- NULL
    diag(w) <- 0
    w <- w / Rfast::rowsums(w)

    wx <- w %*% x[, -1]
    X <- cbind(x, wx)
    if ( a <= 1e-5 ) {
      mod <- Compositional::comp.reg(y, X[, -1], yb = yb)
      bes <- mod$be

    } else {
      aX <- a * X
      ini <- as.vector( solve(crossprod(X), crossprod(X, ya) ) )
      suppressWarnings({
      mod <- minpack.lm::nls.lm( par = ini, fn = reg, ya = ya, ax = aX, a = a, ha = ha, d = d, D = D,
                                 control = minpack.lm::nls.lm.control(maxiter = 2048) )
      })
      bes <- matrix(mod$par, ncol = d)
    }  ## end if (a == 0)

    b <- Bnew
    b[b > i + 1] <- 0
    b[b > 0] <- 1
    wnew <- b * Wnew
    b <- NULL
    wnew <- wnew / Rfast::rowsums(wnew)

    wxnew <- wnew %*% x[, -1]
    Xnew <- cbind(xnew, wxnew)
    yest <- cbind( 1, exp(Xnew %*% bes) )
    yest <- yest/Rfast::rowsums(yest)

    if ( is.null( colnames(x) ) ) {
      p <- dim(x)[2] - 1
      rownames(bes) <- c("constant", paste("X", 1:p, sep = ""), paste("WX", 1:p, sep = "") )
    } else  rownames(bes)  <- c("constant", colnames(x)[-1], paste("W", colnames(x)[-1], sep = "") )

    gama[[ i ]] <- bes[(p + 2) : (2 * p + 1), ]
    be[[ i ]] <- bes[1:(p + 1), ]
    est[[ i ]] <- yest
  }

  runtime <- proc.time() - runtime

  list(runtime = runtime, be = be, gama = gama, est = est)
}



