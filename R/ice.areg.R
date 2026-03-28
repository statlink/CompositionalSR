ice.areg <- function(be, x, ind = 1, frac = 0.2, pos = 0.5) {

  dm <- dim(be)[2] + 1
  n <- dm[1]  ;  d <- dm[2]
  nu <- ceiling(frac * n)
  x <- cbind(1, x)
  xsel <- sort( rangen::Sample(x[, ind + 1], nu) )
  est <- matrix(NA, nu, d)

  for (i in 1:nu) {
    X <- x
    X[, ind + 1] <- xsel[i]
    mu <- cbind( 1, exp(X %*% be) )
    mu <- mu/Rfast::rowsums(mu)
    est[i, ] <- Rfast::colmeans(mu)
  }

  namx <- colnames(x)[ind]
  if ( is.null(namx) )  namx <- paste("X", ind, sep = "")

  ##png(filename = "ice.png", width = 5000, height = 4000, res = 600)

  plot( xsel, est[, 1], type = "l", xlab = namx, ylab = "Fitted proportions",
        xaxt = "n", yaxt = "n", cex.lab = 1.2, cex.axis = 1.2,
        ylim = c( min(est), max(est) ), lwd = 2, bty = "n" )

  v <- seq( min(xsel), max(xsel), length = 10 )
  h <- seq( min(est), max(est), length = 10 )
  mtext( text = round(h, 2), side = 2, at = h, las = 2, font = 2, line = 0.2 )
  mtext( text = round(v, 2), side = 1, at = v, las = 2, font = 2, line = 0.2 )

  for (i in 1:10) {
    lines(rep(v[i], 10), h, col = "lightgrey", lty = 2)
    lines(v, rep(h[i], 10), col = "lightgrey", lty = 2)
  }
  for ( i in 1:d ) lines(xsel, est[, i], col = i, lwd = 2)

  namy <- paste("Y", 1:d, sep = "")
  if ( is.null(namy) )  namy <- paste("Comp. ", 1:d, sep = "")

  legend( x = quantile(x[, ind], pos), y = max(est), legend = namy,
          xpd = TRUE, bty = "n", title = "Components", col = 1:d,
          lwd = rep(2, d), lty = rep(1, d), text.col = 1:d )

  ##dev.off()

}
