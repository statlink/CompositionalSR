aslx.grads <- function(y, x, a, be, gama, coords, k = 10) {

  W <- CompositionalSR::contiguity(coords, k)
  Wx <- W %*% x
  nam <- colnames(x)
  if ( is.null( nam ) ) {
    p <- NCOL(x) 
    nam <- c( paste("X", 1:p, sep = ""), paste("WX", 1:p, sep = "") ) 
  } else  nam <- c(nam, paste("W", nam, sep = "") )
  X <- cbind(x, Wx)
  colnames(X) <- nam
  be <- rbind(be, gama)
  CompositionalSR::ar.grads(y, X, a, be)
}
