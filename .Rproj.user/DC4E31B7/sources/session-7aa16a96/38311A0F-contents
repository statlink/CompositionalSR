alfa.pcreg <- function(y, x, a, k, xnew = NULL, yb = NULL) {
 
  runtime <- proc.time()
  if ( is.null(yb) ) {
    ya <- Compositional::alfa(y, a)$aff
  } else  ya <- yb
  xa <- Compositional::alfa(x, a)$aff
        
  pca <- prcomp(xa, center = FALSE, scale. = FALSE, rank. = k)
  xa <- pca$x
  if ( !is.null(xnew) )  xnew <- Compositional::alfa(xnew, a)$aff %*% pca$rotation
  mod <- CompositionalSR::areg(y, xa, a, xnew = xnew, yb = ya)
  runtime <- proc.time() - runtime  
  mod$runtime <- runtime  
  mod
}

