contiguity <- function(coords, k = 10) {
  #colnames(coords) <- c("latitude", "longitude") 
  cords <- pi * coords / 180  ## from degrees to rads
  a1 <- sin(cords[, 1])
  coord <- cbind( cos(cords[, 1]), a1 * cos(cords[, 2]), a1 * sin(cords[, 2]) )
  W <- Rfast::Dist(coord, square = TRUE)
  b <- Rfast::rowOrder(W)
  b[b > k + 1] <- 0 
  b[b > 0] <- 1
  W <- 1 / W
  W <- b * W
  b <- NULL
  diag(W) <- 0
  W / Rfast::rowsums(W)
}