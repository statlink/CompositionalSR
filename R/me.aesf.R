me.aesf <- function(be, gama, mu, x, x.esf) {
  be <- rbind(be, gama)
  x <- cbind(x, x.esf)
  CompositionalSR::me.ar(be, mu, x)
}
