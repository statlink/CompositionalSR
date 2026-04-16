alfa.wald <- function(y, x, a, R = 499, ncores) {

cl <- parallel::makeCluster(ncores)
n <- dim(x)[1]
colnames(x) <- paste("X", 1:p)
yb <- Compositional::alfa(y, a)$aff
mod <- CompositionalSR::areg(y, x, a, covb = TRUE, yb = yb)
b <- as.vector(mod$be[-1, ])
ep <- grep("Intercept", colnames(mod$covbe))
s <- mod$covbe[-ep, -ep]
stat <- as.vector(b %*% solve(s, b))
parallel::clusterExport(cl, c("x", "y", "yb", "mod", "ep", "s", "stat"), envir = environment())
  statb <- parallel::parLapply(cl, 1:B, function(j) {
  ind <- sample(n, n, replace = TRUE)
  modb <- CompositionalSR::areg(y, x[ind, ], a = a, covb = TRUE, yb = yb)
  s <- mod$covbe[-ep, -ep]
  b <- as.vector(modb$be[-1, ])
  as.numeric(b %*% solve(s, b))
})

statb <- unlist(statb)
( sum(statb >= stat) + 1 ) / R
}
