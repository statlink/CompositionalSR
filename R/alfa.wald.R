alfa.wald <- function(y, x, a, R = 299, ncores) {

  n <- dim(y)[1]
  yb <- Compositional::alfa(y, a)$aff
  mod <- CompositionalSR::areg(y, x, a, covb = TRUE, yb = yb)
  b <- as.vector(mod$be[-1, ])
  ep <- grep("Intercept", colnames( mod$covbe) )
  s <- mod$covbe[-ep, -ep]
  stat <- as.numeric( b %*% solve(s, b) )

  cl <- parallel::makeCluster(ncores)
  parallel::clusterEvalQ(cl, {  # CHANGE 2: UNCOMMENT AND ADD PACKAGES
    library(Compositional)
    library(CompositionalSR)
  })
  parallel::clusterExport(cl, varlist = c("x", "y", "yb", "a", "ep", "n"), envir = environment())
  statb <- parallel::parLapply(cl, 1:R, function(j) {
    ind <- sample(n, n)
    modb <- CompositionalSR::areg(y, x, a = a, covb = TRUE, yb = yb[ind, ])
    s <- modb$covbe[-ep, -ep]
    b <- as.vector(modb$be[-1, ])
    as.numeric( b %*% solve(s, b) )
  })
  parallel::stopCluster(cl)
  statb <- unlist(statb)
  ( sum(statb >= stat) + 1 ) / R
}
