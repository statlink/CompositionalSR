alfa.wald <- function(y, x, a, index = "all", R = 299, ncores) {

  n <- dim(y)[1]
  ya <- Compositional::alfa(y, a)$aff
  mod <- CompositionalSR::areg(y, x, a, covb = TRUE, yb = ya)

  if ( index == "all" ) {
    b <- as.vector(mod$be[-1, ])
    ep <- grep("Intercept", colnames( mod$covbe) )
    s <- mod$covbe[-ep, -ep]
    stat <- as.numeric( b %*% solve(s, b) )

    cl <- parallel::makeCluster(ncores)
    parallel::clusterEvalQ(cl, {  # CHANGE 2: UNCOMMENT AND ADD PACKAGES
      library(Compositional)
      library(CompositionalSR)
    })
    parallel::clusterExport(cl, varlist = c("x", "y", "ya", "a", "ep", "n"), envir = environment())
    statb <- parallel::parLapplyLB(cl, 1:R, function(j) {
      tryCatch({
        ind <- sample(n, n)
        modb <- CompositionalSR::areg(y, x, a = a, covb = TRUE, yb = ya[ind, ])
        s <- modb$covbe[-ep, -ep]
        b <- as.vector(modb$be[-1, ])
        list( value = as.numeric(b %*% solve(s, b)), error = NULL )
      }, error = function(e) {
        list( value = NA_, error = conditionMessage(e) )
      })
    })
    parallel::stopCluster(cl)
    statb <- unlist(statb)
    pvalue <- ( sum(statb >= stat, na.rm = TRUE) + 1 ) / (R - sum( is.na(statb) ) + 1)

  } else {
    b <- as.vector(mod$be[index + 1, ])
    p <- dim(x)[2]
    colnames(x) <- paste("X", 1:p, sep = "")
    ep <- grep(paste("X", index, sep = ""), colnames( mod$covbe) )
    s <- mod$covbe[ep, ep]
    stat <- as.numeric( b %*% solve(s, b) )

    mod <- CompositionalSR::areg(y, x[, -index], a, covb = TRUE, yb = ya, xnew = x[, -index])
    esta <- Compositional::alfa(mod$est, a)$aff
    res <- ya - esta

    cl <- parallel::makeCluster(ncores)
    parallel::clusterEvalQ(cl, {  # CHANGE 2: UNCOMMENT AND ADD PACKAGES
      library(Compositional)
      library(CompositionalSR)
    })
    parallel::clusterExport(cl, varlist = c("x", "y", "a", "res", "ep", "n"), envir = environment())
    statb <- parallel::parLapplyLB(cl, 1:R, function(j) {
      tryCatch({
        ind <- sample(n, n)
        resb <- res[ind, ]
        yb <- esta + resb
        modb <- CompositionalSR::areg(y, x, a = a, covb = TRUE, yb = yb)
        s <- modb$covbe[ep, ep]
        b <- as.vector(modb$be[index + 1, ])
        list( value = as.numeric(b %*% solve(s, b)), error = NULL )
      }, error = function(e) {
        list( value = NA_, error = conditionMessage(e) )
      })
    })
    parallel::stopCluster(cl)
    statb <- unlist(statb)
    pvalue <- ( sum(statb >= stat, na.rm = TRUE) + 1 ) / (R - sum( is.na(statb) ) + 1)

  }
  pvalue
}
