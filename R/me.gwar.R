me.gwar <- function(be, mu, x) {

  n <- dim(mu)[1]  ;  D <- dim(mu)[2]
  d <- D - 1
  p <- NCOL(x)

  mfx_array <- array(0, dim = c(n, D, p))

  for (i in 1:n) {
    beta_loc <- t( matrix( be[i, ], ncol = d) )[, -1, drop = FALSE]  # d x p
    mu_loc <- mu[i, ]        # D

    # Loop over covariates
    for ( k in 1:p ) {
      com <- 0
      for (j in 1:d) {
        com <- com + beta_loc[j, k] * mu_loc[j+1]
      }
      # Reference component
      mfx_array[i, 1, k] <- -mu_loc[1] * com
      # Other components
      for ( comp in 2:D ) {
        mfx_array[i, comp, k ] <- mu_loc[comp] * (beta_loc[comp - 1, k] - com)
      }
    }
  }

  dimnames(mfx_array) <- list( location = paste0("Loc_", 1:n),
    component = paste0("Component_", 1:D), covariate = colnames(x) )
  if ( is.null( dimnames(mfx_array)[[3]] ) ) {
    dimnames(mfx_array)[[ 3 ]] <- paste0("X", 1:p)
  }

  list( me = mfx_array, ame = t( apply(mfx_array, c(2, 3), mean) ) )
}
