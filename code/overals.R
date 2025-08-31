
overals <-
    function(theData,
             sets, 
             knots = knotsQ(theData),
             degrees = 2,
             ordinal = TRUE,
             copies = 1,
             ndim = 2,
             ties = "s",
             missing = "m",
             names = colnames(theData, do.NULL = FALSE),
             active = TRUE,
             itmax = 1000,
             eps = 1e-6,
             seed = 123,
             verbose = FALSE) {
    aname <- deparse(substitute(theData))
    nvars <- ncol(theData)
    nobs <- nrow(theData)
    g <- makeGifi(theData = theData,
               knots = knots,
               degrees = degrees,
               ordinal = reshape(ordinal, nvars),
               ties = reshape(ties, nvars),
               copies = reshape(copies, nvars),
               missing = reshape(missing, nvars),
               active = reshape(active, nvars),
               names = names,
               sets = sets) 
    h <- gifiEngine(
      gifi = g,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    xhat <- h$h
    rhat <- cor (xhat)
    a <- h$a
    y <- xhat
    for (j in 1:ncol(theData)) {
      k <- (1:ndim) + (j - 1) * ndim
      y[, k] <- xhat[, k] %*% a[k, ]
    }
    return (structure (
      list (
        xhat = xhat,
        rhat = rhat,
        objscores = h$x,
        quantifications = y,
        ntel = h$ntel,
        f = h$f
      ),
      class = "overals"
    ))
  }


