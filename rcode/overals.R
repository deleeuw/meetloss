
overals <-
  function(data,
           sets,
           copies,
           knots = knotsQ (data),
           degrees = rep (2, ncol (data)),
           ordinal = rep (TRUE, ncol (data)),
           ndim = 2,
           itmax = 1000,
           eps = 1e-6,
           seed = 123,
           verbose = FALSE)  {
    h <- gifiEngine(
      data = data,
      knots = knots,
      degrees = degrees,
      ordinal = ordinal,
      sets = sets,
      copies = copies,
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
    for (j in 1:ncol(data)) {
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


