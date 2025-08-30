
princals <-
  function(data,
           knots = knotsQ(data),
           degrees = 2,
           ordinal = TRUE,
           copies = 1,
           ndim = 2,
           ties = "s",
           missing = "m",
           names = colnames(data, do.NULL = FALSE),
           active = TRUE,
           itmax = 1000,
           eps = 1e-6,
           seed = 123,
           verbose = FALSE)  {
    aname <- deparse (substitute (data))
    nvars <- ncol (data)
    nobs <- nrow (data)
    g <- makeGifi (
      data = data,
      knots = knots,
      degrees = reshape (degrees, nvars),
      ordinal = reshape (ordinal, nvars),
      sets =  1:nvars,
      copies = reshape (copies, nvars),
      ties = reshape (ties, nvars),
      missing = reshape (missing, nvars),
      active = reshape (active, nvars),
      names = names
    )
    h <- gifiEngine(
      gifi = g,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    a <- v <- z <- d <- y <- o <- as.list (1:nvars)
    dsum <- matrix (0, ndim, ndim)
    for (j in 1:nvars) {
      jgifi <- h$xGifi[[j]][[1]]
      v[[j]] <- jgifi$transform
      a[[j]] <- jgifi$weights
      y[[j]] <- jgifi$scores
      z[[j]] <- jgifi$quantifications
      cy <- crossprod (y[[j]])
      dsum <- dsum + cy
      d[[j]] <- cy
      o[[j]] <- crossprod (h$x, v[[j]])
    }
    return (structure (
      list (
        transform = v,
        rhat = corList (v),
        objectscores = h$x,
        scores = y,
        quantifications = z,
        dmeasures = d,
        lambda = dsum / ncol (data),
        weights = a,
        loadings = o,
        ntel = h$ntel,
        f = h$f
      ),
      class = "princals"
    ))
  }
