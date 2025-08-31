

morals <-
  function(x,
           y,
           xknots = knotsQ(x),
           yknots = knotsQ(y),
           xdegrees = 2,
           ydegrees = 2,
           xordinal = TRUE,
           yordinal = TRUE,
           xties = "s",
           yties = "s",
           xmissing = "m",
           ymissing = "m",
           xnames = colnames (x, do.NULL = FALSE),
           ynames = "Y",
           xactive = TRUE,
           xcopies = 1,
           itmax = 1000,
           eps = 1e-6,
           seed = 123,
           verbose = FALSE) {
    npred <- ncol (x)
    nobs <- nrow (x)
    xdegrees <- reshape (xdegrees, npred)
    xordinal <- reshape (xordinal, npred)
    xties <- reshape (xties, npred)
    xmissing <- reshape (xmissing, npred)
    xactive <- reshape (xactive, npred)
    xcopies <- reshape (xcopies, npred)
    g <- makeGifi (
      theData = cbind (x, y),
      knots = c (xknots, yknots),
      degrees = c (xdegrees, ydegrees),
      ordinal = c (xordinal, yordinal),
      sets =  c (rep(1, npred), 2),
      copies = c (xcopies, 1),
      ties = c (xties, yties),
      missing = c (xmissing, ymissing),
      active = c (xactive, TRUE),
      names = c (xnames, ynames)
    )
    h <- gifiEngine(
      gifi = g,
      ndim = 1,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    xhat <- matrix (0, nobs, 0)
    for (j in 1:npred)
      xhat <- cbind (xhat, h$xGifi[[1]][[j]]$transform)
    yhat <- h$xGifi[[2]][[1]]$transform
    rhat <- cor (cbind (xhat, yhat))
    qxy <- lsRC(xhat, yhat)$solution
    ypred <- xhat %*% qxy
    yres <- yhat - ypred
    smc <- sum (yhat * ypred)
    return (structure (
      list (
        objscores = h$x,
        xhat = xhat,
        yhat = yhat,
        rhat = rhat,
        beta = qxy,
        ypred = ypred,
        yres = yres,
        smc = smc,
        ntel = h$ntel,
        f = h$f
      ),
      class = "morals"
    ))
  }

