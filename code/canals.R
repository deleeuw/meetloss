
canals <-
  function(x,
           y,
           xknots = knotsQ(x),
           yknots = knotsQ(y),
           xdegrees = rep(2, ncol(x)),
           ydegrees = rep(2, ncol(y)),
           xordinal = rep (TRUE, ncol (x)),
           yordinal = rep (TRUE, ncol (y)),
           xcopies = rep (1, ncol (x)),
           ycopies = rep (1, ncol (y)),
           ndim = 2,
           itmax = 1000,
           eps = 1e-6,
           seed = 123,
           verbose = FALSE) {
    h <- gifiEngine(
      data = cbind (x, y),
      knots = c(xknots, yknots),
      degrees = c(xdegrees, ydegrees),
      ordinal = c(xordinal, yordinal),
      sets = c(rep(1, ncol(x)), rep(2, ncol (y))),
      copies = c(xcopies, ycopies),
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    x <- h$h[, 1:sum(xcopies)]
    y <- h$h[, -(1:sum(xcopies))]
    u <- crossprod (x)
    v <- crossprod (y)
    w <- crossprod (x, y)
    a <- solve (chol (u))
    b <- solve (chol (v))
    s <- crossprod (a, w %*% b)
    r <- svd (s)
    xw <- a %*% (r$u)
    yw <- b %*% (r$v)
    xs <- x %*% xw
    ys <- y %*% yw
    xl <- crossprod (x, xs)
    yl <- crossprod (y, ys)
    return (structure (
      list(
        xhat = x,
        yhat = y,
        rhat = cor (cbind (x, y)),
        cancors = r$d,
        xweights = xw,
        yweights = yw,
        xscores = xs,
        yscores = ys,
        xloadings = xl,
        yloadings = yl,
        ntel = h$ntel,
        f = h$f
      ),
      class = "canals"
    ))
  }
