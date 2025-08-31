

criminals <-
  function(x,
           y,
           xknots = knotsQ(x),
           yknots = knotsD(y),
           xdegrees = 2,
           ydegrees = -1,
           xordinal = TRUE,
           yordinal = FALSE,
           xcopies = 1,
           xties = "s",
           yties = "s",
           xmissing = "m",
           ymissing = "m",
           xactive = TRUE,
           xnames = colnames(x, do.NULL = FALSE),
           ynames = "Y",
           ndim = 2,
           itmax = 1000,
           eps = 1e-6,
           seed = 123,
           verbose = FALSE) {
    aname <- deparse(substitute(theData))
    npred <- ncol(x)
    nobs <- nrow(x)
    g <- makeGifi(
      theData = cbind(x, y),
      knots = c(xknots, yknots),
      degrees = c(reshape(xdegrees, npred), ydegrees),
      ordinal = c(reshape(xordinal, npred), yordinal),
      sets =  c(rep(1, npred), 2),
      copies = c(reshape(xcopies, npred), length(unique(y))),
      ties = c(reshape(xties, npred), yties),
      missing = c(reshape(xmissing, npred), ymissing),
      active = c(reshape(xactive, npred), TRUE),
      names = c(xnames, ynames)
    )
    h <- gifiEngine(
      gifi = g,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    x <- matrix(0, nobs, 0)
    for (j in 1:npred)
      x <- cbind(x, h$xGifi[[1]][[j]]$transform)
    y <- as.vector(y)
    g <- ifelse(outer(y, unique(y), "=="), 1, 0)
    d <- colSums(g)
    v <- crossprod(x)
    u <- crossprod(g, x)
    b <- crossprod(u, (1 / d) * u)
    w <- v - b
    e <- eigen(v)
    k <- e$vectors
    l <- sqrt(abs(e$values))
    l <- ifelse(l < 1e-7, 0, 1 / l)
    f <- eigen(outer(l, l) * crossprod(k, b %*% k))
    a <- k %*% (l * f$vectors)
    z <- x %*% a
    u <- (1 / d) * crossprod(g, x)
    return(structure(
      list(
        objectscores = h$x,
        xhat = x,
        loadings = a,
        scores = z,
        groupmeans = u,
        bwratios = f$values,
        ntel = h$ntel,
        f = h$f
      ),
      class = "criminals"
    ))
  }
