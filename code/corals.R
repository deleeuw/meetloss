

corals <-
  function(theData,
           ftype = TRUE,
           xknots = NULL,
           yknots = NULL,
           xdegree = -1,
           ydegree = -1,
           xordinal = FALSE,
           yordinal = FALSE,
           xties = "s",
           yties = "s",
           xmissing = "m",
           ymissing = "m",
           xname = "X",
           yname = "Y",
           ndim = 2,
           itmax = 1000,
           eps = 1e-6,
           seed = 123,
           verbose = FALSE) {
    if (ftype) {
      xy <- preCorals(as.matrix(theData))
      x <- xy[, 1, drop = FALSE]
      y <- xy[, 2, drop = FALSE]
    } else {
      x <- theData[, 1, drop = FALSE]
      y <- theData[, 2, drop = FALSE]
    }
    if (is.null(xknots))
      xknots <- knotsD(x)
    if (is.null(yknots))
      yknots <- knotsD(y)
    g <- makeGifi(
      theData = cbind(x, y),
      knots = c(xknots, yknots),
      degrees = c(xdegree, ydegree),
      ordinal = c(xordinal, yordinal),
      ties = c(xties, yties),
      copies = c(ndim, ndim),
      missing = c(xmissing, ymissing),
      active = c(TRUE, TRUE),
      names = c(xname, yname),
      sets = c(1, 2)
    )
    h <- gifiEngine(
      gifi = g,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    xg <- h$xGifi[[1]][[1]]
    yg <- h$xGifi[[2]][[1]]
    return(structure(
      list(
        burt = crossprod(cbind(g[[1]][[1]]$basis, g[[2]][[1]]$basis)),
        objectscores = h$x,
        xtransform = postCorals(x, xg$transform),
        ytransform = postCorals(y, yg$transform),
        rhat = cor(cbind(xg$transform, yg$transform)),
        xweights = xg$weights,
        yweights = yg$weights,
        xscores = postCorals(x, xg$scores),
        yscores = postCorals(y, yg$scores),
        xdmeasure = crossprod(xg$scores),
        ydmeasure = crossprod(yg$scores),
        xquantifications = xg$quantifications,
        yquantifications = yg$quantifications,
        xloadings = crossprod(xg$transform, h$x),
        yloadings = crossprod(yg$transform, h$x),
        lambda = (crossprod(xg$scores) + crossprod(yg$scores)) / 2,
        ntel = h$ntel,
        f = h$f
      ),
      class = "corals"
    ))
  }
