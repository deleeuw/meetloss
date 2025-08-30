source("gifiStructures.R")
source("gifiEngine.R")
source("gifiUtilities.R")
source("matrix.R")
source("coneRegression.R")
source("splineBasis.R")

homals <-
  function(theData,
           knots = knotsD(theData),
           degrees = -1,
           ordinal = FALSE,
           ndim = 2,
           ties = "s",
           missing = "m",
           names = colnames(theData, do.NULL = FALSE),
           active = TRUE,
           itmax = 1000,
           eps = 1e-6,
           seed = 123,
           verbose = FALSE)  {
    nvars <- ncol(theData)
    g <- makeGifi(
      theData = theData,
      weights = NULL,
      knots = knots,
      degrees = reshape(degrees, nvars),
      ordinal = reshape(ordinal, nvars),
      ties = reshape(ties, nvars),
      copies = rep(ndim, ncol(theData)),
      missing = reshape(missing, nvars),
      active = reshape(active, nvars),
      names = names,
      sets = 1:nvars
    )
    h <- gifiEngine(
      gifi = g,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    a <- v <- z <- d <- y <- o <- as.list(1:ncol(theData))
    dsum <- matrix(0, ndim, ndim)
    nact <- 0
    for (j in 1:nvars) {
      jgifi <- h$xGifi[[j]][[1]]
      v[[j]] <- jgifi$transform
      a[[j]] <- jgifi$weights
      y[[j]] <- jgifi$scores
      z[[j]] <- jgifi$quantifications
      cy <- crossprod(y[[j]])
      if (g[[j]][[1]]$active) {
        dsum <- dsum + cy
        nact <- nact + 1
      }
      d[[j]] <- cy
      o[[j]] <- crossprod(h$x, v[[j]])
    }
    return(structure(
      list(
        transform = v,
        rhat = corList(v),
        objectscores = h$x,
        scores = y,
        quantifications = z,
        dmeasures = d,
        lambda = dsum / nact,
        weights = a,
        loadings = o,
        ntel = h$ntel,
        f = h$f
      ),
      class = "homals"
    ))
  }
