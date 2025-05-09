data(small, package = "homals")
small <- makeNumeric(small)
small <- cbind(NA, small)
small_knots <- c(as.list(NULL), knotsD(small))
small_degrees <- c(-2, -1, -1, -1)
small_ordinal <- c(FALSE, FALSE, FALSE, FALSE)
small_ties <- c("s", "s", "s", "s")
small_copies <- c(2, 2, 2, 2)
small_missing <- c("m", "m", "m", "m")
small_active <- c(TRUE, TRUE, TRUE, TRUE)
small_names <- c("x", "a", "b", "c")
small_sets <- 1:4

h <- makeGifi (
  small,
  small_weights,
  small_knots,
  small_degrees,
  small_ordinal,
  small_ties,
  small_copies,
  small_missing,
  small_active,
  small_names,
  small_sets
)