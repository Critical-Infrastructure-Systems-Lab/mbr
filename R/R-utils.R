#' Scale columns of a matrix
#'
#' Same as `base::scale()` but much faster.
#' @param x A matrix.
#' @param add_attr If TRUE, the column means and standard deviations are returned as attributes. This is consistent with [base::scale()].
#' @return The scaled matrix.
#' @section Reference:
#' This function was adopted from John Muschelli's code on \href{hopstat.wordpress.com/2016/02/23/a-faster-scale-function/}{StackOverflow}, but I changed the underlying functions to calculate mean and standard deviation from `matrixStats` to `Rfast`, which is much faster.
#' @export
colScale <- function(x, add_attr = TRUE) {

  # Get the column means
  # cm <- colmeans(x)
  # csd <- colVars(x, std = TRUE)

  cm <- colMeans(x)
  csd <- colSds(x)
  x <- t((t(x) - cm) / csd)
  if (add_attr) {
    attr(x, "scaled:center") <- cm
    attr(x, "scaled:scale") <- csd
  }
  x
}




#' Scale rows of a Matrix
#'
#' Similar to [colScale]
#' @return The scaled matrix.
#' @export
rowScale <- function(x, add_attr = TRUE) {

  # Get the column means
  # rm <- rowmeans(x)
  # rsd <- rowVars(x, std = TRUE)
  rm <- rowMeans(x)
  rsd <- rowSds(x, std = TRUE)
  x <- (x - rm) / rsd
  if (add_attr) {
    attr(x, "scaled:center") <- rm
    attr(x, "scaled:scale") <- rsd
  }
  x
}

#' Unscale columns of a matrix
#'
#' Backtransform a matrix that was scaled before.
#' @param x A matrix.
#' @param cm A vector of column means
#' @param csd A vector of column standard deviations
#' @param drop If TRUE, drop the attributes
#' @return The unscaled matrix
#' @export
colUnscale <- function(x, cm, csd, drop = TRUE) {
  x <- sweep(x, 2, csd, '*')
  x <- sweep(x, 2, cm, '+')
  if (drop) attributes(x)[2:3] <- NULL
  x
}

#' Unscale rows of a matrix
#'
#' Backtransform a matrix that was scaled before.
#' @param x A matrix.
#' @param rm A vector of row means
#' @param rsd A vector of row standard deviations
#' @param drop If TRUE, drop the attributes
#' @return The unscaled matrix
#' @export
rowUnscale <- function(x, rm, rsd, drop = TRUE) {
  x <- sweep(x, 1, rsd, '*')
  x <- sweep(x, 1, rm, '+')
  if (drop) attributes(x)[2:3] <- NULL
  x
}

