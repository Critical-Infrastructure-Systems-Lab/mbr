#' Scale columns of a matrix
#'
#' Same as `base::scale()` but much faster.
#' @param x A matrix.
#' @param add_attr If TRUE, the column means and standard deviations are returned as attributes. This is consistent with [base::scale()].
#' @return The scaled matrix.
#' @section Reference:
#' This function was adopted from John Muschelli's code on \href{https://hopstat.wordpress.com/2016/02/23/a-faster-scale-function/}{StackOverflow}, but I changed the underlying functions to calculate mean and standard deviation from `matrixStats` to `Rfast`, which is much faster.
colScale <- function(x, add_attr = TRUE) {

  # Get the column means
  cm <- colmeans(x)
  csd <- colVars(x, std = TRUE)

  # cm <- colMeans(x)
  # csd <- colSds(x)
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
#' @inheritParams colScale
#' @return The scaled matrix.
rowScale <- function(x, add_attr = TRUE) {

  # Get the column means
  rm <- rowmeans(x)
  rsd <- rowVars(x, std = TRUE)
  # rm <- rowMeans(x)
  # rsd <- rowSds(x, std = TRUE)
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
#' @return The unscaled matrix
colUnscale <- function(x, cm, csd) {
  t(t(x) * csd + cm)
}

#' Unscale rows of a matrix
#'
#' Backtransform a matrix that was scaled before.
#' @param x A matrix.
#' @param rm A vector of row means
#' @param rsd A vector of row standard deviations
#' @return The unscaled matrix
rowUnscale <- function(x, rm, rsd) {
  x * rsd + rm
}

#' Reconstruction metrics
#'
#' Calculate reconstruction metrics from the instrumental period
#' @param sim A vector of reconstruction output for instrumental period
#' @param obs A vector of all observations
#' @param z A vector of left out indices in cross validation
#' @param norm.fun The function (unquoted name) used to calculate the normalizing constant. Default is `mean()`, but other functions such as `sd()` can also be used. THe function must take a vector as input and return a scalar as output, and must have an argument `na.rm = TRUE`.
#' @return A named vector of performance metrics
#' @examples
#' calculate_metrics(rnorm(100), rnorm(100), z = 1:10)
#' calculate_metrics(rnorm(100), rnorm(100), z = 1:10, norm.fun = sd)
#' @export
calculate_metrics <- function(sim, obs, z, norm.fun = mean) {
  train.obs <- obs[-z]
  obsInd <- which(!is.na(train.obs))
  train.obs <- train.obs[obsInd]
  train.sim <- sim[-z][obsInd]
  val.sim <- sim[z]
  val.obs <- obs[z]

  c(R2    = NSE(train.sim, train.obs), # Use the NSE form of R2
    RE    = RE(val.sim, val.obs, mean(train.obs)),
    CE    = NSE(val.sim, val.obs),
    nRMSE = nRMSE(val.sim, val.obs, norm.fun(obs, na.rm = TRUE)),
    KGE   = KGE(val.sim, val.obs)
  )
}

#' Make cross-validation folds.
#'
#' Make a list of cross-validation folds. Each element of the list is a vector of the cross-validation points for one cross-validation run.
#' @param obs Vector of observations.
#' @param nRuns Number of repetitions.
#' @param frac Fraction of left-out points. For leave-one-out, use `frac = 1`, otherwise use any value less than 1. Default is 0.1 (leave-10%-out).
#' @param contiguous Logical. If `TRUE`, the default, the left-out points are made in contiguous blocks; otherwise, they are scattered randomly.
#' @return A list of cross-validation folds
#' @examples
#' Z <- make_Z(p1Seasonal$Qa, nRuns = 30, frac = 0.25, contiguous = TRUE)
#' @export
make_Z <- function(obs, nRuns = 30, frac = 0.1, contiguous = TRUE) {
  obsInd <- which(!is.na(obs))
  if (frac == 1) {
    split(obsInd, obsInd)
  } else {
    n <- length(obsInd)
    k <- floor(n * frac) # leave-k-out
    if (contiguous) {
      maxInd <- n - k # Highest possible position in of obsInd
      if (maxInd < nRuns) { # Not enough samples, reduce k
        maxInd <- nRuns
        k <- n - nRuns
      }
      lapply(sort(sample(1:maxInd, nRuns)), function(x) obsInd[x:(x + k)])
    } else {
      replicate(nRuns, sort(sample(obsInd, k, replace = FALSE)), simplify = FALSE)
    }
  }
}

#' Nash-Sutcliffe Efficiency
#'
#' @param yhat Model outputs
#' @param y Observations
#' @return NSE
#' @examples
#' NSE(rnorm(100), rnorm(100))
#' @export
NSE <- function(yhat, y) {
  ybar <- mean(y)
  rss <- sum((y - yhat) * (y - yhat))
  tss <- sum((y - ybar) * (y - ybar))
  1 - rss/tss
}

#' Normalized root-mean-square error
#'
#' RMSE is normalized by the normalization constant
#' @param yhat Model outputs
#' @param y Observations
#' @param normConst The normalization constant
#' @return normalized RMSE
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' nRMSE(x, y, sd(y))
#' @export
nRMSE <- function(yhat, y, normConst) {
  rmse <- sqrt(mean((y - yhat) * (y - yhat)))
  rmse / normConst
}

#' Kling-Gupta Efficiency
#'
#' @param yhat Model outputs
#' @param y Observations
#' @return KGE
#' @examples
#' KGE(rnorm(100), rnorm(100))
#' @export
KGE <- function(yhat, y) {
  mu <- mean(y)
  mu_hat <- mean(yhat)
  sigma <- stats::sd(y)
  sigma_hat <- stats::sd(yhat)
  r <- stats::cor(yhat, y)
  alpha <- sigma_hat / sigma
  beta <- mu_hat / mu
  EDsq <- (r - 1) * (r - 1) + (alpha - 1) * (alpha - 1) + (beta - 1) * (beta - 1)
  1 - sqrt(EDsq)
}

#' Reduction of Error
#'
#' @param yhat Model outputs in the validation set
#' @param y Observations in the validation set
#' @param yc_bar Mean observations in the calibration set
#' @return RE
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' yc_bar <- mean(x[1:50])
#' RE(x[51:100], y[51:100], yc_bar)
#' @export
RE <- function(yhat, y, yc_bar) {
  rss <- sum((y - yhat) * (y - yhat))
  tss <- sum((y - yc_bar) * (y - yc_bar))
  1 - rss/tss
}
