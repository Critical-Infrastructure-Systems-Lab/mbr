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
  sigma <- sd(y)
  sigma_hat <- sd(yhat)
  r <- cor(yhat, y)
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
