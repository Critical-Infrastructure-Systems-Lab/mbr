#' Least square with mass balance penalty
#'
#' @param hat A vector of estimated flow in the transformed space.
#' @param obs A vector of observed flow in the transformed space.
#' @param lambda Penalty weight.
#' @param log.seasons A vector containing the indices of the seasons that are log-transformed.
#' @param log.ann TRUE if the annual reconstruction is log-transformed.
#' @param N The number of targets (number of seasons plus one for the annual reconstruction).
#' @param sInd Indices of the seasons, i.e, 1...N-1
#' @return Objective function value: least squares plus a penalty term.
#' @export
lsq_mb <- function(hat, obs, lambda, log.seasons, log.ann, N, sInd) {

  s1 <- sum((hat - obs)^2) # Regression part

  if (lambda == 0) {       # Penalty part
    s2 <- 0
  } else {

    hatBack <- matrix(hat, ncol = N)
    hatBack[, log.seasons] <- exp(hatBack[, log.seasons])

    if (any(is.infinite(hatBack))) {
      s2 <- 1e7 # GA needs finite f value
    } else {
      # Take sum and transform if necessary
      totalSeasonal <- if (log.ann) log(rowSums(hatBack[, sInd])) else rowSums(hatBacks[, sInd])
      s2 <- sum((totalSeasonal - hatBack[, N])^2)
    }
  }
  s1 + lambda * s2
}

#' Objective function from parameters
#'
#' This is a wrapper for `lsq_mb()`. It first calculates `hat`, then calls `lsq_mb()`.
#' This is used in `optim()`, so it returns a scalar.
#' @param beta Parameters
#' @param X Inputs, must have columns of 1 added
#' @param Y Observed Dry, Wet, and Annual log-transformed flows
#' @inheritParams lsq_mb
#' @return Objective function value
#' @export
obj_fun <- function(beta, X, Y, lambda, log.seasons, log.ann, N, sInd) {

  hat <- X %*% beta
  lsq_mb(hat, Y, lambda, log.seasons, log.ann, N, sInd)
}

#' Fit parameters with mass balance criterion
#'
#' @inheritParams obj_fun
#' @return A one-column matrix of beta value
#' @export
mb_fit <- function(X, Y, lambda, log.seasons, log.ann, N, sInd) {

  # Solve the free optimization and use the result as initial value L-BFGS-B search
  # This will speed up the search process
  XTX <- crossprod(X)
  XTY <- crossprod(X, Y)
  betaFree <- solve(XTX, XTY)

  # Solve the constrained optimization (with constraints changed to penalties)
  stats::optim(
    betaFree, obj_fun, method = 'L-BFGS-B',
    X = X, Y = Y, lambda = lambda,
    log.seasons = log.seasons, log.ann = log.ann, N = N, sInd = sInd)$par
}

#' Prepend a column of ones
#'
#' @param x The input matrix
#' @return x with a column of ones prepended, which is named 'Int' for 'intercept'
#' @export
prepend_ones <- function(x) cbind('Int' = rep(1, dim(x)[1]), x)

#' Back-transformation
#'
#' Transform the reconstructed values back to the flow space
#' and convert to data.table
#' @param years A vector of all years in the study period
#' @inheritParams lsq_mb
#' @param season.names A character vector containing the names of the seasons
#' @export
back_trans <- function(hat, years, log.trans, num.targets, season.names) {

  hatBack <- matrix(hat, ncol = num.targets)
  hatBack[, log.trans] <- exp(hatBack[, log.trans])

  data.table(
    Q = c(hatBack),
    season = rep(season.names, each = length(hat) / num.targets),
    year = rep(years, num.targets))
}

#' Mass-balance-adjusted reconstruction
#'
#' @param instQ Instrumental data, in the same order as pc.list. The "season" column must be a factor.
#' @param pc.list List of PC matrices. The first element is for the first season, second element for second season, and so on. The last element is for the annual reconstruction.
#' @param start.year The first year of record
#' @param lambda The penalty weight
#' @export
mb_reconstruction <- function(instQ, pc.list, start.year, lambda = 1, log.trans, num.targets) {

  # Setup
  years     <- start.year:max(instQ$year)
  instInd   <- which(years %in% instQ$year)
  XList     <- lapply(pc.list, prepend_ones)
  XListInst <- lapply(XList, function(x) x[instInd, , drop = FALSE])
  X         <- as.matrix(Matrix::bdiag(XList))
  XTrain    <- as.matrix(Matrix::bdiag(XListInst))
  instQList <- split(instQ, by = 'season')

  Y <- matrix(instQ$Qa, ncol = num.targets)
  Y[, log.trans] <- log(Y[, log.trans])
  Y <- c(Y)

  # Calibration
  log.seasons <- which(log.trans < num.targets)
  log.ann <- max(log.trans) == num.targets
  sInd <- seq_len(num.targets - 1)
  beta <- mb_fit(XTrain, Y, lambda, log.seasons, log.ann, num.targets, sInd)
  # Prediction
  hat <- X %*% beta # Convert to vector
  # Tidy up
  DT <- back_trans(hat, years, log.trans, num.targets, levels(instQ$season))
  DT[, lambda := lambda][]
}


#' Cross-validation
#'
#' @inheritParams mb_reconstruction
#' @param pc.list.inst List of PC matrices
#' @param cv.folds A list containing the cross validation folds
#' @param return.type If 'mb', only the objective function value is returned. If 'metrics', all metrics are returned. If 'Q', all Q predictions are returned.
#' @export
cv_mb <- function(instQ, pc.list, cv.folds, start.year,
                  lambda = 1,
                  log.trans, num.targets,
                  return.type = c('mb', 'metrics', 'Q')) {

  # Setup
  years     <- start.year:max(instQ$year)
  instInd   <- which(years %in% instQ$year)
  XListInst <- lapply(pc.list, function(x) prepend_ones(x[instInd, , drop = FALSE]))
  XTrain    <- as.matrix(Matrix::bdiag(XListInst))
  instQList <- split(instQ, by = 'season')

  Y <- matrix(instQ$Qa, ncol = num.targets)
  Y[, log.trans] <- log(Y[, log.trans])
  Y <- c(Y)
  L <- length(Y)         # Total number of data points

  indMat <- matrix(seq_along(Y), ncol = num.targets)

  log.seasons <- which(log.trans < num.targets)
  log.ann <- max(log.trans) == num.targets
  sInd <- seq_len(num.targets - 1)

  # Cross-validation routine
  one_cv <- function(z) {
    calInd <- c(indMat[-z, ])
    valInd <- c(indMat[ z, ])

    # Calibration

    beta <- mb_fit(XTrain[calInd, ], Y[calInd], lambda, log.seasons, log.ann, num.targets, sInd)

    # Validation
    hat <- XTrain %*% beta # All instrumental period prediction
    fval <- lsq_mb(hat[valInd], Y[valInd], lambda, log.seasons, log.ann, num.targets, sInd)
    if (return.type == 'mb') {
      ans <- fval
    } else {
      Qcv <- merge(
        back_trans(hat, years[instInd], back.funs, num.targets, levels(instQ$season)),
        instQ, by = c('year', 'season'))
      setnames(Qcv, 'V1', 'Q')
      if (return.type == 'Q') {
        ans <- Qcv
      } else {
        metrics <- Qcv[, as.data.table(t(calculate_metrics(Q, Qa, z))), by = season]
        metrics[, fval := fval, by = season]
        ans <- metrics
      }
    }
    ans
  }

  if (return.type == 'mb') { # A vector of fval
    unlist(lapply(cv.folds, one_cv), use.names = FALSE)
  } else { # A data.table of all metrics or all reps
    lapplyrbind(cv.folds, one_cv, id = 'rep')
  }
}

#' Calculate the objective function value in cross-validation
#'
#' Takes the negative value for use in GA as GA maximizes. So the less negative the better.
#' @param choice Binary vector of choices
#' @param pool Data.table listing the sites with significant correlations for each season
#' @param n.site.min Minimum number of sites to be kept for each season
#' @param Xpool The input matrix of all the sites in the pool
#' @export
cv_site_selection <- function(choice, pool, n.site.min = 5,
                              seasons,
                              Xpool, instQ, cv.folds, start.year,
                              lambda = 1,
                              log.trans, num.targets,
                              use.robust.mean = FALSE) {

  poolSubset <- pool[choice == 1L]

  # Constraint: at least n.site.min for each season
  # If not met, return -1e7
  n.sites <- poolSubset[, .N, by = season][, N]

  if (length(n.sites) == num.targets && all(n.sites > n.site.min)) {

    years     <- start.year:max(instQ$year)
    instInd   <- which(years %in% instQ$year)

    pcList <- lapply(seasons, function(s) {
      Qa <- instQ[s]
      X <- Xpool[, poolSubset[s, site]]
      PC <- wPCA(X, use.eigen = FALSE, return.matrix = TRUE)
      sv <- input_selection(
        PC[instInd, , drop = FALSE],
        Qa$Qa, 'leaps backward', nvmax = 8)
      PC[, sv, drop = FALSE]
    })

    cvFval <- cv_mb(instQ, pcList, cv.folds, start.year, lambda, log.trans, num.targets, return.type = 'mb')

    if (use.robust.mean) -dplR::mean(cvFval) else -mean(cvFval)
  } else -1e7
}
