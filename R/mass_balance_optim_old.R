#' Least square with mass balance penalty
#'
#' @param hat A vector of estimated Dry, Wet, and Annual log-transformed flow
#' @param obs A vector of observed Dry, Wet, and Annual log-transformed flow
#' @param lambda Penalty weight
#' @param Jtype Objective function type
#' @return Objective function value: least squares plus a penalty term
#' @export
lsq_mb_old <- function(hat, obs, lambdas, Jtype = 2) {

  s1 <- sum((hat - obs)^2)

  N <- length(obs) / 3
  rowStt <- c(1, N + 1, N * 2 + 1)
  rowEnd <- c(N, N * 2, N * 3)

  Dhat <- hat[rowStt[1]:rowEnd[1]]
  What <- hat[rowStt[2]:rowEnd[2]]
  Qhat <- hat[rowStt[3]:rowEnd[3]]

  D <- obs[rowStt[1]:rowEnd[1]]
  W <- obs[rowStt[2]:rowEnd[2]]
  Q <- obs[rowStt[3]:rowEnd[3]]

  if (length(lambdas) == 1 && lambdas == 0) {
    s2 <- 0
  } else {
    expD <- exp(Dhat)
    expW <- exp(What)
    expQ <- exp(Qhat)
    if (any(is.infinite(expD)) || any(is.infinite(expW)) || any(is.infinite(expQ))) {
      s2 <- rep(Inf, length(lambdas))
    } else {
      if (Jtype == 1) { # Three reconstructions aware of one another's target
        deltaD <- if (any(expQ < expW)) Inf else log(expQ - expW) - D
        deltaW <- if (any(expQ < expD)) Inf else log(expQ - expD) - W
        deltaQ <- log(expD + expW) - Q
        s2 <- sum(deltaD^2) + sum(deltaW^2) + sum(deltaQ^2)
      } else { # Three reconstructions aware of one another
        delta <- log(expD + expW) - Qhat
        s2 <- sum(delta^2)
      }
    }
  }
  fv <- s1 + lambdas * s2
  replace(fv, is.infinite(fv), 1e7)
}

#' Objective function from parameter
#'
#' This is a wrapper for `lsq_mb()`. It first calculates `hat`, then calls `lsq_mb()`.
#' This is used in `optim()`, so it returns a scalar.
#' @param beta Parameters
#' @param X Inputs, must have columns of 1 added
#' @param Y Observed Dry, Wet, and Annual log-transformed flows
#' @param lambda Penalty weight
#' @return Objective function value
#' @export
obj_fun_old <- function(beta, X, Y, lambda, Jtype = 2) {

  hat <- X %*% beta
  lsq_mb_old(hat, Y, lambda, Jtype = Jtype)
}

#' Fit parameters with mass balance criterion
#'
#' @param X Matrix of inputs, constructed from U, V, and Z
#' @param Y Vector of observed log-transformed flow (following Dry, Wet, and Ann order)
#' @param lambdas Vector of weights to search
#' @export
mb_fit_old <- function(X, Y, lambdas = 1, Jtype = 2) {

  # Solve the free optimization
  XTX <- crossprod(X)
  XTY <- crossprod(X, Y)
  betaFree <- solve(XTX, XTY)

  # Solve the constrained optimization (with constraints changed to penalties)
  betaCons <- sapply(lambdas, function(lambda)
    stats::optim(betaFree, obj_fun_old, method = 'L-BFGS-B', X = X, Y = Y, lambda = lambda, Jtype = Jtype)$par)
  colnames(betaCons) <- paste0('l=', lambdas)
  # rownames(betaCons) <- rownames(betaFree)
  betaCons
}

#' Back-transformation
#'
#' Transform the reconstructed values back to the flow space
#' and convert to data.table
backtrans_old <- function(hat, years) {
  N <- dim(hat)[1] / 3
  recCons <- data.table(exp(hat), check.names = FALSE)
  recCons[, season := rep(c('NJ', 'JO', 'WY'), each = N)
        ][, year := rep(years, 3)]
  melt(recCons, id.var = c('season', 'year'), variable.name = 'lambda', value.name = 'Q', variable.factor = FALSE)
}
#' Cross-validation
#'
#' @inheritParams mb_fit
#' @param instQ Instrumental data for all seasons
#' @param z Vector of the indices of the hold-out chunk
#' @param return.type If 'mb', only the objective function value is returned. If 'metrics', all metrics are returned. If 'Q', all Q predictions are returned.
#' @export
cv_mb_old <- function(instQ, pc.dry, pc.wet, pc.ann, cv.folds, start.year,
                  lambdas = 1, Jtype = 2, return.type = c('mb', 'metrics', 'Q')) {

  # Setup
  end.year <- max(instQ$year)
  years <- start.year:end.year
  instInd <- which(years %in% instQ$year)

  U <- prepend_ones(pc.dry[instInd, , drop = FALSE])
  V <- prepend_ones(pc.wet[instInd, , drop = FALSE])
  Z <- prepend_ones(pc.ann[instInd, , drop = FALSE])

  X <- as.matrix(Matrix::bdiag(U, V, Z))
  Y <- instQ[c('NJ', 'JO', 'WY'), log(Qa)]

  N <- length(instInd)
  dryInd <- 1:N
  wetInd <- (N + 1):(N * 2)
  annInd <- (N * 2 + 1):(N * 3)

  # Cross-validation routine
  one_cv <- function(z) {
    calInd <- c(dryInd[-z], wetInd[-z], annInd[-z])
    valInd <- c(dryInd[ z], wetInd[ z], annInd[ z])

    # Calibration
    beta <- mb_fit_old(X[calInd, ], Y[calInd], lambdas, Jtype = Jtype)

    # Validation
    hat <- X %*% beta # All instrumental period prediction
    fval <- lsq_mb_old(hat[valInd], Y[valInd], lambdas, Jtype = Jtype)
    if (return.type == 'mb') fval else {
      Qcv <- merge(backtrans_old(hat, years[instInd]), instQ, by = c('year', 'season'))
      if (return.type == 'Q') Qcv else {
        metrics <- Qcv[, as.data.table(t(calculate_metrics(Q, Qa, z))), by = .(season, lambda)]
        metrics[, fval := fval, by = season]
        metrics[]
      }
    }
  }

  if (return.type == 'mb') { # A vector of fval
    sapply(cv.folds, one_cv)
  } else { # A data.table of all metrics
    lapplyrbind(cv.folds, one_cv, id = 'rep')
  }
}

#' Mass-balance-adjusted reconstruction
#'
#' @param instQ Instrumental data, in the same order as pc.list. The "season" column must be a factor.
#' @param pc.dry,pc.wet,pc.ann List of PC matrices. The first element is for the first season, second element for second season, and so on. The last element is for the annual reconstruction.
#' @param start.year The first year of record
#' @param lambda The penalty weight
#' @export
mb_reconstruction_old <- function(instQ, pc.dry, pc.wet, pc.ann, start.year, lambdas = 1, Jtype = 2) {

  # Setup
  end.year <- max(instQ$year)
  years <- start.year:end.year
  instInd <- which(years %in% instQ$year)

  U      <- prepend_ones(pc.dry)
  V      <- prepend_ones(pc.wet)
  Z      <- prepend_ones(pc.ann)
  X      <- as.matrix(Matrix::bdiag(U, V, Z))
  Xtrain <- as.matrix(Matrix::bdiag(U[instInd, , drop = FALSE],
                                    V[instInd, , drop = FALSE],
                                    Z[instInd, , drop = FALSE]))
  Y <- instQ[c('NJ', 'JO', 'WY'), log(Qa)]

  # Calibration
  beta <- mb_fit_old(Xtrain, Y, lambdas, Jtype = Jtype)
  # Prediction
  hat <- X %*% beta # Convert to vector
  # Tidy up
  backtrans_old(hat, years)
}

#' Calculate the objective function value in cross-validation
#'
#' Takes the negative value for use in GA as GA maximizes. So the less negative the better.
#' @param choice Binary vector of choices
#' @param pool Data.table listing the sites with significant correlations for each season
#' @param n.site.min Minimum number of sites to be kept for each season
#' @param Xpool The input matrix of all the sites in the pool
#' @export
cv_site_selection_old <- function(choice, pool, n.site.min = 5,
                              seasons,
                              Xpool, instQ, cv.folds, start.year,
                              lambda = 1, Jtype = 2,
                              use.robust.mean = FALSE) {

  poolSubset <- pool[choice == 1L]
  # Constraint: at least n.site.min for each season
  # If not met, return -1e7
  n.sites <- poolSubset[, .N, by = season][, N]
  if (length(n.sites) == 3 && all(n.sites > n.site.min)) {
    end.year <- max(instQ$year)
    years <- start.year:end.year
    instInd <- which(years %in% instQ$year)

    pc3seasons <- lapply(seasons, function(s) {
      Qa <- instQ[s]
      X <- Xpool[, poolSubset[s, site]]
      PC <- wPCA(X, use.eigen = FALSE, return.matrix = TRUE)
      sv <- input_selection(PC[instInd, , drop = FALSE],
                            Qa$Qa, 'leaps backward', nvmax = 8)
      PC[, sv, drop = FALSE]
    })
    pcDry <- pc3seasons$NJ
    pcWet <- pc3seasons$JO
    pcAnn <- pc3seasons$WY

    cvFval <- cv_mb_old(instQ, pcDry, pcWet, pcAnn, cv.folds, start.year,
                    lambdas = lambda, Jtype = Jtype, return.type = 'mb')

    if (use.robust.mean) -dplR::mean(cvFval) else -mean(cvFval)
  } else -1e7
}

# memo_cv_site_selection_old <- memoise::memoise(cv_site_selection_old)
