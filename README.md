# Mass Balance Regression

This package implements the Mass-balance-adjusted regression algorithm for sub-annual streamflow reconstruction. The algorithm implements a penalty term to minimize the differences between the total seasonal flow reconstruction and the annual flow reconstruction. Details are presented in Nguyen et al (2020).

Currently the package is only available on GitHub, but it will be available on CRAN soon.

To install the development from GitHub

``` {.r}
install.packages('remotes')
remotes::install_github('ntthung/mbr')
```

The package has two main functions, `mb_reconstruction()` for reconstruction, and `cv_mb()` for cross-validation.

Example reconstruction

``` {.r}
fit <- mb_reconstruction(
  instQ = p1Seasonal,
  pc.list = pc3seasons,
  start.year = 1750,
  lambda = 1,
  log.trans = 1:3
)
```

Example cross-validation

``` {.r}
# Create hold-out chunks
set.seed(24)
cvFolds <- make_Z(
  obs = 1922:2003,
  nRuns = 50, 
  frac = 0.25,
  contiguous = TRUE
)
# Run cross validation
cv <- cv_mb(
  instQ = p1Seasonal,
  pc.list = pc3seasons,
  cv.folds = cvFolds,
  start.year = 1750,
  lambda = 1,
  log.trans = 1:3,
  return.type = 'metric means'
)
# Round up to two decimal places
cv[, (2:6) := lapply(.SD, round, digits = 2), .SDcols = 2:6][]
```

Type `browseVignettes('mbr')` for details.

**References**

Nguyen, H. T. T., Galelli, S., Xu, C., & Buckley, B. (2020). Multi-Proxy, Multi-Season Streamflow Reconstruction with Mass Balance Adjustment. Earth and Space Science Open Archive, 22. <https://doi.org/10.1002/essoar.10504791.1>
