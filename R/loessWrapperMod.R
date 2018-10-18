#' wrapper script for loess modeling
#' @description adapted from bisoreg package
#' @param x predictor values
#' @param y response values
#' @param span.vals values of the tuning parameter to evaluate 
#' using cross validation
#' @param folds number of 'folds' for the cross-validation procedure
#' @examples loessWrapperMod (rnorm(200), rnorm(200), span.vals = 
#' seq(0.25, 1, by = 0.05),folds = 5)
#' @return LOESS model
#' @usage loessWrapperMod(x, y, span.vals = seq(0.25, 1, by = 0.05),
#' folds = 5)
#' @export
loessWrapperMod <- function(x, y, span.vals = seq(0.25, 1, by = 0.05),
folds = 5) {
    mae <- numeric(length(span.vals))
    theta.fit <- function(x, y, span) loess(y ~ x, span = span, 
    surface = "direct")
    theta.predict <- function(fit, x0) predict(fit, newdata = x0)
    ii = 0
    for (span in span.vals) {
        ii <- ii + 1
        y.cv <- bootstrap::crossval(x, y, theta.fit, theta.predict, 
        span = span, 
        ngroup = folds)$cv.fit
        fltr <- !is.na(y.cv)
        mae[ii] <- mean(abs(y[fltr] - y.cv[fltr]))
    }
    span <- span.vals[which.min(mae)]
    out <- loess(y ~ x, span = span, surface = "direct")
    return(out)
}