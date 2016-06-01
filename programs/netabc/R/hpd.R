#' Calculated a weighted highest density interval.
#'
#' @param x sample from distribution
#' @param wt weights of each element in the sample
#' @param conf proportion of probability density contained in interval
#' @return a two-element vector containing upper and lower HPD interval bounds
#' @export
wtd.hpd <- function (x, wt, conf=0.95) {
    cdf <- setDT(wtd.Ecdf(x, wt))
    setkey(cdf, x)
    cdf <- cdf[,tail(.SD, 1), by=x]

    if (is.integer(x)) {
        icdf <- function (p)
            sapply(p, function (pi) cdf[ecdf >= pi, head(x, 1)])
        find.hpd.discrete(icdf, x[which.max(wt)], conf=conf)
    } else {
        icdf <- cdf[,approxfun(ecdf, x)]
        hpd(icdf, conf=conf)
    }
}

#' Find a highest density interval for a discrete-valued distribution.
#'
#' Helper function for wtd.hpd().
#'
#' @param icdf data.table containing inverse CDF values
#' @param conf how much mass the HPD interval should contain
#' @param tol narrowness of steps to search ICDF
find.hpd.discrete <- function (icdf, mode, conf=0.95, tol=1e-3) {
    # interval widths for [0, conf], [tol, conf+tol], ...
    p <- seq(0, 1-conf, tol)
    w <- na.omit(sapply(p, function (pi) icdf(pi + conf) - icdf(pi)))

    # intervals with minimum width 
    mins <- which(w == min(w))
    all.ints <- cbind(icdf(p[mins]), icdf(p[mins]+conf))

    # interval which is closest to being symmetric about the mode    
    keep.int <- which.min(abs((mode - all.ints[,1]) - (all.ints[,2] - mode)))
    c(lower=all.ints[keep.int,1], upper=all.ints[keep.int,2])
}
