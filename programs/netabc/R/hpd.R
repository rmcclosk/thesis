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
    cdf <- unique(cdf)
    icdf <- cdf[,approxfun(ecdf, x)]
    hpd(icdf)
}
