#' Calculated a weighted highest density interval.
#'
#' @param x sample from distribution
#' @param wt weights of each element in the sample
#' @param conf proportion of probability density contained in interval
#' @return a two-element vector containing upper and lower HPD interval bounds
#' @export
wtd.hpd <- function (x, wt, conf=0.95) {
    if (is.integer(x)) {
        message("Integer values passed, calling wtd.hpd.discrete")
        return(wtd.hpd.discrete(x, wt, conf))
    }
    cdf <- setDT(wtd.Ecdf(x, wt))
    setkey(cdf, x)
    cdf <- cdf[,head(.SD, 1), by=x]
    icdf <- cdf[,approxfun(ecdf, x)]
    hpd(icdf, conf=conf)
}

#' Find a highest density interval for a discrete-valued distribution.
#'
#' @param x sample from distribution
#' @param wt weights of each element in the sample
#' @param conf proportion of probability density contained in interval
#' @return a two-element vector containing upper and lower HPD interval bounds
#' @export
wtd.hpd.discrete <- function (x, wt, conf=0.95) {
    if (any(diff(x) != 1)) {
        stop("Not implemented yet for non-consecutive integers.")
    }
    
    # make an empirical CDF and find the mode
    cdf <- setDT(wtd.Ecdf(x, wt))
    mode <- cdf[,unique(x)][cdf[,which.max(diff(ecdf))]]
    
    # record the cdf excluding and including each value
    setkey(cdf, x)
    cdf <- cdf[,list(start=head(ecdf, 1), end=tail(ecdf, 1)), by=x]
    cdf[,start := c(head(start, 1), head(end, -1))]
    
    # find every possible interval (yeah it's a bit hacky, ideally we would
    # stream them or something but I don't want to bother right now because all
    # my use cases are small)
    intervals <- setnames(as.data.table(t(cdf[,combn(x, 2)])), c("start", "end"))
    intervals <- rbind(intervals, cdf[,list(start=x, end=x)])
    setkey(intervals, start, end)
    
    # record the width of each interval and the density they contain
    intervals[,width := end - start + 1]
    intervals[,density := cdf[end,end] - cdf[start,start]]
    
    # pick the minimum width interval containing at least conf of the density
    intervals <- intervals[density >= conf]
    intervals <- intervals[width == min(width)]
    
    # if there are multiple intervals left, pick the one which contains the
    # mode
    intervals <- intervals[mode >= start & mode <= end]
    
    # if there are still multiple left, pick the one with the highest density
    intervals <- intervals[density == min(density)]
    
    # finally, just pick one at random
    intervals <- intervals[sample(1:nrow(intervals), 1)]
    intervals[,c(lower=start, upper=end)]
}
