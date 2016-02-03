#' Merge a set of /sorted/ intervals
#'
#' @param intervals a data.table with columns "start" and "end"
#' @return a reduced data.table with non-overlapping intervals only
#' @export
intervals.merge <- function (dt) {
    if (nrow(dt) <= 1) return (dt)
    dt.copy <- copy(dt)
    i <- 2
    repeat {
        if (i > nrow(dt.copy)) break
        if (dt.copy[i-1,end+1] >= dt.copy[i,start]) {
            dt.copy[i-1,end := dt.copy[i,end]]
            dt.copy <- dt.copy[-i,]
        } else {
            i <- i + 1
        }
    }
    dt.copy
}
