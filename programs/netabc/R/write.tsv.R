#' Write a TSV file
#'
#' This is just a thin wrapper around write.table with different default
#' options.
#' 
#' @param df data.frame to write
#' @param ... other arguments for write.table, in particular the file name
#' @export
write.tsv <- function (df, ...)
{
    write.table(df, ..., col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}
