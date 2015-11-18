#' Operate on a file which may or may not be an archive
#'
#' @param file.name name of the file to operate on
#' @param fun function to apply to the opened file
process.archive <- function (file.name, fun, ...)
{
    if (grepl("\\.bz2?$", file.name, perl=TRUE)) {
        f <- bzfile(file.name, open="r")
    } else if (grepl("\\.gz$", file.name)) {
        f <- gzfile(file.name, open="r")
    } else if (grepl("\\.xz$", file.name)) {
        f <- xzfile(file.name, open="r")
    } else {
        f <- file(file.name, open="r")
    }
    result <- fun(f, ...)
    close(f)
    result
}
