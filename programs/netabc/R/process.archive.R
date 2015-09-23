process.archive <- function (file.name, fun, ...)
{
    if (grepl("\\.bz2?$", file.name, perl=TRUE)) {
        f <- bzfile(file.name)
    } else if (grepl("\\.gz$", file.name)) {
        f <- gzfile(file.name)
    } else if (grepl("\\.xz$", file.name)) {
        f <- xzfile(file.name)
    } else {
        f <- file(file.name)
    }
    result <- fun(f, ...)
    close(f)
    result
}
