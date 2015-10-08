read.mm.symmetric <- function (file)
{
    dim <- scan(file, what=integer(0), n=2, comment.char="%", quiet=TRUE)[1]
    nelem <- dim * (dim + 1) / 2
    data <- scan(file, what=numeric(0), n=2+nelem, comment.char="%", quiet=TRUE)
    data <- data[3:(nelem+2)]

    m <- matrix(nrow=dim, ncol=dim)
    m[lower.tri(m, diag=TRUE)] <- data
    m <- t(m)
    m[lower.tri(m, diag=TRUE)] <- data
    m
}

read.mm.general <- function (file)
{
    dim <- scan(file, what=integer(0), n=2, comment.char="%", quiet=TRUE)
    nelem <- dim[1] * dim[2]
    data <- scan(file, what=numeric(0), n=2+nelem, comment.char="%", quiet=TRUE)
    data <- data[3:(nelem+2)]

    matrix(data, nrow=dim, ncol=dim)
}

# read a matrix in MatrixMarket format
read.mm <- function (file.name)
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

    mm.spec <- readLines(f, n=1)
    result <- NULL
    if (grepl("coordinate", mm.spec))
        result <- readMM(f)
    if (mm.spec == "%%MatrixMarket matrix array real symmetric")
    {
        result <- read.mm.symmetric(f)
    }
    else if (mm.spec == "%%MatrixMarket matrix array real general")
    {
        result <- read.mm.general(f)
    }

    close(f)
    if (is.null(result))
        stop(sprintf("MatrixMarket format '%s' not yet implemented", mm.spec))
    return(result)
}
