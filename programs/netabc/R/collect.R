# Collect data from several files whose names contain metadata in the format
# key1-value1_key2-value2.foo
collect.data <- function (data.files)
{
    names <- basename(data.files) 
    names <- lapply(names, function (x) sub("(.*)\\..*?$", "\\1", x, perl=TRUE))
    names <- sapply(names, strsplit, "_")

    names <- lapply(names, strsplit, "-")
    variables <- lapply(names, sapply, "[[", 1)
    values <- lapply(names, sapply, "[[", 2)
    metadata <- mapply(setNames, values, variables, SIMPLIFY=FALSE)
    metadata <- as.data.frame(do.call(rbind, metadata), stringsAsFactors=FALSE)

    wcon <- textConnection("buf", open="w")
    write.table(metadata, wcon)
    close(wcon)

    rcon <- textConnection(buf, open="r")
    metadata <- read.table(rcon)
    close(rcon)

    data <- lapply(data.files, read.delim)
    data <- suppressWarnings(mapply(cbind, data, by(metadata, 1:nrow(metadata), identity), SIMPLIFY=FALSE))
    do.call(rbind, data)
}
