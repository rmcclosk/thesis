collect.metadata.mm <- function (mm.file)
{
    as.data.frame(yaml.load(substring(process.archive(mm.file, readLines, n=2)[2], 2)))
}

collect.metadata.newick <- function (nwk.file)
{
    as.data.frame(yaml.load(substring(process.archive(nwk.file, readLines, n=1)[1], 2)))
}

collect.metadata <- function (data.files)
{
    data.files <- sort(Sys.glob(data.files))
    names <- basename(data.files)
    names <- sapply(names, function (x) sub("\\.gz", "", x))
    names <- sapply(names, function (x) sub("\\.tar", "", x))
    names <- sapply(names, function (x) sub("\\.bz2?", "", x, perl=TRUE))
    names <- sapply(names, function (x) sub("\\.xz", "", x, perl=TRUE))

    extensions <- lapply(strsplit(names, ".", fixed=TRUE), tail, 1)
    metadata <- do.call(rbind, mapply(function (f, ext) {
        switch(ext, mtx=collect.metadata.mm(f),
                    nwk=collect.metadata.newick(f))
    }, data.files, extensions, SIMPLIFY=FALSE))

    wcon <- textConnection("buf", open="w")
    write.table(metadata, wcon)
    close(wcon)

    rcon <- textConnection(buf, open="r")
    metadata <- read.table(rcon)
    close(rcon)

    metadata
}

# Collect data from several files whose names contain metadata in the format
# key1-value1_key2-value2.foo
collect.data <- function (data.files)
{
    names <- basename(data.files)

    names <- lapply(names, function (x) sub("\\.gz", "", x))
    names <- lapply(names, function (x) sub("\\.tar", "", x))
    names <- lapply(names, function (x) sub("\\.bz2?", "", x, perl=TRUE))
    names <- lapply(names, function (x) sub("\\.xz", "", x, perl=TRUE))
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

    data <- lapply(data.files, process.archive, read.delim)
    data <- suppressWarnings(mapply(cbind, data, by(metadata, 1:nrow(metadata), identity), SIMPLIFY=FALSE))
    do.call(rbind, data)
}
