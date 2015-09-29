collect.metadata.mm <- function (mm.file)
{
    as.data.frame(yaml.load(substring(process.archive(mm.file, readLines, n=2)[2], 2)))
}

collect.metadata.newick <- function (nwk.file)
{
    as.data.frame(yaml.load(substring(process.archive(nwk.file, readLines, n=1)[1], 2)))
}

collect.metadata.tsv <- function (tsv.file)
{
    as.data.frame(yaml.load(substring(process.archive(tsv.file, readLines, n=1)[1], 2)))
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
                    nwk=collect.metadata.newick(f),
                    tsv=collect.metadata.tsv(f))
    }, data.files, extensions, SIMPLIFY=FALSE))

    wcon <- textConnection("buf", open="w")
    write.table(metadata, wcon)
    close(wcon)

    rcon <- textConnection(buf, open="r")
    metadata <- read.table(rcon)
    close(rcon)

    rownames(metadata) <- data.files
    metadata
}

# Collect data and metadata from several files.
collect.data <- function (data.files)
{
    data.files <- sort(Sys.glob(data.files))
    metadata <- collect.metadata(data.files)
    data <- lapply(data.files, process.archive, read.table, header=TRUE)
    data <- suppressWarnings(mapply(cbind, data, 
                                    by(metadata, 1:nrow(metadata), identity),
                                    SIMPLIFY=FALSE))
    do.call(rbind, data)
}
