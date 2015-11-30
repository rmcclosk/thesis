#' Collect YAML-formatted metadata from a MatrixMarket file
#'
#' Metadata is stored in the comment line, which is the second line of the file.
#' It is preceeded by a '%' character.
#'
#' @param mm.file MatrixMarket formatted file
collect.metadata.mm <- function (mm.file)
{
    as.data.frame(yaml.load(substring(process.archive(mm.file, readLines, n=2)[2], 2)))
}

#' Collect YAML-formatted metadata from a Newick file
#'
#' Metadata is stored as a comment (preceded by '#') in the first line of the
#' file.
#'
#' @param nwk.file Newick formatted file
collect.metadata.newick <- function (nwk.file)
{
    as.data.frame(yaml.load(substring(process.archive(nwk.file, readLines, n=1)[1], 2)))
}

#' Collect YAML-formatted metadata from a TSV file
#'
#' Metadata is stored as a comment (preceded by '#') in the first line of the
#' file. This is also used as the default metadata collector for unknown
#' formats.
#'
#' @param tsv.file tab-separated value file
collect.metadata.tsv <- function (tsv.file)
{
    as.data.frame(yaml.load(substring(process.archive(tsv.file, readLines, n=1)[1], 2)))
}

#' Collect YAML-formatted metadata from a GML file
#'
#' Metadata is stored in the "comment" graph attribute. Since igraph doesn't
#' actually load the graph attributes, we have to grep for it in the first few
#' lines of the file.
#'
#' @param tsv.file tab-separated value file
collect.metadata.gml <- function (gml.file)
{
    getter <- function (f) { 
        strsplit(grep("comment", readLines(f, n=20), value=TRUE), '"')[[1]][2]
    }
    as.data.frame(yaml.load(process.archive(gml.file, getter)))
}

#' Collect YAML-formatted metadata from several files.
#'
#' This collates YAML metadata from a list of files or a globbing pattern into
#' a data.frame whose row names are the names of the files. If the files don't
#' all have the same metadata, missing columns are filled in with NA.
#'
#' @param data.files vector of file names, or glob pattern
#' @return A data.frame containing the metadata
#' @seealso \code{\link{collect.data}} to collect data and metadata together
#' @export
collect.metadata <- function (data.files)
{
    # realize the pattern, if there is one
    data.files <- sort(Sys.glob(data.files))

    # get the file names without directories or compression extensions
    names <- basename(data.files)
    names <- sapply(names, function (x) sub("\\.gz", "", x))
    names <- sapply(names, function (x) sub("\\.tar", "", x))
    names <- sapply(names, function (x) sub("\\.xz", "", x))
    names <- sapply(names, function (x) sub("\\.bz2?", "", x, perl=TRUE))

    # collect the metadata, using the files' extensions as a guide
    extensions <- lapply(strsplit(names, ".", fixed=TRUE), tail, 1)
    metadata <- mapply(function (f, ext) {
        switch(ext, mtx=collect.metadata.mm(f),
                    nwk=collect.metadata.newick(f),
                    tsv=collect.metadata.tsv(f),
                    gml=collect.metadata.gml(f),
                    collect.metadata.tsv(f))
    }, data.files, extensions, SIMPLIFY=FALSE)

    # union all the metadata fields
    all.cols <- Reduce(union, lapply(metadata, colnames))
    new.cols <- mapply(setdiff, list(all.cols), lapply(metadata, colnames), SIMPLIFY=FALSE)

    # fill in the fields each file doesn't have with NA 
    metadata <- mapply(function (df, new) {
        if (length(new) > 0) {
            cbind(df, as.data.frame(setNames(as.list(rep(NA, length(new))), new)))
        }
        else {
            df
        }
    }, metadata, new.cols, SIMPLIFY=FALSE)

    # put all together and set row names to files
    metadata <- do.call(rbind, metadata)
    rownames(metadata) <- data.files
    metadata
}

#' Collect data and metadata from several files.
#'
#' This function collects data and YAML-formatted metadata together. It's
#' assumed that the files are some kind of column-oriented thing (CSV or TSV).
#' They are read with read.table, and columns are added for each of the metadata
#' elements.
#' 
#' @param data.files a vector of file names, or a glob pattern
#' @param header passed to read.table, but the default is TRUE
#' @param ... additional arguments for read.table
#' @return a data.frame with collated data and metadata
#' @seealso \code{\link{collect.metadata}} to collect only metadata
#' @export
collect.data <- function (data.files, header=TRUE, ...)
{
    data.files <- sort(Sys.glob(data.files))
    metadata <- collect.metadata(data.files)
    data <- lapply(data.files, process.archive, read.table, header=header, ...)
    data <- suppressWarnings(mapply(cbind, data, 
                                    by(metadata, 1:nrow(metadata), identity),
                                    SIMPLIFY=FALSE))
    do.call(rbind, data)
}
