#' Get all qualifier data from a set of Genbank records.
#'
#' @param records Genbank records
#' @return a list of qualifier data.tables, one for each record
get.genbank.qualifiers <- function (records) {
    ft <- lapply(records, "[[", "GBSeq_feature-table")
    lapply(lapply(ft, "[[", "GBFeature_quals"), rbindlist)
}

#' Parse collection dates from Genbank.
#'
#' @param records Genbank records to parse
#' @return a vector of formatted collection dates with the time unit as a prefix
#' @export
parse.genbank.dates <- function (records) {
    quals <- get.genbank.qualifiers(records)
    dates <- sapply(quals, "[", i=GBQualifier_name == "collection_date", j=GBQualifier_value)
    dates <- unlist(ifelse(sapply(dates, length) == 0, NA, dates))

    days <- as.integer(as.POSIXct(strptime(dates, c("%d-%b-%Y")))) %/% 86400
    dates <- ifelse(is.na(days), dates, paste0("days_", days))

    months <- as.integer(as.POSIXct(strptime(paste0("01-", dates), c("%d-%b-%Y")))) %/% 2628000
    dates <- ifelse(is.na(months), dates, paste0("months_", months))

    years <- as.integer(as.POSIXct(strptime(paste0("01-Jan-", dates), c("%d-%b-%Y")))) %/% 31536000
    dates <- ifelse(is.na(years), dates, paste0("years_", years))
    ifelse(is.na(days) & is.na(months) & is.na(years), "NA_0", dates)
}

#' Parse countries from Genbank.
#'
#' @param records Genbank records to parse
#' @return a vector of countries
#' @export
parse.genbank.countries <- function (records) {
    quals <- get.genbank.qualifiers(records)
    country <- sapply(quals, "[", i=GBQualifier_name == "country", j=GBQualifier_value)
    unlist(ifelse(sapply(country, length) == 0, NA, country))
}

#' Parse genes from Genbank.
#'
#' @param records Genbank records to parse
#' @return a data.frame with 4 columns gbid, start, end, gene
#' @export
parse.genbank.genes <- function (records) {
    seqid <- lapply(records, "[[", "GBSeq_other-seqids")
    gbid <- mapply(grep, list("gi"), seqid, value=TRUE)
    gbid <- sapply(strsplit(gbid, "|", fixed=TRUE), "[[", 2)

    ft <- lapply(records, "[[", "GBSeq_feature-table")
    gene.ft <- lapply(ft, subset, GBFeature_key == "gene")
    gene.qual <- lapply(gene.ft, "[[", "GBFeature_quals")
    gene.qual <- lapply(gene.qual, lapply, subset, GBQualifier_name == "gene")

    intervals <- lapply(gene.ft, "[[", "GBFeature_intervals")
    start <- lapply(intervals, lapply, "[[", "GBInterval_from")
    end <- lapply(intervals, lapply, "[[", "GBInterval_to")
    gene <- lapply(gene.qual, lapply, "[[", "GBQualifier_value")

    keep <- which(sapply(gene, length) > 0)
    gbid <- as.list(gbid[keep])
    gene <- gene[keep]
    start <- start[keep]
    end <- end[keep]

    gene <- mapply(rep, gene, lapply(start, sapply, length), SIMPLIFY=FALSE)
    gbid <- mapply(rep, gbid, lapply(lapply(start, sapply, length), sum), SIMPLIFY=FALSE)
    start <- lapply(start, unlist)
    end <- lapply(end, unlist)

    d <- mapply(data.table, gbid=gbid, start=start, end=end, gene=gene, 
           SIMPLIFY=FALSE)
    d <- rbindlist(d)
    d[,start := as.integer(start)]
    d[,end := as.integer(end)]
    d
}
