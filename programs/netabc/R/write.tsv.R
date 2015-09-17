write.tsv <- function (df, ...)
{
    write.table(df, ..., col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}
