#!/usr/bin/env Rscript

# Calculate the Sackin index for all simulated trees.

library(apTreeshape)

args <- commandArgs(trailingOnly=TRUE)
dirs <- head(args, -1)
tsv.file <- tail(args, 1)

get.sackin <- function (tree.file)
    tryCatch({sackin(as.treeshape(read.tree(tree.file)))},
        warning = function (w) {0}, error = function (e) {0})

tree.files <- list.files(dirs, full.names=TRUE, recursive=TRUE, pattern="*.nwk")

d <- do.call(rbind, strsplit(dirname(tree.files), .Platform$file.sep))
d <- setNames(as.data.frame(d), c("net.type", "T"))
d$sackin <- sapply(tree.files, get.sackin)
write.table(d, tsv.file, quote=FALSE, row.names=FALSE, sep="\t")
