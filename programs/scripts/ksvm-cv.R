#!/usr/bin/env Rscript

# Perform cross-validations of a kernel SVM.

library(kernlab)
library(optparse)
source(file="util.R")

options <- list(
    make_option(c("-n", "--number"), type="integer", default=1000, 
                help="Number of cross-validations"),
    make_option(c("-s", "--seed"), type="integer", default=NULL,
                help="Random seed"),
    make_option(c("-t", "--statistic"), default="accuracy",
                help="Outcome(s) to report"),
    make_option(c("-o", "--overwrite"), type="boolean", default=FALSE,
                help="Overwrite output file if it exists")
)

parser <- OptionParser(usage = "%prog [options] <matrix file> <responses file>",
                       option_list = options)
args <- parse_args(parser, positional_arguments = 3)

matrix.file <- args$args[1]
response.file <- args$args[2]
tsv.file <- args$args[3]

if (file.access(matrix.file) == -1)
    stop(sprintf("Cannot access matrix file %s", matrix.file))
if (file.access(response.file) == -1)
    stop(sprintf("Cannot access response file %s", matrix.file))
if (!args$options$overwrite & file.access(tsv.file) == 0)
    stop(sprintf("Output file %s already exists", tsv.file))

stats <- strsplit(args$options$statistic, ",")[[1]]

set.seed(args$options$seed)
n.cv <- args$options$number
kmat <- read.mm(matrix.file)
y <- read.table(response.file)[,1]

if (is.factor(y)) 
{
    allowed.stats <- c("accuracy")
} else 
{
    allowed.stats <- c()
}

if (!all(stats %in% allowed.stats))
    stop(sprintf("Unrecognized statistics: %s", paste0(stats[!stats %in% allowed.stats], collapse=", ")))

accuracy <- function (pred, y)
    sum(pred == y) / length(y)

# http://stackoverflow.com/questions/1753299/help-using-predict-for-kernlabs-svm-in-r
results <- lapply(1:n.cv, function (i) {
    holdout <- sample.int(nrow(kmat), nrow(kmat)/2)
    
    train <- as.kernelMatrix(kmat[-holdout,-holdout])
    m <- ksvm(train, y[-holdout], kernel="matrix")
    
    test <- as.kernelMatrix(kmat[holdout, -holdout][,SVindex(m), drop=F])
    pred <- predict(m, test)
    sapply(stats, do.call, list(pred, y[holdout]))
})

results <- do.call(rbind, results)
write.table(results, tsv.file, row.names=FALSE, quote=FALSE)
