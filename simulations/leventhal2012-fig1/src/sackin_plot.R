#!/usr/bin/env Rscript

# Make a plot of Sackin index against transmissibility, stratified by network
# type.

library(ggplot2) 

args <- commandArgs(trailingOnly=TRUE)
data.file <- args[1]
pdf.file <- args[2]

q25 <- function (x) quantile (x, 0.25) 
q75 <- function (x) quantile (x, 0.75) 

data <- read.table(data.file, header=TRUE)

ggplot(data, aes(x=T, y=sackin, fill=net.type)) +
    stat_summary(fun.ymin=q25, fun.ymax=q75, geom='ribbon', alpha=0.5) +
    stat_summary(fun.y=median, geom='line', aes(color=net.type)) +
    theme(legend.position=c(1, 0), legend.justification=c(0, 1)) +
    theme_bw()
ggsave(pdf.file)
