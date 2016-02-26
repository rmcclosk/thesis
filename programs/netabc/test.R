#!/usr/bin/env Rscript

library(netabc)

data(iris)

d <- as.matrix(dist(iris[,1:4]))
y <- iris[,5]

netabc::dsvm.cv(d, y)
