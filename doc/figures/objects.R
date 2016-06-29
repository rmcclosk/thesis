#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(igraph))
library(png)

n <- 40
t <- 8
net <- sample_pa(n, m=2, directed=FALSE)
trans.tree <- ladderize(rcoal(t))
phylo <- ladderize(rtree(t))

pdf("objects.pdf", width=7, height=2)
par(mfrow=c(1, 4), mar=c(0, 0, 0, 0))
plot(net, vertex.label=NA, vertex.size=3, vertex.color="black")
plot(trans.tree, show.tip.label=FALSE, direction="down", edge.width=2)
plot(phylo, type="unrooted", show.tip.label=FALSE, edge.width=2)
plot.new()
y <- seq(0.1, 0.9, length.out=5)
segments(0.25, y, 1, y, lwd=2)

nmut <- 12
xmut <- sample(seq(0.25, 1, 0.01), nmut, replace=TRUE)
ymut <- sample(y, nmut, replace=TRUE)
segments(xmut, ymut-0.03, xmut, ymut+0.03, lwd=2)

virus <- readPNG("stock/virus.png")
rasterImage(virus, 0, y-0.08, 0.18, y+0.08)

dev.off()
