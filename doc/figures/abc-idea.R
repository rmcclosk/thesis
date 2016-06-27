#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(netabc))
library(RColorBrewer)

col <- as.list(brewer.pal(3, "Set1"))
n <- 40
t <- 8
br <- list(rep(1, t - 1), c(rep(1, t/2), rep(0.1, t/2)), c(1, rep(0.1, t-1)))
tree <- lapply(mapply(rcoal, list(t), br=br, SIMPLIFY=FALSE), ladderize)
net <- mapply(sample_pa, list(n), power=list(0, 1, 2), m=2, directed=FALSE, SIMPLIFY=FALSE)

pdf("abc-idea.pdf", width=6, height=4)
par(mfrow=c(4, 4), mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.6, "sample", cex=2)
text(0.5, 0.4, "parameters", cex=2)
plot.new()
text(0.5, 0.5, expression(theta[1]), cex=4, col=col[[1]])
plot.new()
text(0.5, 0.5, expression(theta[2]), cex=4, col=col[[2]])
plot.new()
text(0.5, 0.5, expression(theta[3]), cex=4, col=col[[3]])

plot.new()
text(0.5, 0.6, "simulate", cex=2)
text(0.5, 0.4, "networks", cex=2)
mapply(plot, net, edge.color=col, vertex.color=col, vertex.label=NA, 
       vertex.size=6, edge.width=3)
par(mar=c(1, 1, 1, 1))
plot.new()
text(0.5, 0.6, "simulate", cex=2)
text(0.5, 0.35, "trees", cex=2)
mapply(plot, tree, edge.color=col, direction="down", show.tip.label=FALSE,
       edge.width=3)
plot.new()
text(0.5, 0.6, "compare to", cex=2)
text(0.5, 0.35, "observed tree", cex=2)
plot.new()
par(xpd=NA)
text(0.8, 0.5, expression(K(phantom()%.%phantom(),phantom()%.%phantom())),
     cex=3)
plot(ladderize(rcoal(t, br=br[[1]])), show.tip.label=FALSE, direction="down", edge.width=3)
dev.off()
