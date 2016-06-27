#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(netabc))
library(RColorBrewer)

n <- 40
t <- 8
net <- sample_pa(n, m=2, directed=FALSE)
trans.tree <- ladderize(rcoal(t))
phylo <- ladderize(rtree(t))

pdf("objects.pdf", width=7, height=2)
par(mfrow=c(1, 3), mar=c(0, 0, 0, 0))
plot(net, vertex.label=NA, vertex.size=3, vertex.color="black")
plot(trans.tree, show.tip.label=FALSE, direction="down", edge.width=2)
plot(phylo, type="unrooted", show.tip.label=FALSE, edge.width=2)
dev.off()
