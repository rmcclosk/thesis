#!/usr/bin/env Rscript

library(igraph)

n <- 15

set.seed(1)

g1 <- sample_pa(n, m=2, power=0, directed=FALSE)
g2 <- sample_pa(n, m=2, power=2, directed=FALSE)

pdf("nettypes.pdf")
par(mfrow=c(2, 1), mar=c(1, 0, 1, 0))
plot(g1, vertex.label=NA, vertex.color="black", edge.color="black",
     edge.width=2)
plot(g2, vertex.label=NA, vertex.color="black",
     edge.color="black", edge.width=2)
dev.off()
