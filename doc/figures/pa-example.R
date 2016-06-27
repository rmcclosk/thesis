#!/usr/bin/env Rscript

library(igraph)
library(extrafont)

set.seed(1)
g0.small <- sample_pa(50, m=2, power=0, directed=FALSE)
g1.small <- sample_pa(50, m=2, power=1, directed=FALSE)
g2.small <- sample_pa(50, m=2, power=2, directed=FALSE)

pdf("pa-example.pdf", width=6, height=3)
par(mfrow=c(1, 3))
par(mar=c(2, 1, 2, 1) + 0.1)
plot(g0.small, vertex.label=NA, vertex.color="black", vertex.size=5)
plot(g1.small, vertex.label=NA, vertex.color="black", vertex.size=5)
plot(g2.small, vertex.label=NA, vertex.color="black", vertex.size=5)
dev.off()
