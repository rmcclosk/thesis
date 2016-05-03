#!/usr/bin/env Rscript

library(igraph)
library(extrafont)

set.seed(0)
g0.small <- sample_pa(50, m=2, power=0, directed=FALSE)
g2.small <- sample_pa(50, m=2, power=2, directed=FALSE)

pdf("pa-example.pdf", family="Gillius ADF", width=6, height=3)
par(mfrow=c(1, 2))
par(mar=c(2, 2, 2, 2) + 0.1)
plot(g0.small, vertex.label=NA, vertex.color="black", vertex.size=5)
plot(g2.small, vertex.label=NA, vertex.color="black", vertex.size=5)
dev.off()
