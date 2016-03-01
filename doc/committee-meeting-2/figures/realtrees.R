#!/usr/bin/env Rscript

library(ape)
library(extrafont)

bctree <- read.tree(Sys.glob("../../../simulations/abc-pa-bctree/subtree/*"))
cntree <- read.tree(Sys.glob("../../../simulations/abc-pa-cn/tree/*"))

pdf("realtrees.pdf", height=3, family="Gillius ADF")
par(mfrow=c(1, 2), mar=c(0, 0, 2, 0) + 0.1)
plot(bctree, type="fan", show.tip.label=FALSE, main="BC Cluster 0", edge.color="grey40")
plot(cntree, type="fan", show.tip.label=FALSE, main="China CRF07", edge.color="grey40")
dev.off()
