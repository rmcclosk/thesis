#!/usr/bin/env Rscript

library(ape)
library(igraph)

n <- 15

pdf("phylo.pdf", width=0.8, height=0.8)
plot.phylo(rcoal(n), show.tip.label=FALSE, direction="down", edge.width=2,
           no.margin=TRUE)
dev.off()

col <- c('#bdd7e7','#6baed6','#2171b5')
for (i in 1:3) {
    g <- sample_pa(n, m=2, power=(i-1)*0.8, directed=FALSE)
    pdf(paste0("net", i, ".pdf"), width=0.8, height=0.8)
    par(mar=c(0, 0, 0, 0) + 0.1)
    plot(g, vertex.label=NA, vertex.size=8, vertex.color=col[i],
         edge.width=2, vertex.frame.color=col[i], edge.color=col[i])
    dev.off()
    pdf(paste0("tree", i, ".pdf"), width=0.8, height=0.8)
    par(mar=c(0, 0, 0, 0) + 0.1)
    plot.phylo(rcoal(n), show.tip.label=FALSE, direction="down", edge.width=2,
               edge.color=col[i], no.margin=TRUE)
    dev.off()
}
