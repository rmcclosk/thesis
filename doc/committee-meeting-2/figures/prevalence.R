#!/usr/bin/env Rscript

library(igraph)

col <- '#6baed6'
n <- 30

g <- replicate(3, sample_pa(n, m=2, power=1, directed=F), simplify=FALSE)

mst.node <- list(1:(n%/%5), 1:(n%/%2), 1:n)
subg <- mapply(induced_subgraph, g, mst.node, SIMPLIFY=FALSE)
mst <- lapply(subg, minimum.spanning.tree)
mst.edge <- lapply(mst, as_edgelist)
mst.edge <- sapply(lapply(mst.edge, t), c)
mst.edge <- mapply(get.edge.ids, g, mst.edge, SIMPLIFY=FALSE)

edge.col <- mapply(replicate, lapply(g, ecount), col, SIMPLIFY=FALSE)
edge.col <- mapply(function (e, m) ifelse(1:length(e) %in% m, "black", e), 
                   edge.col, mst.edge, SIMPLIFY=FALSE)

vertex.col <- mapply(replicate, lapply(g, vcount), col, SIMPLIFY=FALSE)
vertex.col <- mapply(function (e, m) ifelse(1:length(e) %in% m, "black", e), 
                     vertex.col, mst.node, SIMPLIFY=FALSE)

pdf("prevalence.pdf", height=2)
par(mfrow=c(1, 3), mar=c(0, 0, 0, 0) + 0.1)
mapply(plot, g, edge.color=edge.col, vertex.label=NA, vertex.size=8, 
       edge.width=4, vertex.color=vertex.col, vertex.frame.color=vertex.col)
dev.off()
