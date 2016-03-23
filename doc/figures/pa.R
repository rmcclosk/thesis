#!/usr/bin/env Rscript

library(igraph)

nv <- c(3, 4, 5, 6, 12)

g <- list(NULL)
for (n in nv) {
    g[[length(g)+1]] <- sample_pa(n, m=2, power=1, directed=FALSE, start.graph=tail(g, 1)[[1]])
}
g <- tail(g, length(g)-1)

vcol <- mapply(rep, list(c("black", "white")), 
               as.list(as.data.frame(rbind(nv, max(nv)-nv))), SIMPLIFY=FALSE)

lay <- layout_in_circle(tail(g, 1)[[1]])
pdf("pa.pdf")
mapply(plot, g, vertex.color=vcol, vertex.frame.color=vcol,
       vertex.label=list(NA), vertex.size=list(5), 
       layout=list(lay))
dev.off()
