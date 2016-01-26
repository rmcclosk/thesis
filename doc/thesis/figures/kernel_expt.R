#!/usr/bin/env Rscript

library(igraph)
library(ape)

nnode <- 7
ninf <- 4
set.seed(0)

pow <- c(0.5, 1, 1.5)
col <- c("red", "blue", "green")
inf.col <- "black"

for (i in 1:3) {
    g <- sample_pa(nnode, m=2, power=pow[i], directed=FALSE)
    png(sprintf("tinynet%d.png", i), bg="transparent", width=100, height=100)
    par(mar=c(0, 0, 0, 0) + 0.1)
    plot(g, vertex.size=10, vertex.label=NA, vertex.color=col[i],
         vertex.frame.color=col[i], edge.color=col[i], edge.width=5)
    dev.off()

    png(sprintf("tinyepi%d.png", i), bg="transparent", width=100, height=100)
    par(mar=c(0, 0, 0, 0) + 0.1)
    ecol <- rep(col[i], ecount(g))
    vcol <- rep(c(inf.col, col[i]), c(ninf, vcount(g)-ninf))
    ecol[sample(unique(unlist(g[[1:ninf,1:ninf,edges=TRUE]])), ninf-1)] <- inf.col
    plot(g, vertex.size=10, vertex.label=NA, vertex.color=vcol,
         vertex.frame.color=vcol, edge.color=ecol, edge.width=5)
    dev.off()

    png(sprintf("tinytree%d.png", i), bg="transparent", width=100, height=100)
    par(mar=c(0, 0, 0, 0) + 0.1)
    t <- rcoal(ninf)
    plot(t, direction="down", edge.width=5, edge.color=col[i], show.tip.label=FALSE)
    dev.off()
}
