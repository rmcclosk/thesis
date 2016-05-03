#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(netabc))
library(RColorBrewer)
library(grid)

epidemic.mst <- function (g, nI) {
    t <- mst(g)
    I <- vcount(g)

    while (I > nI) {
        v <- which(degree(t) == 1)[1]
        t <- delete_edges(t, E(t)[inc(v)])
        I <- I - 1
    }

    inodes <- which(degree(t) > 0)
    iedges <- get.edge.ids(g, c(t(as_edgelist(t))))
    list(nodes=inodes, edges=iedges)
}

n <- 40
t <- 12
alpha <- c(0, 2)
m <- 2

br <- list(rep(1, t - 1), c(1, rep(0.1, t-1)))
tree <- lapply(mapply(rcoal, list(t), br=br, SIMPLIFY=FALSE), ladderize)
net <- mapply(sample_pa, list(n), power=as.list(alpha), m=m, directed=FALSE, SIMPLIFY=FALSE)

epi <- lapply(net, epidemic.mst, t)
inodes <- lapply(epi, "[[", "nodes")
iedges <- lapply(epi, "[[", "edges")
ecol <- mapply(rep, list(NA), lapply(net, ecount), SIMPLIFY=FALSE)
ecol <- mapply(function (x, i) { x[i] <- "red"; x}, ecol, iedges, SIMPLIFY=FALSE)
vcol <- mapply(rep, list("black"), lapply(net, ecount), SIMPLIFY=FALSE)
vcol <- mapply(function (x, i) { x[i] <- "red"; x}, vcol, inodes, SIMPLIFY=FALSE)

ewidth <- mapply(rep, list(0), lapply(net, ecount), SIMPLIFY=FALSE)
ewidth <- mapply(function (x, i) { x[i] <- 3; x}, ewidth, iedges, SIMPLIFY=FALSE)

lay <- lapply(net, layout_nicely)

pdf("kernel-idea.pdf", width=6, height=4, family="Gillius ADF")
par(mfcol=c(2, 3), mar=c(0, 2, 0, 2))

# plot networks
mapply(plot, net, edge.color="gray", vertex.color="black", vertex.label=NA,
       vertex.size=6, edge.width=3, layout=lay)

# plot epidemics
mapply(function (g, lay, ecol, ewidth, vcol) {
    plot(g, edge.color="gray", vertex.color="black", vertex.label=NA,
       vertex.size=6, edge.width=3, layout=lay)
    plot(g, edge.color=ecol, vertex.color=vcol, vertex.frame.color=vcol, 
         vertex.label=NA, vertex.size=6, edge.width=ewidth, layout=lay,
         add=TRUE)
}, net, lapply(net, layout_nicely), ecol, ewidth, vcol)

# plot trees
par(mar=c(2, 2, 2, 2))
mapply(plot, tree, edge.color="red", direction="down", show.tip.label=FALSE,
       edge.width=3)

# arrows
grid.lines(x=unit(c(0.3, 0.37), "npc"), y=unit(c(0.25, 0.25), "npc"), 
           arrow=arrow(), gp=gpar(lwd=3))
grid.lines(x=unit(c(0.3, 0.37), "npc"), y=unit(c(0.75, 0.75), "npc"), 
           arrow=arrow(), gp=gpar(lwd=3))
grid.lines(x=unit(c(0.63, 0.7), "npc"), y=unit(c(0.25, 0.25), "npc"), 
           arrow=arrow(), gp=gpar(lwd=3))
grid.lines(x=unit(c(0.63, 0.7), "npc"), y=unit(c(0.75, 0.75), "npc"), 
           arrow=arrow(), gp=gpar(lwd=3))
dev.off()
