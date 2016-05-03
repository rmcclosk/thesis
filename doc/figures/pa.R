#!/usr/bin/env Rscript

library(igraph)
set.seed(4)

nv <- c(3, 4, 5, 6, 12, 12)
nI <- 7

g <- list(NULL)
for (n in nv) {
    g[[length(g)+1]] <- sample_pa(n, m=2, power=1, directed=FALSE, 
                                  start.graph=tail(g, 1)[[1]])
}
g <- tail(g, length(g)-1)
fullg <- tail(g, 1)[[1]]

t <- mst(fullg)
I <- vcount(fullg)
while (I > nI) {
    v <- which(degree(t) == 1)[1]
    t <- delete_edges(t, E(t)[inc(v)])
    I <- I - 1
}

inodes <- which(degree(t) > 0)
iedges <- get.edge.ids(fullg, c(t(as_edgelist(t))))

vcol <- mapply(rep, list(c("black", "white")), 
               as.list(as.data.frame(rbind(nv, max(nv)-nv))), SIMPLIFY=FALSE)
vcol[[length(vcol)]][inodes] <- "red"
ecol <- mapply(rep, list("black"), lapply(g, ecount), SIMPLIFY=FALSE)
ecol[[length(ecol)]][iedges] <- "red"
vsize <- mapply(rep, list(5), as.list(nv), SIMPLIFY=FALSE)
vsize[[length(vsize)]][inodes] <- 10

lay <- layout_in_circle(tail(g, 1)[[1]])
pdf("pa.pdf")
par(mar=c(0, 0, 0, 0) + 0.1)
mapply(plot, g, vertex.color=vcol, vertex.frame.color=vcol,
       vertex.label=list(NA), vertex.size=vsize, edge.color=ecol,
       edge.width=2, layout=list(lay))
dev.off()
