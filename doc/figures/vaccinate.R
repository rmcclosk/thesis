#!/usr/bin/env Rscript

library(igraph)

seed <- 111
n <- 20
nI <- 10
nI2 <- 5

repeat {
    set.seed(seed)
    g <- sample_gnm(n, n)
    if (length(subcomponent(g, 1)) == n) break
    seed <- seed + 1
    print(seed)
}

I <- sample(1:vcount(g), 1)
for (i in 1:nI) {
    I <- c(I, sample(setdiff(as.list(V(g)[nei(I)]), I), 1))
    if (i == nI2) {
        I2 <- unlist(I)
    }
}
I <- unlist(I)

layout <- layout.auto(g)

vcol <- rep("black", vcount(g))
vcol[I] <- "red"

ecol <- rep("black", ecount(g))
ecol[E(g)[I %--% I]] <- "red"

pdf("vaccinate.pdf")
par(mfrow=c(2, 1), mar=c(1, 0, 1, 0))
plot(g, vertex.label=NA, edge.color=ecol, vertex.color=vcol,
     vertex.frame.color=vcol, edge.width=2, layout=layout)

vcol <- rep("black", vcount(g))
vcol[I2] <- "red"

ecol <- rep("black", ecount(g))
ecol[E(g)[I2 %--% I2]] <- "red"

plot(g, vertex.label=NA, edge.color=ecol, vertex.color=vcol, 
     vertex.frame.color=vcol, edge.width=2, layout=layout)
dev.off()
