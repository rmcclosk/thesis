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

npt <- 5
x <- c(rnorm(npt, sd=0.2), rnorm(npt, mean=2, sd=0.2), rnorm(npt, mean=4, sd=0.2))
y <- rnorm(npt*3)
png("tinypca.png", bg="transparent", width=500, height=100)
par(mar=c(0, 0, 0, 0) + 0.1)
plot(x, y, pch=16, cex=3, col=rep(col, each=npt))
box(lwd=4)
dev.off()

x <- c(rnorm(npt, sd=0.2), rnorm(npt, mean=2, sd=0.2), rnorm(npt, mean=4, sd=0.2))
y <- c(rnorm(npt, mean=1, sd=0.2), rnorm(npt, mean=0, sd=0.2), rnorm(npt, mean=1, sd=0.2))
png("tinyksvr.png", bg="transparent", width=100, height=100, type="cairo")
par(mar=c(0, 0, 0, 0) + 0.1)
plot(splinefun(x=c(0, 2, 4), y=c(1, 0, 1)), lwd=4, xlim=c(-0.5, 4.5))
plot(splinefun(x=c(0, 2, 4), y=c(1, 0, 1)), lwd=32, add=TRUE, xlim=c(-0.5, 4.5), 
     col=rgb(0, 0, 0, alpha=0.5))
points(x, y, pch=16, cex=2, col=rep(col, each=npt))
box(lwd=4)
dev.off()

x <- c(rnorm(npt, sd=0.2), rnorm(npt, mean=2, sd=0.2), rnorm(npt, mean=4, sd=0.2))
y <- c(rnorm(npt, mean=0, sd=0.2), rnorm(npt, mean=0.33, sd=0.2), rnorm(npt, mean=1, sd=0.2))
png("tinysvr.png", bg="transparent", width=100, height=100, type="cairo")
par(mar=c(0, 0, 0, 0) + 0.1)
plot(x, y, pch=16, cex=2, col=rep(col, each=npt), type="n")
abline(a=0, b=0.2, lwd=4)
abline(a=0, b=0.2, lwd=32, col=rgb(0, 0, 0, alpha=0.5))
points(x, y, pch=16, cex=2, col=rep(col, each=npt))
box(lwd=4)
dev.off()

x <- c(rnorm(npt, sd=0.2), rnorm(npt, mean=2, sd=0.2), rnorm(npt, mean=4, sd=0.2))
y <- c(rnorm(npt, mean=0, sd=0.2), rnorm(npt, mean=0.5, sd=0.2), rnorm(npt, mean=1, sd=0.2))
png("tinyreg.png", bg="transparent", width=100, height=100)
par(mar=c(0, 0, 0, 0) + 0.1)
plot(x, y, pch=16, cex=2, col=rep(col, each=npt))
abline(a=0, b=0.25, lwd=4)
box(lwd=4)
dev.off()
