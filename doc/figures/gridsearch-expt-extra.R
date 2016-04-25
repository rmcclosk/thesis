#!/usr/bin/env Rscript

library(igraph)
library(ape)
library(extrafont)

set.seed(0)
ninf <- 4
col <- c(rgb(0, 0, 1), rgb(0.4, 0.4, 1), rgb(0.8, 0.8, 1))

for (i in 1:3) {
    png(sprintf("testtree%d.png", i), bg="transparent", width=100, height=100)
    par(mar=c(0, 0, 0, 0) + 0.1)
    t <- rcoal(ninf)
    plot(t, direction="down", edge.width=5, edge.color=col[i], show.tip.label=FALSE)
    dev.off()
}

dens <- density(rnorm(200))
x <- dens$x
y <- dens$y
gradient <- colorRampPalette(c("blue", "white"))(length(x))
true.x <- x[length(x)*0.4]
true.x.col <- gradient[length(x)*0.4]
pdf("tinykscore.pdf", family="Gillius ADF", width=3, height=0.8)
par(mar=c(1, 1, 0, 0) + 0.1, mgp=c(0, 1, 0), cex=1.25)
plot(x, y, type="p", xlim=c(-4, 4), main=NA, pch=16,
     xlab="training value", ylab="score", col=gradient, xaxt="n", yaxt="n")
abline(v=true.x, col=true.x.col)
dev.off()
