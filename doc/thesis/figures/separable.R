#!/usr/bin/env Rscript

library(ape)
library(png)
library(raster)
library(extrafont)

set.seed(4)

tree.col <- c("black", rep(c("blue", "red", "forestgreen"), each=3))
tmp <- tempfile()
tinytree <- as.list(rep(0, 10))
for (i in 1:10) {
    png(tmp, width=50, height=50, bg="transparent")
    plot(rtree(4), show.tip.label=F, direction="down", no.margin=TRUE, edge.width=5,
         edge.color=tree.col[i])
    dev.off()
    tinytree[[i]] <- readPNG(tmp)
}
unlink(tmp)

pdf("separable.pdf", width=6, height=2, family="Gillius ADF")
par(mar=c(0, 0, 0, 0) + 0.1, mfrow=c(1, 3))

plot(NA, xlim=c(0, 1), ylim=c(0, 1), axes=FALSE, xlab=NA, ylab=NA, frame.plot=TRUE)
x <- c(rep(c(0.2, 0.4, 0.8), each=3) + rnorm(9, sd=0.1))
y <- c(rep(c(0.7, 0.4, 0.6), each=3) + rnorm(9, sd=0.1))
rasterImage(tinytree[[1]],0.5,0.5,0.6,0.6)
mapply(rasterImage, tinytree[2:10], as.list(x), as.list(y), as.list(x+0.1), as.list(y+0.1))
text(0.05, 0.95, "A", cex=2)

plot(NA, xlim=c(0, 1), ylim=c(0, 1), axes=FALSE, xlab=NA, ylab=NA, frame.plot=TRUE)
x <- c(rep(c(0.2, 0.4, 0.8), 3) + rnorm(9, sd=0.2))
y <- c(rep(c(0.7, 0.4, 0.6), 3) + rnorm(9, sd=0.2))
rasterImage(tinytree[[1]],0.5,0.5,0.6,0.6)
mapply(rasterImage, tinytree[2:10], as.list(x), as.list(y), as.list(x+0.1), as.list(y+0.1))
text(0.05, 0.95, "B", cex=2)

plot(NA, xlim=c(0, 1), ylim=c(0, 1), axes=FALSE, xlab=NA, ylab=NA, frame.plot=TRUE)
x <- c(rep(0.3, 9) + rnorm(9, sd=0.08))
y <- c(rep(0.3, 9) + rnorm(9, sd=0.08))
rasterImage(tinytree[[1]],0.5,0.5,0.6,0.6)
mapply(rasterImage, tinytree[2:10], as.list(x), as.list(y), as.list(x+0.1), as.list(y+0.1))
text(0.05, 0.95, "C", cex=2)

dev.off()
