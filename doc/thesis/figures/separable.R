#!/usr/bin/env Rscript

library(ape)
library(png)
library(raster)
library(extrafont)

set.seed(0)

tree.col <- c("black", rep(c("blue", "red"), each=3))
tmp <- tempfile()
tinytree <- as.list(rep(0, 7))
for (i in 1:7) {
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
x <- c(0.2, 0.3, 0.1, 0.55, 0.75, 0.65)
y <- c(0.65, 0.75, 0.7, 0.4, 0.45, 0.35)
rasterImage(tinytree[[1]],0.5,0.5,0.6,0.6)
mapply(rasterImage, tinytree[2:7], as.list(x), as.list(y), as.list(x+0.1), as.list(y+0.1))
abline(a=0.1, b=1, lty=2)
text(0.05, 0.95, "A", cex=2)

plot(NA, xlim=c(0, 1), ylim=c(0, 1), axes=FALSE, xlab=NA, ylab=NA, frame.plot=TRUE)
x <- (c(1, 0.5, -0.5, -1, -0.5, 0.5) + 1) * 0.45
y <- (c(0, sqrt(3)/2, sqrt(3)/2, 0, -sqrt(3)/2, -sqrt(3)/2, 0) + 1) * 0.45
x <- x[c(1, 3, 5, 2, 4, 6)]
y <- y[c(1, 3, 5, 2, 4, 6)]
rasterImage(tinytree[[1]],0.5,0.5,0.6,0.6)
mapply(rasterImage, tinytree[2:7], as.list(x), as.list(y), as.list(x+0.1), as.list(y+0.1))
text(0.05, 0.95, "B", cex=2)

plot(NA, xlim=c(0, 1), ylim=c(0, 1), axes=FALSE, xlab=NA, ylab=NA, frame.plot=TRUE)
x <- (c(1, 0.5, -0.5, -1, -0.5, 0.5) + 2) * 0.1
y <- (c(0, sqrt(3)/2, sqrt(3)/2, 0, -sqrt(3)/2, -sqrt(3)/2, 0) + 2) * 0.1
x <- x[c(1, 3, 5, 2, 4, 6)]
y <- y[c(1, 3, 5, 2, 4, 6)]
rasterImage(tinytree[[1]],0.5,0.5,0.6,0.6)
mapply(rasterImage, tinytree[2:7], as.list(x), as.list(y), as.list(x+0.1), as.list(y+0.1))
text(0.05, 0.95, "C", cex=2)

dev.off()
