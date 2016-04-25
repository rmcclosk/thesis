#!/usr/bin/env Rscript

library(ape)
library(extrafont)

t1 <- read.tree(text="(((:1,:1):1,:1):1,:1):0;")
t2 <- read.tree(text="((((:1,:1):1,:1):1,:1):1,:1):0;")

sf1 <- stepfun(c(0, 0.5), c(0, 0.5, 1))
sf2 <- stepfun(c(0, 1/3, 2/3), c(0, 1/3, 2/3, 1))

af1 <- approxfun(c(0, 0.5, 1), c(0, 0.5, 1))
af2 <- approxfun(c(0, 1/3, 2/3, 1), c(0, 1/3, 2/3, 1))

pdf("nltt.pdf", family="Gillius ADF")
par(mfrow=c(2, 2), lwd=3)
old.mar <- par("mar")
par(mar=c(0, 2, 4, 2) + 0.1)
plot(t1, edge.width=3, edge.color="blue", show.tip.label=FALSE, direction="down")
plot(t2, edge.width=3, edge.color="red", show.tip.label=FALSE, direction="down")

par(mar=old.mar)
plot(sf1, col="blue", main="original nLTT", xlab="normalized lineages",
     ylab="normalized time", do.points=FALSE, xlim=c(0, 1), xaxs="i")
lines(sf2, col="red", do.points=FALSE)

plot(af1, col="blue", main="updated nLTT", xlab="normalized lineages",
     ylab="normalized time", xlim=c(0, 1), xaxs="i")
x <- seq(0, 1, 0.1)
lines(x, af2(x) + 0.05, col="red")
dev.off()
