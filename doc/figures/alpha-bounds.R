#!/usr/bin/env Rscript

library(igraph)
library(extrafont)

set.seed(0)

g0.small <- sample_pa(50, m=2, power=0, directed=FALSE)
g2.small <- sample_pa(50, m=2, power=2, directed=FALSE)
g0.large <- sample_pa(5000, m=2, power=0, directed=FALSE)
g2.large <- sample_pa(5000, m=2, power=2, directed=FALSE)

summary(degree(g0.large))
summary(degree(g2.large))

pdf("alpha-bounds.pdf", family="Gillius ADF", width=6, height=6)
par(mfrow=c(2, 2))
par(mar=c(2, 2, 2, 2) + 0.1)

plot(g0.small, vertex.label=NA, vertex.color="black", vertex.size=5)
text(grconvertX(0.05, "ndc", "user"), grconvertY(0.95, "ndc", "user"), "A", cex=2)

plot(g2.small, vertex.label=NA, vertex.color="black", vertex.size=5)
text(grconvertX(0.55, "ndc", "user"), grconvertY(0.95, "ndc", "user"), "B", cex=2)

par(mar=c(4, 4, 2, 2) + 0.1, xpd=NA)

d0 <- density(degree(g0.large))
y0 <- d0$y[d0$x > 0]
x0 <- log(d0$x[d0$x > 0])
plot(x0, y0, type="l", main=NA, xlab="log(degree)", ylab="density")
text(grconvertX(0.05, "ndc", "user"), grconvertY(0.5, "ndc", "user"), "C", cex=2)

d2 <- density(degree(g2.large))
y2 <- d2$y[d2$x > 0]
x2 <- log(d2$x[d2$x > 0])
plot(x2, y2, type="l", main=NA, xlab="log(degree)", ylab="density")
text(grconvertX(0.55, "ndc", "user"), grconvertY(0.5, "ndc", "user"), "D", cex=2)

dev.off()
