#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(netabc))

post <- rnorm(500)
abc <- list(runif(500, -3, 3),
            c(runif(250, -3, 3), rnorm(250)),
            c(runif(100, -3, 3), rnorm(400)))
dens <- lapply(abc, density)
gradient <- colorRampPalette(c("steelblue1", "steelblue1"))(length(dens[[1]]$x)*1.2)

pdf("abc-smc.pdf", family="Gillius ADF", height=4, width=6)
par(mfrow=c(2, 3), mgp=c(0, 0, 0), mar=c(1, 1, 2, 1) + 0.1, xpd=TRUE)
for (i in 1:3) {
    plot(density(post), main=NA, ylab=bquote(f[.(i)](theta*"|"*D)), 
         xlab=expression(theta), tck=0, xaxt="n", yaxt="n", zero.line=FALSE,
         xlim=c(-4, 4) + 0.1)
    points(dens[[i]]$x, dens[[i]]$y, col=gradient, pch=16)
}

true.tree <- tiny.tree(4, ultra=TRUE, color="black")
tree.size <- 0.15
offset <- tree.size/2
epsilon <- c(0.45, 0.35, 0.25)

for (i in 1:3) {
    dist <- sample(dens[[i]]$x, 10, prob=dens[[i]]$y)
    r <- (dist - min(dens[[i]]$x)) / (max(dens[[i]]$x) - min(dens[[i]]$x))
    r <- abs(r - 0.5) * 1.2
    r[r < 0.2] <- rnorm(sum(r < 0.2), mean=0.2, sd=0.05)
    theta <- (1:length(r)) * 2 * pi / length(r)
    coords <- polar2rect(r, theta) + 0.5
    colors <- gradient[as.integer(r * length(gradient))]
    trees <- lapply(colors, function (x) tiny.tree(ntip=4, color=x, ultra=TRUE))
    plot.new()
    for (j in 1:length(trees)) {
        rasterImage(trees[[j]], coords[j,1]-offset, coords[j,2]-offset, 
                    coords[j,1]+offset, coords[j,2]+offset)
        if (r[j] > epsilon[i]) {
            points(coords[j,1], coords[j,2], pch=4, col="red", cex=4)
        }
    }
    rasterImage(true.tree, 0.5-offset, 0.5-offset, 0.5+offset, 0.5+offset)
    symbols(0.5, 0.5, circles=epsilon[i], lty=2, add=T, inches=F)
}

dev.off()
