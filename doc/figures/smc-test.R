#!/usr/bin/env Rscript

library(data.table)
library(extrafont)

trace.file <- system2(c("locate", "check_smc_trace.tsv"), stdout=TRUE)
trace <- fread(trace.file)

trace <- trace[iter == max(iter)]
trace <- trace[sample(1:nrow(trace), prob=weight, replace=TRUE)]
dens <- density(trace$theta)

x <- seq(-4, 4, 0.01)
target <- 0.5 * dnorm(x) + 0.5 * dnorm(x, sd=0.1)

pdf("smc-test.pdf", family="Gillius ADF", width=5, height=5)
plot(x, target, type="n", xlab="theta", ylab="density")
polygon(dens$x, dens$y, col="gray80")
lines(x, target)
dev.off()
