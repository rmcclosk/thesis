#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(netabc))

m <- collect.metadata("../../../simulations/abc-pa-bctree/abc/*.tsv")
bc.file <- rownames(subset(m, m == 2 & rbf_variance == 5))
cn.file <- Sys.glob("../../../simulations/abc-pa-cn/abc/*.tsv")
stopifnot(length(bc.file) == 1)
stopifnot(length(cn.file) == 1)

get.posterior <- function (file) {
    d <- fread(paste("head -n -1", file), select=c("iter", "weight", "attach_power"))
    d <- subset(d, iter == d[,max(iter)])
    idx <- sample(1:nrow(d), prob=d[,weight], replace=TRUE)
    d[idx,]
}

bc.data <- get.posterior(bc.file)
cn.data <- get.posterior(cn.file)

bc.data <- bc.data[,"attach_power",with=FALSE][,data := "BC"]
cn.data <- cn.data[,"attach_power",with=FALSE][,data := "China"]
d <- rbind(bc.data, cn.data)

bc.dens <- bc.data[,density(attach_power)]
cn.dens <- cn.data[,density(attach_power)]

bc.max <- bc.dens$x[which.max(bc.dens$y)]
cn.max <- cn.dens$x[which.max(cn.dens$y)]

lines <- data.frame(attach_power=c(bc.max, cn.max), data=c("BC", "China"),
                    density=c(max(bc.dens$y), max(cn.dens$y)))
lines$label <- round(lines$attach_power, 2)

pdf("bc_vs_cn.pdf", family="Gillius ADF")
ggplot(d, aes(x=attach_power)) + geom_density() + facet_grid(data~.) +
    geom_vline(data=lines, lty=2, aes(xintercept=attach_power)) +
    geom_text(data=lines, aes(x=attach_power - 0.05, y=density + 0.05, label=label), 
              hjust=1, vjust=0) +
    labs(x="preferential attachment power", y="posterior density") +
    theme_bw() +
    theme(text = element_text(size = 20))
dev.off()
