#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(netabc))

n <- 6
g <- SIR.net(sample_pa(n, m=2, power=1, directed=FALSE))
tmp.graph <- tempfile()
tmp.tree <- tempfile()
write.graph(g, tmp.graph, "gml")
system2("nettree", args=tmp.graph, stdout=tmp.tree)
t <- read.tree(tmp.tree)
t$tip.label <- as.character(as.integer(t$tip.label)+1)
t$node.label <- as.character(as.integer(t$node.label)+1)

pdf("nettree.pdf", height=4)
par(mar=c(0, 1, 0, 1) + 0.1, mfrow=c(1, 2))
plot(g, vertex.color='#eff3ff', vertex.label.cex=1.5, vertex.size=20,
     vertex.label.color="black", edge.width=2)
plot(t, show.node.label=TRUE, edge.color="gray50", direction="down", srt=90, 
     label.offset=0.05, cex=1.5, edge.width=2)
dev.off()
