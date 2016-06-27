#!/usr/bin/env Rscript

library(igraph)

seed <- 381

n <- 10
vcol <- rep(c("blue", "magenta"), each=n)

repeat {{
    print(seed)
    set.seed(seed)
    g <- bipartite.random.game(n, n, type="gnm", m=2*n)
    g2 <- rewire(g, each_edge(p = .1, loops = FALSE))
    if (min(degree(g2)) > 0) break
    seed <- seed+1
}} 

seed <- 390
repeat {{
    print(seed)
    set.seed(seed)
    g3 <- bipartite.random.game(n, n, type="gnm", m=2*n)
    g3 <- rewire(g3, each_edge(p = .1, loops = FALSE))
    if (min(degree(g3)) > 0) break
    seed <- seed+1
}} 

layout <-  layout.auto(g2)
pdf("homophily.pdf")
plot(g, vertex.label=NA, vertex.color=vcol, vertex.frame.color=vcol,
     edge.color=NA, edge.width=NA, layout=layout)
plot(g, vertex.label=NA, vertex.color=vcol, vertex.frame.color=vcol,
     edge.color="black", edge.width=2, layout=layout)
plot(g2, vertex.label=NA, vertex.color=vcol, vertex.frame.color=vcol,
     edge.color="black", edge.width=2, layout=layout)
plot(g3, vertex.label=NA, vertex.color=vcol, vertex.frame.color=vcol,
     edge.color="black", edge.width=2)
dev.off()
