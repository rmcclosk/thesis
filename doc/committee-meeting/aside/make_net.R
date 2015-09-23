#!/usr/bin/env Rscript

library(igraph)
library(ape)
g <- barabasi.game(100, m=2, directed=FALSE)
g <- as.directed(g, "mutual")
V(g)$remove <- 0
E(g)$transmit <- 1
write.graph(g, "ba.net", "gml")

clus <- cluster_edge_betweenness(g, directed=FALSE)
clus <- tapply(1:vcount(g), membership(clus), identity, simplify=FALSE)        
subg <- mapply(induced_subgraph, list(g), clus, impl="create_from_scratch",       
               SIMPLIFY=FALSE)                                                   
                                                                                 
# keep only clusters which have non-zero clustering coefficient                  
keep.clus <- clus[head(order(-sapply(subg, transitivity)), 1)]                  

rbw <- rainbow(length(keep.clus))                                                
keep.memb <- mapply(match, list(1:vcount(g)), keep.clus)                       
keep.memb <- unlist(ifelse(!is.na(apply(keep.memb, 1, any)),                     
                           apply(keep.memb, 1, function (x) which(!is.na(x))),   
                           NA))                                                  
keep.edge <- do.call(c, lapply(keep.clus, function (x) E(g)[x %--% x]))        
                                                                                 
edge.col <- rep("gray", ecount(g))                                             
edge.col[keep.edge] <- "black" 

pdf("ba-net.pdf")
par(mar=c(0, 0, 0, 0))
plot(g, 
     edge.arrow.size=0,
     vertex.label=ifelse(is.na(keep.memb), NA, V(g)),
     vertex.size=ifelse(is.na(keep.memb), 3, 10),                                
     vertex.color=rbw[keep.memb],                                                
     edge.color=edge.col)
dev.off()

system("nettree ba.net ba.nwk")

t <- read.tree("ba.nwk")
t$tip.label <- as.character(as.integer(t$tip.label)+1)

tip.col <- rep("black", 100)
tip.col[match(keep.clus[[1]], as.integer(t$tip.label))] <- "red"
tip.cex <- rep(par("cex"), 100)
tip.cex[match(keep.clus[[1]], as.integer(t$tip.label))] <- tip.cex[match(keep.clus[[1]], as.integer(t$tip.label))] * 1.5
tip.edge <- which(as.integer(t$tip.label[t$edge[,2]]) %in% keep.clus[[1]])
edge.col <- rep("black", nrow(t$edge))
edge.col[tip.edge] <- "red"

pdf("ba-tree.pdf")
plot(t, type="fan", tip.color=tip.col, edge.color=edge.col)
dev.off()
