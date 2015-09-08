#!/usr/bin/env Rscript

# Create a random network of the specified type.

suppressPackageStartupMessages(library(igraph))
options(warn=-1)

args <- commandArgs(trailingOnly=TRUE)
seed <- as.integer(args[1])
net.type <- args[2]
nnode <- as.integer(args[3])
mean.degree <- as.numeric(args[4])
remove.rate <- as.numeric(args[5])
transmissibility <- as.numeric(args[6])
gml.file <- args[7]
ws.prob <- as.numeric(args[8])

set.seed(seed)

g <- switch(net.type,
        WS = watts.strogatz.game(1, nnode, mean.degree/2, ws.prob),
        ER = erdos.renyi.game(nnode, mean.degree/nnode),
        BA = barabasi.game(nnode, m=mean.degree/2, directed=FALSE))

V(g)$remove <- remove.rate
E(g)$transmit <- transmissibility / (1-transmissibility) * remove.rate
g <- as.directed(g, 'mutual')
write.graph(g, gml.file, format='gml')
