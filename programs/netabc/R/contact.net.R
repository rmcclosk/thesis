SI.net <- function(net, transmit.rate=1, mode="mutual")
{
    SIR.net(net, transmit.rate=transmit.rate, remove.rate=0, mode=mode)
}

SIR.net <- function(net, transmit.rate=1, remove.rate=0, mode="mutual")
{
    net <- as.directed(net, mode=mode)
    E(net)$transmit <- transmit.rate
    V(net)$remove <- remove.rate
    net
}

add.transmission.clusters <- function (net, cluster.size, num.clusters, 
                                       cluster.rate=2, connect=FALSE)
{
    clusters <- list()
    while (length(clusters) < num.clusters)
    {
        cluster <- sample.int(vcount(net), 1)
        while(cluster %in% unlist(clusters))
        {
            cluster <- cluster + 1
            if (cluster > vcount(net))
                cluster <- 1
        }
        while (length(cluster) < cluster.size)
        {
            nbhd <- table(unlist(adjacent_vertices(g, cluster, mode="out")))
            nbhd <- nbhd[!as.integer(names(nbhd)) %in% cluster]
            nbhd <- nbhd[!as.integer(names(nbhd)) %in% unlist(clusters)]
            if (length(nbhd) == 0) break
            cluster <- c(as.integer(names(nbhd)[which.max(nbhd)]), cluster)
        }
        clusters[[length(clusters)+1]] <- cluster
    }
    
    pairs <- mapply(expand.grid, clusters, clusters, SIMPLIFY=FALSE)
    pairs <- lapply(pairs, subset, Var1 != Var2)
    pairs <- lapply(pairs, function (x) c(t(x)))
    if (connect)
        for (p in pairs) {g <- add.edges(g, p)}

    cluster.edges <- unlist(mapply(get.edge.ids, list(g), pairs, SIMPLIFY=FALSE))
    cluster.edges <- unique(cluster.edges[cluster.edges != 0])
    E(g)[cluster.edges]$transmit <- cluster.rate
    g
}
