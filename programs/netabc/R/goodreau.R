#' Create one of the networks from Goodreau 2004.
#'
#' @param n number of vertices in the graph
#' @param net.type type of network to create, one of "random", "2-weak", 
#'                 "2-strong", "8-weak", "8-strong", "core-weak", "core-strong",
#'                 or "bridges"
#' @return an undirected graph
#' @export
goodreau.net <- function (n, net.type) {
    if (net.type == "random") {
        g <- sample_gnm(n, n)
    } else if (net.type == "2-weak") {
        g <- make_empty_graph(n, directed=FALSE)
        edges <- combn(n, 2)
        intra.edges <- edges[,! xor(edges[1,] <= n/2, edges[2,] <= n/2)]
        inter.edges <- edges[,xor(edges[1,] <= n/2, edges[2,] <= n/2)]
        intra.ties <- c(intra.edges[,sample(ncol(intra.edges), 0.75 * n)])
        inter.ties <- c(inter.edges[,sample(ncol(inter.edges), 0.25 * n)])
        g <- add_edges(g, intra.ties)
        g <- add_edges(g, inter.ties)
    } else if (net.type == "2-strong") {
        g <- make_empty_graph(n, directed=FALSE)
        edges <- combn(n, 2)
        intra.edges <- edges[,! xor(edges[1,] <= n/2, edges[2,] <= n/2)]
        inter.edges <- edges[,xor(edges[1,] <= n/2, edges[2,] <= n/2)]
        intra.ties <- c(intra.edges[,sample(ncol(intra.edges), 0.975 * n)])
        inter.ties <- c(inter.edges[,sample(ncol(inter.edges), 0.025 * n)])
        g <- add_edges(g, intra.ties)
        g <- add_edges(g, inter.ties)
    } else if (net.type == "8-weak") {
        g <- make_empty_graph(n, directed=FALSE)
        edges <- combn(n, 2)
        which.intra <- (edges[1,] - 1) %/% (n/8) == (edges[2,] - 1) %/% (n/8)
        intra.edges <- edges[,which.intra]
        inter.edges <- edges[,!which.intra]
        intra.ties <- c(intra.edges[,sample(ncol(intra.edges), 0.75 * n)])
        inter.ties <- c(inter.edges[,sample(ncol(inter.edges), 0.25 * n)])
        g <- add_edges(g, intra.ties)
        g <- add_edges(g, inter.ties)
    } else if (net.type == "8-strong") {
        g <- make_empty_graph(n, directed=FALSE)
        edges <- combn(n, 2)
        which.intra <- (edges[1,] - 1) %/% (n/8) == (edges[2,] - 1) %/% (n/8)
        intra.edges <- edges[,which.intra]
        inter.edges <- edges[,!which.intra]
        intra.ties <- c(intra.edges[,sample(ncol(intra.edges), 0.975 * n)])
        inter.ties <- c(inter.edges[,sample(ncol(inter.edges), 0.025 * n)])
        g <- add_edges(g, intra.ties)
        g <- add_edges(g, inter.ties)
    } else if (net.type == "core-strong") {
        g <- make_empty_graph(n, directed=FALSE)
        edges <- combn(n, 2)
        intracore.edges <- edges[,edges[1,] <= n/8 & edges[2,] <= n/8]
        intraper.edges <- edges[,edges[1,] > n/8 & edges[2,] > n/8]
        inter.edges <- edges[,xor(edges[1,] <= n/8, edges[2,] <= n/8)]
        intracore.ties <- c(intracore.edges[,sample(ncol(intracore.edges), 0.5 * n)])
        intraper.ties <- c(intraper.edges[,sample(ncol(intraper.edges), 0.475 * n)])
        inter.ties <- c(inter.edges[,sample(ncol(inter.edges), 0.025 * n)])
        g <- add_edges(g, intracore.ties)
        g <- add_edges(g, intraper.ties)
        g <- add_edges(g, inter.ties)
    } else if (net.type == "core-weak") {
        g <- make_empty_graph(n, directed=FALSE)
        edges <- combn(n, 2)
        intracore.edges <- edges[,edges[1,] <= n/8 & edges[2,] <= n/8]
        intraper.edges <- edges[,edges[1,] > n/8 & edges[2,] > n/8]
        inter.edges <- edges[,xor(edges[1,] <= n/8, edges[2,] <= n/8)]
        intracore.ties <- c(intracore.edges[,sample(ncol(intracore.edges), 0.15 * n)])
        intraper.ties <- c(intraper.edges[,sample(ncol(intraper.edges), 0.75 * n)])
        inter.ties <- c(inter.edges[,sample(ncol(inter.edges), 0.1 * n)])
        g <- add_edges(g, intracore.ties)
        g <- add_edges(g, intraper.ties)
        g <- add_edges(g, inter.ties)
    } else if (net.type == "bridges") {
        g <- make_empty_graph(n, directed=FALSE)
        edges <- combn(n, 2)
        csw.edges <- edges[,edges[1,] <= 0.475 * n & edges[2,] >= 0.95 * n]
        marriage.ties <- c(rbind(1:(0.475*n), (1:(0.475*n)) + 0.475*n))
        csw.ties <- c(csw.edges[,sample(ncol(csw.edges), 0.525*n)])
        g <- add_edges(g, marriage.ties)
        g <- add_edges(g, csw.ties)
    } else {
        stop("Unrecognized network type")
    }
    g
}

#' Add a baseline transmission rate for all directed edges not in a contact
#' network.
#'
#' @param net network to which a baseline rate will be added
#' @param rate the baseline transmission rate (should be much lower than the
#' normal transmission rate)
#' @return the modified network
#' @export
add.baseline.rate <- function (net, rate = 0.01)
{
    new.ties <- c(t(as_edgelist(complementer(g))))
    add_edges(net, new.ties, attr=list(transmit=rep(rate, length(new.ties)/2)))
}
