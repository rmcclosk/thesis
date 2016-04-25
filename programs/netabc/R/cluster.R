#' Find the maximum phylogenetic cluster size in a tree.
#'
#' @param tree the tree to cluster
#' @param cutoff genetic distance less than this cutoff indicates phylogenetic
#'               relatedness
#' @return the size of the largest phylogenetic cluster
#' @export
csize.max <- function(tree, cutoff=0.02)
{
    cph <- cophenetic(tree)
    adj <- cph <= cutoff
    g <- graph.adjacency(adj, mode="undirected", diag=FALSE)
    max(clusters(g)$csize)
}
