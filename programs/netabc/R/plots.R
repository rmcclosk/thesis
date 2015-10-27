kpca.plot <- function (kmat, color=NULL, shape=NULL, yaml=NULL)
{
    kmat <- as.kernelMatrix(kmat)
    plot.data <- as.data.frame(rotated(kpca(kmat, features=2)))
    colnames(plot.data) <- c("PC1", "PC2")
    if (!is.null(color))
        plot.data <- cbind(plot.data, color)
    if (!is.null(shape))
        plot.data <- cbind(plot.data, shape)
    plot.aes <- aes_string(x="PC1", y="PC2", color=names(color), shape=names(shape))

    p <- ggplot(plot.data, plot.aes) + geom_point() + theme_bw()
    if (!is.null(yaml))
        p <- p + ggtitle(as.yaml(yaml.load(yaml)))
    p
}

summary.plot <- function (data, x, y, facet.x=".", facet.y=".", group=NULL,
                          fun="mean", x.factor=TRUE, y.factor=FALSE,
                          group.factor=TRUE)
{
    agg.rhs <- paste(x, group, facet.x, facet.y, sep="+")
    agg.formula <- as.formula(paste0(y, "~", agg.rhs))
    plot.data <- aggregate(agg.formula, data, fun)

    if (y.factor)
        plot.data[,y] <- as.factor(plot.data[,y])
    if (x.factor)
        plot.data[,x] <- as.factor(plot.data[,x])
    if (!is.null(group) & group.factor)
        plot.data[,group] <- as.factor(plot.data[,group])

    p <- ggplot(plot.data, aes_string(x=x, y=y, color=group, group=group))
    if (facet.x != "." | facet.y != ".") {
        p <- p + facet_grid(as.formula(paste(facet.x, "~", facet.y)), labeller="label_both")
    }
    p + geom_point() + geom_line() + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

cluster.plot <- function (net, tree, yaml="", palette="Set1", status.only=FALSE, ...)
{
    d <- data.frame(node=as.character(V(net)$id), cluster=as.integer(V(net)$cluster))
    bg <- d$cluster == 0
    if (status.only) {
        d$cluster[!bg] <- 2
    }
    else {
        d$cluster[!bg] <- (d$cluster[!bg] %% 7) + 2
    }
    d$cluster[bg] <- 1
    cluster.colors <- c("#bdbdbd", brewer.pal(max(d$cluster)-1, palette))
    
    labels <- c(tree$tip.label, tree$node.label)
    parents <- labels[tree$edge[,1]]
    children <- labels[tree$edge[,2]]

    parent.cluster <- with(d, cluster[match(parents, node)])
    child.cluster <- with(d, cluster[match(children, node)])
    if (status.only) {
        edge.col <- cluster.colors[pmin(parent.cluster, child.cluster)]
    }
    else {
        edge.col <- ifelse(parent.cluster == child.cluster | parent.cluster == 0 |
                           child.cluster == 0, cluster.colors[parent.cluster], "#000000")
    }

    plot.phylo(tree, edge.color=edge.col, ...)
    title(paste(strwrap(as.yaml(yaml.load(yaml)), 60), collapse="\n"))
}

pcbr.plot <- function (tree, pcbr.out, yaml="", palette="PuBu", ...)
{
    names(pcbr.out) <- c("node", "rate", "cluster")
    pcbr.out$rate <- factor(pcbr.out$rate, levels=sort(unique(pcbr.out$rate)))
    pcbr.out$rate <- as.integer(pcbr.out$rate)
    colors <- brewer.pal(max(pcbr.out$rate), palette)

    labels <- c(tree$tip.label, tree$node.label)
    parents <- labels[tree$edge[,1]]
    parent.rate <- with(pcbr.out, rate[match(parents, node)])
    edge.col <- colors[parent.rate]
    plot.phylo(tree, edge.color=edge.col, ...)
    title(paste(strwrap(as.yaml(yaml.load(yaml)), 60), collapse="\n"))
}
