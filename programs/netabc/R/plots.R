#' Plot a 2-dimensional PCA projection of a kernel-matrix 
#'
#' @param kmat kernel matrix
#' @param color either NULL, or a list of one element to display with color
#' @param shape either NULL, or a list of one element to display with point shape
#' @param yaml either NULL, or a YAML-formatted string to title the plot with
#' @return a ggplot object
#' @export
kpca.plot <- function (kmat, color=NULL, shape=NULL, yaml=NULL)
{
    kmat <- as.kernelMatrix(kmat)
    plot.data <- as.data.frame(rotated(kpca(kmat, features=2)))
    colnames(plot.data) <- c("PC1", "PC2")
    if (!is.null(color)) {
        plot.data <- cbind(plot.data, color)
    }
    if (!is.null(shape)) {
        plot.data <- cbind(plot.data, shape)
    }
    plot.aes <- aes_string(x="PC1", y="PC2", color=names(color), shape=names(shape))

    p <- ggplot(plot.data, plot.aes) + geom_point() + theme_bw()
    if (!is.null(yaml)) {
        p <- p + ggtitle(paste(strwrap(as.yaml(yaml.load(yaml)), 60), collapse="\n"))
    }
    p
}

#' Summarize data by colors, shapes, and faceting
#'
#' The color and shape of the points, as well as which points have lines drawn
#' between them, is dictated by the group parameter.
#'
#' @param data data to summarize
#' @param x name of column to be plotted on the x-axis
#' @param y name of column to be plotted on the y-axis
#' @param facet.x name of column to be faceted in columns
#' @param facet.y name of column to be faceted in rows
#' @param group name of column used to group points, or NULL
#' @param fun summary statistic to display
#' @param x.factor display the x-axis variable as a factor
#' @param y.factor display the y-axis variable as a factor
#' @param group.factor display the grouping variable as a factor
#' @param yaml either NULL, or a YAML-formatted string to title the plot with
#' @return a ggplot object
#' @export
summarize.plot <- function (data, x, y, facet.x=".", facet.y=".", group=NULL,
                          fun="mean", x.factor=TRUE, y.factor=FALSE,
                          group.factor=TRUE, yaml=NULL)
{
    agg.rhs <- paste(x, group, facet.x, facet.y, sep="+")
    agg.formula <- as.formula(paste0(y, "~", agg.rhs))
    plot.data <- aggregate(agg.formula, data, fun)

    if (y.factor) {
        plot.data[,y] <- as.factor(plot.data[,y])
    }
    if (x.factor) {
        plot.data[,x] <- as.factor(plot.data[,x])
    }
    if (!is.null(group) & group.factor) {
        plot.data[,group] <- as.factor(plot.data[,group])
    }

    p <- ggplot(plot.data, aes_string(x=x, y=y, color=group, shape=group, group=group))
    if (facet.x != "." | facet.y != ".") {
        p <- p + facet_grid(as.formula(paste(facet.y, "~", facet.x)), labeller="label_both")
    }
    p <- p + geom_point(size=3) + geom_line() + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    if (!is.null(yaml)) {
        p <- p + ggtitle(paste(strwrap(as.yaml(yaml.load(yaml)), 60), collapse="\n"))
    }
    p
}

#' Plot contact network clusters in a phylogeny
#'
#' This takes a contact network and a phylogeny, where the phylogeny was
#' presumably simulated from the contact network. The node names in the network
#' and phylogeny are matched. The nodes of the contact network should have the
#' "cluster" attribute defined, which indicates which cluster they are part of.
#' These clusters are highlighted with colours in the phylogeny.
#'
#' @param net contact network whose nodes have a "cluster" attribute
#' @param tree phylogeny with tip and node labels from the contact network
#' @param yaml YAML-formatted metadata to title the tree with
#' @param palette color brewer palette to use
#' @param status.only if TRUE, use only one colour to indicate cluster status
#' @param ... other arguments to be passed to plot.phylo
#' @export
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

#' Plot a phylogeny with edges coloured by a label
#'
#' @param tree phylogeny to plot
#' @param labels a matrix or data.frame, where the first two columns are the
#'               node name and the corresponding label
#' @param yaml YAML-formatted metadata to title the tree with
#' @param palette color brewer palette to use
#' @param parent if TRUE, color the edges according to the parent's state,
#'               otherwise the child's
#' @param ... other arguments to be passed to plot.phylo
#' @export
edge.color.plot <- function (tree, labels, yaml="", palette="PuBu", parent=TRUE, ...)
{
    labels <- as.data.frame(labels)
    colnames(labels)[1:2] <- c("node", "rate")
    labels$rate <- as.integer(as.factor(labels$rate))
    colors <- brewer.pal(max(labels$rate), palette)

    tree.labels <- c(tree$tip.label, tree$node.label)
    nodes <- tree.labels[tree$edge[,if (parent) 1 else 2]]
    edge.labels <- with(labels, rate[match(nodes, node)])
    edge.col <- colors[edge.labels]
    plot.phylo(tree, edge.color=edge.col, ...)
    title(paste(strwrap(as.yaml(yaml.load(yaml)), 60), collapse="\n"))
}

#' Plot marginal kernel densities
#'
#' @param d data.frame to plot densities for
#' @param truth true parameter values to indicate on plots
#' @param limits axis limits for plots
#' @export
marginal.plot <- function (d, truth, limits) {                                   
    vary.cols <- colnames(d)[apply(d, 2, function (x) length(unique(x)) > 1)]       
                                                                                 
    if (length(vary.cols) >= 2) {                                                
    # 2D marginals                                                               
        combos <- combn(vary.cols, 2)                                            
        plot.data <- apply(combos, 2, function (c) d[, c, with=FALSE])           
        plot.aes <- apply(combos, 2, function (c) aes_string(x=c[1], y=c[2]))       
        plot.limits <- apply(combos, 2, function (c) limits[c])                  
        plots <- mapply(ggplot, plot.data, plot.aes, SIMPLIFY=FALSE)             
        plots <- mapply(function (p, lim) {                                      
            p + stat_density2d(aes(fill=..level..), geom="polygon") + theme_bw() +
                guides(fill=FALSE) +                                             
                geom_point(data=truth, color="black", size=8) +                  
                geom_point(data=truth, color="white", size=6) +                  
                xlim(lim[[1]][1], lim[[1]][2]) +                                 
                ylim(lim[[2]][1], lim[[2]][2])                                   
        }, plots, plot.limits, SIMPLIFY=FALSE)                                   
        do.call(grid.arrange, c(plots, ncol=ceiling(sqrt(ncol(combos))),         
                                top="2D marginals"))                             
    }                                                                            
                                                                                 
    # 1D marginals                                                               
    plot.aes <- lapply(vary.cols, function (c) aes_string(x=c))                  
    plots <- mapply(ggplot, list(d), plot.aes, SIMPLIFY=FALSE)                   
    plots <- mapply(function (p, col) {                                          
        p + geom_density() + theme_bw() + xlim(limits[[col]][1], limits[[col]][2]) +
            geom_vline(xintercept=truth[,col], linetype=2)                       
    }, plots, vary.cols, SIMPLIFY=FALSE)                                         
    do.call(grid.arrange, c(plots, ncol=ceiling(sqrt(length(vary.cols))),        
                            top="1D marginals"))                                 
}
