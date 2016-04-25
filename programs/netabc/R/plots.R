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

#' Label wrapped facets.
#' http://stackoverflow.com/questions/11979017/changing-facet-label-to-math-formula-in-ggplot2
#' @param gg.plot plot to add facet labels to
#' @param labels labels to add
#' @export
facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1

  g <- ggplotGrob(gg.plot)
  gg <- g$grobs
  strips <- grep("strip_t", names(gg))

  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text",
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }

  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g))
  g
}

#' Make a tiny tree icon.
#'
#' @param ntip number of tips in the tree
#' @param color color for the tree
#' @param ultra if TRUE, make an ultrametric tree
#' @param lwd line width
#' @param ... extra options for plot.phylo
#' @return an image to pass to rasterImage
#' @export
tiny.tree <- function (ntip=4, color="black", ultra=FALSE, lwd=8)
{
    tf <- tempfile()
    t <- ifelse(ultra, rcoal, rtree)(ntip)
    png(tf, bg="transparent", width=100, height=100)
    par(mar=c(0, 0, 0, 0) + 0.1)
    plot(t, direction="down", edge.width=lwd, edge.color=color,
         show.tip.label=FALSE)
    dev.off()
    pic <- readPNG(tf)
    unlink(tf)
    pic
}

#' Convert polar co-ordinates to rectangular.
#'
#' @param r radii of the points
#' @param theta angles of the points
#' @return a two-column matrix with the x and y co-ordinates
#' @export
polar2rect <- function (r, theta) {
    cbind(r*cos(theta), r*sin(theta))
}

#' Plot posterior distributions from a netabc trace for the Barabasi-Albert
#' preferential attachment model.
#'
#' @param trace data.frame output by netabc
#' @export
posterior.plot.pa <- function (trace, I_min=500, I_max=10000, N_min=500, N_max=10000,
                               alpha_min=0, alpha_max=2, m_min=1, m_max=5,
                               show.map=TRUE, true_alpha=NA, true_m=NA, true_I=NA, true_N=NA) {
    setDT(trace)
    d <- trace[iter == max(iter)]
    d <- d[sample(1:nrow(d), prob=weight, replace=TRUE)]
    font.size <- 14
    options(scipen=-1)

    denspoly <- function (x, q=0.95) {
        dens <- density(x)
        dens <- data.frame(x=dens$x, y=dens$y)
        dens <- dens[dens$x >= quantile(x, 1-q) & dens$x <= quantile(x, q),]
        rbind(c(dens[1, "x"], 0), dens, c(dens[nrow(dens), "x"], 0))
    }
    map <- function (x) {
        dens <- density(x)
        dens$x[which.max(dens$y)]
    }

    plot.theme <- theme_bw() +
                  theme(text=element_text(family="Gillius ADF", size=font.size),
                        axis.ticks.y=element_blank(),
                        axis.text.y=element_blank(),
                        legend.position="none")
    palpha <- ggplot(d, aes(x=alpha)) + 
            geom_density() +
            labs(x=expression(alpha), y="") + 
            xlim(alpha_min, alpha_max) +
            geom_polygon(data=denspoly(d[,alpha]), aes(x=x, y=y), fill="black", alpha=0.3) +
            plot.theme
    if (show.map) {
        palpha <- palpha + geom_vline(xintercept=d[,map(alpha)])
    }
    if (!is.na(true_alpha)) {
        palpha <- palpha + geom_vline(xintercept=true_alpha, linetype="dashed")
    }

    pI <- ggplot(d, aes(x=I)) + 
            geom_density() +
            geom_polygon(data=denspoly(d[,I]), aes(x=x, y=y), fill="black", alpha=0.3) +
            labs(x="I", y="") + 
            xlim(I_min, I_max) +
            plot.theme
    if (show.map) {
        pI <- pI + geom_vline(xintercept=d[,map(I)])
    }
    if (!is.na(true_I)) {
        pI <- pI + geom_vline(xintercept=true_I, linetype="dashed")
    }

    pN <- ggplot(d, aes(x=N)) + 
            geom_density() + 
            labs(x="N", y="") +
            geom_polygon(data=denspoly(d[,N]), aes(x=x, y=y), fill="black", alpha=0.3) +
            xlim(N_min, N_max) +
            plot.theme
    if (show.map) {
        pN <- pN + geom_vline(xintercept=d[,map(N)])
    }
    if (!is.na(true_N)) {
        pN <- pN + geom_vline(xintercept=true_N, linetype="dashed")
    }

    d[,shade := m >= quantile(m, 0.05) & m <= quantile(m, 0.95)]
    pm <- ggplot(d, aes(x=floor(m))) +
            geom_histogram(aes(alpha=shade), fill="black", col="black", binwidth=1) +
            labs(x="m", y="") +
            xlim({m_min}-0.5, {m_max}+0.5) +
            scale_alpha_manual(values=c("TRUE"=0.3, "FALSE"=0)) +
            plot.theme
    if (show.map) {
        pm <- pm + geom_vline(xintercept=d[,map(m)])
    }
    if (!is.na(true_m)) {
        pm <- pm + geom_vline(xintercept=true_m, linetype="dashed")
    }

    arrangeGrob(palpha, pI, pm, pN, ncol=2)
}
