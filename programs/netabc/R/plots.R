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

    ggplot(plot.data, aes_string(x=x, y=y, color=group, group=group)) + 
        facet_grid(as.formula(paste(facet.x, "~", facet.y)), labeller="label_both") +
        geom_point() + geom_line() + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}
