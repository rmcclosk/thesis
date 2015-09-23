kpca.plot <- function (kmat, color=NULL, shape=NULL)
{
    kmat <- as.kernelMatrix(kmat)
    plot.data <- as.data.frame(rotated(kpca(kmat, features=2)))
    colnames(plot.data) <- c("PC1", "PC2")
    if (!is.null(color))
        plot.data <- cbind(plot.data, color)
    if (!is.null(shape))
        plot.data <- cbind(plot.data, shape)
    plot.aes <- aes_string(x="PC1", y="PC2", color=names(color), shape=names(shape))
    ggplot(plot.data, plot.aes) + geom_point() + theme_bw()
}

summary.plot <- function (data, x, y, group=NULL, fun="mean", x.factor=TRUE,
                          y.factor=FALSE, group.factor=TRUE)
{
    agg.formula <- as.formula(paste0(y, "~", x, "+", group))
    plot.data <- aggregate(agg.formula, data, fun)

    if (y.factor)
        plot.data[,y] <- as.factor(plot.data[,y])
    if (x.factor)
        plot.data[,x] <- as.factor(plot.data[,x])
    if (!is.null(group) & group.factor)
        plot.data[,group] <- as.factor(plot.data[,group])

    ggplot(plot.data, aes_string(x=x, y=y, color=group, group=group)) + 
        geom_point() + geom_line() + theme_bw()
}
