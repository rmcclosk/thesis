#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(extrafont)

set.seed(0)

dens1 <- density(rnorm(400))
dens2 <- density(runif(400, -3, 3))
dens3 <- density(rnorm(400))
dens4 <- density(runif(400, -3, 3))

d <- data.table(x=c(dens1$x, dens2$x, dens3$x, dens4$x), 
                y=c(dens1$y, dens2$y, dens3$y, dens4$y), 
                precise=rep(c("yes", "no", "yes", "no"), each=length(dens1$x)),
                accurate=rep(c("yes", "no"), each=length(dens1$x)*2))
d[,precise := factor(precise, levels=c("yes", "no"))]
d[,accurate := factor(accurate, levels=c("yes", "no"))]
d[,col := ifelse(precise == "yes", ifelse(accurate == "yes", "green", "red"), "orange")]

lines <- data.frame(
    v=c(0, dens2$x[which.max(dens2$y)], 2, 1),
    precise=c("yes", "no", "yes", "no"),
    accurate=c("yes", "yes", "no", "no")
)

pdf("acc-prec.pdf", height=4)
ggplot(d, aes(x=x, y=y, color=col)) + geom_line() + theme_bw() + xlim(-4, 4) +
    geom_vline(data=lines, aes(xintercept=v), lty=2) +
    facet_grid(precise~accurate, scales="free", labeller="label_both") +
    theme(axis.ticks=element_blank(), legend.position="none",
          panel.grid=element_blank(), axis.text=element_blank()) +
    scale_color_manual(values=c("green", "orange", "red")) +
    labs(x="testing value", y="kernel score")
dev.off()
