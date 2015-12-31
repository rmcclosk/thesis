#' Create one of the group-type from Goodreau 2004.
#'
#' @param n number of vertices in the graph
#' @param num.groups number of equally sized assortative groups, must divide n
#' @param prop.intra proportion of edges which are reserved for intra-group
#'                   contacts
#' @return an undirected graph
#' @export
goodreau.net.groups <- function (n, num.groups, prop.intra)
{
    if (n %% num.groups != 0) {
        stop("Group size must divide number of nodes evenly")
    }

    g <- make_empty_graph(n, directed=FALSE)
    group.size <- n / num.groups

    num.intra.ties <- rep(0, num.groups)
    if (prop.intra * n >= 1) {
        for (i in 1:(prop.intra*n)) {
            group <- sample(num.groups, 1)
            num.intra.ties[group] <- num.intra.ties[group] + 1
        }
    
        intra.pairs <- 1:group.size^2
        head <- (intra.pairs - 1) %/% group.size + 1
        tail <- (intra.pairs - 1) %% group.size + 1
        intra.pairs <- intra.pairs[head < tail]
    
        intra.ties <- c()
        for (i in 1:num.groups) {
            group.start <- (i-1) * group.size + 1
            pairs <- sample(intra.pairs, num.intra.ties[i])
            head <- (pairs - 1) %/% group.size + group.start
            tail <- (pairs - 1) %% group.size + group.start
            intra.ties <- c(intra.ties, c(rbind(head, tail)))
        }
        g <- add_edges(g, intra.ties)
    }

    if ((1-prop.intra) * n >= 1) {
        group.pairs <- 1:num.groups^2
        head <- (group.pairs - 1) %/% num.groups + 1
        tail <- (group.pairs - 1) %% num.groups + 1
        group.pairs <- group.pairs[head < tail]
    
        num.inter.ties <- matrix(0, nrow=num.groups, ncol=num.groups)
        for (i in 1:((1-prop.intra)*n)) {
            if (length(group.pairs) > 1) {
                groups <- sample(group.pairs, 1)
            } else {
                groups <- group.pairs
            }
            head.group <- (groups - 1) %/% num.groups + 1
            tail.group <- (groups - 1) %% num.groups + 1
            num.inter.ties[head.group, tail.group] <- num.inter.ties[head.group, tail.group] + 1
        }
    
        inter.pairs <- 1:group.size^2
        inter.ties <- c()
        for (i in 1:num.groups) {
            head.start <- (i-1)*group.size + 1
            for (j in 1:num.groups) {
                tail.start <- (j-1)*group.size + 1
                pairs <- sample(inter.pairs, num.inter.ties[i,j])
                head <- (pairs - 1) %/% group.size + head.start
                tail <- (pairs - 1) %% group.size + tail.start
                inter.ties <- c(inter.ties, c(rbind(head, tail)))
            }
        }
        g <- add_edges(g, inter.ties)
    }

    V(g)$group <- paste0("group", ((1:n)-1) %/% group.size)
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
    if (rate <= 0) {
        return (g)
    }
    new.ties <- c(t(as_edgelist(complementer(g))))
    add_edges(net, new.ties, attr=list(transmit=rep(rate, length(new.ties)/2)))
}
