make.nltt.depths <- function (depths, stepfun = TRUE)                            
{                                                                                
    # normalize by maximum depth                                                 
    depths <- depths / max(depths)                                               
                                                                                 
    # get number of lineages at each depth                                       
    lin <- cumsum(as.integer(table(depths))) - 1                                 
                                                                                 
    # normalize lineages to be in [0, 1]                                         
    lin <- lin / max(lin)                                                           
                                                                                 
    # make the nLTT functions                                                    
    if (stepfun) {                                                               
        stepfun(unique(depths), c(0, lin))                                       
    }                                                                               
    else {                                                                       
        approxfun(unique(depths), lin)                                           
    }                                                                            
}                                                                                
                                                                                    
make.ltt <- function (t, stepfun = TRUE)                                         
{                                                                                
    # get node depths in sorted order; ignore tips                               
    depths <- node.depth.edgelength(t)                                           
    depths <- sort(depths[(t$Nnode+2):(2*t$Nnode+1)])                            
    depths <- round(depths, 8)                                                      
    make.nltt.depths(depths, stepfun)                                               
}                                                                                   
                                                                                    
#' Calcluate the normalized lineages-through-time statistic.
#'
#' @param t1 first tree
#' @param t2 second tree
#' @param stepfun if FALSE, use a linear interpolation instead of a step function
#' @param ... additional arguments passed to integrate
#' @export
nLTT <- function (t1, t2, stepfun = TRUE, ...)                        
{                                                                                   
    f1 <- make.ltt(t1, stepfun)                                                     
    f2 <- make.ltt(t2, stepfun)                                                  
    f <- function (x) abs(f1(x) - f2(x))                                            
    integrate(f, 0, 1, ...)$value                                                   
}
