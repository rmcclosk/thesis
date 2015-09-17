SI.net <- function(net, transmit.rate=1, mode="mutual")
{
    SIR.net(net, transmit.rate=transmit.rate, remove.rate=0, mode=mode)
}

SIR.net <- function(net, transmit.rate=1, remove.rate=0, mode="mutual")
{
    net <- as.directed(net, mode=mode)
    E(net)$transmit <- transmit.rate
    V(net)$remove <- remove.rate
    net
}
