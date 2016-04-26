#' Plot a ramasr object
#'
#' Plot a ramasr demographic projection
#' 
#' @param x a ramasr object
#' @param ... other plot arguments
#' @author Skipton Woolley
#' @examples 
#' tmat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' matsd <- tmat/10
#' v0 <- c(80,20,0)
#' sim_t10_rep100 <- demo_proj_n(v0=v0,tmat=tmat,matsd=matsd,estdem=TRUE,time=10,nrep=100) 
#' plot(sim_t10_rep100)
#' @export
plot.ramasr <- function (x, ...){
    x <- x$vn
    nrep <- length(x)
    stages <- dim(x[[1]])[1]
    time <- dim(x[[1]])[2]
    par(mfrow = c(1,stages))
    for (i in 1:stages) matplot(0:(time - 1), sapply(x,function(st) st[i, ]), 
                                xlab = "time", ylab = "abundance", 
                                type = 'l', col = rainbow(stages)[i], pch = 1, 
                                main = paste("stage", i),...)
}