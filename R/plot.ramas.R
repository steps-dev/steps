#' Plot a ramas object
#'
#' Plot a ramas demographic projection
#' 
#' @param x a ramas object
#' @param mean_pop logical. If TRUE plots the mean demographic change in population through time. Otherwise it plots the population changes for each stage in the transition matrix.
#' @param ... other plot arguments
#' @author Skipton Woolley
#' @examples 
#' tmat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' matsd <- tmat/10
#' v0 <- c(80,20,0)
#' sim_t20_rep100 <- ramas(v0=v0,tmat=tmat,matsd=matsd,estdem=TRUE,time=20,nrep=100) 
#' plot(sim_t20_rep100,mean_pop=TRUE)
#' plot(sim_t20_rep100,mean_pop=FALSE)
#' @export
plot.ramas <- function(x,mean_pop=TRUE,...){
    x <- x$vn
    nrep <- base::length(x)
    stages <- base::dim(x[[1]])[1]
    time <- base::dim(x[[1]])[2]
    if(mean_pop==TRUE){
    graphics::par(mfrow = c(1,1))
    graphics::matplot(0:(time - 1), (base::sapply(x, function(re) base::apply(re,2, base::sum))),
            type = 'l', xlab = "time", ylab = "abundance",pch = 1,col="#00000030", ...)
    graphics::lines(0:(time - 1), base::apply(base::sapply(x, function(re) base::apply(re,2, sum)), 1, base::mean),
         type = 'l', col='red',lwd=2)  
    } else {
    graphics::par(mfrow = c(1,stages))
    for (i in 1:stages) {graphics::matplot(0:(time - 1), base::sapply(x,function(st) st[i, ]), 
                                          xlab = "time", ylab = "abundance", 
                                          type = 'l', col="#00000030", pch = 1, 
                                          main = base::paste("stage", i),...)
                        graphics::lines(0:(time - 1), base::apply(base::sapply(x,function(st) st[i, ]),1,mean),
                                        type = 'l', col=base::suppressWarnings(RColorBrewer::brewer.pal(base::length(base::unique(stages)),
                                                                                                        "Set1"))[i],lwd=2)
    }
    }
}