#' Plot a metapop_mod object
#'
#' Plot a metapop_mod object.
#' 
#' @param x a metapop_mod object
#' @param ... other plot arguments
#' @author Skipton Woolley
#' @examples 
#' n <- 50
#' meta_data <- data.frame(x1=runif( n, min=-10, max=10), x2=runif( n, min=-10, max=10),
#' area=exp(-seq(.1,10,length.out = n))*10,presence=rbinom(n,1,.8))
#' area <- meta_data$area
#' dist <- as.matrix(with(meta_data, dist(cbind(x1, x2))))
#' presence <- meta_data$presence
#' mp <- metapop_mod(nrep=100, time=50, dist, area, presence, x = 0.42, e = 0.061, y = 1.2)
#' plot(mp)
#' @export
plot.metapop_mod <- function(x,...){
    graphics::matplot(0:x$time, base::sapply(x$mp, function(zz) base::apply(zz,2, base::sum)),
                      type = 'l', lty = 1, xlab = "time", ylab = "population size",pch = 1,col="#00000030", ...)
    graphics::lines(0:(x$time), x$sim_p_obs,type = 'l', col='red',lwd=2)  
}