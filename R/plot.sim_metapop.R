#' Plot a sim_metapop object
#'
#' Plot a sim_metapop object.
#' 
#' @param x a sim_metapop object
#' @param ... other plot arguments
#' @author Skipton Woolley
#' @examples 
#' n <- 50
#' meta_data <- data.frame(x1=runif( n, min=-10, max=10), x2=runif( n, min=-10, max=10),
#' area=exp(-seq(.1,10,length.out = n))*10,presence=rbinom(n,1,.8))
#' area <- meta_data$area
#' dist <- as.matrix(with(meta_data, dist(cbind(x1, x2))))
#' presence <- meta_data$presence
#' mp <- sim_metapop(nrep=100, time=50, dist, area, presence, x = 0.42, e = 0.061, y = 15.2)
#' plot(mp)
#' @export
plot.sim_metapop <- function(x,...){
    # graphics::par(mfrow = c(1,2))
    # M <- base::as.matrix(exp(-x$alpha * x$dist))
    # base::diag(M) <- 0
    # M <- M * base::outer(x$area, x$area)
    # tmp <- base::eigen(M)
    # y <- tmp$vectors[, 1]^2
    # cl <-  (1 - y/max(y))*254 + 1
    # dia <- sqrt(max(y));
    # plot(locations, asp=1, cex=sqrt(x$area)/dia, xlab="", ylab="", pch=21,col="blue", bg=cl)
    graphics::matplot(0:x$time, base::sapply(x$mp, function(zz) base::apply(zz,2, base::sum)),
                      type = 'l', xlab = "time", ylab = "abundance",pch = 1,col="#00000030", ...)
    graphics::lines(0:(x$time), x$sim_p_obs,type = 'l', col='red',lwd=2)  
}