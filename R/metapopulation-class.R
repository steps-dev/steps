#' @title metapopulation objects
#' @name metapopulation
#' @rdname metapopulation
#' @param nrep Number of simulations.
#' @param time Number of time steps
#' @param habitat object
#' @param y1 incidence function parameters
#' @param x1 incidence function parameters
#' @param e1 Minimum area of patches
#' @param alpha Exponential decay rate of patch connectivity (dispersion parameter)
#' @param beta double parameter that represents the shape of the dispersal kernel.
#' @param disp_fun char Which dispersal function to use. 'H' uses hanski(1994), if 'S' uses shaw(1995).
#' @param coln_fun char which colonisation function to use. 'H' = Hanski (1994), 'M'= Moilanen (2004), 'O'= Ovaskainen (2002). See decription for details.
#' @param locations NULL or NumericMatrix Longitudes and latitudes of coordinates of the patches
#' @param c double colonisation scale parameter see Ovaskainen 2002.
#' @export
#' @author Skipton Woolley
#' @seealso \code{link{metapop_n_cpp}}
#' @description do some metapopulation modelling in R
#' @examples 
#' habitat <- as.habitat(list(coordinates = data.frame(x=runif( 10, min=-10, max=10),
#'                                                     y=runif( 10, min=-10, max=10)),
#'                                area = as.data.frame(exp(-seq(.1,10,length.out = 10))*10),
#'                                population = as.population(t(rmultinom(10, 
#'                                size = 100, prob = c(0.8,0.2,0.01)))),
#'                                features = data.frame(temperature = 10)))
#' mp <- metapopulation(nrep=10, time=100, habitat,
#'                   x1 = 0.42, e1 = 0.0061, y1 = 1.2)

setGeneric("metapopulation",
           function(nrep=10, time=20, habitat,
                    y1 = 1, x1 = 1, e1=1, alpha = 1, beta = 1,
                    disp_fun = 'H',locations = NULL,
                    c=1,coln_fun='H') {
             if(!is.habitat(habitat))stop('metapopulation needs a habitat')
             dist <- as.matrix(distance(habitat))
             area <- as.numeric(unlist(area(habitat)))
             locations <- NULL
             pop <- population(habitat)# need to integrate demographic and presence correctly.
             presence <- as.vector(ifelse(pop[,ncol(pop)]>0,1,0))
             # call c++ function that does this loop.
             mp <- dlmpr::metapop_n_cpp(nrep=nrep, time=time, dist=dist, area=area, presence=presence,
                                        y = y1, x = x1, e = e1, alpha = alpha, beta = beta,
                                        disp_fun = disp_fun,locations = locations,
                                        c=c,coln_fun = coln_fun)
             sim_p_obs <- base::apply(plyr::ldply(base::lapply(mp,base::colSums),function(x)x),2,base::mean)
             sim_i_obs <- base::apply(plyr::ldply(base::lapply(base::lapply(mp,function(x)x[,-1]),function(x)base::rowSums(x)/time),function(x)x),2,base::mean)
             results <- base::list(mp = mp, sim_p_obs = sim_p_obs, sim_i_obs = sim_i_obs,
                             nrep=nrep, time = time, dist = dist, area = area, y = y1, x = x1, e = e1, 
                             alpha = alpha, beta=beta,locations = locations)
             base::class(results) <- "metapopulation"
             return(results)
           }
)

#' @rdname metapopulation
#' @name plot.metapopulation
#' @method plot metapopulation
#' @param x a metapopulation object
#' @param ... other plot arguments
#' @author Skipton Woolley
#' @examples 
#' n <- 50
#' meta_data <- data.frame(x1=runif( n, min=-10, max=10), x2=runif( n, min=-10, max=10),
#' area=exp(-seq(.1,10,length.out = n))*10,presence=rbinom(n,1,.8))
#' area <- meta_data$area
#' dist <- as.matrix(with(meta_data, dist(cbind(x1, x2))))
#' presence <- meta_data$presence
#' mp <- metapopulation(nrep=100, time=50, dist, area, presence, x = 0.42, e = 0.061, y = 1.2)
#' plot(mp)
#' @export
plot.metapopulation <- function(x,...){
  graphics::matplot(0:x$time, base::sapply(x$mp, function(zz) base::apply(zz,2, base::sum)),
                    type = 'l', lty = 1, xlab = "time", ylab = "population size",pch = 1,col="#00000030", ...)
  graphics::lines(0:(x$time), x$sim_p_obs,type = 'l', col='red',lwd=2)  
}