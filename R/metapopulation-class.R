#' @title metapopulation objects
#' @name metapopulation
#' @rdname metapopulation
#' @param nrep Number of simulations.
#' @param time Number of time steps
#' @param habitat object
#' @param y1 incidence function parameters
#' @param x1 incidence function parameters
#' @param e1 Minimum area of patches
#' @param dispersal object see \link[bbgdm]{dispersal} an object that contains information on:
#' alpha, the decay rate of patch connectivity (dispersion parameter); beta, parameter that controls
#'  the shape of the dispersal kernel; and disp_fun, which selects the disperal function: 'H' uses hanski(1994)
# and 'S' uses shaw(1995).
#' @param coln_fun char which colonisation function to use. 'H' = Hanski (1994), 'M'= Moilanen (2004), 'O'= Ovaskainen (2002). See decription for details.
#' @param locations NULL or NumericMatrix Longitudes and latitudes of coordinates of the patches
#' @param c double colonisation scale parameter see Ovaskainen 2002.
#' @export
#' @author Skipton Woolley
#' @seealso \code{link{metapop_n_cpp}}
#' @description do some metapopulation modelling in R
#' @examples 
#' habitat <- as.habitat(list(coordinates = data.frame(x=runif( 100, min=-20, max=20),
#'                                                     y=runif( 100, min=-20, max=20)),
#'                                area = as.data.frame(exp(-seq(.1,10,length.out = 100))*10),
#'                                population = as.population(t(rmultinom(100, 
#'                                size = 100, prob = c(0.8,0.2,0.01)))),
#'                                features = data.frame(temperature = 10)))
#' params <- list(alpha=1,beta=1,disp_fun="H")
#' adult.dispersal <- dispersal(params) 
#' mp <- metapopulation(nrep=10, time=100, habitat=habitat, dispersal=adult.dispersal,
#'                   x1 = 0.42, e1 = 0.061, y1 = 1.2)
#'                   
#' library(raster)
#' set.seed(42)
#' xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
#' Dd <- as.matrix(dist(xy))
#' w <- exp(-1/nrow(xy) * Dd)
#' Ww <- chol(w)
#' xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
#' coordinates(xy) <- ~x+y
#' r <- rasterize(xy, raster(points2grid(xy)), 'z')
#' r2 <- raster(r)
#' res(r2) <- 0.01
#' r2 <- resample(r, r2)
#' proj4string(r2) <- '+init=epsg:4283'                 
#' habitat <- as.habitat(list(r2,population = as.population(t(rmultinom(1, 
#'                                size = 100, prob = c(0.8,0.2,0.1))))))
#' mp1 <- metapopulation(nrep=10, time=100, habitat=habitat, dispersal=adult.dispersal,
#'                   x1 = 0.42, e1 = 0.061, y1 = .2)
#'                   
#' habitat <- as.habitat(list(r2))
#' mp2 <- metapopulation(nrep=10, time=100, habitat=habitat, dispersal=adult.dispersal,
#'                   x1 = 0.42, e1 = 0.0061, y1 = .2)

setGeneric("metapopulation",
           function(nrep=10, time=20, habitat, dispersal, y1 = 1, x1 = 1, e1=1,
                    locations = NULL, c=1,coln_fun='H') {
             if(!is.habitat(habitat))stop('metapopulation needs a habitat')
             dist <- as.matrix(distance(habitat))
             area <- as.numeric(unlist(area(habitat)))
             locations <- as.matrix(coordinates(habitat))
             pop <- population(habitat)# need to integrate demographic and presence correctly.
             presence <- as.vector(ifelse(pop[,ncol(pop)]>0,1,0))
             dist <- dispersal$dist(habitat)
             # call c++ function that does this loop.
             mp <- dlmpr::metapop_n_cpp(nrep=nrep, time=time, dist=dist, area=area, presence=presence,
                                        y = y1, x = x1, e = e1, alpha = dispersal$alpha, beta = dispersal$beta,
                                        disp_fun = dispersal$disp_fun,locations = locations,
                                        c=c,coln_fun = coln_fun)
             sim_p_obs <- base::apply(plyr::ldply(base::lapply(mp,base::colSums),function(x)x),2,base::mean)
             sim_i_obs <- base::apply(plyr::ldply(base::lapply(base::lapply(mp,function(x)x[,-1]),function(x)base::rowSums(x)/time),function(x)x),2,base::mean)
             results <- base::list(mp = mp, sim_p_obs = sim_p_obs, sim_i_obs = sim_i_obs,
                             nrep=nrep, time = time, dist = dist, area = area, y = y1, x = x1, e = e1, 
                             dispersal=dispersal,locations = locations)
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
  graphics::plot(0:x$time, apply(base::sapply(x$mp, function(zz) base::apply(zz,2, base::sum)),1,max),
                 ylim=c(0,max(base::sapply(x$mp, function(zz) base::apply(zz,2, base::sum)))+1),
                 type = 'n', lty = 1, xlab = "time", ylab = "population size", ...)
  ci <- base::apply(base::sapply(x$mp, function(zz) base::apply(zz,2, base::sum)),1,function(x)quantile(x,c(0.025,0.975)))
  polygon(c(0:x$time,rev(0:x$time)),c(ci[1,],rev(ci[2,])),col="grey80",border=NA)
  graphics::lines(0:(x$time), x$sim_p_obs,type = 'l', col='black',lwd=2)  
}