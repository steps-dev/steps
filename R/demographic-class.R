#' @useDynLib dhmpr
#' @importFrom Rcpp sourceCpp
NULL

#' @title demographic objects
#' @name demographic
#' @rdname demographic
#' @param pop vector Starting abundances at each timestep.
#' @param tmat matrix. Transition matrix
#' @param matsd matrix. Matrix with the standard deviation of the probabilities in tmat. 
#' @param stmat martix. Matrix indicating for each transition probability in mat which part (i.e. which proportion) should be considered resulting from fecundity.
#' @param estamb Logical. Should environmental stochasticity be considered to projet the dynamics of the population?
#' @param estdem Logical. Should demographic stochasticity be employed to project the dynamics of the population?
#' @param equalsign Logical. Should the environmental deviations have all the same sign and magnitude? 
#' @param tmat_fecundity Logical. Should you use the first row of tmat as fecundities?
#' @param nrep int number of simulations to run.
#' @param time int length of the demographic trajectory.
#' @return demographic_mod object a list of demographic projections
#' @author Skipton Woolley
#' @seealso \code{link{demo_proj_n_cpp}}
#' @description do some demographic modelling in R
#' 
#' @examples 
#' tmat <- as.transition(matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),
#' nrow = 3,ncol = 3,byrow = TRUE))
#' matsd <- tmat$stage_matrix/10
#' pop <- as.population(c(80,20,0))
#' proj <- demographic(pop=pop,tmat=tmat,matsd=matsd,time=100) 
#' @export
#' 
setGeneric("demographic",
  function (pop, tmat, matsd = NULL, stmat = NULL,
            estamb = FALSE, estdem = FALSE, 
            equalsign = TRUE, tmat_fecundity = FALSE, nrep = 10, 
            time = 10) {
  if(!is.transition(tmat)) stop("The transition matrix is not of transition class")
    tmat <- tmat$stage_matrix
    vn <- NULL
      for (i in 1:nrep) {
          vn[[i]] <- base::cbind(pop, pop)
      }
  # call c++ function that does this loop.
  v <- dhmpr::demo_proj_n_cpp(vn, tmat, matsd = matsd, estamb = estamb, estdem = estdem, 
                       equalsign = equalsign, stmat = stmat, tmat_fecundity = tmat_fecundity,
                       nrep = nrep, time = time)
  vn <- base::lapply(v, function(x) x[,-1])
  di <- dim(tmat)[1]
  m.names <- dimnames(tmat)[[1]]
  if(base::is.null(m.names)) m.names <- base::paste0("stage.",1:di)
  ea<- base::eigen(tmat)
  lambda <-base::abs(ea$values[1])
  ssd <- base::abs(ea$vectors[,1]/sum(ea$vectors[,1]) ) 
  ae <- base::eigen(base::t(tmat))
  vr <- base::abs(ae$vectors[,1]/ae$vectors[1,1] )
  sensitivity <- (vr %*% t(ssd))/(base::t(vr) %*% ssd)[1,1]
  elasticity <- sensitivity*tmat/lambda
  vn <- base::list(vn = vn, tmat = tmat,lambda=lambda, stable.stage.distribution = ssd,
                   reproductive.value =vr, sensitivity = sensitivity,
                   elasticity=elasticity,m.names= m.names)
  class(vn) <- "demographic"
  return(vn)
}
)

#' @rdname demographic
#' @method plot demographic
#' @param x demographic model object
#' @param mean_pop logical. If TRUE plots the mean demographic change in population through time. Otherwise it plots the population changes for each stage in the transition matrix.
#' @param ... other plot arguments
#' @author Skipton Woolley
#' @examples 
#' tmat <- as.transition(matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),
#' nrow = 3,ncol = 3,byrow = TRUE))
#' matsd <- tmat$stage_matrix/10
#' pop <- as.population(c(80,20,0))
#' proj <- demographic(pop=pop,tmat=tmat,matsd=matsd,estdem=TRUE,time=20,nrep=100) 
#' plot(proj,mean_pop=TRUE)
#' plot(proj,mean_pop=FALSE)
#' @export
plot.demographic <- function(x,mean_pop=TRUE,...){
  x <- x$vn
  nrep <- base::length(x)
  stages <- base::dim(x[[1]])[1]
  time <- base::dim(x[[1]])[2]
  if(mean_pop==TRUE){
    graphics::par(mfrow = c(1,1))
    graphics::matplot(0:(time - 1), (base::sapply(x, function(re) base::apply(re,2, base::sum))),
                      type = 'n', xlab = "time", ylab = "abundance",pch = 1,col="#00000030", ...)
    ci <- base::apply(base::sapply(x, function(re) base::apply(re,2, base::sum)),1,function(x)quantile(x,c(0.025,0.975)))
    polygon(c(0:(time-1),rev(0:(time-1))),c(ci[1,],rev(ci[2,])),col="grey80",border=NA)
    graphics::lines(0:(time - 1), base::apply(base::sapply(x, function(re) base::apply(re,2, sum)), 1, base::mean),
                    type = 'l', col='black',lwd=2)  
  } else {
    graphics::par(mfrow = c(1,stages))
    for (i in 1:stages) {
      graphics::matplot(0:(time - 1), base::sapply(x,function(st) st[i, ]), 
                                           xlab = "time", ylab = "abundance", 
                                           type = 'n', col="#00000030", pch = 1, 
                                           main = base::paste("stage", i),...)
      ci <- base::apply( base::sapply(x,function(st) st[i, ]),1,function(x)quantile(x,c(0.025,0.975)))
      polygon(c(0:(time-1),rev(0:(time-1))),c(ci[1,],rev(ci[2,])),col="grey80",border=NA)
      
      graphics::lines(0:(time - 1), base::apply(base::sapply(x,function(st) st[i, ]),1,mean),
                      type = 'l', col=base::suppressWarnings(RColorBrewer::brewer.pal(base::length(base::unique(stages)),
                                                                                      "Set1"))[i],lwd=2)
    }
  }
}