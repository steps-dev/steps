#' calculate carrying capacity from occupancy
#' 
#' @param x a raster of species habitit suitability (occupancy)
#' @param a=6
#' @param b=3
#' @param c=0.5
#' @param cc_mod model form for converting hs to carrying capacity. 
#' @author Skipton Woolley

cc_fun<-function(x,a=6,b=3,c=0.5,cc_mod=c('exp','logit','linear')){
  type <- match.arg(cc_mod)
  switch(type,
         exp = exp((a*x)-b),
         linear = a*(x)-b,
         logit = a+(1/(1+exp(-b*x+c))))
}
  
