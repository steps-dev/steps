#' @title demography objects 
#' @name as.demography
#' @rdname demography
#' @description demography is one of the main functions for dhmpr, it helps you construct and process stage based matricies.
#' Once a demography object is created \code{summary} will return:
#'  The finite rate of increase "lambda",
#'  the stable stage distribution,
#'  the reproductive value,
#'  and the sensitivities and elasticities matrices.
#'  \code{plot} will return a graph object that plots the stages between and amoungest each stage of the population matrix.
#' @param x For as.demography, x is a square matrix, that has stage stages between population stages.
#' @param method can be 'global' or 'local' stage matrices. 'global' is the default and a single stage-based stage matrix 
#' will be used for all populations. if 'local' method is called a cell specific stage matrix will be setup for all non-NA cells
#'  in the underlying habitat_suitability.
#' @return An object of class demography, i.e, resulting from as.demography.
#' @author Skipton Woolley
#' @export
#' @examples 
#' mat <- matrix(c(.53,0,.52,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult') 
#' demo <- as.demography(mat)

as.demography <- function(x, type='global', habsuit=NULL, stage_matrix_sd=0){
  #object <- list(...)
  #if(length(object)==0) object <- list(1)
  if(type=='local'){
    if(is.null(habsuit)){
      stop('You must include a "habitat_suitability" spatial grid if you are using method "local".')
    } 
    if(!sapply(list(habsuit),inherits,c("RasterLayer","RasterBrick","RasterStack"))){
      stop('The "habitat_suitability" input must be a spatial grid if you are using method "local".')
    }
    if(type=='local' & sapply(list(habsuit),inherits,c("RasterBrick","RasterStack"))){
      habsuit <- habsuit[[1]] 
    }
  }
  transition <- switch(type,
                       global = global_stage_matrix(x),      
                       local = local_stage_matrices(x, habsuit))

# if(length(object)>1){
#   transition[["env_stoch_matrix"]] <- object[[length(object)+1]]
# }else{
#   transition[["env_stoch_matrix"]] <- object[[length(object)]]
# }
  transition[["env_stoch_matrix"]] <- stage_matrix_sd
  
  structure(list(transition),class='demography')
  return(transition)
}

global_stage_matrix <- function(x){
  x <- base::as.matrix(x)
  if(base::diff(base::dim(x)) !=0) stop("Needs to be a square matrix with stage probabilities between each stage.")
  stage_matrixCheck(x) 
  di <- base::dim(x)[[1]]
  m.names <- base::dimnames(x)[[1]]
  if(base::is.null(m.names)) m.names <- base::paste0("stage.",1:di)
  base::dimnames(x) <- base::list(m.names, m.names)
  transition <- structure(list(global_stage_matrix=x),class='demography')
  return(transition)
}


local_stage_matrices <- function(x,habsuit){
  
  #check that the original stage-based matrix is of the right dimensions etc.
  x <- base::as.matrix(x)
  if(base::diff(base::dim(x)) !=0) stop("Needs to be a square matrix with stage probabilities between each stage.")
  stage_matrixCheck(x) 
  di <- base::dim(x)[[1]]
  m.names <- base::dimnames(x)[[1]]
  if(base::is.null(m.names)) m.names <- base::paste0("stage.",1:di)
  
  #get the number of cells that are available in the landscape/habitat_suitability object - this will be all non-NA cells.
  stopifnot(inherits(habsuit,c("RasterLayer","RasterBrick","RasterStack")))
  npops <- length(habsuit[!is.na(habsuit[])])
  
  # set up a matrix of (ncell)(rows)*(nstage*nstage)(cols) and fill with original stage matrix
  all_stage_matrices <- matrix(rep(c(x),npops),npops,length(c(x)),byrow = TRUE)
  colnames(all_stage_matrices)<-paste0(rep(m.names,each=di),1:di)
  
  transition <- structure(list(global_stage_matrix=x,local_stage_matrices=all_stage_matrices),class='demography')
  return(transition)
}

#' @rdname demography
#' @param object an object of \code{demography} class
#' @export
#' @examples
#' summary(demo)
#' 
summary.demography <- function(x,...){
    demographyCheck(x)
    x <- x$global_stage_matrix
    name.mat<-deparse(substitute(x))
    x <- x
    di <- base::dim(x)[1]
    m.names <- base::dimnames(x)[[1]] 
    ea<- base::eigen(x)
    lambda <- base::abs(ea$values[1]) 
    ssd <- base::abs(ea$vectors[,1]/base::sum(ea$vectors[,1]) ) 
    ae <- base::eigen(base::t(x))
    vr <- base::abs(ae$vectors[,1]/ae$vectors[1,1] )
    sensitivity <-  (vr  %*%  base::t(ssd))  / (base::t(vr) %*% ssd)[1,1]
    elasticity <- sensitivity * x / lambda
    
    #### Add stochasticity summary...
    result<- list(lambda=lambda, stable.stage.distribution = ssd,
                  reproductive.value =vr, sensitivity = sensitivity,
                  elasticity=elasticity,name.mat=name.mat,m.names= m.names)
    return (result)
  }

#' @rdname demography
#' @name stages
#' @param demography a demography object
#' @export
#' @examples
#' # get component stages
#' stages(demo)

stages <- function (demography) {
  # given a list of demographys, extract all of the mentioned stages
  stages <- colnames(demography$global_stage_matrix)
  return (stages)
}

#' @rdname demography
#' @importFrom igraph graph.adjacency
#' @export
#' @author Nick Golding
#' @examples 
#' plot(demo)
#' 

plot.demography <- function (x, ...) {
  # plot a dynamic using igraph
  
  # extract the stage matrix & create an igraph graph x
  graphics::par(mar=c(2,4,4,2))
  demographyCheck(x)
  x <- x$global_stage_matrix
  textmat <- base::t(x)
  textmat[textmat>0]<-base::paste0('p(',base::as.character(textmat[textmat>0]),')')
  textmat[textmat=='0'] <-''
  linkmat <- textmat != ''
  g <- igraph::graph.adjacency(linkmat, weighted = TRUE)
  
  # extract edge labels
  labels <- textmat[igraph::get.edges(g, base::seq_len(base::sum(linkmat)))]
  
  # vertex plotting details
  igraph::V(g)$color <- grDevices::grey(0.9)
  igraph::V(g)$label.color <- grDevices::grey(0.4)
  igraph::V(g)$label.family <- 'sans'
  igraph::V(g)$size <- 50
  igraph::V(g)$frame.color <- NA
  
  # edge plotting details
  igraph::E(g)$color <- grDevices::grey(0.5)
  igraph::E(g)$curved <- igraph::curve_multiple(g, 0.1)
  igraph::E(g)$arrow.size <- 0.5
  igraph::E(g)$label <- labels
  igraph::E(g)$loop.angle <- 4
  igraph::E(g)$label.color <- grDevices::grey(0.4)
  
  graphics::plot(g)
  
  # return the igraph x
  return (base::invisible(g))
  
}

#' @rdname demography
#' @name is.demography
#' @export
#' @examples
#' is.demography(demo)
is.demography <- function (x) {
  inherits(x, 'demography')
}


stage_matrixCheck <- function (x) {
  stopifnot(ncol(x) == nrow(x))
  stopifnot(is.matrix(x))
  stopifnot(all(is.finite(x)))
}

demographyCheck <- function (x) {
  stopifnot(inherits(x,'demography'))
  stopifnot(is.list(x))
}

# NOT SURE IF THIS IS USEFUL...YET.
#dots <- function(...) {
#  eval(substitute(alist(...)))
#}


#### message for Casey: I realised I messed up the demographic projections.
#### I had the matrix multipication around the wrong way.
#### I had vec_pops%*%stage_matrix, where is should have been stage_matrix%*%vec_pops.
#### I'm going to write a function which shouild do all the demographic projections 
#### and hopefully sort out any issues.
#' @rdname demography
#' @name estimate_demography
#' @param demographic_params a list of parameters which can be used to manipulate demographic projections.
#'\itemize{
#'\item{stage_matrix_sd}{Matrix with the standard deviation of the probabilities in \code{mat}.}
#'}
#' @export

#estimate_demography <- function(demography_object, habitat_object, time_step, stage_matrix_sd=0, seed=NULL){
estimate_demography <- function(demography_object, habitat_object, time_step, seed=NULL){
   
  pop_vec <- lapply(populations(habitat_object),function(x)c(x[]))
  pop_mat <- do.call(cbind,pop_vec)
  if(any(is.na(pop_mat))){
    stop(paste0("\nNA values detected in iteration ", time_step))
    #cat(paste0("\nNA values detected in iteration ", time_step))
  }
  
  if(inherits(demography_object$env_stoch_matrix,'numeric')){
    stage_matrix_sd <- demography_object$global_stage_matrix
    stage_matrix_sd[] <- demography_object$env_stoch_matrix
  }else{
    stage_matrix_sd <- demography_object$env_stoch_matrix
  }
 
  ns <- length(stages(demography_object))
  # if(!all(dim(demography_object$global_stage_matrix)==ns)){
  #   stop("Check dimensions of stage matrix.")
  # }
  if('local_stage_matrices' %in% attributes(demography_object)$names){
    if(time_step==1){cat("\nUsing a local demographic stage matricies; one per cell\n")}
    pop_mat_new <- t(sapply(1:nrow(pop_mat),function(x) estdemo(pop_mat[x,],matrix(demography_object$local_stage_matrices[x,],ns,ns), stage_matrix_sd, seed)))

  }else{
    if(time_step==1){cat("\nUsing a global demographic stage matrix for all cells\n")}
    pop_mat_new <- t(sapply(1:nrow(pop_mat),function(x)estdemo(pop_mat[x,],demography_object$global_stage_matrix, stage_matrix_sd, seed))) 
  }

  r <- populations(habitat_object)[[1]]
  
  # add function for ceiling density dependence - based on total individuals in all stages  
  # if(inherits(carrying_capacity(habitat_object),c("RasterStack","RasterBrick"))){
  #   if(sum(unlist(lapply(populations(habitat_object), function(x) cellStats(x,sum)))) > cellStats(carrying_capacity(habitat_object)[[time_step]], sum)){
  #     rescale <- cellStats(carrying_capacity(habitat_object)[[time_step]], sum)/sum(unlist(lapply(populations(habitat_object), function(x) cellStats(x,sum))))
  #     pop_mat_dd <- pop_mat_new*rescale
  #     pops_updated <- lapply(split(pop_mat_dd, rep(1:ncol(pop_mat_dd), each = nrow(pop_mat_dd))),function(x){r[]<-x;return(r)})
  #   }else{
  #     pops_updated <- lapply(split(pop_mat_new, rep(1:ncol(pop_mat_new), each = nrow(pop_mat_new))),function(x){r[]<-x;return(r)})
  #   }
  # }
  
  if(inherits(carrying_capacity(habitat_object),c("RasterLayer"))){
    if(sum(unlist(lapply(populations(habitat_object), function(x) cellStats(x,sum)))) > cellStats(carrying_capacity(habitat_object), sum)){
      rescale <- cellStats(carrying_capacity(habitat_object), sum)/sum(unlist(lapply(populations(habitat_object), function(x) cellStats(x,sum))))
      pop_mat_dd <- pop_mat_new*rescale
      pops_updated <- lapply(split(pop_mat_dd, rep(1:ncol(pop_mat_dd), each = nrow(pop_mat_dd))),function(x){r[]<-x;return(r)})
    }else{
      pops_updated <- lapply(split(pop_mat_new, rep(1:ncol(pop_mat_new), each = nrow(pop_mat_new))),function(x){r[]<-x;return(r)})
    }    
  }

  pops <- lapply(pops_updated, `attr<-`, "habitat", "populations")
  return(pops)
}


# message for Casey: This function will do the internal projections for demographic processes, so we can add in more complicated function as we need. 
estdemo <- function(popvec, stage_matrix, stage_matrix_sd=0, seed=NULL){
  
    # no stochastisity
    #popvec_new <- stage_matrix%*%popvec
    #return(popvec_new)
    
    # with env stochasticity by supplying standard deviation for random normal draw (truncated)
    if(!is.null(seed)) set.seed(seed)
    
    # if(!inherits(stage_matrix_sd,"matrix")){
    #   popvec_es <- structure(sapply(stage_matrix, function(x) if(x!=0){pmax(rnorm(1,x,stage_matrix_sd),0)}else{0}), dim=dim(stage_matrix))%*%popvec
    # }else{
      popvec_es <- structure(sapply(stage_matrix, function(x) pmax(rnorm(1,x,stage_matrix_sd),0)), dim=dim(stage_matrix))%*%popvec
    # }
    #if(any(is.na(popvec_es))){
    #  print(i)
    #}
    # add in demographic stochasticity.....
    #return(popvec_es)
    
    if(!is.null(seed)) set.seed(seed)
    
    popvec_ds <- as.matrix(
      rowSums(
        structure(
          c(rpois(length(popvec_es), stage_matrix[1, ] * popvec_es),
            apply(stage_matrix[-1,],1,function(x) rbinom(length(popvec_es), round(as.vector(popvec_es),0), x))
            ),
          dim=dim(stage_matrix))
        ),
      cols=1)

    return(popvec_ds)
}

#rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
#  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
#}

#stage_matrix %*% popvec_es

# n <- length(popvec_es)
# 
# set.seed(1)
# ans <- rep(NA, n)
# for (j in seq_len(n)) {
#   print(stage_matrix[, j]* popvec_es[j])
#   ans[j] <- sum(rpois(n, stage_matrix[, j]* popvec_es[j]))
# }
# 
# stage_matrix[, j]* popvec_es[j]
# 
# # to do the whole loop efficiently
# lambdas <- t(sweep(stage_matrix, 2, popvec_es, "*"))
# counts <- lambdas
# counts[] <- rpois(length(lambdas), lambdas)
# colSums(counts)
# 
# # to do just the first row
# counts <- rpois(n, stage_matrix[1, ] * popvec_es)
# 
# # for a subset of cells
# which_mat <- matrix(c("S", "F", "F", "F",
#                       "G", "S", "", "",
#                       "", "G", "S", "",
#                       "", "", "G", "S"),
#                     ncol = n, nrow = n, byrow = TRUE)
# 
# # find fecundities
# which_F <- which(which_mat == "F", arr.ind = TRUE)
# lambdas <- stage_matrix[which_F] * popvec_es[which_F[, 2]]
# 
# ans
# set.seed(1)
# rpois(length(popvec_es), stage_matrix %*% popvec_es)
# 
# 
# t( apply( stage_matrix , 1 , `*` , popvec_es ) )
# 
# t(apply(t(stage_matrix)*as.vector(popvec_es),1,sum))
# 
# 
# apply(stage_matrix[-1,],1,function(x) rbinom(length(popvec_es), round(as.vector(popvec_es),0), x))


#           Stage_0-1 Stage_1-2 Stage_2-3 Stage_3+
# Stage_0-1      0.00     0.000     0.302    0.302
# Stage_1-2      0.94     0.000     0.000    0.000
# Stage_2-3      0.00     0.884     0.000    0.000
# Stage_3+       0.00     0.000     0.793    0.793

#            [,1]
# [1,]   0.000000
# [2,]   8.618288
# [3,]  36.548778
# [4,] 140.850154

#                 [,1]
# Stage_0-1  53.574477
# Stage_1-2   0.000000
# Stage_2-3   7.618566
# Stage_3+  140.677353
