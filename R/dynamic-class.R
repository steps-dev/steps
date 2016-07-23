#' @title dynamic object
#' @rdname dynamic 
#' @name 
#' @description Very similar to \link[pop]{dynamic}, \code{dynamic} is an object which stores all the relevent 
#' stage-matrices, transitions, population(s), habitat and dispersal for using the \link[dhmpr]{simulation}.
#' The main difference between \link[pop]{dynamic} and \code{dynamic} is that dhmpr requires the input of set objects to run.

dynamic <- function(transition,population,...){

  objects <- dots(...)
  stopifnot(is.transition(transition))
  stopifnot(is.population(population))
  stopifnot(is.habitat(habitat))
  stopifnot(is.dispersal(dispersal))
  atts_t <- attributes(transition)
  atts_p <- attributes(population)
  atts_h <- attributes(habitat)
  atts_d <- attributes(dispersal)
  
  
  return(objects)
}

#' @rdname dynamic
#' @export
is.dynamic <- function (x) inherits(x, 'dynamic')

as.dynamic <- function (x) {
  if (!is.dynamic(x)) {
    class(x) <- c('dynamic', class(x))
  }
  return (x)
}

# params <- list(alpha=list('larvae'=2,'juvenile'=0,'adult'=3),
#                probability=list('larvae'=0.2,'juvenile'=0,'adult'=0.6))
# 
# disps  <- dispersal(params)
# 
# disp_probs <- probdisp(disps,habitat)
# 
# stage_disp_mat <- stage_disperal_matrix(disp_probs,p,states,M_sub)
# 
# dynamic_mat <- Matrix::t(P) %*% stage_disp_mat %*% P %*% stage_demog_mat %*%  c(rmultinom(10,200,c(.8,.3,.1))) #intial_populations

as.matrix.dynamic <- function(dynamic,...){
  
  # if patches = 1
  p <- patches(habitat)
  
  #get the stage demography matrix
  stage_demog_mat <- stage_demographic_matrix(dynamic)
  
  if(p==1){
    dynamic_mat <- stage_demog_mat
  } else {
  
  stage_disp_mat <- stage_disperal_matrix(dynamic)    
  P <- make.P(dynamic)
  
  #dynamic matrix
  dynamic_mat <-  Matrix::t(P) %*% stage_disp_mat %*% P %*% stage_demog_mat
  }
  
  return(dynamic_mat)
}

stage_demographic_matrix <- function(dynamic){
  # this function sets up the stage demographic matrix for the vec-matrix metapopulation matrix.
  # it takes a dynamic object which contains patch, and stage matrix.
  
  #get transition object
  trans <- transition(dynamic)
  
  #get the stage matrix
  mat <- trans$stage_matrix
  
  #get the habitat from dynamic object
  habitat <- habitat(dynamic)
  
  #get the number of patches
  p <- patches(habitat)
  
  #set up the kronecker mat as a sparse matrix
  stage_demog_mat <- Matrix::kronecker(Matrix::Diagonal(p),mat)
  
  return(stage_demog_mat)
}


stage_disperal_matrix <- function(dynamic){#disp_probs,states,M_sub){
  
  #get transition object
  trans <- transition(dynamic)
  
  #get the habitat from dynamic object
  habitat <- habitat(dynamic)
  
  # get the number of patches
  p <- patches(habitat)
  
  # get the states
  states <- states(tran)
  
  # get the number of states
  s <- length(states)
  
  # get dispersal 
  disp_probs <- probdisp(disps,habitat)
  
  # get the states that have dispersal matricies
  disp_states <- rep(0,s)
  disp_states[match(names(disp_probs),states)] <- match(names(disp_probs),states)
  
  # set up the sparse matrix
  dim_mat <-s*s*p^2
  stage_disp_mat <- Matrix::sparseMatrix(dims = c(1,dim_mat), i={}, j={},x=rep(0.0,dim_mat))
  
  #get the M-sub mat as a vec to save memory.
  dims<-p^2
  index<-0:(p-1)
  Index<-1+index*p+index
  M_sub<-sparseMatrix(i = rep(1,p), j = Index, x= rep(1.0,p), dims=c(1,dims))
  
  #Calculate the indices which have to be added for all p^2 indices of M_ii
  indx1 <- seq(1, p * p);
  xp <- (indx1 - 1) - ((indx1 -1) %/% p) * p
  yp <- (indx1 - 1) %/% p;
  indx2 <- xp + yp * s * p
  
  # add data to dispersal matrix
  z <- 1
  m <- 1
  for (i in seq_along(states)){
    if(i==disp_states[i]){   #check which stages have dispersal
      stage_disp_mat[z+indx2] <- disp_probs[[m]][indx1]
      m <- m + 1
    } else {
      stage_disp_mat[z+indx2] <- M_sub[indx1] # if no dispersal make it the M_sub mat.
    }
    z <- 1+(1+i-1)*(s*p^2+p)
  }
  dim(stage_disp_mat)<-rep(s*p,2)
  
  #check the resulting matrix
  # image(stage_disp_mat)
  
  return(stage_disp_mat)
  
}

make.P <- function(dynamic){
  
  # get patches and stages from dynamic object
  habitat <- habita(dynamic)
  p <- patches(habitat) #(or/cells)
  s <- stages(dynamic) #stages
  
  rw<-cl<-rep(0, s*p)
  
  index<-0
  for( i in 1:s){
    for( j in 1:p){
      index<-index+1
      cl[index]<- (j-1)* s + i
      rw[index]<- (i-1)* p + j
    }
  }
  
  P <- Matrix::sparseMatrix(i = rw, j = cl, x= rep(1,length(rw)),dims=rep(s*p,2))
  return(P)
}
