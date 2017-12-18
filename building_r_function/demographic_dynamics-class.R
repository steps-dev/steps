#' @title demography_dynamics objects
#' @name demography_dynamics
#' @rdname demography_dynamics
#' @description demography_dynamics are functions for altering the underlying demographic process in the \link[code]{experiment}. 
#' Demography_dynamics are functions which either directly affect the stage-based demographic matrix or other demographic processes excluding dispersal, any functions which alter dispersal are called from \link[dhmpr]{dispersal} or \link[dhmpr]{dispersal_dynamics}.
#' \code{demography_dynamics} sets up internal or custom functions to work with \code{demography} and \code{experiment} objects.
#'  
#' @export
#' @examples 
#' ## Create population
#' library(raster)
#' library(dhmpr)
#' #set a demography matrix
#' mat <- matrix(c(.53,0,.62,0.15,0.87,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult') 
#' demog <- as.demography(mat)
#' n_stages <- length(states(demog))

alter_adult_survival <- function(demography, stage='adult', survial=0.9){
                                #first check that demography is a demographic class
                                if(!is.demography(demography))stop('check that demography is a demography object')
                                if(length(which(stages(demography)==stage))==0)stop('check names of the stage and make sure you have the right name')
                                idx <-which(stages(demography)==stage)
                                demography$stage_matrix[idx,idx]
                                  
 }


library(raster)
library(dhmpr)
# set a demography matrix
mat <- matrix(c(.53,0,.62,0.15,0.87,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult') 
demog <- as.demography(mat)
n_stages <- length(states(demog))
print(demog)

survival <- calc(habsuit,function(x){x[x<0.5] <- runif(1,0,.5) ;return(x)})
c(adult_survival[])

all_cells_stage_matricies <- matrix(rep(c(demog$stage_matrix),npops),npops,length(c(demog$stage_matrix)),byrow = TRUE)

stages_all <- 1:9
stages_lar <- 1:3
stages_juv <- 4:6
stages_adl <- 7:9

# survival prob affects which stage?
# let's try with adult
all_cells_stage_matricies[,stages_adl]<-c(survival[])*all_cells_stage_matricies[,stages_adl]

pop_mat_new <- t(sapply(1:npops,function(i)pop_mat[i,]%*%matrix(all_cells_stage_matricies[i,],dim(demog$stage_matrix)[1],dim(demog$stage_matrix)[2])))

r <- habitat_suitability(habitat)
pops_updated <- lapply(split(pop_mat_new, rep(1:ncol(pop_mat_new), each = nrow(pop_mat_new))),function(x){r[]<-x;return(r)})

# looks like adult pops are going down in less suitable areas - maholla.
plot(stack(c(pops_updated,habsuit)))
