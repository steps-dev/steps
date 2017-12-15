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
#' 
alter_adult_survival <- function(demography,survial=0.9){
                                #first check that demography is a demographic class
                                if(!is.demography(demography))stop('check that demography is a demography object')
                                  
    
 
 
 }
#' ##Define a function for manipulating habitat.
#' ##This can be a custom function for manipulating rasters or existing functions. 
#' ##Create a named list with corresponding parameters and values
#' params = list(habitat,
#'              fire_start_location = sample(ncell(habitat_suitability(habitat)),10),
#'              prob = 0.24,
#'              continue_to_burn_prob = 0.01)
#'               
#'## check habitat_suitability_dynamics produces expected output
#' fire_affected_habitat_suitability <- as.habitat_dynamics(fire_module,params,check=TRUE)                   
#'
#'## If it does? create the habitat_dynamics (yay!).
#'fire_habitat_suitability_dynamics <- as.habitat_dynamics(fun,params) 

library(raster)
library(dhmpr)
# set a demography matrix
mat <- matrix(c(.53,0,.62,0.15,0.87,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult') 
demog <- as.demography(mat)
n_stages <- length(states(demog))
print(demog)

#set up starting populations
set.seed(42)
xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
Dd <- as.matrix(dist(xy))
w <- exp(-1/nrow(xy) * Dd)
Ww <- chol(w)
xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
coordinates(xy) <- ~x+y
r <- rasterize(xy, raster(points2grid(xy)), 'z')
proj4string(r) <- '+init=epsg:4283'
r[] <- scales::rescale(r[],to=c(0,1))
random_populations <- sampleRandom(r, size=50, na.rm=TRUE, sp=TRUE)
random_populations@data <- as.data.frame(t(rmultinom(50, size = 50, prob = c(0.8,0.2,0.1))))

#set up habitat
features <- list('habitat_suitability_map'=as.habitat_suitability(r),
                 'population'=as.populations(random_populations),
                 'carrying_capacity'=as.carrying_capacity(100))
habitat <- as.habitat(features)
print(habitat)

## just playing with idea of make cell specific demography matricies which could be altered.
pop_vec <- lapply(populations(habitat),function(x)c(x[]))
pop_mat <- do.call(cbind,pop_vec)
habsuit <- habitat_suitability(habitat)

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
