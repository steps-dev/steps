library(raster)
library(dhmpr)
# set a transition matrix
mat <- matrix(c(.53,0,.52,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult') 
trans <- as.transition(mat)
n_stages <- length(states(trans))

#set up starting populations
random_populations <- sampleRandom(r, size=50, na.rm=TRUE, sp=TRUE)
random_populations@data <- as.data.frame(t(rmultinom(50, size = 100, prob = c(0.8,0.2,0.1))))

#set up habitat
features <- list('habitat_suitability_map'=as.habitat_suitability(r),
                 'population'=as.populations(random_populations),
                 'carrying_capacity'=as.carrying_capacity(100))
habitat <- as.habitat(features)

# set up dispersal parameters.
dispersal_params <- as.dispersal(list(dispersal_distance=list('larvae'=3,'juvenile'=0,'adult'=3),
                dispersal_kernel=list('larvae'=exp(-c(0:2)),'juvenile'=0,'adult'=exp(-c(0:2)*.2)),
                dispersal_proportion=list('larvae'=0.1,'juvenile'=0,'adult'=0.3)))  

## dispersal using cellular automata.                                 
dispersed_populations <- dispersal(dispersal_params, habitat, method='ca')
populations(habitat) <- dispersed_populations

## now I've set up a custom function which manipulates rasters (based on a dodgy fire simualtion)
##This can be a custom function for manipulating rasters or existing functions. 
fun <- fire_spread

##Create a named list with corresponding parameters and values
params = list(habitat_suitability(habitat),
              fire_start_location = sample(ncell(hs),10),
              prob = 0.24,
              continue_to_burn_prob = 0.01)
               
## check module produces expected output
fire_module <- as.module(fun,params,check=TRUE)       

# set up a density dependence fucntion. This isn't working as before :(
ddfun <- function (pop) {
  adult_density <- pop 
  .9 * exp(-adult_density/80)* pop
}

# apply dispersal to the population (need to run this separately for each stage)
pops <- list(as.matrix(r1),as.matrix(r2),as.matrix(r3))
pops_time<-list('t1'=pops)
n_time_steps <- 10
for(j in 1:n_time_steps){
  pops_new <- list()
    # iterate through life histories
    for(i in 1:n_stages){
      pops_new[[i]] <- dispersalFFT(popmat = pops_time[[j]][[i]], fs = fs) # add in stage specific despersal function. 
      }
  pop_vec <- lapply(pops_new,c)
  pop_mat <- do.call(cbind,pop_vec)
  
  #update populations
  pops_n <- pop_mat %*% (trans$stage_matrix)
  
  ## density dependence bit working. 
  pops_n[,3] <- ddfun(pops_n[,3])
  pops_time[[j+1]]<-lapply(split(pops_n, rep(1:ncol(pops_n), each = nrow(pops_n))),function(x)matrix(x,nrow = nrow(pops_new[[1]])))
  cat(j,'\n')
}
  
N <- lapply(pops_time,sum)

larvae <- lapply(pops_time, "[[", 1)
juve <- lapply(pops_time, "[[", 2)
adult <- lapply(pops_time, "[[", 3)

N <- lapply(adult,mean)


larv <- lapply(larvae, function(x) raster(matrix(x,nrow = nrow(r2))))
st.l <- stack(larv)

juv <- lapply(juve, function(x) raster(matrix(x,nrow = nrow(r2))))
st.j <- stack(juv)

adl <- lapply(adult, function(x) raster(matrix(x,nrow = nrow(r2))))
st.a <- stack(adl)
