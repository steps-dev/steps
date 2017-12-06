library(raster)
library(dhmpr)
# set a transition matrix
mat <- matrix(c(.53,0,.72,0.25,0.87,0,0,0.32,0.9),nrow = 3,ncol = 3,byrow = TRUE)
colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult') 
trans <- as.transition(mat)
n_stages <- length(states(trans))
print(trans)

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
random_populations@data <- as.data.frame(t(rmultinom(50, size = 100, prob = c(0.8,0.2,0.1))))

#set up habitat
features <- list('habitat_suitability_map'=as.habitat_suitability(r),
                 'population'=as.populations(random_populations),
                 'carrying_capacity'=as.carrying_capacity(100))
habitat <- as.habitat(features)
print(habitat)

# set up dispersal parameters.
dispersal_params <- as.dispersal(list(dispersal_distance=list('larvae'=3,'juvenile'=0,'adult'=6),
                dispersal_kernel=list('larvae'=exp(-c(0:2)),'juvenile'=0,'adult'=exp(-c(0:5)*.2)),
                dispersal_proportion=list('larvae'=0.1,'juvenile'=0,'adult'=0.6)))  
print(dispersal_params)

## dispersal using cellular automata.                                 
dispersed_populations <- dispersal(dispersal_params, habitat, method='ca')
populations(habitat) <- dispersed_populations

## now I've set up a custom function which manipulates rasters (based on a dodgy fire simualtion)
## This will be used as a module to manipulate the habitat suitability between time steps. 
## create a named list with corresponding parameters and values
params = list(habitat,
              fire_start_location = sample(ncell(habitat_suitability(habitat)),10),
              prob = 0.24,
              continue_to_burn_prob = 0.01)
               
## check module produces expected output
fires <- as.module(fire_module,params)       

## altenative you can run the fire module as follows:
RUN_FIRE_screaming <- fire_module(habitat,sample(ncell(habitat_suitability(habitat)),10),
                                  prob = 0.24,
                                  continue_to_burn_prob = 0.01)

# set up a density dependence fucntion. This isn't working as before :( But it should give you some ideas.
ddfun <- function (pop) {
  adult_density <- pop 
  if(adult_density>10) adult_density <- 10 * runif(1)
  adult_density
}

#now we are going to try and set up an experiment using the above functions.
pops_at_time_step <- dispersal_at_time_step <- fire_history_df_list <- fire_at_time_step <- habitat_suit_at_time_step <- list()
pops_at_time_step[[1]] <- populations(habitat) 
habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat) 
fire_at_time_step[[1]] <- fire_module(habitat,sample(ncell(habitat_suitability(habitat)),10),
  prob = 0.24,
  continue_to_burn_prob = 0.01)
n_time_steps <- 11

## begin experiment through time
for(i in 1:n_time_steps){
  
    #update the landscape based on last fire. 
    fire_prob <- fire_at_time_step[[i]]/max(fire_at_time_step[[i]][])   
    
    #now you could have some equation on how fire effects underlying habitat suitability. I'm just going to make     # something up :) = habitat_suitability*(1-fire_prob) 
    new_hs <- habitat_suitability(habitat)*(1-fire_prob)
    attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheridence.
    habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs 
    
    # let's disperse people
    dispersed_populations <- dispersal(dispersal_params, habitat, method='ca')
    dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
    
    # get a vector of populations per life history
    pop_vec <- lapply(populations(habitat),function(x)c(x[]))
    pop_mat <- do.call(cbind,pop_vec)
  
    # update populations (this could be done much more nicely with pop)
    pops_n <- pop_mat %*% (trans$stage_matrix)# * matrix(runif(length(trans$stage_matrix)),dim(trans$stage_matrix)[1],dim(trans$stage_matrix)[2]))
    
    ## update density dependence for adult populations. 
    pops_n[,3] <- ddfun(pops_n[,3])
    
    #now update the populations - ideally this could be turned into a function (update_pops())
    r <- habitat_suitability(habitat)
    pops_updated <- lapply(split(pops_n, rep(1:ncol(pops_n), each = nrow(pops_n))),function(x){r[]<-x;return(r)})
    pops <- lapply(pops_updated, `attr<-`, "habitat", "populations")
    pops_at_time_step[[i+1]] <- populations(habitat) <- pops
    
    #now let's start another fire! Mwhahahhaha.
    fire_at_time_step[[i+1]] <- fire_module(habitat,sample(ncell(habitat_suitability(habitat)),10),
      prob = 0.24,
      continue_to_burn_prob = 0.01,
      fire_history = )
}

#lets look at populations
#spatially
larvae <- stack(lapply(pops_at_time_step, "[[", 1))
juve <- stack(lapply(pops_at_time_step, "[[", 2))
adult <- stack(lapply(pops_at_time_step, "[[", 3))

pops_after_experiment <- stack(larvae,juve,adult)
names(pops_after_experiment) <- c('larvae','juveniles','adults')
# par(mfrow=c(3,1))
plot(pops_after_experiment,nr=1)


tot_abn <- list()
for(i in 1:12){
  tot_abn[[i]] <- sum(stack(pops_at_time_step[[i]]))
}

rb <- brick(tot_abn)
library(viridis)
animation::saveGIF(animate(rb,pause=.001,main='dispersal',n=1,col=rev(viridis::viridis(100))), movie.name = "pop_growth_t24.gif",loop=TRUE)

#through time
par(mfrow=c(1,3))

larvae_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
plot(larvae_n,type='l',ylab='Total Larvae Abundance',xlab="Time (years)",lwd=2,col='springgreen')

juve_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')

adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')


#let's look at fire through time
#spatially
fire_history <- stack(fire_at_time_step)
plot(fire_history)
