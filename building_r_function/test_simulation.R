## now we could run the whole thing manytimes to simulate populations with uncertainty.
library(raster)
library(dhmpr)
# set a transition matrix
mat <- matrix(c(.53,0,.62,0.20,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
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


# here I create starting populations based on habitat suitability and stable starting state x 100 for good measure.
starting_populations <- lapply(1:3,function(x)r*100*summary(trans)$stable.stage.distribution[x])
start_pops <- stack(starting_populations)

dispersal_params <- as.dispersal(list(dispersal_distance=list('larvae'=3,'juvenile'=0,'adult'=6),
                                      dispersal_kernel=list('larvae'=exp(-c(0:2)),'juvenile'=0,'adult'=exp(-c(0:5)*.2)),
                                      dispersal_proportion=list('larvae'=0.1,'juvenile'=0,'adult'=0.6)))  
print(dispersal_params)


#set up habitat
simulations <- 10
sim_results <- list()
for(j in 1:simulations){
  features <- list('habitat_suitability_map'=as.habitat_suitability(r),
    'population'=as.populations(start_pops),
    'carrying_capacity'=as.carrying_capacity(100))
  habitat <- as.habitat(features)
  # print(habitat)
  
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
  
  ## altenative you can run the fire module as follows:
  RUN_FIRE_screaming <- fire_module(habitat,sample(ncell(habitat_suitability(habitat)),10),
    prob = 0.24,
    continue_to_burn_prob = 0.01)
  
  # set up a density dependence fucntion. This isn't working as before :( But it should give you some ideas.
  ddfun <- function (pop) {
    adult_density <- pop 
    .9 * exp(-adult_density/80)* pop
  }
  
  #now we are going to try and set up an experiment using the above functions.
  pops_at_time_step <- dispersal_at_time_step <- fire_at_time_step <- habitat_suit_at_time_step <- list()
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
    pops_n <- pop_mat %*% (trans$stage_matrix * matrix(runif(length(trans$stage_matrix),max=2.5),dim(trans$stage_matrix)[1],dim(trans$stage_matrix)[2])) 
    
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
      continue_to_burn_prob = 0.01)
  }
  sim_results[[j]] <- pops_at_time_step
  cat(j,"\n")
}

# spatial
mean_lar <- all_lar <- list()
mean_juv <- all_juv <- list()
mean_adl <- all_adl <- list()

for(i in 1:12){
  all_lar <- list()
  all_juv <- list()
  all_adl <- list()
  
   for(j in 1:10){
      all_lar[[j]] <- sim_results[[j]][[i]][[1]]
      all_juv[[j]] <- sim_results[[j]][[i]][[2]]
      all_adl[[j]] <- sim_results[[j]][[i]][[3]]
      }
      mean_lar[[i]] <- mean(stack(all_lar))
      mean_juv[[i]] <- mean(stack(all_juv))
      mean_adl[[i]] <- mean(stack(all_adl))
}

plot(stack(mean_lar))
plot(stack(mean_juv))
plot(stack(mean_adl))



#through time
par(mfrow=c(1,3))
larvaes_ns <- matrix(20*12,20,12)
juveniles_ns <- matrix(20*12,20,12)
adults_ns <- matrix(20*12,20,12)
for(i in 1:12){ 
  larvaes_ns[,i] <- sapply(lapply(lapply(sim_results, "[[", i),"[[",1), function(x)mean(x[]))
  juveniles_ns[,i] <- sapply(lapply(lapply(sim_results, "[[", i),"[[",2), function(x)mean(x[]))
  adults_ns[,i] <- sapply(lapply(lapply(sim_results, "[[", i),"[[",3), function(x)mean(x[]))
}

ln <- (apply(larvaes_ns, 2, quantile, c(0.1, 0.90)))
jn <- (apply(juveniles_ns, 2, quantile, c(0.1, 0.90)))
an <- (apply(adults_ns, 2, quantile, c(0.1, 0.90)))

lnm <- (apply(larvaes_ns, 2, mean))
jnm <- (apply(juveniles_ns, 2, mean))
anm <- (apply(adults_ns, 2, mean))


par(mfrow=c(1,3))
x <- 1:12
y <- 100
plot(x,seq(0,100,length.out = 12),type = 'n', ylab = 'population',xlab = 'time (years)',axes=F)
polygon(x = c(x, rev(x)),
  y = c(ln[1, ], rev(ln[2,])),
  col = grey(0.8),
  border = NA)
lines(1:12,lnm, lwd = 1,type="b",pch=16)
text(max(x),y-10,"Mean larvae population\n at each time step \n with variance", adj=1, family="serif")
axis(1, at=x, label=x, family="serif")
axis(2, at=seq(0,100,10),label=seq(0,100,10), family="serif",las=2)

x <- 1:12
y <- 100
plot(x,seq(0,100,length.out = 12),type = 'n', ylab = 'population',xlab = 'time (years)',axes=F)
polygon(x = c(x, rev(x)),
  y = c(jn[1, ], rev(jn[2,])),
  col = grey(0.8),
  border = NA)
lines(1:12,jnm, lwd = 1,type="b",pch=16)
text(max(x),y-10,"Mean juvenile population\n at each time step \n with variance", adj=1, family="serif")
axis(1, at=x, label=x, family="serif")
axis(2, at=seq(0,100,10),label=seq(0,100,10), family="serif",las=2)

x <- 1:12
y <- 100
plot(x,seq(0,100,length.out = 12),type = 'n', ylab = 'population',xlab = 'time (years)',axes=F)
polygon(x = c(x, rev(x)),
  y = c(an[1, ], rev(an[2,])),
  col = grey(0.8),
  border = NA)
lines(1:12,anm, lwd = 1,type="b",pch=16)
text(max(x),y-10,"Mean adult population\n at each time step \n with variance", adj=1, family="serif")
axis(1, at=x, label=x, family="serif")
axis(2, at=seq(0,100,10),label=seq(0,100,10), family="serif",las=2)


