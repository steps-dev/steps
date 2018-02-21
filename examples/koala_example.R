#library(dhmpr)
library(raster)
library(foreach)

#system('rm /home/casey/Research/Github/dhmpr/src/*.so /home/casey/Research/Github/dhmpr/src/*.o')

# Questions for Nat/Skip:
# 1) Best way to input populations - as multiple (based on life stages) or single components? - improve warning message for now
# 2) Where should density dependence be incorporated? In the carrying capacity attributes? function that uses the carrying capacity and populations in age classes (structure) that contribute to check -> adjusts survival/fecundity in demo or population in habitat
# 3) How to incorporate stochasticity? In the habitat or demographic dynamics objects?  look into past code from skip - inquire about demographic stochasticity
# 4) Reset habitat layers after disturbance based on years of effect

# develop all data inputs:

.pardefault <- par()
par(mar=c(0.5,0.5,0.5,0.5))

#Note, this is an age-structured matrix and is not generic
transition_matrix <- matrix(c(0.000,0.000,0.302,0.302,
                              0.940,0.000,0.000,0.000,
                              0.000,0.884,0.000,0.000,
                              0.000,0.000,0.793,0.793),
                            nrow = 4, ncol = 4, byrow = TRUE)
colnames(transition_matrix) <- rownames(transition_matrix) <- c('Stage_0-1','Stage_1-2','Stage_2-3','Stage_3+')

transition_matrix_es <- matrix(c(0.000,0.000,1,1,
                                 1,0.000,0.000,0.000,
                                 0.000,1,0.000,0.000,
                                 0.000,0.000,1,1),
                               nrow = 4, ncol = 4, byrow = TRUE)
colnames(transition_matrix_es) <- rownames(transition_matrix) <- c('Stage_0-1','Stage_1-2','Stage_2-3','Stage_3+')

############ MODEL INPUTS ###############

#intended timesteps
n <- 20

#### 1 ####
koala.hab.suit <- raster("data/Koala/Habitat/HS_crop_aggregate.tif") # read in spatial habitat suitability raster
plot(koala.hab.suit, box = FALSE, axes = FALSE)

#### 2 ####
vec.suit <- c(seq(.9,.3,length.out = n/2),seq(.3,.9,length.out = n/2))
koala.hab.suit.s <- stack(unlist(foreach(i=1:20) %do% {koala.hab.suit * vec.suit[i]}))

#### 3 ####
koala.hab.k <- koala.hab.suit*0.866
plot(koala.hab.k, box = FALSE, axes = FALSE)

#### 4 #### - No longer required
#vec.k <- c(seq(1,0.1,length.out = n/2),seq(1,0.1,length.out = n/2))
#koala.hab.k.s <- stack(unlist(foreach(i=1:20) %do% {koala.hab.k * vec.k[i]}))
#koala.hab.k.s <- stack(koala.hab.k,koala.hab.k,koala.hab.k,koala.hab.k,koala.hab.k,
#                       koala.hab.k,koala.hab.k,koala.hab.k,koala.hab.k,koala.hab.k,
#                       koala.hab.k,koala.hab.k,koala.hab.k,koala.hab.k,koala.hab.k,
#                       koala.hab.k,koala.hab.k,koala.hab.k,koala.hab.k,koala.hab.k)

#### 5 #### - No longer required
#koala.hab.k.n <- 20

#### 6 ####
koala.hab.k.func <- function(x) x*0.866

#### 7 ####
koala.hab.pop <- (koala.hab.suit/0.866)*0.03 # create new raster to hold initial abundance values
plot(koala.hab.pop, box = FALSE, axes = FALSE)

#### 8 ####
vec.pop <- c(seq(.9,.7,length.out = n/2),seq(.9,.8,length.out = n/2))
koala.hab.pop.s <- stack(unlist(foreach(i=1:20) %do% {koala.hab.pop * vec.pop[i]}))

#### 5 ####
koala.hab.pop.n <- 1

#### 10 ####
koala.demo.glob <- as.demography(transition_matrix, type="global", transition_matrix_es)
#n_stages <- length(stages(koala.demo.glob))
print(koala.demo.glob)
summary(koala.demo.glob)

#### 11 ####
koala.demo.loc <-as.demography(transition_matrix, type="local", koala.hab.suit, transition_matrix_es)
#n_stages <- length(stages(koala.demo.loc))
print(koala.demo.loc)
summary(koala.demo.loc)

#### 12 ####
koala.disp.bar <- koala.hab.suit*0
koala.disp.bar[cellFromRow(koala.disp.bar,nrow(koala.disp.bar)/2)] <- 1

#### 13 ####
koala.disp.bar.s <- stack(mget(rep("koala.disp.bar",n)))

#### 14 & 15 ####
koala.disp.param <- as.dispersal(list(dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
                                      dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
                                      dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0)
                                      )
                                 )  
print(koala.disp.param)

#### 16 ####
koala.dist.fire <- stack(list.files("data/Koala/Fire", full = TRUE, pattern = '*agg'))[[18]]

#### 17 ####
koala.dist.fire.s <- stack(rep(list.files("data/Koala/Fire", full = TRUE, pattern = '*agg'), length.out = n))

#### 18 ####
# koala.dist.fire.func.ran <- function (x, n) {
#   x[cellFromRowCol(x,sample(c(1:nrow(x)),n),sample(c(1:ncol(x)),n))] <- 1
#   return(x)
# }

koala.dist.fire.func <- function (x,habitat) {
  return(x*habitat)
}

#### 19 ####
koala.dist.fire.func.param <- koala.dist.fire.s
koala.dist.fire.func.param2 <- koala.dist.fire

n_time_steps <- 20

####### Permutation 1 ########

features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit),
                 'population'=as.populations(stack(koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[1],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[2],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[3],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[4]
                                                   )
                                             ),
                 'carrying_capacity'=as.carrying_capacity(koala.hab.k))
habitat <- as.habitat(features)
print(habitat)

demo.g <- as.demography(transition_matrix, type="global", transition_matrix_es)
demo.l <- as.demography(transition_matrix, type="local", koala.hab.suit, transition_matrix_es)

hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)

pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()

pops_at_time_step[[1]] <- populations(habitat)
habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat)

#system.time(
  pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
  for(i in 1:n_time_steps){

    if(any(i == seq(4,n_time_steps,3))){ # reset back to orginal or modified hs assuming time period of fire effect (3 years)
      r <- habitat_suit_at_time_step[[1]]
      r[] <- 1
      disturbance_at_time_step[[i]] <- r
      habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- habitat_suit_at_time_step[[1]]
    }else{ # change habitat suitability based on disturbance
      disturbance_at_time_step[[i]] <- hab_dyn[[2]][[i]]
      new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),time_step=i)
      attr(new_hs,"habitat") <- "habitat_suitability"
      habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
    }    

    dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca',time_step=i) # disperse populations in landscape
    dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
    
    pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i) # update populations based on transition matrix (incorporates stochasticity and density dependence)

    setTxtProgressBar(pb, i)
  }
  close(pb)
#) ##### approx 30 seconds to do 20 timesteps 

Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))

plot(Stage_0_1)
plot(Stage_1_2)
plot(Stage_2_3)
plot(Stage_3_)

#through time
par(.pardefault)
par(mfrow=c(1,4))

sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')

adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')

sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')

juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')

######################################

####### Permutation 2 ########

features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit.s),
                 'population'=as.populations(stack(koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[1],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[2],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[3],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[4]
                 )
                 ),
                 'carrying_capacity'=as.carrying_capacity(koala.hab.k))
habitat <- as.habitat(features)
print(habitat)

hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)

pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()

pops_at_time_step[[1]] <- populations(habitat)
habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat,1)

  pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
  for(i in 1:n_time_steps){
    
    disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
    
    new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat,time_step=i),time_step=i)
    attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
    habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat,time_step=i) <- new_hs
    
    dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca',time_step=i)
    dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
    
    pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)

Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))

plot(Stage_0_1)
plot(Stage_1_2)
plot(Stage_2_3)
plot(Stage_3_)

#through time
par(.pardefault)
par(mfrow=c(1,4))

sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')

adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')

sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')

juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')

######################################

####### Permutation 3 ########

# features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit),
#                  'population'=as.populations(stack(koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[1],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[2],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[3],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[4]
#                  )
#                  ),
#                  'carrying_capacity'=as.carrying_capacity(koala.hab.k.s))
# habitat <- as.habitat(features)
# print(habitat)
# 
# hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)
# 
# pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()
# 
# pops_at_time_step[[1]] <- populations(habitat)
# habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat,1)
# 
#   pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
#   for(i in 1:n_time_steps){
#     
#     disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
#     
#     new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),time_step=i)
#     attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
#     habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
#     
#     dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca',time_step=i)
#     dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
#     
#     pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob,habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
#     
#     setTxtProgressBar(pb, i)
#   }
#   close(pb)
# 
# Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
# Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
# Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
# Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))
# 
# plot(Stage_0_1)
# plot(Stage_1_2)
# plot(Stage_2_3)
# plot(Stage_3_)
# 
# #through time
# par(.pardefault)
# par(mfrow=c(1,4))
# 
# sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
# plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')
# 
# adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
# plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')
# 
# sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
# plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')
# 
# juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
# plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')
# 
# ######################################
# 
# ####### Permutation 4 ########
# 
# features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit),
#                  'population'=as.populations(stack(koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[1],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[2],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[3],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[4]
#                  )
#                  ),
#                  'carrying_capacity'=as.carrying_capacity(koala.hab.k.n))
# habitat <- as.habitat(features)
# print(habitat)
# 
# hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)
# 
# pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()
# 
# pops_at_time_step[[1]] <- populations(habitat)
# habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat,1)
# 
#   pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
#   for(i in 1:n_time_steps){
#     
#     disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
#     
#     new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),time_step=i)
#     attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
#     habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
#     
#     dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca')
#     dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
#     
#     pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
#     
#     setTxtProgressBar(pb, i)
#   }
#   close(pb)
# 
# Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
# Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
# Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
# Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))
# 
# plot(Stage_0_1)
# plot(Stage_1_2)
# plot(Stage_2_3)
# plot(Stage_3_)
# 
# #through time
# par(.pardefault)
# par(mfrow=c(1,4))
# 
# sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
# plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')
# 
# adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
# plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')
# 
# sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
# plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')
# 
# juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
# plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')
# 
# ######################################
# 
# ####### Permutation 5 ########
# 
# features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit.s),
#                  'population'=as.populations(stack(koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[1],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[2],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[3],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[4]
#                  )
#                  ),
#                  'carrying_capacity'=as.carrying_capacity(koala.hab.k.n))
# habitat <- as.habitat(features)
# print(habitat)
# 
# hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)
# 
# pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()
# 
# pops_at_time_step[[1]] <- populations(habitat)
# habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat,1)
# 
#   pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
#   for(i in 1:n_time_steps){
#     
#     disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
#     
#     new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat,time_step=i),time_step=i)
#     attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
#     habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat,time_step=i) <- new_hs
#     
#     dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca',time_step=i)
#     dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
#     
#     pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
#     
#     setTxtProgressBar(pb, i)
#   }
#   close(pb)
# 
# Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
# Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
# Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
# Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))
# 
# plot(Stage_0_1)
# plot(Stage_1_2)
# plot(Stage_2_3)
# plot(Stage_3_)
# 
# #through time
# par(.pardefault)
# par(mfrow=c(1,4))
# 
# sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
# plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')
# 
# adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
# plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')
# 
# sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
# plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')
# 
# juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
# plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')
# 
# ######################################
# 
# ####### Permutation 6 ########
# 
# features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit.s),
#                  'population'=as.populations(stack(koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[1],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[2],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[3],
#                                                    koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[4]
#                  )
#                  ),
#                  'carrying_capacity'=as.carrying_capacity(koala.hab.k.s))
# habitat <- as.habitat(features)
# print(habitat)
# 
# hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)
# 
# pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()
# 
# pops_at_time_step[[1]] <- populations(habitat)
# habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat,1)
# 
#   pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
#   for(i in 1:n_time_steps){
#     
#     disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
#     
#     new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat,time_step=i),time_step=i)
#     attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
#     habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat,time_step=i) <- new_hs
#     
#     dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca',time_step=i)
#     dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
#     
#     pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
#     
#     setTxtProgressBar(pb, i)
#   }
#   close(pb)
# 
# Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
# Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
# Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
# Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))
# 
# plot(Stage_0_1)
# plot(Stage_1_2)
# plot(Stage_2_3)
# plot(Stage_3_)
# 
# #through time
# par(.pardefault)
# par(mfrow=c(1,4))
# 
# sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
# plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')
# 
# adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
# plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')
# 
# sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
# plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')
# 
# juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
# plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')
# 
# ######################################

####### Permutation 7 ########

features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit),
                 'population'=as.populations(koala.hab.pop),
                 'carrying_capacity'=as.carrying_capacity(koala.hab.k))

habitat <- as.habitat(features)
print(habitat)

hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)

pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()

pops_at_time_step[[1]] <- populations(habitat)
habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat)

  pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
  for(i in 1:n_time_steps){
    
    disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
    
    new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),time_step=i)
    attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
    habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
    
    dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca',time_step=i)
    dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
    
    pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)

Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))

plot(Stage_0_1)
plot(Stage_1_2)
plot(Stage_2_3)
plot(Stage_3_)

#through time
par(.pardefault)
par(mfrow=c(1,4))

sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')

adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')

sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')

juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')

######################################

####### Permutation 8 ########

features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit),
                 'population'=as.populations(koala.hab.pop.n),
                 'carrying_capacity'=as.carrying_capacity(koala.hab.k))
habitat <- as.habitat(features)
print(habitat)

hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)

pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()

pops_at_time_step[[1]] <- populations(habitat)
habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat)

  pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
  for(i in 1:n_time_steps){
    
    disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
    
    new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),time_step=i)
    attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
    habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
    
    dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca',time_step=i)
    dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
    
    pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)

Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))

plot(Stage_0_1)
plot(Stage_1_2)
plot(Stage_2_3)
plot(Stage_3_)

#through time
par(.pardefault)
par(mfrow=c(1,4))

sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')

adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')

sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')

juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')

######################################

####### Permutation 9 ########

features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit),
                 'population'=as.populations(stack(koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[1],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[2],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[3],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[4]
                 )
                 ),
                 'carrying_capacity'=as.carrying_capacity(koala.hab.k))
habitat <- as.habitat(features)
print(habitat)

koala.disp.param2 <- as.dispersal(c(koala.disp.param,barrier_type=1,barriers_map=koala.disp.bar,use_barriers=TRUE))

hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)

pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()

pops_at_time_step[[1]] <- populations(habitat)
habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat)

  pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
  for(i in 1:n_time_steps){
    
    disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
    
    new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),time_step=i)
    attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
    habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
    
    dispersed_populations <- dispersal(koala.disp.param2,habitat,method='ca',time_step=i)
    dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
    
    pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)

Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))

plot(Stage_0_1)
plot(Stage_1_2)
plot(Stage_2_3)
plot(Stage_3_)

#through time
par(.pardefault)
par(mfrow=c(1,4))

sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')

adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')

sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')

juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')

######################################

####### Permutation 10 ########

features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit),
                 'population'=as.populations(stack(koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[1],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[2],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[3],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[4]
                 )
                 ),
                 'carrying_capacity'=as.carrying_capacity(koala.hab.k))
habitat <- as.habitat(features)
print(habitat)

koala.disp.param2 <- as.dispersal(c(koala.disp.param,barrier_type=1,barriers_map=koala.disp.bar.s,use_barriers=TRUE))

hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param)

pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()

pops_at_time_step[[1]] <- populations(habitat)
habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat)

  pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
  for(i in 1:n_time_steps){
    
    disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
    
    new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),time_step=i)
    attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
    habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
    
    dispersed_populations <- dispersal(koala.disp.param2,habitat,method='ca',time_step=i)
    dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
    
    pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)

Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))

plot(Stage_0_1)
plot(Stage_1_2)
plot(Stage_2_3)
plot(Stage_3_)

#through time
par(.pardefault)
par(mfrow=c(1,4))

sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')

adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')

sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')

juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')

######################################

####### Permutation 11 ########

features <- list('habitat_suitability_map'=as.habitat_suitability(koala.hab.suit),
                 'population'=as.populations(stack(koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[1],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[2],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[3],
                                                   koala.hab.pop*summary(koala.demo.glob)$stable.stage.distribution[4]
                 )
                 ),
                 'carrying_capacity'=as.carrying_capacity(koala.hab.k))
habitat <- as.habitat(features)
print(habitat)

hab_dyn <- as.habitat_dynamics(koala.dist.fire.func,koala.dist.fire.func.param2)

pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()

pops_at_time_step[[1]] <- populations(habitat)
habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat)

  pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
  for(i in 1:n_time_steps){
    
    disturbance_at_time_step[[i]] <- koala.dist.fire
    
    new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),time_step=i)
    attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
    habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
    
    dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca',time_step=i)
    dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
    
    pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=1, seed=NULL, time_step=i)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)

Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))

plot(Stage_0_1)
plot(Stage_1_2)
plot(Stage_2_3)
plot(Stage_3_)

#through time
par(.pardefault)
par(mfrow=c(1,4))

sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')

adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')

sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')

juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')

######################################


#Combining all into a single function

pop_sim <- function(hab.suit,hab.pop,hab.k,demo.obj,demo.sd,dist.func,dist.param,n_time_steps){
  features <- list('habitat_suitability_map'=as.habitat_suitability(hab.suit),
                   'population'=as.populations(stack(hab.pop*summary(demo.obj)$stable.stage.distribution[1],
                                                     hab.pop*summary(demo.obj)$stable.stage.distribution[2],
                                                     hab.pop*summary(demo.obj)$stable.stage.distribution[3],
                                                     hab.pop*summary(demo.obj)$stable.stage.distribution[4]
                   )
                   ),
                   'carrying_capacity'=as.carrying_capacity(hab.k))
  habitat <- as.habitat(features)
  print(habitat)
  
  hab_dyn <- as.habitat_dynamics(dist.func,dist.param)
  
  pops_at_time_step <- dispersal_at_time_step <- disturbance_at_time_step <- habitat_suit_at_time_step <- list()
  
  pops_at_time_step[[1]] <- populations(habitat)
  habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat)
  
  #system.time(
  pb <- txtProgressBar(min = 0, max = n_time_steps, style = 3)
  for(i in 1:n_time_steps){
    
    if(any(i == seq(4,n_time_steps,3))){ # reset back to orginal or modified hs assuming time period of fire effect (3 years)
      r <- habitat_suit_at_time_step[[1]]
      r[] <- 1
      disturbance_at_time_step[[i]] <- r
      habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- habitat_suit_at_time_step[[1]]
    }else{ # change habitat suitability based on disturbance
      disturbance_at_time_step[[i]] <- hab_dyn[[2]][[i]]
      new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),time_step=i)
      attr(new_hs,"habitat") <- "habitat_suitability"
      habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
    }    
    
    dispersed_populations <- dispersal(koala.disp.param,habitat,method='ca',time_step=i) # disperse populations in landscape
    dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
    
    pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob, habitat, stage_matrix_sd=demo.sd, seed=NULL, time_step=i) # update populations based on transition matrix (incorporates stochasticity and density dependence)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  #) ##### approx 30 seconds to do 20 timesteps 
  
  Stage_0_1 <- stack(lapply(pops_at_time_step, "[[", 1))
  Stage_1_2 <- stack(lapply(pops_at_time_step, "[[", 2))
  Stage_2_3 <- stack(lapply(pops_at_time_step, "[[", 3))
  Stage_3_ <- stack(lapply(pops_at_time_step, "[[", 4))

  #par(.pardefault)
  par(mfrow=c(1,4))
  
  sup_adult_n <- sapply(lapply(pops_at_time_step, "[[", 4)[], function(x)sum(x[]))
  plot(sup_adult_n,type='l',ylab='Total Super-Adult Abundance',xlab="Time (years)",lwd=2,col='gray')
  
  adult_n <- sapply(lapply(pops_at_time_step, "[[", 3)[], function(x)sum(x[]))
  plot(adult_n,type='l',ylab='Total Adult Abundance',xlab="Time (years)",lwd=2,col='tomato')
  
  sub_adult_n <- sapply(lapply(pops_at_time_step, "[[", 2)[], function(x)sum(x[]))
  plot(sub_adult_n,type='l',ylab='Total Sub-Adult Abundance',xlab="Time (years)",lwd=2,col='dodgerblue')
  
  juve_n <- sapply(lapply(pops_at_time_step, "[[", 1)[], function(x)sum(x[]))
  plot(juve_n,type='l',ylab='Total Juvenile Abundance',xlab="Time (years)",lwd=2,col='darkgreen')
  
}

pop_sim(koala.hab.suit,koala.hab.pop*0.5,koala.hab.k*.05,koala.demo.glob,0.5,koala.dist.fire.func,koala.dist.fire.func.param,20)




##### pop function #####

# surv_juvenile <- tr(juvenile ~ juvenile, p(0.5))
# surv_subadult <- tr(subadult ~ subadult, p(0.7))
# surv_adult <- tr(adult ~ adult, p(0.8))
# surv_superadult <- tr(superadult ~ superadult, p(0.8))
# 
# at_foot <- tr(subadult ~ juvenile, p(0.5))
# mature <- tr(adult ~ subadult, p(0.8))
# territorial <- tr(superadult ~ adult, p(0.8))
# 
# prob_birth <- tr(juvenile ~ superadult, p(0.5))
# fecundity <- tr(juvenile ~ superadult, r(1))
# 
# 
# survival <- pop::dynamic(surv_juvenile,
#                     surv_subadult,
#                     surv_adult,
#                     surv_superadult)
# 
# growth <- pop::dynamic(at_foot,
#                   mature,
#                   territorial)
# 
# recruitment <- pop::dynamic(prob_birth * fecundity)
# 
# par(mfrow = c(1, 4))
# plot(survival); title(main = 'survival')
# plot(growth); title(main = 'growth')
# plot(recruitment); title(main = 'recruitment')
# 
# all <- pop::dynamic(survival,
#                growth,
#                recruitment)
# 
# par(mfrow = c(1, 2))
# plot(all)
# 
# (A <- as.matrix(all))
# 
# popbio::lambda(A)
# 
# (ss <- popbio::stable.stage(A))
# 
# population <- round(ss * 100)
# 
# plot(popdemo::project(A, population, time = 50))
# 
# sim <- simulation(dynamic = all,
#                   population = population,
#                   timesteps = 50,
#                   replicates = 30)
# 
# plot(sim)

#koala.pop[koala.pop<quantile(koala.hab,probs = c(0.999))] <- 0 # set abundance to zero below a threshold
#koala.pop[koala.pop!=0] <- 100 # set abundance (100 in this case) for remaining values

#random_pop <- sampleRandom(koala.pop, size=50, na.rm=TRUE, sp=TRUE)
#random_pop@data <- as.data.frame(t(rmultinom(50, size = 1, prob = summary(trans)$stable.stage.distribution)))

#koala.k <-  calc(koala.hab, function(x){ifelse( x >= 0.154, x/cellStats(koala.hab,max), 0)}) # create raster of carrying capacity values (k)


# system.time(dispersed_populations <- dispersal(dispersal_params, habitat, method='ca'))
# 
# ##### note, all dispersed populations the same?!?
# plot(dispersed_populations[[1]], box = FALSE, axes = FALSE)
# plot(dispersed_populations[[2]], box = FALSE, axes = FALSE)
# plot(dispersed_populations[[3]], box = FALSE, axes = FALSE)
# plot(dispersed_populations[[4]], box = FALSE, axes = FALSE)
# 
# populations(habitat) <- dispersed_populations # only showing how to update populations - happens inside experiment after each time step

####  THIS IS A TEST FOR A MODULE WITH A CUSTOM FUNCTION ####
# params = list(habitat,
#               fire_start_location = sample(ncell(habitat_suitability(habitat)),10),
#               prob = 0.24,
#               continue_to_burn_prob = 0.01)
# 
# ## check module produces expected output
# fires <- as.module(fire_module,params)       
# 
# fires <- as.module(fire.list)

## altenative you can run the fire module as follows:
# RUN_FIRE_screaming <- fire_module(habitat,sample(ncell(habitat_suitability(habitat)),10),
#                                   prob = 0.24,
#                                   continue_to_burn_prob = 0.01)
##############################

# set up a density dependence function. This isn't working as before :( But it should give you some ideas.
# ddfun <- function (pop, K) {
#   K*pop/(pop+(K-pop)*exp(-0.5)) ##### changed to Beverton-Holt
# }

#now we are going to try and set up an experiment using the above functions.
# experiment <- function(transition, habitat, dispersal, module, dd){} <- ##### WHAT?
# pops_at_time_step <- dispersal_at_time_step <- fire_at_time_step <- habitat_suit_at_time_step <- list()
# 
# ##### this shows how lists are populated...initial values??? no dispersal???
# pops_at_time_step[[1]] <- populations(habitat) 
# habitat_suit_at_time_step[[1]] <- habitat_suitability(habitat) 
# fire_at_time_step[[1]] <- fire_module(habitat,sample(ncell(habitat_suitability(habitat)),10),
#                                       prob = 0.24,
#                                       continue_to_burn_prob = 0.01)

# ddfun <- function (pop, K) {
#   K*pop/(pop+(K-pop)*exp(-0.5)) ##### changed to Beverton-Holt
# }


# #update the landscape based on last fire ##### shouldn't this indexing be different?
# #fire_prob <- fire_at_time_step[[i]]/max(fire_at_time_step[[i]][])   
# 
# #now you could have some equation on how fire effects underlying habitat suitability. I'm just going to make     # something up :) = habitat_suitability*(1-fire_prob) 
# #new_hs <- habitat_suitability(habitat)*(1-fire_prob)
# disturbance_at_time_step[[i]] <- koala.dist.fire.s[[i]]
# new_hs <- run_habitat_dynamics(hab_dyn,habitat_suitability(habitat),i)
# attr(new_hs,"habitat") <- "habitat_suitability" # i need to update attribute inheritance.
# habitat_suit_at_time_step[[i+1]] <- habitat_suitability(habitat) <- new_hs
# 
# 
# # let's disperse people
# dispersed_populations <- dispersal(koala.disp.param, habitat, method='ca')
# dispersal_at_time_step[[i]] <- populations(habitat) <- dispersed_populations
# 
# # get a vector of populations per life history
# #pop_vec <- lapply(populations(habitat),function(x)c(x[]))
# #pop_mat <- do.call(cbind,pop_vec)
# 
# # update populations (this could be done much more nicely with pop)
# #pops_n <- t(koala.demo.glob$global_stage_matrix%*%t(pop_mat))
# pops_at_time_step[[i+1]] <- populations(habitat) <- estimate_demography(koala.demo.glob,habitat)
# 
# ## update density dependence for adult populations.  <- ##### why adults???
# #pops_n[,4] <- ddfun(pops_n[,4], cellStats(koala.hab.k, max))
# 
# #now update the populations - ideally this could be turned into a function (update_pops())
# #r <- habitat_suitability(habitat)
# #pops_updated <- lapply(split(pops_n, rep(1:ncol(pops_n), each = nrow(pops_n))),function(x){r[]<-x;return(r)})
# #pops <- lapply(pops_updated, `attr<-`, "habitat", "populations")
# #pops_at_time_step[[i+1]] <- populations(habitat) <- pops
# 
# #now let's start another fire! Mwhahahhaha.
# # fire_at_time_step[[i+1]] <- fire_module(habitat,sample(ncell(habitat_suitability(habitat)),10),
# #                                         prob = 0.24,
# #                                         continue_to_burn_prob = 0.01)