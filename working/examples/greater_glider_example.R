library(steps)
library(raster)
library(viridis)
library(doMC)
library(foreach)
library(future)
plan(multiprocess)
plan(sequential)

############ MODEL INPUTS ###############

#Note, this is an age-structured matrix and is not generic
gg_trans_mat <- matrix(c(0.00,0.00,0.50,
                         0.50,0.00,0.00,
                         0.00,0.85,0.85),
                         nrow = 3, ncol = 3, byrow = TRUE)
colnames(gg_trans_mat) <- rownames(gg_trans_mat) <- c('Newborn','Juvenile','Adult')

gg_stable_states <- abs( eigen(gg_trans_mat)$vectors[,1] / base::sum(eigen(gg_trans_mat)$vectors[,1]) ) 

#gg_hab_suit <- raster("working/misc/Petauroides_volans_SDM.GH_PM.tif") / 1000
gg_hab_suit <- aggregate(raster("working/misc/Petauroides_volans_SDM.GH_PM.tif") / 1000, fact=10, fun=mean)
names(gg_hab_suit) <- "Habitat"
#plot(gg_hab_suit, box = FALSE, axes = FALSE, col = viridis(100))

gg_hab_k <- floor(gg_hab_suit*30)
gg_hab_k[is.na(gg_hab_k)] <- 0
names(gg_hab_k) <- "Carrying Capacity"
#plot(gg_hab_k, box = FALSE, axes = FALSE, col = viridis(3))

gg_popN <- stack(replicate(ncol(gg_trans_mat), floor(gg_hab_suit*30)))
gg_popN <- gg_popN*gg_stable_states

registerDoMC(cores=4)
gg_pop <- stack(
  foreach(i = 1:nlayers(gg_popN)) %dopar% {
    m <- ceiling(cellStats(gg_popN[[i]], max, na.rm=T))
    calc(gg_popN[[i]], fun = function(x) rbinom(prob = (x/m), size = m, n = 1))
  })
names(gg_pop) <- colnames(gg_trans_mat)
#plot(gg_pop, box = FALSE, axes = FALSE, col = viridis(100))

TotpopN <- sum(cellStats(gg_pop, 'sum', na.rm=T)) # Get total population size to check sensible
init_pop_size <- sum(gg_pop)
#plot(init_pop_size, box = FALSE, axes = FALSE, col = viridis(25))

gg_disp_bar <- gg_hab_suit*0
gg_disp_bar[cellFromRow(gg_disp_bar,c(nrow(gg_disp_bar)/2,(nrow(gg_disp_bar)/2)+1))] <- 1
#plot(gg_disp_bar, box = FALSE, axes = FALSE, col = viridis(100))

# gg_pop_source <- gg_pop[[3]]
# gg_pop_source[] <- 0
# gg_pop_source[sample(which(getValues(gg_pop[[3]]) >= 3), 25)] <- 1
# plot(gg_pop_source, box = FALSE, axes = FALSE)
# 
# gg_pop_sink <- gg_pop[[3]]
# gg_pop_sink[] <- 0
# gg_pop_sink[sample(which(getValues(gg_pop[[1]]) == 1 |
#                             getValues(gg_pop[[2]]) == 1 |
#                             getValues(gg_pop[[3]]) == 1),
#                       cellStats(gg_pop_source, sum))] <- 1
# plot(gg_pop_sink, box = FALSE, axes = FALSE)
# 
# gg_surv <- list(stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F03R+')),
#                    stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F01+')),
#                    stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F02NR+')),
#                    stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F03NR+')))
# 
# gg_fec <- list(NULL,
#                   NULL,
#                   stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F02R+')),
#                   stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F03R+')))



gg_landscape <- landscape(population = gg_pop,
                          suitability = gg_hab_suit,
                          carrying_capacity = gg_hab_k)

gg_pop_dynamics <- population_dynamics(change = growth(transition_matrix = gg_trans_mat,
                                                       global_stochasticity = 0.1),
                                       dispersal = fast_dispersal(stages = 3,
                                                                  dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.5),
                                                                  dispersal_proportion = 0.5),
                                       modification = NULL,
                                       density_dependence = population_cap(stages = 3))

gg_results <- simulation(landscape = gg_landscape,
                         population_dynamics = gg_pop_dynamics,
                         habitat_dynamics = NULL,
                         timesteps = 50,
                         replicates = 3)

plot(gg_results)

plot(gg_results[1], type = "raster", stages = 0, animate = TRUE, timesteps = c(5,10,15,20,25,30,35,40,45,50))


