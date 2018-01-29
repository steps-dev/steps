#######################################################################
## Collate koala data for pop model example                        ####
#######################################################################

# Step 1 - crop habitat suitability, survival, fecundity & fire layers

## CROP RASTER CODE - just including for my own reference ##########################################
library(rgdal)

wd<-"C:/Users/nbriscoe/ownCloud/PopDynModels_BrigalowBeltSth" # Location of data
wd_out<-"C:/Users/nbriscoe/ownCloud/Brigalow_Belt_Sth/Koala_ex" #Save output here

## Habitat suitability 
HS<-raster(paste0(wd,"/Rasters/Maxent_HS_model/Pred_mod11_mask_Alb_NM_Nov.tif"))
plot(HS) # Can update later...

# Crop to smaller rectangle that has some variation (large + smaller patches/high and low quality)
Crop.shp<-readOGR(dsn=paste0(wd,"/Shapefiles&Tables"), layer="Koala_crop_ex")

HS_crop<-crop(HS,Crop.shp)
plot(HS_crop)

writeRaster(HS_crop, paste0(wd_out,"/Habitat Suitability/HS_crop.tif"))

## Survival of age classes - crop
F1<-stack(list.files("C:/Users/nbriscoe/ownCloud/PopDynModels_BrigalowBeltSth/Rasters/UpIntake_F1/TestProj", full=TRUE,pattern="\\.tif$"))
FNR2<-stack(list.files("C:/Users/nbriscoe/ownCloud/PopDynModels_BrigalowBeltSth/Rasters/UpIntake_F2NR/TestProj", full=TRUE,pattern="\\.tif$"))
FR2<-stack(list.files("C:/Users/nbriscoe/ownCloud/PopDynModels_BrigalowBeltSth/Rasters/UpIntake_F2R/TestProj", full=TRUE,pattern="\\.tif$"))
FNR3<-stack(list.files("C:/Users/nbriscoe/ownCloud/PopDynModels_BrigalowBeltSth/Rasters/UpIntake_F3NR/TestProj", full=TRUE,pattern="\\.tif$"))
FR3<-stack(list.files("C:/Users/nbriscoe/ownCloud/PopDynModels_BrigalowBeltSth/Rasters/UpIntake_F3R/TestProj", full=TRUE,pattern="\\.tif$"))

F1<-crop(F1, Crop.shp)
FNR2<-crop(FNR2,Crop.shp)
FR2<-crop(FR2,Crop.shp)
FNR3<-crop(FNR3,Crop.shp)
FR3<-crop(FR3,Crop.shp)

writeRaster(F1, paste0(wd_out,"/Dyn Sur & Fec/Sur_F1/Sur_F1_.tif"), bylayer=TRUE)
writeRaster(FR2, paste0(wd_out,"/Dyn Sur & Fec/Rep_F2/Sur_F2R_.tif"), bylayer=TRUE)
writeRaster(FNR2, paste0(wd_out,"/Dyn Sur & Fec/Sur_F2/Sur_F2NR_.tif"), bylayer=TRUE)
writeRaster(FR3, paste0(wd_out,"/Dyn Sur & Fec/Rep_F3 & Sur_F0/Sur_F3R_.tif"), bylayer=TRUE)
writeRaster(FNR3, paste0(wd_out,"/Dyn Sur & Fec/Sur_F3/Sur_F3NR_.tif"), bylayer=TRUE)

## Note that NR = multiplier on survival of stage, R = multiplier on reproduction, survival of F0 = survival of F3R 

# Fire
Fire<-stack(list.files("C:/Users/nbriscoe/ownCloud/PopDynModels_BrigalowBeltSth/Rasters/Fire_Albers", full=TRUE, pattern="\\.tif$"))
plot(Fire)
Fire<-crop(Fire, Crop.shp)
writeRaster(Fire, paste0(wd_out,"/Fire/Fire_.tif"),bylayer=TRUE)

#########################################################################

## Step 2 - Set up pop dynamics using pop  ##################################

library(pop)
library(popdemo)

## Set up the transitions
# For now use 'advance' function to get around limits on transitions 0< p >1
advance <- as.transfun(function (landscape) {1},
                       param = list(),
                       type = "probability")

# probability of surviving - in this case all transition to next life stage, except S3 
survival_S0 <- tr(S0 ~ S0, p(0.94))
survival_S1 <- tr(S1 ~ S1, p(0.884))
survival_S2 <- tr(S2 ~ S2, p(0.793))
survival_S3 <- tr(S3 ~ S3, p(0.793))

# probability of moving to the next state
M0_1 <- tr(S1 ~ S0, advance) # Augustine 1998
M1_2 <- tr(S2 ~ S1, advance) # Dique et al., 2003)
M2_3 <- tr(S3 ~ S2, advance) # # Kavenaugh et al 2007 
M3_3 <-tr(S3 ~ S3, advance) # Kavenaugh et al 2007

survival <- dynamic(survival_S0,survival_S1,survival_S2, survival_S3)
#survival <- dynamic(survival_S3)

growth <- dynamic(M0_1,
                  M1_2,
                  M2_3)
# probability of rep and number of offspring
prob_repS2 <- tr(S0 ~ S2, p(0.7619))
prob_repS3<-tr(S0 ~ S3, p(0.7619)) # 
fecundityS2 <- tr(S0 ~ S2, r(0.5)) # Female offspring only
fecundityS3 <- tr(S0 ~ S3, r(0.5)) # Female offspring only

# we want the product of these two things, which we might estimate separately
recruitment <- dynamic(prob_repS2 * fecundityS2,
                       prob_repS3 * fecundityS3)

# note that the order of transitions in the dynamic determines the order in
# which transitions occur, and may influence the dynamics. Can't make them grow (transition to next stage) before they have a chance to reproduct

all <- dynamic(survival,
               recruitment,
               growth)

(A<- as.matrix(all))

popbio::lambda(A)

(ss <- popbio::stable.stage(A))

# get an initial population
population <- round(ss * 100)

# plot deterministic trajectory from this population
par(mfrow = c(1, 1))
plot(popdemo::project(A, population, time = 50))

sim <- simulation(dynamic = all,
                  population = population,
                  timesteps = 50,
                  replicates = 30)

# plot abundance of the three life stages
par(mfrow = c(1, 4))
plot(sim)

## Density dependence ## ##NEEDS FIXING ##
# density-dependent adult survival function (patch is a required argument)
landscape <- as.landscape(NULL)
ddfun <- function (landscape) {
  all_density <- (population(landscape)$S3 + population(landscape)$S2 + population(landscape)$S1) /area(landscape)$area
  param$p * exp(-all_density/param$area)
}
# Make density dependent on all stages except S0

##What should dd be? Change to match data...in best quality habitat 2 koala per ha (1 female)

# turn it into a transfun object
dd <- as.transfun(ddfun,
                  param = list(p = 0.9,
                               area = 100),
                  type = 'probability')

# use it in a transition and update the dynamic
survival_S1_dd <- tr(S1 ~ S1, dd)
survival_S2_dd <- tr(S2 ~ S2, dd)
survival_S3_dd <- tr(S3 ~ S3, dd)

# updated survival transitions
survival_dd<-dynamic(survival_S0,survival_S1_dd,survival_S2_dd, survival_S3_dd)
all_dd <- dynamic(survival_dd,
                  recruitment,
                  growth)

# run the simulation (a little longer and with more simulations this time)
sim_dd <- simulation(dynamic = all_dd,
                     population = population,
                     timesteps = 100,
                     replicates = 100)

# and plot it
par(mfrow = c(1, 4))
plot(sim_dd)

proj_dd <- projection(dynamic = all_dd,
                      population = population,
                      timesteps = 100)
# and plot it
par(mfrow = c(1, 4))
plot(proj_dd)


## Metapopulation dyanmics ##
?landscape
patches <- as.landscape(list(population = population(landscape(all_dd)),
                             area = data.frame(area = c(0.5, 2, 5)),
                             coordinates = data.frame(x = c(-10, 3, -2),
                                                      y = c(-10, -1, 1)),
                             features = features(landscape(all_dd))))
# plot it
par(mfrow = c(1, 2))
symbols(x = patches[, 1],
        y = patches[, 2],
        circles = sqrt(patches[, 3] / pi),
        xlim = c(-10, 4),
        ylim = c(-10, 3),
        fg = grey(0.4),
        bg = grey(0.9),
        xlab = 'x',
        ylab = 'y',
        inches = FALSE,
        asp = 1)
text(x = patches[, 1],
     y = patches[, 2],
     labels = 1:3)

##
all_metapop <- dynamic(all_dd) # Does this mean dispersal happens last?
landscape(all_metapop) <- patches

sim_metapop <- simulation(dynamic = all_metapop,
                          population = population,
                          timesteps = 100,
                          replicates = 100)

# and plot it
par(mfrow = c(3, 3))
plot(sim_metapop, 'S1', patches = 1:3)
plot(sim_metapop, 'S2', patches = 1:3)
plot(sim_metapop, 'S3', patches = 1:3)

# Add dispersal
#Dispersal for F1 (0.35*exp(-Dij^1/3.36)) - for now just use below
Dij<-seq(0,11,by=0.2)

Disp_F1<-0.35*exp(-Dij*(1/3.36))
Disp_F2<-0.25*exp(-Dij*(1/3.36))
par(mfrow=c(1,2))
plot(Disp_F1~Dij)
plot(Disp_F2~Dij)

S1_dispersal <- tr(S1 ~ S1, p(0.35) * d(1/3.36))
S2_dispersal <- tr(S2 ~ S2, p(0.25) * d(1/3.36))

# add this to the dynamic and assign the spatial structure
all_metapop_disp <- dynamic(all_dd, S1_dispersal,S2_dispersal) # Does this mean dispersal happens last?
landscape(all_metapop_disp) <- patches

sim_metapop_disp <- simulation(dynamic = all_metapop_disp,
                          population = population,
                          timesteps = 100,
                          replicates = 100)
# and plot it
par(mfrow = c(3, 3))
plot(sim_metapop_disp, 'S1', patches = 1:3)
plot(sim_metapop_disp, 'S2', patches = 1:3)
plot(sim_metapop_disp, 'S3', patches = 1:3)

## What happens to individuals that can't disperse?

distance(patches) ## Patch 1 too far away to get dispersal - move closer and see what happens

patches <- as.landscape(list(population = population(landscape(all_dd)),
                             area = data.frame(area = c(0.5, 2, 5)),
                             coordinates = data.frame(x = c(1, 3, -2),
                                                      y = c(-2, -1, 1)),
                             features = features(landscape(all_dd))))
landscape(all_metapop_disp) <- patches # update patches
sim_metapop_disp <- simulation(dynamic = all_metapop_disp,
                               population = population,
                               timesteps = 100,
                               replicates = 100)

# and plot it
par(mfrow = c(3, 3))
plot(sim_metapop_disp, 'S1', patches = 1:3)
plot(sim_metapop_disp, 'S2', patches = 1:3)
plot(sim_metapop_disp, 'S3', patches = 1:3)

# See how things are stored...$simulations[replicates][timestep,stage*patch]
range(unlist(sim_metapop_disp$simulations[[1]][,1])) # min & max population size for the first replicate for first life stage in first patch
