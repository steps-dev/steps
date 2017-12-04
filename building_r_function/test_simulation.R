library(raster)
#devtools::install_github('skiptoniam/dhmpr')
library(dhmpr)
set.seed(42)
xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
Dd <- as.matrix(dist(xy))
w <- exp(-1/nrow(xy) * Dd)
Ww <- chol(w)
xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
coordinates(xy) <- ~x+y
r <- rasterize(xy, raster(points2grid(xy)), 'z')
r2 <- raster(r)
res(r2) <- 0.01

r2 <- resample(r, r2)
proj4string(r2) <- '+proj=longlat'
rthr <- r2>quantile(r2[],.8)

## now the fun begins
n_patches <- ncol(rthr)*ncol(rthr)
n_stages <- length(states(trans))

#need the pops bit first mat %*% pop
range01 <- function(x){(x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))}

probs <- range01(r2[])
pops <- matrix(NA,n_patches,3)
for(i in seq_len(n_patches))pops[i,] <- rmultinom(1,100*probs[i],c(.80,.40,.10))

pop_s1 <- pops[,1]
pop_s2 <- pops[,2]
pop_s3 <- pops[,3]

r1 <- r2 <- r3 <- rthr
r1[]<-pop_s1
r2[]<-pop_s2
r3[]<-pop_s3

y <- seq_len(ncol(r1))
x <- seq_len(ncol(r1))

# dispersal function acting on distance matrix
# cut-off dispersal at the minimum dimension of the grid 
f <- function (d, cutoff = 2) {
  ifelse (d > cutoff, 0, exp(-d))
}

# f <- function (d) exp(-d)

# initial population on grid (one stage)
# pop <- matrix(rpois(length(x) * length(y),10),length(y), length(x))

# setup for the fft approach (run this once, before the simulation)
fs <- setupFFT(x = x, y = y, f = f, factor = 2)

#stage based dispersal 
# fs[[1]]<-larval_dispersal_fun



ddfun <- function (pop) {
  adult_density <- pop 
  .9 * exp(-adult_density/80)* pop
}
#this is an example of threshold 
ddfun2 <- function (pop) {
  adult_density <- pop 
  ifelse(pop>400,400,pop)
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
  
  ## density dependence bit working. 
  pops_n <- pop_mat %*% (trans$stage_matrix)
  pops_n[,3] <- ddfun2(pops_n[,3])
  
  
  pops_time[[j+1]]<-lapply(split(pops_n, rep(1:ncol(pops_n), each = nrow(pops_n))),function(x)matrix(x,nrow = nrow(pops_new[[1]])))
  cat(j,'\n')
}
  
N <- lapply(pops_time,sum)

larvae <- lapply(pops_time, "[[", 1)
juve <- lapply(pops_time, "[[", 2)
adult <- lapply(pops_time, "[[", 3)

N <- lapply(adult,mean)
# 
# # turn it into a transfun object
# dd <- dhmpr::as.customfun(ddfun,
#                        param = list(p = 0.9,
#                                area = 100),
#                         type = 'probability')

# r[]<-ddfun(N,250)

larv <- lapply(larvae, function(x) raster(matrix(x,nrow = nrow(r2))))
st.l <- stack(larv)

juv <- lapply(juve, function(x) raster(matrix(x,nrow = nrow(r2))))
st.j <- stack(juv)

adl <- lapply(adult, function(x) raster(matrix(x,nrow = nrow(r2))))
st.a <- stack(adl)

nl <- 10
prefix <- paste(sample(letters, 6, replace=TRUE), collapse='')
grDevices::png(sprintf('%s/%s_%%0%sd.png', tempdir(), prefix, nchar(nl)),
               type='cairo', width=width, height=height)
plot_t <- function(i) {
  p <- rasterVis::levelplot(st.a[[i]], margin=FALSE)
  print(p)
}

mapply(plot_t, seq_len(nl)) 
grDevices::dev.off()

width = 800; height=800
oopt <- animation::ani.options(
  ani.width=width, ani.height=height, interval=0.05, ani.dev='png', 
  ani.type='png', nmax=nl, autobrowse=FALSE)
on.exit(animation::ani.options(oopt), add=TRUE)
ff <- list.files(tempdir(), sprintf('^%s.*\\.png$', prefix), full.names=TRUE)
animation::im.convert(ff,output = 'c:/Users/swoolley/Desktop/test.gif', extra.opts='-dispose Background', 
                      clean=TRUE)
