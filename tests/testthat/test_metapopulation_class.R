context('metapopulation-class')

n <- 50
meta_data <- data.frame(x1=runif( n, min=-10, max=10), x2=runif( n, min=-10, max=10),
area=exp(-seq(.1,10,length.out = n))*10,presence=rbinom(n,1,.8))
area <- meta_data$area
dist <- as.matrix(with(meta_data, dist(cbind(x1, x2))))
presence <- meta_data$presence
locations <- meta_data[,c('x1','x2')]

test_that('metapopulation classes work', {

  mp <- metapopulation(nrep=1, time=10, dist=dist, area=area, presence=presence,locations=locations,
                       x1 = 0.42, e1 = 0.0061, y1 = 1.2)
  mp1 <- metapopulation(nrep=10, time=10, dist=dist, area=area, presence=presence,locations=locations,
                       x1 = 0.42, e1 = 0.0061, y1 = 1.2)

})