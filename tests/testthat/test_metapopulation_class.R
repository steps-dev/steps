context('metapopulation-class')

habitat <- as.habitat(list(coordinates = data.frame(x=runif( 10, min=-10, max=10),
                                                     y=runif( 10, min=-10, max=10)),
                                area = as.data.frame(exp(-seq(.1,10,length.out = 10))*10),
                                population = as.population(t(rmultinom(10, 
                                size = 100, prob = c(0.8,0.2,0.01)))),
                                features = data.frame(temperature = 10)))
params <- list(alpha=1,beta=1,disp_fun="H")
adult.disperal <- dispersal(params) 

test_that('metapopulation classes work', {

  mp  <-  metapopulation(nrep=1, time=10, habitat=habitat, dispersal=adult.disperal,
                       x1 = 0.42, e1 = 0.0061, y1 = 1.2)
  mp1 <- metapopulation(nrep=100, time=100, habitat, dispersal=adult.disperal,
                        x1 = 0.2, e1 = 0.0061, y1 = 1)

})