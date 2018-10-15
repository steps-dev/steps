fire_fun <- function (ncol = 30, nrow = 30) {
  r <- raster(ncol = ncol, nrow = nrow)
  n <- ncol * nrow
  vals <- rbinom(n, 1, 0.01) * runif(n)
  r[] <- vals / max(vals)
  r
}

one_year <- function (fire, hs) {
  hs * (1 - fire)
}

rescale <- function (x) {
  x <- x - min(x)
  x / max(x)
}

# regeneration over
lag_year <- function (timestep = 1,
                      fire,
                      hs,
                      lag = 3,
                      regeneration_function = function (time) {-time}) {

  # lags  
  years_since_fire <- c(0, seq_len(lag))
  
  # fire weights
  relative_regeneration <- regeneration_function(years_since_fire)
  regeneration <- rescale(relative_regeneration)
  
  lags <- timestep - years_since_fire
  valid <- lags > 0
  lags <- lags[valid]
  
  weights <- regeneration[valid]
  
  # apply the regeneration weights to previous fires
  fires_weighted <- fire[[lags]] * weights
  
  # get habitat reduction in each previous years, downweighted by regeneration
  # time
  annual_impact <- 1 - fires_weighted
  
  # get the cumulative impact
  impact <- prod(annual_impact)
  
  # apply the habitat reduction
  hs * impact
  
}

library (raster)
set.seed(123)
n_years <- 10

fire_severity <- stack(replicate(n_years, fire_fun()))

hs_base <- fire_fun()
hs_base[] <- runif(length(hs_base), 0.9, 1)

plot(fire_severity)
plot(hs_base, zlim = c(0, 1))

hs_annual <- lapply(seq_len(n_years),
       lag_year,
       fire_severity,
       hs_base)


hs_annual_stack <- stack(hs_annual)

idx <- 5

plot(fire_severity[[idx]])
plot(hs_annual_stack[[idx]], zlim = c(0, 1))

plot(hs_annual_stack, zlim = c(0, 1))


plot(one_year(fire_severity, hs))
