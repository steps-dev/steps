% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_results-class.R
\name{plot.simulation_results}
\alias{plot.simulation_results}
\title{Plot the results of a simulation}
\usage{
\method{plot}{simulation_results}(x, object = "population",
  type = "graph", stages = NULL, animate = FALSE,
  timesteps = c(1:3), panels = c(3, 3), emp = FALSE, ...)
}
\arguments{
\item{x}{a simulation_results object}

\item{object}{the \code{simulation_results} object to plot - can be 'population'
(default), 'suitability' or 'carrying_capacity'}

\item{type}{the plot type - 'graph' (default) or 'raster'}

\item{stages}{life-stages to plot - must be specified for 'raster' plot
types; by default all life-stages will be plotted}

\item{animate}{if plotting type 'raster' would you like to animate the
rasters?}

\item{timesteps}{timesteps to display when plotting rasters}

\item{panels}{the number of columns and rows to use when plotting raster
timeseries - default is 3 x 3 (e.g. c(3,3) )}

\item{emp}{should the expected minimum population of the simulation be
plotted?}

\item{...}{further arguments passed to/from other methods}
}
\description{
Methods to visually inspect the results of a simulation. Both linear graphs
and spatial-explicit grids can be generated to illustrate population changes through
time and space.
}
\examples{

\dontrun{
ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)

pd <- population_dynamics(change = growth(egk_mat),
                          dispersal = kernel_dispersal(max_distance = c(0, 0, 2000),
                          dispersal_kernel = exponential_dispersal_kernel(distance_decay = 1000)),
                          density_dependence = ceiling_density())

sim <- simulation(landscape = ls,
                  population_dynamics = pd,
                  habitat_dynamics = NULL,
                  timesteps = 20)

# Plot the population trajectories by life-stage
plot(sim)

# Plot the spatial distributions of total populations for first nine timesteps
plot(sim, type = "raster", stages = 0, timesteps = 1:9) 
}
}