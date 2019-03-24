#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix na_matrix(int nr, int nc){
  NumericMatrix m(nr, nc);
  std::fill(m.begin(), m.end(), NumericVector::get_na());
  return m;
}


// [[Rcpp::export]]
IntegerVector shuffle_vec(int min, int max){
  IntegerVector vec = seq(min, max);
  std::random_shuffle(vec.begin(), vec.end());
  return vec;
}


// [[Rcpp::export]]
bool barrier_to_dispersal(int source_x,
                          int source_y,
                          int sink_x,
                          int sink_y,
                          NumericMatrix barriers_map){
  
  /*
   ** Initialise parameters used in the function
   */
  
  int  dist_x, dist_y, i, pixel_x, pixel_y, distance_max;
  bool barrier_found = false;
  int barrier_counter = 0;
  
  /*
   ** Calculate the distance in both dimensions between the source and sink
   ** pixels and take the largest of the two.
   */
  
  dist_x = source_x - sink_x;
  dist_y = source_y - sink_y;
  if (abs(dist_x) >= abs(dist_y)){
    distance_max = abs(dist_x);
  } else {
    distance_max = abs(dist_y);
  }

  /*
   ** Check the possible paths (5 variants) from source to sink and see if there is a path
   ** without barriers. If any barrier cell is found, the barrier counter is incremented.
   */
  
  for (i = 1; i <= distance_max; i++){
    pixel_x = round(sink_x + (1.0 * i / distance_max * dist_x));
    pixel_y = round(sink_y + (1.0 * i / distance_max * dist_y));
    if (R::rbinom(1, barriers_map(pixel_x, pixel_y)) == 1){
      barrier_counter++;
      break;
    }
  }
  
  for (i = 1; i <= distance_max; i++){
    pixel_x = round(sink_x - 0.49 + (((i - 1.0) / distance_max * dist_x) +
      ((1.0 / distance_max * dist_x) / 2.0)));
    pixel_y = round(sink_y - 0.49 + (((i - 1.0) / distance_max * dist_y) +
      ((1.0 / distance_max * dist_y) / 2.0)));
    if (R::rbinom(1, barriers_map(pixel_x, pixel_y)) == 1){
      barrier_counter++;
      break;
    }
  }
  
  for (i = 1; i <= distance_max; i++){
    pixel_x = round(sink_x + 0.49 + (((i - 1.0) / distance_max * dist_x) +
      ((1.0 / distance_max * dist_x) / 2.0)));
    pixel_y = round(sink_y - 0.49 + (((i - 1.0) / distance_max * dist_y) +
      ((1.0 / distance_max * dist_y) / 2.0)));
    if (R::rbinom(1, barriers_map(pixel_x, pixel_y)) == 1){
      barrier_counter++;
      break;
    }
  }
  
  for (i = 1; i <= distance_max; i++){
    pixel_x = round(sink_x - 0.49 + (((i - 1.0) / distance_max * dist_x) +
      ((1.0 / distance_max * dist_x) / 2.0)));
    pixel_y = round(sink_y + 0.49 + (((i - 1.0) / distance_max * dist_y) +
      ((1.0 / distance_max * dist_y) / 2.0)));
    if (R::rbinom(1, barriers_map(pixel_x, pixel_y)) == 1){
      barrier_counter++;
      break;
    }
  }
  
  for (i = 1; i <= distance_max; i++){
    pixel_x = round(sink_x + 0.49 + (((i - 1.0) / distance_max * dist_x) +
      ((1.0 / distance_max * dist_x) / 2.0)));
    pixel_y = round(sink_y + 0.49 + (((i - 1.0) / distance_max * dist_y) +
      ((1.0 / distance_max * dist_y) / 2.0)));
    if (R::rbinom(1, barriers_map(pixel_x, pixel_y)) == 1){
      barrier_counter++;
      break;
    }
  }
  
  /*
   ** If any barriers cells have been found (dependent on permeability)
   ** return barrier_found, otherwise, return no barrier_found
   */
  
  if (barrier_counter > 0){
    barrier_found = true;
  }
  
  return(barrier_found);
  
}


// [[Rcpp::export]]
IntegerVector can_source_cell_disperse(int source_x,
                                       int source_y,
                                       NumericMatrix carrying_capacity_available, 
                                       NumericMatrix habitat_suitability_map,
                                       NumericMatrix barriers_map,
                                       int dispersal_distance,
                                       NumericVector dispersal_kernel){
  
  /*
   ** Initialise parameters used in the function
   */
  
  int ncols = carrying_capacity_available.ncol();
  int nrows = carrying_capacity_available.nrow();
  int sink_x, sink_y, real_distance;
  IntegerVector sink_x_vec = shuffle_vec(source_x - dispersal_distance, source_x + dispersal_distance);
  IntegerVector sink_y_vec = shuffle_vec(source_y - dispersal_distance, source_y + dispersal_distance);
  double prob_colonisation, rnd;
  IntegerVector sink_found(2, -1);
  bool barrier;
  
  for (IntegerVector::iterator sink_x = sink_x_vec.begin(); sink_x != sink_x_vec.end(); sink_x++){
    for (IntegerVector::iterator sink_y = sink_y_vec.begin(); sink_y != sink_y_vec.end(); sink_y++){
      
      /*
       ** 1. Test of basic conditions to see if a pixel could be a potential sink cell:
       **  - The pixel must be within the limits of the matrix's extent.
       **  - The pixel must have available carrying capacity to allow recruitment.
       **  - The pixel must not be NA.
       */
      
      if ((*sink_x >= 0) && (*sink_x < nrows) && (*sink_y >= 0) && (*sink_y < ncols)){
        if ((carrying_capacity_available(*sink_x, *sink_y) > 0) && !R_IsNA(carrying_capacity_available(*sink_x, *sink_y))){
          if(!R_IsNA(habitat_suitability_map(*sink_x, *sink_y))){
            
            /*
             ** 2. Compute the distance between potential sink and source pixel
             ** and check if it is less than or equal to maximum dispersal distance.
             ** The distance is computed in pixel units.
             */
            
            real_distance = round(sqrt((*sink_x - source_x) * (*sink_x - source_x) + (*sink_y - source_y) * (*sink_y - source_y)));
            
            if ((real_distance > 0) && (real_distance <= dispersal_distance)){
              
              /*
               ** 3. Compute the probability of colonisation of the sink pixel.
               ** This probability depends on two factors:
               **  - Distance between source and sink cells.
               **  - Habitat suitability.
               */
              
              prob_colonisation = dispersal_kernel[real_distance - 1] * (habitat_suitability_map(*sink_x, *sink_y));
              
              /*
               ** Generate random number
               */
              
              rnd = as<double>(Rcpp::runif(1));
              
              /*
               ** Add stochasticity to colonisation event
               */
              
              if (rnd < prob_colonisation){
                
                /*
                 ** 4. When we reach this stage, the last thing we need to check for
                 ** is whether there is a "barrier" obstacle between the source
                 ** and sink pixel. We check this last as it requires significant
                 ** computing time.
                 */
                
                barrier = barrier_to_dispersal(source_x, source_y, *sink_x, *sink_y, barriers_map);
                
                if (!barrier){
                  sink_found[0] = *sink_x;
                  sink_found[1] = *sink_y;
                }
                
              }
            }
          }
        }
      }
    }
  }
  
  return(sink_found);
  
}


// [[Rcpp::export]]
int proportion_of_population_to_disperse(int source_x,
                                         int source_y,
                                         int sink_x,
                                         int sink_y, 
                                         NumericMatrix starting_population_state,
                                         NumericMatrix current_carrying_capacity, 
                                         double dispersal_proportion){
  
  /*
   ** Initialise parameters used in the function
   */
  
  int source_pop;
  int source_pop_dispersed = 0;
  
  /*
   ** Get intitial population
   */
  
  source_pop = starting_population_state(source_x, source_y);
  
  /*
   ** Only calculate proportion of population to disperse if
   ** there is population in the source and if there is available
   ** space in the target sink cell
   */
  
  if(!R_IsNA(source_pop) && (source_pop >= 1)){
    
    source_pop_dispersed = R::rbinom(source_pop, dispersal_proportion);
    
    if (current_carrying_capacity(sink_x, sink_y) <= source_pop_dispersed){
      source_pop_dispersed = 0;
    }
    
  }
  
  return(source_pop_dispersed);
  
}

 
// [[Rcpp::export]]
NumericMatrix rcpp_dispersal(NumericMatrix starting_population_state,
                    NumericMatrix potential_carrying_capacity,
                    NumericMatrix habitat_suitability_map,
                    NumericMatrix barriers_map,
                    int dispersal_steps,
                    int dispersal_distance,
                    NumericVector dispersal_kernel,
                    double dispersal_proportion){
  
  /*
   ** Initialise parameters used in the function
   */
  
  int ncols = starting_population_state.ncol();
  int nrows = starting_population_state.nrow();
  NumericMatrix carrying_capacity_available = na_matrix(nrows, ncols);
  NumericMatrix future_population_state;
  int dispersal_step, i, j;
  IntegerVector source_x_vec = shuffle_vec(0, (ncols - 1));
  IntegerVector source_y_vec = shuffle_vec(0, (nrows - 1));

  /*
   ** fill future population matrix with zeros
   */
  
  for(i = 0; i < nrows; i++){
    for(j = 0; j < ncols; j++){
      future_population_state(i, j) = 0;
    }
  }
  
  /*
   ** create a carrying capacity available matrix by subtracting
   ** current population state from carrying capacity of landscape
   */
  
  for(i = 0; i < nrows; i++){
    for(j = 0; j < ncols; j++){
      if(!R_IsNA(starting_population_state(i, j))){
        carrying_capacity_available(i, j) = potential_carrying_capacity(i, j) - starting_population_state(i, j);
        if(carrying_capacity_available(i, j) < 0) carrying_capacity_available(i, j) = 0;
        if(barriers_map(i, j) == 1) carrying_capacity_available(i, j) = 0;
      }
    }
  }
  
  /* *********************** */
  /* Dispersal starts here.  */
  /* *********************** */
  
  for(dispersal_step = 0; dispersal_step < dispersal_steps; dispersal_step++){
    
    for (IntegerVector::iterator source_x = source_x_vec.begin(); source_x != source_x_vec.end(); source_x++){
      for (IntegerVector::iterator source_y = source_y_vec.begin(); source_y != source_y_vec.end(); source_y++){
        
        if(!R_IsNA(starting_population_state(*source_x, *source_y))){
          
          /*
           ** Sink cell search: Can the source pixel disperse? There are four
           ** conditions to be met for a source pixel to disperse:
           **  1. Sink pixel is currently suitable.
           **  2. Sink pixel has available carrying capacity.
           **  3. Sink pixel is within dispersal distance of source cell.
           **  4. There is no obstacle (barrier) between the pixel to be colonised
           **  (sink pixel) and the pixel that is already colonised (source pixel).
           */
          
          IntegerVector cell_in_dispersal_distance = can_source_cell_disperse(*source_x,
                                                                              *source_y,
                                                                              carrying_capacity_available,
                                                                              habitat_suitability_map,
                                                                              barriers_map,
                                                                              dispersal_distance,
                                                                              dispersal_kernel);
          
          if(cell_in_dispersal_distance[0] >= 0 && cell_in_dispersal_distance[1] >= 0){
            
            int sink_x = cell_in_dispersal_distance[0];
            int sink_y = cell_in_dispersal_distance[1];
            
            int source_pop_dispersed = proportion_of_population_to_disperse(*source_x,
                                                                            *source_y,
                                                                            sink_x,
                                                                            sink_y,
                                                                            starting_population_state,
                                                                            carrying_capacity_available,
                                                                            dispersal_proportion);
            
            starting_population_state(*source_x, *source_y) = starting_population_state(*source_x, *source_y) - source_pop_dispersed;
            
            future_population_state(sink_x, sink_y) = source_pop_dispersed;
            
          }
        }
      }
    }
    
    for(i = 0; i < nrows; i++){
      for(j = 0; j < ncols; j++){
        
        future_population_state(i, j) = future_population_state(i, j) + starting_population_state(i, j);
        
      }
    }
    
  }
  
  return(future_population_state);
}