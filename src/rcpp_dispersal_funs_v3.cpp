#include <Rcpp.h>
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm - from http://gallery.rcpp.org/articles/stl-random-shuffle/
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
IntegerVector shuffle_vec(int min, int max){
  IntegerVector vec = seq(min, max);
  std::random_shuffle(vec.begin(), vec.end(), randWrapper);
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
    pixel_x = round( (double) sink_x + (1.0 * i / distance_max * dist_x) );
    pixel_y = round( (double) sink_y + (1.0 * i / distance_max * dist_y) );
    if (R::rbinom(1, barriers_map(pixel_y, pixel_x)) == 1){
      barrier_counter++;
      break;
    }
  }
  
  for (i = 1; i <= distance_max; i++){
    pixel_x = round( (double) sink_x - 0.49 + (((i - 1.0) / distance_max * dist_x ) +
      ((1.0 / distance_max * dist_x) / 2.0)));
    pixel_y = round( (double) sink_y - 0.49 + (((i - 1.0) / distance_max * dist_y ) +
      ((1.0 / distance_max * dist_y) / 2.0)));
    if (R::rbinom(1, barriers_map(pixel_y, pixel_x)) == 1){
      barrier_counter++;
      break;
    }
  }
  
  for (i = 1; i <= distance_max; i++){
    pixel_x = round( (double) sink_x + 0.49 + (((i - 1.0) / distance_max * dist_x ) +
      ((1.0 / distance_max * dist_x) / 2.0)));
    pixel_y = round( (double) sink_y - 0.49 + (((i - 1.0) / distance_max * dist_y ) +
      ((1.0 / distance_max * dist_y) / 2.0)));
    if (R::rbinom(1, barriers_map(pixel_y, pixel_x)) == 1){
      barrier_counter++;
      break;
    }
  }
  
  for (i = 1; i <= distance_max; i++){
    pixel_x = round( (double) sink_x - 0.49 + (((i - 1.0) / distance_max * dist_x ) +
      ((1.0 / distance_max * dist_x) / 2.0)));
    pixel_y = round( (double) sink_y + 0.49 + (((i - 1.0) / distance_max * dist_y ) +
      ((1.0 / distance_max * dist_y) / 2.0)));
    if (R::rbinom(1, barriers_map(pixel_y, pixel_x)) == 1){
      barrier_counter++;
      break;
    }
  }
  
  for (i = 1; i <= distance_max; i++){
    pixel_x = round( (double) sink_x + 0.49 + (((i - 1.0) / distance_max * dist_x ) +
      ((1.0 / distance_max * dist_x) / 2.0)));
    pixel_y = round( (double) sink_y + 0.49 + (((i - 1.0) / distance_max * dist_y ) +
      ((1.0 / distance_max * dist_y) / 2.0)));
    if (R::rbinom(1, barriers_map(pixel_y, pixel_x)) == 1){
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
                                       NumericMatrix iterative_population_state,
                                       NumericMatrix future_population_state,
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
  int i, j, real_distance;
  IntegerVector sink_x_vec = shuffle_vec(source_x - dispersal_distance, source_x + dispersal_distance);
  IntegerVector sink_y_vec = shuffle_vec(source_y - dispersal_distance, source_y + dispersal_distance);
  int nx = sink_x_vec.size();
  int ny = sink_y_vec.size();
  double prob_colonisation, rnd;
  IntegerVector sink_found(2, -1);
  bool barrier;
  
  for (i = 0; i < nx; i++){
    for (j = 0; j < ny; j++){
      
      /*
       ** Assign sink coordinates from shuffled vector of integers
       */
      
      int sink_x = sink_x_vec[i];
      int sink_y = sink_y_vec[j];
      
      /*
       ** 1. Test of basic conditions to see if a pixel could be a potential sink cell:
       **  - The pixel must be within the limits of the matrix's extent.
       **  - The pixel must have available carrying capacity to allow recruitment.
       **  - The pixel must not be NA.
       */
      
      if ((sink_x >= 0) && (sink_x < nrows) && (sink_y >= 0) && (sink_y < ncols)){
        
        int sink_carrying_cap = carrying_capacity_available(sink_y, sink_x) - (iterative_population_state(sink_y, sink_x) + future_population_state(sink_y, sink_x));
        
        if ((sink_carrying_cap > 0) && !R_IsNA(sink_carrying_cap)){
          if(!R_IsNA(habitat_suitability_map(sink_y, sink_x))){
            
            /*
             ** 2. Compute the distance between potential sink and source pixel
             ** and check if it is less than or equal to maximum dispersal distance.
             ** The distance is computed in pixel units.
             */
            
            real_distance = round( (double) sqrt( (double) (sink_x - source_x) * (sink_x - source_x) + (sink_y - source_y) * (sink_y - source_y) ) );
            
            if ((real_distance > 0) && (real_distance <= dispersal_distance)){
              
              /*
               ** 3. Compute the probability of colonisation of the sink pixel.
               ** This probability depends on two factors:
               **  - Distance between source and sink cells.
               **  - Habitat suitability.
               */
              
              prob_colonisation = dispersal_kernel[real_distance - 1] * (habitat_suitability_map(sink_y, sink_x));
              
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
                
                barrier = barrier_to_dispersal(source_x, source_y, sink_x, sink_y, barriers_map);
                
                if (!barrier){
                  sink_found[0] = sink_x;
                  sink_found[1] = sink_y;
                  break;
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
List rcpp_dispersal(NumericMatrix starting_population_state,
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
  NumericMatrix carrying_capacity_available(nrows, ncols);
  NumericMatrix iterative_population_state(nrows, ncols);
  NumericMatrix dispersers(nrows, ncols);
  NumericMatrix failed_dispersers(nrows, ncols);
  NumericMatrix future_population_state(nrows, ncols);
  int dispersal_step, i, j, individual;
  IntegerVector source_x_vec = shuffle_vec(0, (ncols - 1));
  IntegerVector source_y_vec = shuffle_vec(0, (nrows - 1));
  int nx = source_x_vec.size();
  int ny = source_y_vec.size();
  
  /*
   ** copy values from source to initialised matrices and set zeros
   ** in carrying capacity available matrix where barriers are present
   */
  
  for(i = 0; i < nrows; i++){
    for(j = 0; j < ncols; j++){
      
      iterative_population_state(i, j) = starting_population_state(i, j);
            
      carrying_capacity_available(i, j) = potential_carrying_capacity(i, j);
      if(barriers_map(i, j) == 1) carrying_capacity_available(i, j) = 0;

    }
  }
  
  /* *********************** */
  /* Dispersal starts here.  */
  /* *********************** */
  
  for(dispersal_step = 0; dispersal_step < dispersal_steps; dispersal_step++){
    
    for (i = 0; i < nx; i++){
      for (j = 0; j < ny; j++){
        
        /*
         ** Assign source coordinates from shuffled vector of integers
         */
        
        int source_x = source_x_vec[i];
        int source_y = source_y_vec[j];        
        
        /*
         ** Verify if there is population in the cell
         */
        
        if(!R_IsNA(iterative_population_state(source_y, source_x)) && iterative_population_state(source_y, source_x) > 0){
          
          /*
           ** Realised number of individuals that can disperse - based on dispersal proportion
           */
          
          int source_population = iterative_population_state(source_y, source_x);
          int source_pop_dispersed = R::rbinom(source_population, dispersal_proportion);
          
          dispersers(source_y, source_x) = source_pop_dispersed;
          
          /*
           ** Loop through individuals that can disperse
           */
          
          for (individual = 0; individual < source_pop_dispersed; individual++){
            
            /*
             ** Sink cell search: Can the source pixel disperse? There are four
             ** conditions to be met for a source pixel to disperse:
             **  1. Sink pixel is currently suitable.
             **  2. Sink pixel has available carrying capacity.
             **  3. Sink pixel is within dispersal distance of source cell.
             **  4. There is no obstacle (barrier) between the pixel to be colonised
             **  (sink pixel) and the pixel that is already colonised (source pixel).
             */
            
            IntegerVector cell_in_dispersal_distance = can_source_cell_disperse(source_x,
                                                                                source_y,
                                                                                iterative_population_state,
                                                                                future_population_state,
                                                                                carrying_capacity_available,
                                                                                habitat_suitability_map,
                                                                                barriers_map,
                                                                                dispersal_distance,
                                                                                dispersal_kernel);
            
            if(cell_in_dispersal_distance[0] >= 0 && cell_in_dispersal_distance[1] >= 0){
              
              int sink_x = cell_in_dispersal_distance[0];
              int sink_y = cell_in_dispersal_distance[1];
 
              iterative_population_state(source_y, source_x) = iterative_population_state(source_y, source_x) - 1;
              
              future_population_state(sink_y, sink_x) = future_population_state(sink_y, sink_x) + 1;
              
            } else {
              
              failed_dispersers(source_y, source_x) = failed_dispersers(source_y, source_x) + 1;
            
            }
          }
        }
      }
    }
    
    for(i = 0; i < nrows; i++){
      for(j = 0; j < ncols; j++){
        
        future_population_state(i, j) = future_population_state(i, j) + iterative_population_state(i, j);
        
      }
    }
    
  }
  
  return(List::create(Named("future_population") = future_population_state,
                      Named("dispersed") = dispersers,
                      Named("failed") = failed_dispersers));
}