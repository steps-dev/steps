#include <Rcpp.h>
#include <random>
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm - from http://gallery.rcpp.org/articles/stl-random-shuffle/
//inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
IntegerVector shuffle_vec(int min, int max){
  IntegerVector vec = seq(min, max);
  std::shuffle(vec.begin(), vec.end(), std::mt19937(std::random_device()()));
  return vec;
}

// [[Rcpp::export]]
IntegerVector can_source_cell_disperse(int source_y,
                                       int source_x,
                                       NumericMatrix iterative_population_state,
                                       NumericMatrix future_population_state,
                                       NumericMatrix carrying_capacity_available, 
                                       NumericMatrix permeability_map,
                                       int max_cells,
                                       int min_cells){
  
  /*
   ** Initialise parameters used in the function
   */
  
  int ny = carrying_capacity_available.nrow();
  int nx = carrying_capacity_available.ncol();
  int cell;
  IntegerVector sink_found(2, -1);
  int dest_y;
  int dest_x;
  int ndir = 4; // Rook's case
  IntegerVector dest_y_vec(ndir);
  IntegerVector dest_x_vec(ndir);
  IntegerVector possible(ndir);
  NumericVector prob_vec(ndir);
  int direction, sink_carrying_cap;
  IntegerVector selected_direction(1);
  
  /*
   ** Assign initial destination coordinates from source cell
   */
  
  dest_y = source_y;
  dest_x = source_x;
  
  for (cell = 0; cell < max_cells; cell++){ // Increment cell movements up to a maximum
    
    /*
     ** Fill vectors with surrounding cell locations (rook's case)
     */
    
    dest_y_vec[0] = dest_y + 1;
    dest_y_vec[1] = dest_y;
    dest_y_vec[2] = dest_y - 1;
    dest_y_vec[3] = dest_y;
    
    dest_x_vec[0] = dest_x;
    dest_x_vec[1] = dest_x + 1;
    dest_x_vec[2] = dest_x;
    dest_x_vec[3] = dest_x - 1;
    
    /*
     ** Reset parameters used in the loop
     */
    
    possible[0] = 0;
    possible[1] = 0;
    possible[2] = 0;
    possible[3] = 0;
    
    for (direction = 0; direction < ndir; direction++){
      
      /*
       ** Which directions are possible?
       */
      
      if(dest_y_vec[direction] >= 0 &&
         dest_y_vec[direction] < ny &&
         dest_x_vec[direction] >= 0 &&
         dest_x_vec[direction] < nx &&
         !R_IsNA(carrying_capacity_available(dest_y_vec[direction], dest_x_vec[direction])) &&
         !R_IsNaN(carrying_capacity_available(dest_y_vec[direction], dest_x_vec[direction])) &&
         !R_IsNA(permeability_map(dest_y_vec[direction], dest_x_vec[direction])) &&
         !R_IsNaN(permeability_map(dest_y_vec[direction], dest_x_vec[direction])) &&
         permeability_map(dest_y_vec[direction], dest_x_vec[direction]) > 0){
        
        possible[direction] = 1;
        
      }
      
    }
    
    if (is_true(all(possible == 0))){
      
      /*
       ** If no directions are possible, dispersal fails
       */
      
      break;
      
    } else {
      
      /*
       ** Reset parameters used in the loop
       */
      
      prob_vec[0] = 0;
      prob_vec[1] = 0;
      prob_vec[2] = 0;
      prob_vec[3] = 0;
      
      for (direction = 0; direction < ndir; direction++){
        
        /*
         ** Of the directions that are possible, which are most likely?
         */
        
        if(possible[direction] == 1){
          
          prob_vec[direction] = permeability_map(dest_y_vec[direction], dest_x_vec[direction]);
          
        }
        
      }
      
    }
    
    /*
     ** Choose direction based on habitat suitability and landscape barriers
     ** (permeability above)
     */
    
    selected_direction = sample(ndir, 1, false, prob_vec, false);
    
    /*
     ** Update destination coordinates with most likely sink cell
     */
    
    dest_y = dest_y_vec[selected_direction[0]];
    dest_x = dest_x_vec[selected_direction[0]];
    
    /*
     ** Check for carrying capacity after minimum number of cells is reached and,
     ** if available, return destination cell
     */
    
    if (cell > min_cells) {
      sink_carrying_cap = 0;
      
      if(!R_IsNA(carrying_capacity_available(dest_y, dest_x)) &&
         !R_IsNaN(carrying_capacity_available(dest_y, dest_x)) &&
         !R_IsNA(iterative_population_state(dest_y, dest_x)) &&
         !R_IsNaN(iterative_population_state(dest_y, dest_x)) &&
         !R_IsNA(future_population_state(dest_y, dest_x)) &&
         !R_IsNaN(future_population_state(dest_y, dest_x))){
         
         sink_carrying_cap = carrying_capacity_available(dest_y, dest_x) - (iterative_population_state(dest_y, dest_x) + future_population_state(dest_y, dest_x));
        
      }
      
      if(sink_carrying_cap >= 1){
        sink_found[0] = dest_y;
        sink_found[1] = dest_x;
        break;
      }
      
    }
  }
  return(sink_found);
}


// [[Rcpp::export]]
List rcpp_dispersal(NumericMatrix starting_population_state,
                    NumericMatrix potential_carrying_capacity,
                    NumericMatrix permeability_map,
                    int max_cells,
                    int min_cells,
                    double dispersal_proportion){
  
  /*
   ** Initialise parameters used in the function
   */
  
  int ny = starting_population_state.nrow();
  int nx = starting_population_state.ncol();
  NumericMatrix carrying_capacity_available(ny, nx);
  NumericMatrix iterative_population_state(ny, nx);
  NumericMatrix future_population_state(ny, nx);
  NumericMatrix dispersers(ny, nx);
  NumericMatrix failed_dispersers(ny, nx);
  int y, x, individual;
  IntegerVector source_y_vec = shuffle_vec(0, (ny - 1));
  IntegerVector source_x_vec = shuffle_vec(0, (nx - 1));
  
  /*
   ** copy values from source to initialised matrices and set zeros
   ** in carrying capacity available matrix where barriers are present
   */
  
  for(y = 0; y < ny; y++){
    for(x = 0; x < nx; x++){
      
      iterative_population_state(y, x) = starting_population_state(y, x);
      
      carrying_capacity_available(y, x) = potential_carrying_capacity(y, x);
      if(permeability_map(y, x) == 0) carrying_capacity_available(y, x) = 0;
      
    }
  }
  
  /* *********************** */
  /* Dispersal starts here.  */
  /* *********************** */
  
  for (y = 0; y < ny; y++){
    for (x = 0; x < nx; x++){
      
      /*
       ** Assign source coordinates from shuffled vector of integers
       */
      
      int source_y = source_y_vec[y];        
      int source_x = source_x_vec[x];
      
      /*
       ** Verify if there is population in the cell
       */
      
      if(!R_IsNA(iterative_population_state(source_y, source_x)) &&
         !R_IsNaN(iterative_population_state(source_y, source_x)) &&
         iterative_population_state(source_y, source_x) > 0){
        
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
           ** Sink cell search: Can the source pixel disperse? Each individual will
           ** search the landscape within the boundaries until it finds a sink pixel
           ** with available carrying capacity.
           */
          
          IntegerVector eligible_sink_cell_yx = can_source_cell_disperse(source_y,
                                                                         source_x,
                                                                         iterative_population_state,
                                                                         future_population_state,
                                                                         carrying_capacity_available,
                                                                         permeability_map,
                                                                         max_cells,
                                                                         min_cells);
          
          if(eligible_sink_cell_yx[0] >= 0 && eligible_sink_cell_yx[1] >= 0){
            
            int sink_y = eligible_sink_cell_yx[0];
            int sink_x = eligible_sink_cell_yx[1];
            
            iterative_population_state(source_y, source_x) = iterative_population_state(source_y, source_x) - 1;
            
            future_population_state(sink_y, sink_x) = future_population_state(sink_y, sink_x) + 1;
            
          } else {
            
            failed_dispersers(source_y, source_x) = failed_dispersers(source_y, source_x) + 1;
            
          }
        }
      }
    }
  }
  
  for(y = 0; y < ny; y++){
    for(x = 0; x < nx; x++){
      
      future_population_state(y, x) = future_population_state(y, x) + iterative_population_state(y, x);
      
    }
  }
  
  return(List::create(Named("future_population") = future_population_state,
                      Named("dispersed") = dispersers,
                      Named("failed") = failed_dispersers));
}
