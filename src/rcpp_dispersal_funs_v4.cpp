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
IntegerVector can_source_cell_disperse(int source_x,
                                       int source_y,
                                       NumericMatrix iterative_population_state,
                                       NumericMatrix future_population_state,
                                       NumericMatrix carrying_capacity_available, 
                                       NumericMatrix permeability_map,
                                       int max_cells){
  
  /*
   ** Initialise parameters used in the function
   */
  
  int ncols = carrying_capacity_available.ncol();
  int nrows = carrying_capacity_available.nrow();
  int cell;
  IntegerVector sink_found(2, -1);
  int dest_x;
  int dest_y;
  int ndir = 4; // Rook's case
  IntegerVector dest_x_vec(ndir);
  IntegerVector dest_y_vec(ndir);
  IntegerVector possible(ndir);
  NumericVector prob_vec(ndir);
  int direction, sink_carrying_cap;
  IntegerVector selected_direction(1);
  
  /*
   ** Assign initial destination coordinates from source cell
   */
  
  dest_x = source_x;
  dest_y = source_y;
  
  for (cell = 0; cell < max_cells; cell++){ // Increment cell movements up to a maximum
    
    /*
     ** Fill vectors with surrounding cell locations (rook's case)
     */
    
    dest_x_vec[0] = dest_x;
    dest_x_vec[1] = dest_x + 1;
    dest_x_vec[2] = dest_x;
    dest_x_vec[3] = dest_x - 1;
    
    dest_y_vec[0] = dest_y + 1;
    dest_y_vec[1] = dest_y;
    dest_y_vec[2] = dest_y - 1;
    dest_y_vec[3] = dest_y;
    
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
      
      if(dest_x_vec[direction] >= 0 &&
         dest_x_vec[direction] < nrows &&
         dest_y_vec[direction] >= 0 &&
         dest_y_vec[direction] < ncols &&
         !R_IsNA(carrying_capacity_available(dest_y_vec[direction], dest_x_vec[direction]))){
         
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
    
    dest_x = dest_x_vec[selected_direction[0]];
    dest_y = dest_y_vec[selected_direction[0]];
    
    /*
     ** Check for carrying capacity and, if available, return destination cell
     */
    
    sink_carrying_cap = carrying_capacity_available(dest_y, dest_x) - (iterative_population_state(dest_y, dest_x) + future_population_state(dest_y, dest_x));
    
    if(sink_carrying_cap >= 1){
      sink_found[0] = dest_x;
      sink_found[1] = dest_y;
      break;
    }
    
  }
  
  return(sink_found);
  
}


// [[Rcpp::export]]
List rcpp_dispersal(NumericMatrix starting_population_state,
                    NumericMatrix potential_carrying_capacity,
                    NumericMatrix permeability_map,
                    int max_cells,
                    double dispersal_proportion){
  
  /*
   ** Initialise parameters used in the function
   */
  
  int ncols = starting_population_state.ncol();
  int nrows = starting_population_state.nrow();
  NumericMatrix carrying_capacity_available(nrows, ncols);
  NumericMatrix iterative_population_state(nrows, ncols);
  NumericMatrix future_population_state(nrows, ncols);
  NumericMatrix dispersers(nrows, ncols);
  NumericMatrix failed_dispersers(nrows, ncols);
  int i, j, individual;
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
      if(permeability_map(i, j) == 1) carrying_capacity_available(i, j) = 0;
      
    }
  }
  
  /* *********************** */
  /* Dispersal starts here.  */
  /* *********************** */
  
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
           ** Sink cell search: Can the source pixel disperse? The each individual will
           ** search the landscape within the boundaries until it finds a sink pixel
           ** with available carrying capacity.
           */
          
          IntegerVector eligible_sink_cell_xy = can_source_cell_disperse(source_x,
                                                                         source_y,
                                                                         iterative_population_state,
                                                                         future_population_state,
                                                                         carrying_capacity_available,
                                                                         permeability_map,
                                                                         max_cells);
          
          if(eligible_sink_cell_xy[0] >= 0 && eligible_sink_cell_xy[1] >= 0){
            
            int sink_x = eligible_sink_cell_xy[0];
            int sink_y = eligible_sink_cell_xy[1];
            
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
  
  return(List::create(Named("future_population") = future_population_state,
                      Named("dispersed") = dispersers,
                      Named("failed") = failed_dispersers));
}