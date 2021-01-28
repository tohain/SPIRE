#include "surface_projection.hpp"
#include "distance_transform.hpp"
#include "img_out.hpp"
#include "homotopic_thinning.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "percolation_analysis.hpp"

#include <png.h>

void print_menue(){

  std::vector<std::string> options = {"SET_NTUCS", "SET_A", "SET_TYPE", "SET_LEVEL_SET", "EXIT"};

  for(unsigned int ii=0; ii<options.size(); ii++){
    std::cout << ii << ":   " << options.at(ii) << std::endl;;
  }
    

}




int main( int argc, char* argv[] ){

  double progress;
  std::string status;
  
  /*
   * create a projection first
   */
  surface_projection sp( progress, status );
  sp.set_n_points_x( 150 );
  sp.set_n_points_z( 75 );
  sp.update_geometry();
  sp.update_containers();

  sp.compute_projection();


  return EXIT_SUCCESS;
}
