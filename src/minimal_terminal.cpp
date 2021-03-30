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

void print_menue(){

  std::vector<std::string> options = {"SET_NTUCS", "SET_A", "SET_TYPE", "SET_LEVEL_SET", "EXIT"};

  for(unsigned int ii=0; ii<options.size(); ii++){
    std::cout << ii << ":   " << options.at(ii) << std::endl;;
  }
    

}


int main( int argc, char* argv[] ){

  
  return EXIT_SUCCESS;
}
