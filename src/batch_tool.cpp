/* SPIRE - Structure Projection Image Recognition Environment
 * Copyright (C) 2021 Tobias Hain
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see https://www.gnu.org/licenses.
 */


#include "surface_projection.hpp"
#include "distance_transform.hpp"
#include "img_out.hpp"
#include "homotopic_thinning.hpp"

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>

#include "percolation_analysis.hpp"

#include "auxiliary.hpp"
#include <png.h>

#include "batch_lib.hpp"


typedef struct {
  std::vector<double> membranes;
  std::vector<double> filled_channels;
  std::vector<double> resolution;
  std::string fn_prefix;

  std::string mode = "ordered";
  unsigned long N;
  bool quadratic = false;
} cmd_options;

void print_help(){

  std::cout << "batch tool; create a large batch of projections" << std::endl
	    << std::endl
	    << "fixed parameters for the entire sweep:" << std::endl
	    << "   --membranes       : comma separated pairs of pos,width" << std::endl
            << "                       ALWAYS includes a membrane of 0,0.3." << std::endl
	    << "                       can be disabled using filled_channels" << std::endl
	    << "   --filled_channels : comma separated number of 1 or 0." << std::endl
	    << "                       1 is a filled channel, 0 can disable membranes" << std::endl
	    << "   --resolution      : two comma separated integers. Resolution in" << std::endl
	    << "                       x,y dimension" << std::endl
	    << "   --fn_prefix       : prefix to attach to all output files" << std::endl
	    << "   --mode            : (ordered) iterate through all paramter combinations" << std::endl
    	    << "                     : (random)  pick random parameters" << std::endl
	    << "   --N               : number of projections to be computed (random mode only)" << std::endl
	    << "   --quadratic       : only allows quadratic projections, height is set to width (random mode only)" << std::endl
    
	    << std::endl

	    << "sweep parameters, each parameter can be assign a comma separated list," << std::endl
	    << "or a range given as start:stop:step" << std::endl
	    << "in case of random mode, provide upper and lower bound for the interval"
	    << "in which samples are drawn uniformly" << std::endl
	    << std::endl
	    << "last given parameter is increased most frequently, first least frequently" << std::endl
	    << "recommending using miller indeces as last ones, otherwise configurations" << std::endl
	    << "might be skipped" << std::endl;
  

  for( auto it : surface_projection::parameter_names ){
    std::cout << "   --" << it << std::endl;
  }
  
  std::cout << std::endl << std::endl;
  
}





/** 
 * check the most critical parameters for consistency and apply them
 */
bool check_and_set_fixed_pars( cmd_options &ops, surface_projection &sp ){

  if( ops.membranes.size() % 2 == 1 ){
    throw std::string ("odd number of parameters in membrane array found");
    return false;
  }
  
  // check if 0 membrane is given
  bool found_zero = false;
  for( unsigned int ii=0; ii<ops.membranes.size(); ii+=2){
    if( ops.membranes[ii] == 0 ){
      found_zero = true;
      break;
    }
  }
  if( !found_zero ){
    throw std::string( "couldn't find main membrane" );
    return false;
  }

  // check if channel filled number matches number of channels
  unsigned int nr_channels = ops.membranes.size() + 1;
  if( ops.filled_channels.size() != nr_channels ){
    throw std::string("nr of channels mismatches number of filled_channels parameters" );
    return false;
  }


  // check resolution
  if( ops.resolution.size() != 2 ){
    throw std::string("error resolution; provide pair of x,z");
    return false;
  }

  
  // set membranes and filled channels, not sure if this is needed
  // each time something changes  
  sp.set_membranes( ops.membranes );

  std::vector<int> fills (ops.filled_channels.size(), 0 );
  for( unsigned int ii=0; ii < ops.filled_channels.size(); ii++){
    fills[ii] =static_cast<int>( ops.filled_channels[ii] );
  }
  sp.set_channel_fill( fills );

  sp.set_n_points_x( static_cast<int>( ops.resolution[0] ) );
  sp.set_n_points_z( static_cast<int>( ops.resolution[1] ) );		     
  
  return true;

}


std::vector<double> split_and_map( std::string in ){

  std::vector<std::string> split = my_utility::str_split( in, ',' );
  std::vector<double> out ( split.size(), 0 );
  for( size_t ii=0; ii<split.size(); ii++){
    out[ii] = std::stod( split[ii] );
  }
  return out;
}


void parse_cmd_options( int argc, char* argv[], batch_creation &bc, cmd_options &ops ){
  
  char *found_it = NULL;
  
  for( int ii=1; ii<argc; ii++) {    
    for( unsigned int jj=0; jj<surface_projection::parameter_names.size(); jj++){

      // check all the iterable parameters
      found_it = strstr( argv[ii], surface_projection::parameter_names[jj].c_str() );
      
      if( found_it != NULL ){
	bc.add_parameter( surface_projection::parameter_names[jj], std::string( argv[ii+1] ) );
      }
    }


    /*
     * now check the "fixed" parameters
     */
    if( strcmp( argv[ii], "--membranes" ) == 0 ){
      ops.membranes = split_and_map( std::string( argv[ii+1] ) );
    }
    if( strcmp( argv[ii], "--filled_channels" ) == 0 ){
      ops.filled_channels = split_and_map ( std::string( argv[ii+1] ) );
    }

    if( strcmp( argv[ii], "--resolution" ) == 0 ){
      ops.resolution = split_and_map ( std::string( argv[ii+1] ) );
    }

    if( strcmp( argv[ii], "--fn_prefix" ) == 0 ){
      ops.fn_prefix = std::string( argv[ii+1] );
    }

    if( strcmp( argv[ii], "--mode" ) == 0 ){
      ops.mode = std::string( argv[ii+1] );
    }

    if( strcmp( argv[ii], "--N" ) == 0 ){
      ops.N = strtoul( argv[ii+1], NULL, 10 );
    }

    if( strcmp( argv[ii], "--quadratic") == 0 ){
      ops.quadratic = true;
    }
    
  }
    
}




class cmd_callback : public sp_callback {

public:
  cmd_callback( surface_projection &sp,
		batch_creation &bc_,
		std::string prefix ) :
    sp_callback( sp ),
    fn_prefix( prefix ),
    bc( bc_ ),
    counter (0)
  {
    summary.open( fn_prefix + "_summary.txt" );
  }

  ~cmd_callback(){
    summary.close();
  }
  
  bool operator()( std::vector< std::vector<double>::iterator > pars ){

    std::string fn = fn_prefix + "_" + std::to_string( counter ) + ".png";

    sp.validate_uc_scaling();
    sp.update_geometry();
    sp.compute_projection();
    sp.save_to_png( fn, invert, "LIN" );

    summary << fn << "    ";

    for( auto it : surface_projection::parameter_names ){
      summary << it << "=" << sp.get_parameter( it ) << " ";
    }

    summary << std::endl;
    
      
    counter++;
    return true;
  }

  long counter;
  
  batch_creation &bc;
  std::string fn_prefix;
  std::ofstream summary;
};


int main( int argc, char* argv[] ){

  if( argc == 1 ){
    print_help();
  }

  global_settings gs ( "global_settins.conf" );
  cmd_options ops;
  surface_projection sp ( gs );
  batch_creation bc ( sp, time(NULL) );
  
  try {
    parse_cmd_options( argc, argv, bc, ops );
    if( !check_and_set_fixed_pars( ops, sp ) ){
      return EXIT_FAILURE;
    }
  } catch ( std::string s ){
    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // get the functor

  cmd_callback cmdcb ( sp, bc, ops.fn_prefix );
  cmdcb.set_img_mode( true, "LIN" );
  
  if( ops.mode == "ordered" ){    
    long elements_created = bc.do_loop( cmdcb  );
    std::cout << "made " << elements_created << " projections" << std::endl;
  } else if ( ops.mode == "random" ){

    std::cout << "creating " << ops.N << " random projections" << std::endl;

    long long counter = 0;
    std::ofstream f_out ( ops.fn_prefix + "_faulty.txt" );
    
    std::vector<std::vector<double>::iterator > dummy_it;    
    for( unsigned long ii=0; ii<ops.N; ii++ ){

      if( !bc.set_random_parameters( ops.quadratic ) ){
	f_out << counter << " " << " caught invalid_parameter_exception" << std::endl;
      }

      counter++;
      cmdcb( dummy_it );
    }

    f_out.close();
    
  }
  
  return EXIT_SUCCESS;
}


