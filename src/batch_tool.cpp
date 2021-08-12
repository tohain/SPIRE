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
#include "img_manip.hpp"

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
  std::string sum_format;

  std::string mode = "ordered";
  unsigned long N;
  bool quadratic = false;

  double gaussian_noise = -1;
  double gaussian_blur = -1;
  std::vector<double> grains = std::vector<double> (0,0);
  std::vector<std::string> post_processing = std::vector<std::string> (0,"");


  std::vector<std::string> to_skip;
  std::vector<int> strides;
  
} cmd_options;

void print_help(){

  std::cout << "batch tool; create a large batch of projections" << std::endl
	    << std::endl
	    << "fixed parameters for the entire sweep:" << std::endl
	    << "   --membranes       : comma separated pairs of pos,width" << std::endl
            << "                       ALWAYS includes a membrane of 0,0.03" << std::endl
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
	    << "   --summary_format  : the format of the summary file (human,csv)" << std::endl
	    << "   --gaussian_noise  : adds gaussian noise to the image. provide magnitude as parameter" << std::endl
	    << "   --guassian_blur   : adds a gaussian blur to the image. Provide kernel size as parameter" << std::endl
	    << "   --add_grains      : adds grains to the image; comma separated parameters:" << std::endl
	    << "                     : grain_size_mu,grain_size_sigma,grain_number_mu,grain_number_sigam,intensity,blur,mix" << std::endl
	    << "   --skip            : comma separated list of parameters and strides (parameter0,n0,parameter1,n1,...) " << std::endl
	    << "                     : to keep constant for a number of computed projections" << std::endl


	    << "   --post_processing : post_processing mode, only applies specified filters to space" << std::endl
	    << "                       separated list of files. Must come last or will cause errors!" << std::endl
	    << "                       OVERWRITES FILES!!!!" << std::endl
    
	    << std::endl

	    << "sweep parameters, each parameter can be assign a comma separated list," << std::endl
	    << "or a range given as start:stop:step" << std::endl
	    << "Random mode: provide comma separated tupel. First character indicates"
	    << " random distribution, further floats are distribution parameters:" << std::endl
	    << "0: uniform, with parameters lower and upper bound: 0,lo,hi" << std::endl
	    << "1: normal, with parameters mu and sigma: 1,mu,sigma" << std::endl    
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

  size_t pos_found = std::string::npos;
  
  for( int ii=1; ii<argc; ii++) {
    
    for( unsigned int jj=0; jj<surface_projection::parameter_names.size(); jj++){

      // check all the iterable parameters
      pos_found = std::string( argv[ii] ).find( surface_projection::parameter_names[jj] );
            
      if( pos_found != std::string::npos && (pos_found == 2 || pos_found == 1 ) ) {
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

    if( strcmp( argv[ii], "--summary_format" ) == 0 ){
      ops.sum_format = std::string( argv[ii+1] );
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


    // post processing
    if( strcmp( argv[ii], "--gaussian_noise" ) == 0 ){
      ops.gaussian_noise = std::stod( argv[ii+1] );
    }
    if( strcmp( argv[ii], "--gaussian_blur" ) == 0 ){
      ops.gaussian_blur = std::stod( argv[ii+1] );
    }

    if( strcmp( argv[ii], "--add_grains" ) == 0 ){
      auto split = my_utility::str_split( std::string( argv[ii+1] ), ',' );
      if( split.size() == 7 ){

	ops.grains.resize( split.size() );
	for( unsigned int i=0; i<split.size(); i++){
	  ops.grains[i] = std::stod( split[i] );
	}

      } else {
	std::cerr << "--add_grains: wrong number of parameters" << std::endl;
      }
    }


    if( strcmp( argv[ii], "--skip" ) == 0 ){
      auto split = my_utility::str_split( std::string( argv[ii+1] ), ',' );      

      if( split.size() % 2 == 0 ){
	for( unsigned int ii=0; ii<split.size(); ii++){
	  if( ii%2==0 ){
	    ops.to_skip.push_back( split[ii] );
	  } else {
	    ops.strides.push_back( std::stoi( split[ii] ) );
	  }
	}	
      } else {
	std::cerr << "incorrect format of --skip parameter" << std::endl;
      }
      
    }

    
    if( strcmp( argv[ii], "--post_processing" ) == 0 ){
      for( int jj=ii+1; jj<argc; jj++){
	ops.post_processing.push_back( std::string( argv[jj] ) );
      }
    }    
    
  }
    
}




class cmd_callback : public sp_callback {

public:
  cmd_callback( surface_projection &sp,
		batch_creation &bc_,
		std::string prefix,
		std::string summary_format_,
		double gaussian_noise_,
		double gaussian_blur_,
		std::vector<double> grains_ ) :
    sp_callback( sp ),
    fn_prefix( prefix ),
    bc( bc_ ),
    counter (0),
    summary_format( summary_format_ ),
    gaussian_noise( gaussian_noise_ ),
    gaussian_blur( gaussian_blur_ ),
    grains( grains_ )
  {
    summary.open( fn_prefix + "_summary.txt" );
    if( summary_format != "human" && summary_format != "csv" ){
      std::cerr << "unkown summary format, using \"csv\"" << std::endl;
      summary_format = "csv";
    }
    if( summary_format == "csv" ){
      summary << "filename";
      for( unsigned int ii=0; ii<surface_projection::parameter_names.size(); ii++){
	summary << "," << surface_projection::parameter_names.at(ii);
      }
      summary << std::endl;
    }
  }

  ~cmd_callback(){
    summary.close();
  }
  
  bool operator()( std::vector< std::vector<double>::iterator > pars ){

    std::stringstream ssfn;
    ssfn << fn_prefix << "_" << std::setfill('0') << std::setw(5) << counter << ".png";
    std::string fn = ssfn.str();//fn_prefix + "_" + std::to_string( counter ) + ".png";

    sp.validate_uc_scaling();
    sp.update_geometry();
    sp.compute_projection();


    unsigned char *img = sp.get_image( invert, "LIN" );

    if( grains.size() > 0 ){

      unsigned char *noise = new unsigned char[sp.get_width()*sp.get_height()]();
      noise = image_manipulation::create_grains( sp.get_width(), sp.get_height(), grains[0], grains[1], grains[2], grains[3], grains[4] );
      if( grains[5] > 0 ){
	image_manipulation::gaussian_blur( noise, sp.get_width(), sp.get_height(), grains[5] );
      }
      image_manipulation::add_images( img, noise, 1.0-grains[6], grains[6], sp.get_width(), sp.get_height(), img );
      delete[]( noise );
      
    }    

    if( gaussian_blur > 0 ){
      image_manipulation::gaussian_blur( img, sp.get_width(), sp.get_height(), gaussian_blur );
    }

    if( gaussian_noise > 0 ){
      image_manipulation::gaussian_noise( img, sp.get_width(), sp.get_height(), gaussian_noise );
    }

    image_manipulation::write_png( fn, img, sp.get_width(), sp.get_height() );
    delete[] (img);
    
    // write the summary file entry

    summary << fn;;
    if( summary_format == "human" ){
      summary << "    ";
      for( auto it : surface_projection::parameter_names ){
	summary << it << "=" << sp.get_parameter( it ) << " ";
      }
    } else if ( summary_format == "csv" ){
      for( unsigned int ii=0; ii<surface_projection::parameter_names.size(); ii++){
	summary << "," << sp.get_parameter( surface_projection::parameter_names.at(ii) );
      }
    }
    summary << std::endl;
    
      
    counter++;
    return true;
  }

  long counter;
  
  batch_creation &bc;
  std::string fn_prefix;
  std::ofstream summary;
  std::string summary_format;
  double gaussian_blur;
  double gaussian_noise;
  std::vector<double> grains;
};


int main( int argc, char* argv[] ){

  if( argc == 1 ){
    print_help();
    return EXIT_SUCCESS;
  }

  
  global_settings gs ( "global_settins.conf" );
  cmd_options ops;
  surface_projection sp ( gs );
  batch_creation bc ( sp, time(NULL) );
  
  try {
    parse_cmd_options( argc, argv, bc, ops );
    if( ops.post_processing.size() == 0 && !check_and_set_fixed_pars( ops, sp ) ){
      return EXIT_FAILURE;
    }
  } catch ( std::string s ){
    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }




  // check if we're only post-processing than nothing needs to be
  // initialized
  if( ops.post_processing.size() > 0 ){
    
    // get an array of filenames
    for( unsigned int ii=0; ii<ops.post_processing.size(); ii++){
      unsigned char *img;
      int width, height;
      image_manipulation::read_png( ops.post_processing[ii], &img, &width, &height );

      if( ops.grains.size() > 0 ){
	unsigned char *noise = new unsigned char[width*height]();
	noise = image_manipulation::create_grains( width, height, ops.grains[0], ops.grains[1], ops.grains[2], ops.grains[3], ops.grains[4] );
	if( ops.grains[5] > 0 ){
	  image_manipulation::gaussian_blur( noise, width, height, ops.grains[5] );
	}
	image_manipulation::add_images( img, noise, 1.0-ops.grains[6], ops.grains[6], width, height, img );
	delete[]( noise );	
      }      

      if( ops.gaussian_blur > 0 ){
	image_manipulation::gaussian_blur( img, width, height, ops.gaussian_blur );
      }
      
      if( ops.gaussian_noise > 0 ){
	image_manipulation::gaussian_noise( img, width, height, ops.gaussian_noise );
      }
      image_manipulation::write_png( ops.post_processing[ii], img, width, height );
      delete[] (img);
    }
    return EXIT_SUCCESS;
  }
  

  // get the functor

  cmd_callback cmdcb ( sp, bc, ops.fn_prefix, ops.sum_format, ops.gaussian_noise, ops.gaussian_blur, ops.grains );
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


      
      // create the skip string to pass to set_random_parameters
      std::vector<std::string> skip;
      for( unsigned int kk=0; kk<ops.strides.size(); kk++){
	
	if( counter % ops.strides[kk] && counter>0 ){
	  skip.push_back( ops.to_skip[kk] );
	}
	
      }
            
      if( !bc.set_random_parameters( ops.quadratic, skip ) ){
	f_out << counter << " " << " caught invalid_parameter_exception" << std::endl;
      }

      counter++;
      cmdcb( dummy_it );
    }

    f_out.close();
    
  }
  
  return EXIT_SUCCESS;
}


