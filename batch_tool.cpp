#include "surface_projection.hpp"
#include "distance_transform.hpp"
#include "img_out.hpp"
#include "homotopic_thinning.hpp"

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include <cstring>

#include "percolation_analysis.hpp"

#include "auxiliary.hpp"

#include <png.h>



typedef struct {

  std::vector<std::string> parameters;
  std::vector< std::vector<double> > values;

  std::vector<double> membranes;
  std::vector<double> filled_channels;
  std::vector<double> resolution;

  std::string fn_prefix;  
} options;


std::vector<std::string> iter_par_names = { "struct_types", "uc_scale_ab", "uc_scale_c", "surface_level", "slice_thickness", "slice_height", "slice_width", "slice_position", "miller_h", "miller_k", "miller_l" };


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

	    << std::endl

	    << "sweep parameters, each parameter can be assign a comma separated list," << std::endl
	    << "or a range given as start:stop:step" << std::endl;

  for( auto it : iter_par_names ){
    std::cout << "   --" << it << std::endl;
  }
  
  std::cout << std::endl << std::endl;
  
}

template <class T>
std::vector<T> parse_string( std::string in ){

  auto parts_range = my_utility::str_split( in, ':' );
  auto parts_singles = my_utility::str_split( in, ',' );  

  std::vector<T> out;
  
  if( parts_range.size() > 1 ){

    double start = std::stod(parts_range[0]);
    double end = std::stod(parts_range[1]);
    double step = std::stod(parts_range[2]);    
    if( step == 0 ){
      throw std::string ("invalid range in " + in );
    }   
    if( start < end && step < 0 ){
      throw std::string ("invalid range in " + in );
    }
    if( start > end && step > 0 ){
      throw std::string ("invalid range in " + in );
    }    

    if ( start < end ){
      for( double ii=start; ii<=end; ii+=step ){
	out.push_back( ii );
      }
    } else {
      for( double ii=start; ii>=end; ii+=step ){
	out.push_back( ii );
      }
    }
    
  } else if ( parts_singles.size() >= 1 ){
    for( auto it : parts_singles ){
      out.push_back( static_cast<T>( std::stod( it ) ) );
    }
  }
  
  return out;
}


void parse_cmd_options( int argc, char* argv[], options &ops ){
  
  char *found_it = NULL;
  
  for( int ii=1; ii<argc; ii++) {
    
    for( unsigned int jj=0; jj<iter_par_names.size(); jj++){


      /*
       * check all the iterable parameters
       */
      found_it = strstr( argv[ii], iter_par_names[jj].c_str() );
      
      if( found_it != NULL ){
	ops.parameters.push_back( iter_par_names[jj] );
	ops.values.push_back( parse_string<double> ( std::string( argv[ii+1] ) ) );
      }
    }


    /*
     * now check the "fixed" parameters
     */
    if( strcmp( argv[ii], "--membranes" ) == 0 ){
      ops.membranes = parse_string<double> ( std::string( argv[ii+1] ) );
    }
    if( strcmp( argv[ii], "--filled_channels" ) == 0 ){
      ops.filled_channels = parse_string<double> ( std::string( argv[ii+1] ) );
    }

    if( strcmp( argv[ii], "--resolution" ) == 0 ){
      ops.resolution = parse_string<double> ( std::string( argv[ii+1] ) );
    }

    if( strcmp( argv[ii], "--fn_prefix" ) == 0 ){
      ops.fn_prefix = std::string( argv[ii+1] );
    }    
    
  }
    
}


/**
 * sets the according parameter of the surface projection module. Just
 * a wrapper around the setters/getters of \ref surface_projection
 */
void set_parameter( surface_projection &sp, std::string par, double val ){

  if( par == iter_par_names[0] ){
    sp.set_type( static_cast<int>( val ) );
  }
  
  if( par == iter_par_names[1] ){
    sp.set_uc_scale_ab( val );
  }
  if( par == iter_par_names[2] ){
    sp.set_uc_scale_c( val );
  }  
  if( par == iter_par_names[3] ){
    sp.set_surface_level( val );
  }


  if( par == iter_par_names[4] ){
    sp.set_slice_thickness( val );
  }
  if( par == iter_par_names[5] ){
    sp.set_slice_height( val );
  }
  if( par == iter_par_names[6] ){
    sp.set_slice_width( val );
  }
  if( par == iter_par_names[7] ){
    sp.set_slice_position( val );
  }  

  if( par == iter_par_names[8] ){
    sp.set_h( static_cast<int>( val ) );
  }
  if( par == iter_par_names[9] ){
    sp.set_k( static_cast<int>( val ) );
  }
  if( par == iter_par_names[10] ){
    sp.set_l( static_cast<int>( val ) );
  }
  
}


std::string create_filename( std::string prefix, std::string suffix, std::vector< std::vector<double>::iterator > pars, std::string delimiter = "_" ){

  std::stringstream filename;

  if( prefix != "" ){
    filename << prefix;
  }
  
  for( auto it : pars ){
    if( prefix != "" ){
      filename << delimiter;
    }
    filename << *it;
  }

  filename << suffix;
  
  return filename.str();
}


/** 
 * this function checks, if the direction of the vector hkl has been
 * visited already
 * 
 * returns true, if we already checked that direction
 */
template <class T>
bool check_hkl_doubles( T h, T k, T l, std::set< std::vector<int> > &visited ){

  // the absolute value only works for cubic symmetry!!!!!
  int h_int = static_cast<int> ( h );
  int k_int = static_cast<int> ( k );
  int l_int = static_cast<int> ( l );  

  int gcd = std::abs( my_utility::gcd_euclid( std::vector<int> {std::abs(h_int),
								std::abs(k_int),
								std::abs(l_int)} ) );

  std::vector<int> reduced_vec = { h_int / gcd, k_int / gcd, l_int / gcd };
  
  if( visited.find( reduced_vec ) == visited.end() ){
    // did not find the vector, new direction!
    visited.insert( reduced_vec );
    return false;
  } else {

    return true;
  }
}


void print_cur_parameters( options &ops, std::vector< std::vector<double>::iterator > its, std::ostream &out = std::cout){
  for(unsigned int ii=0; ii<ops.parameters.size(); ii ++){
    out << ops.parameters[ii] << "=" << *(its[ii]) << " ";
  }
}

/**
 * this function loops over all possible combinations of all iterators
 * in the \ref parameters vector in \ref ops.
 *
 * It keeps track of which hkl direction already has been computed and
 * avoids multiple combination of these directions
 *
 * returns the number of processed elements
 */
int do_loop( options &ops, surface_projection &sp, void (*callback)(surface_projection &sp, std::vector< std::vector<double>::iterator > &pars, options &ops) ){

  // get a set of iterators
  std::vector< std::vector<double>::iterator > its;
  for( unsigned int ii=0; ii<ops.values.size(); ii++ ){
    its.push_back( ops.values[ii].begin() );
  }

  // to measure progress, get total amount of elements
  long total_elements = 1;
  for( unsigned int ii=0; ii<ops.parameters.size(); ii++){
    total_elements *= ops.values[ii].size();
  }
  std::cout << "processing a total of " << total_elements << " configurations" << std::endl;

  long counter = 0;
  long computed = 0;
  
  std::set<std::vector<int> > hkl_visited;
  
  // start looping. The last element is incremented most frequently.
  if( its.size() > 0 ){
  
    while( true ){
      
      // check if any iterator is at its end, starting at last
      for( unsigned int ii=its.size()-1; ii>0; ii--){ 
	
	if( its[ii] == ops.values[ii].end() ){
	  its[ii] = ops.values[ii].begin();
	  its[ii-1]++;
	}
	
      }

      //if first element is at its end, we're done
      if( its[0] == ops.values[0].end() ){
	break;
      }
      
      // we can now set the state of the sp object by calling something like this
      bool caught_exception = false;
      for( unsigned int ii=0; ii<ops.values.size(); ii++ ){
	try{
	  set_parameter( sp, ops.parameters[ii], *(its[ii]) );
	} catch ( invalid_parameter_exception e ){
	  print_cur_parameters( ops, its );
	  std::cout << ": ";
	  caught_exception = true;
	  std::cout << e.msg << std::endl;
	} catch ( std::string s ){
	  caught_exception = true;
	  std::cerr << s << std::endl;
	}
      }


      // check the hkl direction
      auto h_it = std::find( ops.parameters.begin(), ops.parameters.end(), iter_par_names[8] );
      auto k_it = std::find( ops.parameters.begin(), ops.parameters.end(), iter_par_names[9] );
      auto l_it = std::find( ops.parameters.begin(), ops.parameters.end(), iter_par_names[10] );      
      auto h_pos = std::distance( ops.parameters.begin(), h_it );
      auto k_pos = std::distance( ops.parameters.begin(), k_it );
      auto l_pos = std::distance( ops.parameters.begin(), l_it );
      bool visited = true; // in case we have hkl=(0,0,0) we won't
			   // compute it, additional safety
      if( !(*its[h_pos]==0 && *its[k_pos]==0 && *its[l_pos] ==0 ) ){
	visited = check_hkl_doubles<double>( *its[h_pos], *its[k_pos], *its[l_pos], hkl_visited );
      }

      
      // call the callback only if we didn't run into any issues given
      // the current set of parameters
      if( !caught_exception && !visited ){
	callback( sp, its, ops );
	computed++;
      }
      
      // increment the last element
      its[ its.size() - 1 ]++;

      counter++;
      // might not want to call it each time we increment the counter,
      // but since computing the projection is the heavy task I'd
      // argue we won't create a bottle neck here
      if( counter % 1 == 0){
	std::cout << "computing " << counter << "/" << total_elements
		  << " (" << std::setw(5) << std::setprecision(4)
		  << ((float)counter / total_elements) * 100.0
		  << "\%)\r" << std::flush;
      }
    }
    
  }  

  std::cout << std::endl << "done!" << std::endl;
  return computed;
}


/** 
 * check the most critical parameters for consistency and apply them
 */
bool check_and_set_fixed_pars( options &ops, surface_projection &sp ){

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

/**
 * the function which is called after each parameter update
 */
void update_and_compute( surface_projection &sp,
			 std::vector< std::vector<double>::iterator > &pars,
			 options &ops ){

  sp.update_geometry();
  sp.compute_projection();
  sp.write_png( create_filename( ops.fn_prefix, ".png", pars ) );
  
}

int main( int argc, char* argv[] ){

  if( argc == 1 ){
    print_help();
  }
  
  double progress;
  std::string status;
  surface_projection sp( progress, status );


  options ops;
  try {
    parse_cmd_options( argc, argv, ops );
    if( !check_and_set_fixed_pars( ops, sp ) ){
      return EXIT_FAILURE;
    }
  } catch ( std::string s ){
    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }


  long elements_created = do_loop( ops, sp, update_and_compute );
  std::cout << "made " << elements_created << " projections" << std::endl;
  
  return EXIT_SUCCESS;
}
