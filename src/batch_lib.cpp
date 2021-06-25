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

#include "batch_lib.hpp"

/**
 * Constructor
 */
batch_creation::batch_creation( surface_projection &sp_, int seed ) :
  sp( sp_ ), generator( std::mt19937() ) {

  generator.seed( seed );
}


/**
 * Parses a string of comma separated values or a range into a vector
 * of doubles
 */ 
void batch_creation::add_parameter( std::string parameter, std::string in ){

  auto parts_range = my_utility::str_split( in, ':' );
  auto parts_singles = my_utility::str_split( in, ',' );  

  ops.parameters.push_back( parameter );
  ops.values.push_back( std::vector<double> (0, 0) );
  
  // processing a range
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

    int size = std::fabs( (end-start) / step );
    ops.values[ ops.values.size() -1 ].reserve( size + 1 );


    
    if ( start < end ){
      for( double ii=start; ii<=end; ii+=step ){
	ops.values[ ops.values.size() -1 ].push_back( ii );
      }
    } else {
      for( double ii=start; ii>=end; ii-=step ){
	ops.values[ ops.values.size() -1 ].push_back( ii );
      }
    }

    
  } else if ( parts_singles.size() >= 1 ){
    // processing a list of single values
    ops.values[ ops.values.size() -1 ].resize(parts_singles.size(), 0);

    for( size_t ii=0; ii<parts_singles.size(); ii++ ){
      ops.values[ ops.values.size() -1 ][ii] = std::stod( parts_singles[ii] );
    }
  }
  
}


/**
 * Adds a set of values to the loop
 */
void batch_creation::add_parameter( std::string parameter, std::vector<double> vals ){
  ops.parameters.push_back( parameter );
  ops.values.push_back( vals );
}
    


std::string batch_creation::create_filename( std::vector< std::vector<double>::iterator > pars, std::string delimiter ){

  std::stringstream filename;

  /*
  if( prefix != "" ){
    filename << prefix;
  }

  if( prefix != "" ){
    filename << delimiter;
  }

  for( auto it : ops.membranes ){
    if( it >= 0 )
      filename << "+";
    filename << it;
  }
  filename << delimiter;
  
  for( auto it : ops.filled_channels ){
    filename << it;
  }
  filename << delimiter;
  
  for( auto it : pars ){
    filename << *it << delimiter;
  }

  filename << suffix;
  */  

  return filename.str();
}


/** 
 * this function checks, if the direction of the vector hkl has been
 * visited already
 * 
 * returns true, if we already checked that direction
 */
template <class T>
bool batch_creation::check_hkl_doubles( T h, T k, T l, std::set< std::vector<int> > &visited ){

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


void batch_creation::print_cur_parameters( std::vector< std::vector<double>::iterator > its, std::ostream &out ){
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
long batch_creation::do_loop( sp_callback &callback ){

  // the position of the miller indeces in values array
  unsigned int h_pos=std::numeric_limits<unsigned int>::max();
  unsigned int k_pos=std::numeric_limits<unsigned int>::max();
  unsigned int l_pos=std::numeric_limits<unsigned int>::max();  
  unsigned int nr_miller_provided = 0;
  
  // get a set of iterators
  std::vector< std::vector<double>::iterator > its;
  for( unsigned int ii=0; ii<ops.values.size(); ii++ ){
    // find the position of the miller indeces, will need them later
    if( ops.parameters[ii] == surface_projection::parameter_names[9] ){
      h_pos = ii;
      nr_miller_provided++;
    }
    if( ops.parameters[ii] == surface_projection::parameter_names[10] ){
      k_pos = ii;
      nr_miller_provided++;
    }
    if( ops.parameters[ii] == surface_projection::parameter_names[11] ){
      l_pos = ii;
      nr_miller_provided++;
    }    
    its.push_back( ops.values[ii].begin() );
  }

  // find the "top-level" miller index, if this is resetted, we mus
  // reset the visited directions, too!
  // top level is the first appearing in values array
  unsigned int hkl_first = std::min( h_pos, std::min(k_pos, l_pos) );
  
  long computed = 0;
  
  std::set<std::vector<int> > hkl_visited;
  
  // start looping. The last element is incremented most frequently.
  if( its.size() > 0 ){
  
    while( true ){
      
      // check if any iterator is at its end, starting at last
      bool stop = false;
      for( int ii=its.size()-1; ii>=0; ii--){ 
	
	if( its[ii] == ops.values[ii].end() ){
	  its[ii] = ops.values[ii].begin();

	  if( ii == 0 ){
	    stop = true;
	  } else {	  
	    its[ii-1]++;
	  }

	  // reset visited directions
	  if( ii == hkl_first ){
	    hkl_visited.clear();
	  }
	}

	// if no miller index is provided, reset the visited map after
	// each parameter change
	if( nr_miller_provided == 0 ){
	  hkl_visited.clear();
	}
	
      }

      //if first element is at its end, we're done
      if( stop ){
	break;
      }
      
      // we can now set the state of the sp object by calling something like this
      bool caught_exception = false;
      for( unsigned int ii=0; ii<ops.values.size(); ii++ ){
	try{
	  sp.set_parameter( ops.parameters[ii], *(its[ii]) );
	} catch ( invalid_parameter_exception e ){
	  print_cur_parameters( its );
	  std::cout << ": ";
	  caught_exception = true;
	  std::cout << e.msg << std::endl;
	} catch ( std::string s ){
	  caught_exception = true;
	  std::cerr << s << std::endl;
	}
      }
      
      // call the callback only if we didn't run into any issues given
      // the current set of parameters
      if( !caught_exception ){
	// check if this direction was computed already	
	if( !check_hkl_doubles<double>( sp.get_h(),
					sp.get_k(),
					sp.get_l(),
					hkl_visited ) ){
	  if ( !callback( its ) ){
	    // if callback returns false quit and stop
	    break;	    
	  }
	  computed++;
	}
      }
      
      // increment the last element
      its[ its.size() - 1 ]++;

      callback.inc_counter();
    }
    
  }  

  return computed;
}

/**
 * clears all parameters
 */
void batch_creation::reset_parameters(){
  ops.parameters.clear();
  ops.values.clear();
}

/**
 * computes the total number of possible combinations. Some might be
 * skipped though, due to orientation or other criteria
 */
long batch_creation::total_combinations(){

  // to measure progress, get total amount of elements
  long total_elements = 1;
  for( unsigned int ii=0; ii<ops.parameters.size(); ii++){
    total_elements *= ops.values[ii].size();
  }
  return total_elements;
}

/**
 * Draws and sets random parameters to the sp object. Only values for
 * parameters present in the ops.parameters array are generated. The
 * first two values in the corresponding values array are taken as
 * lower and upper boundary of uniformly distributed samples
 *
 * \param[in] qudratic if set to true, only qudartic projections are sampled (height set to width)
 * \param[out] Returns true if all parameters set without errors, false otherwise
 */
bool batch_creation::set_random_parameters( bool quadratic ){

  bool caught_exception = false;


  
  for( unsigned int ii=0; ii<ops.parameters.size(); ii++){

    double val;

    // get the distribution first
    if( int( ops.values[ii][0] ) == 0 ){   // uniform distribtion      
      std::uniform_real_distribution<double> dist (ops.values[ii][1], ops.values[ii][2] );
      val = dist( generator );            
    } else if ( int( ops.values[ii][0] == 1 ) ) {      
      std::normal_distribution<double> dist ( ops.values[ii][1], ops.values[ii][2] );      
      val = dist( generator );
    }

      
      try {
	sp.set_parameter( ops.parameters[ii], val );
      } catch (invalid_parameter_exception e ){
	// we probably tried to set (hkl)=(000). In this case the
	// value is set to 001, so nothing inherently "bad" happens,
	// fixing it is not worth the effort right now, but keep it in
	// mind.
	caught_exception = true;
      }
  }

  if( quadratic ) sp.set_slice_height( sp.get_slice_width() ); // set the slice height to slice width

  return !caught_exception;
}
