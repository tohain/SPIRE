#include "batch_lib.hpp"

/**
 * Constructor
 */
batch_creation::batch_creation( surface_projection &sp_ ) :
  sp( sp_ ) {
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
      for( double ii=start; ii>=end; ii+=step ){      
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
int batch_creation::do_loop( sp_callback &callback ){

  // the position of the miller indeces in values array
  unsigned int h_pos=std::numeric_limits<unsigned int>::max();
  unsigned int k_pos=std::numeric_limits<unsigned int>::max();
  unsigned int l_pos=std::numeric_limits<unsigned int>::max();  
  
  // get a set of iterators
  std::vector< std::vector<double>::iterator > its;
  for( unsigned int ii=0; ii<ops.values.size(); ii++ ){
    // find the position of the miller indeces, will need them later
    if( ops.parameters[ii] == surface_projection::parameter_names[8] ){
      h_pos = ii;
    }
    if( ops.parameters[ii] == surface_projection::parameter_names[9] ){
      k_pos = ii;
    }
    if( ops.parameters[ii] == surface_projection::parameter_names[10] ){
      l_pos = ii;
    }    
    its.push_back( ops.values[ii].begin() );
  }

  // find the "top-level" miller index, if this is resetted, we mus
  // reset the visited directions, too!
  // top level is the first appearing in values array
  unsigned int hkl_first = std::min( h_pos, std::min(k_pos, l_pos) );
  
  // to measure progress, get total amount of elements
  long total_elements = 1;
  for( unsigned int ii=0; ii<ops.parameters.size(); ii++){
    total_elements *= ops.values[ii].size();
  }
  std::cout << "processing a total of " << total_elements << " possible configurations" << std::endl;

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

	  // reset visited directions
	  if( ii == hkl_first ){
	    hkl_visited.clear();
	  }
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
	  callback( its );
	  computed++;
	}
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






