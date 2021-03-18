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

#ifndef BATCH_TOOL_LIB
#define BATCH_TOOL_LIB

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


  /// define a call-back functor, implemented by the caller
class sp_callback {
public:
  sp_callback( surface_projection &sp_ ) : sp( sp_ ), external_counter (0) {}

  /**
   * purely virtual function implemented in descending classes. The
   * iterators provide the current state of the loop
   */ 
  virtual bool operator()( std::vector< std::vector<double>::iterator > pars ) = 0;

  void inc_counter(){
    external_counter++;
  }

  void set_img_mode ( bool invert_, std::string scaling_ ){
    invert = invert_;
    scaling = scaling_;
  }
  
  
  surface_projection &sp;
  // keeps track of how many parameter configurations has been
  // considered, including the skipped ones
  unsigned external_counter;
  bool invert;
  std::string scaling;
};

class batch_creation {

public:
  
  typedef struct {

    std::vector<std::string> parameters;
    std::vector< std::vector<double> > values;
    
    std::string fn_prefix;  
  } batch_options;

  /// Constructor
  batch_creation( surface_projection &sp, int seed );
  
  /// Parses a string of comma separated values or a range into a
  /// vector of doubles
  void add_parameter( std::string parameter, std::string in );

  /// Adds a list of values to the loop
  void add_parameter( std::string parameter, std::vector<double> vals );

  /// Deletes all parameters to iterate over
  void reset_parameters();

  /// gets the total amount of possible combinations
  long total_combinations();
  
  /// creates a filename from the parameters
  std::string create_filename( std::vector< std::vector<double>::iterator > pars,
			       std::string delimiter = "_" );

  /// checks if equvialent hkl values has already been computed
  template <class T>
  bool check_hkl_doubles( T h, T k, T l, std::set< std::vector<int> > &visited );

  /// executes the loop; updates parameters and calls the callback after each update
  long do_loop( sp_callback &callback );
  
  /// print the currently set parameters to the ostream
  void print_cur_parameters( std::vector< std::vector<double>::iterator > its,
			     std::ostream &out = std::cout);

  /// get random parameters
  bool set_random_parameters( bool quadratic );

  batch_options ops;
  surface_projection &sp;

  boost_rand rng;
};





#endif
