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
  sp_callback( surface_projection &sp_ ) : sp( sp_ ){}

  /**
   * purely virtual function implemented in descending classes. The
   * iterators provide the current state of the loop
   */ 
  virtual void operator()( std::vector< std::vector<double>::iterator > pars ) = 0;
  
  surface_projection &sp;
};

class batch_creation {

public:
  
  typedef struct {

    std::vector<std::string> parameters;
    std::vector< std::vector<double> > values;
    
    std::string fn_prefix;  
  } batch_options;

  /// Constructor
  batch_creation( surface_projection &sp );

  /// Parses a string of comma separated values or a range into a
  /// vector of doubles
  void add_parameter( std::string parameter, std::string in );

  /// Adds a list of values to the loop
  void add_parameter( std::string parameter, std::vector<double> vals );


  /// creates a filename from the parameters
  std::string create_filename( std::vector< std::vector<double>::iterator > pars,
			       std::string delimiter = "_" );

  /// checks if equvialent hkl values has already been computed
  template <class T>
  bool check_hkl_doubles( T h, T k, T l, std::set< std::vector<int> > &visited );

  /// executes the loop; updates parameters and calls the callback after each update
  int do_loop( sp_callback &callback );
  
  /// print the currently set parameters to the ostream
  void print_cur_parameters( std::vector< std::vector<double>::iterator > its,
			     std::ostream &out = std::cout);


private:
  batch_options ops;
  surface_projection &sp;
  
};





#endif
