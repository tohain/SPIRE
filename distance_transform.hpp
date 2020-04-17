/*
 * This file is part of the Surface Projection Tool
 *
 * Author: T. Hain
 * 07/2019
 *
 * Licence: TBA
 */


#ifndef SP_DISTTRANS_I
#define SP_DISTTRANS_I

#include <limits>
#include <vector>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>

/** \brief This class computes the euclidean distance map of a structure
 *
 * This class computes the distance transformation of a picture. It is
 * based on an algorithm published by Pedro F. Felzenszwalb and Daniel
 * P. Huttenlocher, THEORY OF COMPUTING , Volume 8 (2012), pp. 415–428.
 * 10.4086/toc.2012.v008a019
 *
 * It is based on a one-dimensional distance transform. For higher
 * dimensional cases, the 1D transform is applied repeatetly. See the
 * publication for more details.
 *
 * This class comes as a template, so distance maps of ints, floats,
 * doubles,[...] can be computed without first converting the
 * vector. Internally it will be converted to a double in the \ref
 * eval_grid_function function. So make sure it can be converted to a
 * floating point number
 * 
 */
template <class T, class M = float>
class distance_transform {

public:

  /// Constructor doing nothing. This is very dangerous to use I
  /// think, please use \ref set_parameters immediately after!
  distance_transform();
  
  /// Constructor for a 1D transform
  distance_transform( std::vector<T> data, int n, double pixsize = 1);

  /// Constructor for a 2D transform
  distance_transform( std::vector<T> data, int n, int k, double pixsize_x = 1, double pixsize_y = 1 );

  /// Constructor for a 3D transform
  distance_transform( std::vector<T> data, int n, int k, int l, double pixsize_x = 1, double pixsize_y = 1, double pixsize_z = 1 );  

  /// Destructor, freeing memory
  ~distance_transform();
  
  /// Performs the distance transformation
  void compute_distance_map();

  /// get distance map
  std::vector<M> get_distance_map() const ;

  /// Update the parameters
  void set_parameters( std::vector<T> data,
		       std::vector<unsigned int> dim,
		       std::vector<double> pixsize );
  
  /// outputs the map to console
  void print_map() const;
  
private:

  /// 1D distance transform
  std::vector<M> do_distance_transform( std::vector<M> &grid, double pixsize = 1 );

  /// evaluates the function on the given data
  template <class U>
  std::vector<M> eval_grid_function( std::vector<U> &data, double max );
    
  /// The image to transform
  std::vector<T> img;
  
  /// The array holding the distance map
  std::vector<M> map;

  /// The dimension of the object to transform
  unsigned int dim;
  /// The number of pixels of the object in each direction (x,y,z)
  unsigned int size[3];
  /// The size of a pixel in length units for each direction (x,y,z)
  double pix_size[3];
};



#endif
