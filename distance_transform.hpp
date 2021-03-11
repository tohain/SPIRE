/* Projection tool - compute planar projection of triply periodic
 * minimal surfaces 
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


#ifndef SP_DISTTRANS_I
#define SP_DISTTRANS_I

#include <limits>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <map>
#include <unordered_set>


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
 * ToDo: remove implicit integral type <-> floating point type conversions
 *
 *
 */
template <class T, class M = float>
class distance_transform {

public:

  /// Constructor doing nothing. This is very dangerous to use I
  /// think, please use \ref set_parameters immediately after!
  distance_transform();
  
  /// Constructor for a 1D transform
  distance_transform( std::vector<T> data, int n, bool periodic, double pixsize = 1.0);

  /// Constructor for a 2D transform
  distance_transform( std::vector<T> data, int n, int k, bool periodic, double pixsize_x = 1.0, double pixsize_y = 1.0 );

  /// Constructor for a 3D transform
  distance_transform( std::vector<T> data, int n, int k, int l, bool periodic, double pixsize_x = 1.0, double pixsize_y = 1.0, double pixsize_z = 1.0 );  

  /// Destructor, freeing memory
  ~distance_transform();
  
  /// Performs the distance transformation
  void compute_distance_map();

  /// computes the largest radius covering transform
  void compute_max_radius_covering();
  
  /// get distance map
  std::vector<M> get_distance_map() const ;

  /// returns the max. radius cover transform
  std::vector<M> get_max_radius_covering() const ;

  /// Update the parameters
  void set_parameters( std::vector<T> data,
		       std::vector<unsigned int> dim,
		       std::vector<double> pixsize,
		       bool periodic );
  
  /// outputs the map to console
  void print_map( std::ostream &out = std::cout ) const;
  
private:

  /// 1D distance transform
  std::vector<M> do_distance_transform( std::vector<M> &grid, bool periodic, double pixsize = 1 );

  /// evaluates the function on the given data
  template <class U>
  std::vector<M> eval_grid_function( std::vector<U> &data, double max );

  /// returns all pixels which are within a radius r of the provided pixel
  std::unordered_set<unsigned int> get_voxels_in_ball( double r, int id );

  
  /// The image to transform
  std::vector<T> img;
  
  /// The array holding the distance map
  std::vector<M> map;

  /// Container for maximum radius covering transform
  std::vector<M> mrct;
  std::vector<bool> marked;
  
  /// The dimension of the object to transform
  unsigned int dim;
  /// The number of pixels of the object in each direction (x,y,z)
  unsigned int size[3];
  /// The size of a pixel in length units for each direction (x,y,z)
  double pix_size[3];

  /// Is the data periodic?
  bool is_periodic;
};



#endif
