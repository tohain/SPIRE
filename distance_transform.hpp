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


/** \brief This class computes the euclidean distance map of a structure
 *
 * This class computes the distance transformation of a picture. It is
 * based on an algorithm published by Pedro F. Felzenszwalb and Daniel
 * P. Huttenlocher, T HEORY OF C OMPUTING , Volume 8 (2012), pp. 415–428.
 * 10.4086/toc.2012.v008a019
 *
 * It is based on a one-dimensional distance transform. For higher
 * dimensional cases, the 1D transform is applied repeatetly. See the
 * publication for more details
 */
class distance_transform {

public:
  
  /// Constructor for a 1D transform
  distance_transform( std::vector<double> data, int n, double pixsize = 1);

  /// Constructor for a 2D transform
  distance_transform( std::vector<double> data, int n, int k, double pixsize_x = 1, double pixsize_y = 1 );

  /// Constructor for a 3D transform
  distance_transform( std::vector<double> data, int n, int k, int l, double pixsize_x = 1, double pixsize_y = 1, double pixsize_z = 1 );  

  /// Destructor, freeing memory
  ~distance_transform();
  
  /// Performs the distance transformation
  void compute_distance_map();


  /// get distance map
  std::vector<double> get_distance_map() const ;
  
private:

  /// 1D distance transform
  std::vector<double> do_distance_transform( std::vector<double> &grid, double pixsize = 1 );

  /// evaluates the function on the given data
  std::vector<double> eval_grid_function( std::vector<double> &data );
  
private:

  /// The image to transform
  std::vector<double> img;
  
  /// The array holding the distance map
  std::vector<double> map;

  /// The dimension of the object to transform
  unsigned int dim;
  /// The number of pixels of the object in each direction (x,y,z)
  unsigned int size[3];
  /// The size of a pixel in length units for each direction (x,y,z)
  double pix_size[3];
};


#endif
