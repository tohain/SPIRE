#ifndef HOMOTOPIC_THINNING_H
#define HOMOTOPIC_THINNING_H


#include "iterable_voxel.hpp"

#include <unordered_set>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cmath>
#include <cstdlib>

/*
 * This class provides functionality to compute the skeleton of a 3D
 * Object. It is an implementation of the algorithms described in
 *
 * Pudney (1998): https://doi.org/10.1006/cviu.1998.0680
 * Bertrand, Malandain (1994): https://doi.org/10.1016/0167-8655(94)90046-9
 * Vincent (1991) : https://doi.org/10.1117/12.45227
 */
template <class T>
class homotopic_thinning {


public:
  /// Constructor
  homotopic_thinning ( int nx, int ny, int nz, std::vector<T> image_, std::vector<float> dmap_ );  

  /// Computes and returns the skeleton of the channel
  std::unordered_set<int> find_channel_skeleton( T ch_id );
  

private:
  /// Returns all points on the interface of the channel
  std::unordered_set<int> get_exterior_points( T channel_id, int m );

  /// Checks if a point can be deleted without chaning the topology
  bool point_deletable( int id, int m, int n );

  /// Computes all connected components around a voxel
  std::vector< std::unordered_set<int> > get_connected_components( int center, int m, bool img );
  
  /// The image to reduce
  std::vector<T> image;

  /// the distance map
  std::vector<float> dmap;

  /// Dimension of the image
  int n_points_x, n_points_y, n_points_z;

};


#endif
