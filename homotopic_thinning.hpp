#ifndef HOMOTOPIC_THINNING_H
#define HOMOTOPIC_THINNING_H


#include "iterable_voxel.hpp"

#include <unordered_set>
#include <algorithm>
#include <map>
#include <unordered_map>

class homotopic_thinning {


public:

  homotopic_thinning ( int nx, int ny, int nz, std::vector<int> image_, std::vector<double> dmap_ );  

  std::unordered_set<int> find_channel_skeleton( int ch_id );
  


private:

  std::unordered_set<int> get_exterior_points( int channel_id, int m );

  bool point_deletable( int id, int m, int n );
  
  std::vector< std::unordered_set<int> > get_connected_components( int center, int m, bool img );
  
  /// the image
  std::vector<int> image;

  /// the distance map
  std::vector<double> dmap;
  
  int n_points_x, n_points_y, n_points_z;


};


#endif
