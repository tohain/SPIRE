#ifndef ITERABLE_VOXEL_H
#define ITERABLE_VOXEL_H


#include <unordered_set>
#include <vector>
#include <string>

/*
 * A quick & dirt implementation of a pixel, which can iterate thorugh
 * a 3D image. Accounts for periodic boundary conditions. The notation
 * chosen in this class is +/-x = f/b, +/-y = l/r, +/-z=u/d, where
 * directions are given as the stringtripel (xyz), e.g. (111) means
 * f->l->u;
 */
class iterable_voxel {
  
public:
  /// Constructor
  iterable_voxel( int id_, int x_, int y_, int z_ ) : id(id_), x(x_), y(y_), z(z_) {}

  /// Overload () operator for easy access
  int operator()(){
    return id;
  }
  
  iterable_voxel r();
  iterable_voxel l();
  iterable_voxel u();
  iterable_voxel d();
  iterable_voxel f();
  iterable_voxel b();  
  
  
  /// Returns all neighbours with 6 connectivity (only main axes)
  std::unordered_set<int> get_6_neighbors();
  
  /// Returns all neighbours with 18 connectivity (main axis, and in plane diagonals)
  std::unordered_set<int> get_18_neighbors();
  
  /// Returns all neighbours with 26 connectivity (all diagonals and main axes)
  std::unordered_set<int> get_26_neighbors();
  
  
  /// Returns all neighbours and their directions with 6 connectivity (only
  /// main axes)
  std::vector< std::pair<int, std::string> > get_6_neighbors_dir();
  
  /// Returns all neighbours and their directions with 18 connectivity (only
  /// main axes)  
  std::vector< std::pair<int, std::string> > get_18_neighbors_dir();
  
  /// Returns all neighbours and their directions with 26 connectivity (only
  /// main axes)
  std::vector< std::pair<int, std::string> > get_26_neighbors_dir();
  
  /// Returns the particles in a circle around the current particle
  /// perpendicular to dir
  std::vector<int> get_circle( std::string dir );
  
private:
  /// the id of the currently considered particles
  int id;
  /// int the dimensions of the grid
  int x, y, z;
};




#endif
