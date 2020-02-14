#ifndef ITERABLE_VOXEL_H
#define ITERABLE_VOXEL_H


#include <unordered_set>
#include <vector>

class iterable_voxel {
  
public:
  
  iterable_voxel( int id_, int x_, int y_, int z_ ) : id(id_), x(x_), y(y_), z(z_) {}
  
  int operator()(){
    return id;
  }
  
  iterable_voxel r();
  iterable_voxel l();
  iterable_voxel u();
  iterable_voxel d();
  iterable_voxel f();
  iterable_voxel b();  
  
  
  /// +/-x = f/b, +/-y = l/r,
  /// +/-z=u/d
  std::unordered_set<int> get_6_neighbors();
  
  /// +/-x = f/b, +/-y = l/r,
  /// +/-z=u/d
  std::unordered_set<int> get_18_neighbors();
  
  /// +/-x = f/b, +/-y = l/r,
  /// +/-z=u/d
  std::unordered_set<int> get_26_neighbors();
  
  
  
  /// +/-x = f/b, +/-y = l/r,
  /// +/-z=u/d
  std::vector< std::pair<int, std::string> > get_6_neighbors_dir();
  
  /// +/-x = f/b, +/-y = l/r,
  /// +/-z=u/d
  std::vector< std::pair<int, std::string> > get_18_neighbors_dir();
  
  /// +/-x = f/b, +/-y = l/r,
  /// +/-z=u/d
  std::vector< std::pair<int, std::string> > get_26_neighbors_dir();
  
    
  std::vector<int> get_circle( std::string dir );
  
private:
  /// the id of the currently considered particles
  int id;
  /// int the dimensions of the grid
  int x, y, z;
};




#endif
