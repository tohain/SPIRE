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
  inline int operator()(){
    return id;
  }


  inline void set( int id_ ){
    id = id_;
  }

  

  /// one voxel to the right
  inline iterable_voxel r(){
    int m = id % (x*z);
    if( (x*z) - z <= m &&
	m < (x*z) ){
      return iterable_voxel( id - (( x - 1) * z), x, y, z );
    } else {
      return iterable_voxel( id + z, x, y, z );
    }
  }

  /// one voxel to the left
  inline iterable_voxel l(){
    int m = id % (x*z);
    if( 0 <= m && m < z ){
      return iterable_voxel( id + ((x - 1) * z), x, y, z );
    } else {
      return iterable_voxel(id - z, x, y, z);
    }
  }
  

  /// one voxel up
  inline iterable_voxel u(){
    if(id % z == z - 1){
      return iterable_voxel(id-(z-1), x, y, z);
    } else {
      return iterable_voxel(id+1, x, y, z);
    }  
  }

  /// one voxel down
  inline iterable_voxel d(){
    if( id % z == 0 ){
      return iterable_voxel( id+(z-1), x, y, z);
    } else {
      return iterable_voxel( id - 1, x, y, z);
    }
  }
  

  /// one voxel backwards
  inline iterable_voxel b(){
    if( id < x * z ){
      return iterable_voxel( id + ( x * z * (y-1) ), x, y, z);
    } else {
      return iterable_voxel( id - (x * z), x, y, z );
    }
  }

  /// one voxel forward
  inline iterable_voxel f(){
    if( id >= (y-1) * x * z ){
      return iterable_voxel( id - ( x * z * (y-1) ), x, y, z );
    } else {
      return iterable_voxel(id + (x * z), x, y, z);
    }
  }


    /// one voxel to the right
  inline void rr(){
    int m = id % (x*z);
    if( (x*z) - z <= m &&
	m < (x*z) ){
      id -= (( x - 1) * z);
    } else {
      id += z;
    }
  }

  /// one voxel to the left
  inline void ll(){
    int m = id % (x*z);
    if( 0 <= m && m < z ){
      id += ((x - 1) * z);
    } else {
      id -= z;
    }
  }
  

  /// one voxel up
  inline void uu(){
    if(id % z == z - 1){
      id-=(z-1);
    } else {
      id+=1;
    }  
  }

  /// one voxel down
  inline void dd(){
    if( id % z == 0 ){
      id+=(z-1);
    } else {
      id -= 1, x, y, z;
    }
  }
  

  /// one voxel backwards
  inline void bb(){
    if( id < x * z ){
      id+= x * z * (y-1);
    } else {
      id-= (x * z);
    }
  }

  /// one voxel forward
  inline void ff(){
    if( id >= (y-1) * x * z ){
      id -= x * z * (y-1);
    } else {
      id += x * z;
    }
  }


  /// just the next voxel, moves to next col, row automatically
  
  
  inline void get_6_neighbors( std::vector<int> &nbs ){
    
    nbs.resize (6, 0);
    
    nbs[0] = u()();
    nbs[1] = d()();
  
    nbs[2] = l()();
    nbs[3] = r()();
    
    nbs[4] = f()();
    nbs[5] = b()();    
  }
  
  
  inline void get_18_neighbors(std::vector<int> &nbs){
    
    get_6_neighbors( nbs );
    nbs.resize (18, 0);
    
    nbs[6] = f().l()();
    nbs[7] = f().r()();
    nbs[8] = b().l()();
    nbs[9] = b().r()();
    
    nbs[10] = l().u()();
    nbs[11] = l().d()();
    nbs[12] = r().u()();
    nbs[13] = r().d()();
    
    nbs[14] = f().u()();
    nbs[15] = f().d()();
    nbs[16] = b().u()();
    nbs[17] = b().d()();    
  }    
  
  
  
  inline void get_26_neighbors( std::vector<int> &nbs){
    
    get_18_neighbors( nbs );
    nbs.resize(26, 0);
    
    nbs[18] = u().f().l()();
    nbs[19] = u().f().r()();
    nbs[20] = u().b().l()();
    nbs[21] = u().b().r()();
    
    nbs[22] = d().f().l()();
    nbs[23] = d().f().r()();
    nbs[24] = d().b().l()();
    nbs[25] = d().b().r()();
    
  }
  
  /*
    iterable_voxel r();
    iterable_voxel l();
    iterable_voxel u();
  iterable_voxel d();
  iterable_voxel f();
  iterable_voxel b();  
  
  
  /// Returns all neighbours with 6 connectivity (only main axes)
  inline std::unordered_set<int> get_6_neighbors();
  
  /// Returns all neighbours with 18 connectivity (main axis, and in plane diagonals)
  inline std::unordered_set<int> get_18_neighbors();
  
  /// Returns all neighbours with 26 connectivity (all diagonals and main axes)
  inline std::unordered_set<int> get_26_neighbors();
    */

    
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
