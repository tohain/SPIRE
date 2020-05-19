#include <iostream>
#include <fstream>
#include <algorithm>

#include "surface_projection.hpp"
#include "img_out.hpp"



int main(int argc, char* argv[])
{


  //get the surface projection done
  double progress;// = new double;
  std::string status;// = new char[20];


  surface_projection sp (progress, status);

  //the entire cube
  sp.set_slice_width( 1 );
  sp.set_slice_height( 0 );


  int res = 150;
  
  sp.set_n_points_x( res );
  sp.set_n_points_y( res );
  sp.set_n_points_z( res );
  
  sp.set_type( 2 );
  sp.set_surface_level( 0.0 ) ;

  sp.update_geometry();
  sp.update_containers();


  
  std::ofstream volumes_out ("vols_p.dat");



  for( double level_set = -3.1; level_set <= 3.1; level_set+= 0.1 ){

    sp.set_surface_level( level_set );
    sp.compute_projection();
    sp.compute_volume();
    auto vols = sp.get_channel_volumes();

    volumes_out << level_set;
    for( auto it : vols ){
      volumes_out << " " << it;
    }
    volumes_out << std::endl;
  }

  

    
}
