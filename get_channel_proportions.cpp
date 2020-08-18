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
  sp.set_slice_length( 1 );
  sp.set_slice_height( 1 );
  sp.set_slice_width( 1 );
  sp.set_slice_position( 0 );

  // unit unit cell
  sp.set_uc_scale_ab( 1 );
  sp.set_uc_scale_c( 1 );
  sp.update_a();

  sp.set_h( 0 );
  sp.set_k( 0 );
  sp.set_l( 1 );
  sp.set_orientation_from_hkl();
  sp.compute_uc_dim_in_orientation();
  auto dim = sp.get_uc_dim_in_orientation();
  
  sp.set_slice_length( dim[0] );
  sp.set_slice_height( dim[1] );
  sp.set_slice_width( dim[2] );  
  
  // at lower resolution the membrane isn't that thick.
  int res = 200;
  
  sp.set_n_points_x( res );
  sp.set_n_points_y_to_unitcell();
  sp.set_n_points_z_to_unitcell();

  // iterate over all available surfaces
  for(unsigned int ii=0; ii < sp.get_surface_choices().size(); ii++){

    // open the outputfile
    std::ofstream volumes_out ("vols_"+sp.get_surface_choices()[ii]+".dat");
    
    std::cout << sp.get_surface_choices()[ii] << std::endl;

    // set up surface
    sp.set_type( ii );
    sp.set_surface_level( 0.0 ) ;
    
    sp.update_geometry();
    sp.update_containers();
  


    // this is somewhat arbitrary and depends on the surface where to start and end
    double begin = -3.5, end=3.52, step=0.01;

    std::vector< std::vector<double> > vol_container;
    std::vector<double> levels;

    // get volumes
    for( double level_set = begin; level_set <= end; level_set+= step ){
      
      sp.set_surface_level( level_set );
      sp.compute_projection();
      
      sp.compute_volume();

      vol_container.push_back( sp.get_channel_volumes() );
      levels.push_back( level_set );
    }
      

      
    

    /*
     * output
     */
    
    // human readable first
    for( unsigned int jj=0; jj<levels.size(); jj++){
      volumes_out << levels[jj] << " ";
      for( auto it : vol_container[jj] ){
	volumes_out << " " << it;
      }
      volumes_out << std::endl;
    }
    
    
    
    // header code
    
    volumes_out << std::endl << std::endl << std::endl <<"const std::map<double, double> "
		<< sp.get_surface_choices()[ii] << "_SURFACE_VOL = {";
    
    
    for( unsigned int jj=0; jj<levels.size(); jj++){
      volumes_out << "{" << vol_container[jj][0] + (vol_container[jj][1]*0.5)
		  << ", " << levels[jj] << "}";

      if( jj+1 < levels.size() ){
	volumes_out << ", ";
      }

    }


    volumes_out << "};" << std::endl << std::endl << std::endl;
    

    volumes_out << "const std::map<double, double> " 
	      << sp.get_surface_choices()[ii] << "_SURFACE_VOL_INV = {";  

    for( unsigned int jj=0; jj<levels.size(); jj++){
      
      volumes_out << "{" << levels[jj] << ", "
		  << vol_container[jj][0] + (vol_container[jj][1]*0.5) << "}"; 
      if( jj+1 < levels.size() ){
	volumes_out << ", ";
      }
    }
    
    volumes_out << "};";
    
  }

}


