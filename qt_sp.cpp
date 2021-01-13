/* Projection tool - compute planar projection of triply periodic
 * minimal surfaces 
 * Copyright (C) 2020 Tobias Hain
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


#include "qt_sp.hpp"


sp_qt::sp_qt( double &p, std::string &s) : surface_projection(p, s){

  connect( this, SIGNAL( geometry_changed() ), this, SLOT( update_geometry_() ) );
  
}

sp_qt::~sp_qt(){

}



void sp_qt::update_geometry_(){

  try {
    update_geometry();
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
    set_orientation_from_hkl();
    //compute_uc_dim_in_orientation();
    compute_smallest_uc();

}

void sp_qt::change_surface_type( int ind ){

  //backup channel proportion
  double vol_prop = get_channel_prop();

  set_type( ind );

  set_channel_vol_prop( vol_prop );

  emit geometry_changed();
  emit parameter_changed();
}


void sp_qt::change_z_points( int val ){
  try {
    set_n_points_z( val );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_xy_points( int val ){
  try {
    set_n_points_x( val );
    set_n_points_y_to_unitcell( );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
  emit geometry_changed();
  emit parameter_changed();
}



void sp_qt::change_uc_scale_ab( double ab ){
  
  set_uc_scale_ab( ab );

  // if we have membranes of cubic symmetries, change c scaling as
  // well, since it needs to be cubic
  if( surface_choices[type].substr(0,8) != "Wurtzite" ){
    set_uc_scale_c( ab );
  }

  update_a();  
  
  emit geometry_changed();
  emit parameter_changed();
}


void sp_qt::change_uc_scale_c( double c ){

  //only allow if we have other than cubic symmetries

  if( surface_choices[type].substr(0,8) == "Wurtzite" ){
  
    set_uc_scale_c( c );
    update_a();
  
    emit geometry_changed();
    emit parameter_changed();
  }
}



void sp_qt::change_vol_prop( double val ){
  set_channel_vol_prop( val );
  emit parameter_changed();
}

void sp_qt::change_lvl_set( double val ){
  set_surface_level( val );
  emit parameter_changed();
}



void sp_qt::change_hkl( int _h, int _k, int _l ){

  if( !(_h == 0 && _k == 0 && _l == 0 ) ){
    set_h( _h ); set_k( _k ); set_l( _l );
    emit geometry_changed();
    emit parameter_changed();
  } else {
    emit send_message( "At least one parameter of hkl must be non-zero", 1 );
    emit parameter_changed();
  } 
}


void sp_qt::change_slice_thickness( double val ){
  try {
    set_slice_thickness( val );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_slice_width( double val ){
  try {
    set_slice_width( val );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_slice_height( double val ){
  try {
    set_slice_height( val );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }    
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_slice_position( double val ){
  set_slice_position( val );
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_membranes( std::vector<double> val ){
  set_membranes( val );
  emit parameter_changed();
}

void sp_qt::compute_projection(){

  update_containers(); 
  
  //get the points in the slice
  emit set_status( 0, 1 );
  //emit send_message( "Busy", 0 );
  set_up_points();
  
  //reset grid
  memset( grid.data(), 0, sizeof(short) * grid.size() );
  
  //get grid
  set_grid();

  // apply channel colors
  for( unsigned int ii=0; ii < channel_filled.size(); ii++){

    if( ii % 2 == 0 ){
      //channel, if it is zero, do nothing, otherwise mark it
      if( channel_filled[ii] != 0 ){
	set_channel_color( ii + 1, 1 );
      }
    } else {
      //membrane, if it is one, do nothing, otherwise mark it
      if( channel_filled[ii] != 1 ){
	set_channel_color( ii + 1, 0 );
      }      
    }

  }
  
  //get projection
  //emit send_message( "Computing Projection", 0 );
  memset( projection.data(), 0, sizeof(float) * projection.size() );
  project_grid();
  
  emit projection_changed();
  emit set_status( 0, 0 );
  //emit send_message( "Ready", 0 );
  
}


void sp_qt::do_something(){

  std::cout << "doing something" << std::endl;

  dt.set_parameters( grid, std::vector<unsigned int> {n_points_x, n_points_y, n_points_z},
		     std::vector<double> {dx, dy, dz}, true);  
  dt.compute_distance_map();
  dt.compute_max_radius_covering();
  

  //print_channel_width_distribution( "channel_width_dist.dat" );
  print_max_rad_transform_dist( "mrct.dat" );

  std::vector<float> mrct = dt.get_max_radius_covering();
  std::vector<float> dmap = dt.get_distance_map();  
  
  // print out slices
  double min = 0;
  double mrct_max = *std::max_element( mrct.begin(), mrct.end() );
  double dmap_max = *std::max_element( dmap.begin(), dmap.end() );  

  unsigned char *img_data = new unsigned char[ n_points_x * n_points_y ];
  unsigned char *dmap_data = new unsigned char[ n_points_x * n_points_y ];  
  
  for( int zz = 0; zz < n_points_z; zz++){
    std::ofstream data_out ("mrct_" + std::to_string(zz) + ".dat" );
    std::ofstream dmap_out ("dmap_" + std::to_string(zz) + ".dat" );
    std::ofstream grid_out ("grid_" + std::to_string(zz) + ".dat" );
    for( int yy = 0; yy < n_points_y; yy++){      
      for( int xx = 0; xx < n_points_x; xx++ ){

	int ind = zz + xx * n_points_z + yy*n_points_z*n_points_x;

	double mrct_val = (mrct[ind] / mrct_max) * 255.0;
	unsigned char val_c_mrct = static_cast<unsigned char> (mrct_val);

	double dmap_val = (dmap[ind] / dmap_max) * 255.0;
	unsigned char val_c_dmap = static_cast<unsigned char> (dmap_val);
	
	int img_ind = xx + yy*n_points_x;
	
	img_data[img_ind] = val_c_mrct;
	dmap_data[img_ind] = val_c_dmap;

	data_out << std::setprecision(3) << std::setw(6) << mrct[ind] << " ";
	dmap_out << std::setprecision(3) << std::setw(6) << sqrt(dmap[ind]) << " ";
	grid_out << std::setw(1) << std::setprecision(1) << grid[ind] << " ";
      }
      data_out << std::endl;
      dmap_out << std::endl;
      grid_out << std::endl;
    }

    data_out.close();
    dmap_out.close();
    
    write_image( "layer_" + std::to_string(zz) + ".ppm", img_data, n_points_x, n_points_y );
    write_image( "dmap_" + std::to_string(zz) + ".ppm", dmap_data, n_points_x, n_points_y );    
    
    
  }

  std::cout << "done" << std::endl;
  
}



void sp_qt::copy_parameters( sp_qt *source ){

  // orientation
  set_h( source->get_h() );
  set_k( source->get_k() );
  set_l( source->get_l() );  

  // update orientation
  set_orientation_from_hkl();  

  // slice
  set_slice_thickness( source->get_slice_thickness() );
  set_slice_height( source->get_slice_height() );
  set_slice_width( source->get_slice_width() );
  set_slice_position( source->get_slice_position() );

  // points
  set_n_points_x( source->get_width() );
  set_n_points_y( source->get_height() );
  set_n_points_z( source->get_depth() );
  
  // uc scale
  set_uc_scale_ab( source->get_uc_scale_ab() );
  set_uc_scale_c( source->get_uc_scale_c() );

  update_a();
  
  // structure
  set_surface_level( source->get_surface_level() );
  set_type( source->get_type() );
  set_membranes( source->get_membranes() );

  try {
    update_geometry();
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }

  
  emit geometry_changed();
  emit parameter_changed();
}




void sp_qt::update_measurements( QString what ){  

  // set "measurement" (first digit) to "busy" (second digits)
  emit set_status( 1, 1 );

  try {
    update_geometry();
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }  

  emit send_message( "Computing projection", 1 );
  compute_projection();

  emit send_message( "Computing " + what, 1 );  
  
  if( what == "Volumes" ){
    compute_volume();
  }

  if( what == "Areas" ){
    compute_surface_area();
  }

  if( what == "Networks" ){
    compute_channel_network();
  }

  if( what == "Percolation" ){
    compute_percolation_threshold();
  }

  emit send_message("Measurement Done");
  emit set_status( 1, 0 );

  emit measurements_updated();
}


void sp_qt::save_grid( QString fn ){

  print_grid( fn.toStdString() );

  emit send_message( "Saved grid!", 1);
}

void sp_qt::save_topological_network( int id, QString fn ){

  print_topological_network( id, fn.toStdString() );

  emit send_message( "Saved network!", 1);
}

void sp_qt::save_surface_points( int id, QString fn ){

  print_channel_surface_points( id, fn.toStdString() );

  emit send_message( "Saved membrane!", 1);
}

void sp_qt::set_slice_dim_to_uc(){

  set_slice_to_uc( 0.1 );

  /*
  if( uc_dim_in_orientation[0] > 0 )
    set_slice_width( uc_dim_in_orientation[0] );
  if( uc_dim_in_orientation[1] > 0 )
    set_slice_height( uc_dim_in_orientation[1] );
  if( uc_dim_in_orientation[2] > 0 )
    set_slice_thickness( uc_dim_in_orientation[2] );
  */

  emit geometry_changed();
  emit parameter_changed();

}


void sp_qt::change_channel_color( int id, int val ){
  set_channel_color( id, val );
}
