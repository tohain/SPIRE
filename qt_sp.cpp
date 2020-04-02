#include "qt_sp.hpp"


sp_qt::sp_qt( double &p, std::string &s) : surface_projection(p, s){


  connect( this, SIGNAL( parameter_changed() ), this, SLOT( compute_projection() ) );

  connect( this, SIGNAL( geometry_changed() ), this, SLOT( update_geometry_() ) );
  
}

sp_qt::~sp_qt(){

}



void sp_qt::update_geometry_(){
  update_geometry();
  update_containers();
}

void sp_qt::change_surface_type( int ind ){
  set_type( ind );
  emit parameter_changed();
}


void sp_qt::change_ntucs( int val ){
  set_ntucs( val );
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_z_points( int val ){
  set_n_points_z( val );
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_xy_points( int val ){
  set_n_points_x( val );
  set_n_points_y( val );  
  emit geometry_changed();
  emit parameter_changed();
}


void sp_qt::change_uc_size( double val ){
  set_a( val );
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_vol_prop( double val ){
  set_channel_vol_prop( val );
  emit parameter_changed();
}



void sp_qt::compute_projection(){
    
  //get the points in the slice  
  emit status_updated( "Computing the points" );
  set_up_points();
  
  //reset grid
  memset( grid.data(), 0, sizeof(int) * grid.size() );
  
  //get grid
  emit status_updated( "Setting up the grid" );
  set_grid();
  
  //get projection
  emit status_updated( "Computing Projection" );
  memset( projection.data(), 0, sizeof(double) * projection.size() );
  project_grid();

  emit projection_changed();
  emit status_updated( "Ready" );
}
