#include "qt_sp.hpp"


sp_qt::sp_qt( double &p, std::string &s) : surface_projection(p, s){

  connect( this, SIGNAL( geometry_changed() ), this, SLOT( update_geometry_() ) );
  
}

sp_qt::~sp_qt(){

}



void sp_qt::update_geometry_(){
  set_orientation_from_hkl();
  update_periodicity_length();  
  update_geometry();
  update_containers();
}

void sp_qt::change_surface_type( int ind ){

  //backup channel proportion
  double vol_prop = get_channel_prop();

  set_type( ind );

  set_channel_vol_prop( vol_prop );
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


void sp_qt::change_slice_width( double val ){
  set_slice_width( val );
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_slice_position( double val ){
  set_slice_height( val );
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_membranes( std::vector<double> val ){
  set_membranes( val );
  emit parameter_changed();
}

void sp_qt::compute_projection(){
  
  //get the points in the slice  
  emit send_message( "Busy", 0 );
  set_up_points();
  
  //reset grid
  memset( grid.data(), 0, sizeof(int) * grid.size() );
  
  //get grid
  set_grid();
  
  //get projection
  //emit send_message( "Computing Projection", 0 );
  memset( projection.data(), 0, sizeof(double) * projection.size() );
  project_grid();

  emit projection_changed();
  emit send_message( "Ready", 0 );
  
}


void sp_qt::do_something(){
  std::cout << "Doing something" << std::endl;
}


void sp_qt::copy_parameters( sp_qt *source ){


  set_h( source->get_h() );
  set_k( source->get_k() );
  set_l( source->get_l() );  

  set_ntucs( source->get_ntucs() );
  set_a( source->get_a() );
  set_membranes( source->get_membranes() );
  set_channel_vol_prop( source->get_channel_prop() );

  set_slice_width( source->get_slice_width() );
  set_slice_height( source->get_slice_height() );

  set_type( source->get_type() );

  emit geometry_changed();
  emit parameter_changed();
}




void sp_qt::update_measurements(){

  emit send_message( "Busy", 3 );
  
  update_geometry();

  emit send_message( "Computing projection", 1 );
  
  compute_projection();

  emit send_message( "Computing volumes", 1 );
  
  compute_volume();

  emit send_message( "Computing areas", 1 );
  
  compute_surface_area();

  emit send_message( "Computing networks", 1 );
  
  compute_channel_network();
  
  emit send_message( "Finished computing measurements", 1 );

  emit send_message( "Ready", 3 );

  emit measurements_updated();
}



void sp_qt::save_grid( QString fn ){

  print_grid( fn.toStdString() );
  
}

void sp_qt::save_topological_network( int id, QString fn ){

  print_topological_network( id, fn.toStdString() );
  
}

void sp_qt::save_surface_points( int id, QString fn ){

  print_channel_surface_points( id, fn.toStdString() );
  
}
