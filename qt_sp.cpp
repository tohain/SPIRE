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
  set_n_points_y_to_unitcell( );
  std::cout << "n_points=(" << n_points_x << "," << n_points_y <<","<< n_points_z << ")" << std::endl;
  emit geometry_changed();
  emit parameter_changed();
}


void sp_qt::change_uc_scale_ab( double ab ){
  
  set_uc_scale_ab( ab );

  // if we have membranes of cubic symmetries, change c scaling as
  // well, since it needs to be cubic
  if( surface_choices[type] != "Wurtzite" ){
    set_uc_scale_c( ab );
  }

  update_a();  
  
  emit geometry_changed();
  emit parameter_changed();
}


void sp_qt::change_uc_scale_c( double c ){

  //only allow if we have other than cubic symmetries

  if( surface_choices[type] == "Wurtzite" ){
  
    set_uc_scale_c( c );
    update_a();
  
    emit geometry_changed();
    emit parameter_changed();
  }
}


  /*
void sp_qt::change_uc_size_c( double val ){
  std::vector<double> tmp_size = get_a();
  set_a( tmp_size[0], tmp_size[1], val );  
  emit geometry_changed();
  emit parameter_changed();
}
  */

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

void sp_qt::change_slice_length( double val ){
  set_slice_length( val );
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_slice_height( double val ){
  set_slice_height( val );
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

  compute_slice_size();
  
  //get the points in the slice
  emit set_status( 0, 1 );
  //emit send_message( "Busy", 0 );
  set_up_points();
  
  //reset grid
  memset( grid.data(), 0, sizeof(short) * grid.size() );
  
  //get grid
  set_grid();
  
  //get projection
  //emit send_message( "Computing Projection", 0 );
  memset( projection.data(), 0, sizeof(float) * projection.size() );
  project_grid();
  
  emit projection_changed();
  emit set_status( 0, 0 );
  //emit send_message( "Ready", 0 );
  
}


void sp_qt::do_something(){
  std::cout << "Doing something" << std::endl;
}


void sp_qt::copy_parameters( sp_qt *source ){


  set_h( source->get_h() );
  set_k( source->get_k() );
  set_l( source->get_l() );  

  set_ntucs( source->get_ntucs() );


  set_uc_scale_ab( source->get_uc_scale_ab() );
  set_uc_scale_c( source->get_uc_scale_c() );
  update_a();
  
  set_membranes( source->get_membranes() );
  set_channel_vol_prop( source->get_channel_prop() );

  set_slice_width( source->get_slice_width() );
  set_slice_position( source->get_slice_position() );

  set_type( source->get_type() );

  emit geometry_changed();
  emit parameter_changed();
}




void sp_qt::update_measurements(){
  
  emit set_status( 1, 1 );
  
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
