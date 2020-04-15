#include "qt_gui.hpp"

void GUI::set_up_ui(){

  // the widgets for the tabs
  controls = new QTabWidget( this );
  controls_basic = new QWidget();
  controls_save = new QWidget();
  controls_measurement = new QWidget();

  controls->addTab( controls_basic, "Basic" );
  controls->addTab( controls_measurement, "Measurements" );
  controls->addTab( controls_save, "Save Image" );

  // stats tab
  detailled_stats = new QLabel( controls_save );

  // save tab
  save_grid_control = new QPushButton( "Save grid", controls_save );
  save_surface_points_control = new QPushButton( "Save membrane points", controls_save );
  save_topological_network_control = new QPushButton( "Save network", controls_save );
  path_prefix_control = new QLineEdit( controls_save );
  
  // buttons for controls_basic
  button_quit = new QPushButton("Quit", controls_basic);
  button_save = new QPushButton("Save", controls_basic);
  button_render = new QPushButton("Compute", controls_basic);
  button_measure = new QPushButton("Measure", controls_basic);  

  /*
   * structure control
   */ 
  ntucs_control = new QT_labeled_obj<QSpinBox> ( "Unit cells", controls_basic );
  ntucs_control->object()->setRange(1, 50);
  
  uc_size_control = new QT_labeled_obj<QDoubleSpinBox> ( "Unit cell size", controls_basic );
  uc_size_control->object()->setRange(0.01, 100);
  uc_size_control->object()->setSingleStep(0.01);    
  
  channel_prop_control = new QT_labeled_obj<QDoubleSpinBox> ( "Volume proportions", controls_basic );
  channel_prop_control->object()->setRange(0, 1);
  channel_prop_control->object()->setSingleStep(0.01);  
  
  surface_type_control = new QT_labeled_obj<QComboBox> ( "Surface type", controls_basic );
  std::vector<std::string> sfc_types = sp->get_surface_choices();
  for(unsigned int ii=0; ii<sfc_types.size(); ii++){
    surface_type_control->object()->insertItem( ii, QString( sfc_types.at(ii).c_str() ) );
  }

  /*
   * Resolution control
   */
  xy_points_control = new QT_labeled_obj<QSpinBox> ( "XY resolution", controls_basic );
  xy_points_control->object()->setRange(1, 3000);
 
  z_points_control = new QT_labeled_obj<QSpinBox> ( "Z resolution", controls_basic );  
  z_points_control->object()->setRange(1, 3000);

  pix_size_indicator = new QLabel("");
  pix_size_indicator->setWordWrap( true );
  pix_size_indicator->setAlignment(Qt::AlignCenter);
  pix_size_indicator->setTextFormat(Qt::PlainText);

  invert_control = new QCheckBox( controls_basic );
  invert_control->setText( "Invert colors" );

  autoupdate_control = new QCheckBox( controls_basic );
  autoupdate_control->setText( "Autoupdate" );  
  
  
  //Set up the draw area
  draw_area = new QLabel();
  draw_area->setAlignment(Qt::AlignVCenter | Qt::AlignHCenter);
  draw_area->setMinimumSize(500,500);
  draw_area->setSizePolicy(QSizePolicy::MinimumExpanding,
                     QSizePolicy::MinimumExpanding);



  /*
   * Slice control 
   */

  slice_width_control = new QT_labeled_obj<QDoubleSpinBox>( "Slice width", controls_basic);
  slice_width_control->object()->setRange(0, 1);
  slice_width_control->object()->setSingleStep(0.01);
  
  slice_position_control = new QT_labeled_obj<QDoubleSpinBox>( "Slice height", controls_basic);    
  slice_position_control->object()->setRange(0, 1);
  slice_position_control->object()->setSingleStep(0.01);  
  
  miller_h_control = new QT_labeled_obj<QSpinBox> ("h", controls_basic );
  miller_h_control->object()->setRange(-100, 100);
  
  miller_k_control = new QT_labeled_obj<QSpinBox> ("k", controls_basic );
  miller_k_control->object()->setRange(-100, 100);
  
  miller_l_control = new QT_labeled_obj<QSpinBox> ("l", controls_basic );  
  miller_l_control->object()->setRange(-100, 100);

  /*
   * membrane control
   */

  membranes_control = new QTableWidget( controls_basic );
  membranes_control->insertColumn( 0 );
  membranes_control->insertColumn( 0 );

  QStringList table_header;
  table_header << "Distance" << "Width";
  membranes_control->setHorizontalHeaderLabels( table_header );
  
  add_membrane_control = new QPushButton ( "Add" );
  rm_membrane_control = new QPushButton( "Remove" );

  membranes_label = new QLabel( "Membranes", controls_basic );
  
  
  /*
   * Layout
   */
  
  //set up spacers
  v_spacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
  h_spacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);


  // status bar
  status_bar = new QStatusBar();
  status_bar_status_m = new QLabel();
  status_bar_status_m->setAlignment(Qt::AlignCenter);
  status_bar_status_p = new QLabel();
  status_bar_status_p->setAlignment(Qt::AlignCenter);  
  
  status_bar_vols = new QLabel();
  status_bar_areas = new QLabel();  
  status_bar_pixs = new QLabel();
  status_bar_mins = new QLabel();
  status_bar->addPermanentWidget( status_bar_mins );  
  status_bar->addPermanentWidget( status_bar_vols );
  status_bar->addPermanentWidget( status_bar_areas );  
  status_bar->addPermanentWidget( status_bar_pixs );  
  status_bar->addPermanentWidget( status_bar_status_m );
  status_bar->addPermanentWidget( status_bar_status_p );  

  
  //the main layout of the form
  main_layout = new QVBoxLayout( this );
  
  sub_main_layout = new QHBoxLayout();
  sub_main_layout->addWidget( controls);
  sub_main_layout->addWidget( draw_area );

  status_bar_layout = new QHBoxLayout();
  status_bar_layout->addWidget( status_bar );
  
  main_layout->addLayout( sub_main_layout );
  main_layout->addLayout( status_bar_layout );
  
  // the layout of the buttons
  buttons_layout = new QHBoxLayout();
  // the layout of the tabs
  controls_basic_layout = new QVBoxLayout( controls_basic );
  controls_save_layout = new QVBoxLayout( controls_save );
  controls_measurement_layout = new QVBoxLayout( controls_measurement );  



  controls_measurement_layout->addWidget( detailled_stats );

  controls_save_layout->addWidget( path_prefix_control );
  controls_save_layout->addWidget( save_grid_control );
  controls_save_layout->addWidget( save_surface_points_control );
  controls_save_layout->addWidget( save_topological_network_control );
    

  
  
  structure_settings = new QHBoxLayout();
  resolution_settings = new QHBoxLayout();
  slice_settings = new QHBoxLayout();
  membrane_settings = new QHBoxLayout();
  membrane_buttons_layout = new QVBoxLayout();
  
  structure_settings->addLayout( ntucs_control->layout() );
  structure_settings->addLayout( uc_size_control->layout() );
  structure_settings->addLayout( channel_prop_control->layout() );
  structure_settings->addLayout( surface_type_control->layout() );  

  resolution_settings->addLayout( xy_points_control->layout() );
  resolution_settings->addLayout( z_points_control->layout() );
  resolution_settings->addWidget( pix_size_indicator );
  resolution_settings->addWidget( invert_control );
  resolution_settings->addWidget( autoupdate_control );

  slice_settings->addLayout( slice_width_control->layout() );
  slice_settings->addLayout( slice_position_control->layout() );
  slice_settings->addItem( h_spacer );
  slice_settings->addLayout( miller_h_control->layout() );
  slice_settings->addLayout( miller_k_control->layout() );
  slice_settings->addLayout( miller_l_control->layout() );    

  membrane_buttons_layout->addWidget( membranes_label );
  membrane_buttons_layout->addWidget( add_membrane_control );
  membrane_buttons_layout->addWidget( rm_membrane_control );
  membrane_settings->addLayout( membrane_buttons_layout );
  membrane_settings->addWidget( membranes_control );
  
  controls_basic_layout->addLayout( structure_settings );
  controls_basic_layout->addItem( v_spacer );
  controls_basic_layout->addLayout( slice_settings );
  controls_basic_layout->addItem( v_spacer );
  controls_basic_layout->addLayout( membrane_settings );
  controls_basic_layout->addItem( v_spacer );  
  controls_basic_layout->addLayout( resolution_settings );
  controls_basic_layout->addItem( v_spacer );
  controls_basic_layout->addLayout( buttons_layout );
    
  buttons_layout->addWidget( button_render );
  buttons_layout->addWidget( button_measure );
  buttons_layout->addWidget( button_save );
  buttons_layout->addWidget( button_quit );





  
}

void GUI::set_up_signals_and_slots(){

  
  /*
   * update parameters
   */

  // structure type
  connect(surface_type_control->object(), QOverload<int>::of(&QComboBox::currentIndexChanged),
	  sp, &sp_qt::change_surface_type);
  // xy resolution
  connect(xy_points_control->object(), QOverload<int>::of(&QSpinBox::valueChanged),
	  sp, &sp_qt::change_xy_points );
  // z resolution
  connect(z_points_control->object(), QOverload<int>::of(&QSpinBox::valueChanged),
	  sp, &sp_qt::change_z_points );
  // number of unit cells
  connect( ntucs_control->object(), QOverload<int>::of(&QSpinBox::valueChanged),
	   sp, &sp_qt::change_ntucs );
  // unit cell size
  connect( uc_size_control->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_uc_size );
  // channel proportion
  connect( channel_prop_control->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_vol_prop );

  // slice width
  connect( slice_width_control->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_slice_width );
  // slice position
  connect( slice_position_control->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_slice_position );
  // miller hkl

  connect( miller_h_control->object(), QOverload<int>::of(&QSpinBox::valueChanged),
	   this, &GUI::change_orientation );
  connect( miller_k_control->object(), QOverload<int>::of(&QSpinBox::valueChanged),
  	   this, &GUI::change_orientation );
  connect( miller_l_control->object(), QOverload<int>::of(&QSpinBox::valueChanged),
	   this, &GUI::change_orientation );



  /*
   * gui stuff
   */

  //buttons
  connect( button_render, SIGNAL( clicked() ), sp, SLOT( compute_projection() ) );
  //connect( button_render, SIGNAL( clicked() ), this, SLOT( compute_stats() ) );
  
  connect( button_save, SIGNAL( clicked() ), this, SLOT( save_image_to_file() ) );

  connect( button_measure, SIGNAL( clicked() ), this, SLOT( measure() ) );  
  connect( button_quit, SIGNAL( clicked() ), this, SLOT( quit_app() ) );  

  connect( add_membrane_control, SIGNAL( clicked() ), this, SLOT( add_membrane() ) );
  connect( rm_membrane_control, SIGNAL( clicked() ), this, SLOT( rm_membrane() ) );  

  connect( membranes_control, &QTableWidget::cellChanged, this, &GUI::write_membranes );
  //connect( sp, &sp_qt::parameter_changed, this, &GUI::read_membranes );
  
  // redraw picture when surface_projection class updated the projection
  connect( sp, &sp_qt::projection_changed, this, &GUI::update_view );
  // update the gui if surface_projection updated its status
  //connect( sp, &sp_qt::status_updated , this, &GUI::update_status  );
  // read the updated parameters
  connect( sp, &sp_qt::parameter_changed, this, &GUI::update_gui_from_sp );

  // actiave autoupdate
  connect( autoupdate_control, SIGNAL( stateChanged(int) ), this, SLOT( change_autoupdate(int) ) );

  //invert or not
  connect( invert_control, SIGNAL( stateChanged(int) ), this, SLOT( update_view() ) );  
  
  // get messages from subclass
  connect( sp, &sp_qt::send_message, this, &GUI::output_message );
  connect( sp_stats, &sp_qt::send_message, this, &GUI::output_message );  

  //projection
  connect( this, &GUI::call_compute_projection, sp, &sp_qt::compute_projection );



  // measurements
  connect( this, &GUI::call_update_stats, sp_stats, &sp_qt::update_measurements );
  connect( sp_stats, &sp_qt::measurements_updated, this, &GUI::update_stats );
  connect( sp_stats, &sp_qt::measurements_updated, this, &GUI::update_detailled_stats );  


  // saving
  connect( save_grid_control, &QPushButton::clicked, this, &GUI::save_grid );
  connect( save_surface_points_control, &QPushButton::clicked, this, &GUI::save_surface_points );
  connect( save_topological_network_control, &QPushButton::clicked, this, &GUI::save_network );

  connect( this, &GUI::call_save_grid, sp, &sp_qt::save_grid );
  connect( this, &GUI::call_save_surface, sp, &sp_qt::save_surface_points );
  connect( this, &GUI::call_save_network, sp_stats, &sp_qt::save_topological_network );    
  
  
}



GUI::GUI( QApplication *_app, QWidget *parent ) : QWidget( parent ), app(_app){

  //initialize surface projection
  sp = new sp_qt( progress, status );
  sp->set_n_points_x( 200 );
  sp->set_n_points_y( 200 );
  sp->set_n_points_z( 150 );  
  sp->update_geometry();
  sp->update_containers();

  //initialize surface projection
  
  sp_stats = new sp_qt( progress_stats, status_stats );
  sp_stats->set_n_points_x( 50 );
  sp_stats->set_n_points_y( 50 );
  sp_stats->set_n_points_z( 50 );  
  sp_stats->update_geometry();
  sp_stats->update_containers();


  
  set_up_ui();


  update_gui_from_sp();
  read_membranes();
  
  //set a black background iamge
  uchar *blank_bckgrnd = new uchar[10000]();
  image = new QImage( blank_bckgrnd, 100, 100, 100, QImage::Format_Grayscale8 );
  img_pix = new QPixmap();
  img_pix->convertFromImage( *image );
  

  set_up_signals_and_slots();
  
 
  thread = new QThread();
  sp->moveToThread(thread);
  thread->start();

  t_stats = new QThread();
  sp_stats->moveToThread( t_stats );
  t_stats->start();

  measure();
  emit call_compute_projection();  
}




void GUI::quit_app(){
  thread->quit();
  t_stats->quit();
  app->quit();
}


void GUI::update_view(){
  

  unsigned char* img_data = sp->get_image( invert_control->isChecked() );

  image = new QImage( img_data, sp->get_width(), sp->get_height(), sp->get_width(), QImage::Format_Grayscale8 );

  image->save( "img.jpg" );


  img_pix->convertFromImage( *image );
  draw_area->setPixmap( *img_pix);

  delete( img_data );
  
}


void GUI::paintEvent( QPaintEvent * event ){

   if( img_pix != NULL ){

     //get label dimensions
    int w = draw_area->width();
    int h = draw_area->height();
    
    //set a scaled pixmap to a w x h window keeping its aspect ratio 
    draw_area->setPixmap(img_pix->scaled(w,h,Qt::KeepAspectRatio));
  }

  

}




void GUI::update_status( QString s ){
  output_message( s, 0 );
}





void GUI::update_gui_from_sp(){

  ntucs_control->object()->setValue( sp->get_ntucs() );
  uc_size_control->object()->setValue( sp->get_a() );
  channel_prop_control->object()->setValue( sp->get_channel_prop() );
  surface_type_control->object()->setCurrentIndex( sp->get_type() );


  z_points_control->object()->setValue( sp->get_depth() );
  xy_points_control->object()->setValue( sp->get_width() );

  slice_width_control->object()->setValue( sp->get_slice_width() );
  slice_position_control->object()->setValue( sp->get_slice_height() );  
  
  miller_h_control->object()->setValue( sp->get_h() );
  miller_k_control->object()->setValue( sp->get_k() );
  miller_l_control->object()->setValue( sp->get_l() );  
  
  update_stats( );

  read_membranes();
  
}


void GUI::write_membranes(int row, int col){
  
  //get the data from the table
  int rows = membranes_control->rowCount();
  int cols = membranes_control->columnCount();

  std::vector<double> new_membranes( rows*2, 0 );

  for(unsigned int rr=0; rr<rows; rr++){

    double dist, width;
    try {
      dist = std::stod( membranes_control->item( rr, 0 )->text().toStdString() );
      width = std::stod( membranes_control->item( rr, 1 )->text().toStdString() );
    } catch ( std::invalid_argument e ) {
      output_message( "Invalid argument, resetting to last good configuration" );
      read_membranes();
      return;
    }

    if( rr== 0 ){
      if( dist != 0 || width < 0.02 ){
	dist = 0;
	width = 0.02;
	output_message( "Innermost membrane must be at d=0 and a minimum width of 0.02" );
      }  
    }

    new_membranes[2*rr] = dist;
    new_membranes[2*rr+1] = width;
  }

  sp->change_membranes( new_membranes );
  
}

void GUI::read_membranes(){

  membranes_control->setRowCount( 0 );
  
  std::vector<double> membranes = sp->get_membranes();
  for(unsigned int ii=0; ii<membranes.size(); ii+=2){
    add_membrane( membranes[ii], membranes[ii+1] );
  }
  

}


void GUI::add_membrane( double first, double second ){

  //temporarliry disconnect signals
  disconnect( membranes_control, &QTableWidget::cellChanged, this, &GUI::write_membranes );
  
  int curRow = membranes_control->rowCount();
  membranes_control->insertRow( curRow  );

  membranes_control->setItem( curRow, 0, new QTableWidgetItem() );
  membranes_control->setItem( curRow, 1, new QTableWidgetItem() );  

  membranes_control->item( curRow, 0)->setText( QString::number( first ) );
  membranes_control->item( curRow, 1)->setText( QString::number( second ) );

  //reconnect them
  connect( membranes_control, &QTableWidget::cellChanged, this, &GUI::write_membranes );
}


void GUI::rm_membrane(){
  int ind = membranes_control->currentRow();
  if( ind > 0 ){
    membranes_control->removeRow( ind );
    write_membranes(0, 0);
  }
  else
    output_message( "Innermost membrane can't be removed" );
}

void GUI::save_image_to_file(){
  image->save("image.png");
}


/*
 * types:
 * 0 = status_p
 * 1 = message
 * 2 = stats
 * 3 = status_m
 */
void GUI::output_message( QString msg, int type ){

  std::stringstream text;
  
  if( type == 0 ){
    text.str("");
    text << "Projection" << std::endl << std::endl << msg.toStdString();
    status_bar_status_p->setText( QString( text.str().c_str() )  );
  } else if ( type == 1 ){
    status_bar->showMessage( msg );
  } else if( type == 3 ){
    text.str("");
    text << "Measurement:" << std::endl << std::endl << msg.toStdString();
    status_bar_status_m->setText( text.str().c_str() );
  }

}


void GUI::change_orientation( int val ){

  int h = miller_h_control->object()->value();
  int k = miller_k_control->object()->value();
  int l = miller_l_control->object()->value();  

  sp->change_hkl( h, k, l );
  
}

void GUI::update_stats(){

  std::stringstream pix_size;  
  pix_size << "  XY Pixel size: " << std::setprecision(3) << sp->get_dx() << std::endl
	   << "  Z Pixel size: " << std::setprecision(3) << sp->get_dz()
	   << "          ";
  //pix_size_indicator->setText( QString( pix_size.str().c_str() ) );
  
  status_bar_pixs->setText( QString( pix_size.str().c_str() ) );


  std::stringstream vols, areas, mins;
  vols << "Volumes" << std::endl;
  std::vector<double> vol_data = sp_stats->get_channel_volumes();
  double tmp1 = (vol_data.size()>0) ? vol_data[0] : -1;
  double tmp2 = (vol_data.size()>0) ? vol_data[vol_data.size() - 1] : -1;
  vols << tmp1 << std::endl;
  vols << tmp2;



  areas << "Areas" << std::endl;
  std::vector<double> area_data = sp_stats->get_membrane_surface_area();
  tmp1 = (area_data.size()>0) ? area_data[0] : -1;
  tmp2 = (area_data.size()>0) ? area_data[area_data.size() - 1] : -1;
  areas << tmp1 << std::endl;
  areas << tmp2;  


  // channel diameters
  double inner_min = sp_stats->get_minimal_channel_diameter( 0 );
  double outer_min = sp_stats->get_minimal_channel_diameter( 1 );  
  mins << "Min. channel diameter" << std::endl;
  mins << "Inner membrane: " << inner_min << std::endl;
  mins << "Outer membrane: " << outer_min;


  status_bar_vols->setText( QString( vols.str().c_str() ) );
  status_bar_areas->setText( QString( areas.str().c_str() ) );
  status_bar_mins->setText( QString( mins.str().c_str() ) );
  
}



void GUI::update_detailled_stats(){

  std::stringstream out;


  out << "Volumes:" << std::endl << std::endl;
  std::vector<double> vol_data = sp_stats->get_channel_volumes();
  for( unsigned int ii=0; ii<vol_data.size(); ii+=2){
    out << "Channel " << (ii+2)/2 << ": " << vol_data.at(ii) << std::endl;
  }
  out << std::endl;
  for( unsigned int ii=1; ii<vol_data.size(); ii+=2){
    out << "Membrane " << (ii+1)/2 << ": " << vol_data.at(ii) << std::endl;
  }  
  

  out << std::endl << std::endl << "Areas:" << std::endl << std::endl;
  std::vector<double> area_data = sp_stats->get_membrane_surface_area();
  for( unsigned int ii=0; ii < area_data.size(); ii++){
    out << "Membrane " << ii+1 << ": " << area_data.at(ii) << std::endl;
  }


  // channel diameters
  double inner_min = sp_stats->get_minimal_channel_diameter( 0 );
  double outer_min = sp_stats->get_minimal_channel_diameter( 1 );  
  out << std::endl << std::endl;
  out << "Min. channel diameter" << std::endl << std::endl;
  out << "Inner membrane: " << inner_min << std::endl;
  out << "Outer membrane: " << outer_min;
  

  detailled_stats->setText( out.str().c_str() );
  
}


void GUI::change_autoupdate( int state ){
  
  if( state ){
    connect( sp, &sp_qt::parameter_changed, sp, &sp_qt::compute_projection );
    output_message( "Actiavted autoupdate" );
  } else {
    disconnect( sp, &sp_qt::parameter_changed, sp, &sp_qt::compute_projection );
    output_message( "Deactiavted autoupdate" );    
  }

}


void GUI::measure(){
    sp_stats->copy_parameters( sp );
  sp_stats->set_n_points_x( 76 );
  sp_stats->set_n_points_y( 76 );
  sp_stats->set_n_points_z( 76 );  

  emit call_update_stats();   
}


std::string GUI::get_prefix(){

  std::string tmp = path_prefix_control->text().toStdString();
  if( tmp == "" ){
    tmp += "structure";
  }
  return tmp;

}

void GUI::save_grid(){

  std::string fn = get_prefix();
  fn += "_grid.dat";
  emit call_save_grid( QString( fn.c_str() ) );
  
}

void GUI::save_network(){

  std::string fn = get_prefix();
  fn += "_network_inner.dat";
  emit call_save_network( 0, QString( fn.c_str() ) );

  fn = get_prefix();
  fn += "_network_outer.dat";
  emit call_save_network( 1, QString( fn.c_str() ) );
  
}

void GUI::save_surface_points(){

  int nr_membranes = membranes_control->rowCount();

  std::cout << nr_membranes << std::endl;
  for( unsigned int i=0; i < nr_membranes; i++){
    std::string fn = get_prefix();
    fn+= "_membrane_" + std::to_string( i ) + ".dat";
    //remember that channel counts are started at 1 not 0
    emit call_save_surface( i+1, QString( fn.c_str() ) );
  }

  
}
