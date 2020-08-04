#include "qt_gui.hpp"

void GUI::set_up_ui(){

  // the widgets for the tabs
  controls = new QTabWidget( this );
  controls_basic = new QWidget();
  controls_save = new QWidget();
  controls_measurement = new QWidget();
  manual_widget = new QWidget();

  controls->addTab( controls_basic, "Parameters" );
  controls->addTab( controls_measurement, "Measurements" );
  controls->addTab( controls_save, "Export" );
  controls->addTab( manual_widget, "Manual" );


  //manual tab
  manual = new QPlainTextEdit( manual_widget );
  manual->setPlainText( QString( ttips.manual.c_str() ) );
  
  // stats tab
  detailled_stats = new QLabel( controls_save );

  // save tab
  choose_path_prefix = new QPushButton( "Choose location and prefix", controls_save );
  save_grid_control = new QPushButton( "Save grid", controls_save );
  save_surface_points_control = new QPushButton( "Save membrane points", controls_save );
  save_topological_network_control = new QPushButton( "Save network", controls_save );
  path_prefix_control = new QT_v_labeled_obj<QLineEdit>( "Path and prefix", controls_save );
  
  // buttons for controls_basic
  button_quit = new QPushButton("Quit", controls_basic);
  button_save = new QPushButton("Save Image", controls_basic);
  button_render = new QPushButton("Compute Projection", controls_basic);
  button_measure = new QPushButton("Measure", controls_basic);  
  
  /*
   * structure control
   */ 
  ntucs_control = new QT_v_labeled_obj<QSpinBox> ( "Unit cells", controls_basic );
  ntucs_control->object()->setRange(1, 50);
  
  uc_size_control_a = new QT_v_labeled_obj<QDoubleSpinBox> ( "Unit cell size (a)", controls_basic );
  uc_size_control_a->object()->setRange(0.01, 999);
  uc_size_control_a->object()->setSingleStep(0.01);
  uc_size_control_c = new QT_v_labeled_obj<QDoubleSpinBox> ( "Unit cell size (c)", controls_basic );
  uc_size_control_c->object()->setRange(0.01, 999);
  uc_size_control_c->object()->setSingleStep(0.01); 
  
  channel_prop_control = new QT_v_labeled_obj<QDoubleSpinBox> ( "Volume proportions", controls_basic );
  //channel_prop_control->object()->setRange(0, 1);
  channel_prop_control->object()->setSingleStep(0.01);  
  
  surface_type_control = new QT_v_labeled_obj<QComboBox> ( "Surface type", controls_basic );
  std::vector<std::string> sfc_types = sp->get_surface_choices();
  for(unsigned int ii=0; ii<sfc_types.size(); ii++){
    surface_type_control->object()->insertItem( ii, QString( sfc_types.at(ii).c_str() ) );
  }

  /*
   * Resolution control
   */
  xy_points_control = new QT_v_labeled_obj<QSpinBox> ( "XY resolution", controls_basic );
  xy_points_control->object()->setRange(1, 3000);
 
  z_points_control = new QT_v_labeled_obj<QSpinBox> ( "Z resolution", controls_basic );  
  z_points_control->object()->setRange(1, 3000);

  pix_size_indicator = new QLabel("");
  pix_size_indicator->setWordWrap( true );
  pix_size_indicator->setAlignment(Qt::AlignCenter);
  pix_size_indicator->setTextFormat(Qt::PlainText);

  invert_control = new QCheckBox( controls_basic );
  invert_control->setText( "Invert colors" );

  autoupdate_control = new QCheckBox( controls_basic );
  autoupdate_control->setText( "Autoupdate" );  

  image_scaling_control = new QT_v_labeled_obj<QComboBox>( "Scaling", controls_basic );
  std::vector<std::string> imgs_types = sp->get_img_scaling_choices();
  for(unsigned int ii=0; ii<imgs_types.size(); ii++){
    image_scaling_control->object()->insertItem( ii, QString( imgs_types.at(ii).c_str() ) );
  }
  
  //Set up the draw area
  draw_area = new QLabel();
  draw_area->setAlignment(Qt::AlignVCenter | Qt::AlignHCenter);
  draw_area->setMinimumSize(500,500);
  draw_area->setSizePolicy(QSizePolicy::MinimumExpanding,
                     QSizePolicy::MinimumExpanding);



  /*
   * Slice control 
   */

  slice_width_control = new QT_h_labeled_obj<QDoubleSpinBox>( "Slice width", controls_basic);
  slice_width_control->object()->setRange(0, 250);
  slice_width_control->object()->setSingleStep(0.01);

  slice_length_control = new QT_h_labeled_obj<QDoubleSpinBox>( "Slice length", controls_basic);
  slice_length_control->object()->setRange(0, 250);
  slice_length_control->object()->setSingleStep(0.01);

  slice_height_control = new QT_h_labeled_obj<QDoubleSpinBox>( "Slice height", controls_basic);
  slice_height_control->object()->setRange(0, 250);
  slice_height_control->object()->setSingleStep(0.01);  
  
  slice_position_control = new QT_h_labeled_obj<QDoubleSpinBox>( "Slice position", controls_basic);    
  slice_position_control->object()->setRange(0, 1);
  slice_position_control->object()->setSingleStep(0.01);  
  
  miller_h_control = new QT_h_labeled_obj<QSpinBox> ("h", controls_basic );
  miller_h_control->object()->setRange(-100, 100);
  
  miller_k_control = new QT_h_labeled_obj<QSpinBox> ("k", controls_basic );
  miller_k_control->object()->setRange(-100, 100);
  
  miller_l_control = new QT_h_labeled_obj<QSpinBox> ("l", controls_basic );  
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
  manual_widget_layout = new QVBoxLayout( manual_widget );


  controls_measurement_layout->addWidget( detailled_stats );
  
  controls_save_layout->addLayout( path_prefix_control->layout() );
  controls_save_layout->addWidget( choose_path_prefix );
  controls_save_layout->addItem( v_spacer );  
  controls_save_layout->addWidget( save_grid_control );
  controls_save_layout->addWidget( save_surface_points_control );
  controls_save_layout->addWidget( save_topological_network_control );
    

  manual_widget_layout->addWidget( manual );
  
  structure_settings = new QHBoxLayout();
  resolution_settings = new QHBoxLayout();
  slice_orientation_layout = new QVBoxLayout();
  slice_dimension_layout = new QVBoxLayout();
  slice_settings = new QHBoxLayout();  
  membrane_settings = new QHBoxLayout();
  membrane_buttons_layout = new QVBoxLayout();
  
  structure_settings->addLayout( ntucs_control->layout() );
  structure_settings->addLayout( uc_size_control_a->layout() );
  structure_settings->addLayout( uc_size_control_c->layout() );  
  structure_settings->addLayout( channel_prop_control->layout() );
  structure_settings->addLayout( surface_type_control->layout() );  

  resolution_settings->addLayout( xy_points_control->layout() );
  resolution_settings->addLayout( z_points_control->layout() );
  resolution_settings->addLayout( image_scaling_control->layout() );
  resolution_settings->addWidget( invert_control );
  resolution_settings->addWidget( autoupdate_control );




  
  slice_dimension_layout->addLayout( slice_width_control->layout() );
  slice_dimension_layout->addLayout( slice_height_control->layout() );
  slice_dimension_layout->addLayout( slice_length_control->layout() );  
  slice_dimension_layout->addLayout( slice_position_control->layout() );

  slice_settings->addLayout( slice_dimension_layout );


  slice_settings->addItem( h_spacer );
  
  slice_orientation_layout->addLayout( miller_h_control->layout() );
  slice_orientation_layout->addLayout( miller_k_control->layout() );
  slice_orientation_layout->addLayout( miller_l_control->layout() );    

  slice_settings->addLayout( slice_orientation_layout );
  
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



void GUI::set_up_tooltips(){

   ntucs_control->object()->setToolTip( QString( ttips.ntucs_tooltip.c_str() ) );
   uc_size_control_a->object()->setToolTip( QString( ttips.aa_tooltip.c_str() ) );
   uc_size_control_c->object()->setToolTip( QString( ttips.ac_tooltip.c_str() ) );
   channel_prop_control->object()->setToolTip( QString( ttips.level_set_tooltip.c_str() ) );
   surface_type_control->object()->setToolTip( QString( ttips.type_tooltip.c_str() ) );

   slice_width_control->object()->setToolTip( QString( ttips.slicewidth_tooltip.c_str() ) );
   //slice_length_control->object()->setToolTip( QString( ttips.slicewidth_tooltip.c_str() ) );
   //slice_height_control->object()->setToolTip( QString( ttips.slicewidth_tooltip.c_str() ) );   
   slice_position_control->object()->setToolTip( QString( ttips.sliceheight_tooltip.c_str() ) );

   miller_h_control->object()->setToolTip( QString( ttips.hkl_tooltip.c_str() ) );
   miller_k_control->object()->setToolTip( QString( ttips.hkl_tooltip.c_str() ) );
   miller_l_control->object()->setToolTip( QString( ttips.hkl_tooltip.c_str() ) );

   xy_points_control->object()->setToolTip( QString( ttips.nr_points_xy_tooltip.c_str() ) );
   z_points_control->object()->setToolTip( QString( ttips.nr_points_z_tooltip.c_str() ) );

   invert_control->setToolTip( QString( ttips.invert_tooltip.c_str() ) );
   autoupdate_control->setToolTip( QString( ttips.autoupdate_tooltip.c_str() ) );
   image_scaling_control->object()->setToolTip( QString( ttips.image_scaling_tooltip.c_str() ) );

   membranes_control->setToolTip ( QString( ttips.membranes_settings_tooltip.c_str() ) );
   add_membrane_control->setToolTip ( QString( ttips.membranes_add_tooltip.c_str() ) );
   rm_membrane_control->setToolTip ( QString( ttips.membranes_remove_tooltip.c_str() ) );   


   button_quit->setToolTip( QString( ttips.button_quit.c_str() ) );
   button_measure->setToolTip( QString( ttips.button_measure.c_str() ) );
   button_render->setToolTip( QString( ttips.button_render.c_str() ) );
   button_save->setToolTip( QString( ttips.button_save.c_str() ) );   

   save_grid_control->setToolTip( QString( ttips.button_save_grid.c_str() ) );
   save_surface_points_control->setToolTip( QString( ttips.button_save_membranes.c_str() ) );
   save_topological_network_control->setToolTip( QString( ttips.button_save_network.c_str() ) );   

   path_prefix_control->object()->setToolTip( QString( ttips.save_prefix_tooltip.c_str() ) );
   
   choose_path_prefix->setToolTip( QString( ttips.choose_prefix_path.c_str() ) );
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
  connect( uc_size_control_a->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_uc_scale_ab );
  connect( uc_size_control_c->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_uc_scale_c );


  // channel proportion
  connect( channel_prop_control->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_vol_prop );

  // slice width
  connect( slice_width_control->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_slice_width );
  connect( slice_height_control->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_slice_height );
  connect( slice_length_control->object(), QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	   sp, &sp_qt::change_slice_length );  
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

  // scaling
  connect(image_scaling_control->object(), QOverload<int>::of(&QComboBox::currentIndexChanged),
	  this, &GUI::update_view);  
  
  // get messages from subclass
  connect( sp, &sp_qt::send_message, this, &GUI::output_message );
  connect( sp, &sp_qt::set_status, this, &GUI::set_state );
  connect( sp_stats, &sp_qt::send_message, this, &GUI::output_message );  
  connect( sp_stats, &sp_qt::set_status, this, &GUI::set_state );
  
  //projection
  connect( this, &GUI::call_compute_projection, sp, &sp_qt::compute_projection );



  // measurements
  connect( this, &GUI::call_update_stats, sp_stats, &sp_qt::update_measurements );
  connect( sp_stats, &sp_qt::measurements_updated, this, &GUI::update_stats );
  connect( sp_stats, &sp_qt::measurements_updated, this, &GUI::update_detailled_stats );  

  connect( this, &GUI::call_set_measurement_status, this, &GUI::set_measurements_status );

  // saving
  connect( choose_path_prefix, &QPushButton::clicked, this, &GUI::choose_export_prefix );
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
  set_up_tooltips();

  update_gui_from_sp();
  read_membranes();
  
  //set a black background iamge
  img_data = new uchar[10000]();
  image = new QImage( img_data, 100, 100, 100, QImage::Format_Grayscale8 );
  img_pix = new QPixmap();
  img_pix->convertFromImage( *image );
  


  thread = new QThread();
  sp->moveToThread(thread);
  thread->start();

  t_stats = new QThread();
  sp_stats->moveToThread( t_stats );
  t_stats->start();

  
  set_up_signals_and_slots();
  
  //measure();
  emit call_compute_projection();  
}




void GUI::quit_app(){

  auto reply = QMessageBox::question( this, "Really quit?", "Are you sure you want to quit?", QMessageBox::Yes|QMessageBox::No );
  
  if( reply == QMessageBox::Yes ){
    thread->quit();
    t_stats->quit();
    app->quit();
  } else {
    // do nothing, just close dialog box
  }
  
}


void GUI::update_view(){

  // delete old img_data if existant.
  if( img_data != NULL ){
    delete( img_data );
  }
  img_data = sp->get_image( invert_control->isChecked(), image_scaling_control->object()->currentText().toStdString() );

  // this does *NOT* seem to copy the img_data into its own object, so
  // keep that img_data array around!
  image = new QImage( img_data, sp->get_width(), sp->get_height(), sp->get_width(), QImage::Format_Grayscale8 );

  img_pix->convertFromImage( *image );
  draw_area->setPixmap( *img_pix);  
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

  uc_size_control_a->object()->setValue( sp->get_uc_scale_ab() );
  uc_size_control_c->object()->setValue( sp->get_uc_scale_c() );  

  channel_prop_control->object()->setValue( sp->get_channel_prop() );
  surface_type_control->object()->setCurrentIndex( sp->get_type() );

  z_points_control->object()->setValue( sp->get_depth() );
  xy_points_control->object()->setValue( sp->get_width() );

  slice_width_control->object()->setValue( sp->get_slice_width() );
  slice_position_control->object()->setValue( sp->get_slice_position() );  
  
  miller_h_control->object()->setValue( sp->get_h() );
  miller_k_control->object()->setValue( sp->get_k() );
  miller_l_control->object()->setValue( sp->get_l() );  

  sp->update_periodicity_length();
  
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

  QString filename = QFileDialog::getSaveFileName( this, "Select File", "./", tr("*.png"), NULL, QFileDialog::DontUseNativeDialog );
  
  if( !filename.endsWith( ".png", Qt::CaseInsensitive ) ){
    filename.append(".png");
  } 
  
  image->save( filename );
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

  pix_size << " box size X/Y:  "  << sp->get_L()[0] << "/" << sp->get_L()[1]  << std::endl;
  pix_size << " box heigth :  ";
  pix_size << sp->get_slice_width() << "    " << std::endl;
  pix_size << "periodicity length : ";
  pix_size << sp->get_periodicity_length() << std::endl;

  
  
  status_bar_pixs->setText( QString( pix_size.str().c_str() ) );


  
  std::stringstream vols, areas, mins;


  /*
  vols << "Volumes" << std::endl;
  std::vector<double> vol_data = sp_stats->get_channel_volumes();
  double tmp1 = (vol_data.size()>0) ? vol_data[0] : -1;
  double tmp2 = (vol_data.size()>0) ? vol_data[vol_data.size() - 1] : -1;
  vols << tmp1 << std::endl;
  vols << tmp2;
  */

  /*

  areas << "Areas" << std::endl;
  std::vector<double> area_data = sp_stats->get_membrane_surface_area();
  tmp1 = (area_data.size()>0) ? area_data[0] : -1;
  tmp2 = (area_data.size()>0) ? area_data[area_data.size() - 1] : -1;
  areas << tmp1 << std::endl;
  areas << tmp2;  
  */

  areas << "Pixelinfo" << std::endl;
  areas << "x: " << sp->get_width() << " (" << sp->get_dx() << ")" << std::endl;
  areas << "y: " << sp->get_height() << " (" << sp->get_dy() << ")" << std::endl;
  areas << "z: " << sp->get_depth() << " (" << sp->get_dz() << ")" << std::endl;  


  vols << "uc info" << std::endl;
  vols << "x: " << sp->get_a()[0] << std::endl;// << "(" << sp->inv_a[0] << ")" << std::endl;
  vols << "y: " << sp->get_a()[1] << std::endl;//"(" << sp->inv_a[1] << ")" << std::endl;
  vols << "z: " << sp->get_a()[2] << std::endl;//"(" << sp->inv_a[2] << ")" << std::endl;  
  
  /*
  
  // channel diameters
  double inner_min = sp_stats->get_minimal_channel_diameter( 0 );
  double outer_min = sp_stats->get_minimal_channel_diameter( 1 );  
  mins << "Min. channel diameter" << std::endl;
  mins << "Inner membrane: " << inner_min << std::endl;
  mins << "Outer membrane: " << outer_min;

  */

  mins << "orientation" << std::endl;
  mins << "phi=" << sp->get_phi() << std::endl;
  mins << " theta=" << sp->get_theta() << std::endl;

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


void GUI::choose_export_prefix(){
  QString prefix = QFileDialog::getSaveFileName( this, "Select Prefix", "./", tr("*"), NULL, QFileDialog::DontUseNativeDialog );

  path_prefix_control->object()->setText( prefix );
}


std::string GUI::get_prefix(){
  
  std::string tmp = path_prefix_control->object()->text().toStdString();
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

  for( unsigned int i=0; i < nr_membranes; i++){
    std::string fn = get_prefix();
    fn+= "_membrane_" + std::to_string( i ) + ".dat";
    //remember that channel counts are started at 1 not 0
    emit call_save_surface( i+1, QString( fn.c_str() ) );
  }

  
}


void GUI::update_unitcell_size(){
  emit call_change_uc_size( uc_size_control_a->object()->value(), uc_size_control_c->object()->value() );
}

/**
 * switches the state indicator of the GUI
 *
 *  \param[in] what 0: projection, 1: measurement
 * \param[in] state 0: ready, 1: busy
 */ 
void GUI::set_state( int what, int state ){

  std::stringstream text;
  text.str("");
  QLabel *display;

  if( what == 0 ){
    text << "Projection" << std::endl << std::endl;
    display = status_bar_status_p;
  } else if ( what == 1 ){
    text << "Measurement:" << std::endl << std::endl;
    display = status_bar_status_m;    
  } 


  std::string style;
  
  if( state == 0 ){    
    text << "Ready";
    style = "QLabel { background-color : green; color : black; }";
  } else if ( state == 1){
    text << "Busy";
    style = "QLabel { background-color : red; color : black; }";
  }
  display->setStyleSheet( QString( style.c_str() ) );
  display->setText( text.str().c_str() );

}


/**
 * Changes the background color of the measurement labels, depending
 * if they are up to date
 *
 * \param[in] state 0: up to date, 1: outdated
 */
void GUI::set_measurements_status( int state ){

  std::string style;
  
  if( state == 0 ){
    style = "QLabel { background-color : red; color : black; }";
  } else if ( state == 1 ){
    style = "QLabel { background-color : green; color : black; }";
  }
  
  status_bar_pixs->setStyleSheet( QString( style.c_str() ) );
  status_bar_vols->setStyleSheet( QString( style.c_str() ) );
  status_bar_areas->setStyleSheet( QString( style.c_str() ) );
  status_bar_mins->setStyleSheet( QString( style.c_str() ) );

}



void GUI::do_something(){  
  std::cout << "done something\n"; 
  
}
