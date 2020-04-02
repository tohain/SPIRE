#include "qt_gui.hpp"

void GUI::set_up_ui(){

  // the widgets for the tabs
  controls = new QTabWidget( this );
  controls_basic = new QWidget();
  controls_advanced = new QWidget();
  controls_all = new QWidget();

  controls->addTab( controls_basic, "Basic" );
  controls->addTab( controls_advanced, "Advanced" );
  controls->addTab( controls_all, "All" );  

  // buttons for controls_basic
  button_quit = new QPushButton("Quit", controls_basic);
  button_save = new QPushButton("Save", controls_basic);
  button_render = new QPushButton("Compute", controls_basic);
  button_update_view = new QPushButton("Update", controls_basic);  

  /*
   * structure control
   */ 
  ntucs_control = new QT_labeled_obj<QSpinBox> ( "Unit cells", controls_basic );
  ntucs_control->object()->setRange(0, 50);
  
  uc_size_control = new QT_labeled_obj<QDoubleSpinBox> ( "Unit cell size", controls_basic );
  uc_size_control->object()->setRange(0, 100);
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
  
  slice_position_control = new QT_labeled_obj<QDoubleSpinBox>( "Slice width", controls_basic);    
  slice_position_control->object()->setRange(0, 1);
  slice_position_control->object()->setSingleStep(0.01);  
  
  miller_h_control = new QT_labeled_obj<QSpinBox> ("h", controls_basic );
  miller_h_control->object()->setRange(-100, 100);
  
  miller_k_control = new QT_labeled_obj<QSpinBox> ("k", controls_basic );
  miller_k_control->object()->setRange(-100, 100);
  
  miller_l_control = new QT_labeled_obj<QSpinBox> ("l", controls_basic );  
  miller_l_control->object()->setRange(-100, 100);


  /*
   * Layout
   */
  
  //set up spacers
  v_spacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
  h_spacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

  
  //the main layout of the form
  main_layout = new QHBoxLayout( this );
  main_layout->addWidget( controls);
  main_layout->addWidget( draw_area );
  
  // the layout of the buttons
  buttons_layout = new QHBoxLayout();
  // the layout of the tabs
  controls_basic_layout = new QVBoxLayout( controls_basic );
  controls_advanced_layout = new QVBoxLayout( controls_advanced );
  controls_all_layout = new QVBoxLayout( controls_all );  

  structure_settings = new QHBoxLayout();
  resolution_settings = new QHBoxLayout();
  slice_settings = new QHBoxLayout();
  
  structure_settings->addLayout( ntucs_control->layout() );
  structure_settings->addLayout( uc_size_control->layout() );
  structure_settings->addLayout( channel_prop_control->layout() );
  structure_settings->addLayout( surface_type_control->layout() );  

  resolution_settings->addLayout( xy_points_control->layout() );
  resolution_settings->addLayout( z_points_control->layout() );
  resolution_settings->addWidget( pix_size_indicator );
  resolution_settings->addWidget( invert_control );

  slice_settings->addLayout( slice_width_control->layout() );
  slice_settings->addLayout( slice_position_control->layout() );
  slice_settings->addItem( h_spacer );
  slice_settings->addLayout( miller_h_control->layout() );
  slice_settings->addLayout( miller_k_control->layout() );
  slice_settings->addLayout( miller_l_control->layout() );    
  
  controls_basic_layout->addLayout( structure_settings );
  controls_basic_layout->addItem( v_spacer );
  controls_basic_layout->addLayout( slice_settings );
  controls_basic_layout->addItem( v_spacer );
  controls_basic_layout->addLayout( resolution_settings );
  controls_basic_layout->addItem( v_spacer );
  controls_basic_layout->addLayout( buttons_layout );
    
  buttons_layout->addWidget( button_render );
  buttons_layout->addWidget( button_update_view );
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
	   sp, &sp_qt::change_h );
  connect( miller_k_control->object(), QOverload<int>::of(&QSpinBox::valueChanged),
	   sp, &sp_qt::change_k );
  connect( miller_l_control->object(), QOverload<int>::of(&QSpinBox::valueChanged),
	   sp, &sp_qt::change_l );  



  /*
   * gui stuff
   */

  //buttons
  connect( button_render, SIGNAL( clicked() ), sp, SLOT( compute_projection() ) );
  connect( button_save, SIGNAL( clicked() ), this, SLOT( save_image_to_file() ) );
  connect( button_quit, SIGNAL( clicked() ), this, SLOT( quit_app() ) );
  connect( button_update_view, SIGNAL( clicked() ), this, SLOT( update_view() ) );

  // redraw picture when surface_projection class updated the projection
  connect( sp, &sp_qt::projection_changed, this, &GUI::update_view );
  // update the gui if surface_projection updated its status
  connect( sp, &sp_qt::status_updated , this, &GUI::update_status  );
  // read the updated parameters
  connect( sp, &sp_qt::parameter_changed, this, &GUI::update_gui_from_sp );

}



GUI::GUI( QApplication *_app, QWidget *parent ) : QWidget( parent ), app(_app){

  //initialize surface projection
  sp = new sp_qt( progress, status );
  sp->update_geometry();
  sp->update_containers();
  
  set_up_ui();


  update_gui_from_sp();
  
  //set a black background iamge
  uchar *blank_bckgrnd = new uchar[10000]();
  image = new QImage( blank_bckgrnd, 100, 100, 100, QImage::Format_Grayscale8 );
  img_pix = new QPixmap();
  img_pix->convertFromImage( *image );
  

 
  thread = new QThread();
  sp->moveToThread(thread);
  thread->start();
  
  set_up_signals_and_slots();
  
}






void GUI::quit_app(){
  app->quit();
}


void GUI::update_view(){
  

  unsigned char* img_data = sp->get_image( invert_control->isChecked() );

  image = new QImage( img_data, sp->get_width(), sp->get_height(), sp->get_width(), QImage::Format_Grayscale8 );

  image->save( "img.jpg" );


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
  std::cout << s.toStdString() << std::endl;
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
  
  std::stringstream pix_size;
  pix_size << "XY Pixel size: " << std::setprecision(3) << sp->get_dx() << std::endl
	   << "Z Pixel size: " << std::setprecision(3) << sp->get_dz();
  pix_size_indicator->setText( QString( pix_size.str().c_str() ) );

}



void GUI::save_image_to_file(){


  image->save("image.png");

}
