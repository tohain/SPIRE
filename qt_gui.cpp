#include "qt_gui.hpp"


GUI::GUI( QApplication *_app, QWidget *parent ) : QWidget( parent ), app(_app){

  //set a black background iamge
  uchar *blank_bckgrnd = new uchar[10000]();
  image = new QImage( blank_bckgrnd, 100, 100, 100, QImage::Format_Grayscale8 );
  img_pix = new QPixmap();
  img_pix->convertFromImage( *image );
  
  //initialize surface projection
  sp = new worker( progress, status );
  sp->set_n_points_x( 175 );
  sp->set_n_points_y( 175 );
  sp->set_n_points_z( 150 );  
  sp->update_geometry();
  sp->update_containers();

  button_quit = new QPushButton("Quit", this);
  button_save = new QPushButton("Compute", this);
  button_render = new QPushButton("Show", this);

  draw_area = new QLabel();
  draw_area->setAlignment(Qt::AlignVCenter | Qt::AlignHCenter);
  draw_area->setMinimumSize(100,100);
  
  main_layout = new QHBoxLayout( this );
  buttons_layout = new QVBoxLayout( this );
  
  buttons_layout->addWidget( button_quit );
  buttons_layout->addWidget( button_save );
  buttons_layout->addWidget( button_render );

  main_layout->addLayout( buttons_layout);
  main_layout->addWidget( draw_area );


  thread = new QThread();
  sp->moveToThread(thread);

  connect(sp, SIGNAL(workRequested()), thread, SLOT(start()));
  connect(thread, SIGNAL(started()), sp, SLOT(doWork()));
  connect(sp, SIGNAL(finished()), thread, SLOT(quit()));

  //connect signals and slots

  //connect( button_quit, SIGNAL (clicked(int) ), this, SLOT (printme(int)));
  connect( button_render, SIGNAL( clicked() ), this, SLOT( pushed_button() ) );
  connect( button_save, SIGNAL( clicked() ), this, SLOT( pushed_button_2() ) );
  connect( button_quit, SIGNAL( clicked() ), this, SLOT( pushed_button_3() ) );    
  
  
}


void worker::doWork(){

  compute_projection();
  std::cout << "done!" << std::endl;
  emit finished();
}



void GUI::pushed_button_3(){
  app->quit();

}

void GUI::pushed_button_2(){

  sp->requestWork();

}

void GUI::pushed_button(){
  

  unsigned char* img_data = sp->get_image();

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



void worker::requestWork(){
  emit workRequested();
}

