/* Projection tool - compute planar projection of triply periodic
 * minimal surfaces 
 * Copyright (C) 2021 Tobias Hain
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

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <sstream>

#include <QWidget>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QPainter>
#include <QImage>
#include <QLabel>
#include <QPixmap>
#include <QApplication>
#include <QTableWidget>
#include <QHeaderView>
#include <QStatusBar>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QSpacerItem>
#include <QCheckBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QScrollArea>
#include <QTextStream>
#include <QFile>
#include <QProgressBar>
#include <QRegularExpression>
#include <QtSvg/QSvgWidget>

#include <QThread>
#include <QMutex>

#include "qt_sp.hpp"
#include "sp_gui_tooltips.h"
#include "slice_orientation_visualisation.hpp"
#include "vec_mat_math.hpp"

#include "global_settings.hpp"
#include "batch_lib.hpp"


/**
 * a short qt wrapper around the batch_lib to it has signals and slots
 */
class qt_bc : public QObject, public batch_creation {

  Q_OBJECT;

public:
  
  qt_bc ( surface_projection &sp_ ) : batch_creation( sp_ ) {
  }

  void set_callback( sp_callback *callback ){
    cb = callback;
  }

public slots:
  void start_loop(){
    this->do_loop( *cb );    
    emit finished_batch_loop();
  }
  
signals:
  void finished_batch_loop();
  
private:
  sp_callback *cb;

};

/** 
 * The callback functor for batch creation
 */
class gui_callback : public QObject, public sp_callback {

  Q_OBJECT;
  
public:
  gui_callback( surface_projection &sp,
		batch_creation &bc_,
		std::string prefix ) :
    sp_callback( sp ),
    fn_prefix( prefix ),
    bc( bc_ ),
    internal_counter (0),
    keep_going ( true )
  {

  }

  ~gui_callback(){
  }

  void set_prefix( std::string new_pre ){
    fn_prefix = new_pre;   
  }

  void init(){
    summary.open( fn_prefix + "_summary.txt" );
    internal_counter=0;
    external_counter=0;
  }

  void finalize(){
    summary.close();
  }



  /**
   * this function computes the projection and adds the relevant
   * parameters to the png
   *
   * the parameters are added in a margin below the image. The size of
   * the text is tried to kept constant at a height of
   * pixsize_fractional of the image height, but minium of
   * 25pixels. If the image resolution is too low, the image is
   * upscaled so it can accomodate the text in the desired resolution.
   *
   * The text is arranged on a grid, with a cell length of the longest
   * word size. If a shorter word fits the line, it is squeezed in,
   * however, only at the start of a grid cell
   */
  bool operator()( std::vector< std::vector<double>::iterator > pars ){


    /*
     * compute the projection and get the image
     */

    std::string fn = fn_prefix + "_" +  std::to_string( internal_counter ) + ".png";       
    sp.update_geometry();
    sp.compute_projection();
    sp.write_png( fn + ".png" );
    
    unsigned char* img = sp.get_image( invert, scaling );
    int img_x = sp.get_width(); int img_y = sp.get_height();


    /*
     * creating the legend
     */ 

    // constants
    double pixsize_fractional = 0.045; // size of the font in pixels
    double margin = 1.2;
    unsigned int min_pixsize = 25;

        
    //
    // would be pixel size if image is not scaled
    //
    unsigned int pixsize = (pixsize_fractional * img_y);

    // if pix size gets too small, it is not rendered readable, so we
    // need to upscale the image
    double scale = 1.0;
    if( pixsize < min_pixsize ){
      scale = (float) min_pixsize / pixsize;
    }

    // setting the pixel size; this makes the point size not used
    QFont font = QFont();
    font.setCapitalization( QFont::AllLowercase ); // not sure if we need this;
      
    // set the current font size, so the longest word length can be
    // measured in the non-scaled image
    font.setPixelSize( pixsize );
    
    //
    // get longest string in pixels to compute grid sizing (cell length)
    //
    QFontMetrics fm ( font );
    
    // store the text
    std::vector<std::string> parameter_words;

    // the length of the longest word in unscaled pixels
    unsigned int max_width = 0;
    for( unsigned int ii=0; ii<bc.ops.parameters.size(); ii++){

      std::stringstream ss;

      // get the short parameter name
      size_t ind = std::distance( surface_projection::parameter_names.begin(),
				  std::find( surface_projection::parameter_names.begin(),
					     surface_projection::parameter_names.end(),
					     bc.ops.parameters[ii] ) );


      // construct the string containing name and parameter
      ss << surface_projection::parameter_names_short[ind] << "=";
      auto types = sp.get_surface_choices();

      // special treatment for surface type so string is printed instead of index
      if( bc.ops.parameters[ind] == surface_projection::parameter_names[0] ){
	ss << types[ int( sp.get_parameter( surface_projection::parameter_names[0] )  )  ] << std::endl;
      } else {      
	//ss << *(pars[ii]) << std::endl;
	std::cout << surface_projection::parameter_names[ind]
		  << sp.get_parameter( surface_projection::parameter_names[ind] ) << std::endl;
	ss << sp.get_parameter( surface_projection::parameter_names[ind] ) << std::endl;
      }

      // store words
      parameter_words.push_back( ss.str() );

      // compute length
      unsigned int word_width = fm.horizontalAdvance( ss.str().c_str() );      
      if( word_width > max_width ){
	max_width = word_width;
      }      
    }

    
    //
    // compute the number of columns and rows we need, not really
    // needed for the grid but we need to estimate how many pixels we
    // need to add for the margin text
    //
    unsigned int nr_cols = int(img_x / max_width);
    if( nr_cols == 0 ) nr_cols++;
    unsigned int nr_rows = int(pars.size() / nr_cols)+1;

    // now add a few rows as writing area and realloc the memory
    unsigned int new_img_y = img_y + int( margin * nr_rows * pixsize);    
    img = (unsigned char*) realloc (img, sizeof(unsigned char) * img_x * new_img_y);

    // make the background white/black
    memset( img + sizeof(unsigned char) * img_x * img_y,
	    invert ? 0 : 255,
	    sizeof(unsigned char)*img_x*(new_img_y-img_y) );
    
    // create the qimage object where we can write on
    QImage qimg = QImage( img, img_x, new_img_y, img_x, QImage::Format_Grayscale8 );
    // scale it, if image is large enough hopefully this does nothing
    qimg = qimg.scaledToHeight( int(scale * new_img_y) );

    // set real font size after image scaling
    font.setPixelSize( int(scale * pixsize) );    
    
    // font color
    QPen pen = QPen();
    pen.setColor( invert ? "#ffffff" : "#000000" );

    // setup painter
    QPainter poet ( &qimg );
    poet.setFont( font );
    poet.setPen( pen );    

    // update the new font size
    fm = QFontMetrics( font );

    // put the words below the text. Only start at grid locations, but
    // if a cell is cutoff by the image, check if we can squeeze a
    // shorter word in
    int ind=0, x=0, y = scale * (img_y + pixsize - 0.5*(1-margin)*pixsize);
    int word_width_pix;
    int c=0; //current col
    while( ind < parameter_words.size() ){

      // get word length in scale pixels
      word_width_pix = fm.horizontalAdvance( parameter_words[ind].c_str() );

      // check if word matches this line
      if( x + word_width_pix > qimg.width() ){
	// line does not fit that word, advance to next one
	x = 0;
	y += scale * margin * pixsize;
	c=0;
      }
      
      // draw it
      poet.drawText( x, y, parameter_words[ind].c_str() );
      // increase the col
      c++;
      // set x to the next grid start
      x = c * scale * max_width;
            
      ind++;
    }
    
    // write the image
    qimg.save( QString( fn.c_str() ) );

    // free memory
    free( img );
    
    summary << fn << "    ";

    for( auto it : surface_projection::parameter_names ){
      summary << it << "=" << sp.get_parameter( it ) << " ";
    }

    summary << std::endl;
    internal_counter++;
    emit updated_batch_progress( external_counter );

    return keep_going;
  }

  // keeps track of how many projections has been computed, used for
  // filenames
  unsigned int internal_counter;
  
  batch_creation &bc;
  std::string fn_prefix;
  std::ofstream summary;

  bool keep_going;
		       
public slots:
  void update_status ( bool go_on ){
    keep_going = go_on;
  }
  
signals:
  void updated_batch_progress( int p );
  
};



/**
 * short wrapper for an object including a label combined with an
 * object
 */
template <class QT_O>
class QT_labeled_obj {

public:
  /**
   * Constructor orientation is a char array of length 2. First
   * character: v/h horizontal or vertical laytou, second character
   * l/o label or object first
   *
   * \param[in] orientation Orientation and Layout
   * \param[in] label_str Text of the label
   * param[in] parent The QT parent
   */
  QT_labeled_obj( std::string orientation, std::string label_str, QWidget *parent = NULL ){
    obj = new QT_O ( parent );
    lbl = new QLabel( QString(label_str.c_str()), parent );

    if( orientation[0] == 'v' ){
      lyt = new QVBoxLayout();
    } else {
      lyt = new QHBoxLayout();
    }

    if( orientation[1] == 'l' ){      
      lyt->addWidget( lbl );
      lyt->addWidget( obj );
    } else {
      lyt->addWidget( obj );
      lyt->addWidget( lbl );
    }
  }
  
  ~QT_labeled_obj(){
    delete( obj );
    delete( lbl );
    delete( lyt );
  }

  QBoxLayout* layout(){
    return lyt;
  }
  QT_O* object(){
    return obj;
  }
  QLabel* label(){
    return lbl;
  }  
  
private:
  QT_O *obj;
  QLabel *lbl;
  QBoxLayout *lyt;
};


class GUI : public QWidget {

  Q_OBJECT

public:
  explicit GUI( QApplication *app, QLocale *def_locale_, QWidget *parent = 0 );


  ~GUI();
    
protected:
  
  void paintEvent( QPaintEvent * event );
  
private:

  int sp_state, sp_stat_state;

  const unsigned int space_between_items = 20;
  
  // locale
  QLocale *def_locale;
  
  //tooltips
  tooltips ttips;

  // functions
  void set_up_ui();
  void set_up_signals_and_slots();
  void set_up_tooltips();
  
  // ui elements
  QLabel *draw_area;
  QSvgWidget *orientation_visualisation;

  QTabWidget *controls;
  QWidget *parameters_widget;
  QWidget *measurement_widget;

  QWidget *batch_widget;
  QScrollArea *batch_scroll_area;
  QWidget *batch_scroll_subwidget;
  
  QWidget *save_widget;
  QWidget *about_widget;
  QWidget *license_widget;  

  // manual tab
  QPlainTextEdit *manual;

  // measurement tab
  QT_labeled_obj<QTableWidget> *measurements_slice;
  QT_labeled_obj<QTableWidget> *measurements_uc;  

  // save tab
  QPushButton *choose_path_prefix;
  QPushButton *save_grid_control;
  QPushButton *save_surface_points_control;

  QT_labeled_obj<QLineEdit> *path_prefix_control;
  
  // structure control
  QT_labeled_obj<QDoubleSpinBox> *uc_size_control_a;
  QT_labeled_obj<QDoubleSpinBox> *uc_size_control_c;
  QT_labeled_obj<QDoubleSpinBox> *channel_prop_control;
  QT_labeled_obj<QComboBox> *level_par_type;
  QT_labeled_obj<QComboBox> *surface_type_control;
  

  // resolution control
  QT_labeled_obj<QSpinBox> *z_points_control;
  QT_labeled_obj<QSpinBox> *x_points_control;
  QT_labeled_obj<QComboBox> *image_scaling_control;  
  QT_labeled_obj<QCheckBox> *invert_control;  
  QCheckBox *autoupdate_control;  


  // slice parameters control
  QT_labeled_obj<QDoubleSpinBox> *slice_thickness_control;
  QT_labeled_obj<QDoubleSpinBox> *slice_width_control;
  QT_labeled_obj<QDoubleSpinBox> *slice_height_control;  
  QT_labeled_obj<QDoubleSpinBox> *slice_position_control;

  QPushButton *button_set_to_uc_dim;
  QCheckBox *draw_uc_control;
  QVBoxLayout *auto_uc_layout;
  
  QT_labeled_obj<QSpinBox> *miller_h_control;
  QT_labeled_obj<QSpinBox> *miller_k_control;
  QT_labeled_obj<QSpinBox> *miller_l_control;  
  
  //membranes
  QLabel *membranes_label;
  QTableWidget *membranes_control;
  QPushButton *add_membrane_control;
  QPushButton *rm_membrane_control;

  QLabel *fill_channels_label;
  QScrollArea *fill_channels_control_container;
  QWidget *fill_channels_control_content;
  QVBoxLayout *fill_channels_container_layout;
  std::vector<QCheckBox*> fill_channels;
  std::vector<bool> channel_states;
  
  //basic buttons
  QPushButton *button_render;
  QPushButton *button_save;

  QPushButton *button_write_pars;
  QPushButton *button_read_pars;
  
  QPushButton *button_measure_vol_area;
  QPushButton *button_measure_percthres;

  // batch creation
  QLabel *batch_instructions;
  std::vector<QT_labeled_obj<QLineEdit>*> batch_values;
  QT_labeled_obj<QLineEdit> *batch_output_name;
  QT_labeled_obj<QPushButton> *batch_choose_folder;
  QProgressBar *batch_progress;
  QPushButton *batch_compute_start;
  QPushButton *batch_compute_stop;
  
  
  //status bar
  QStatusBar *status_bar;
  QLabel *status_bar_status_m;
  QLabel *status_bar_status_p;  
  QLabel *status_bar_uco;   // --> uc in orientation
  QLabel *status_bar_or; //--> orientation
  QLabel *status_bar_pixs; //--> pixel info
  QLabel *status_bar_uc; //-->> uc_dim


  //about
  QLabel *qt_logo, *qt_text;
#ifdef USE_CGAL
  QLabel *cgal_logo, *cgal_text;
  QHBoxLayout *about_cgal_layout;
#endif
  QLabel *about_us;
  QLabel *refs_ack;

  // license
  QTextEdit *licenses;

  
  // layouts

  // the main layout, holds status bar and sub_main_layout
  QHBoxLayout *main_layout;
  
  // slice properties
  QHBoxLayout *slice_settings;
  QVBoxLayout *slice_orientation_layout;
  QVBoxLayout *slice_dimension_layout;
  
  // hold draw area and status bar
  QVBoxLayout *draw_and_status_layout;
  
  QHBoxLayout *buttons_layout;
  QVBoxLayout *buttons_savewrite;
  QVBoxLayout *buttons_projection;
  
  QHBoxLayout *structure_settings;
  QHBoxLayout *surface_level_settings;
  QHBoxLayout *resolution_settings;

  QHBoxLayout *membrane_settings;
  QHBoxLayout *membrane_buttons_layout;
  
  QVBoxLayout *parameters_widget_layout;
  QVBoxLayout *measurement_widget_layout;
  QHBoxLayout *measurement_widget_buttons_layout;
  QVBoxLayout *save_widget_layout;

  
  QVBoxLayout *batch_widget_layout;
  QHBoxLayout *batch_widget_buttons_layout;
  QVBoxLayout *batch_widget_parameters_layout;
  QHBoxLayout *batch_widget_output_layout;
  
  QVBoxLayout *about_widget_layout;
  QHBoxLayout *about_qt_layout;

  QVBoxLayout *license_widget_layout;
  
  std::vector< QSpacerItem* > v_spacer;
  std::vector< QSpacerItem* > h_spacer;  
  QFrame *h_line_1, *h_line_2, *h_line_3;

  
  //
  // background members
  //
  
  const QImage *image;
  QPixmap *img_pix;
  unsigned char* img_data;

  QPainter *uc_artist;
  
  // the thread holding the projection class
  QThread *thread;

  // another copy of the projection class, but used for measurements
  QThread *t_stats;

  // a thread to handle the batch computation
  QThread *bc_thread;
  QThread *cb_thread;
  
  // the surface projection class
  sp_qt *sp;

  // surface projection class to compute measurements, not sure if
  // that's needed after all
  sp_qt *sp_stats;

  // the batch_creation controller class
  qt_bc *bc;

  // the callback functor to call the projection thingy
  gui_callback *cb;
  
  // a pointer to the application executing this form
  QApplication *app;
  
signals:


  void call_change_uc_size( double ab, double c );
  void call_change_hkl( int h, int k, int l );

  void call_compute_projection();
  void call_update_stats( QString what );


  void call_save_grid( QString fn );
  void call_save_surface( int id, QString fn );
  void call_save_network( int id, QString fn );

  void call_set_measurement_status( int state );

  void call_change_channel_color( int id, int val );

  void call_show_dialog_box();

  void request_parameter_change( double value );

  void emit_start_batch();
  void emit_stop_batch( bool go_on );			 
			 
public slots:

  void change_surface_par_type( int index );
  
  void check_channel_color();
  
  void update_unitcell_size();
  
  void update_gui_from_sp();
  
  
  void update_view();
  void draw_unitcell();

  
  void quit_app();
  void save_image_to_file();

  void write_membranes(int row, int col);
  void read_membranes();
  void add_membrane( double first = 0, double second = 0);
  void rm_membrane();

  // argument needed to match signal of QSpinBox
  void change_orientation( int val );

  void output_message( QString msg, int type = 1);

  void update_stats();
  void update_measurements();

  void change_autoupdate( int state );
  void request_compute_projection();

  void measure_vol_area();
  void measure_network();
  void measure_percolation();

  void choose_export_prefix();
  std::string get_prefix();
  void save_grid();
  void save_network();
  void save_surface_points();

  void set_state( int what, int state );

  void set_measurements_status( int state );

  void read_parameters();
  void write_parameters();
  
  void update_fill_channels();
  
  void set_parameter();


  bool check_batch_string_consistency( QString inp );
  std::vector<double> check_struct_type_string_consistency( QString inp );  
  void set_batch_output();
  void start_batch_computing();
  void stop_batch_computing();
  void finalize_batch_loop();
  
  void do_something();
  
};



