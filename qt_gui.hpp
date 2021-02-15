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

#include <iostream>
#include <iomanip>
#include <string>
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

#include <QtSvg/QSvgWidget>

#include <QThread>
#include <QMutex>

#include "qt_sp.hpp"
#include "sp_gui_tooltips.h"
#include "slice_orientation_visualisation.hpp"
#include "vec_mat_math.hpp"

#include "global_settings.hpp"

template <class QT_O>
class QT_v_labeled_obj {

public:
  QT_v_labeled_obj( std::string label_str, QWidget *parent = NULL ){
    obj = new QT_O ( parent );
    lbl = new QLabel( QString(label_str.c_str()), parent );
    lyt = new QVBoxLayout();

    lyt->addWidget( lbl );
    lyt->addWidget( obj );
  }
  
  ~QT_v_labeled_obj(){
    delete( obj );
    delete( lbl );
    delete( lyt );
  }

  QVBoxLayout* layout(){
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
  QVBoxLayout *lyt;
};

template <class QT_O>
class QT_h_labeled_obj {

public:
  QT_h_labeled_obj( std::string label_str, QWidget *parent = NULL, bool _left = true ){
    obj = new QT_O ( parent );
    lbl = new QLabel( QString(label_str.c_str()), parent );
    lyt = new QHBoxLayout();
    left = _left;

    if( left ){
      lyt->addWidget( lbl );
      lyt->addWidget( obj );
    } else {
      lyt->addWidget( obj );
      lyt->addWidget( lbl );
    }
  }
  
  ~QT_h_labeled_obj(){
    delete( obj );
    delete( lbl );
    delete( lyt );
  }

  QHBoxLayout* layout(){
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

  bool left; // is label left or right?
  
  QHBoxLayout *lyt;  
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
  QWidget *save_widget;
  QWidget *about_widget;
  QWidget *license_widget;  

  // manual tab
  QPlainTextEdit *manual;

  // stats tab
  //QLabel *detailled_stats;
  QTableWidget *detailled_stats;

  // save tab
  QPushButton *choose_path_prefix;
  QPushButton *save_grid_control;
  QPushButton *save_surface_points_control;

  QT_v_labeled_obj<QLineEdit> *path_prefix_control;
  
  // structure control
  QT_v_labeled_obj<QDoubleSpinBox> *uc_size_control_a;
  QT_v_labeled_obj<QDoubleSpinBox> *uc_size_control_c;
  QT_v_labeled_obj<QDoubleSpinBox> *channel_prop_control;
  QT_v_labeled_obj<QComboBox> *level_par_type;
  QT_v_labeled_obj<QComboBox> *surface_type_control;
  

  // resolution control
  QT_v_labeled_obj<QSpinBox> *z_points_control;
  QT_v_labeled_obj<QSpinBox> *x_points_control;
  QT_v_labeled_obj<QComboBox> *image_scaling_control;  
  QT_v_labeled_obj<QCheckBox> *invert_control;  
  QCheckBox *autoupdate_control;  


  // slice parameters control
  QT_h_labeled_obj<QDoubleSpinBox> *slice_thickness_control;
  QT_h_labeled_obj<QDoubleSpinBox> *slice_width_control;
  QT_h_labeled_obj<QDoubleSpinBox> *slice_height_control;  
  QT_h_labeled_obj<QDoubleSpinBox> *slice_position_control;

  QPushButton *button_set_to_uc_dim;
  QCheckBox *draw_uc_control;
  QVBoxLayout *auto_uc_layout;
  
  QT_h_labeled_obj<QSpinBox> *miller_h_control;
  QT_h_labeled_obj<QSpinBox> *miller_k_control;
  QT_h_labeled_obj<QSpinBox> *miller_l_control;  
  
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

  //status bar
  QStatusBar *status_bar;
  QLabel *status_bar_status_m;
  QLabel *status_bar_status_p;  
  //QLabel *status_bar_pixs;
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

  
  // the surface projection class
  sp_qt *sp;

  sp_qt *sp_stats;

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
  void update_detailled_stats();

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

  void do_something();
  
};



