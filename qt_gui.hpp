
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

#include <QtSvg/QSvgWidget>

#include <QThread>
#include <QMutex>

#include "qt_sp.hpp"
#include "sp_gui_tooltips.h"
#include "slice_orientation_visualisation.hpp"




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


  double progress, progress_stats;
  std::string status, status_stats;
  
protected:
  
  void paintEvent( QPaintEvent * event );
  
private:

  int sp_state, sp_stat_state;
  
  
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
  QWidget *controls_basic;
  QWidget *controls_measurement;
  QWidget *controls_save;
  QWidget *manual_widget;

  // manual tab
  QPlainTextEdit *manual;

  // stats tab
  QLabel *detailled_stats;

  // save tab
  QPushButton *choose_path_prefix;
  QPushButton *save_grid_control;
  QPushButton *save_surface_points_control;
  QPushButton *save_topological_network_control;

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
  QLabel *pix_size_indicator;
  QCheckBox *invert_control;  
  QCheckBox *autoupdate_control;  


  // slice parameters control
  QT_h_labeled_obj<QDoubleSpinBox> *slice_thickness_control;
  QT_h_labeled_obj<QDoubleSpinBox> *slice_width_control;
  QT_h_labeled_obj<QDoubleSpinBox> *slice_height_control;  
  QT_h_labeled_obj<QDoubleSpinBox> *slice_position_control;

  QPushButton *button_set_to_uc_dim;
  
  QT_h_labeled_obj<QSpinBox> *miller_h_control;
  QT_h_labeled_obj<QSpinBox> *miller_k_control;
  QT_h_labeled_obj<QSpinBox> *miller_l_control;  
  
  //membranes
  QLabel *membranes_label;
  QTableWidget *membranes_control;
  QPushButton *add_membrane_control;
  QPushButton *rm_membrane_control;

  QScrollArea *fill_channels_control_container;
  QWidget *fill_channels_control_content;
  QVBoxLayout *fill_channels_container_layout;
  std::vector<QCheckBox*> fill_channels;
  std::vector<bool> channel_states;
  
  //basic buttons
  QPushButton *button_quit;
  QPushButton *button_render;
  QPushButton *button_save;

  QPushButton *button_write_pars;
  QPushButton *button_read_pars;
  
  QPushButton *button_measure_vol_area;
  QPushButton *button_measure_network;

  QPushButton *button_measure_max_rad_dist;
  QPushButton *button_measure_channel_width_dist;  

  //status bar
  QStatusBar *status_bar;
  QLabel *status_bar_status_m;
  QLabel *status_bar_status_p;  
  QLabel *status_bar_pixs;
  QLabel *status_bar_vols;
  QLabel *status_bar_areas;
  QLabel *status_bar_mins;

  // layouts

  // the main layout, holds status bar and sub_main_layout
  QVBoxLayout *main_layout;

  // the status bar at the bottom of the window
  QHBoxLayout *status_bar_layout;
  
  // slice properties
  QHBoxLayout *slice_settings;
  QVBoxLayout *slice_orientation_layout;
  QVBoxLayout *slice_dimension_layout;
  
  // holds the tab bar and the drawing area
  QHBoxLayout *sub_main_layout;


  
  QHBoxLayout *buttons_layout;
  
  QHBoxLayout *structure_settings;
  QHBoxLayout *resolution_settings;

  QHBoxLayout *membrane_settings;
  QVBoxLayout *membrane_buttons_layout;
  
  QVBoxLayout *controls_basic_layout;
  QVBoxLayout *controls_measurement_layout;
  QHBoxLayout *controls_measurement_buttons_layout;
  QVBoxLayout *controls_save_layout;

  QVBoxLayout *manual_widget_layout;

  QSpacerItem *v_spacer;
  QSpacerItem *h_spacer;  


  //
  // background members
  //
  
  const QImage *image;
  QPixmap *img_pix;
  unsigned char* img_data;
  
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



