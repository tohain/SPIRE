
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

#include <QThread>
#include <QMutex>

#include "qt_sp.hpp"
#include "sp_gui_tooltips.h"


template <class QT_O>
class QT_labeled_obj {

public:
  QT_labeled_obj( std::string label_str, QWidget *parent = NULL ){
    obj = new QT_O ( parent );
    lbl = new QLabel( QString(label_str.c_str()), parent );
    lyt = new QVBoxLayout();

    lyt->addWidget( lbl );
    lyt->addWidget( obj );
  }
  
  ~QT_labeled_obj(){
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





class GUI : public QWidget {

  Q_OBJECT

public:
  explicit GUI( QApplication *app, QWidget *parent = 0 );


  double progress, progress_stats;
  std::string status, status_stats;
  
protected:
  
  void paintEvent( QPaintEvent * event );
  
private:

  //tooltips
  tooltips ttips;

  // functions
  void set_up_ui();
  void set_up_signals_and_slots();
  void set_up_tooltips();
  
  // ui elements
  QLabel *draw_area;

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

  QT_labeled_obj<QLineEdit> *path_prefix_control;
  
  // structure control
  QT_labeled_obj<QSpinBox> *ntucs_control;
  QT_labeled_obj<QDoubleSpinBox> *uc_size_control_a;
  QT_labeled_obj<QDoubleSpinBox> *uc_size_control_c;
  QT_labeled_obj<QDoubleSpinBox> *channel_prop_control;
  QT_labeled_obj<QComboBox> *surface_type_control;
  

  // resolution control
  QT_labeled_obj<QSpinBox> *z_points_control;
  QT_labeled_obj<QSpinBox> *xy_points_control;
  QT_labeled_obj<QComboBox> *image_scaling_control;  
  QLabel *pix_size_indicator;
  QCheckBox *invert_control;  
  QCheckBox *autoupdate_control;  


  // slice parameters control
  QT_labeled_obj<QDoubleSpinBox> *slice_width_control;
  QT_labeled_obj<QDoubleSpinBox> *slice_position_control;

  QT_labeled_obj<QSpinBox> *miller_h_control;
  QT_labeled_obj<QSpinBox> *miller_k_control;
  QT_labeled_obj<QSpinBox> *miller_l_control;  
  
  //membranes
  QLabel *membranes_label;
  QTableWidget *membranes_control;
  QPushButton *add_membrane_control;
  QPushButton *rm_membrane_control;

  
  //basic buttons
  QPushButton *button_quit;
  QPushButton *button_measure;
  QPushButton *button_render;
  QPushButton *button_save;


  //status bar
  QStatusBar *status_bar;
  QLabel *status_bar_status_m;
  QLabel *status_bar_status_p;  
  QLabel *status_bar_pixs;
  QLabel *status_bar_vols;
  QLabel *status_bar_areas;
  QLabel *status_bar_mins;

  // layouts


  
  QVBoxLayout *main_layout;

  QHBoxLayout *sub_main_layout;
  QHBoxLayout *status_bar_layout;
  
  QHBoxLayout *buttons_layout;
  
  QHBoxLayout *structure_settings;
  QHBoxLayout *resolution_settings;
  QHBoxLayout *slice_settings;
  QHBoxLayout *membrane_settings;
  QVBoxLayout *membrane_buttons_layout;
  
  QVBoxLayout *controls_basic_layout;
  QVBoxLayout *controls_measurement_layout;
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

  void call_compute_projection();
  void call_update_stats();


  void call_save_grid( QString fn );
  void call_save_surface( int id, QString fn );
  void call_save_network( int id, QString fn );

  void call_set_measurement_status( int state );


  void call_show_dialog_box();
					       
public slots:



  
  void update_gui_from_sp();
  
  void update_status( QString s);
  
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

  void measure();

  void choose_export_prefix();
  std::string get_prefix();
  void save_grid();
  void save_network();
  void save_surface_points();

  void set_state( int what, int state );

  void set_measurements_status( int state );

  
  void do_something();
  
};



