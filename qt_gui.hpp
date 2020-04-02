
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
#include <QTabWidget>

#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QSpacerItem>

#include <QThread>
#include <QMutex>

#include "qt_sp.hpp"

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


  double progress;
  std::string status;
  
protected:
  
  void paintEvent( QPaintEvent * event );
  
private:


  // functions
  void set_up_ui();
  void set_up_signals_and_slots();
  
  // ui elements
  QLabel *draw_area;

  QTabWidget *controls;
  QWidget *controls_basic;
  QWidget *controls_advanced;
  QWidget *controls_all;


  QT_labeled_obj<QSpinBox> *ntucs_control;
  QT_labeled_obj<QDoubleSpinBox> *uc_size_control;
  QT_labeled_obj<QDoubleSpinBox> *channel_prop_control;
  QT_labeled_obj<QComboBox> *surface_type_control;
  

  QT_labeled_obj<QSpinBox> *z_points_control;
  QT_labeled_obj<QSpinBox> *xy_points_control;
  QLabel *pix_size_indicator;
  
  QPushButton *button_quit;
  QPushButton *button_save;
  QPushButton *button_render;

  QHBoxLayout *main_layout;
  QHBoxLayout *buttons_layout;
  
  QHBoxLayout *structure_settings;

  QHBoxLayout *resolution_settings;
  
  QVBoxLayout *controls_basic_layout;
  QVBoxLayout *controls_advanced_layout;
  QVBoxLayout *controls_all_layout;    


  QSpacerItem *spacer;
  

  //
  // background members
  //
  
  const QImage *image;
  QPixmap *img_pix;
  
  // the thread holding the projection class
  QThread *thread;

  // the surface projection class
  sp_qt *sp;

  // a pointer to the application executing this form
  QApplication *app;
  
signals:

    
public slots:

  void update_gui_from_sp();
  
  void update_status( QString s);
  
  void update_view();
  void pushed_button_2();
  void pushed_button_3();    
  
};
