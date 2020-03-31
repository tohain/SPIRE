
#include <iostream>
#include <string>
#include <vector>

#include <QWidget>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QPainter>
#include <QImage>
#include <QLabel>
#include <QPixmap>
#include <QApplication>


#include <QThread>
#include <QMutex>

#include "surface_projection.hpp"


class worker : public QObject, public surface_projection {

  Q_OBJECT;
  
public:
  worker (double &p, std::string &s) : surface_projection(p, s){
  }
  ~worker(){
    std::cout << "destructor called" << std::endl;
  }
  void requestWork();
  void abort();

private:


signals:
  void workRequested();
  void finished();

public slots:
  void doWork();  

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

  void setup_sp_object();

  
  const QImage *image;
  QPixmap *img_pix;
  
  QLabel *draw_area;
  
  QPushButton *button_quit;
  QPushButton *button_save;
  QPushButton *button_render;

  QPushButton *tmp;
  
  QHBoxLayout *main_layout;

  QVBoxLayout *buttons_layout;


  QThread *thread;


  worker *sp;


  QApplication *app;
  
signals:

    
public slots:


  void pushed_button();
  void pushed_button_2();
  void pushed_button_3();    
  
};
