#include <QApplication>


#include "qt_gui.hpp"

int main( int argc, char *argv[]){

  QApplication app( argc, argv );

  GUI gui( &app );
  gui.show();
  
  return app.exec();

}
