#include <QApplication>


#include "qt_gui.hpp"

int main( int argc, char *argv[]){

  QApplication app( argc, argv );

  QLocale default_locale = QLocale();
		       
  app.setWindowIcon(QIcon(":/resources/icon/icon.ico"));
  
  GUI gui( &app, &default_locale );
  gui.show();
  
  return app.exec();

}
