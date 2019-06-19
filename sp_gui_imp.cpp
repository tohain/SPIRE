#include "sp_gui_imp.h"

#include <cstdlib>
#include <ctime>

sp_gui_imp::sp_gui_imp( wxWindow* parent )
:
  sp_gui( parent )
{
  srand( time(NULL) );
}

void sp_gui_imp::compute( wxCommandEvent& event )
{
  drawImg();
}

void sp_gui_imp::close_app( wxCommandEvent& event )
{
  Close( true );
}



void sp_gui_imp::drawImg(){

  int width = 100;
  int height = 100;
  
  unsigned char* imgArray = (unsigned char*) malloc( sizeof(unsigned char) * width * height * 3 );

  char val = rand() % 255;

  button_compute->SetLabel( std::to_string(val) );
  
  memset( imgArray, val, sizeof(char) * width * height * 3 );

  wxImage img (  width, height, imgArray );
  
  preview = new wxBitmap( img );

  m_bitmap1->SetBitmap( *preview );
  
}
