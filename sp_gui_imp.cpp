#include "sp_gui_imp.h"
#include "img_out.hpp"

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

void sp_gui_imp::repaint( wxCommandEvent& event )
{
  int sw, sh;
  m_bitmap1->GetSize(&sw, &sh);
  std::cout << sw << ", " << sh << std::endl;

  int new_size = std::min( sw, sh);
  img->Rescale(new_size, new_size);
  
  preview = new wxBitmap( *img );
  
  m_bitmap1->SetBitmap( *preview );
}



void sp_gui_imp::drawImg(){



  
  surface_projection sp;
  sp.compute_projection();
  unsigned char *proj = sp.get_image();

  int width = sp.get_width();
  int height = sp.get_height();
  
  unsigned char* imgArray = (unsigned char*) malloc( sizeof(unsigned char) * width * height * 3 );
  for(int ii=0; ii<width; ii++){
    for(int jj=0; jj<height; jj++){    
      int ind = ii*width + jj;

      imgArray[3*ind]=proj[ind];
      imgArray[3*ind+1]=proj[ind];
      imgArray[3*ind+2]=proj[ind];      
    }
  }


  img = new wxImage(  width, height, imgArray );

  int sw, sh;
  m_bitmap1->GetSize(&sw, &sh);
  std::cout << sw << ", " << sh << std::endl;

  int new_size = std::min( sw, sh);
  img->Rescale(new_size, new_size);
  
  preview = new wxBitmap( *img );

  m_bitmap1->SetBitmap( *preview );
  


}
