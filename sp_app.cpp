
#include "sp_gui_imp.h"

#include <cstring>


class sp_app : public wxApp {

public:
  virtual bool OnInit();
};


//main method and stuff
wxIMPLEMENT_APP(sp_app);

bool sp_app::OnInit(){
  sp_gui_imp *frame = new sp_gui_imp(NULL);
  frame->Show( true );
  return true;
}






