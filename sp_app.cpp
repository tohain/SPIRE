#include "sp_app.h"


bool sp_app::OnInit(){
  sp_gui_imp *frame = new sp_gui_imp(NULL);
  frame->Show( true );
  return true;
}



wxIMPLEMENT_APP(sp_app);

