#include <wx/wx.h>
#include "sp_gui_imp.h"


class sp_gui_app : public wxApp{
public:
  virtual bool OnInit();
};

bool sp_gui_app::OnInit(){
  sp_gui_imp *frame = new sp_gui_imp(NULL);
  frame->Show( true );
  return true;
}


wxIMPLEMENT_APP(sp_gui_app);
