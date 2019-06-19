#include "sp_gui_imp.h"

sp_gui_imp::sp_gui_imp( wxWindow* parent )
:
sp_gui( parent )
{

}

void sp_gui_imp::compute( wxCommandEvent& event )
{
  wxMessageBox( "Running Computations",
		"Doing the computions", wxOK | wxICON_INFORMATION );
}

void sp_gui_imp::close_app( wxCommandEvent& event )
{
  Close( true );
}

