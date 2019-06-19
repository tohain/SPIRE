///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Jun 18 2019)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#pragma once

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/statbmp.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/string.h>
#include <wx/sizer.h>
#include <wx/button.h>
#include <wx/frame.h>
#include <wx/wx.h>

///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// Class sp_gui
///////////////////////////////////////////////////////////////////////////////
class sp_gui : public wxFrame
{
	private:

	protected:
		wxStaticBitmap* m_bitmap1;
		wxButton* button_compute;
		wxButton* button_close;

		// Virtual event handlers, overide them in your derived class
		virtual void compute( wxCommandEvent& event ) = 0;
		virtual void close_app( wxCommandEvent& event ) = 0;


	public:

		sp_gui( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Surface Projection GUI"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 500,300 ), long style = wxDEFAULT_FRAME_STYLE|wxTAB_TRAVERSAL );

		~sp_gui();

};

