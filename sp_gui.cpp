///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Jun 18 2019)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#include "sp_gui.h"

///////////////////////////////////////////////////////////////////////////

sp_gui::sp_gui( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxFrame( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	wxBoxSizer* bSizer2;
	bSizer2 = new wxBoxSizer( wxHORIZONTAL );

	wxBoxSizer* bSizer3;
	bSizer3 = new wxBoxSizer( wxVERTICAL );

	m_bitmap1 = new wxStaticBitmap( this, wxID_ANY, wxNullBitmap, wxDefaultPosition, wxDefaultSize, 0 );
	bSizer3->Add( m_bitmap1, 0, wxALL, 5 );


	bSizer2->Add( bSizer3, 1, wxEXPAND, 5 );

	wxBoxSizer* bSizer4;
	bSizer4 = new wxBoxSizer( wxVERTICAL );

	button_compute = new wxButton( this, wxID_ANY, wxT("Compute"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer4->Add( button_compute, 0, wxALL, 5 );

	button_close = new wxButton( this, wxID_ANY, wxT("Close"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer4->Add( button_close, 0, wxALL, 5 );


	bSizer2->Add( bSizer4, 1, wxEXPAND, 5 );


	this->SetSizer( bSizer2 );
	this->Layout();

	this->Centre( wxBOTH );

	// Connect Events
	button_compute->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::compute ), NULL, this );
	button_close->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::close_app ), NULL, this );
}

sp_gui::~sp_gui()
{
	// Disconnect Events
	button_compute->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::compute ), NULL, this );
	button_close->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::close_app ), NULL, this );

}
