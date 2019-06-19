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

	wxGridSizer* gSizer3;
	gSizer3 = new wxGridSizer( 1, 1, 0, 0 );

	gSizer3->SetMinSize( wxSize( 100,100 ) );
	m_bitmap1 = new wxStaticBitmap( this, wxID_ANY, wxNullBitmap, wxDefaultPosition, wxSize( -1,-1 ), 0 );
	gSizer3->Add( m_bitmap1, 0, wxALL|wxEXPAND, 5 );


	bSizer2->Add( gSizer3, 1, wxEXPAND, 5 );

	wxBoxSizer* bSizer4;
	bSizer4 = new wxBoxSizer( wxVERTICAL );

	bSizer4->SetMinSize( wxSize( 100,100 ) );
	button_compute = new wxButton( this, wxID_ANY, wxT("Compute"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer4->Add( button_compute, 0, wxALL, 5 );

	button_repaint = new wxButton( this, wxID_ANY, wxT("Repaint"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer4->Add( button_repaint, 0, wxALL, 5 );

	button_close = new wxButton( this, wxID_ANY, wxT("Close"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer4->Add( button_close, 0, wxALL, 5 );


	bSizer2->Add( bSizer4, 1, wxEXPAND, 5 );


	this->SetSizer( bSizer2 );
	this->Layout();
	bSizer2->Fit( this );

	this->Centre( wxBOTH );

	// Connect Events
	button_compute->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::compute ), NULL, this );
	button_repaint->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::repaint ), NULL, this );
	button_close->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::close_app ), NULL, this );
}

sp_gui::~sp_gui()
{
	// Disconnect Events
	button_compute->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::compute ), NULL, this );
	button_repaint->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::repaint ), NULL, this );
	button_close->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::close_app ), NULL, this );

}
