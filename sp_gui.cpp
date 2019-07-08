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
	this->SetSizeHints( wxSize( 750,500 ), wxDefaultSize );

	wxBoxSizer* root;
	root = new wxBoxSizer( wxHORIZONTAL );

	wxGridSizer* img_preview;
	img_preview = new wxGridSizer( 1, 1, 0, 0 );

	img_preview->SetMinSize( wxSize( 300,300 ) );
	m_bitmap1 = new wxStaticBitmap( this, wxID_ANY, wxNullBitmap, wxDefaultPosition, wxSize( -1,-1 ), 0 );
	img_preview->Add( m_bitmap1, 0, wxALL|wxEXPAND, 5 );


	root->Add( img_preview, 4, wxEXPAND, 1 );

	wxBoxSizer* controls;
	controls = new wxBoxSizer( wxVERTICAL );

	m_scrolledWindow1 = new wxScrolledWindow( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxHSCROLL|wxVSCROLL );
	m_scrolledWindow1->SetScrollRate( 5, 5 );
	wxBoxSizer* bSizer17;
	bSizer17 = new wxBoxSizer( wxVERTICAL );

	geometry = new wxPanel( m_scrolledWindow1, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxBoxSizer* bSizer27;
	bSizer27 = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* bSizer29;
	bSizer29 = new wxBoxSizer( wxHORIZONTAL );

	m_staticText15 = new wxStaticText( geometry, wxID_ANY, wxT("Geometry"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText15->Wrap( -1 );
	bSizer29->Add( m_staticText15, 0, wxALL, 5 );

	m_staticline2 = new wxStaticLine( geometry, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer29->Add( m_staticline2, 1, wxEXPAND | wxALL, 5 );


	bSizer27->Add( bSizer29, 0, wxEXPAND, 5 );

	wxBoxSizer* bSizer101;
	bSizer101 = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* dimension;
	dimension = new wxBoxSizer( wxHORIZONTAL );


	dimension->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText8 = new wxStaticText( geometry, wxID_ANY, wxT("Nr. of Unit Cells"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
	m_staticText8->Wrap( -1 );
	dimension->Add( m_staticText8, 0, wxALL, 5 );

	ntucs_ctl = new wxSpinCtrl( geometry, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 999, 1 );
	dimension->Add( ntucs_ctl, 0, wxALL, 5 );


	dimension->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText9 = new wxStaticText( geometry, wxID_ANY, wxT("Unit Cell Size"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
	m_staticText9->Wrap( -1 );
	dimension->Add( m_staticText9, 0, wxALL, 5 );

	a_ctl = new wxSpinCtrlDouble( geometry, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 999, 1, 0.01 );
	a_ctl->SetDigits( 2 );
	dimension->Add( a_ctl, 0, wxALL, 5 );


	dimension->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer101->Add( dimension, 1, wxEXPAND, 5 );

	wxBoxSizer* surface_type;
	surface_type = new wxBoxSizer( wxHORIZONTAL );


	surface_type->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText11 = new wxStaticText( geometry, wxID_ANY, wxT("Structure Type"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
	m_staticText11->Wrap( -1 );
	surface_type->Add( m_staticText11, 0, wxALL, 5 );

	wxArrayString type_ctlChoices;
	type_ctl = new wxChoice( geometry, wxID_ANY, wxDefaultPosition, wxDefaultSize, type_ctlChoices, 0 );
	type_ctl->SetSelection( 0 );
	surface_type->Add( type_ctl, 0, wxALL, 5 );


	surface_type->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText27 = new wxStaticText( geometry, wxID_ANY, wxT("Channel Proportion"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText27->Wrap( -1 );
	surface_type->Add( m_staticText27, 0, wxALL, 5 );

	level_set_ctrl = new wxSpinCtrlDouble( geometry, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, -999, 999, 0, 0.001 );
	level_set_ctrl->SetDigits( 4 );
	surface_type->Add( level_set_ctrl, 0, wxALL, 5 );


	surface_type->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer101->Add( surface_type, 1, wxEXPAND, 5 );


	bSizer27->Add( bSizer101, 0, wxEXPAND, 5 );

	wxBoxSizer* bSizer381;
	bSizer381 = new wxBoxSizer( wxHORIZONTAL );

	m_staticText26 = new wxStaticText( geometry, wxID_ANY, wxT("MyLabeMembranesl"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText26->Wrap( -1 );
	bSizer381->Add( m_staticText26, 0, wxALL, 5 );

	m_staticline8 = new wxStaticLine( geometry, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer381->Add( m_staticline8, 1, wxEXPAND | wxALL, 5 );


	bSizer27->Add( bSizer381, 1, wxEXPAND, 5 );

	wxBoxSizer* bSizer351;
	bSizer351 = new wxBoxSizer( wxHORIZONTAL );

	wxBoxSizer* bSizer371;
	bSizer371 = new wxBoxSizer( wxVERTICAL );

	membranes_ctrl = new wxListCtrl( geometry, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLC_REPORT|wxLC_SINGLE_SEL|wxLC_SORT_ASCENDING );
	bSizer371->Add( membranes_ctrl, 1, wxALL, 5 );

	wxBoxSizer* bSizer363;
	bSizer363 = new wxBoxSizer( wxHORIZONTAL );

	membrane_idx = new wxSpinCtrl( geometry, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 10, 0 );
	membrane_idx->Enable( false );

	bSizer363->Add( membrane_idx, 0, wxALL, 5 );

	membrane_dist = new wxSpinCtrlDouble( geometry, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, -100, 100, 0, 0.01 );
	membrane_dist->SetDigits( 4 );
	bSizer363->Add( membrane_dist, 0, wxALL, 5 );

	membrane_width = new wxSpinCtrlDouble( geometry, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 100, 0, 0.01 );
	membrane_width->SetDigits( 4 );
	bSizer363->Add( membrane_width, 0, wxALL, 5 );


	bSizer371->Add( bSizer363, 0, wxEXPAND, 5 );


	bSizer351->Add( bSizer371, 0, wxEXPAND, 5 );

	wxBoxSizer* bSizer362;
	bSizer362 = new wxBoxSizer( wxVERTICAL );

	m_button5 = new wxButton( geometry, wxID_ANY, wxT("Add"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer362->Add( m_button5, 0, wxALL, 5 );

	m_button6 = new wxButton( geometry, wxID_ANY, wxT("Delete"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer362->Add( m_button6, 0, wxALL, 5 );


	bSizer351->Add( bSizer362, 0, wxEXPAND, 5 );


	bSizer351->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer27->Add( bSizer351, 0, wxEXPAND, 5 );


	geometry->SetSizer( bSizer27 );
	geometry->Layout();
	bSizer27->Fit( geometry );
	bSizer17->Add( geometry, 0, wxEXPAND | wxALL, 5 );

	orientation = new wxPanel( m_scrolledWindow1, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxBoxSizer* bSizer6;
	bSizer6 = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* bSizer30;
	bSizer30 = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* bSizer31;
	bSizer31 = new wxBoxSizer( wxHORIZONTAL );

	m_staticText16 = new wxStaticText( orientation, wxID_ANY, wxT("Orientation"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText16->Wrap( -1 );
	bSizer31->Add( m_staticText16, 0, wxALL, 5 );

	m_staticline3 = new wxStaticLine( orientation, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer31->Add( m_staticline3, 1, wxEXPAND | wxALL, 5 );


	bSizer30->Add( bSizer31, 0, wxEXPAND, 5 );


	bSizer6->Add( bSizer30, 0, wxEXPAND, 5 );

	wxBoxSizer* bSizer312;
	bSizer312 = new wxBoxSizer( wxHORIZONTAL );

	ortype_ang_ctl = new wxRadioButton( orientation, wxID_ANY, wxT("Angles"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer312->Add( ortype_ang_ctl, 0, wxALL, 5 );


	bSizer312->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer6->Add( bSizer312, 1, wxEXPAND, 5 );

	wxBoxSizer* angles;
	angles = new wxBoxSizer( wxHORIZONTAL );


	angles->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText4 = new wxStaticText( orientation, wxID_ANY, wxT("Theta"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText4->Wrap( -1 );
	angles->Add( m_staticText4, 0, wxALL, 5 );

	theta_ctl = new wxSpinCtrlDouble( orientation, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxSP_WRAP, 0, 100, 0, 0.01 );
	theta_ctl->SetDigits( 7 );
	theta_ctl->Enable( false );

	angles->Add( theta_ctl, 0, wxALL, 5 );


	angles->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText5 = new wxStaticText( orientation, wxID_ANY, wxT("Phi"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText5->Wrap( -1 );
	angles->Add( m_staticText5, 0, wxALL, 5 );

	phi_ctl = new wxSpinCtrlDouble( orientation, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxSP_WRAP, 0, 100, 0, 0.01 );
	phi_ctl->SetDigits( 7 );
	phi_ctl->Enable( false );

	angles->Add( phi_ctl, 0, wxALL, 5 );


	angles->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer6->Add( angles, 0, wxEXPAND, 5 );

	wxBoxSizer* bSizer322;
	bSizer322 = new wxBoxSizer( wxHORIZONTAL );

	otype_miller_ctl = new wxRadioButton( orientation, wxID_ANY, wxT("Miller Indeces"), wxDefaultPosition, wxDefaultSize, 0 );
	otype_miller_ctl->SetValue( true );
	bSizer322->Add( otype_miller_ctl, 0, wxALL, 5 );


	bSizer322->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer6->Add( bSizer322, 1, wxEXPAND, 5 );

	wxBoxSizer* hkl;
	hkl = new wxBoxSizer( wxHORIZONTAL );


	hkl->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText6 = new wxStaticText( orientation, wxID_ANY, wxT("h"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText6->Wrap( -1 );
	hkl->Add( m_staticText6, 0, wxALL, 5 );

	h_ctl = new wxSpinCtrl( orientation, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxSP_WRAP, -999, 999, 0 );
	hkl->Add( h_ctl, 0, wxALL, 5 );


	hkl->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText7 = new wxStaticText( orientation, wxID_ANY, wxT("k"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText7->Wrap( -1 );
	hkl->Add( m_staticText7, 0, wxALL, 5 );

	k_ctl = new wxSpinCtrl( orientation, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxSP_WRAP, -999, 999, 0 );
	hkl->Add( k_ctl, 0, wxALL, 5 );


	hkl->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText81 = new wxStaticText( orientation, wxID_ANY, wxT("l"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText81->Wrap( -1 );
	hkl->Add( m_staticText81, 0, wxALL, 5 );

	l_ctl = new wxSpinCtrl( orientation, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxSP_WRAP, -999, 999, 1 );
	hkl->Add( l_ctl, 0, wxALL, 5 );


	hkl->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer6->Add( hkl, 0, wxEXPAND, 5 );

	wxBoxSizer* bSizer3311;
	bSizer3311 = new wxBoxSizer( wxHORIZONTAL );

	m_staticText241 = new wxStaticText( orientation, wxID_ANY, wxT("Periodicity Length/Unit cell size in current orientation"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText241->Wrap( -1 );
	bSizer3311->Add( m_staticText241, 0, wxALL, 5 );

	m_staticline61 = new wxStaticLine( orientation, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer3311->Add( m_staticline61, 1, wxEXPAND | wxALL, 5 );


	bSizer6->Add( bSizer3311, 1, wxEXPAND, 5 );

	wxBoxSizer* bSizer361;
	bSizer361 = new wxBoxSizer( wxHORIZONTAL );

	text_periodicitylength = new wxTextCtrl( orientation, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxTE_CENTER );
	text_periodicitylength->Enable( false );

	bSizer361->Add( text_periodicitylength, 1, wxALL|wxEXPAND, 5 );


	bSizer6->Add( bSizer361, 0, wxEXPAND, 5 );

	wxBoxSizer* bSizer331;
	bSizer331 = new wxBoxSizer( wxHORIZONTAL );

	m_staticline6 = new wxStaticLine( orientation, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer331->Add( m_staticline6, 1, wxEXPAND | wxALL, 5 );


	bSizer6->Add( bSizer331, 1, wxEXPAND, 5 );

	wxBoxSizer* bSizer35;
	bSizer35 = new wxBoxSizer( wxHORIZONTAL );


	bSizer35->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText18 = new wxStaticText( orientation, wxID_ANY, wxT("Slice width"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText18->Wrap( -1 );
	bSizer35->Add( m_staticText18, 0, wxALL, 5 );

	slicewidth_ctl = new wxSpinCtrlDouble( orientation, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 999, 0.2, 0.01 );
	slicewidth_ctl->SetDigits( 3 );
	bSizer35->Add( slicewidth_ctl, 0, wxALL, 5 );


	bSizer35->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText19 = new wxStaticText( orientation, wxID_ANY, wxT("Slice position"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText19->Wrap( -1 );
	bSizer35->Add( m_staticText19, 0, wxALL, 5 );

	sliceheight_ctl = new wxSpinCtrlDouble( orientation, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 1, 0, 0.01 );
	sliceheight_ctl->SetDigits( 3 );
	bSizer35->Add( sliceheight_ctl, 0, wxALL, 5 );


	bSizer35->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer6->Add( bSizer35, 0, wxEXPAND, 5 );


	orientation->SetSizer( bSizer6 );
	orientation->Layout();
	bSizer6->Fit( orientation );
	bSizer17->Add( orientation, 0, wxEXPAND | wxALL, 5 );

	resolution = new wxPanel( m_scrolledWindow1, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxBoxSizer* bSizer10;
	bSizer10 = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* bSizer32;
	bSizer32 = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* bSizer33;
	bSizer33 = new wxBoxSizer( wxHORIZONTAL );

	m_staticText17 = new wxStaticText( resolution, wxID_ANY, wxT("Resolution"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText17->Wrap( -1 );
	bSizer33->Add( m_staticText17, 0, wxALL, 5 );

	m_staticline4 = new wxStaticLine( resolution, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer33->Add( m_staticline4, 1, wxEXPAND | wxALL, 5 );


	bSizer32->Add( bSizer33, 1, wxEXPAND, 5 );


	bSizer10->Add( bSizer32, 0, wxEXPAND, 5 );

	wxBoxSizer* n_points;
	n_points = new wxBoxSizer( wxHORIZONTAL );


	n_points->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText91 = new wxStaticText( resolution, wxID_ANY, wxT("Nr. of points xy"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText91->Wrap( -1 );
	n_points->Add( m_staticText91, 0, wxALL, 5 );

	n_points_xy_ctl = new wxSpinCtrl( resolution, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxSP_WRAP, 1, 1000000, 150 );
	n_points->Add( n_points_xy_ctl, 0, wxALL, 5 );


	n_points->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText10 = new wxStaticText( resolution, wxID_ANY, wxT("Nr. of points z"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText10->Wrap( -1 );
	n_points->Add( m_staticText10, 0, wxALL, 5 );

	n_points_z_ctl = new wxSpinCtrl( resolution, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxSP_WRAP, 1, 1000000, 50 );
	n_points->Add( n_points_z_ctl, 0, wxALL, 5 );


	n_points->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer10->Add( n_points, 1, wxEXPAND, 5 );

	wxBoxSizer* distance;
	distance = new wxBoxSizer( wxHORIZONTAL );


	distance->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText111 = new wxStaticText( resolution, wxID_ANY, wxT("dx"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText111->Wrap( -1 );
	distance->Add( m_staticText111, 0, wxALL, 5 );

	dx_ctl = new wxSpinCtrlDouble( resolution, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 100, 0, 1 );
	dx_ctl->SetDigits( 3 );
	dx_ctl->Enable( false );

	distance->Add( dx_ctl, 0, wxALL, 5 );


	distance->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText12 = new wxStaticText( resolution, wxID_ANY, wxT("dy"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText12->Wrap( -1 );
	distance->Add( m_staticText12, 0, wxALL, 5 );

	dy_ctl = new wxSpinCtrlDouble( resolution, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 100, 0, 1 );
	dy_ctl->SetDigits( 3 );
	dy_ctl->Enable( false );

	distance->Add( dy_ctl, 0, wxALL, 5 );


	distance->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText13 = new wxStaticText( resolution, wxID_ANY, wxT("dz"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText13->Wrap( -1 );
	distance->Add( m_staticText13, 0, wxALL, 5 );

	dz_ctl = new wxSpinCtrlDouble( resolution, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 100, 0, 1 );
	dz_ctl->SetDigits( 3 );
	dz_ctl->Enable( false );

	distance->Add( dz_ctl, 0, wxALL, 5 );


	distance->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer10->Add( distance, 1, wxEXPAND, 5 );


	resolution->SetSizer( bSizer10 );
	resolution->Layout();
	bSizer10->Fit( resolution );
	bSizer17->Add( resolution, 0, wxEXPAND | wxALL, 5 );

	options = new wxPanel( m_scrolledWindow1, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxBoxSizer* bSizer291;
	bSizer291 = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* bSizer311;
	bSizer311 = new wxBoxSizer( wxHORIZONTAL );

	m_staticText211 = new wxStaticText( options, wxID_ANY, wxT("Options"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText211->Wrap( -1 );
	bSizer311->Add( m_staticText211, 0, wxALL, 5 );

	m_staticline51 = new wxStaticLine( options, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer311->Add( m_staticline51, 1, wxEXPAND | wxALL, 5 );


	bSizer291->Add( bSizer311, 0, wxEXPAND, 5 );

	wxBoxSizer* bSizer321;
	bSizer321 = new wxBoxSizer( wxHORIZONTAL );


	bSizer321->Add( 0, 0, 1, wxEXPAND, 5 );

	invert_ctl = new wxCheckBox( options, wxID_ANY, wxT("Invert Colors"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer321->Add( invert_ctl, 0, wxALL, 5 );


	bSizer321->Add( 0, 0, 1, wxEXPAND, 5 );

	m_staticText23 = new wxStaticText( options, wxID_ANY, wxT("Auto update upto"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText23->Wrap( -1 );
	bSizer321->Add( m_staticText23, 0, wxALL, 5 );

	max_prev_points_ctl = new wxSpinCtrl( options, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 999, 100 );
	bSizer321->Add( max_prev_points_ctl, 0, wxALL, 5 );

	m_staticText24 = new wxStaticText( options, wxID_ANY, wxT("Points"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText24->Wrap( -1 );
	bSizer321->Add( m_staticText24, 0, wxALL, 5 );


	bSizer321->Add( 0, 0, 1, wxEXPAND, 5 );


	bSizer291->Add( bSizer321, 0, wxEXPAND, 5 );


	options->SetSizer( bSizer291 );
	options->Layout();
	bSizer291->Fit( options );
	bSizer17->Add( options, 1, wxEXPAND | wxALL, 5 );


	bSizer17->Add( 0, 0, 1, wxEXPAND, 5 );


	m_scrolledWindow1->SetSizer( bSizer17 );
	m_scrolledWindow1->Layout();
	bSizer17->Fit( m_scrolledWindow1 );
	controls->Add( m_scrolledWindow1, 1, wxEXPAND | wxALL, 5 );

	control = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxBoxSizer* bSizer36;
	bSizer36 = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* bSizer37;
	bSizer37 = new wxBoxSizer( wxHORIZONTAL );

	m_staticText21 = new wxStaticText( control, wxID_ANY, wxT("Control"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText21->Wrap( -1 );
	bSizer37->Add( m_staticText21, 0, wxALL, 5 );

	m_staticline5 = new wxStaticLine( control, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer37->Add( m_staticline5, 1, wxEXPAND | wxALL, 5 );


	bSizer36->Add( bSizer37, 0, wxEXPAND, 5 );

	wxBoxSizer* bSizer38;
	bSizer38 = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* bSizer40;
	bSizer40 = new wxBoxSizer( wxHORIZONTAL );

	m_filePicker2 = new wxFilePickerCtrl( control, wxID_ANY, wxEmptyString, wxT("Select a file"), wxT("*.*"), wxDefaultPosition, wxDefaultSize, wxFLP_SAVE|wxFLP_USE_TEXTCTRL );
	bSizer40->Add( m_filePicker2, 1, wxALL|wxEXPAND, 5 );

	b_render = new wxButton( control, wxID_ANY, wxT("Render"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer40->Add( b_render, 0, wxALL, 5 );

	b_save = new wxButton( control, wxID_ANY, wxT("Save"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer40->Add( b_save, 0, wxALL, 5 );


	bSizer38->Add( bSizer40, 1, wxEXPAND, 5 );

	wxBoxSizer* bSizer26;
	bSizer26 = new wxBoxSizer( wxHORIZONTAL );

	m_gauge1 = new wxGauge( control, wxID_ANY, 100, wxDefaultPosition, wxDefaultSize, wxGA_HORIZONTAL );
	m_gauge1->SetValue( 0 );
	bSizer26->Add( m_gauge1, 1, wxALL, 5 );


	bSizer38->Add( bSizer26, 1, wxEXPAND, 5 );

	wxBoxSizer* bSizer41;
	bSizer41 = new wxBoxSizer( wxHORIZONTAL );

	b_quit = new wxButton( control, wxID_ANY, wxT("Close"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer41->Add( b_quit, 1, wxALL|wxEXPAND, 5 );


	bSizer38->Add( bSizer41, 1, wxEXPAND, 5 );


	bSizer36->Add( bSizer38, 1, wxEXPAND, 5 );


	control->SetSizer( bSizer36 );
	control->Layout();
	bSizer36->Fit( control );
	controls->Add( control, 0, wxEXPAND | wxALL, 5 );


	root->Add( controls, 3, wxEXPAND, 5 );

	wxGridSizer* gSizer2;
	gSizer2 = new wxGridSizer( 1, 1, 0, 0 );

	text_help = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxTE_BESTWRAP|wxTE_MULTILINE|wxTE_NO_VSCROLL|wxTE_WORDWRAP );
	text_help->Enable( false );

	gSizer2->Add( text_help, 0, wxALL|wxEXPAND, 5 );


	root->Add( gSizer2, 1, wxEXPAND, 5 );


	this->SetSizer( root );
	this->Layout();
	m_menubar1 = new wxMenuBar( 0 );
	m_menu1 = new wxMenu();
	wxMenuItem* m_menuItem2;
	m_menuItem2 = new wxMenuItem( m_menu1, wxID_ANY, wxString( wxT("Quit") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu1->Append( m_menuItem2 );

	m_menubar1->Append( m_menu1, wxT("File") );

	m_menu2 = new wxMenu();
	wxMenuItem* m_menuItem1;
	m_menuItem1 = new wxMenuItem( m_menu2, wxID_ANY, wxString( wxT("About") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu2->Append( m_menuItem1 );

	m_menubar1->Append( m_menu2, wxT("Help") );

	this->SetMenuBar( m_menubar1 );


	this->Centre( wxBOTH );

	// Connect Events
	ntucs_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_ntucs_ctl ), NULL, this );
	ntucs_ctl->Connect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::ntucs_change ), NULL, this );
	a_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_a_ctl ), NULL, this );
	a_ctl->Connect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::a_change ), NULL, this );
	type_ctl->Connect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( sp_gui::surface_change ), NULL, this );
	type_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_type_ctl ), NULL, this );
	level_set_ctrl->Connect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::level_change ), NULL, this );
	membranes_ctrl->Connect( wxEVT_COMMAND_LIST_ITEM_SELECTED, wxListEventHandler( sp_gui::membrane_selected ), NULL, this );
	membrane_dist->Connect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::mem_change_d ), NULL, this );
	membrane_width->Connect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::mem_change_w ), NULL, this );
	m_button5->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::membrane_add ), NULL, this );
	m_button6->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::membrane_delete ), NULL, this );
	ortype_ang_ctl->Connect( wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( sp_gui::selection_angles ), NULL, this );
	theta_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_theta_ctl ), NULL, this );
	theta_ctl->Connect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::theta_change ), NULL, this );
	phi_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_phi_ctl ), NULL, this );
	phi_ctl->Connect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::phi_change ), NULL, this );
	otype_miller_ctl->Connect( wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( sp_gui::selection_miller ), NULL, this );
	h_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_h_ctl ), NULL, this );
	h_ctl->Connect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::h_change ), NULL, this );
	k_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_k_ctl ), NULL, this );
	k_ctl->Connect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::k_change ), NULL, this );
	l_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_l_ctl ), NULL, this );
	l_ctl->Connect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::l_change ), NULL, this );
	slicewidth_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_slicewidth_ctl ), NULL, this );
	slicewidth_ctl->Connect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::slicewidth_change ), NULL, this );
	sliceheight_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_sliceposition_ctl ), NULL, this );
	sliceheight_ctl->Connect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::sliceheight_change ), NULL, this );
	n_points_xy_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_nrpointsxy_ctl ), NULL, this );
	n_points_xy_ctl->Connect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::n_points_xy_change ), NULL, this );
	n_points_z_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_npointsz_ctl ), NULL, this );
	n_points_z_ctl->Connect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::n_points_z_change ), NULL, this );
	invert_ctl->Connect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( sp_gui::invert_change ), NULL, this );
	invert_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_invert_ctl ), NULL, this );
	max_prev_points_ctl->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_updateprev_ctl ), NULL, this );
	m_filePicker2->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_filename_ctl ), NULL, this );
	b_render->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::button_render ), NULL, this );
	b_save->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::button_save ), NULL, this );
	b_quit->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::button_quit ), NULL, this );
}

sp_gui::~sp_gui()
{
	// Disconnect Events
	ntucs_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_ntucs_ctl ), NULL, this );
	ntucs_ctl->Disconnect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::ntucs_change ), NULL, this );
	a_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_a_ctl ), NULL, this );
	a_ctl->Disconnect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::a_change ), NULL, this );
	type_ctl->Disconnect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( sp_gui::surface_change ), NULL, this );
	type_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_type_ctl ), NULL, this );
	level_set_ctrl->Disconnect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::level_change ), NULL, this );
	membranes_ctrl->Disconnect( wxEVT_COMMAND_LIST_ITEM_SELECTED, wxListEventHandler( sp_gui::membrane_selected ), NULL, this );
	membrane_dist->Disconnect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::mem_change_d ), NULL, this );
	membrane_width->Disconnect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::mem_change_w ), NULL, this );
	m_button5->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::membrane_add ), NULL, this );
	m_button6->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::membrane_delete ), NULL, this );
	ortype_ang_ctl->Disconnect( wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( sp_gui::selection_angles ), NULL, this );
	theta_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_theta_ctl ), NULL, this );
	theta_ctl->Disconnect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::theta_change ), NULL, this );
	phi_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_phi_ctl ), NULL, this );
	phi_ctl->Disconnect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::phi_change ), NULL, this );
	otype_miller_ctl->Disconnect( wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( sp_gui::selection_miller ), NULL, this );
	h_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_h_ctl ), NULL, this );
	h_ctl->Disconnect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::h_change ), NULL, this );
	k_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_k_ctl ), NULL, this );
	k_ctl->Disconnect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::k_change ), NULL, this );
	l_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_l_ctl ), NULL, this );
	l_ctl->Disconnect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::l_change ), NULL, this );
	slicewidth_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_slicewidth_ctl ), NULL, this );
	slicewidth_ctl->Disconnect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::slicewidth_change ), NULL, this );
	sliceheight_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_sliceposition_ctl ), NULL, this );
	sliceheight_ctl->Disconnect( wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler( sp_gui::sliceheight_change ), NULL, this );
	n_points_xy_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_nrpointsxy_ctl ), NULL, this );
	n_points_xy_ctl->Disconnect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::n_points_xy_change ), NULL, this );
	n_points_z_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_npointsz_ctl ), NULL, this );
	n_points_z_ctl->Disconnect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler( sp_gui::n_points_z_change ), NULL, this );
	invert_ctl->Disconnect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( sp_gui::invert_change ), NULL, this );
	invert_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_invert_ctl ), NULL, this );
	max_prev_points_ctl->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_updateprev_ctl ), NULL, this );
	m_filePicker2->Disconnect( wxEVT_SET_FOCUS, wxFocusEventHandler( sp_gui::focus_filename_ctl ), NULL, this );
	b_render->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::button_render ), NULL, this );
	b_save->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::button_save ), NULL, this );
	b_quit->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( sp_gui::button_quit ), NULL, this );

}
