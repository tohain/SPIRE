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
#include <wx/stattext.h>
#include <wx/statline.h>
#include <wx/spinctrl.h>
#include <wx/choice.h>
#include <wx/panel.h>
#include <wx/radiobut.h>
#include <wx/textctrl.h>
#include <wx/checkbox.h>
#include <wx/scrolwin.h>
#include <wx/filepicker.h>
#include <wx/button.h>
#include <wx/gauge.h>
#include <wx/frame.h>

///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// Class sp_gui
///////////////////////////////////////////////////////////////////////////////
class sp_gui : public wxFrame
{
	private:

	protected:
		wxStaticBitmap* m_bitmap1;
		wxScrolledWindow* m_scrolledWindow1;
		wxPanel* geometry;
		wxStaticText* m_staticText15;
		wxStaticLine* m_staticline2;
		wxStaticText* m_staticText8;
		wxSpinCtrl* ntucs_ctl;
		wxStaticText* m_staticText9;
		wxSpinCtrlDouble* a_ctl;
		wxStaticText* m_staticText11;
		wxChoice* type_ctl;
		wxStaticText* m_staticText20;
		wxSpinCtrlDouble* d_ctl;
		wxPanel* orientation;
		wxStaticText* m_staticText16;
		wxStaticLine* m_staticline3;
		wxRadioButton* ortype_ang_ctl;
		wxStaticText* m_staticText4;
		wxSpinCtrlDouble* theta_ctl;
		wxStaticText* m_staticText5;
		wxSpinCtrlDouble* phi_ctl;
		wxRadioButton* otype_miller_ctl;
		wxStaticText* m_staticText6;
		wxSpinCtrl* h_ctl;
		wxStaticText* m_staticText7;
		wxSpinCtrl* k_ctl;
		wxStaticText* m_staticText81;
		wxSpinCtrl* l_ctl;
		wxStaticText* m_staticText241;
		wxStaticLine* m_staticline61;
		wxTextCtrl* text_periodicitylength;
		wxStaticLine* m_staticline6;
		wxStaticText* m_staticText18;
		wxSpinCtrlDouble* slicewidth_ctl;
		wxStaticText* m_staticText19;
		wxSpinCtrlDouble* sliceheight_ctl;
		wxPanel* resolution;
		wxStaticText* m_staticText17;
		wxStaticLine* m_staticline4;
		wxStaticText* m_staticText91;
		wxSpinCtrl* n_points_xy_ctl;
		wxStaticText* m_staticText10;
		wxSpinCtrl* n_points_z_ctl;
		wxStaticText* m_staticText111;
		wxSpinCtrlDouble* dx_ctl;
		wxStaticText* m_staticText12;
		wxSpinCtrlDouble* dy_ctl;
		wxStaticText* m_staticText13;
		wxSpinCtrlDouble* dz_ctl;
		wxPanel* options;
		wxStaticText* m_staticText211;
		wxStaticLine* m_staticline51;
		wxCheckBox* invert_ctl;
		wxStaticText* m_staticText23;
		wxSpinCtrl* max_prev_points_ctl;
		wxStaticText* m_staticText24;
		wxPanel* control;
		wxStaticText* m_staticText21;
		wxStaticLine* m_staticline5;
		wxFilePickerCtrl* m_filePicker2;
		wxButton* b_render;
		wxButton* b_save;
		wxGauge* m_gauge1;
		wxButton* b_quit;
		wxTextCtrl* text_help;

		// Virtual event handlers, overide them in your derived class
		virtual void focus_ntucs_ctl( wxFocusEvent& event ) = 0;
		virtual void ntucs_change( wxSpinEvent& event ) = 0;
		virtual void focus_a_ctl( wxFocusEvent& event ) = 0;
		virtual void a_change( wxSpinDoubleEvent& event ) = 0;
		virtual void surface_change( wxCommandEvent& event ) = 0;
		virtual void focus_type_ctl( wxFocusEvent& event ) = 0;
		virtual void focus_d_ctl( wxFocusEvent& event ) = 0;
		virtual void d_change( wxSpinDoubleEvent& event ) = 0;
		virtual void selection_angles( wxCommandEvent& event ) = 0;
		virtual void focus_theta_ctl( wxFocusEvent& event ) = 0;
		virtual void theta_change( wxSpinDoubleEvent& event ) = 0;
		virtual void focus_phi_ctl( wxFocusEvent& event ) = 0;
		virtual void phi_change( wxSpinDoubleEvent& event ) = 0;
		virtual void selection_miller( wxCommandEvent& event ) = 0;
		virtual void focus_h_ctl( wxFocusEvent& event ) = 0;
		virtual void h_change( wxSpinEvent& event ) = 0;
		virtual void focus_k_ctl( wxFocusEvent& event ) = 0;
		virtual void k_change( wxSpinEvent& event ) = 0;
		virtual void focus_l_ctl( wxFocusEvent& event ) = 0;
		virtual void l_change( wxSpinEvent& event ) = 0;
		virtual void focus_slicewidth_ctl( wxFocusEvent& event ) = 0;
		virtual void slicewidth_change( wxSpinDoubleEvent& event ) = 0;
		virtual void focus_sliceposition_ctl( wxFocusEvent& event ) = 0;
		virtual void sliceheight_change( wxSpinDoubleEvent& event ) = 0;
		virtual void focus_nrpointsxy_ctl( wxFocusEvent& event ) = 0;
		virtual void n_points_xy_change( wxSpinEvent& event ) = 0;
		virtual void focus_npointsz_ctl( wxFocusEvent& event ) = 0;
		virtual void n_points_z_change( wxSpinEvent& event ) = 0;
		virtual void invert_change( wxCommandEvent& event ) = 0;
		virtual void focus_invert_ctl( wxFocusEvent& event ) = 0;
		virtual void focus_updateprev_ctl( wxFocusEvent& event ) = 0;
		virtual void button_render( wxCommandEvent& event ) = 0;
		virtual void button_save( wxCommandEvent& event ) = 0;
		virtual void button_quit( wxCommandEvent& event ) = 0;


	public:

		sp_gui( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Surface Projection GUI"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 1094,569 ), long style = wxDEFAULT_FRAME_STYLE|wxRESIZE_BORDER|wxTAB_TRAVERSAL );

		~sp_gui();

};

