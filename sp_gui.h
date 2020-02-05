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
#include <wx/listctrl.h>
#include <wx/button.h>
#include <wx/panel.h>
#include <wx/radiobut.h>
#include <wx/textctrl.h>
#include <wx/checkbox.h>
#include <wx/scrolwin.h>
#include <wx/filepicker.h>
#include <wx/gauge.h>
#include <wx/menu.h>
#include <wx/statusbr.h>
#include <wx/timer.h>
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
		wxStaticText* m_staticText27;
		wxSpinCtrlDouble* level_set_ctrl;
		wxStaticText* m_staticText26;
		wxStaticLine* m_staticline8;
		wxListCtrl* membranes_ctrl;
		wxSpinCtrl* membrane_idx;
		wxSpinCtrlDouble* membrane_dist;
		wxSpinCtrlDouble* membrane_width;
		wxButton* m_button5;
		wxButton* m_button6;
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
		wxMenuBar* m_menubar1;
		wxMenu* m_menu1;
		wxMenu* m_menu2;
		wxStatusBar* status;
		wxTimer timer;

		// Virtual event handlers, overide them in your derived class
		virtual void focus_ntucs_ctl( wxFocusEvent& event ) = 0;
		virtual void ntucs_change( wxSpinEvent& event ) = 0;
		virtual void focus_a_ctl( wxFocusEvent& event ) = 0;
		virtual void a_change( wxSpinDoubleEvent& event ) = 0;
		virtual void surface_change( wxCommandEvent& event ) = 0;
		virtual void focus_type_ctl( wxFocusEvent& event ) = 0;
		virtual void focus_level_set( wxFocusEvent& event ) = 0;
		virtual void level_change( wxSpinDoubleEvent& event ) = 0;
		virtual void membrane_selected( wxListEvent& event ) = 0;
		virtual void focus_membrane_list( wxFocusEvent& event ) = 0;
		virtual void focus_mem_dst( wxFocusEvent& event ) = 0;
		virtual void mem_change_d( wxSpinDoubleEvent& event ) = 0;
		virtual void focus_mem_width( wxFocusEvent& event ) = 0;
		virtual void mem_change_w( wxSpinDoubleEvent& event ) = 0;
		virtual void membrane_add( wxCommandEvent& event ) = 0;
		virtual void focus_button_add( wxFocusEvent& event ) = 0;
		virtual void membrane_delete( wxCommandEvent& event ) = 0;
		virtual void focus_button_delete( wxFocusEvent& event ) = 0;
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
		virtual void focus_filename_ctl( wxFocusEvent& event ) = 0;
		virtual void button_render( wxCommandEvent& event ) = 0;
		virtual void button_save( wxCommandEvent& event ) = 0;
		virtual void button_quit( wxCommandEvent& event ) = 0;
		virtual void timer_tick( wxTimerEvent& event ) = 0;


	public:

		sp_gui( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Surface Projection GUI"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 1281,624 ), long style = wxDEFAULT_FRAME_STYLE|wxRESIZE_BORDER|wxTAB_TRAVERSAL );

		~sp_gui();

};

