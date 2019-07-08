﻿#ifndef __sp_gui_imp__
#define __sp_gui_imp__

/**
@file
Subclass of sp_gui, which is partly generated by wxFormBuilder.
*/

#include "sp_gui.h"

//// end generated include

#include <cstdlib>
#include <ctime>
#include "surface_projection.hpp"
#include "img_out.hpp"
#include <wx/wx.h>
#include "sp_gui_tooltips.h"

/** \brief The class is implementing the GUI of this software. 
 * 
 * While it does not hold the code to generate the forms, it
 * implements all the events and hold some more code to handle the
 * parameter input.
 * 
 */
class sp_gui_imp : public sp_gui
{
protected:
  // Handlers for sp_gui events.
  void focus_ntucs_ctl( wxFocusEvent& event );		
  void focus_a_ctl( wxFocusEvent& event );
  void focus_type_ctl( wxFocusEvent& event );

  void focus_theta_ctl( wxFocusEvent& event );
  void focus_phi_ctl( wxFocusEvent& event );
  void focus_h_ctl( wxFocusEvent& event );
  void focus_k_ctl( wxFocusEvent& event );
  void focus_l_ctl( wxFocusEvent& event );
  void focus_slicewidth_ctl( wxFocusEvent& event );
  void focus_sliceposition_ctl( wxFocusEvent& event );
  void focus_nrpointsxy_ctl( wxFocusEvent& event );
  void focus_npointsz_ctl( wxFocusEvent& event );
  void focus_updateprev_ctl( wxFocusEvent& event );
  void focus_invert_ctl( wxFocusEvent& event );
  void focus_filename_ctl( wxFocusEvent& event );
  void ntucs_change( wxSpinEvent& event );
  void a_change( wxSpinDoubleEvent& event );
  void surface_change( wxCommandEvent& event );
  void d_change( wxSpinDoubleEvent& event );
  void selection_angles( wxCommandEvent& event );
  void theta_change( wxSpinDoubleEvent& event );
  void phi_change( wxSpinDoubleEvent& event );
  void selection_miller( wxCommandEvent& event );
  void h_change( wxSpinEvent& event );
  void k_change( wxSpinEvent& event );
  void l_change( wxSpinEvent& event );
  void slicewidth_change( wxSpinDoubleEvent& event );
  void sliceheight_change( wxSpinDoubleEvent& event );
  void n_points_xy_change( wxSpinEvent& event );
  void n_points_z_change( wxSpinEvent& event );
  void invert_change( wxCommandEvent& event );
  void button_render( wxCommandEvent& event );
  void button_save( wxCommandEvent& event );  
  void button_quit( wxCommandEvent& event );
  void membrane_selected( wxListEvent& event );
  void mem_change_d( wxSpinDoubleEvent& event );
  void mem_change_w( wxSpinDoubleEvent& event );
  void membrane_add( wxCommandEvent& event );
  void membrane_delete( wxCommandEvent& event );  
public:
  /** Constructor */
  sp_gui_imp( wxWindow* parent );
  //// end generated class members
  
private:

  /// saves the currently selected membrane
  void membrane_save();
  
  /// Parse membranes into list view
  void read_membranes();
  
  /// Transfers all parameters from the GUI into the \ref sp obejct
  void write_parameters();

  /// Tansfers all parameters from the sp class to gui
  void read_parameters();

  /// Sets ranges of controls affected by periodicity
  void update_controls_periodicity();
  
  /// Redraws the currently loaded image after resizeing the window
  void redraw();

  /// Updates the orientation from miller indeces
  void update_orientation_from_hkl();  
  
  /// Computes the projection, scales it to the wxStaticBitmap size
  /// and draws it on screen
  void compute_and_draw();

  /// Draw a preview if resolution is low enough
  void draw_preview();

  //the bitmap to show
  wxBitmap *preview;
  //image object to convert 2D array to a bitmap and scale
  wxImage *img;

  //The object doing the calculations
  surface_projection *sp;

  //A class storing the help strings
  tooltips *help;
  
};

#endif // __sp_gui_imp__
