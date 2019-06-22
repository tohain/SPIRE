#include "sp_gui_imp.h"

sp_gui_imp::sp_gui_imp( wxWindow* parent )
:
sp_gui( parent )
{

  //initialize surface_projection object
  sp = new surface_projection();
  
  //set up gui

  //set range of the spinners
  phi_ctl->SetRange(0, 2*M_PI );
  theta_ctl->SetRange(0, M_PI );  

  //read availabe surfaces
  for( auto it : sp->get_surface_choices() ){
    type_ctl->Append( it );
  }
  type_ctl->SetSelection( 0 );

  //set up the gui with the parameters from sp object
  read_parameters();

  //do the first image
  compute_and_draw();
}


/**
 *  Writes all parameters to to sp class and calls necessary methods
 *  to update it
 */
void sp_gui_imp::write_parameters(){

  //ntucs
  sp->set_ntucs( ntucs_ctl->GetValue() );
  //a (unit cell length);
  sp->set_a( a_ctl->GetValue() );
  //type
  sp->set_type( type_ctl->GetSelection() );
  //mem_width
  sp->set_mem_width( d_ctl->GetValue() );

  if( ortype_ang_ctl->IsEnabled() ){
    
    //theta
    sp->set_theta( theta_ctl->GetValue() );
    //phi
    sp->set_phi( phi_ctl->GetValue() );

    sp->set_periodicity_length( -1 );
  } else {
    update_orientation_from_hkl();
  }
 
  //slice
  sp->set_slice_width( slicewidth_ctl->GetValue() );
  sp->set_slice_height( sliceheight_ctl->GetValue() );

  //points
  sp->set_n_points_x( n_points_xy_ctl->GetValue() );
  sp->set_n_points_y( n_points_xy_ctl->GetValue() );
  sp->set_n_points_z( n_points_z_ctl->GetValue() );
  
  //update gui
  dx_ctl->SetValue( sp->get_dx() );
  dy_ctl->SetValue( sp->get_dy() );
  dz_ctl->SetValue( sp->get_dz() );  
}


/**
 *  reads all parameters from the sp obejct and set all gui objects
 */
void sp_gui_imp::read_parameters(){
  //geometry
  ntucs_ctl->SetValue( sp->get_ntucs() );
  a_ctl->SetValue( sp->get_a() );
  type_ctl->SetSelection( sp->get_type() );
  d_ctl->SetValue( sp->get_mem_width() );

  theta_ctl->SetValue( sp->get_theta() );
  phi_ctl->SetValue( sp->get_phi() );
  h_ctl->SetValue( sp->get_h() );
  k_ctl->SetValue( sp->get_k() );
  l_ctl->SetValue( sp->get_l() );


  sliceheight_ctl->SetValue( sp->get_slice_height() );
  slicewidth_ctl->SetValue( sp->get_slice_width() );
  //check for ranges
  update_controls_periodicity();
 
    
  n_points_xy_ctl->SetValue( sp->get_width() );
  n_points_z_ctl->SetValue( sp->get_depth() );

  dx_ctl->SetValue( sp->get_dx() );
  dy_ctl->SetValue( sp->get_dy() );
  dz_ctl->SetValue( sp->get_dz() );
}

/**
 * This function sets the ranges of all controls, which are affected
 * by periodicty or aperiodicity of the chosen orientation
 */
void sp_gui_imp::update_controls_periodicity(){
  //we need to check, if we are periodic or not. So wether this value
  //is fractional and in [0,1] or unlimited
  if( sp->get_periodicity_length() == -1 ){
    //aperiodic. No limits :)
    slicewidth_ctl->SetRange(-9999, 9999);
  } else {
    //periodic. fractional
    slicewidth_ctl->SetRange(0, 1);
  }

  if( sp->get_periodicity_length() == -1 ){
    sliceheight_ctl->SetRange(-9999, 9999);
  } else {
    sliceheight_ctl->SetRange(0, 1);
  }
  
}


/*
 *  These methods control the GUI. Their puropse is primarly to update
 *  the sp class with their new parameters, update other GUI elements
 *  and start computations
 */

void sp_gui_imp::ntucs_change( wxSpinEvent& event )
{
  //update parameter
  sp->set_ntucs( ntucs_ctl->GetValue() );
  //update geometry
  sp->update_geometry();

  //update gui
  dx_ctl->SetValue( sp->get_dx() );
  dy_ctl->SetValue( sp->get_dy() );
  dz_ctl->SetValue( sp->get_dz() );  

  //periodicity
  sp->update_periodicity_length();
  update_controls_periodicity();
  
  //redraw  
  draw_preview();
}

void sp_gui_imp::a_change( wxSpinDoubleEvent& event )
{
  //update parameter
  sp->set_a( a_ctl->GetValue() );

  //update geometry
  sp->update_geometry();

  //update periodicity
  sp->update_periodicity_length();
  update_controls_periodicity();
  
  //update gui
  dx_ctl->SetValue( sp->get_dx() );
  dy_ctl->SetValue( sp->get_dy() );
  dz_ctl->SetValue( sp->get_dz() );  

  //redraw
  draw_preview();  
}

void sp_gui_imp::surface_change( wxCommandEvent& event )
{
  //update parameter
  sp->set_type( type_ctl->GetSelection() );

  //redraw
  draw_preview();
}

void sp_gui_imp::d_change( wxSpinDoubleEvent& event )
{
  //update parameter
  sp->set_mem_width( d_ctl->GetValue() );
  
  //redraw
  draw_preview();  
}

void sp_gui_imp::theta_change( wxSpinDoubleEvent& event )
{
  //update parameter
  sp->set_theta( theta_ctl->GetValue() );

  //update periodicity
  sp->update_periodicity_length();
  update_controls_periodicity();

  //redraw
  draw_preview();  
}

void sp_gui_imp::phi_change( wxSpinDoubleEvent& event )
{
  //update parameter
  sp->set_phi( phi_ctl->GetValue() );

  //update periodicity
  sp->update_periodicity_length();
  update_controls_periodicity();
	       
  //redraw
  draw_preview();  
}

void sp_gui_imp::selection_angles(  wxCommandEvent& event ){
  h_ctl->Enable(false);
  k_ctl->Enable(false);
  l_ctl->Enable(false);

  theta_ctl->Enable(true);
  phi_ctl->Enable(true);    
}

void sp_gui_imp::selection_miller(  wxCommandEvent& event ){
  theta_ctl->Enable(false);
  phi_ctl->Enable(false);

  h_ctl->Enable(true);
  k_ctl->Enable(true);
  l_ctl->Enable(true);
}

void sp_gui_imp::h_change( wxSpinEvent& event )
{
  if( h_ctl->GetValue() != 0 || k_ctl->GetValue() != 0 || l_ctl->GetValue() ){
    
    //update parameter
    sp->set_h( h_ctl->GetValue() );

    //update the angles
    sp->set_orientation_from_hkl();

    //get new periodicity
    sp->update_periodicity_length();

    //update the gui
    update_controls_periodicity();
    theta_ctl->SetValue( sp->get_theta() );
    phi_ctl->SetValue( sp->get_phi() );    
               
    //redraw
    draw_preview();
    
  } else {
    wxMessageBox( wxT("Invalid orientation. At least one indeces must be > 0"), wxT("Error"), wxICON_INFORMATION);
    h_ctl->SetValue(1);
 } 
}

void sp_gui_imp::k_change( wxSpinEvent& event )
{
  if( h_ctl->GetValue() != 0 || k_ctl->GetValue() != 0 || l_ctl->GetValue() ){
    
    //update parameter
    sp->set_k( k_ctl->GetValue() );

    //update the angles
    sp->set_orientation_from_hkl();

    //get new periodicity
    sp->update_periodicity_length();

    //update the gui
    update_controls_periodicity();
    theta_ctl->SetValue( sp->get_theta() );
    phi_ctl->SetValue( sp->get_phi() );    
               
    //redraw
    draw_preview();
    
  } else {
    wxMessageBox( wxT("Invalid orientation. At least one indeces must be > 0"), wxT("Error"), wxICON_INFORMATION);
    k_ctl->SetValue(1);
 }
}

void sp_gui_imp::l_change( wxSpinEvent& event )
{
  if( h_ctl->GetValue() != 0 || k_ctl->GetValue() != 0 || l_ctl->GetValue() ){
    
    //update parameter
    sp->set_l( l_ctl->GetValue() );

    //update the angles
    sp->set_orientation_from_hkl();

    //get new periodicity
    sp->update_periodicity_length();

    //update the gui
    update_controls_periodicity();
    theta_ctl->SetValue( sp->get_theta() );
    phi_ctl->SetValue( sp->get_phi() );
               
    //redraw
    draw_preview();
    
  } else {
    wxMessageBox( wxT("Invalid orientation. At least one indeces must be > 0"), wxT("Error"), wxICON_INFORMATION);
    l_ctl->SetValue(1);
 } 
}

void sp_gui_imp::slicewidth_change( wxSpinDoubleEvent& event )
{
  //update parameter
  sp->set_slice_width( slicewidth_ctl->GetValue() );

  //update geometry
  sp->update_geometry();

  //update gui
  dz_ctl->SetValue( sp->get_dz() );
  
  //redraw
  draw_preview();  
}

void sp_gui_imp::sliceheight_change( wxSpinDoubleEvent& event )
{
  //update parameter
  sp->set_slice_height( sliceheight_ctl->GetValue() );

  //redraw
  draw_preview();    
}

void sp_gui_imp::n_points_xy_change( wxSpinEvent& event ){
  //update parameter
  sp->set_n_points_x( n_points_xy_ctl->GetValue() );
  sp->set_n_points_y( n_points_xy_ctl->GetValue() );  

  //update_geometry
  sp->update_geometry();

  //update container size
  sp->update_containers();
  
  //update gui
  dx_ctl->SetValue( sp->get_dx() );
  dy_ctl->SetValue( sp->get_dy() );

  //redraw
  draw_preview();
}


void sp_gui_imp::n_points_z_change( wxSpinEvent& event ){

  //update parameter
  sp->set_n_points_z( n_points_z_ctl->GetValue() );  

  //update_geometry
  sp->update_geometry();

  //update container size
  sp->update_containers();

  //update gui
  dz_ctl->SetValue( sp->get_dz() );

  //redraw
  draw_preview();
}



void sp_gui_imp::invert_change( wxCommandEvent& event )
{
  //redraw
  draw_preview();    
}

void sp_gui_imp::button_render( wxCommandEvent& event )
{
  compute_and_draw();
}

void sp_gui_imp::button_save( wxCommandEvent& event )
{
  std::string fn (m_filePicker2->GetPath().c_str());
  if( fn == "" ){
    wxMessageBox( wxT("Please select a valid filename"),
		  wxT("Nothing to see here"),
		  wxICON_INFORMATION);
  } else {
    write_image( fn, sp->get_projection(), sp->get_width(), sp->get_height(), invert_ctl->GetValue() );
  }
}

void sp_gui_imp::button_quit( wxCommandEvent& event )
{
  Close( true );
}

void sp_gui_imp::update_orientation_from_hkl(){
  //get indeces
  int h = h_ctl->GetValue();
  int k = k_ctl->GetValue();
  int l = l_ctl->GetValue();  

  //check for valid orientation
  if( h > 0 || k > 0 || l > 0 ){
    //set indeces
    sp->set_h( h ); sp->set_k( k ); sp->set_l( l );

    //convert indeces to angles
    sp->set_orientation_from_hkl();

    //update gui
    theta_ctl->SetValue( sp->get_theta() );
    phi_ctl->SetValue( sp->get_phi() );
  } else {
    wxMessageBox( wxT("Invalid orientation. At least one indeces must be > 0"), wxT("Error"), wxICON_INFORMATION);
    h_ctl->SetValue(1);
  }
}


/**
 * Does the same as \ref compute_and_draw(), but checks first, if a
 * live preview is still wanted or the number of points are too high
 */
void sp_gui_imp::draw_preview(){
  if (max_prev_points_ctl->GetValue() > n_points_xy_ctl->GetValue() ){
    compute_and_draw();
  }
}


void sp_gui_imp::compute_and_draw(){

  //update parameters
  //  write_parameters();
  //update geometry and initialize points
  //  sp->update_geometry();

  //update periodiciy
  //  sp->update_periodicity_length();
  
  //  sp->update_containers();

  //compute the projection
  sp->compute_projection();
  unsigned char *proj = sp->get_image( invert_ctl->GetValue() );
  
  int w = sp->get_width();
  int h = sp->get_height();
  
  //convert the greyscale image array from sp to rgb array
  unsigned char* rgb_img = (unsigned char*) malloc( sizeof(unsigned char) * w * h * 3 );
  for(int ii=0; ii<w; ii++){
    for(int jj=0; jj<h; jj++){    
      int ind = ii*w + jj;      
            
      rgb_img[3*ind]=proj[ind];
      rgb_img[3*ind+1]=proj[ind];
      rgb_img[3*ind+2]=proj[ind];
      
    }
  }

  //now create an wxImage object, since that can be loaded into
  //awxBitmap, which then can be displayed in wxStaticBitmap. wxImage
  //also allows for scaling

  //delete old one if existant
  if( img != NULL ){
    delete( img );
  }

  //create wximage from array
  img = new wxImage( w, h, rgb_img );

  //rescale the image and draw it
  redraw();
}


void sp_gui_imp::redraw()
{
  //get the current size of the staticBitmap, which is the size we
  //want to scale to
  int sw, sh;
  m_bitmap1->GetSize(&sw, &sh);

  //scale to the smaller edge to keept aspect ratio
  int new_size = std::min( sw, sh);
  img->Rescale(new_size, new_size);

  //delete old memory
  if( preview != NULL ){
    delete( preview );
  }

  //create new bitmap
  preview = new wxBitmap( *img );
  //load that into the form
  m_bitmap1->SetBitmap( *preview );
}
