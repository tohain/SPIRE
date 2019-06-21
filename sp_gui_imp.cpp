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

 
  //set up the first image
  update_parameters();
  update_orientation_from_hkl();
  
  compute_and_draw();
}


void sp_gui_imp::update_parameters(){

  //ntucs
  sp->set_ntucs( ntucs_ctl->GetValue() );
  //a
  sp->set_a( a_ctl->GetValue() );
  //type
  sp->set_type( std::string( type_ctl->GetString( type_ctl->GetSelection() ).c_str() ) );
  //mem_width
  sp->set_mem_width( d_ctl->GetValue() );
  //theta
  sp->set_theta( theta_ctl->GetValue() );
  //theta
  sp->set_theta( theta_ctl->GetValue() );
  //phi
  sp->set_phi( phi_ctl->GetValue() );
  //skip hkl for now

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

void sp_gui_imp::ntucs_change( wxSpinEvent& event )
{
  draw_preview();
}

void sp_gui_imp::a_change( wxSpinDoubleEvent& event )
{
  draw_preview();  
}

void sp_gui_imp::surface_change( wxCommandEvent& event )
{
  draw_preview();  
}

void sp_gui_imp::d_change( wxSpinDoubleEvent& event )
{
  draw_preview();  
}

void sp_gui_imp::theta_change( wxSpinDoubleEvent& event )
{
  draw_preview();  
}

void sp_gui_imp::phi_change( wxSpinDoubleEvent& event )
{
  draw_preview();  
}

void sp_gui_imp::h_change( wxSpinEvent& event )
{
  draw_preview();  
  update_orientation_from_hkl();
}

void sp_gui_imp::k_change( wxSpinEvent& event )
{
  update_orientation_from_hkl();
  draw_preview();    
}

void sp_gui_imp::l_change( wxSpinEvent& event )
{
  update_orientation_from_hkl();
  draw_preview();    
}

void sp_gui_imp::slicewidth_change( wxSpinDoubleEvent& event )
{
  draw_preview();  
}

void sp_gui_imp::sliceheight_change( wxSpinDoubleEvent& event )
{
  draw_preview();    
}

void sp_gui_imp::invert_change( wxCommandEvent& event )
{
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
  if( h != 0 || k != 0 || l != 0 ){
    //convert indeces to angles
    sp->set_orientation_from_hkl( h, k, l );

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
  update_parameters();
  //update geometry and initialize points
  sp->update_geometry();

  //update periodiciy
  sp->update_periodicity_length();
  
  sp->update_containers();

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
