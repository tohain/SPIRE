#include "sp_gui_imp.h"

/** \brief Default constructor
 *
 * Sets a few element properties which are read from the \ref
 * surface_projection class and draws a first preview
 */
sp_gui_imp::sp_gui_imp( wxWindow* parent ) : sp_gui( parent ) {
  
  //set up status
  progress = 1;
  status_string = "Ready";

  //initialize surface_projection object
  sp = new surface_projection(progress, status_string);

  //set img empty, will be initialized later on
  img = NULL;
  preview = NULL;
  
  //set up gui

  //status bar
  status->SetFieldsCount( 2 );
  status->SetStatusText( status_string, 0 );
  status->SetStatusText( "100%", 1 );  
  
  //set up columns of list view
  membranes_ctrl->AppendColumn( "ID" );
  membranes_ctrl->AppendColumn( "DIST" );
  membranes_ctrl->AppendColumn( "WIDTH" );  
  read_membranes();
  membranes_ctrl->SetItemState(0, wxLIST_STATE_SELECTED, wxLIST_STATE_SELECTED);
  
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

  //add handlers to save image in differen formats
  img->AddHandler( new wxPNGHandler );
  img->AddHandler( new wxTIFFHandler );

  help = new tooltips();
}


/** 
 * parses the membranes from the membranes vector in
 * the surface projection object and views the data
 * in the list
 */
void sp_gui_imp::read_membranes(){

  //clear all items
  membranes_ctrl->DeleteAllItems();
  
  for( unsigned int ii=0; ii<sp->get_membranes().size(); ii+=2 ){
    membranes_ctrl->InsertItem( int(ii/2.), std::to_string(int(ii/2.)) );
    membranes_ctrl->SetItem( int(ii/2.), 1, std::to_string(sp->get_membranes()[ii]) );
    membranes_ctrl->SetItem( int(ii/2.), 2, std::to_string(sp->get_membranes()[ii+1]) );
  }
}



/**
 *  Writes all parameters from the GUI into the \ref sp class and
 *  calls necessary methods to update the GUI and the \ref sp class
 */
void sp_gui_imp::write_parameters(){

  //ntucs
  sp->set_ntucs( ntucs_ctl->GetValue() );
  //a (unit cell length);
  sp->set_a( a_ctl->GetValue() );
  //type
  sp->set_type( type_ctl->GetSelection() );


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
 *  reads all parameters from the sp obejct and updates GUI controls
 */
void sp_gui_imp::read_parameters(){
  //geometry
  ntucs_ctl->SetValue( sp->get_ntucs() );
  a_ctl->SetValue( sp->get_a() );
  type_ctl->SetSelection( sp->get_type() );

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
    slicewidth_ctl->SetRange(0, 9999);
    sliceheight_ctl->SetRange(-9999, 9999);
    
    //write to indicator
    text_periodicitylength->SetValue("Not periodic");
  } else {
    //periodic. fractional
    slicewidth_ctl->SetRange(0, 1);
    sliceheight_ctl->SetRange(0, 1);

    //write to indicator
    text_periodicitylength->SetValue("");    
    *text_periodicitylength << sp->get_a()*sp->get_periodicity_length();
  }
  
}


void sp_gui_imp::focus_level_set( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->level_set_tooltip;
}

void sp_gui_imp::focus_membrane_list( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->membrane_list_tooltip;
}
  
void sp_gui_imp::focus_mem_dst( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->mem_dst_tooltip;
}

void sp_gui_imp::focus_mem_width( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->mem_width_tooltip;
}

void sp_gui_imp::focus_button_delete( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->mem_delete_tooltip;
}

void sp_gui_imp::focus_button_add( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->mem_add_tooltip;
}

void sp_gui_imp::timer_tick( wxTimerEvent& event ){

  //update status GUI elemnts
  status->SetStatusText( status_string, 0 );
  
  std::cout << "tick\n";

}


void sp_gui_imp::focus_ntucs_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  if(sp->get_periodicity_length() == -1 ){
    *text_help << help->ntucs_tooltip;
  } else {
    *text_help << help->ntucs_p_tooltip;
  }  
}

void sp_gui_imp::focus_filename_ctl( wxFocusEvent& event ){
 text_help->SetValue("");
 *text_help << help->filename_tooltip;
}


void sp_gui_imp::focus_a_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->a_tooltip;
}

void sp_gui_imp::focus_theta_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->theta_tooltip;  
}

void sp_gui_imp::focus_phi_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->phi_tooltip;  
}

void sp_gui_imp::focus_h_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->hkl_tooltip;  
}

void sp_gui_imp::focus_k_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->hkl_tooltip;  
}

void sp_gui_imp::focus_l_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->hkl_tooltip;  
}

void sp_gui_imp::focus_slicewidth_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  if(sp->get_periodicity_length() == -1 ){
    *text_help << help->slicewidth_tooltip;
  } else {
    *text_help << help->slicewidth_p_tooltip;
  }
}

void sp_gui_imp::focus_sliceposition_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  if(sp->get_periodicity_length() == -1 ){
    *text_help << help->sliceheight_tooltip;
  } else {
    *text_help << help->sliceheight_p_tooltip;
  }
}

void sp_gui_imp::focus_nrpointsxy_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->nr_points_xy_tooltip;  
}

void sp_gui_imp::focus_npointsz_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->nr_points_z_tooltip;  
}

void sp_gui_imp::focus_invert_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->invert_tooltip;  
}

void sp_gui_imp::focus_updateprev_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->update_prev_tooltip;  
}

void sp_gui_imp::focus_type_ctl( wxFocusEvent& event ){
  text_help->SetValue("");
  *text_help << help->type_tooltip;  
}




/*
 *  These methods control the GUI. Their puropse is primarly to update
 *  the sp class with their new parameters, update other GUI elements
 *  and start computations
 */

/**
 * Handles a change of number in unit cells. Needs to update geometry,
 * periodicity and GUI
 */
void sp_gui_imp::ntucs_change( wxSpinEvent& event )
{
  //update parameter
  try {
    sp->set_ntucs( ntucs_ctl->GetValue() );
  } catch( invalid_parameter_exception e ){
    wxMessageBox( e.details(),
		  e.what(),
		  wxICON_INFORMATION);
    //get corrected value
    ntucs_ctl->SetValue( sp->get_ntucs() );
  }
  
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
  try {
    sp->set_a( a_ctl->GetValue() );
  } catch ( invalid_parameter_exception e ){
    wxMessageBox( e.details(),
		  e.what(),
		  wxICON_INFORMATION);
    //get corrected value
    a_ctl->SetValue( sp->get_a() );
  }


  //update periodicity
  sp->update_periodicity_length();
  update_controls_periodicity();

  //update geometry
  sp->update_geometry();
  
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
  try {
    sp->set_type( type_ctl->GetSelection() );
  } catch ( invalid_parameter_exception e ){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION );
    type_ctl->SetSelection( sp->get_type() );
  }

  //redraw
  draw_preview();
}


void sp_gui_imp::theta_change( wxSpinDoubleEvent& event )
{
  //update parameter
  try {
    sp->set_theta( theta_ctl->GetValue() );
  } catch( invalid_parameter_exception e ){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION);
    theta_ctl->SetValue( sp->get_theta() );
  }

  //update periodicity
  sp->update_periodicity_length();
  update_controls_periodicity();

  //redraw
  draw_preview();  
}

void sp_gui_imp::phi_change( wxSpinDoubleEvent& event )
{
  //update parameter
  try {
    sp->set_phi( phi_ctl->GetValue() );
  } catch( invalid_parameter_exception e ){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION);
    phi_ctl->SetValue( sp->get_phi() );
  } 

  //update periodicity
  sp->update_periodicity_length();
  update_controls_periodicity();
	       
  //redraw
  draw_preview();  
}

/**
 * Disables the hkl controls, enables angle controls
 */
void sp_gui_imp::selection_angles(  wxCommandEvent& event ){
  h_ctl->Enable(false);
  k_ctl->Enable(false);
  l_ctl->Enable(false);

  theta_ctl->Enable(true);
  phi_ctl->Enable(true);    
}

/**
 * Enables the hkl controls, disables angle controls
 */
void sp_gui_imp::selection_miller(  wxCommandEvent& event ){
  theta_ctl->Enable(false);
  phi_ctl->Enable(false);

  h_ctl->Enable(true);
  k_ctl->Enable(true);
  l_ctl->Enable(true);
}

void sp_gui_imp::h_change( wxSpinEvent& event )
{
   
  //update parameter
  try {
    sp->set_h( h_ctl->GetValue() );
  } catch (invalid_parameter_exception e){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION);
    h_ctl->SetValue( sp->get_h() );
  }
  
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
    

}

void sp_gui_imp::k_change( wxSpinEvent& event )
{

  //update parameter
  try {
    sp->set_k( k_ctl->GetValue() );
  } catch (invalid_parameter_exception e){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION);
    k_ctl->SetValue( sp->get_k() );
  }
  
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
}

void sp_gui_imp::l_change( wxSpinEvent& event )
{
  //update parameter
  try {
    sp->set_l( l_ctl->GetValue() );
  } catch (invalid_parameter_exception e){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION);
    l_ctl->SetValue( sp->get_l() );
  }

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
}

void sp_gui_imp::slicewidth_change( wxSpinDoubleEvent& event )
{
  //update parameter
  try {
    sp->set_slice_width( slicewidth_ctl->GetValue() );
  } catch ( invalid_parameter_exception e ){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION );
    slicewidth_ctl->SetValue( sp->get_slice_width() );
  }

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
  try {
    sp->set_slice_height( sliceheight_ctl->GetValue() );
  } catch ( invalid_parameter_exception e ){
    wxMessageBox( e.details(), e.what(), wxICON_ERROR);
    sliceheight_ctl->SetValue( sp->get_slice_height() );
  }

  //redraw
  draw_preview();    
}

void sp_gui_imp::n_points_xy_change( wxSpinEvent& event ){
  //update parameter
  try {
    sp->set_n_points_x( n_points_xy_ctl->GetValue() );
    sp->set_n_points_y( n_points_xy_ctl->GetValue() );
  } catch ( invalid_parameter_exception e ){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION );
    n_points_xy_ctl->SetValue( sp->get_width() );
  }

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
  try {
    sp->set_n_points_z( n_points_z_ctl->GetValue() );
  } catch ( invalid_parameter_exception e ){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION );
    n_points_z_ctl->SetValue( sp->get_depth() );
  }

  //update_geometry
  sp->update_geometry();

  //update container size
  sp->update_containers();

  //update gui
  dz_ctl->SetValue( sp->get_dz() );

  //redraw
  draw_preview();
}


/** 
 * a membrane from the list view was selecetd. Load props into
 * textboxes 
 */
void sp_gui_imp::membrane_selected( wxListEvent& event ){

  //get index;
  int i_idx = event.GetIndex();
    
  membrane_idx->SetValue( std::to_string(i_idx) );
  membrane_width->SetValue( membranes_ctrl->GetItemText( i_idx, 2 ) );
  membrane_dist->SetValue( membranes_ctrl->GetItemText( i_idx, 1 ) );
  
}


void sp_gui_imp::membrane_save(){


  double d,w;
  long i;
  i = membrane_idx->GetValue();
  d = membrane_dist->GetValue();
  w = membrane_width->GetValue();
  
  //save new data to sp object
  try {
    sp->edit_membrane( i, d, w );
  } catch ( invalid_parameter_exception e ){
    wxMessageBox( e.details(),
		  e.what(),
		  wxICON_ERROR);
  }

  //update values
  read_membranes();
  
  //redraw
  draw_preview();
}


void sp_gui_imp::mem_change_w( wxSpinDoubleEvent& event ){
  membrane_save();
}

void sp_gui_imp::mem_change_d( wxSpinDoubleEvent& event ){
   membrane_save();
}

void sp_gui_imp::membrane_add( wxCommandEvent& event ){

  double d,w;
  long i;
  i = membrane_idx->GetValue();
  d = membrane_dist->GetValue();
  w = membrane_width->GetValue();

  sp->add_membrane( d, w );
  read_membranes();

  //draw preview
  draw_preview();
}


void sp_gui_imp::membrane_delete( wxCommandEvent& event ){

  // delete membrane
  try {
    sp->delete_membrane( membrane_idx->GetValue() );
  } catch ( invalid_parameter_exception e ){
    wxMessageBox( e.details(),
		  e.what(),
		  wxICON_INFORMATION);
  }

  //update gui
  read_membranes();
  draw_preview();
}


void sp_gui_imp::level_change(  wxSpinDoubleEvent& event ){

  try {
    sp->set_channel_vol_prop( level_set_ctrl->GetValue() );
  } catch ( invalid_parameter_exception e ){
    wxMessageBox( e.details(),
		  e.what(),
		  wxICON_INFORMATION);
  }

  //read possibly (exception handling) new values
  level_set_ctrl->SetValue( sp->get_channel_prop() );

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


  sp->print_grid("grid.dat");

  /*
  
  distance_transform<int> dt ( sp->get_grid(), sp->get_width(), sp->get_height(), sp->get_depth(), sp->get_dx(), sp->get_dy(), sp->get_dz() );  
  dt.compute_distance_map();
  std::vector<double> dmap = dt.get_distance_map();  
  

  homotopic_thinning compute_network ( sp->get_width(), sp->get_height(), sp->get_depth(), sp->get_channel(), dmap );

  auto network = compute_network.find_channel_skeleton( 1 );

  auto coords = sp->get_points();
  
  std::ofstream out ("network_1.dat");
  for( auto it : network ){
    out << coords[3*it] << " " << coords[3*it+1] << " " << coords[3*it+2] << std::endl;
  }
  out.close();

  network = compute_network.find_channel_skeleton( 2 );
  
  out.open("network_2.dat");
  for( auto it : network ){
    out << coords[3*it] << " " << coords[3*it+1] << " " << coords[3*it+2] << std::endl;
  }
  out.close();  


  network = compute_network.find_channel_skeleton( 3 );
  
  out.open("network_3.dat");
  for( auto it : network ){
    out << coords[3*it] << " " << coords[3*it+1] << " " << coords[3*it+2] << std::endl;
  }
  out.close();  


  auto sur = sp->get_surface_points( 1 );
  out.open("interface_1.dat");
  for( auto it : sur ){
    out << coords[3*it] << " " << coords[3*it+1] << " " << coords[3*it+2] << std::endl;
  }
  out.close();

  sur = sp->get_surface_points( 2 );
  out.open("interface_2.dat");
  for( auto it : sur ){
    out << coords[3*it] << " " << coords[3*it+1] << " " << coords[3*it+2] << std::endl;
  }
  out.close();  

  sur = sp->get_surface_points( 3 );
  out.open("interface_3.dat");
  for( auto it : sur ){
    out << coords[3*it] << " " << coords[3*it+1] << " " << coords[3*it+2] << std::endl;
  }
  out.close();  
  
  
  //std::string fn (m_filePicker2->GetPath().c_str());

  //sp->find_channel_graph( 1, true );
  //sp->find_channel_skeleton( 1 );
  
    //sp->print_grid( fn );
  
  /*
  if( fn == "" ){
    wxMessageBox( wxT("Please select a valid filename"),
		  wxT("Nothing to see here"),
		  wxICON_INFORMATION);
  } else {

    //chck for image type to be saved to
    std::string extension = fn.substr( fn.size() - 4, 4 );
    if( extension[0] == '.' ){
      extension = extension.substr(1, 3 );
    }

    if( extension == "png" || extension == "tiff" ){
      img->SaveFile( fn );
    } else {
      wxMessageBox( wxT("Image type not supported"),
		    wxT("Error"),
		    wxICON_ERROR);
    }
  }
  */  
}

void sp_gui_imp::button_quit( wxCommandEvent& event )
{
  Close( true );
}


/**
 * This function calls the neccessary function from the underlying
 * \ref surface_projection class to update the slice orientation
 * according to the Miller indeces and update the GUI
 */
void sp_gui_imp::update_orientation_from_hkl(){

  //set parameters
  try {
    sp->set_h( h_ctl->GetValue() ); sp->set_k( k_ctl->GetValue() ); sp->set_l( l_ctl->GetValue() );
  } catch (invalid_parameter_exception e){
    wxMessageBox( e.details(), e.what(), wxICON_INFORMATION);
    h_ctl->SetValue( sp->get_h() );k_ctl->SetValue( sp->get_k() );l_ctl->SetValue( sp->get_l() );
  }

  //convert indeces to angles
  sp->set_orientation_from_hkl();
  
  //periodicity
  sp->update_periodicity_length();
  update_controls_periodicity();
    
  //update gui
  theta_ctl->SetValue( sp->get_theta() );
  phi_ctl->SetValue( sp->get_phi() );
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


/*
 * This method updates the status and progress bar while the
 * computation of the surface projection is performed. It extis when a
 * progress of 1 is reached
 */
void sp_gui_imp::update_status_bar(){

  std::stringstream ps ("");
  
  while( true ){
    std::string tmp_string ( status_string );
    status->SetStatusText( tmp_string, 0 );
    ps.str("");
    ps << std::setprecision(4) << progress * 100 << "%";
    status->SetStatusText( ps.str(), 1 );
    if( progress == 1 && tmp_string == "Ready"){
      return;
    }
  }

  
}

/**
 * This method calls the neccessary functions from \ref surface
 * projection to compute the projection with the current parameters
 * and draw it onto the GUI
 */
void sp_gui_imp::compute_and_draw(){

  
  
  //compute the projection in a separate thread
  std::thread computation ( &surface_projection::compute_projection, sp );
  //keep the interface up to date in a separate thread
  std::thread progress_status ( &sp_gui_imp::update_status_bar, this );

  //wait for the threads to join
  computation.join();
  progress_status.join();

  //m_progressbar->SetValue(0);
  
  /*
  sp->compute_volume();

  try {
    sp->compute_surface_area();
  } catch(std::out_of_range e){
    std::cout << e.what() << std::endl;
  }

  
  std::cout << "volumes: ";
  for( auto it : sp->get_channel_volumes() ){
    std::cout << it << " ";
  }
  std::cout << std::endl;

  std::cout << "surface_areas: ";
  for( auto it : sp->get_membrane_surface_area() ){
    std::cout << it << " ";
  }
  std::cout << std::endl;
  
  */

  unsigned char *proj = sp->get_image( invert_ctl->GetValue() );

  //width and height of image
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


/**
 * Redraws (without recomputing) the current image and especially
 * rescales the current image to the current size of the
 * wxStaticBitmap
 */
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
