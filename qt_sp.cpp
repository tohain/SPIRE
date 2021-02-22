/* Projection tool - compute planar projections of triply periodic
 * minimal surfaces 
 * Copyright (C) 2021 Tobias Hain
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see https://www.gnu.org/licenses.
 */


#include "qt_sp.hpp"


sp_qt::sp_qt() : surface_projection(){

  connect( this, SIGNAL( geometry_changed() ), this, SLOT( update_geometry_() ) );
  
}

sp_qt::~sp_qt(){

}



void sp_qt::update_geometry_(){

  try {
    update_geometry();
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
  
}

void sp_qt::change_surface_type( int ind ){

  //backup channel proportion
  double vol_prop = get_channel_vol_prop();

  set_type( ind );

  set_channel_vol_prop( vol_prop );

  emit geometry_changed();
  emit parameter_changed();
}


void sp_qt::change_z_points( int val ){
  try {
    set_n_points_z( val );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_xy_points( int val ){
  try {
    set_n_points_x( val );
    set_n_points_y_to_unitcell( );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
  emit geometry_changed();
  emit parameter_changed();
}



void sp_qt::change_uc_scale_ab( double ab ){
  
  set_uc_scale_ab( ab );

  // if we have membranes of cubic symmetries, change c scaling as
  // well, since it needs to be cubic
  if( surface_choices[type].substr(0,8) != "Wurtzite" ){
    set_uc_scale_c( ab );
  }

  update_a();  
  
  emit geometry_changed();
  emit parameter_changed();
}


void sp_qt::change_uc_scale_c( double c ){

  //only allow if we have other than cubic symmetries

  if( surface_choices[type].substr(0,8) == "Wurtzite" ){
  
    set_uc_scale_c( c );
    update_a();
  
    emit geometry_changed();
    emit parameter_changed();
  }
}



void sp_qt::change_vol_prop( double val ){
  set_channel_vol_prop( val );
  emit parameter_changed();
}

void sp_qt::change_lvl_set( double val ){
  set_surface_level( val );
  emit parameter_changed();
}



void sp_qt::change_hkl( int _h, int _k, int _l ){

  if( !(_h == 0 && _k == 0 && _l == 0 ) ){
    set_h( _h ); set_k( _k ); set_l( _l );
    emit geometry_changed();
    emit parameter_changed();
  } else {
    emit send_message( "At least one parameter of hkl must be non-zero", 1 );
    emit parameter_changed();
  } 
}


void sp_qt::change_slice_thickness( double val ){
  try {
    set_slice_thickness( val );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_slice_width( double val ){
  try {
    set_slice_width( val );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_slice_height( double val ){
  try {
    set_slice_height( val );
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }    
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_slice_position( double val ){
  set_slice_position( val );
  emit geometry_changed();
  emit parameter_changed();
}

void sp_qt::change_membranes( std::vector<double> val ){
  set_membranes( val );
}

void sp_qt::compute_projection(){

  update_containers(); 
  
  //get the points in the slice
  emit set_status( 0, 1 );
  set_up_points();
  
  //reset grid
  memset( grid.data(), 0, sizeof(short) * grid.size() );
  
  //get grid
  set_grid();

  // apply channel colors
  for( unsigned int ii=0; ii < channel_filled.size(); ii++){

    if( ii % 2 == 0 ){
      //channel, if it is zero, do nothing, otherwise mark it
      if( channel_filled[ii] != 0 ){
	set_channel_color( ii + 1, 1 );
      }
    } else {
      //membrane, if it is one, do nothing, otherwise mark it
      if( channel_filled[ii] != 1 ){
	set_channel_color( ii + 1, 0 );
      }      
    }

  }
  
  //get projection
  memset( projection.data(), 0, sizeof(float) * projection.size() );
  project_grid();
  
  emit projection_changed();
  emit set_status( 0, 0 );
}


void sp_qt::do_something(){

  std::cout << "doing something" << std::endl;

  dt.set_parameters( grid, std::vector<unsigned int> {n_points_x, n_points_y, n_points_z},
		     std::vector<double> {dx, dy, dz}, true);  
  dt.compute_distance_map();
  dt.compute_max_radius_covering();
  

  //print_channel_width_distribution( "channel_width_dist.dat" );
  print_max_rad_transform_dist( "mrct.dat" );

  std::vector<float> mrct = dt.get_max_radius_covering();
  std::vector<float> dmap = dt.get_distance_map();  
  
  // print out slices
  double min = 0;
  double mrct_max = *std::max_element( mrct.begin(), mrct.end() );
  double dmap_max = *std::max_element( dmap.begin(), dmap.end() );  

  unsigned char *img_data = new unsigned char[ n_points_x * n_points_y ];
  unsigned char *dmap_data = new unsigned char[ n_points_x * n_points_y ];  
  
  for( int zz = 0; zz < n_points_z; zz++){
    std::ofstream data_out ("mrct_" + std::to_string(zz) + ".dat" );
    std::ofstream dmap_out ("dmap_" + std::to_string(zz) + ".dat" );
    std::ofstream grid_out ("grid_" + std::to_string(zz) + ".dat" );
    for( int yy = 0; yy < n_points_y; yy++){      
      for( int xx = 0; xx < n_points_x; xx++ ){

	int ind = zz + xx * n_points_z + yy*n_points_z*n_points_x;

	double mrct_val = (mrct[ind] / mrct_max) * 255.0;
	unsigned char val_c_mrct = static_cast<unsigned char> (mrct_val);

	double dmap_val = (dmap[ind] / dmap_max) * 255.0;
	unsigned char val_c_dmap = static_cast<unsigned char> (dmap_val);
	
	int img_ind = xx + yy*n_points_x;
	
	img_data[img_ind] = val_c_mrct;
	dmap_data[img_ind] = val_c_dmap;

	data_out << std::setprecision(3) << std::setw(6) << mrct[ind] << " ";
	dmap_out << std::setprecision(3) << std::setw(6) << sqrt(dmap[ind]) << " ";
	grid_out << std::setw(1) << std::setprecision(1) << grid[ind] << " ";
      }
      data_out << std::endl;
      dmap_out << std::endl;
      grid_out << std::endl;
    }

    data_out.close();
    dmap_out.close();
    
    write_image( "layer_" + std::to_string(zz) + ".ppm", img_data, n_points_x, n_points_y );
    write_image( "dmap_" + std::to_string(zz) + ".ppm", dmap_data, n_points_x, n_points_y );    
    
    
  }

  std::cout << "done" << std::endl;
  
}



void sp_qt::copy_parameters( sp_qt *source ){

  // orientation
  set_h( source->get_h() );
  set_k( source->get_k() );
  set_l( source->get_l() );  

  // slice
  set_slice_thickness( source->get_slice_thickness() );
  set_slice_height( source->get_slice_height() );
  set_slice_width( source->get_slice_width() );
  set_slice_position( source->get_slice_position() );

  // points
  set_n_points_x( source->get_width() );
  set_n_points_y( source->get_height() );
  set_n_points_z( source->get_depth() );
  
  // uc scale
  set_uc_scale_ab( source->get_uc_scale_ab() );
  set_uc_scale_c( source->get_uc_scale_c() );
  
  // structure
  set_surface_level( source->get_surface_level() );
  set_type( source->get_type() );
  set_membranes( source->get_membranes() );
  set_channel_fill( source->get_channel_fill() );
    
  try {
    update_geometry();
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }

  
  emit geometry_changed();
  emit parameter_changed();
}



/**
 * Saves a png of the current projection with a legend of all provided
 * parameters. If empty, all parameters are added to the image
 *
 * the parameters are added in a margin below the image. The size of
 * the text is tried to kept constant at a height of
 * pixsize_fractional of the image height, but minium of
 * 25pixels. If the image resolution is too low, the image is
 * upscaled so it can accomodate the text in the desired resolution.
 *
 * The text is arranged on a grid, with a cell length of the longest
 * word size. If a shorter word fits the line, it is squeezed in,
 * however, only at the start of a grid cell
 *
 * reallocates the memory of image, therefore returns a pointer to the
 * reallocated memory
 *
 * renders all parameters, if none are given
 */
unsigned char* sp_qt::save_png_legend( unsigned char* img,
			     unsigned int img_x,
			     unsigned int img_y,
			     std::string textcolor,
			     std::string backcolor,
			     std::string fn,
			     std::vector<std::string> parameters ){

  
  if( parameters.size() == 0 ){
    parameters = std::vector<std::string> ( surface_projection::parameter_names );
  }
  
  
  /*
   * creating the legend
   */ 
  
  // constants determining the layout
  double pixsize_fractional = 0.045; // size of the font (height) in
				     // fractions of the height of the
				     // image
  double margin = 1.2; // one line will have the height of
		       // margin*fontsize(in pixels)
  unsigned int min_pixsize = 25;  // the minimum font height(in
				  // pixels)

  // would be pixel size if image is not scaled
  unsigned int pixsize = (pixsize_fractional * img_y);  

  // if pix size gets too small, it is not rendered readable, so we
  // need to upscale the image
  double scale = 1.0;
  if( pixsize < min_pixsize ){
    scale = (float) min_pixsize / pixsize;
  }
  
  
  // setting the pixel size; this makes the point size not used
  QFont font = QFont();
  font.setCapitalization( QFont::AllLowercase ); // not sure if we want this;  
  
  // set the current font size, so the longest word length can be
  // measured in the non-scaled image
  font.setPixelSize( pixsize );
  
  // get longest string in pixels to compute grid sizing (cell length)
  QFontMetrics fm ( font );  

  // store the text
  std::vector<std::string> parameter_words;


  // the length of the longest word in unscaled pixels
  unsigned int max_width = 0;

  // create all the parameter words
  for( unsigned int ii=0; ii<parameters.size(); ii++){
    std::stringstream ss;
    
    // get the position in the parameter_names array of the ii-th
    // parameter in the batch loop
    size_t ind = std::distance( surface_projection::parameter_names.begin(),
				std::find( surface_projection::parameter_names.begin(),
					   surface_projection::parameter_names.end(),
					   parameters[ii] ) );
    
    // so now surface_projection::parameter_names[ind] == bc.ops.parameters[ii]

    // construct the string containing name and parameter
    ss << surface_projection::parameter_names_short[ind] << "=";

    // special treatment for surface type so string is printed instead of index
    if( parameters[ii] == surface_projection::parameter_names[0] ){
      ss << surface_choices[ int( get_parameter( parameters[ii] )  )  ] << std::endl;
    } else {      
      ss << std::setprecision(3) << get_parameter( parameters[ii] ) << std::endl;
    }
    
    // store words
    parameter_words.push_back( ss.str() );
    
    // compute length
    unsigned int word_width = fm.horizontalAdvance( ss.str().c_str() );      
    if( word_width > max_width ){
      max_width = word_width;
    }      
  }
  

  //
  // compute the number of columns and rows we need, not really
  // needed for the grid but we need to estimate how many pixels we
  // need to add for the margin text
  //
  unsigned int nr_cols = int(img_x / max_width);
  if( nr_cols == 0 ) nr_cols++;
  unsigned int nr_rows = int(parameters.size() / nr_cols)+1;
  
  // now add a few rows as writing area and realloc the memory
  unsigned int new_img_y = img_y + int( margin * nr_rows * pixsize);    


  /*
   * This is very risky (but efficient!) since img originally was
   * allocated with new, but is reallocated with realloc. Seems to run
   * fine with gcc 10 though, but keep that in mind!
   */
  img = (unsigned char*) realloc (img, sizeof(unsigned char) * img_x * new_img_y);
  
  // make the background white/black
  memset( img + sizeof(unsigned char) * img_x * img_y,
	  0,
	  sizeof(unsigned char)*img_x*(new_img_y-img_y) );
  

  // create the qimage object where we can write on
  QImage qimg = QImage( img, img_x, new_img_y, img_x, QImage::Format_Grayscale8 );
  // scale it, if image is large enough hopefully this does nothing
  qimg = qimg.scaledToHeight( int(scale * new_img_y) );

  // set real font size after image scaling
  font.setPixelSize( int(scale * pixsize) );    
  
  // font color
  QPen pen = QPen();
  pen.setColor( textcolor.c_str() );

  // setup painter
  QPainter poet ( &qimg );
  poet.setFont( font );
  poet.setPen( pen );
  
  // draw the background
  poet.fillRect( 0, scale*img_y, qimg.width(), qimg.height()-scale*img_y, backcolor.c_str() );
  
  
  // update the new font size
  fm = QFontMetrics( font );
  
  
  // put the words below the text. Only start at grid locations, but
  // if a cell is cutoff by the image, check if we can squeeze a
  // shorter word in
  int ind=0, x=0, y = scale * img_y;
  int word_width_pix;
  int c=0; //current col
  while( ind < parameter_words.size() ){
    
    // get word length in scale pixels
    word_width_pix = fm.horizontalAdvance( parameter_words[ind].c_str() );
    
    // check if word matches this line
    if( x + word_width_pix > qimg.width() ){
      // line does not fit that word, advance to next one
      x = 0;
      y += scale * margin * pixsize;
      c=0;
    }
    
    // draw it
    poet.drawText( x, y, word_width_pix, scale * margin * pixsize, Qt::AlignCenter,
		   parameter_words[ind].c_str() );
    // increase the col
    c++;
    // set x to the next grid start
    x = c * scale * max_width;
    
    ind++;
  }

  qimg.save( fn.c_str() );

  return img;
}




void sp_qt::update_measurements( QString what ){  

  // set "measurement" (first digit) to "busy" (second digits)
  emit set_status( 1, 1 );

  try {
    update_geometry();
  } catch ( invalid_parameter_exception e ){
    emit send_message( e.what() );
  }  

  emit send_message( "Computing projection", 1 );
  compute_projection();

  // the distance map of the current object is not up to date, since
  // the "fill channel" option is perforemd after the last change of
  // the grid, so keep that in mind!


  
  emit send_message( "Computing " + what, 1 );  
  
  if( what == "Volumes" ){
    compute_volume();
  }

  if( what == "Areas" ){
    compute_surface_area();
  }

  if( what == "Networks" ){
    compute_channel_network();
  }

  // assuming a periodic structure
  if( what == "Percolation_P" ){
    compute_percolation_threshold( true );
  }

  // assuming a aperiodic structure
  if( what == "Percolation_A" ){
    compute_percolation_threshold( false );
  }  

  emit send_message("Measurement Done");
  emit set_status( 1, 0 );

  emit measurements_updated( QString( what ) );
}


void sp_qt::save_grid( QString fn ){

  print_grid( fn.toStdString() );

  emit send_message( "Saved grid!", 1);
}

void sp_qt::save_topological_network( int id, QString fn ){

  print_topological_network( id, fn.toStdString() );

  emit send_message( "Saved network!", 1);
}

void sp_qt::save_surface_points( int id, QString fn ){

  print_channel_surface_points( id, fn.toStdString() );

  emit send_message( "Saved membrane!", 1);
}

void sp_qt::set_slice_dim_to_uc(){

  set_slice_to_uc( 0.1 );

  emit geometry_changed();
  emit parameter_changed();

}


void sp_qt::change_channel_color( int id, int val ){
  set_channel_color( id, val );
}
