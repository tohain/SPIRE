#include "surface_projection.hpp"


invalid_parameter_exception::invalid_parameter_exception( std::string _msg ) : msg(_msg){
}

const char* invalid_parameter_exception::what() const throw() {
    return "Invalid parameter choice!";
}

const std::string invalid_parameter_exception::details() const {
  return msg;
}




/** This function prints the grid to a file. The file has the form 
 * x y z color
 * This function does NOT compute the grid, so make sure it is
 * up to date before calling this function
 */
void surface_projection::print_grid( std::string fn ){

  std::ofstream out ( fn );
  out << "#x y z color channel" << std::endl;
  for( unsigned int ii=0; ii<points.size(); ii+=3 ){
    out << points[ii] << " " << points[ii+1] << " ";
    out << points[ii+2] << " " << grid[ii/3] << " ";
    if( channel[ii/3.] < 0 ){
      out << -channel[ii/3.];
    } else {
      out << channel[ii/3.];
    }
    out << std::endl;
  }

  out.close();
}


/** Standard constructor initialize with the standard values and
 *  derive some more quantities
 */
surface_projection::surface_projection( double *p, char* stat) : ntucs(1), slice_width(0.1), slice_height(0.5), a(1), inv_a(2*M_PI), n_points_x(90), n_points_y(90), n_points_z(50), type(2), h(0), k(0), l(1), surface_level( 0.0f ), progress(p), status( stat ) {

  //update geometry
  //sets dx,dy,dz, L
  update_geometry();

  //resize containers
  update_containers();

  //update orientation from hkl
  //sets theta,phi
  set_orientation_from_hkl();
  
  //periodicity
  //sets periodicty length
  update_periodicity_length();

  
  //add main membrane. This can't be deleted
  membranes.push_back( 0 );
  membranes.push_back( 0.03 );

}


surface_projection::~surface_projection(){

}


double surface_projection::dot_prod ( std::vector<double> v, std::vector<double> w ){
  return v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
}



std::vector<double> surface_projection::dot_prod( Matrix m, std::vector<double> v ){
  std::vector<double> r (3, 0);
  r[0] = dot_prod( m.v, v);
  r[1] = dot_prod( m.w, v);
  r[2] = dot_prod( m.z, v);
  return r;
}


/**
 * Returns a rotation matrix. The rotation is by ang against the
 * clock around the x-Axis
 */
Matrix surface_projection::get_x_rot_m (double ang) const {
  Matrix R;
  R.v = { 1, 0,           0          };
  R.w = { 0, cos(ang), -sin(ang) };
  R.z = { 0, sin(ang),  cos(ang) };
  return R;
}

/**
 * Returns a rotation matrix. The rotation is by ang against the
 * clock around the z-Axis
 */
Matrix surface_projection::get_z_rot_m (double ang) const {
  Matrix R;
  R.v = { cos(ang), -sin(ang), 0 };
  R.w = { sin(ang),  cos(ang), 0 };
  R.z = { 0,         0,        1 };  
  return R;
}

 
/**
 * Nodal representation of a Gyroid Ia\bar(3)d surface. From
 *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
 * (1991) 83: 407. https://doi.org/10.1007/BF01313411
 */
double surface_projection::level_set_gyroid( double x, double y, double z, double inv_a) {
  return sin(inv_a*x)*cos(inv_a*y) + sin(inv_a*y)*cos(inv_a*z) + cos(inv_a*x)*sin(inv_a*z);
}

/**
 * Nodal representation of a Diamond Pn\bar(3)m surface. From
 *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
 * (1991) 83: 407. https://doi.org/10.1007/BF01313411
 */
double surface_projection::level_set_diamond( double x, double y, double z, double inv_a) {
  return cos(inv_a*x)*cos(inv_a*y)*cos(inv_a*z) - sin(inv_a*x)*sin(inv_a*y)*sin(inv_a*z);
}
 
/**
 * Nodal representation of a Primitive Im\bar(3)m surface. From
 *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
 * (1991) 83: 407. https://doi.org/10.1007/BF01313411
 */
double surface_projection::level_set_primitive( double x, double y, double z, double inv_a) {
  return cos(inv_a*x)+cos(inv_a*y)+cos(inv_a*z);
}

/**
 * For debugging purposes: just a simple layer
 *
 */
double surface_projection::level_set_layer( double x, double y, double z, double inv_a) {
  return mod(z, 1./a);
}

/**
 * For debugging purposes: just a simple sphere
 *
 */
double surface_projection::level_set_sphere( double x, double y, double z, double inv_a) {
  return (x*x+y*y+z*z);
}

/**
 * Computes and sets the theta and phi angles to match the orientation
 * given by the Miller indeces
 */
void surface_projection::set_orientation_from_hkl(){

  // Normal vector of the plane given by Miller indeces
  std::vector<double> n = {double(h), double(k), double(l)};

  //make it unit length
  double norm_scale = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
  n[0]/=norm_scale;n[1]/=norm_scale;n[2]/=norm_scale;
 
  //get angle for rotation in xy plane
  double _phi;
  if( k == 0 ){
    if( h < 0)
        _phi = 3*M_PI/2.;
    else
        _phi = M_PI/2.;
  } else {
    if( h > 0){
        if (k > 0)
            _phi = atan2( n[0], n[1] ) ;
        else 
            _phi = 2*M_PI - atan2( n[0], -n[1] ) ;
    } else{
        if (k > 0)
            _phi = M_PI - atan2( -n[0], n[1] ) ;
        else 
            _phi = M_PI + atan2( -n[0], -n[1] ) ;
    }
  }


  //get length of projection of normal vector in xy plane
  double n_xy = sqrt(n[0]*n[0]+n[1]*n[1]);
   
  //get angle to rotate around x axis to align vector with {0, 0, 1}  
  double _theta;
  if( l == 0 ){
    _theta = M_PI/2.;
  } else {
      if(l > 0)
        _theta = atan2( n_xy, n[2] );
      else
        _theta = M_PI - atan2( n_xy, -n[2] );
  }

  //set the angles
  theta = _theta; phi = _phi;
}



void surface_projection::set_up_points(){
  
  //vectors of slice
  std::vector<double> nx = {1, 0, 0};
  std::vector<double> ny = {0, 1, 0};
  std::vector<double> nz = {0, 0, 1};  

  //rotation matrices
  //substract the angle from 2pi since we're rotating in mathematical
  //negative orientation (with the clock
  Matrix Rx = get_x_rot_m(2*M_PI - theta), Rz = get_z_rot_m(2*M_PI - phi);

  //rotate
  nx = dot_prod( Rz, dot_prod( Rx, nx));
  ny = dot_prod( Rz, dot_prod( Rx, ny));
  nz = dot_prod( Rz, dot_prod( Rx, nz));  
  
  //the total index of that particle in the array
  int ind = 0;
  double iy = -L_2;
  //a single loop would be nicer, but i find this less confusing  
  for(unsigned int ii=0; ii<n_points_y; ii++){ //height, the row index (vertical)
    double jx = -L_2;
    for(unsigned int jj=0; jj<n_points_x; jj++){ //width, the column index! (horizontal)

      double kz;

      //check if we are periodic. If so, slice_height will be handeled differently
      if( periodicity_length == -1 ){
	
	kz = slice_height - 0.5*slice_width;// kz = ((kk*dz)-L_2.) + (slice_height+L_2.) - 0.5*slice_width;
	//aperiodic. slice_height is an absolute length
      } else {
	
        kz = periodicity_length*(slice_height - 0.5*slice_width);// kz = ((kk*dz)-L_2.) + ((periodicity_length*slice_height)+L_2.) - 0.5*periodicity_length*slice_width;
	//periodic. slice_height is the fraction of the periodicity length
      }
      
      for(unsigned int kk=0; kk<n_points_z; kk++){ //depth
	
	
	//now compute the absolute position of all voxels
	  double x = jx*nx[0] + iy*ny[0] +kz*nz[0];
	  
	  double y = jx*nx[1] + iy*ny[1] +kz*nz[1];
	  
	  double z = jx*nx[2] + iy*ny[2] +kz*nz[2];
	  
	  //assign the position
	  points[ind]=x;
	  points[ind+1]=y;
	  points[ind+2]=z;
	  
	  kz += dz; // next z-pixel
	  ind += 3; // running index
      }
      jx += dx; // next x-pixel
    }
    iy += dy; // next y-pixel
  }
  
}


/**
 * Evalutes the voxels in the level set. Sets the "colors" of the
 * voxels to 1 if they are within the membrane. Also updates the
 * channel array
 */
void surface_projection::set_grid(){

  // First compute the distance map of the current grid

  //We need to set an interval where points are getting markes. Since
  //the surface is samples it's highly unlikely that a point has
  //exactly the wanted level set constraint
  double min_width = 0.025;
  
  // the level set value of the present points
  double level = 0; 
  
  for( unsigned int ii=0; ii<points.size(); ii+=3 ){

    //compute level set
    
    if( type == 0 ){ //gyroid
      level = level_set_gyroid( points[ii], points[ii+1], points[ii+2], inv_a );
    } else if( type == 1 ){ //diamon
      level = level_set_diamond( points[ii], points[ii+1], points[ii+2], inv_a );
    } else if( type == 2 ){ //primitive
      level = level_set_primitive( points[ii], points[ii+1], points[ii+2], inv_a );
    } else if( type == 3 ){ //layer
      level = level_set_layer( points[ii], points[ii+1], points[ii+2], inv_a );
    } else if( type == 4 ){ //sphere
      level = level_set_sphere( points[ii], points[ii+1], points[ii+2], inv_a );
    }
    else {
      throw std::string("type not supported");
    }
    
    //mark appropriate voxels '1'. Also check if we are in the membrane, outside or inside
    if( level <= (surface_level + min_width) && level > ( surface_level - min_width ) ){
      grid[int(ii/3)] = 1;
    }



    if( level > surface_level ){
      //outside
      channel[int(ii/3.)] = -1;
    } else if ( level <= surface_level ) {
      //inside
      channel[int(ii/3.)] = 1;
    } 
	
  }


  // next step is to compute the distance transform. We need to label
  // all points according to in which channel they are also. We do
  // that by comparing the distances to the main membrane with the
  // distance map
  
  // compute the distance transform of the grid
  //distance_transform<int> dt ( grid, n_points_x, n_points_y, n_points_z, dx, dy, dz );
  distance_transform<int> dt ( grid, n_points_x, n_points_y, n_points_z, dx, dy, dz);  
  dt.compute_distance_map();
  std::vector<double> dmap = dt.get_distance_map();  


  // now assign channel number
  
  // get all the membranes
  std::vector<double> mem_pos ( membranes.size(), 0);
  for(unsigned int ii=0; ii<membranes.size(); ii+=2){
    mem_pos[ii] = membranes[ii] + membranes[ii+1]/2.;
    mem_pos[ii+1] = membranes[ii] - membranes[ii+1]/2.;
  }
  mem_pos.push_back( -std::numeric_limits<double>::max() );
  mem_pos.push_back( std::numeric_limits<double>::max() );  
  //sort
  std::sort( mem_pos.begin(), mem_pos.end() );

  
  for( unsigned int ii=0; ii<channel.size(); ii++){
    //start with the inner most shell
    int ch_id=0;

    //find the right channel
    for( ; ch_id < mem_pos.size()-1; ch_id++ ){     
      if( channel[ii] * sqrt(dmap[ii]) >= mem_pos[ch_id] &&
	  channel[ii] * sqrt(dmap[ii]) < mem_pos[ch_id+1])
	  break;
    }    
    
    //found the appropriate membrane, save the channel number
    //the sign still marks, if in or outside of main membrane
    channel[ii] = (ch_id+1)*channel[ii];
  }

  // now also mark each point, which is within the desired thickness of the membrane
  for( unsigned int ii=0; ii<grid.size(); ii++){

    bool mark_point = false;
    double real_distance = sqrt( dmap[ii] );
    
    //check if this point is within the "main" membrane
    // this will extend symmetrically into both channels
    if( real_distance < 0.5*membranes[1] ){
      mark_point = true;
    }

    //check if this point is within several other membranes
    for( unsigned int jj=2; jj<membranes.size(); jj+=2){

      if( membranes[jj] >= 0 ){
	//outside membrane
	
	if( real_distance > membranes[jj] - membranes[jj+1]/2. &&
	    real_distance < membranes[jj] + membranes[jj+1]/2. &&
	    channel[ii] > 0 ){
	  
	  mark_point = true;
	}
      } else {
	//inside

	if( real_distance > -membranes[jj] - membranes[jj+1]/2. &&
	    real_distance < -membranes[jj] + membranes[jj+1]/2. &&
	    channel[ii] < 0 ){
	  
	  mark_point = true;
	}	
      }

    }//end cycle membranes
    
    if( mark_point ){
      grid[ii] = 1;
    }
  }
  
}


void surface_projection::project_grid (){
  
  
  unsigned int ind_grid = 0;
  unsigned int ind_proj = 0;
  for(unsigned int ii=0; ii<n_points_y; ii++){ //y direcion=height (vertical)
    for(unsigned int jj=0; jj<n_points_x; jj++){//x direction=width (horizontal)
      for(unsigned int kk=0; kk<n_points_z; kk++){
        

	if(grid[ind_grid] == 1){
	  projection[ind_proj] += 1;
	}

    ind_grid++; // run 3d-index
	
      }
    ind_proj++;// run 2d-index
    }
  }

}

/**
 * This function resizes the containers for the voxel coordinates,
 * gird and projection after a geometry change
 */
void surface_projection::update_containers(){
  //resize and reset arrays
  points.resize( n_points_x * n_points_y * n_points_z * 3, 0 );
  //memset( points.data(), 0, sizeof(double) * points.size() );
  
  /// Array holding the color (electron density) of the voxels
  grid.resize( n_points_x * n_points_y * n_points_z, 0 );
  //memset( grid.data(), 0, sizeof(int) * grid.size() );

  /// The channel
  channel.resize( n_points_x * n_points_y * n_points_z, 0 );
  
  /// The 2D projection
  projection.resize( n_points_x * n_points_y, 0 );
  memset( projection.data(), 0, sizeof(double) * projection.size());
}

/**
 * This function updates the geometry, i.e. recomputes all values
 * depending on the unit cell size and number of unit cells.
 */
void surface_projection::update_geometry(){
  //update box
  L = ntucs * a;

  L_2 = L/2.;

  //update points number
  dx = L / n_points_x;
  dy = L / n_points_y;
  dz = slice_width / n_points_z;  
}

/**
 * This function computes the projection. It is basically just a
 * wrapper around the functions performing the single steps. It might
 * be more efficient to do that in a single loop, but with modern
 * compiler optimizations and loop unrolling this should not make too
 * much of a difference
 */
void surface_projection::compute_projection( ){
    
  //get the points in the slice    
  status = "Setting up grid";
  set_up_points();
  
  //reset grid
  memset( grid.data(), 0, sizeof(int) * grid.size() );
  
  //get grid
  status = "Setting up the grid";
  set_grid();
  
  //get projection
  status = "Computing Projection";
  memset( projection.data(), 0, sizeof(double) * projection.size() );
  project_grid();
}


/** 
 * Computes the normal vector of the given orientation. In principle
 * this is equivalent to computing hkl, except we're not making even
 * numbers
 */
std::vector<double> surface_projection::get_normal() {

  std::vector<double> n = {0, 0, 1};  

  //rotation matrices
  //substract the angle from 2pi since we're rotating in mathematical
  //negative orientation (with the clock
  Matrix Rx = get_x_rot_m(2*M_PI - theta), Rz = get_z_rot_m(2*M_PI - phi);

  //rotate
  n = dot_prod( Rz, dot_prod( Rx, n));  

  double scale = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2] ); //necessary?
  n[0]/=scale;n[1]/=scale;n[2]/=scale;
  
  return n;
}

/**
 * Returns the modulo of two floating point numbers. That's nmerically
 * not a nice thing, since we have to take rounding errors into
 * account.
 */
inline double surface_projection::mod( double lhs, double rhs ){
    
  double divisor = lhs / rhs;

  //let's check if we're very close to the next integer
  double divisor_round  = round(divisor);

  //if we are, then let's round to that so we're exact again
  int divisor_int;
  if( std::fabs(divisor_round - divisor) < tolerance ){    
    divisor_int = int(divisor_round);
  } else {
    divisor_int = int(divisor);
  }
  
  double ret = lhs - (divisor_int * rhs);
  return ret;
}

/**
 * Computes the periodicity length. That means, how far the plane must
 * be moved perpendicular to its normal vector so it will reach
 * another chrystallographic equivalent position. This is implemented
 * by looking for a scaling factor for the normal vector so it will
 * point to a point periodically equivalent to the origin (where we
 * started)
 */
void surface_projection::update_periodicity_length(){

  //get the normal of the current orientation
  std::vector<double> n = get_normal();

  // the step size is chosen, so that one step will take at least one
  // direction from one periodic box to the next. Everything in
  // between can't be periodic anyways
  double step_size;
  if( std::fabs( n[0] ) > tolerance ){
    step_size = a / n[0];
  } else if ( std::fabs( n[1] ) > tolerance ){
    step_size = a / n[1];
  } else {
    step_size = a / n[2];
  }

  if (step_size < 0) step_size  = -step_size ;

  //stop after some steps, rather arbitrary number
  long long max_steps = 100;

  //count how many steps we need
  long long step_count = 1;


  bool more_steps = true;

  //increase steps, scale normal vector and see if we are at a
  //periodic origin again
  while( step_count < max_steps && more_steps){

    double xx = step_count * step_size * n[0];
    double yy = step_count * step_size * n[1];
    double zz = step_count * step_size * n[2];

    double ddxx = mod(xx, a);
    double ddyy = mod(yy, a);
    double ddzz = mod(zz, a);
    
    if( std::fabs( ddxx ) < tolerance &&
	std::fabs( ddyy ) < tolerance &&
	std::fabs( ddzz ) < tolerance  ){
      more_steps = false;
      break;
    }

    //not there yet, keep on goind
    step_count++;
  }

  
  if( step_count == max_steps ){
    //no succes
    periodicity_length = -1;
  } else {
    //success
    periodicity_length =step_count * step_size;
  }
}


/**
 * This function basically only converts the projection array to a
 * rescaled array, where the minimal value is set to 0 and the maximum
 * to 255. Also inverts if necessary. Memory management could be
 * imporved I guess
 */
unsigned char* surface_projection::get_image(bool invert){

  //new image array
  unsigned char* img = (unsigned char*) malloc( sizeof(unsigned char) * projection.size() );
  
  //find max value
  double max= *std::max_element( projection.begin(), projection.end() );
  double min= *std::min_element( projection.begin(), projection.end() );

  
  //write image data
  for(unsigned int ii=0; ii<projection.size(); ii++){
      //scale
      int u_scaled = static_cast<int>((projection[ii] - min)/(max - min)*255);
      //invert
      if(invert)
	u_scaled = 255 - u_scaled;
      //write in array
      img[ii] = u_scaled;
  }  

  return img;  
}

/** 
 * adds a membrane to the membrane array
 *
 *\param[in] dist The distance to the "original" membrane
 *\param[in] width The width of the membrane
 */
void surface_projection::add_membrane(double dist, double width ){

  membranes.push_back( dist );
  membranes.push_back( width );
}



/**
 * This function computes the volumes of the channels. It is just
 * adding up the volumes of the appropriate voxels
 *
 */
void surface_projection::compute_volume(){

  volumes = std::vector<double> (membranes.size() + 1, 0 );
  
  double vox_vol = dx*dy*dz;
  int ind;
  
  for( unsigned int ii=0; ii<channel.size(); ii++){
    
    ind = channel[ii];
    if( ind < 0 ) ind*=-1;
    volumes[ind-1]+=vox_vol;
  }

}



/**
 * returns the id of the pixel to the right
 */
int surface_projection::p_right( int val ){
  int m = val % (n_points_x*n_points_z);
  if( (n_points_x*n_points_z) - n_points_z <= m &&
      m < (n_points_x*n_points_z) ){
    return val - (( n_points_x - 1) * n_points_z);
  } else {
    return val + n_points_z;
  }
}

/**
 * returns the id of the pixel to the left
 */
int surface_projection::p_left( int val ){
  int m = val % (n_points_x*n_points_z);
  if( 0 <= m && m < n_points_z ){
    return val + ((n_points_x - 1) * n_points_z);
  } else {
    return val - n_points_z;
  }
}

/**
 * returns the id of the pixel above
 */
int surface_projection::p_up( int val ){
  if(val % n_points_z == n_points_z - 1){
    return val-(n_points_z-1);
  } else {
    return val+1;
  }  
}

/**
 * returns the id of the pixel below
 */
int surface_projection::p_down( int val ){
  if( val % n_points_z == 0 ){
    return val+(n_points_z-1);
  } else {
    return val - 1;
  }
}


/**
 * returns the id of the pixel forward
 */
int surface_projection::p_back( int val ){
  if( val < n_points_x * n_points_z ){
    return val + ( n_points_x * n_points_z * (n_points_y-1) );
  } else {
    return val - (n_points_x * n_points_z);
  }
}

/**
 * returns the id of the pixel backwards
 */
int surface_projection::p_for( int val ){
  if( val >= (n_points_y-1) * n_points_x * n_points_z ){
    return val - ( n_points_x * n_points_z * (n_points_y-1) );
  } else {
    return val + (n_points_x * n_points_z);
  }
}


/**
 * computes the surface area of the membranes. Outside and inside
 * separately
 */
void surface_projection::compute_surface_area(){

  surface_area = std::vector<double> ( membranes.size(), 0 );

  // iterate over all cells, and check if they are adjacent to a cell
  // of different channel. Always check front, up, and right
  for( unsigned int ii=0; ii<grid.size(); ii++){

    int u = channel.at( p_up(ii)  );
    int f = channel.at( p_for(ii) );
    int r = channel.at( p_right(ii) );
    int d = channel.at( p_down(ii)  );
    int b = channel.at( p_back(ii) );
    int l = channel.at( p_left(ii) );    
    int ch_id = channel.at( ii );
    
    if( u < 0) u *= -1;
    if( f < 0) f *= -1;
    if( r < 0) r *= -1;
    if( d < 0) d *= -1;
    if( b < 0) b *= -1;
    if( l < 0) l *= -1;           
    if( ch_id < 0 ) ch_id*=-1;


    
    if( ch_id > u ){
      //std::cout << "pix_id=" << ii << " ";
      //std::cout << "met up. chid=" << ch_id << " up=" << u << std::endl;
      surface_area.at( ch_id - 2 ) += dx*dy;
    }
    if( ch_id > f ){
      //std::cout << "pix_id=" << ii << " ";      
      //std::cout << "met for. chid=" << ch_id << " f=" << f << std::endl;      
      surface_area.at( ch_id - 2 ) += dx*dz;
    }
    if( ch_id > r ){
      //std::cout << "pix_id=" << ii << " ";      
      //std::cout << "met right chid=" << ch_id << " r=" << r << std::endl;      
      surface_area.at( ch_id - 2 ) += dy*dz;
    }

    
    if( ch_id > d ){
      //std::cout << "pix_id=" << ii << " ";      
      //std::cout << "met down. chid=" << ch_id << " do=" << d << std::endl;      
      surface_area.at( ch_id - 2 ) += dx*dy;
    }
    if( ch_id > b ){
      //std::cout << "pix_id=" << ii << " ";      
      //std::cout << "met ba. chid=" << ch_id << " back=" << b << std::endl;      
      surface_area.at( ch_id - 2 ) += dx*dz;
    }
    if( ch_id > l ){
      //std::cout << "pix_id=" << ii << " ";      
      //std::cout << "met left. chid=" << ch_id << " lo=" << l << std::endl;      
      surface_area.at( ch_id - 2 ) += dy*dz;
    }            
    
  }

}



/*************************
 *  Getters
 *************************/

std::vector<double> surface_projection::get_channel_volumes() const {
  return volumes;
}

std::vector<double> surface_projection::get_membrane_surface_area() const {
  return surface_area;
}

std::vector<double> surface_projection::get_membranes() const {
  return membranes;
}

std::vector<double> surface_projection::get_projection() const {
  return projection;
}

int surface_projection::get_width() const {
  return n_points_x;
}

int surface_projection::get_height() const {
  return n_points_y;
}

int surface_projection::get_depth() const {
  return n_points_z;
}

const std::vector<std::string> surface_projection::get_surface_choices(){
  return surface_choices;
}

int surface_projection::get_h() const{
  return h;
}

int surface_projection::get_k() const{
  return k;
}

int surface_projection::get_l() const{
  return l;
}

int surface_projection::get_ntucs() const{
  return ntucs;
}

int surface_projection::get_type() const {
  return type;
}

double surface_projection::get_a() const{
  return a;
}

double surface_projection::get_surface_level() const {
  return surface_level;
}

double surface_projection::get_slice_width() const {
  return slice_width;
}

double surface_projection::get_slice_height() const {
  return slice_height;
}

double surface_projection::get_theta() const{
  return theta;
}

double surface_projection::get_phi() const{
  return phi;
}

double surface_projection::get_dx() const{
  return dx;
}

double surface_projection::get_dy() const{
  return dy;
}

double surface_projection::get_dz() const{
  return dz;
}

double surface_projection::get_periodicity_length() const{
  return periodicity_length;
}

std::vector<int> surface_projection::get_grid() const {
  return grid;
}


/*************************
 * Setters
 *************************/

void surface_projection::edit_membrane( int id, double dist, double width ){

  if( 2*id >= membranes.size() ){
    throw invalid_parameter_exception( "Membrane not found" );
  } else if( width < 0 ){    
    throw invalid_parameter_exception( "Membrane width can't be negative" );
  } else {
    if( id == 0 && dist != 0 ){
      throw invalid_parameter_exception( "Main membrane distance must be 0" );
    }
    membranes[2*id] = dist;
    membranes[2*id + 1] = width;
  }

}

void surface_projection::delete_membrane( int id ){
  
  if( 2*id >= membranes.size() ){
    throw invalid_parameter_exception( "Membrane not found" );
  } else if ( id == 0 ){
    throw invalid_parameter_exception( "Main membrane can't be deleted" );
  } else {
    membranes.erase( membranes.begin() + 2*id, membranes.begin()+(2*id+2) );
  }

}

void surface_projection::set_theta( double ang ){
  if( ang > M_PI ){
    theta  = M_PI - tolerance;
    throw invalid_parameter_exception( "Angle too large" );
  } else if( ang < 0 ) {
    theta = 0;
    throw invalid_parameter_exception( "No negative angles allowed" );
  } else {
    theta = ang;
  }
}

void surface_projection::set_phi( double ang ){
  if( ang > 2*M_PI ){
    phi  = 2*M_PI - tolerance;
    throw invalid_parameter_exception( "Angle too large" );
  } else if( ang < 0 ) {
    phi = 0;
    throw invalid_parameter_exception( "No negative angles allowed" );
  } else {
    phi = ang;
  }
}

void surface_projection::set_type(int val ){
  if( val < 0 || val >= surface_choices.size() ){
    val = 0;
    throw invalid_parameter_exception( "Surface not available" );
  } else {
    type = val;
  }
}

void surface_projection::set_ntucs( int val ){
  if( val <= 0 ){
    ntucs = 1;
    throw invalid_parameter_exception("At least one unit cell must be projected");
  } else {
    ntucs = val;
  }
  
  //update geometry and arrrays
  update_geometry();
}

void surface_projection::set_slice_width ( double val ){
  if( val < 0 ){
    slice_width = 0.1;
    throw invalid_parameter_exception("Slice width can't be negative");
  } else {
    slice_width = val;
  }
}

void surface_projection::set_slice_height ( double val ){
  //check for invalid parameter
  if( periodicity_length > 0 &&
      (slice_height < 0 || slice_height > 1 ) ){
    //set valid value
    slice_height = 0;
    throw invalid_parameter_exception("periodic orientation, slice height must be in [0, 1]");
  }
  //set value
  slice_height = val;
}


void surface_projection::set_surface_level( double val ){
  surface_level = val;
}

void surface_projection::set_a ( double val ){
  if( val < tolerance ){
    a = 1;
    inv_a = 2*M_PI/(a); // period for nodal representations     
    throw invalid_parameter_exception("Unit cell can't be smaller than 0!");
  } else {
    a = val;
    inv_a = 2*M_PI/(a); // period for nodal representations     
  }


}


void surface_projection::set_n_points_x( int val ){
  if( val <= 0){
    n_points_x = 50;
    //recompute the resolution
    dx = L / n_points_x;
    throw invalid_parameter_exception("Must have a minimum of 1 point");
  } else {
    n_points_x = val;
    //recompute the resolution
    dx = L / n_points_x;
  }
}

void surface_projection::set_n_points_y( int val ){
  if( val <= 0){
    n_points_y = 50;
    //recompute the resolution
    dy = L / n_points_y;      
    throw invalid_parameter_exception("Must have a minimum of 1 point");
  } else {    
    n_points_y = val;
    //recompute the resolution
    dy = L / n_points_y;      
  }   
}

void surface_projection::set_n_points_z( int val ){
  if( val <= 0){
    n_points_z = 50;
    //recompute the resolution
    dz = slice_width / n_points_z;      
    throw invalid_parameter_exception("Must have a minimum of 1 point");
  } else {
    n_points_z = val;
    //recompute the resolution
    dz = slice_width / n_points_z;      
  }
}


void surface_projection::set_h( int val ){
  if( val == 0 && k == 0 && l == 0){
    h = 1;
    throw invalid_parameter_exception("At least one of hkl must be non-zero");
  }
  h = val;
}

void surface_projection::set_k( int val ){
  if( h == 0 && val == 0 && l == 0){
    k = 1;
    throw invalid_parameter_exception("At least one of hkl must be non-zero");
  }
  k = val;
}

void surface_projection::set_l( int val ){
  if( h == 0 && k == 0 && val == 0){
    l = 1;
    throw invalid_parameter_exception("At least one of hkl must be non-zero");
  }
  l = val;
}

void surface_projection::set_periodicity_length( double val ){
  periodicity_length = val;
}
