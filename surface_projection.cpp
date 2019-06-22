#include "surface_projection.hpp"


// eport as pgm image
#include "img_out.hpp"

/// Standard constructor initialize with the standard values and derive some more quantities
surface_projection::surface_projection() : ntucs(1), slice_width(0.1), slice_height(0.5), mem_width(0.2), a(1), n_points_x(50), n_points_y(50), n_points_z(50), type(2), h(0), k(0), l(1) {

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

 

double surface_projection::level_set_gyroid( double x, double y, double z, double a) {
  double s = (2*M_PI)/(a);
  return sin(s*x)*cos(s*y) + sin(s*y)*cos(s*z) + cos(s*x)*sin(s*z);
}

double surface_projection::level_set_diamond( double x, double y, double z, double a) {
  double s = (2*M_PI)/(a);
  return cos(s*x)*cos(s*y)*cos(s*z) - sin(s*x)*sin(s*y)*sin(s*z);
}

double surface_projection::level_set_primitive( double x, double y, double z, double a) {
  double s = (2*M_PI)/(a);
  return cos(s*x)+cos(s*y)+cos(s*z);
}

/**
 * Computes the rotation angles needed to get the orientation given by
 * the Miller indeces hkl. Only works if all indeces are positive!
 */
void surface_projection::set_orientation_from_hkl(){

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


  //rotate vector, so it will aligned with yz plane
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


  //apply the angles
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
  
  //a single loop would be nicer, but i find this less confusing  
  for(unsigned int ii=0; ii<n_points_y; ii++){ //height, the row index (vertical)
    for(unsigned int jj=0; jj<n_points_x; jj++){ //width, the column index! (horizontal)
      for(unsigned int kk=0; kk<n_points_z; kk++){ //depth

	//the total index of that particle in the array
	int ind=ii*(n_points_x*n_points_z) + jj*n_points_z + kk; 

	//now compute the absolute position of all voxels
	double x,y,z;

	//check if we are periodic. If so, slice_height will be handeled differently


	if( periodicity_length == -1 ){
	  //aperiodic. slice_height is an absolute length
	  x = ((jj*dx)-L/2.)*nx[0] + ((ii*dy)-L/2.)*ny[0] + ((kk*dz)-L/2)*nz[0] + (slice_height+L/2.)*nz[0] - 0.5*slice_width*nz[0];
	  y = ((jj*dx)-L/2.)*nx[1] + ((ii*dy)-L/2.)*ny[1] + ((kk*dz)-L/2)*nz[1] + (slice_height+L/2.)*nz[1] - 0.5*slice_width*nz[1];
	  z = ((jj*dx)-L/2.)*nx[2] + ((ii*dy)-L/2.)*ny[2] + ((kk*dz)-L/2)*nz[2] + (slice_height+L/2.)*nz[2] - 0.5*slice_width*nz[2];


	} else {
	  //periodic. slice_height is the fraction of the periodicity length
	  x = ((jj*dx)-L/2.)*nx[0] + ((ii*dy)-L/2.)*ny[0] + ((kk*dz)-L/2)*nz[0] + ((periodicity_length*slice_height)+L/2.)*nz[0] - 0.5*periodicity_length*slice_width*nz[0];
	  y = ((jj*dx)-L/2.)*nx[1] + ((ii*dy)-L/2.)*ny[1] + ((kk*dz)-L/2)*nz[1] + ((periodicity_length*slice_height)+L/2.)*nz[1] - 0.5*periodicity_length*slice_width*nz[1];
	  z = ((jj*dx)-L/2.)*nx[2] + ((ii*dy)-L/2.)*ny[2] + ((kk*dz)-L/2)*nz[2] + ((periodicity_length*slice_height)+L/2.)*nz[2] - 0.5*periodicity_length*slice_width*nz[2];
	}

	
	//assign the position
	points[3*ind]=x;
	points[(3*ind)+1]=y;
	points[(3*ind)+2]=z;
      }
    }
  }
  
}


void surface_projection::set_grid (){
  
  for( unsigned int ii=0; ii<points.size(); ii+=3 ){

    double level = 0;
    if( type == 0 ){ //gyroid
      level = level_set_gyroid( points[ii], points[ii+1], points[ii+2], a );
    } else if( type == 1 ){ //diamon
      level = level_set_diamond( points[ii], points[ii+1], points[ii+2], a );
    } else if( type == 2 ){ //primitive
      level = level_set_primitive( points[ii], points[ii+1], points[ii+2], a );
    }
    
    if( level < mem_width && level > -mem_width ){
      grid[int(ii/3)] = 1;
    }

  }

}


void surface_projection::project_grid (){


  for(unsigned int ii=0; ii<n_points_y; ii++){ //y direcion=height (vertical)
    for(unsigned int jj=0; jj<n_points_x; jj++){//x direction=width (horizontal)
      for(unsigned int kk=0; kk<n_points_z; kk++){
        
	int ind = ii *( n_points_x * n_points_z ) + jj * n_points_z + kk;

	if(grid[ind] == 1){
	  projection[(ii*n_points_x)+jj] += 1;
	}
	
      }
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
  memset( points.data(), 0, sizeof(double) * points.size() );
  
  /// Array holding the color (electron density) of the voxels
  grid.resize( n_points_x * n_points_y * n_points_z, 0 );
  memset( grid.data(), 0, sizeof(int) * grid.size() );
  
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

  //update points number
  dx = L / n_points_x;
  dy = L / n_points_y;
  dz = slice_width / n_points_z;  
}


void surface_projection::compute_projection( ){

  /*
  std::cout << "computing with:\n";
  std::cout << "theta=" << theta << " phi=" << phi << " L=" << L << std::endl;
  std::cout << "slice_width=" << slice_width << " slice_height=" << slice_height << std::endl;
  std::cout << "mem_width=" << mem_width <<  " a=" << a << std::endl;
  std::cout << "type=" << type << std::endl;
  */
    //get the points in the slice    
    set_up_points();
    
    //reset grid
    memset( grid.data(), 0, sizeof(int) * grid.size() );

    //get grid
    set_grid();
    
    //get projection
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

  //stop after some steps
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




/*************************
 *  Getters
 *************************/

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

double surface_projection::get_mem_width() const {
  return mem_width;
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

/*************************
 * Setters
 *************************/


void surface_projection::set_theta( double ang ){
  if( ang > M_PI )
   theta  = M_PI;
  else if( ang < 0 )
    theta = 0;
  else
    theta = ang;
}

void surface_projection::set_phi( double ang ){
  if( ang > 2*M_PI )
    phi  = 2*M_PI;
  else if( ang < 0 )
    phi = 0;
  else
    phi = ang;
}

void surface_projection::set_type(int val ){
  type = val;
}

void surface_projection::set_ntucs( int val ){
  if( val < 0 )
    ntucs = 1;
  else
    ntucs = val;

  //update geometry and arrrays
  update_geometry();
}

void surface_projection::set_slice_width ( double val ){
  if( val < 0 )
    slice_width = 0;
  else
    slice_width = val;
}

void surface_projection::set_slice_height ( double val ){
  slice_height = val;
}


void surface_projection::set_mem_width ( double val ){
  if( val < 0 )
    mem_width = 0;
  else
    mem_width = val;
}

void surface_projection::set_a ( double val ){
  if( val < 0 )
    a = 1.0;
  else
    a = val;
}


void surface_projection::set_n_points_x( int val ){
  if( val < 0)
    n_points_x = 1;
  else
    n_points_x = val;

  //recompute the resolution
  dx = L / n_points_x;
}

void surface_projection::set_n_points_y( int val ){
  if( val < 0)
    n_points_y = 1;
  else
    n_points_y = val;

  //recompute the resolution
  dy = L / n_points_y;  
}

void surface_projection::set_n_points_z( int val ){
  if( val < 0)
    n_points_z = 1;
  else
    n_points_z = val;

  //recompute the resolution
  dz = slice_width / n_points_z;  
}


void surface_projection::set_h( int val ){
  if( val < 0 )
    h = 0;
  else
    h = val;  
}

void surface_projection::set_k( int val ){
  if( val < 0 )
    k = 0;
  else
    k = val;  
}

void surface_projection::set_l( int val ){
  if( val < 0 )
    l = 0;
  else
    l = val;  
}

void surface_projection::set_periodicity_length( double val ){
  periodicity_length = val;
}
