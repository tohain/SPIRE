#include "surface_projection.hpp"


// eport as pgm image
#include "img_out.hpp"

/// Standard constructor initialize with the standard values and derive some more quantities
surface_projection::surface_projection() : ntucs(1), slice_width(0.2), slice_height(0.5), mem_width(0.2), a(1), theta(0), phi(0), n_points_x(50), n_points_y(50), n_points_z(50), type("Primitive") {

  //update geometry
  update_geometry();
  //resize containers
  update_containers();

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
void surface_projection::set_orientation_from_hkl( int h, int k, int l ){

  std::vector<double> n = {double(h), double(k), double(l)};

  //make it unit length
  double norm_scale = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
  n[0]/=norm_scale;n[1]/=norm_scale;n[2]/=norm_scale;
 
  //get angle for rotation in xy plane
  double _phi;
  if( n[1] == 0 ){
    _phi = M_PI/2.;
  } else {
    _phi = atan2( n[0], n[1] ) ;
  }

  //rotate vector, so it will aligned with yz plane
  Matrix Rz = get_z_rot_m( _phi );  
  n = dot_prod( Rz, n );

  //get angle to rotate around x axis to align vector with {0, 0, 1}  
  double _theta;
  if( n[2] == 0 ){
    _theta = M_PI/2.;
  } else {
    _theta = atan2( n[1], n[2] );
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

	//compute the position of this voxel
	double x=((jj*dx)-L/2.)*nx[0] + ((ii*dy)-L/2.)*ny[0] + ((kk*dz)-L/2)*nz[0] + (slice_height+L/2.)*nz[0] - 0.5*slice_width*nz[0];
	double y=((jj*dx)-L/2.)*nx[1] + ((ii*dy)-L/2.)*ny[1] + ((kk*dz)-L/2)*nz[1] + (slice_height+L/2.)*nz[1] - 0.5*slice_width*nz[1];
	double z=((jj*dx)-L/2.)*nx[2] + ((ii*dy)-L/2.)*ny[2] + ((kk*dz)-L/2)*nz[2] + (slice_height+L/2.)*nz[2] - 0.5*slice_width*nz[2];

	
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
    if( type == "Gyroid" ){
      level = level_set_gyroid( points[ii], points[ii+1], points[ii+2], a );
    } else if( type == "Diamond" ){
      level = level_set_diamond( points[ii], points[ii+1], points[ii+2], a );
    } else if( type == "Primitive" ){
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

void surface_projection::set_type(std::string val ){
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
