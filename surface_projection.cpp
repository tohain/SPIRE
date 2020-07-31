#include "surface_projection.hpp"



/*
 * homeotopic thinning
 * skeletsination
 * medial axis
 * network generation
 * deformation retract
 *
 * Vanessa robins
 * morse graph/ morse complex
 */



#ifdef USE_CGAL


/**
 * This struct is copy & paste from a CGAL example for surface
 * triangulations. Not quite sure what it does, but I think it defines
 * something like the niehgborhood of points
 */
struct Perimeter {
  double bound;
  
  Perimeter(double bound): bound(bound){
  }


  // overload the () operator of this class making it a Functor
  template <typename AdvancingFront, typename Cell_handle>
  double operator() (const AdvancingFront& adv, Cell_handle& c,
                     const int& index) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return adv.smallest_radius_delaunay_sphere (c, index);
    }
    
    // If perimeter > bound, return infinity so that facet is not used
    double d  = 0;
    d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                              c->vertex((index+2)%4)->point()));
    if(d>bound) return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();
    // Otherwise, return usual priority value: smallest radius of
    // delaunay sphere
    return adv.smallest_radius_delaunay_sphere (c, index);
  }
};

// Some typedefs for CGAL
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3  Point_3;
typedef K::Vector_3 Vector_3;
typedef std::array<std::size_t,3> Facet;




// overload the output operator for facets 
namespace std {
  std::ostream&
  operator<<(std::ostream& os, const Facet& f)
  {
    os << "3 " << f[0] << " " << f[1] << " " << f[2];
    return os;
  }
}


#endif

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


/** This function prints all points on the surface of the given cahnnel to a
 * gnuplottable file
 */
void surface_projection::print_channel_surface_points( int mem_id, std::string fn ){

  std::vector<int> ch_s_points = get_surface_points( mem_id );
  
  std::ofstream out ( fn );
  out << " # x y z" << std::endl;

  for( auto it : ch_s_points ){
    out << points[3*it] << " " << points[3*it+1]
	<< " " << points[3*it+2] << std::endl;
  }


  out.close();
}

/** This function saves all points making up the topological network
 * of the inner or outer channel of the structure
 */
void surface_projection::print_topological_network( int which, std::string fn ){

  std::vector< std::unordered_set<int> > net_points = get_channel_network();

  std::ofstream out ( fn );
  out << " # x y z " << std::endl;

  for( auto it : net_points.at( which ) ){
    out << points[3*it] << " " << points[3*it+1]
	<< " " << points[3*it+2] << std::endl;
  }

  out.close();  
}


/** Standard constructor initialize with the standard values and
 *  derive some more quantities
 */
surface_projection::surface_projection( double &p, std::string &stat) : ntucs(1), slice_width(1), slice_height(0), a(3,1), inv_a(3,2*M_PI), L(3,1), L_2(3,0.5), n_points_x(76), n_points_y(76), n_points_z(76), type(2), h(0), k(0), l(1), surface_level( 0.0f ), progress(p), status( stat ), s_tables() {


  //set the channel proportion to 0.5
  set_channel_vol_prop( 0.5 );
  
  //update orientation from hkl
  //sets theta,phi
  set_orientation_from_hkl();
  
  //periodicity
  //sets periodicty length
  update_periodicity_length();

  //update geometry
  //sets dx,dy,dz, L
  update_geometry();

  //resize containers
  update_containers();
  
  //add main membrane. This can't be deleted
  membranes.push_back( 0 );
  membranes.push_back( 0.02 );

}


surface_projection::~surface_projection(){

}

/**
 * Quick and dirty dot product between two vectors
 */
double surface_projection::dot_prod ( std::vector<double> v, std::vector<double> w ){
  assert( v.size() == 3 && w.size() == 3 );
  return v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
}


/**
 * Quick and dirty matrix-vector product
 */
std::vector<double> surface_projection::dot_prod( Matrix m, std::vector<double> v ){
  assert( v.size() == 3 );
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
double surface_projection::level_set_gyroid( double x, double y, double z, std::vector<double> inv_a) {
  return sin(inv_a[0]*x)*cos(inv_a[1]*y) + sin(inv_a[1]*y)*cos(inv_a[2]*z) + cos(inv_a[0]*x)*sin(inv_a[2]*z);
}

/**
 * Nodal representation of a Diamond Pn\bar(3)m surface. From
 *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
 * (1991) 83: 407. https://doi.org/10.1007/BF01313411
 */
double surface_projection::level_set_diamond( double x, double y, double z, std::vector<double> inv_a) {
  return cos(inv_a[0]*x)*cos(inv_a[1]*y)*cos(inv_a[2]*z) - sin(inv_a[0]*x)*sin(inv_a[1]*y)*sin(inv_a[2]*z);
}
 
/**
 * Nodal representation of a Primitive Im\bar(3)m surface. From
 *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
 * (1991) 83: 407. https://doi.org/10.1007/BF01313411
 */
double surface_projection::level_set_primitive( double x, double y, double z, std::vector<double> inv_a) {
  return cos(inv_a[0]*x)+cos(inv_a[1]*y)+cos(inv_a[2]*z);
}

/**
 * For debugging purposes: just a simple layer
 *
 */
double surface_projection::level_set_layer( double x, double y, double z, std::vector<double> inv_a) {
  return mod(z, 1./a[2]);
}

/**
 * For debugging purposes: just a simple sphere
 *
 */
double surface_projection::level_set_sphere( double x, double y, double z, std::vector<double> inv_a) {
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
        _phi = M_PI;
    else
        _phi = 0;
  } else {
    if( h == 0 ){
      if( k >= 0 ){
	_phi = 0.5*M_PI;
      } else {
	_phi = 1.5*M_PI;
      }
    } else if( h > 0){
        if (k >= 0)
            _phi = atan2( n[1], n[0] ) ;
        else 
            _phi = 2*M_PI - atan2( -n[1], n[0] ) ;
    } else{
        if (k > 0)
            _phi = M_PI - atan2( n[1], -n[0] ) ;
        else 
            _phi = M_PI + atan2( -n[1], -n[0] ) ;
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

  long max_points = n_points_x*n_points_y*n_points_z;
  
  //the total index of that particle in the array
  int ind = 0;
  double iy = -L_2[1];
  
  //a single loop would be nicer, but i find this less confusing  
  for(unsigned int ii=0; ii<n_points_y; ii++){ //height, the row index (vertical)
    double jx = -L_2[0];
    for(unsigned int jj=0; jj<n_points_x; jj++){ //width, the column index! (horizontal)

      double kz;

      //check if we are periodic. If so, slice_height will be handeled differently
      if( periodicity_length == -1 ){
	//aperiodic, slice_height is an absolute length

	kz = (slice_height - 0.5*slice_width );

      } else {
	//periodic. slice_height is the fraction of the periodicity length 

	//periodicity length is in fractions of a
        kz = periodicity_length*(slice_height - 0.5*slice_width);

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
	  progress = ind / double(max_points);
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
  //the surface is sampled it's highly unlikely that a point has
  //exactly the wanted level set constraint. This is a bit iffy, since
  //the minimal thickness of the membrane is not given in length
  //units, but this level set value. However, this is the only way to do it
  double min_width = 0.02; 
  
  // the level set value of the present points
  double level = 0; 

  status = "Computing membrane";
  
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

    //mark if the pixel is in the "outside" or "inside" channel
    //respective to the main membrane
    if( level > surface_level ){
      //outside
      channel[int(ii/3.)] = -1;
    } else if ( level <= surface_level ) {
      //inside
      channel[int(ii/3.)] = 1;
    } 

    progress = ii/(double(points.size()));
  }


  // next step is to compute the distance transform. We need to label
  // all points according to which channel they are in also. We do
  // that by comparing the distances to the main membrane with the
  // distance map

  status = "Computing distance map";
  progress = 0;
  // compute the distance transform of the grid
  dt.set_parameters( grid, std::vector<unsigned int> {n_points_x, n_points_y, n_points_z},
		     std::vector<double> {dx, dy, dz});  
  dt.compute_distance_map();
  std::vector<float> dmap = dt.get_distance_map();  

  // now assign channel number

  status = "Computing channels";
  
  // get the boundaries of all the membranes
  std::vector<double> mem_pos ( membranes.size(), 0);
  for(unsigned int ii=0; ii<membranes.size(); ii+=2){
    mem_pos[ii] = (membranes[ii] + membranes[ii+1]/2.);
    mem_pos[ii+1] = (membranes[ii] - membranes[ii+1]/2.);
  }
  mem_pos.push_back( -std::numeric_limits<double>::max() );
  mem_pos.push_back( std::numeric_limits<double>::max() );  

  //sort them in ascending order
  std::sort( mem_pos.begin(), mem_pos.end() );

  //now iterate over all pixels and assign each pixel a channel number
  for( unsigned int ii=0; ii<channel.size(); ii++){

    //check if the distance of the current pixel is within the
    //boundaries of the current membrane. The channel array already
    //holds the sign of the distance value to the main membrane

    //start with the inner most shell
    int ch_id=0;

    //iterate over all the channels
    for( ; ch_id < mem_pos.size()-1; ch_id++ ){     
      if( channel[ii] * sqrt(dmap[ii]) >= mem_pos[ch_id] &&
	  channel[ii] * sqrt(dmap[ii]) < mem_pos[ch_id+1])
	  break;
    }    
    
    //found the appropriate membrane, save the channel number
    //the sign still marks, if in or outside of main membrane
    channel[ii] = (ch_id+1)*channel[ii];

    progress = ii/(0.5*channel.size());
  }

  //now iterate over all pixels again and mark the ones, which are
  //within a secondary membrane
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
	
	if( real_distance > (membranes[jj] - membranes[jj+1]/2.) &&
	    real_distance < (membranes[jj] + membranes[jj+1]/2.) &&
	    channel[ii] > 0 ){
	  
	  mark_point = true;
	}
      } else {
	//inside

	if( real_distance > (-membranes[jj] - membranes[jj+1]/2.) &&
	    real_distance < (-membranes[jj] + membranes[jj+1]/2.) &&
	    channel[ii] < 0 ){
	  
	  mark_point = true;
	}	
      }

    }//end cycle membranes

    //mark the point if applicable
    if( mark_point ){
      grid[ii] = 1;
    }

    progress = 0.5 + (ii/double(grid.size()));
  }
  
}

/**
 * This function projects the 3D grid into a 2D picture. Basically
 * only adds up all pixel values in each stack
 */
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
  memset( projection.data(), 0, sizeof(float) * projection.size());
}

/**
 * This function updates the geometry, i.e. recomputes all values
 * depending on the unit cell size and number of unit cells.
 */
void surface_projection::update_geometry(){

  //update box
  L[0] = ntucs * a[0];
  L[1] = ntucs * a[1];

  L_2[0] = L[0]/2.0;
  L_2[1] = L[1]/2.0;

  //update points number
  dx = L[0] / n_points_x;
  dy = L[1] / n_points_y;
  if( periodicity_length == -1 ){
    dz = (slice_width) / n_points_z;
  } else {
    dz = (periodicity_length*slice_width) / n_points_z;
  }

}

/**
 * This function computes the projection. It is basically just a
 * wrapper around the functions performing the single steps. It might
 * be more efficient to do that in a single loop, but with modern
 * compiler optimizations and loop unrolling this should not make too
 * much of a difference
 */
void surface_projection::compute_projection( ){

  progress = 0;
  
  //get the points in the slice  
  status = "Computing the points";
  set_up_points();

  progress = 0.3;
  
  //reset grid
  memset( grid.data(), 0, sizeof(int) * grid.size() );
  
  //get grid
  status = "Setting up the grid";
  set_grid();

  progress = 0.6;
  
  //get projection
  status = "Computing Projection";
  memset( projection.data(), 0, sizeof(float) * projection.size() );
  project_grid();
  progress = 1.0;
  status = "Ready";
  progress = 1;
}




/**
 * This function computes the topological network of a channel by
 * using homotopic thinning.
 *
 * Though in theory you could compute the network of the channels
 * inbetween the membranes, in most cases that would not make too much
 * sense, so so far only the innermost and outermost channels are
 * computed
 *
 * IMORTANT: make sure \ref set_grid() has been called after any
 * parameter change affecting the grid before calling this function to
 * make sure the distance map is up to date!
 */
void surface_projection::compute_channel_network(){

  
  
  // bring in the homotpoic thinning code
  homotopic_thinning<short> ht ( n_points_x, n_points_y, n_points_z, channel, dt.get_distance_map() );
  
  // resize container
  topological_network.resize( 2 );

  // do the work
  topological_network.at(0) = ht.find_channel_skeleton( 1 );
  topological_network.at(1) = ht.find_channel_skeleton( membranes.size()+1 );
}



/**
 * Computes the minimal diameter of the provided channel. Make sure
 * the distance map is up to date!
 *
 * \param[in] channel_id Choices: 0->inner channel, 1->outer channel
 */
double surface_projection::get_minimal_channel_diameter( int channel_id ){
  
  if( channel_id >= topological_network.size() ){
    return -1;
  } else {
    
    auto dmap = dt.get_distance_map();
    
    double min = std::numeric_limits<double>::max();
    for( auto it : topological_network.at( channel_id ) ) {
      
      if( sqrt( dmap[it] ) < min ){
	min = sqrt( dmap[it] );
      }
      
    }
    return min;
  }
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
  
  return n;
}

/**
 * Returns the modulo of two floating point numbers. That's
 * numerically not a nice thing, since we have to take rounding errors
 * into account.
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
 * Computes the periodicity length in units of a. That means, how
 * far the plane must be moved perpendicular to its normal vector so
 * it will reach another chrystallographic equivalent position. This
 * is implemented by looking for a scaling factor for the normal
 * vector so it will point to a point periodically equivalent to the
 * origin (where we started)
 */
void surface_projection::update_periodicity_length(){

  //get the normal of the current orientation
  std::vector<double> n = get_normal();

  // the step size is chosen, so that one step will take at least one
  // direction from one periodic box to the next. Everything in
  // between can't be periodic anyways
  double step_size;
  if( std::fabs( n[0] ) > tolerance ){
    step_size = a[0] / n[0];
  } else if ( std::fabs( n[1] ) > tolerance ){
    step_size = a[1] / n[1];
  } else {
    step_size = a[2] / n[2];
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

    double ddxx = mod(xx, a[0]);
    double ddyy = mod(yy, a[1]);
    double ddzz = mod(zz, a[2]);
    
    if( std::fabs( ddxx ) < tolerance &&
	std::fabs( ddyy ) < tolerance &&
	std::fabs( ddzz ) < tolerance  ){
      more_steps = false;
      break;
    }

    //not there yet, keep on going
    step_count++;
  }

  
  if( step_count == max_steps ){
    //no succes
    periodicity_length = -1;
  } else {
    //success
    periodicity_length = (step_count * step_size);
  }
}


/**
 * This function basically only converts the projection array to a
 * rescaled array, where the minimal value is set to 0 and the maximum
 * to 255. Also inverts if necessary. Memory management could be
 * imporved I guess
 */
unsigned char* surface_projection::get_image(bool invert, std::string scaling){

  //new image array
  unsigned char* img = new unsigned char[ projection.size() ]();
  //(unsigned char*) malloc( sizeof(unsigned char) * projection.size() );
  
  //find min and max value
  double max= *std::max_element( projection.begin(), projection.end() );
  double min= *std::min_element( projection.begin(), projection.end() );

  if( scaling == "LOG" ){
    max = log( max + 1 );
    min = log( min + 1 );    
  }

  //write image data
  for(unsigned int ii=0; ii<projection.size(); ii++){
    //scale
    float pixel_val;
    short u_scaled;

    if( scaling == "LIN" ){
      pixel_val = ((projection[ii] - min)/(max - min)*255);
    } else if( scaling == "LOG" ){
      pixel_val = log(projection[ii]+1);
      pixel_val = ((pixel_val - min)/(max - min))*255;
    }
    u_scaled = static_cast<short>( pixel_val );
      //invert
    if(invert)
      u_scaled = 255 - u_scaled;
    //write in array
    img[ii] = u_scaled;
  }  

  return img;  
}


/**
 * This function calls \ref get_image to generate an image from the
 * \ref projection array, but additionally adds a scale
 */
unsigned char* surface_projection::get_image_with_scale(std::string loc, bool invert ){

  // the scale bar should have a constant length in terms of the
  // picture size. NOT PIXELS, meaning not matter how many pixels
  // there are, the scale is always, say 10% of the edge length of the
  // picture and oriented towards the x direction (horizontal)

  double bar_length = 0.2;

  //compute the pixel size
  int bar_pix_l = bar_length * get_width();
  int bar_pix_w = 0.1 * bar_length * get_height();

  //find the positions
  int margin_x, margin_y;
  if( loc[0] == 't' ){
    margin_y = 0.1 * get_height();
  } else {
    margin_y = get_height() - (0.1 * get_height()) - bar_pix_w;
  }

  if( loc[1] == 'l' ){
    margin_x = 0.1 * get_width();
  } else {
    margin_x = get_width() - (0.1 * get_width()) - bar_pix_l;
  }
  


  //get the image
  unsigned char* img = get_image( invert );

  //make the bar white or black
  for( unsigned int yy=margin_y; yy < margin_y + bar_pix_w; yy++ ){
    for( unsigned int xx=margin_x; xx < margin_x + bar_pix_l; xx++ ){
      if( invert )
	img[ get_width()*yy + xx ] = 0;
      else
	img[ get_width()*yy + xx ] = 255;
    }
  }


  // now the important part: compute the actual length of the bar. We
  // should take the orientation of the slice into account
  

  
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

  status = "Computing channel volumes";
  progress = 0.0;
  
  volumes = std::vector<double> (membranes.size() + 1, 0 );
  
  double vox_vol = dx*dy*dz;
  int ind;
  
  for( unsigned int ii=0; ii<channel.size(); ii++){
    
    ind = channel[ii];
    if( ind < 0 ) ind*=-1;
    volumes[ind-1]+=vox_vol;

    progress = (double) ii / channel.size();
  }


  progress = 1.0;
  status = "Ready";
  
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
 * computes the surface area of the membranes. Uses CGAL, in CASE
 * that's not found, don't do anything. This works best - or rather
 * works only if the voxels are cubes, i.e. are having the same length
 * in all three spatial dimensions.
 *
 * In future I will probably implement a lookup table, if CGAL is not
 * available
 *
 * TODO: find a good choice for the peremiter parameter!
 *
 */
void surface_projection::compute_surface_area(){


  status = "Computing surface area";
  
  surface_area = std::vector<double> ( int( 0.5*membranes.size() ), 0 );

#ifdef USE_CGAL
  
  for( unsigned int ii=0; ii<int(0.5*membranes.size()); ii++){

    // get all the points making  up the membrane
    auto membrane_points = get_surface_points( (2*ii)+1, 6 );


    /*
    std::ofstream p_out ( "surface_points_" + std::to_string( ii ) + ".dat" );
    for( auto it : membrane_points ){
      p_out << points[3*it] << " " << points[3*it+1] << " " << points[3*it+2] << std::endl;
    }
    p_out.close();
    */
    
    
    // copy and paste the points into a CGAL compatible array
    std::vector<Point_3> CGAL_membrane_points;
    for(int jj=0; jj<membrane_points.size(); jj++){
      double x = points.at( 3 * membrane_points.at(jj) );
      double y = points.at( 3 * membrane_points.at(jj) + 1 );
      double z = points.at( 3 * membrane_points.at(jj) + 2 );
      CGAL_membrane_points.push_back( Point_3( x, y, z) );
    }


    // get the surface triangulation
    std::vector<Facet> facets;
    Perimeter perimeter( 0.5 );
    CGAL::advancing_front_surface_reconstruction(CGAL_membrane_points.begin(),
						 CGAL_membrane_points.end(),
						 std::back_inserter(facets),
						 perimeter);
    

    //compute surface area from triangulation
    double total_area = 0;
    
    for( auto it : facets ){
      Vector_3 ab = CGAL_membrane_points.at( it[0] ) - CGAL_membrane_points.at( it[1] );
      Vector_3 ac = CGAL_membrane_points.at( it[0] ) - CGAL_membrane_points.at( it[2] );    
      Vector_3 cross = CGAL::cross_product( ab, ac );
      double area = 0.5 * sqrt( cross.squared_length() );
      total_area+=area;
    }

    //safe area
    surface_area[ii] = total_area;


    
    // for debugging: output surfaces
    std::ofstream triang_out ( "membrane_" + std::to_string(ii) + ".off" );
    triang_out << "OFF\n" << CGAL_membrane_points.size() << " " << facets.size() << " 0\n";
    std::copy(CGAL_membrane_points.begin(),
	      CGAL_membrane_points.end(),
	      std::ostream_iterator<Point_3>(triang_out, "\n"));
    std::copy(facets.begin(),
	      facets.end(),
	      std::ostream_iterator<Facet>(triang_out, "\n"));
    triang_out.close();
    
  }



  status = "Ready";
  progress = 1;

#endif

  
}

/*
 * This function finds all the points which make up the interface to
 * the next, outer membrane
 */
std::vector<int> surface_projection::get_surface_points( int ch_id, int n ){

  std::vector<int> interface_points;

  for( unsigned int ii=0; ii<channel.size(); ii++){

    if( std::abs(channel[ii]) == std::abs( ch_id ) ){
    
      iterable_voxel point ( ii, n_points_x, n_points_y, n_points_z );
      
      std::unordered_set<int> nbs;
      if( n == 6 ){
	nbs = point.get_6_neighbors();
      } else if( n == 18 ){
	nbs = point.get_18_neighbors();
      } else if ( n == 26 ){
	nbs = point.get_26_neighbors();
      } else {
	throw std::string("no adjancy model found for n=" + std::to_string(n) );
      }

      bool found_background = false;
      for( auto it : nbs ){
	if( std::abs( channel[it] ) > std::abs( ch_id ) ){
	  found_background = true;
	  break;
	}
      }
      
      
      if( found_background ){
	interface_points.push_back( point() );
      }

    }
  }

  return interface_points;
  
}






/*************************
 *  Getters
 *************************/

std::vector<short> surface_projection::get_channel() const {
  return channel;
}

std::vector<double> surface_projection::get_channel_volumes() const {
  return volumes;
}

std::vector< std::unordered_set<int> > surface_projection::get_channel_network() const {
  return topological_network;
}

std::vector<double> surface_projection::get_membrane_surface_area() const {
  return surface_area;
}

std::vector<float> surface_projection::get_distance_map() const {
  return dt.get_distance_map();
}

std::vector<double> surface_projection::get_membranes() const {
  return membranes;
}

std::vector<float> surface_projection::get_projection() const {
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

const std::vector<std::string> surface_projection::get_img_scaling_choices(){
  return img_scaling_choices;
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

std::vector<double> surface_projection::get_a() const{
  return a;
}

double surface_projection::get_surface_level() const {
  return surface_level;
}

double surface_projection::get_channel_prop() const {
  double val =  s_tables.get_prop( surface_choices[type], surface_level );
  return val;
}

double surface_projection::get_slice_width() const {
  return slice_width;
}

double surface_projection::get_slice_height() const {
  return slice_height;
}

std::vector<double> surface_projection::get_L() const {
  return L;
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

std::vector<short> surface_projection::get_grid() const {
  return grid;
}

std::vector<float> surface_projection::get_points() const {
  return points;
}


/*************************
 * Setters
 *************************/

void surface_projection::edit_membrane( int id, double dist, double width ){

  if( 2*id >= membranes.size() ){
    throw invalid_parameter_exception( "Membrane not found" );
  } else if( width < 0 ){    
    throw invalid_parameter_exception( "Membrane width can't be negative" );
  } else if ( id == 0 && width < 0.02 ) {
    membranes[2*id + 1] = 0.02;
    throw invalid_parameter_exception( "Width of main membrane can't be less than 0.02" );
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

void surface_projection::set_membranes( std::vector<double> mems ){
  membranes = mems;
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


/*
 * this function is a wrapper for the surface_level parameter. Since
 * this is a rather technical/mathematical parameter, this function
 * takes as an argument the porortion of the *inner* channel and looks
 * up/interpolates the appropriate surface_level value in a table
 */
void surface_projection::set_channel_vol_prop( double vol ){
    surface_level = s_tables.get_level( surface_choices[ type ], vol );    
}



void surface_projection::set_surface_level( double val ){
  surface_level = val;
}

void surface_projection::set_a ( double val ){
  if( val < tolerance ){
    a[0] = 1;
    a[1] = 1;
    a[2] = 1;
    inv_a[1] = 2*M_PI/(a[0]); // period for nodal representations
    inv_a[2] = 2*M_PI/(a[1]);
    inv_a[3] = 2*M_PI/(a[2]);
    throw invalid_parameter_exception("Unit cell can't be smaller than 0!");
  } else {
    a[0] = val;
    a[1] = val;
    a[2] = val;
    inv_a[0] = 2*M_PI/(a[0]); // period for nodal representations
    inv_a[1] = 2*M_PI/(a[1]); // period for nodal representations        
    inv_a[2] = 2*M_PI/(a[2]); // period for nodal representations
  }
}


void surface_projection::set_a ( double ax, double ay, double az ){
  if( ax < tolerance || ay < tolerance || az < tolerance ){
    a[0] = 1;
    a[1] = 1;
    a[2] = 1;
    inv_a[0] = 2*M_PI/(a[0]); // period for nodal representations
    inv_a[1] = 2*M_PI/(a[1]); // period for nodal representations
    inv_a[2] = 2*M_PI/(a[2]); // period for nodal representations
    throw invalid_parameter_exception("Unit cell can't be smaller than 0!");
  } else {
    a[0] = ax;
    a[1] = ay;
    a[2] = az;
    inv_a[0] = 2*M_PI/(ax); // period for nodal representations
    inv_a[1] = 2*M_PI/(ay); // period for nodal representations        
    inv_a[2] = 2*M_PI/(az); // period for nodal representations
  }
}

void surface_projection::set_n_points_x( int val ){
  if( val <= 0){
    n_points_x = 50;
    //recompute the resolution
    dx = L[0] / n_points_x;
    throw invalid_parameter_exception("Must have a minimum of 1 point");
  } else {
    n_points_x = val;
    //recompute the resolution
    dx = L[0] / n_points_x;
  }
}

void surface_projection::set_n_points_y( int val ){
  if( val <= 0){
    n_points_y = 50;
    //recompute the resolution
    dy = L[1] / n_points_y;      
    throw invalid_parameter_exception("Must have a minimum of 1 point");
  } else {    
    n_points_y = val;
    //recompute the resolution
    dy = L[1] / n_points_y;      
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


