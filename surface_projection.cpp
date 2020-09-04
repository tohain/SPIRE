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
  return msg.c_str();//"Invalid parameter choice!";
}

const std::string invalid_parameter_exception::details() const {
  return msg;
}

/**
 * writes all parameters to an ASCI file, so certain settings can be
 * reloaded 
 */
void surface_projection::write_parameters( std::string outfile ){

  std::ofstream out ( outfile );

  out << " # parameters" << std::endl;;


  out << "Lx=" << L[0] << std::endl
      << "Ly=" << L[1] << std::endl
      << "Lz=" << L[2] << std::endl;

  out << "slice_position=" << slice_position << std::endl;
  
  out << "ax=" << a[0] << std::endl
      << "ay=" << a[1] << std::endl
      << "az=" << a[2] << std::endl;  

  out << "uc_scale_ab=" << uc_scale_ab << std::endl
      << "uc_scale_c=" << uc_scale_c << std::endl;

  out << "surface_level=" << surface_level << std::endl;

  out << "<membranes>" << std::endl;
  for( unsigned int ii=0; ii < membranes.size(); ii+=2 ){
    out << "pos=" << membranes[ii] << ";width=" << membranes[ii+1] << std::endl;
  }
  out << "</membranes>" << std::endl;

  out << "channel_filled=";
  for(unsigned int ii=0; ii<channel_filled.size(); ii++){
    out << channel_filled[ii];
    if( ii < channel_filled.size()-1 )
      out<<",";
  }
  out << std::endl;

  out << "h=" << h << std::endl
      << "k=" << k << std::endl
      << "l=" << l << std::endl;

  out << "n_points_x=" << n_points_x << std::endl
      << "n_points_z=" << n_points_z << std::endl;

  out << "type=" << type << std::endl;
  
  out.close();
}

/**
 * reads all parameters from an ASCI file
 */
void surface_projection::read_parameters( std::string infile ){

  std::ifstream in ( infile );

  // read file line by line
  std::string line;  
  while( !in.eof() ){

    std::getline( in, line );

    //check if line starts with a comment
    int char_pos = 0;
    while( line[char_pos] == ' ' )
      char_pos++;
    if( char_pos < line.size() && line[char_pos] != '#' ){
      // get content

      auto line_split = my_utility::str_split( line, '=' );

      if( line_split[0] == "Lx" ){
	L[0] = std::stod( line_split[1] );
      }
      if( line_split[0] == "Ly" ){
	L[1] = std::stod( line_split[1] );
      }
      if( line_split[0] == "Lz" ){
	L[2] = std::stod( line_split[1] );
      }

      if( line_split[0] == "slice_position" ){
	slice_position = std::stod( line_split[1] );
      }

      if( line_split[0] == "ax" ){
	a[0] = std::stod( line_split[1] );
      }
      if( line_split[0] == "ay" ){
	a[1] = std::stod( line_split[1] );
      }
      if( line_split[0] == "az" ){
	a[2] = std::stod( line_split[1] );
      }


      if( line_split[0] == "uc_scale_ab" ){
	uc_scale_ab = std::stod( line_split[1] );
      }
      if( line_split[0] == "uc_scale_c" ){
	uc_scale_c = std::stod( line_split[1] );
      }

      if( line_split[0] == "surface_level" ){
	surface_level = std::stod( line_split[1] );
      }

      /*
       * membranes
       */
      if( line_split[0] == "<membranes>" ){

	membranes.clear();
	
	// scan until all membranes are read
	std::getline( in, line );
	do{
	  auto mem = my_utility::str_split(line, ';' );

	  membranes.push_back( std::stod( my_utility::str_split(mem[0], '=' )[1] ) );
	  membranes.push_back( std::stod( my_utility::str_split(mem[1], '=' )[1] ) );

	  std::getline( in, line );
	} while( line != "</membranes>" );

      }

      /*
       * channel_filled
       */

      if( line_split[0] == "channel_filled" ){
	channel_filled.clear();
	auto vals = my_utility::str_split( line_split[1], ',' );
	for( unsigned int ii=0; ii<vals.size(); ii++ ){
	  channel_filled.push_back( std::stoi( vals[ii] ) );
	} 
      }


      if( line_split[0] == "h" ){
	h = std::stoi( line_split[1] );
      }
      if( line_split[0] == "k" ){
	k = std::stoi( line_split[1] );
      }
      if( line_split[0] == "l" ){
	l = std::stoi( line_split[1] );
      }

      if( line_split[0] == "n_points_x" ){
	n_points_x = std::stoi( line_split[1] );
      }
      if( line_split[0] == "n_points_z" ){
	n_points_z = std::stoi( line_split[1] );
      }

      if( line_split[0] == "type" ){
	type = std::stoi( line_split[1] );
      }
      

    }
    
  }


  int nr_channel_soll = membranes.size() + 1;

  if( nr_channel_soll != channel_filled.size() ){
    throw std::string("There was an error reading the parameter file; parameters not completly restored");
    
  }
  
  update_geometry();

  compute_projection();

  in.close();
  
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
surface_projection::surface_projection( double &p, std::string &stat) : slice_position(0), a(3,1), inv_a(3,2*M_PI), uc_scale_ab( 1.0 ), uc_scale_c( 1.0 ), L(3,1), L_2(3,0.5), n_points_x(76), n_points_y(76), n_points_z(76), type(2), h(0), k(0), l(1), uc_dim_in_orientation(3, 1), surface_level( 0.0f ), progress(p), status( stat ), s_tables() {

  //set the channel proportion to 0.5
  set_channel_vol_prop( 0.5 );
  
  //update orientation from hkl
  //sets theta,phi
  set_orientation_from_hkl();
  
  //periodicity
  //sets periodicty length
  compute_uc_dim_in_orientation();

  //update geometry
  //sets dx,dy,dz, L
  update_geometry();

  //resize containers
  update_containers();
  
  //add main membrane. This can't be deleted
  membranes.push_back( 0 );
  membranes.push_back( 0.02 );

  // add the first two channels
  channel_filled.push_back( 0 ); // inner channel
  channel_filled.push_back( 1 ); // main membrane
  channel_filled.push_back( 0 ); // outer channel
}


surface_projection::~surface_projection(){

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

double surface_projection::level_set_wurtzite_0_05( double x, double y, double z, std::vector<double> inv_a) {

  double freq_x = inv_a[0];
  double freq_y = inv_a[1];
  double freq_z = inv_a[2];

  return  + 0.000277634 - 0.362363*cos(2*+freq_z*z) + 0.0532718*cos(6*+freq_z*z) + 0.223887*cos(2*+freq_y*y) + 0.225075*sin(2*+freq_y*y)*sin(+freq_z*z) + 0.0848054*cos(2*+freq_y*y)*cos(2*+freq_z*z) - 0.0526959*sin(2*+freq_y*y)*sin(3*+freq_z*z) - 0.0578834*sin(2*+freq_y*y)*sin(5*+freq_z*z) + 0.0688335*cos(4*+freq_y*y)*cos(4*+freq_z*z) - 0.0614875*cos(6*+freq_y*y)*cos(2*+freq_z*z) + 0.447937*cos(2*freq_x*x)*cos(+freq_y*y) - 0.45065*cos(2*freq_x*x)*sin(+freq_y*y)*sin(+freq_z*z) + 0.170029*cos(2*freq_x*x)*cos(+freq_y*y)*cos(2*+freq_z*z) + 0.106884*cos(2*freq_x*x)*sin(+freq_y*y)*sin(3*+freq_z*z) - 0.0802985*cos(2*freq_x*x)*cos(+freq_y*y)*cos(4*+freq_z*z) + 0.117166*cos(2*freq_x*x)*sin(+freq_y*y)*sin(5*+freq_z*z) - 0.065981*cos(2*freq_x*x)*cos(+freq_y*y)*cos(6*+freq_z*z) - 0.0589097*cos(2*freq_x*x)*sin(+freq_y*y)*sin(7*+freq_z*z) + 0.215803*cos(2*freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z) - 0.0570703*cos(2*freq_x*x)*cos(5*+freq_y*y) - 0.0687129*cos(2*freq_x*x)*sin(5*+freq_y*y)*sin(+freq_z*z) - 0.0620785*cos(2*freq_x*x)*cos(7*+freq_y*y)*cos(4*+freq_z*z) + 0.108815*cos(4*freq_x*x)*cos(2*+freq_z*z) - 0.0992979*cos(4*freq_x*x)*sin(2*+freq_y*y)*sin(3*+freq_z*z) + 0.138495*cos(4*freq_x*x)*cos(2*+freq_y*y)*cos(4*+freq_z*z) - 0.0568124*cos(4*freq_x*x)*cos(4*+freq_y*y) + 0.0686862*cos(4*freq_x*x)*sin(4*+freq_y*y)*sin(+freq_z*z) - 0.0736126*cos(4*freq_x*x)*cos(6*+freq_y*y)*cos(2*+freq_z*z) - 0.0564133*cos(6*freq_x*x)*cos(+freq_y*y) + 0.067969*cos(6*freq_x*x)*sin(+freq_y*y)*sin(+freq_z*z) - 0.12331*cos(6*freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z) - 0.0627811*cos(6*freq_x*x)*cos(5*+freq_y*y)*cos(4*+freq_z*z) - 0.0620622*cos(8*freq_x*x)*cos(2*+freq_y*y)*cos(4*+freq_z*z);

}

double surface_projection::level_set_wurtzite_0_075( double x, double y, double z, std::vector<double> inv_a) {

  double freq_x = inv_a[0];
  double freq_y = inv_a[1];
  double freq_z = inv_a[2];  

  return  + 0.000277634 - 0.362363*cos(2*+freq_z*z) + 0.223887*cos(2*+freq_y*y) + 0.225075*sin(2*+freq_y*y)*sin(+freq_z*z) + 0.0848054*cos(2*+freq_y*y)*cos(2*+freq_z*z) + 0.447937*cos(2*freq_x*x)*cos(+freq_y*y) - 0.45065*cos(2*freq_x*x)*sin(+freq_y*y)*sin(+freq_z*z) + 0.170029*cos(2*freq_x*x)*cos(+freq_y*y)*cos(2*+freq_z*z) + 0.106884*cos(2*freq_x*x)*sin(+freq_y*y)*sin(3*+freq_z*z) - 0.0802985*cos(2*freq_x*x)*cos(+freq_y*y)*cos(4*+freq_z*z) + 0.117166*cos(2*freq_x*x)*sin(+freq_y*y)*sin(5*+freq_z*z) + 0.215803*cos(2*freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z) + 0.108815*cos(4*freq_x*x)*cos(2*+freq_z*z) - 0.0992979*cos(4*freq_x*x)*sin(2*+freq_y*y)*sin(3*+freq_z*z) + 0.138495*cos(4*freq_x*x)*cos(2*+freq_y*y)*cos(4*+freq_z*z) - 0.12331*cos(6*freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z);

}


double surface_projection::level_set_wurtzite_0_1( double x, double y, double z, std::vector<double> inv_a) {

  double freq_x = inv_a[0];
  double freq_y = inv_a[1];
  double freq_z = inv_a[2];

  return  + 0.000277634 - 0.362363*cos(2*+freq_z*z) + 0.223887*cos(2*+freq_y*y) + 0.225075*sin(2*+freq_y*y)*sin(+freq_z*z) + 0.447937*cos(2*freq_x*x)*cos(+freq_y*y) - 0.45065*cos(2*freq_x*x)*sin(+freq_y*y)*sin(+freq_z*z) + 0.170029*cos(2*freq_x*x)*cos(+freq_y*y)*cos(2*+freq_z*z) + 0.106884*cos(2*freq_x*x)*sin(+freq_y*y)*sin(3*+freq_z*z) + 0.117166*cos(2*freq_x*x)*sin(+freq_y*y)*sin(5*+freq_z*z) + 0.215803*cos(2*freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z) + 0.108815*cos(4*freq_x*x)*cos(2*+freq_z*z) + 0.138495*cos(4*freq_x*x)*cos(2*+freq_y*y)*cos(4*+freq_z*z) - 0.12331*cos(6*freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z);

}

double surface_projection::level_set_wurtzite_0_2( double x, double y, double z, std::vector<double> inv_a) {

  double freq_x = inv_a[0];
  double freq_y = inv_a[1];
  double freq_z = inv_a[2];  

  return  + 0.000277634 - 0.362363*cos(2*+freq_z*z) + 0.223887*cos(2*+freq_y*y) + 0.225075*sin(2*+freq_y*y)*sin(+freq_z*z) + 0.447937*cos(2*freq_x*x)*cos(+freq_y*y) - 0.45065*cos(2*freq_x*x)*sin(+freq_y*y)*sin(+freq_z*z) + 0.215803*cos(2*freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z);

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
  //Matrix Ry = get_y_rot_m(2*M_PI - theta), Rz = get_z_rot_m(2*M_PI - phi);
  Matrix Ry = VEC_MAT_MATH::get_y_rot_m(theta), Rz = VEC_MAT_MATH::get_z_rot_m(phi);

  //rotate
  nx = VEC_MAT_MATH::dot_prod( Rz, VEC_MAT_MATH::dot_prod( Ry, nx));
  ny = VEC_MAT_MATH::dot_prod( Rz, VEC_MAT_MATH::dot_prod( Ry, ny));
  nz = VEC_MAT_MATH::dot_prod( Rz, VEC_MAT_MATH::dot_prod( Ry, nz));  

  std::cout << "nz=(" << nz[0] << "," << nz[1] << "," << nz[2] << ")" << std::endl;
  
  long max_points = n_points_x*n_points_y*n_points_z;
  
  //the total index of that particle in the array
  int ind = 0;
  double iy = -L_2[1];
  
  //a single loop would be nicer, but i find this less confusing  
  for(unsigned int ii=0; ii<n_points_y; ii++){ //height, the row index (vertical)
    double jx = -L_2[0];
    for(unsigned int jj=0; jj<n_points_x; jj++){ //width, the column index! (horizontal)

      double kz;

      // the start of the slice in absolute length units
      kz = slice_position - 0.5*L[2];
      
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


/*
 * This function computes the dimension of the unit cell of the
 * current structure in the current orientation. Basically just a
 * wrapper around the \ref compute_periodicity function.
 */
void surface_projection::compute_uc_dim_in_orientation(){
 
  // base vectors of slice
  std::vector<double> nx = {1, 0, 0};
  std::vector<double> ny = {0, 1, 0};
  std::vector<double> nz = {0, 0, 1};  
  
  //rotation matrices
  //substract the angle from 2pi since we're rotating in mathematical
  //negative orientation (with the clock
  //Matrix Ry = get_y_rot_m(2*M_PI - theta), Rz = get_z_rot_m(2*M_PI - phi);
  Matrix Ry = VEC_MAT_MATH::get_y_rot_m(theta), Rz = VEC_MAT_MATH::get_z_rot_m(phi);

  //rotate
  nx = VEC_MAT_MATH::dot_prod( Rz, VEC_MAT_MATH::dot_prod( Ry, nx));
  ny = VEC_MAT_MATH::dot_prod( Rz, VEC_MAT_MATH::dot_prod( Ry, ny));
  nz = VEC_MAT_MATH::dot_prod( Rz, VEC_MAT_MATH::dot_prod( Ry, nz));  

  // get the periodicities along the rotate base vectors
  uc_dim_in_orientation[0] = compute_periodicity( nx );
  uc_dim_in_orientation[1] = compute_periodicity( ny );
  uc_dim_in_orientation[2] = compute_periodicity( nz );

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
    } else if( type == 3 ){ //wurtzite_0.05
      level = level_set_wurtzite_0_05( points[ii], points[ii+1], points[ii+2], inv_a );
    } else if( type == 4 ){ //wurtzite_0.075
      level = level_set_wurtzite_0_075( points[ii], points[ii+1], points[ii+2], inv_a );
    } else if( type == 5 ){ //wurtzite_0.1
      level = level_set_wurtzite_0_1( points[ii], points[ii+1], points[ii+2], inv_a );
    } else if( type == 6 ){ //wurtzite_0.2
      level = level_set_wurtzite_0_2( points[ii], points[ii+1], points[ii+2], inv_a );            
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

  set_orientation_from_hkl();
  compute_uc_dim_in_orientation();
  update_a();
  
  L_2[0] = L[0]/2.0;
  L_2[1] = L[1]/2.0;
  L_2[2] = L[2]/2.0; //not really needed, just for completeness

  //update nr. of pixels
  set_n_points_y_to_unitcell();
  
  //update points number
  dx = L[0] / n_points_x;
  dy = L[1] / n_points_y;
  dz = L[2] / n_points_z;  
}

/**
 * This function computes the projection. It is basically just a
 * wrapper around the functions performing the single steps. It might
 * be more efficient to do that in a single loop, but with modern
 * compiler optimizations and loop unrolling this should not make too
 * much of a difference
 */
void surface_projection::compute_projection( ){


  update_containers();
  
  progress = 0;
  
  //get the points in the slice  
  status = "Computing the points";
  set_up_points();



  progress = 0.3;
  
  //reset grid
  memset( grid.data(), 0, sizeof(short) * grid.size() );
  
  //get grid
  status = "Setting up the grid";
  set_grid();

  progress = 0.6;

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
  Matrix Ry = VEC_MAT_MATH::get_y_rot_m(2*M_PI - theta), Rz = VEC_MAT_MATH::get_z_rot_m(phi);

  //rotate
  n = VEC_MAT_MATH::dot_prod( Rz, VEC_MAT_MATH::dot_prod( Ry, n));  
  
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


/*
 * This function computes the periodicity length in a given
 * direction. That means, how far must the current plane (given by the
 * normal vector n) be moven along its normal vector n so it will
 * reach another chrystallographic equivalent layer. This is
 * implemented by looking for a scaling factor for the normal vector
 * so it will point to a point periodically equivalent to the origin
 * (where we started)
 */
double surface_projection::compute_periodicity( std::vector<double> n ){

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
    return -1;
  } else {
    //success
    return (step_count * step_size);
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
 * sets the "color" of a membrane. This means all voxels belonging to
 * that membrane are set to "0" or "1". This way it can be chosen, if
 * the channel contributes towards the projection
 *
 * \param[in] mem_id the id of the channel to color (starting at 1)
 * \param[in] val the value to set the channel voxels to (0,1)
 */
void surface_projection::set_channel_color( int mem_id, int val ){

  channel_filled[mem_id-1] = val;
  
  for( unsigned int ii=0; ii<channel.size(); ii++ ){
    if( std::abs( channel[ii]) == std::abs( mem_id ) ){
      grid[ii] = val;
    }
  }  
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

int surface_projection::get_type() const {
  return type;
}

std::vector<double> surface_projection::get_a() const{
  return a;
}

double surface_projection::get_uc_scale_ab() const {
  return uc_scale_ab;
}

double surface_projection::get_uc_scale_c() const {
  return uc_scale_c;
}

double surface_projection::get_surface_level() const {
  return surface_level;
}

double surface_projection::get_channel_prop() const {
  double val =  s_tables.get_prop( surface_choices[type], surface_level );
  return val;
}

double surface_projection::get_slice_width() const {
  return L[2];
}

double surface_projection::get_slice_height() const {
  return L[1];
}

double surface_projection::get_slice_length() const {
  return L[0];
}

double surface_projection::get_slice_position() const {
  return slice_position;
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

std::vector<double> surface_projection::get_uc_dim_in_orientation() const{
  return uc_dim_in_orientation;
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
  update_channel_fill_container();
}

void surface_projection::update_channel_fill_container(){

  // clear previous data
  channel_filled.clear();
  channel_filled.resize( membranes.size() + 1, 0 );

  for( unsigned int ii=0; ii<channel_filled.size(); ii++){
    if( ii % 2 == 1 )
      channel_filled[ii] = 1;
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
  } else {
    type = val;
  }
}

void surface_projection::set_slice_width ( double val ){
  if( val < 0 ){
    L[2] = 0.1;
    throw invalid_parameter_exception("Slice width can't be negative");
  } else {
    L[2] = val;
  }
}


void surface_projection::set_slice_length ( double val ){
  L[0] = val;
}


void surface_projection::set_slice_height ( double val ){
  L[1] = val;
}

void surface_projection::set_slice_position ( double val ){
  slice_position = val;
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

void surface_projection::set_uc_scale_ab( double val ){
  uc_scale_ab = val;
}

void surface_projection::set_uc_scale_c( double val ){
  uc_scale_c = val;
}

void surface_projection::update_a(){
    a[0] = uc_scale_ab * unitcell_dim[type][0];
    a[1] = uc_scale_ab * unitcell_dim[type][1];
    a[2] = uc_scale_c * unitcell_dim[type][2];
    inv_a[0] = 2*M_PI/(a[0]); // period for nodal representations
    inv_a[1] = 2*M_PI/(a[1]); // period for nodal representations        
    inv_a[2] = 2*M_PI/(a[2]); // period for nodal representations  
}


void surface_projection::set_n_points_x( int val ){
  if( val <= 0){
    n_points_x = 50;
    //recompute the resolution
    dx = L[0] / n_points_x;
    throw invalid_parameter_exception("n_points_x: Must have a minimum of 1 point");
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
    throw invalid_parameter_exception("n_points_y: Must have a minimum of 1 point");
  } else {    
    n_points_y = val;
    //recompute the resolution
    dy = L[1] / n_points_y;      
  }   
}

void surface_projection::set_n_points_y_to_unitcell(){

  int val = n_points_x * ( L[1] / L[0] );

  if( val <= 0){
    n_points_y = 50;
    //recompute the resolution
    dy = L[1] / n_points_y;      
    throw invalid_parameter_exception("n_points_y: Must have a minimum of 1 point");
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
    dz = L[2] / n_points_z;      
    throw invalid_parameter_exception("n_points_z: Must have a minimum of 1 point");
  } else {
    n_points_z = val;
    //recompute the resolution
    dz = L[2] / n_points_z;      
  }
}

void surface_projection::set_n_points_z_to_unitcell(){

  int val = n_points_x * ( L[2] / L[0] );

  if( val <= 0){
    n_points_y = 50;
    //recompute the resolution
    dz = L[2] / n_points_z;
    throw invalid_parameter_exception("n_points_z: Must have a minimum of 1 point");
  } else {    
    n_points_z = val;
    //recompute the resolution
    dy = L[2] / n_points_z;
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


std::vector<int> surface_projection::get_channel_fill() const{
  return channel_filled;
}
