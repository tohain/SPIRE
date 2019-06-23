
#ifndef SP_PROJ_I
#define SP_PROJ_I

#include <cassert>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstring>
#include <cstdio>

  /// Quick and dirty 3x3 Matrix object, consisting of 3 rows
struct Matrix {
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> z;  
};


///enum to label surfaces

class surface_projection {
public:

  //constructor
  surface_projection();
  ~surface_projection();

  /// Updates the projection with new coordinates
  void update_geometry();
  /// Updates container capacities
  void update_containers();

  /// Computes the projection
  void compute_projection();

  /// Computes the periodicity length of the given orientation (should
  /// only use, if orientation given by hkl?)
  void update_periodicity_length();
  

  std::vector<double> get_projection() const;

  unsigned char* get_image(bool invert = false);

  //setters and getters
  int get_width() const;
  int get_height() const;
  int get_depth() const;
  int get_h() const;
  int get_k() const;
  int get_l() const;
  int get_ntucs() const;
  int get_type() const;
  
  double get_a() const;
  double get_mem_width() const;
  double get_slice_width() const;
  double get_slice_height() const;
  double get_theta() const;
  double get_phi() const;
  double get_dx() const;
  double get_dy() const;
  double get_dz() const;  
  double get_periodicity_length() const;  

  
  void set_orientation_from_hkl();
  
  const std::vector<std::string> get_surface_choices();
  
  void set_theta( double ang );
  void set_phi( double ang );  

  void set_type( int val );

  void set_ntucs( int val );
  void set_slice_width( double val );
  void set_slice_height( double val );

  void set_mem_width( double val );
  void set_a( double val );

  void set_n_points_x( int val );
  void set_n_points_y( int val );
  void set_n_points_z( int val );
  
  void set_h( int val );
  void set_k( int val );
  void set_l( int val );  

  void set_periodicity_length( double val );
  
  //private:

  //numerical tolerance
  double tolerance = 1e-5;
  
  /// Dot product of two vectors
  double dot_prod ( std::vector<double> v, std::vector<double> w );

  /// Quick and dirty 3x3 matric by vector multiplication
  std::vector<double> dot_prod( Matrix m, std::vector<double> v );

  /// Quick and dirty modulo between two doubles
  inline double mod( double lhs, double rhs );
  
  /// Returns the normal vector of unit length for the given
  /// orientation
  std::vector<double> get_normal();
  
  /// Creates and returns a rotation matrix around x axis
  Matrix get_x_rot_m( double ang ) const;

  /// Creates and returns a rotation matrix around z axis
  Matrix get_z_rot_m( double ang ) const;   
  
  ///Nodal approximation of the level set function of a g surface
  double level_set_gyroid( double x, double y, double z, double a);

  ///Nodal approximation of the level set function of a d surface
  double level_set_diamond( double x, double y, double z, double a);

  ///Nodal approximation of the level set function of a p surface
  double level_set_primitive( double x, double y, double z, double a);  
  
  /// Initializes the voxels in the slice
  void set_up_points();


  
  /// Evaluates all voxels with the chosen level sets. Sets colors of voxel
  void set_grid();


  
  /// Computes the 2D projection
  void project_grid ();

  

  /////////////////////////////////////////////
  // parameter
  //////////////////////////////////

  /*
   * box parameters
   */

  //how many unitcells in the box?
  int ntucs;

  /// Edge length of cubic "sim"box.
  double L;
  
  double dx;
  double dy;
  double dz;  

  //slice dimension
  double slice_width;
  /// Slice position in range [0, 1]. In relation to the periodicty,
  /// so we'll always move in one unit cell in the given orientation
  double slice_height;

  //membrane thickness
  double mem_width;
  //gyroid unit cell length
  double a;
  
  //slice orientation
  double theta;
  double phi;

  //miller indeces
  int h;
  int k;
  int l;

  /// The length of the unitbox in the current orientation. -1 if there
  ///is none, or it is too large to handle
  double periodicity_length;
  
  //get number of points
  unsigned int n_points_x;
  unsigned int n_points_y;
  unsigned int n_points_z;

  /// Type of surface to project
  int type;


  /// Available surfaces
  const std::vector<std::string> surface_choices = {"Gyroid",
						    "Diamond",
						    "Primitive"};
  
  /******************************************************************
   * data, do not access manually unless you know what you're doing
   ******************************************************************/

  /// Array holding coordinates of the voxels
  std::vector<double> points;
  /// Array holding the color (electron density) of the voxels
  std::vector<int> grid;

  /// The 2D projection
  std::vector<double> projection;
};



#endif
