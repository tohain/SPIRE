
#ifndef SP_PROJ_I
#define SP_PROJ_I


#include <cassert>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstring>
#include <cstdio>
#include <exception>


/** \brief Quick and dirty implementation of an invalid parameter exception
 *
 * This exception is thrown, if a parameter choice would cause
 * faulty maths, unphysical behavior or just makes no sense at all
 */
class invalid_parameter_exception : public std::exception {

public:
  /// Constructor, initializes \ref msg
  invalid_parameter_exception( std::string _msg );

  /// RETURNS a short string what this exception is
  virtual const char* what() const throw();
  /// Returns a string conatining more details
  const std::string details() const;

  /// String holding details on the exception
  std::string msg;
};

/// Quick and dirty 3x3 Matrix object, consisting of 3 rows
struct Matrix {
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> z;  
};



/** \brief The computational core class of this code. Computes the projection calculations and stores al * l parameters.
 *
 * This class contains all neccessary function to compute the 2D
 * projection of a surface. It also stores all parameters. The GUI is
 * basically only a front end to use the code provided in this
 * class. A terminal interface should be easy to implement if wanted.
 */
class surface_projection {
public:

  /// Default Constructor
  surface_projection();
  /// Destructor
  ~surface_projection();

  /// Recomputes the geometric slice properties
  void update_geometry();

  /// Updates container capacities
  void update_containers();

  /// Computes the projection, wrapper around several functions
  void compute_projection();

  /// Computes the periodicity length of the given orientation
  void update_periodicity_length();
  
  /// Computes and sets theta and phi from the Miller indeces
  void set_orientation_from_hkl();

  /// Converts the \ref projection array in a rescaled image array
  unsigned char* get_image(bool invert = false);

  /// Returns the number of points in the slice in x direction
  int get_width() const;
  /// Returns the number of points in the slice in y direction  
  int get_height() const;
  /// Returns the number of points in the slice in z direction  
  int get_depth() const;

  /// Returns h of the Miller indeces of the slice orientation
  int get_h() const;
  /// Returns k of the Miller indeces of the slice orientation  
  int get_k() const;
  /// Returns l of the Miller indeces of the slice orientation  
  int get_l() const;
  /// Returns the size of the slice in nr. of unit cells
  int get_ntucs() const;
  /// Returns the currently set type of the surface
  int get_type() const;
  /// Returns the projection array
  std::vector<double> get_projection() const;

  /// Returns the periodicity (unit cell size) of the surface
  double get_a() const;
  /// Returns the width of the membrane
  double get_mem_width() const;
  /// Returns the width of the slice
  double get_slice_width() const;
  /// Returns the position of the slice along its normal vector
  double get_slice_height() const;
  /// Returns theta of the current orientation
  double get_theta() const;
  /// Returns phi of the current orientation  
  double get_phi() const;
  /// Returns the distance between two points in x direction
  double get_dx() const;
  /// Returns the distance between two points in y direction  
  double get_dy() const;
  /// Returns the distance between two points in z direction  
  double get_dz() const;
  /// Returns the unit cell size (periodicity length) in the current orientation
  double get_periodicity_length() const;  
  /// Returns the available surfaces
  const std::vector<std::string> get_surface_choices();


  /// Sets the theta angle of the orientation of the slice
  void set_theta( double ang );
  /// Sets phi angle of the orientation of the slice
  void set_phi( double ang );
  /// Sets the type of the surface to be projected
  void set_type( int val );
  /// Sets the number of unit cells of the slice
  void set_ntucs( int val );
  /// Sets the slice width
  void set_slice_width( double val );
  /// Sets the position of the slice along ints normal vector
  void set_slice_height( double val );
  /// Sets the membrane width
  void set_mem_width( double val );
  /// Sets the unit cell size of the surface
  void set_a( double val );
  /// Sets the number of points in the slice in x direction
  void set_n_points_x( int val );
  /// Sets the number of points in the slice in y direction  
  void set_n_points_y( int val );
  /// Sets the number of points in the slice in z direction  
  void set_n_points_z( int val );
  /// Sets the \ref h Miller index
  void set_h( int val );
  /// Sets the \ref k Miller index  
  void set_k( int val );
  /// Sets the \ref l Miller index  
  void set_l( int val );  
  /// Sets the \ref periodicity_length
  void set_periodicity_length( double val );


  
private:

  /// Numerical tolerance when two floating point numbers are considered equal
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
  
  /// Nodal approximation of the level set function of a g surface
  double level_set_gyroid( double x, double y, double z, double a);

  /// Nodal approximation of the level set function of a d surface
  double level_set_diamond( double x, double y, double z, double a);

  /// Nodal approximation of the level set function of a p surface
  double level_set_primitive( double x, double y, double z, double a);  
  
  /// Computes the position of the voxels in the slice
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

  /// The number of unit cells the slice covers
  int ntucs;

  /// Edge length of the slice
  double L;
  double L_2; // L_2 = L/2
  
  /// Distance between to points in x direction and voxel size in x direction
  double dx;
  /// Distance between to points in y direction and voxel size in y direction  
  double dy;
  /// Distance between to points in z direction and voxel size in z direction  
  double dz;  

  /** \brief Width of the slice 
   *
   * This parameter is within [0,1] if the current orientation is
   * periodic (i.e. periodicity_length>0) and is then interpreted as a
   * fraction of the periodicity_length. For aperiodic orientations it
   * is given in absolute length units
   */
  double slice_width;
  /** \brief Slice position along its normal vector
   *
   * This parameter is within [0,1] if the current orientation is
   * periodic (i.e. periodicity_length>0) and is then interpreted as a
   * fraction of the periodicity_length. For aperiodic orientations it
   * is given in absolute length units
   */
  double slice_height;

  /** \brief Membrane width
   *
   * This gives the width of the membrane, however, I'm not sure in
   * which units. It is the value the level set conditions are
   * compared to, to see if it is withing the level set conditions
   */
  double mem_width;
  /// Unit cell size of the surface in {100},{010},{001} direciton
  double a;
  //gyroid period
  double inv_a;
  
  /// Theta angle of the orientation of the slice
  double theta;
  /// Phi angle of the orientation of the slice  
  double phi;

  /// h Miller index
  int h;
  /// k Miller index  
  int k;
  /// l Miller index  
  int l;

  /** \brief The length of the unitbox in the current orientation.
   *
   * -1 if there if the structure is aperiodic within reasonable
   * lengthscales
   */
  double periodicity_length;
  
  /// The number of points in the slice in x direction
  unsigned int n_points_x;
  /// The number of points in the slice in y direction  
  unsigned int n_points_y;
  /// The number of points in the slice in z direction  
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