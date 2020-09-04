
#ifndef SP_PROJ_I
#define SP_PROJ_I


#include <iostream>
#include <cassert>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstring>
#include <cstdio>
#include <exception>
#include "auxiliary.hpp"
#include "vec_mat_math.hpp"

#ifdef USE_CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/tuple.h>
#include <CGAL/Vector_3.h>
#endif

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



#include "img_out.hpp"
#include "distance_transform.hpp"
#include "homotopic_thinning.hpp"
#include "surface_tables.hpp"


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
  surface_projection(double &progress, std::string &status);
  /// Destructor
  ~surface_projection();

  /// Recomputes the geometric slice properties
  void update_geometry();

  /// Updates container capacities
  void update_containers();

  /// Computes the projection, wrapper around several functions
  void compute_projection();

  /// Computes the periodicity length in the given direction
  double compute_periodicity( std::vector<double> n );

  /// Computes the dimension of the unitcell in the current orientation
  void compute_uc_dim_in_orientation();
  
  /// Computes and sets theta and phi from the Miller indeces
  void set_orientation_from_hkl();

  /// Adds a membrane
  void add_membrane( double dist, double width );

  /// sets the "color" of a membrane
  void set_channel_color( int mem_id, int val );

  /// updates the size and resets the filled channels control
  void update_channel_fill_container();
  
  /// Computes the volumes of the channels
  void compute_volume();

  /// Computes the surface area of the membranes
  void compute_surface_area();

  /// Computes the minimal channel diameter of the given channel
  void compute_channel_network();

  /// computes the minimal diamter of a channel.
  double get_minimal_channel_diameter( int channel_id );
  
  /// Returns all the points which are on the surface of the channel
  std::vector<int> get_surface_points( int ch_id, int n = 26 );

  /// Returns the distance map. Be careful if it is up to date. If in
  /// doubt, call \ref set_grid() to update it
  std::vector<float> get_distance_map() const;
  
  /// Converts the \ref projection array in a rescaled image array
  unsigned char* get_image(bool invert = false, std::string scaling = "lin");

  /// Converts the \ref projection array in a rescaled image array and
  /// adds a scale
  unsigned char* get_image_with_scale(std::string loc = "bl" , bool invert = false );
  
  /// Outputs the grid
  void print_grid( std::string fn );
  
  /// Outputs surface points of the given membrane
  void print_channel_surface_points( int mem_id, std::string fn );

  /// Output the points making up the minimal topological network
  void print_topological_network( int which, std::string fn );

  /// write parameters to an ASCI file
  void write_parameters( std::string outfile );

  /// reads parameters to an ASCI file
  void read_parameters( std::string outfile );  
  
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
  /// Returns the currently set type of the surface
  int get_type() const;
  /// Returns the projection array
  std::vector<float> get_projection() const;
  /// Returns the grid
  std::vector<short> get_grid() const;
  /// returns a copy of the channel array
  std::vector<short> get_channel() const;  
  /// Returns the points (=positions) of the points
  std::vector<float> get_points() const;  
  /// returns the channel_filled array
  std::vector<int> get_channel_fill() const;  
  /// Returns the unit cell size the surface in [001] direction
  std::vector<double> get_a() const;
  /// returns unit cell scale in ab direction
  double get_uc_scale_ab() const;
  /// returns unit cell scale in c direction
  double get_uc_scale_c() const;
  /// Returns the surface level
  double get_surface_level() const;
  /// Returns the ratio of the volumes of the two channels
  double get_channel_prop() const;
  /// Returns the membranes
  std::vector<double> get_membranes() const;
  /// Return channel volumes
  std::vector<double> get_channel_volumes() const;
  /// Return the topological network of the channels
  std::vector< std::unordered_set<int> > get_channel_network() const;
  /// Return surface areas of the membranes
  std::vector<double> get_membrane_surface_area() const;
  /// Returns the width of the slice
  double get_slice_width() const;
  /// Returns the height of the slice
  double get_slice_height() const;
  /// Returns the length of the slice
  double get_slice_length() const;  
  /// Returns the position of the slice along its normal vector
  double get_slice_position() const;
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
  std::vector<double> get_uc_dim_in_orientation() const;
  /// Returns the slice dimensions
  std::vector<double> get_L() const;
  /// Returns the available surfaces
  const std::vector<std::string> get_surface_choices();
  /// Returns the available image scaling methods
  const std::vector<std::string> get_img_scaling_choices();

  /// Sets the theta angle of the orientation of the slice
  void set_theta( double ang );
  /// Sets phi angle of the orientation of the slice
  void set_phi( double ang );
  /// Sets the type of the surface to be projected
  void set_type( int val );
  /// Sets the slice width // z 
  void set_slice_width( double val );
  /// Sets the slice heigth // y
  void set_slice_height( double val );
  /// Sets the slice length // x
  void set_slice_length( double val );


  /// Sets the position of the slice along ints normal vector
  void set_slice_position( double val );



  
  /// set membrane
  void edit_membrane( int id, double dist, double width );
  /// deletes a meberane
  void delete_membrane( int id );
  /// set membranes
  void set_membranes( std::vector<double> mems );
  
  /// Sets the surface_level
  void set_surface_level( double val );
  /// Sets the channel proportion
  void set_channel_vol_prop( double val );
  /// recomputes the unitcell dimensions from \ref uc_scale_ab , \ref uc_scape_c and \ref unitcell_dim
  void update_a( );
  /// Sets the unit cell scaling in c direction
  void set_uc_scale_ab ( double val );
  /// Sets the unit cell scaling in c direction
  void set_uc_scale_c ( double val );  
  /// Sets the number of points in the slice in x direction
  void set_n_points_x( int val );
  /// Sets the number of points in the slice in y direction
  void set_n_points_y( int val );
  /// Sets the number of points in the slice in y direction according
  /// to the dimension of the unitcell
  void set_n_points_y_to_unitcell();  
  /// Sets the number of points in the slice in z direction  
  void set_n_points_z( int val );
  /// Sets the number of points in the slice in z direction according
  /// to the dimension of the unitcell
  void set_n_points_z_to_unitcell();    
  /// Sets the \ref h Miller index
  void set_h( int val );
  /// Sets the \ref k Miller index  
  void set_k( int val );
  /// Sets the \ref l Miller index  
  void set_l( int val );  

  
protected:

  /// Numerical tolerance when two floating point numbers are considered equal
  double tolerance = 1e-5;
  
  /// Quick and dirty modulo between two doubles
  inline double mod( double lhs, double rhs );
  
  /// Returns the normal vector of unit length for the given
  /// orientation
  std::vector<double> get_normal();  
  
  /// Nodal approximation of the level set function of a g surface
  double level_set_gyroid( double x, double y, double z, std::vector<double> a);

  /// Nodal approximation of the level set function of a d surface
  double level_set_diamond( double x, double y, double z, std::vector<double> a);

  /// Nodal approximation of the level set function of a p surface
  double level_set_primitive( double x, double y, double z, std::vector<double> a);

  /// Nodal approximation of the level set function of a w surface with tolerance 0.05
  double level_set_wurtzite_0_05( double x, double y, double z, std::vector<double> a);
  /// Nodal approximation of the level set function of a w surface with tolerance 0.075
  double level_set_wurtzite_0_075( double x, double y, double z, std::vector<double> a);
  /// Nodal approximation of the level set function of a w surface with tolerance 0.1
  double level_set_wurtzite_0_1( double x, double y, double z, std::vector<double> a);
  /// Nodal approximation of the level set function of a w surface with tolerance 0.2
  double level_set_wurtzite_0_2( double x, double y, double z, std::vector<double> a);    

  /// for debugging: just a single layer
  double level_set_layer( double x, double y, double z, std::vector<double> a);

  /// for debugging: just a single sphere
  double level_set_sphere( double x, double y, double z, std::vector<double> a);  
  
  /// Computes the position of the voxels in the slice
  void set_up_points();
  
  /// Evaluates all voxels with the chosen level sets. Sets colors of voxel
  void set_grid();
  
  /// Computes the 2D projection
  void project_grid ();

  /// Returns the id of a pixel up
  inline int p_up( int val );
  inline int p_down( int val );
  inline int p_right( int val );
  inline int p_left( int val );
  inline int p_for( int val );
  inline int p_back( int val );  


  
  
  /////////////////////////////////////////////
  // parameters                           ////
  ////////////////////////////////////////////

  /*
   * box parameters
   */

  /// Edge length of the slice
  std::vector<double> L;
  std::vector<double> L_2; // L_2 = L/2
  
  /// Distance between to points in x direction and voxel size in x direction
  double dx;
  /// Distance between to points in y direction and voxel size in y direction  
  double dy;
  /// Distance between to points in z direction and voxel size in z direction  
  double dz;  

  /** \brief Slice position along its normal vector in length units
   *
   * It has the same units as the unit cell dimension
   */
  double slice_position;

  /// Rectangular unit cell size
  std::vector<double> a;
  /// Surface period
  std::vector<double> inv_a;

  /// Scaling of unit cell (ab)
  double uc_scale_ab;
  /// Scaling of unit cell (c) for rectangular/hexagonal symmetries (wurtzite so far)
  double uc_scale_c;

  /// This parameter controlls the proportions of the two different
  /// channels. 0 means equal volume in both of them
  double surface_level;


  /// Holds information for the membranes in pairs of (distance,
  /// width). For n membranes this container has then 2*n entries
  std::vector<double> membranes;

  /// Holds information about the channel color (if it is "filled"
  /// thus contributing to the projection). Membranes are by default
  /// "filled", where as space in between channels are not. Same order
  /// as the membranes vector.
  std::vector<int> channel_filled;
  
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

  /** \brief The length of the unitcell in the given orientation of the slice.
   *
   * In [001] direction, this equals the unitcell_dim. -1 if there if
   * the structure is aperiodic within reasonable lengthscales
   */
  std::vector<double> uc_dim_in_orientation;
  
  /// The number of points in the slice in x direction
  unsigned int n_points_x;
  /// The number of points in the slice in y direction  
  unsigned int n_points_y;
  /// The number of points in the slice in z direction  
  unsigned int n_points_z;

  /// Type of surface to project
  int type;



  /// the distance transform of the grid. This is updated in the \ref
  /// set_grid() function and then can be accessed using the getter
  distance_transform<short> dt;

  /// The volumes of the channels
  std::vector<double> volumes;

  /// The surface area of the membranes
  std::vector<double> surface_area;
  
  /// Available surfaces
  const std::vector<std::string> surface_choices = {"Gyroid",
						    "Diamond",
						    "Primitive",
						    "Wurtzite_0.05",
						    "Wurtzite_0.075",
						    "Wurtzite_0.1",
						    "Wurtzite_0.2",						    
  };

  /// The dimension of the unitcell to keep the symmetry. Now the
  /// unitcell dimension is set using the parameter a above, however,
  /// some membranes are non-cubic with a fixed a/b a/c ratio. This
  /// ratio needs to be supplied here so 
  const std::vector< std::vector<double> > unitcell_dim= { {1.0,1.0,1.0},
							   {1.0,1.0,1.0},
							   {1.0,1.0,1.0},
							   {2.0, sqrt(3.0), 1.732692},
							   {2.0, sqrt(3.0), 1.732692},
							   {2.0, sqrt(3.0), 1.732692},
							   {2.0, sqrt(3.0), 1.732692},							   
  };

  /// Available image scalings
  const std::vector<std::string> img_scaling_choices = {"LIN",
							"LOG"};

  /// Lookup tables for surface properties
  SURFACE_TABLES s_tables;
  
  /******************************************************************
   * data, do not access manually unless you know what you're doing
   ******************************************************************/

  /// Array holding coordinates of the voxels
  std::vector<float> points;
  /// Array holding the color (electron density) of the voxels
  std::vector<short> grid;
  /// Array denoting the channel a voxel is in. Membranes are
  /// considered channels. 1 is the innermost channel, 2 the
  /// innermoste membrane, ... (except for the sphere surface, in this
  /// case the order is inverted for whatever reason - I guess wierd
  /// level set values.
  /// negative/positive values denotes if the channel
  /// is inside or outside of the main "membrane" (=the membrane set
  /// by the level set constraint)
  std::vector<short> channel;

  /// This vector holds the minimal topological networks of the
  /// channels
  std::vector< std::unordered_set<int> > topological_network;
  
  /// A variable where the progress of the computations are stored in
  double &progress;
  /// A string the current status of the code is written into. Max len=200
  std::string &status;
  
  /// The 2D projection
  std::vector<float> projection;
};


#endif
