
#ifndef SP_PROJ_I
#define SP_PROJ_I

#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstring>

  /// Quick and dirty 3x3 Matrix object, consisting of 3 rows
struct Matrix {
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> z;  
};


class surface_projection {
public:

  //constructor
  surface_projection();
  ~surface_projection();

  void compute_projection();

  double* get_projection();

  unsigned char* get_image();

  int get_width();
  int get_height();
  int get_depth();
  
private:

  /// Dot product of two vectors
  double dot_prod ( std::vector<double> v, std::vector<double> w );

  /// Quick and dirty 3x3 matric by vector multiplication
  std::vector<double> dot_prod( Matrix m, std::vector<double> v );
  
  ///Nodal approximation of the level set function of a g surface
  double level_set_gyroid( double x, double y, double z, double a);

  ///Nodal approximation of the level set function of a d surface
  double level_set_diamond( double x, double y, double z, double a);

  ///Nodal approximation of the level set function of a p surface
  double level_set_primitive( double x, double y, double z, double a);  

  /// Quick and dirty auxiliary function 
  template <class T>
  void set_to_zero( std::vector< std::vector<T> > &data );

  template <class T>
  void set_to_zero( std::vector<T> &data );
  
  /// Initializes the voxels in the slice
  void set_up_points();


  
  /// Evaluates all voxels with the chosen level sets. Sets colors of voxel
  void set_grid (std::string type, double d);


  
  /// Computes the 2D projection
  void project_grid ();

  

  /////////////////////////////////////////////
  // parameter
  //////////////////////////////////

    //set up parameters
  double L = 10;
  double dx = .02;
  double dy = .02;
  double dz = .02;  

  //slice dimension
  double slice_width = 1;
  double slice_height = -1;

  //membrane thickness
  double mem_width = 0.5;
  //gyroid unit cell length
  double a = 0.1;
  
  //slice orientation
  double theta = 0*M_PI;
  double phi = 0.06*M_PI;

  //get number of points
  unsigned int n_points_x = int(L/dx);
  unsigned int n_points_y = int(L/dy);
  unsigned int n_points_z = int(slice_width/dz);


  ////////////////
  // data
  ///////////////////

  /// Array holding coordinates of the voxels
  std::vector<double> points;
  /// Array holding the color (electron density) of the voxels
  std::vector<int> grid;

  /// The 2D projection
  std::vector<double> projection;
};



#endif
