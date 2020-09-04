


#ifndef VEC_MAT_MATH_HPP
#define VEC_MAT_MATH_HPP

#include <vector>
#include <cmath>

/// Quick and dirty 3x3 Matrix object, consisting of 3 rows
typedef struct {
  std::vector<double> v;
  std::vector<double> w;
  std::vector<double> z;  
} Matrix;


class VEC_MAT_MATH {

public:

  /// Dot product of two vectors
  static double dot_prod ( std::vector<double> v, std::vector<double> w );

  /// Quick and dirty 3x3 matric by vector multiplication
  static std::vector<double> dot_prod( Matrix m, std::vector<double> v );

  /// Creates and returns a rotation matrix around x axis
  static Matrix get_x_rot_m( double ang );

  /// Creates and returns a rotation matrix around y axis
  static Matrix get_y_rot_m( double ang );  

  /// Creates and returns a rotation matrix around z axis
  static Matrix get_z_rot_m( double ang );

};


#endif
