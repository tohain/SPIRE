


#ifndef VEC_MAT_MATH_HPP
#define VEC_MAT_MATH_HPP

#include <vector>
#include <cmath>
#include <iostream>

/// Quick and dirty 3x3 Matrix object, consisting of 3 rows
class Matrix {

public:  

  Matrix(){};
  Matrix( std::vector<double> x, std::vector<double> y , std::vector<double> z );
  

  std::vector<double> v;
  std::vector<double> w;
  std::vector<double> z;

  void print();
  
};


class VEC_MAT_MATH {

public:

  /// Dot product of two vectors
  static double dot_prod ( std::vector<double> v, std::vector<double> w );

  /// Quick and dirty 3x3 matric by vector multiplication
  static std::vector<double> dot_prod( Matrix m, std::vector<double> v );
  static std::vector<double> dot_prod( std::vector<double> v, Matrix m );  

  /// Quick and dirty cross produckt
  static std::vector<double> cross_prod( std::vector<double> a, std::vector<double> b );
  
  /// Quick and dirty matrix product
  static Matrix dot_prod( Matrix m, Matrix n );
  
  /// Creates and returns a rotation matrix around x axis
  static Matrix get_x_rot_m( double ang );

  /// Creates and returns a rotation matrix around y axis
  static Matrix get_y_rot_m( double ang );  

  /// Creates and returns a rotation matrix around z axis
  static Matrix get_z_rot_m( double ang );

  /// returns two vectors, creating a orthogonal base
  static std::vector< std::vector<double> > get_orthogonal_base( std::vector<double> n );

  /// returns unit length vector
  static std::vector<double> get_unit( std::vector<double> v );

  
};


#endif
