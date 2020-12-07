/* Projection tool - compute planar projection of triply periodic
 * minimal surfaces 
 * Copyright (C) 2020 Tobias Hain
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see https://www.gnu.org/licenses.
 */

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
