/* Projection tool - compute planar projections of triply periodic
 * minimal surfaces 
 * Copyright (C) 2021 Tobias Hain
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
template<class T = double>
class Matrix {

public:  

  Matrix(){};
  Matrix( std::vector<T> x, std::vector<T> y , std::vector<T> z ) : v(x), w(y), z(z) {};
  
  void print(){
    std::cout << "/" << v[0] << " " << v[1] << " " << v[2] << "\\" << std::endl;
    std::cout << "|" << w[0] << " " << w[1] << " " << w[2] << "|" << std::endl;
    std::cout << "\\" << z[0] << " " << z[1] << " " << z[2] << "/" << std::endl;
  }


  std::vector<T> v;
  std::vector<T> w;
  std::vector<T> z;
};


class VEC_MAT_MATH {

public:

  
  /// Dot product of two vectors
  template<class T>
  static T dot_prod ( std::vector<T> v, std::vector<T> w ){
    return v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
  }
  /// Dot product of two vectors
  template<class T, class M>
  static double dot_prod ( std::vector<T> v, std::vector<M> w ){
    return v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
  }

  /// scalar product of scalar and vector
  template<class T>  
  static std::vector<double> s_prod ( double s, std::vector<T> w ){
    return std::vector<double> { s*w[0], s*w[1], s*w[2] };  
  }
  
  /// sum of two vector
  template<class T>
  static std::vector<T> vec_add( std::vector<T> v, std::vector<T> w ){
    return std::vector<T> { v[0]+w[0], v[1]+w[1], v[2]+w[2] };  
  }
  
  /// Quick and dirty 3x3 matric by vector multiplication
  template<class T>
  static std::vector<T> dot_prod( Matrix<T> m, std::vector<T> v ){
    std::vector<T> r (3, 0);
    r[0] = dot_prod( m.v, v);
    r[1] = dot_prod( m.w, v);
    r[2] = dot_prod( m.z, v);
    return r;
  }

  /// Quick and dirty 3x3 matric by vector multiplication
  template<class T, class M>
  static std::vector<double> dot_prod( Matrix<T> m, std::vector<M> v ){
    std::vector<double> r (3, 0);
    r[0] = dot_prod( m.v, v);
    r[1] = dot_prod( m.w, v);
    r[2] = dot_prod( m.z, v);
    return r;
  }


  
  template<class T>
  static std::vector<T> dot_prod( std::vector<T> v, Matrix<T> m ){
    
    std::vector<T> n_col_1 = {m.v[0], m.w[0], m.z[0]};
    std::vector<T> n_col_2 = {m.v[1], m.w[1], m.z[1]};
    std::vector<T> n_col_3 = {m.v[2], m.w[2], m.z[2]};  
    
    std::vector<T> r (3, 0);
    r[0] = dot_prod( v, n_col_1 );
    r[1] = dot_prod( v, n_col_2 );
    r[2] = dot_prod( v, n_col_3);
    return r;
  }

  /// Quick and dirty cross produckt
  template<class T>
  static std::vector<T> cross_prod( std::vector<T> a, std::vector<T> b ){

    std::vector<T> out (3,0);
    
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
    
    return out;
  }

    /// Quick and dirty cross produckt
  template<class T, class M>
  static std::vector<double> cross_prod( std::vector<T> a, std::vector<M> b ){

    std::vector<double> out (3,0);
    
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
    
    return out;
  }

  
  /// Quick and dirty matrix product
  template<class T>
  static Matrix<T> dot_prod( Matrix<T> m, Matrix<T> n ){
    
    std::vector<T> n_col_1 = {n.v[0], n.w[0], n.z[0]};
    std::vector<T> n_col_2 = {n.v[1], n.w[1], n.z[1]};
    std::vector<T> n_col_3 = {n.v[2], n.w[2], n.z[2]};  
    
    std::vector<T> prod_row_1 = { dot_prod( m.v, n_col_1 ),
      dot_prod( m.v, n_col_2 ),
      dot_prod( m.v, n_col_3 ) };
    std::vector<T> prod_row_2 = { dot_prod( m.w, n_col_1 ),
      dot_prod( m.w, n_col_2 ),
      dot_prod( m.w, n_col_3 ) };
    std::vector<T> prod_row_3 = { dot_prod( m.z, n_col_1 ),
      dot_prod( m.z, n_col_2 ),
      dot_prod( m.z, n_col_3 ) };
  
    return Matrix<T>( prod_row_1, prod_row_2, prod_row_3 );
  }
  
  

  /** \brief Creates and returns a rotation matrix around x axis
   * Returns a rotation matrix. The rotation is by ang against the
   * clock around the x-Axis
   */
  static Matrix<double> get_x_rot_m( double ang ){
    Matrix<double> R;
    R.v = { 1, 0,           0          };
    R.w = { 0, cos(ang), -sin(ang) };
    R.z = { 0, sin(ang),  cos(ang) };
    return R;
  }


  /**   \brief Creates and returns a rotation matrix around y axis
   * Returns a rotation matrix. The rotation is by ang against the
   * clock around the y-Axis
   */
  static Matrix<double> get_y_rot_m( double ang ){
    Matrix<double> R;
    R.v = {  cos(ang), 0, sin(ang) };
    R.w = {     0    , 1,    0     };
    R.z = { -sin(ang), 0, cos(ang) };
    return R;
  }


  /** \brief Creates and returns a rotation matrix around z axis
   * Returns a rotation matrix. The rotation is by ang against the
   * clock around the z-Axis
   */
  static Matrix<double> get_z_rot_m( double ang ){
    Matrix<double> R;
    R.v = { cos(ang), -sin(ang), 0 };
    R.w = { sin(ang),  cos(ang), 0 };
    R.z = { 0,         0,        1 };  
    return R;
  }
  
  /// Creates and returns a rotation matrix around an arbitrary axis
  /// from bronstein
  template<class T>  
  static Matrix<double> get_rot_n(std::vector<T> n, double ang ){
    
    double cc = 1-cos(ang);
    double ss = sin(ang);
    
    // row 1
    double R_11 = cos(ang) + n[0]*n[0]*cc;
    double R_12 = n[0]*n[1]*cc - n[2]*ss;
    double R_13 = n[0]*n[2]*cc + n[1]*ss;
    
    // row 2
    double R_21 = n[1]*n[0]*cc + n[2]*ss;
    double R_22 = cos(ang) + n[1]*n[1]*cc;
    double R_23 = n[1]*n[2]*cc-n[0]*ss;
    
    // row 3
    double R_31 = n[2]*n[0]*cc - n[1]*ss;
    double R_32 = n[2]*n[1]*cc + n[0]*ss;
    double R_33 = cos(ang) + n[2]*n[2]*cc;
    
    Matrix<double> R;
    R.v = { R_11, R_12, R_13 };
    R.w = { R_21, R_22, R_23 };
    R.z = { R_31, R_32, R_33 };
    return R;
  }
  
  /// returns two vectors, creating a orthogonal base
  static std::vector< std::vector<double> > get_orthogonal_base( std::vector<double> n ){
    
    auto unit_n = get_unit( n );
    
    double precision = 1e-8; // arbitrary choice
    
    // find the component which is non zero
    int non_zero;
    
    for(unsigned int ii=0; ii<3; ii++){
      if( std::fabs( unit_n[ii] ) > precision ){
	non_zero = ii;
	break;
      }
    }
    
    std::vector<double> a (3, 0), b(3, 0);
    
    if( non_zero == 0 ){
      a[0] = -unit_n[1]/unit_n[0]; a[1] = 1;
    } else if( non_zero == 1 ){
      a[0] = 1; a[1] = -unit_n[0]/unit_n[1];
    } else {
      a[0] = 1; a[2] = -unit_n[0]/unit_n[2];
    }
    
    // get a third orthogonal vector
    a = get_unit(a);
    b = cross_prod( unit_n, a );
    b = get_unit(b);
    
    std::vector< std::vector<double> > ret = { a, b };
    
    return ret;
  }
  
  /// returns unit length vector
  template<class T>    
  static std::vector<double> get_unit( std::vector<T> v ){
    double ll = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    std::vector<double> ret = { v[0]/ll, v[1]/ll, v[2]/ll };
    return ret;
  }

  /// check if the two vectors are linear dependent
  template<class T>  
  static bool check_lin_dep( std::vector<T> v, std::vector<T> w ){
    
    double scalar_prod = dot_prod( v, w );
    
    if( scalar_prod == 0 ){
      return false;
    }
    
    double l_v = sqrt( dot_prod(v, v) );
    double l_w = sqrt( dot_prod(w, w) );
    
    if( scalar_prod < l_v * l_w ){
      return false;
    } else if ( scalar_prod < l_v * l_w == 0 ){
      return true;
    }
    
  }
  
  
};


#endif
