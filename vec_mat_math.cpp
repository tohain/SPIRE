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


#include "vec_mat_math.hpp"


Matrix::Matrix( std::vector<double> x, std::vector<double> y , std::vector<double> z ) : v(x), w(y), z(z) {

}


void Matrix::print(){

  std::cout << "/" << v[0] << " " << v[1] << " " << v[2] << "\\" << std::endl;
  std::cout << "|" << w[0] << " " << w[1] << " " << w[2] << "|" << std::endl;
  std::cout << "\\" << z[0] << " " << z[1] << " " << z[2] << "/" << std::endl;

}

/**
 * Quick and dirty scalar product between scalar and vector
 */
std::vector<double> VEC_MAT_MATH::s_prod ( double s, std::vector<double> w ){
  return std::vector<double> { s*w[0], s*w[1], s*w[2] };  
}

/**
 * Quick and dirty sum of two vectors
 */
std::vector<double> VEC_MAT_MATH::vec_add ( std::vector<double> v, std::vector<double> w ){
  return std::vector<double> { v[0]+w[0], v[1]+w[1], v[2]+w[2] };  
}

/**
 *  Quick and dirty dot product between two vectors
 */
double VEC_MAT_MATH::dot_prod ( std::vector<double> v, std::vector<double> w ){

  return v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
}


/**
 * Quick and dirty matrix-vector product
 */
std::vector<double> VEC_MAT_MATH::dot_prod( Matrix m, std::vector<double> v ){

  std::vector<double> r (3, 0);
  r[0] = dot_prod( m.v, v);
  r[1] = dot_prod( m.w, v);
  r[2] = dot_prod( m.z, v);
  return r;
}


std::vector<double> VEC_MAT_MATH::dot_prod( std::vector<double> v, Matrix m ){

  std::vector<double> n_col_1 = {m.v[0], m.w[0], m.z[0]};
  std::vector<double> n_col_2 = {m.v[1], m.w[1], m.z[1]};
  std::vector<double> n_col_3 = {m.v[2], m.w[2], m.z[2]};  
  
  std::vector<double> r (3, 0);
  r[0] = dot_prod( v, n_col_1 );
  r[1] = dot_prod( v, n_col_2 );
  r[2] = dot_prod( v, n_col_3);
  return r;
}

Matrix VEC_MAT_MATH::dot_prod( Matrix m, Matrix n ){

  std::vector<double> n_col_1 = {n.v[0], n.w[0], n.z[0]};
  std::vector<double> n_col_2 = {n.v[1], n.w[1], n.z[1]};
  std::vector<double> n_col_3 = {n.v[2], n.w[2], n.z[2]};  

  std::vector<double> prod_row_1 = { dot_prod( m.v, n_col_1 ), dot_prod( m.v, n_col_2 ), dot_prod( m.v, n_col_3 ) };
  std::vector<double> prod_row_2 = { dot_prod( m.w, n_col_1 ), dot_prod( m.w, n_col_2 ), dot_prod( m.w, n_col_3 ) };
  std::vector<double> prod_row_3 = { dot_prod( m.z, n_col_1 ), dot_prod( m.z, n_col_2 ), dot_prod( m.z, n_col_3 ) };
  
  
  return Matrix( prod_row_1, prod_row_2, prod_row_3 );
}

/**
 * Returns a rotation matrix. The rotation is by ang against the
 * clock around the x-Axis
 */
Matrix VEC_MAT_MATH::get_x_rot_m (double ang) {
  Matrix R;
  R.v = { 1, 0,           0          };
  R.w = { 0, cos(ang), -sin(ang) };
  R.z = { 0, sin(ang),  cos(ang) };
  return R;
}


/**
 * Returns a rotation matrix. The rotation is by ang against the
 * clock around the y-Axis
 */
Matrix VEC_MAT_MATH::get_y_rot_m (double ang) {
  Matrix R;
  R.v = {  cos(ang), 0, sin(ang) };
  R.w = {     0    , 1,    0     };
  R.z = { -sin(ang), 0, cos(ang) };
  return R;
}

/**
 * Returns a rotation matrix. The rotation is by ang against the
 * clock around the z-Axis
 */
Matrix VEC_MAT_MATH::get_z_rot_m (double ang) {
  Matrix R;
  R.v = { cos(ang), -sin(ang), 0 };
  R.w = { sin(ang),  cos(ang), 0 };
  R.z = { 0,         0,        1 };  
  return R;
}

std::vector<double> VEC_MAT_MATH::get_unit( std::vector<double> v ){

  double ll = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  std::vector<double> ret = { v[0]/ll, v[1]/ll, v[2]/ll };

  return ret;
}


std::vector<double> VEC_MAT_MATH::cross_prod( std::vector<double> a, std::vector<double> b ){

  std::vector<double> out (3,0);

  out[0] = a[1]*b[2] - a[2]*b[1];
  out[1] = a[2]*b[0] - a[0]*b[2];
  out[2] = a[0]*b[1] - a[1]*b[0];

  return out;
}


/**
 * returns two vectors which are perpendicular to the input vector
 */
std::vector< std::vector<double> > VEC_MAT_MATH::get_orthogonal_base( std::vector<double> n ){

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
