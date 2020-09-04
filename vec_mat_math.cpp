#include "vec_mat_math.hpp"


/**
 * Quick and dirty dot product between two vectors
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

