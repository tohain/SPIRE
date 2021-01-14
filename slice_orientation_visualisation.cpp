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


#include "slice_orientation_visualisation.hpp"

/**
 * draws a visualization of the slice and the orientation of the viewing angle inside it
 * \param[in] base A 1d array containing the three base vectors in column-major order. First vector is x direction, last one is normal onto viewing plane
 */
std::string draw_slice_visualization( std::vector<double> base, std::vector<double> L ){


  
  //get an orthonormal base
  std::vector<double> x = {base[0], base[1], base[2]};
  std::vector<double> y = {base[3], base[4], base[5]};
  std::vector<double> z = {base[6], base[7], base[8]};  

  x = VEC_MAT_MATH::get_unit( x );
  z = VEC_MAT_MATH::get_unit( z );
  std::vector<double> y_rect = VEC_MAT_MATH::cross_prod( z, x );

  double max_L = std::max( L[0], std::max( L[1], L[2] ) );
  double Lx = L[0] / max_L;
  double Ly = L[1] / max_L;
  double Lz = L[2] / max_L;

  double min_L = std::min( Lx, std::min( Ly, Lz ) );
  
  // the orientation of the projection plane
  double phi_proj = 5.061455, theta_proj = 3.49066;

  sketch main_sketch;

  // the edgelength of the cube
  double a = 50;

  // draw the hkl plane centered around origin
  point r1 = VEC_MAT_MATH::vec_add( VEC_MAT_MATH::s_prod(  min_L*a, x ),
				    VEC_MAT_MATH::s_prod(  min_L*a, y_rect ) );
  point r2 = VEC_MAT_MATH::vec_add( VEC_MAT_MATH::s_prod( -min_L*a, x ),
				    VEC_MAT_MATH::s_prod(  min_L*a, y_rect ) );
  point r3 = VEC_MAT_MATH::vec_add( VEC_MAT_MATH::s_prod( -min_L*a, x ),
				    VEC_MAT_MATH::s_prod( -min_L*a, y_rect ) );  
  point r4 = VEC_MAT_MATH::vec_add( VEC_MAT_MATH::s_prod(  min_L*a, x ),
				    VEC_MAT_MATH::s_prod( -min_L*a, y_rect ) );  
   // draw a cube  
  // cube centered around origin
  point c1 {-Lx*a,-Ly*a,-Lz*a}, c2{-Lx*a, Ly*a, -Lz*a}, c3{Lx*a, Ly*a, -Lz*a}, c4{Lx*a,-Ly*a,-Lz*a};
  point c5 {-Lx*a,-Ly*a,Lz*a}, c6{-Lx*a, Ly*a, Lz*a}, c7{Lx*a, Ly*a, Lz*a}, c8{Lx*a,-Ly*a,Lz*a};  

  // draw the background of the cube
  main_sketch.add_polygon( polygon( std::vector<point> {c1, c4, c8, c5}, "black", "grey", 2, .3 ) );
  main_sketch.add_polygon( polygon( std::vector<point> {c1, c2, c6, c5}, "black", "grey", 2, .3 ) );
  main_sketch.add_polygon( polygon( std::vector<point> {c4, c3, c7, c8}, "black", "grey", 2, .3 ) );
  main_sketch.add_polygon( polygon( std::vector<point> {c1, c2, c3, c4}, "black", "grey", 2, .3 ) );
    
  // draw the plane on top of the background
  main_sketch.add_polygon( polygon( std::vector<point> {r1, r2, r3, r4}, "blue", "blue", 2, 0.3 ) );

  // last, draw the foreground of the cube
  main_sketch.add_polygon( polygon( std::vector<point> {c2, c3, c7, c6}, "black", "grey", 2, .3 ) );
  main_sketch.add_polygon( polygon( std::vector<point> {c5, c6, c7, c8}, "black", "grey", 2, .3 ) );  
  
  // get the vector indicating the angle of view of the skecth
  Matrix<double> Ry_proj = VEC_MAT_MATH::get_y_rot_m(theta_proj), Rz_proj = VEC_MAT_MATH::get_z_rot_m(phi_proj);
  point n_proj {1,0,0};
  n_proj = VEC_MAT_MATH::dot_prod( Rz_proj, VEC_MAT_MATH::dot_prod( Ry_proj, n_proj) );
  
  // draw the normal vector
  main_sketch.add_line( line( point {0,0,0}, VEC_MAT_MATH::s_prod(a, z), "red", 3 ) );
  // draw the two in-plane vectors
  main_sketch.add_line( line( point {0,0,0}, VEC_MAT_MATH::s_prod(0.5*a, x), "blue", 1 ) );  
  main_sketch.add_line( line( point {0,0,0}, VEC_MAT_MATH::s_prod(0.5*a, VEC_MAT_MATH::get_unit(y)), "yellow", 1 ) );
  
  // create canvas and get 2D projection
  canvas can (n_proj, 75, 75, 150, 150, point {0,0,0});
  can.project_sketch( main_sketch );

  return can.draw_projection_svg();

  
}
