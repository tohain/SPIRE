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

std::string draw_slice_orientation( int h, int k, int l, double theta, double phi, double Lx, double Ly, double Lz ){


  double max_L = std::max( Lx, std::max( Ly, Lz ) );
  Lx /= max_L;
  Ly /= max_L;
  Lz /= max_L;

  double min_L = std::min( Lx, std::min( Ly, Lz ) );
  
  // the orientation of the projection plane
  double phi_proj = 5.061455, theta_proj = 3.49066;

  sketch main_sketch;

  // the edgelength of the cube
  double a = 50;

  // cuttin plane
  point n {static_cast<double>(h), static_cast<double>(k), static_cast<double>(l)};
  auto base = VEC_MAT_MATH::get_orthogonal_base( n );
  auto unit_normal = VEC_MAT_MATH::get_unit( n );

  
  // initialize the plane, centered about origin
  point r1{-min_L*a, -min_L*a, 0}, r2{-min_L*a, min_L*a, 0}, r3{min_L*a, min_L*a, 0}, r4{min_L*a, -min_L*a, 0};

  // rotate the initial plane to match the plane indicated by miller indeces
  Matrix R1 = VEC_MAT_MATH::get_y_rot_m(theta), R2 = VEC_MAT_MATH::get_z_rot_m(phi);
  r1 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r1) );
  r2 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r2) );
  r3 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r3) );
  r4 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r4) );

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
  Matrix Ry_proj = VEC_MAT_MATH::get_y_rot_m(theta_proj), Rz_proj = VEC_MAT_MATH::get_z_rot_m(phi_proj);
  point n_proj {1,0,0};
  n_proj = VEC_MAT_MATH::dot_prod( Rz_proj, VEC_MAT_MATH::dot_prod( Ry_proj, n_proj) );
  
  // draw the normal vector
  main_sketch.add_line( line( point {0,0,0}, point {a*unit_normal[0], a*unit_normal[1], a*unit_normal[2] }, "red", 2 ) );
  
  // create canvas and get 2D projection
  canvas can (n_proj, 75, 75, 150, 150, point {0,0,0});
  can.project_sketch( main_sketch );

  return can.draw_projection_svg();
  
}
