#include "slice_orientation_visualisation.hpp"

std::string draw_slice_orientation( int h, int k, int l, double theta, double phi ){

  // the orientation of the projection plane
  double phi_proj = 5.061455, theta_proj = 3.49066;

  sketch main_sketch;

  // the edgelength of the cube
  double a = 100;
    
  // draw a cube  
  point c1 {0,0,0}, c2{0, a, 0}, c3{a, a, 0}, c4{a,0,0};
  point c5 {0,0,a}, c6{0, a, a}, c7{a, a, a}, c8{a,0,a};  
    
  // cube
  main_sketch.add_polygon( polygon( std::vector<point> {c1, c4, c8, c5}, "black", "grey", 2, .3 ) );
  main_sketch.add_polygon( polygon( std::vector<point> {c1, c2, c6, c5}, "black", "grey", 2, .3 ) );
  main_sketch.add_polygon( polygon( std::vector<point> {c4, c3, c7, c8}, "black", "grey", 2, .3 ) );
  main_sketch.add_polygon( polygon( std::vector<point> {c2, c3, c7, c6}, "black", "grey", 2, .3 ) );
  main_sketch.add_polygon( polygon( std::vector<point> {c1, c2, c3, c4}, "black", "grey", 2, .3 ) );
  main_sketch.add_polygon( polygon( std::vector<point> {c5, c6, c7, c8}, "black", "grey", 2, .3 ) );  
  // highlight coordinate axes    
  main_sketch.add_line( line( point {0,0,0}, point{a,0,0}, "red", 5 ) );
  main_sketch.add_line( line( point {0,0,0}, point{0,a,0}, "green", 5 ) );
  main_sketch.add_line( line( point {0,0,0}, point{0,0,a}, "yellow", 5 ) );
  
  // cuttin plane
  point n {h, k, l};
  auto base = VEC_MAT_MATH::get_orthogonal_base( n );

  // initialize the plane, centered about origin, so there won't be a
  // translation when rotating
  point r1{-a/2.0, -a/2.0, 0}, r2{-a/2.0, a/2.0, 0}, r3{a/2.0, a/2.0, 0}, r4{a/2.0, -a/2.0, 0};

  // rotate the initial plane to match the plane indicated by miller indeces
  Matrix R1 = VEC_MAT_MATH::get_y_rot_m(theta), R2 = VEC_MAT_MATH::get_z_rot_m(phi);
  r1 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r1) );
  r2 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r2) );
  r3 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r3) );
  r4 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r4) );

  // offset plane to be centered around midpoint of cube
  r1[0]+=a/2.0; r1[1]+=a/2.0; r1[2]+=a/2.0;
  r2[0]+=a/2.0; r2[1]+=a/2.0; r2[2]+=a/2.0;
  r3[0]+=a/2.0; r3[1]+=a/2.0; r3[2]+=a/2.0;
  r4[0]+=a/2.0; r4[1]+=a/2.0; r4[2]+=a/2.0;        

  // draw the plane
  main_sketch.add_polygon( polygon( std::vector<point> {r1, r2, r3, r4}, "green", "green", 2, 0.3 ) );
    
  Matrix Ry_proj = VEC_MAT_MATH::get_y_rot_m(theta_proj), Rz_proj = VEC_MAT_MATH::get_z_rot_m(phi_proj);
  point n_proj {1,0,0};
  n_proj = VEC_MAT_MATH::dot_prod( Rz_proj, VEC_MAT_MATH::dot_prod( Ry_proj, n_proj) );
  
  canvas can (n_proj, 250, 250, 75, 150);
  can.project_sketch( main_sketch );
  return can.draw_projection_svg();
  
}
