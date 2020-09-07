#include "slice_orientation_visualisation.hpp"

std::string draw_slice_orientation( int h, int k, int l, double theta, double phi ){

  // the orientation of the projection plane
  double phi_proj = 5.061455, theta_proj = 3.49066;

  sketch main_sketch;

  // the edgelength of the cube
  double a = 50;

  // cuttin plane
  point n {h, k, l};
  auto base = VEC_MAT_MATH::get_orthogonal_base( n );
  auto unit_normal = VEC_MAT_MATH::get_unit( n );

  
  // initialize the plane, centered about origin
  point r1{-a, -a, 0}, r2{-a, a, 0}, r3{a, a, 0}, r4{a, -a, 0};

  // rotate the initial plane to match the plane indicated by miller indeces
  Matrix R1 = VEC_MAT_MATH::get_y_rot_m(theta), R2 = VEC_MAT_MATH::get_z_rot_m(phi);
  r1 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r1) );
  r2 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r2) );
  r3 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r3) );
  r4 = VEC_MAT_MATH::dot_prod( R2, VEC_MAT_MATH::dot_prod( R1, r4) );

   // draw a cube  
  // cube centered around origin
  point c1 {-a,-a,-a}, c2{-a, a, -a}, c3{a, a, -a}, c4{a,-a,-a};
  point c5 {-a,-a,a}, c6{-a, a, a}, c7{a, a, a}, c8{a,-a,a};  


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
