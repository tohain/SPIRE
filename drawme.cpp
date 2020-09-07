#include "drawme.hpp"

/**
 * class sketch
 */ 

sketch::sketch() {
}


void sketch::add_line( line l ){
  lines.push_back( l );
}

void sketch::add_polygon( polygon p ){
  polys.push_back( p );
}

std::vector<line> sketch::get_lines(){
  return lines;
}

std::vector<polygon> sketch::get_polys(){
  return polys;
}



/**
 * class line
 */ 

line::line( point start_, point end_, std::string color_, double width_ ) : start(start_), end(end_), color(color_), width(width_) {

}




/**
 * class polygon
 */

polygon::polygon( std::vector<point> points_, std::string stroke_color_, std::string fill_color_, double stroke_width_, double alpha_ ) : points( points_), stroke_color( stroke_color_ ), fill_color( fill_color_ ), stroke_width( stroke_width_), alpha( alpha_) {

}



/**
 * class canvas
 */

canvas::canvas( std::vector<double> n, int res_x_, int res_y_, point mid_ ) : normal(n), res_x(res_x_), res_y(res_y_), mid_point( mid_ ){

  size_x = res_x;
  size_y = res_y;
  
  pixvals.resize( res_x * res_y, 0 );
  
}

canvas::canvas( std::vector<double> n, int res_x_, int res_y_ ) : normal(n), res_x(res_x_), res_y(res_y_){

  size_x = res_x;
  size_y = res_y;
  
  mid_point = point { 0, 0, 0};
  
  pixvals.resize( res_x * res_y, 0 );
  
}

canvas::canvas( std::vector<double> n, int res_x_, int res_y_, double size_x_, double size_y_, point mid_ ) : normal(n), res_x(res_x_), res_y(res_y_), mid_point( mid_ ), size_x( size_x_), size_y( size_y_ ){

  pixvals.resize( res_x * res_y, 0 );  
}

canvas::canvas( std::vector<double> n, int res_x_, int res_y_, double size_x_, double size_y_ ) : normal(n), res_x(res_x_), res_y(res_y_), size_x( size_x_ ), size_y( size_y_ ){
  mid_point = point { 0, 0, 0};
  
  pixvals.resize( res_x * res_y, 0 );
  
}


/**
 * This function updates the projection matrix from the normal vector
 * of the projection plane
 */
void canvas::update_projection_matrix(){

  auto base = VEC_MAT_MATH::get_orthogonal_base( normal );
  projection_matrix = Matrix(base[0], base[1], point {0,0,0});
}


/**
 * This function projects the sketch onto the canvas
 */
void canvas::project_sketch( sketch &sk ){

  // clear former sketch
  lines.clear();
  polys.clear();
  
  // get the projection matrix
  update_projection_matrix();

  // get the projection of the midpoint
  mid_point_proj = VEC_MAT_MATH::dot_prod( projection_matrix, mid_point );
  
  // get lines and polygons to project
  std::vector<line> dd_lines = sk.get_lines();
  std::vector<polygon> dd_polys = sk.get_polys();

  // project each line
  for( auto ll : dd_lines ){
    
    std::vector<double> start_flat = VEC_MAT_MATH::dot_prod( projection_matrix, ll.start );
    std::vector<double> end_flat = VEC_MAT_MATH::dot_prod( projection_matrix, ll.end );

    // copy objects so color and stroke will be copied
    line flat_line = ll;
    // updated with projected coordinates
    flat_line.start = start_flat;flat_line.end = end_flat;

    lines.push_back(flat_line);
  }

  // same as for the lines
  for( auto poly : dd_polys ){
    
    std::vector<point> flat_points;
    for( unsigned int ii = 0; ii < poly.points.size(); ii++){
      flat_points.push_back( VEC_MAT_MATH::dot_prod( projection_matrix, poly.points[ii] ) );     
    }

    polygon flat_pp = poly;
    flat_pp.points = flat_points;

    polys.push_back( flat_pp );
       
  }
  
}


/**
 * Outputs the projection to an svg image, using my own svg library
 */
std::string canvas::draw_projection_svg(){

  svg_canvas can( res_x, res_y );

  
  double scale_x = (double) res_x / size_x;
  double scale_y = (double) res_y / size_y;

  
  for( auto pp : polys ){

    svg_polygon svg_pp ( pp.stroke_color, pp.stroke_width, pp.fill_color, pp.alpha );

    for( unsigned int ii = 0; ii < pp.points.size(); ii++){
      svg_pp.add_point( (scale_x * (-mid_point_proj[0] + pp.points[ii][0])) + res_x/2.0,
			(scale_y * (-mid_point_proj[1] + pp.points[ii][1])) + res_y/2.0 );
    }

    can.add_object( &svg_pp );
  }

  for( auto ll : lines ){

    svg_line tmp ( (scale_x * (-mid_point_proj[0] + ll.start[0])) + res_x/2.0,
		   scale_y * (-mid_point_proj[1] + ll.start[1]) + res_y/2.0,
		   scale_x * (-mid_point_proj[0] + ll.end[0]) + res_x/2.0,
		   scale_y * (-mid_point_proj[1] + ll.end[1]) + res_y/2.0,
		   ll.color, ll.width  );
    
    can.add_object( &tmp );
  }
  
  return can.to_string();
  
}


void canvas::print_lines(){
  for(auto it : lines ){

    std::cout << "(" << it.start[0] << " " << it.start[1] << " " << it.start[2] << ") --> ("
	      << it.end[0] << " " << it.end[1] << " " << it.end[2] << ")" << std::endl;  
  }
}
