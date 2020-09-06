

#ifndef DRAWME_HPP
#define DRAWME_HPP

#include <string>
#include <iostream>
#include "vec_mat_math.hpp"

#include "svglib.hpp"

typedef std::vector<double> point;


/**
 * \brief A simple line
 */
class line {

public:
  line( point start, point end, std::string color = "black", double width = 1 );
  
  point start;
  point end;

  std::string color;
  double width;  
};


class polygon {

public:

  polygon( std::vector<point> points, std::string stroke_color = "black", std::string fill_color = "transparent", double stroke_width = 1, double alpha = 1 );

  std::vector<point> points;
  std::string stroke_color;
  std::string fill_color;
  double stroke_width;
  double alpha;
  
};

class sketch {

public:

  sketch();

  void add_line( line l );
  void add_polygon( polygon p );
  
  std::vector<line> get_lines();
  std::vector<polygon> get_polys();
private:


  std::vector<line> lines;
  std::vector<polygon> polys;
  
};


class canvas {
public:

  canvas( std::vector<double> n, int res_x, int res_y  );
  canvas( std::vector<double> n, int res_x, int res_y, double mid_x, double mid_y  );  

  void project_sketch( sketch &sk );

  std::string draw_projection_svg(  );
  
  void print_lines();
  
private:

  // objects to draw
  std::vector<line> lines;
  std::vector<polygon> polys;

  
  std::vector<double> normal;

  void update_projection_matrix();
  Matrix projection_matrix;
  
  std::vector<int> pixvals;

  // the resolution of the output pixmap
  int res_x, res_y;

  // the point in coordinates to put in the middle
  double mid_x, mid_y;
  
  // the corners of the image
  double lo_x, lo_y;
  double hi_x, hi_y;
  
};


#endif
