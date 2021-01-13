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
  canvas( std::vector<double> n, int res_x, int res_y, double size_x, double size_y  );
  canvas( std::vector<double> n, int res_x, int res_y, point mid_  );
  canvas( std::vector<double> n, int res_x, int res_y, double size_x, double size_y, point mid_  );

  void project_sketch( sketch &sk );

  std::string draw_projection_svg();
  
  void print_lines();
  
private:

  // objects to draw
  std::vector<line> lines;
  std::vector<polygon> polys;

  
  std::vector<double> normal;

  void update_projection_matrix();
  Matrix<double> projection_matrix;
  
  std::vector<int> pixvals;

  // the resolution of the output pixmap
  int res_x, res_y;
  int size_x, size_y;
  
  // the point in coordinates to put in the middle
  point mid_point, mid_point_proj;
    
};


#endif
