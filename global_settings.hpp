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


/*
 * this file just holds some globally used constants
 */

#ifndef GLOBAL_SETTINGS_H
#define GLOBAL_SETTINGS_H

#include <string>
#include <fstream>
#include "auxiliary.hpp"

class global_settings {

public:

  global_settings ( std::string filename ){
    
    std::ifstream read_config ( filename );
    
    if( read_config.good() && !read_config.fail() ){

      // could read the file
      std::string line;

      
      while( !read_config.eof() ){

	std::getline( read_config, line );

	my_utility::strip( line, ' ' );
	my_utility::strip( line, '=' );	
	
	if( line != "" && line[0] != '#' ){
	  
	  auto pair = my_utility::str_split( line, '=' );

	  if( pair.size() == 2 ){
	    my_utility::strip( pair[0] );
	    my_utility::strip( pair[1] );	    

 	    if( pair[0] == "g_color_b1" ){
	      g_color_b1 = pair[1];
	    }

	    if( pair[0] == "g_color_b2" ){
	      g_color_b2 = pair[1];
	    }

	    if( pair[0] == "g_color_n" ){
	      g_color_n = pair[1];
	    }

	    if( pair[0] == "res_measure_vol" ){
	      res_measure_vol = std::stoi( pair[1] );
	    }
	    if( pair[0] == "res_measure_area" ){
	      res_measure_area = std::stoi( pair[1] );
	    }
	    if( pair[0] == "res_measure_perc" ){
	      res_measure_perc = std::stoi( pair[1] );
	    }
	    if( pair[0] == "level_min_width" ){
	      level_min_width = std::stod( pair[1] );
	    }
	  }
	}	
      }      
    }    
  }

  void print(){

    std::cout << "g_color_b1=" << g_color_b1 << std::endl;
    std::cout << "g_color_b2=" << g_color_b2 << std::endl;
    std::cout << "g_color_n=" << g_color_n << std::endl;    
    std::cout << "res_measure_vol=" << res_measure_vol << std::endl;
    std::cout << "res_measure_area=" << res_measure_area << std::endl;
    std::cout << "res_measure_perc=" << res_measure_perc << std::endl;    
  }

  //
  // the thickness in "level-set units" of the initial main membrane
  // at f(x,y,z)=0
  //
  double level_min_width = 0.05;

  
  //
  // the color of the base vectors
  //
  std::string g_color_b1 = "#5aa02c";
  std::string g_color_b2 = "#ffcc00";
  std::string g_color_n  = "#c83737";

  //
  // the resolution used to measure the volume of the channels.
  // -1 means the current resolution from the GUI is used
  //
  int res_measure_vol = 76;


  //
  // the resolution used to measure the volume of the channels.
  // -1 means the current resolution from the GUI is used
  //
  int res_measure_area = 76;
  
  //
  // the resolution used to measure the percolation thresdhold
  // -1 means the current resolution from the GUI is used
  //
  int res_measure_perc = 76;    

};




#endif
