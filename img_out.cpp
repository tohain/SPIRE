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

#include "img_out.hpp"

void write_image( std::string fn, unsigned char *data, int width, int height, bool invert ){
  
  //open filestream
  std::ofstream out ( fn );

  //identifier pgm
  out << "P2";

  //width
  out << " " << width;

  //height
  out << " " << height;

  //max value
  out << " " << 255 << std::endl;

  //write image data
  for(unsigned int ii=0; ii<height; ii++){ //rows  = height (vertical)
    for(unsigned int jj=0; jj<width; jj++){ //cols = width  (horizontal)

      int ind = ii*width + jj;
      
      
      out << int(data[ind]) << " ";
    }
    out << std::endl;
  }  

  out.close();
  
}
