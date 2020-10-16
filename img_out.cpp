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
