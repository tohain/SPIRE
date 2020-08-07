#include "img_out.hpp"



void write_image( std::string fn, unsigned char *data, int width, int height, bool invert ){

  //find max value
  unsigned char min=std::numeric_limits<char>::max(), max=std::numeric_limits<char>::min();
  for( unsigned int ii=0; ii<width*height; ii++){

    if( data[ii] < min ){
      min = data[ii];
    }

    if( data[ii] > max ){
      max = data[ii];
    }
    
  }
  
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
      
      //scale
      int u_scaled = static_cast<int>((data[ind] - min)/(max - min)*255);
      //invert if wanted
      if( invert )
	u_scaled = 255 - u_scaled;

      out << u_scaled << " ";
    }
    out << std::endl;
  }  

  out.close();
  
}
