#include "img_out.hpp"



void write_image( std::string fn, std::vector< std::vector<double> > data ){

  //find max value
  double max=std::numeric_limits<double>::min(),min=std::numeric_limits<double>::max();
  for(unsigned int ii=0; ii<data.size(); ii++){
    for(unsigned int jj=0; jj<data[ii].size(); jj++){
      if( data[ii][jj] > max )
	max = data[ii][jj];
      if( data[ii][jj] < min )
	min = data[ii][jj];
    }
  }

  //open filestream
  std::ofstream out ( fn );

  //identifier pgm
  out << "P2";

  //width
  out << " " << data[0].size();

  //height
  out << " " << data.size();

  //max value
  out << " " << 255 << std::endl;

  //write image data
  for(unsigned int ii=0; ii<data.size(); ii++){
    for(unsigned int jj=0; jj<data[ii].size(); jj++){
      //scale
      int u_scaled = static_cast<int>((data[ii][jj] - min)/(max - min)*255);
      out << 255-u_scaled << " ";
    }
    out << std::endl;
  }  

  out.close();
  
}
