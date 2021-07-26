#include "img_manip.hpp"



image_manipulation::image_manipulation(){

}


double symm_gaussian_2D( double x, double y, double sigma ){
  
  return (1.0/( 2 * M_PI * sigma * sigma ) ) * exp( -(x*x + y*y) / (2*sigma*sigma)   );
  
}


void image_manipulation::gaussian_blur( unsigned char* img, unsigned int width, unsigned height, unsigned int kernel_size ){

  // create a new intermediate image

  unsigned char* soft_img = new unsigned char[width*height]();

  
  // first compute the convolution function (gaussian)

  // choose a 3sigma cutoff
  int conv_size = 6*kernel_size;
  if( conv_size % 2 == 0 ) conv_size++;
  int cv2 = int(conv_size*0.5); // half convolution size

  double* conv_func = new double[conv_size*conv_size]();

  double sum = 0;
  for( int ii=0; ii < conv_size; ii++ ){
    for( int jj=0; jj < conv_size; jj++ ){
      conv_func[ii+conv_size*jj] = symm_gaussian_2D( (int(0.5*conv_size) - ii), (int(0.5*conv_size) - jj), kernel_size );
      sum+=conv_func[ii+conv_size*jj];
    }
  }

  // normalize
  for( int ii=0; ii < conv_size*conv_size; ii++ ){
    conv_func[ii] /= sum;
  }
  

  // apply the gaussian filter in the width dim first

  // iterate over all pixles
  for( unsigned int yy = 0; yy < height; yy++){
    for( unsigned int xx = 0; xx < width; xx++){

      double pix_val = 0;      

      // do the convolution
      for( unsigned cx = 0; cx < conv_size; cx++){
	for( unsigned cy = 0; cy < conv_size; cy++){

	  int ind_x = xx - cv2 + cx;
	  int ind_y = yy - cv2 + cy;
	  // pbc
	  if( ind_x < 0 )
	    ind_x += width;
	  if( ind_x >= width )
	    ind_x -= width;

	  if( ind_y < 0 )
	    ind_y += height;
	  if( ind_y >= height )
	    ind_y -+ height;

	  // total index in 1d array
	  int ind = ind_y * width + ind_x;
	  
	  pix_val += conv_func[cx*conv_size + cy] * img[ind];
	}
      }
      soft_img[ yy*width + xx ] = (unsigned char ) pix_val;      
    }
  }

  // blend images, for now it's just replaced, could get fancier in future
  for( unsigned int i=0; i < width*height; i++){
    img[i] = soft_img[i];
  }

  delete[] (soft_img);
}


void image_manipulation::gaussian_noise( unsigned char* img, unsigned int width, unsigned int height, double magnitude ){
  
  // bring in the random numbers
  std::mt19937 generator;
  std::normal_distribution<double> dist (0, 1.0/6.0);

  // create a new image
  short *noisy_img = new short[width*height]();
  
  for( unsigned int ii=0; ii< width*height; ii++){
    noisy_img[ii] = img[ii] + static_cast<short>( dist(generator) * magnitude);
  }

  //rescale
  short min=0, max=0;
  for( unsigned int i=0; i < width*height; i++){
    if( noisy_img[i] > max )
      max = noisy_img[i];
    if( noisy_img[i] < min )
      min = noisy_img[i];
  }

  for( unsigned int i=0; i < width*height; i++){
    img[i] = (unsigned char) (((float)noisy_img[i] - min) / (max - min) * 255.0);    
  }

  delete[] (noisy_img);
}

  


