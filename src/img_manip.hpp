#ifndef IMG_MANIP_I
#define IMG_MANIP_I

#include <cmath>
#include <iostream>
#include <random>

class image_manipulation {

public:
  
  image_manipulation();



  /**
   * \brief applies a gaussian blur to the input image
   *
   * This function smoothens the image by applying a gaussian blur to
   * the image. This is achieved by convolving the image with a 2D
   * gaussian. Replace the image at the input pointer. 
   *
   *
   * \param[in] img the image to apply the filter to
   * \param[in] width The width of the image in pixels
   * \param[in] height The height of the image in pixels
   * \param[in] The standard deviation of the 2D gaussian in x and y direction
   * \param[in] magnitude the height of the 2D gaussian
   */
  static void gaussian_blur( unsigned char* img, unsigned int width, unsigned int height, unsigned int kernel_size );

  static void gaussian_noise( unsigned char* img, unsigned int width, unsigned int height, double magnitude );
  
};



#endif
