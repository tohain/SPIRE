#ifndef IMG_MANIP_I
#define IMG_MANIP_I

#include <cmath>
#include <iostream>
#include <random>
#include <cstring>
#include <ctime>
#include <chrono>

#ifdef HAVE_PNG
#include <png.h>
#include <string>
#include <cstdio>
#endif

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

  /**
   * \brief adds gaussian noise to the image
   *
   * This function adds noise to an image by multiplying the magnitude
   * with normal-distributed random numbers of sigma=1/6 and mu=0.
   * The image is rescaled after the noise is added
   *
   * \param[in] img The image to add noise to
   * \param[in] width The width of the image in pixels
   * \param[in] height The height of the image in pixels
   * \param[in] magnitude The magnitude of the noise
   */
  static void gaussian_noise( unsigned char* img, unsigned int width, unsigned int height, double magnitude );



  /**
   * \brief Adds grains to the image
   *
   * Adds grains to the image, i.e. grainy gaussian noise.
   *
   * \param[in] img The image to add noise to
   * \param[in] width The width of the image in pixels
   * \param[in] height The height of the image in pixels
   * \param[in] grain_size_center The center of the normal distribution for grain_size
   * \param[in] grain_size_width The standard deviation of the normal distribution for grain_size
   * \param[in] grain_number_center The center of the normal distribution for grain_number
   * \param[in] grain_number_width The standard deviation of the normal distribution for grain_number
   * \param[in] magnitude The intensity of the grains
   */
  static void add_grains( unsigned char* img, unsigned int width, unsigned int height,
			  int grain_size_center, int grain_size_width,
			  int grain_number_center, int grain_number_width,
			  double magnitude );
  

#ifdef HAVE_PNG
  /**
   * \brief Reads an png image from a file into memory as a byte array.
   *
   * This function uses pnglib to read an image into a byte array
   * \param[in] file The png to open
   * \param[out] img The pointer to the byte array. Memory is allocated in the function.
   * \param[out] width The width of the image in pixels
   * \param[out] height The height of the image in pixels
   */
  static void read_png( std::string file, unsigned char** img, int *width, int *height );


  /**
   * \brief Writes a monochromatic png image from a byte array to a file
   *
   * This function uses pnglib to read an image into a byte array
   * \param[in] file The filename to write the image to
   * \param[in] img The pointer to the byte array. 
   * \param[in] width The width of the image in pixels
   * \param[in] height The height of the image in pixels
   */
  static void write_png( std::string file, unsigned char* img, unsigned int width, unsigned int height );  
#endif
  
};



#endif
