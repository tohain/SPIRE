/*
 * This file is part of the Surface Projection Tool
 *
 * Author: T. Hain
 * 07/2019
 *
 * Licence: TBA
 */


#include "distance_transform.hpp"

/**
 * Allocates necessary memory for the transformed map and initializes
 * class parameters
 */
distance_transform::distance_transform( double *data, int n ){
  //initialize parameters
  img = data;
  dim = 1;  
  size[0]=n; size[1]=-1; size[2]=-1;
  
  //allocate memory
  map = std::vector<double> (n, 0);
}


/**
 * Allocates necessary memory for the transformed map and initializes
 * class parameters. Image must be provided as a concatenated list of
 * rows
 *
 * \param[in] data The array holding the image data. 
 * \param[in] n the number of columns in the image
 * \param[in] k the number of rows in the image
 */
distance_transform::distance_transform( double *data, int n, int k ){
  //initialize parameters
  img = data;
  dim = 2;
  size[0]=n; size[1]=k; size[2]=-1;
  
  //allocate memory
  map = std::vector<double> (n*k, 0);
}


/**
 * Allocates necessary memory for the transformed map and initializes
 * class parameters
 *
 * \param[in] data The array holding the image data
 * \param[in] n the number of columns in the image
 * \param[in] k the number of rows in the image
 * \param[in] l the number of layers(?) in the image
 */
distance_transform::distance_transform( double *data, int n, int k, int l ){
  //initialize parameters
  img = data;
  dim = 3;
  size[0]=n; size[1]=k; size[2]=l;
  
  //allocate memory
  map = std::vector<double> (n*k*l, 0);
}


/**
 * Destructor, just frees memory
 */
distance_transform::~distance_transform(){
}



/**
 * This function performs the actual 1D distance transform. grid is
 * the function evaluated on the grid. In case of the "standard" EDM
 * it is infinity for background pixels and 0 for object pixels.
 *
 * This algorithm looks for the lowest envelope of parabolas situated
 * at poins (i, grid[i])
 */
std::vector<double> distance_transform::do_distance_transform( std::vector<double> &grid ){

  std::vector<double> transform ( grid.size(), 0 );
  
  //Index of rightmost parabola in lower envelope
  unsigned int k = 0;
  //Locations of parabolas in lower envelope
  std::vector<int> v = {0}; 
  //Locations of boundaries between parabolas
  std::vector<double> z = {-std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};

  //find the envelope
  for( unsigned int ii=1; ii<grid.size(); ii++){
    //intersection of two parabolas
    double s = ( ( grid[ii] + (ii*ii) ) - (grid[ v[k] ] + (v[k]*v[k])) ) / (2*(ii-v[k]));

    //check where the intersection is. If it is left of the previous
    //parabola, the current one is the new lowest envelope
    while( s <= z[k] ){
      k--;
      s = ( ( grid[ii] + (ii*ii) ) - (grid[ v[k] ] + (v[k]*v[k])) ) / (2*(ii-v[k]));
    }

    //found the previous lowest envelope
    k++;
    if( k < v.size() ){
      v.at(k) = ii;
    } else {
      v.push_back( ii );
    }

    if( k < z.size() ){
      z[k] = s;
    } else {
      z.push_back( s );
    }

    if( k+1 < z.size() ){
      z[k+1] = std::numeric_limits<double>::max();
    } else {
      z.push_back( std::numeric_limits<double>::max() );
    }    
  }


  //fill in the transform by evaluating the parabola envelope at the
  //grid
  k=0;
  for(unsigned int ii=0; ii<grid.size(); ii++){
    while( z[k+1] < ii ){
      k++;
    }
    transform[ii] = ( (ii-v[k])*(ii-v[k]) ) + grid[ v[k] ];
  }

  return transform;
}


/**
 * computes the cost function of the given image data and writes it to grid
 */
std::vector<double> distance_transform::eval_grid_function( std::vector<double> &data ){

  //initialize memory and set it to zero
  std::vector<double> grid (data.size(), 0);

  for(unsigned int ii=0; ii<data.size(); ii++){
    if( data[ii] == 0 ){
      grid[ii] = std::numeric_limits<double>::max();
    } else {
     grid[ii] = 0;
    }
  }

  return grid;
}



/**
 * This method computes the distance map. Calls the base function
 * couple of times, depending on the dimension. The result is stored
 * in \ref map
 */
void distance_transform::compute_distance_map(){

  if( dim == 1 ){
    //one dimensional case

    //cpy img data into a vector
    std::vector<double> data (size[0], 0);
    memcpy( data.data(), img, size[0] * sizeof(double ) );
    
    //get grid
    std::vector<double> grid = eval_grid_function( data );

    //do the distance transform
    map = do_distance_transform( grid );
  }

  else if ( dim == 2 ){
    //2d case

    std::vector<double> tmp_map ( map.size(), 0 );
    
    //do the distance transform for each column first
    for( unsigned int ii=0; ii < size[0]; ii++ ){

      //get the ii-th column
      std::vector<double> col (size[1], 0);
      for( unsigned int jj=0; jj<size[1]; jj++){
	col.at(jj) = img[ii + jj * size[0]];
      }

      //get the grid
      std::vector<double> col_grid = eval_grid_function( col );

      //get distance map of column
      std::vector<double> col_map = do_distance_transform( col_grid );

      //cpy preliminary map in to map vector
      memcpy( tmp_map.data()+(ii*size[1]), col_map.data(), sizeof(double) * col_map.size() );
      
    }

    /*
     * careful, map is now not row based, but column based. We need to
     * reverse that in the following step
     */

    //now iterate over all rows and compute the distance map row wise
    for( unsigned int ii=0; ii<size[1]; ii++){

      //get the iith row
      std::vector<double> row (size[0], 0);
      for( unsigned int jj=0; jj<size[0]; jj++){
	row[jj] = tmp_map[ ii + jj*size[1] ];
      }

      //now get the transform
      std::vector<double> row_map = do_distance_transform( row );

      //cpy the final map data to map object
      memcpy( map.data()+(ii*size[0]), row_map.data(), sizeof(double) * row_map.size() );
    }
       

  }
  

}











/*
 * getters
 */
std::vector<double> distance_transform::get_distance_map() const {
  return map;
}
