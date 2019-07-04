/*
 * This file is part of the Surface Projection Tool
 *
 * Author: T. Hain
 * 07/2019
 *
 * Licence: TBA
 */


#include "distance_transform.hpp"


// Start of with a few template instantiations
template class distance_transform<int>;
template class distance_transform<float>;
template class distance_transform<double>;


/**
 * Allocates necessary memory for the transformed map and initializes
 * class parameters
 *
 * \param[in] data The data to transform
 * \param[in] n The length of the data - not necessary, but for consistency
 * \param[in] pixsize The size of a pixel in length units
 */
template <class T>
distance_transform<T>::distance_transform( std::vector<T> data, int n, double pixsize ){
  //initialize parameters
  img = data;
  dim = 1;  
  size[0]=n; size[1]=-1; size[2]=-1;
  pix_size[0]=pixsize; pix_size[1]=1; pix_size[2]=1;
  
  //allocate memory
  map = std::vector<double> (n, 0);
}


/**
 * Allocates necessary memory for the transformed map and initializes
 * class parameters. The data must be provided as followed: with
 * increasing index, first the rows gets filled, then the columns
 *
 * \param[in] data The array holding the image data. 
 * \param[in] n the width of the grid (=nr of cols)
 * \param[in] k the height of the grid (=nr of rows)
 * \param[in] pixsize_x The size of a pixel in x(n) direction in length units
 * \param[in] pixsize_y The size of a pixel in y(k) direction in length units
 */
template <class T>
distance_transform<T>::distance_transform( std::vector<T> data, int n, int k, double pixsize_x, double pixsize_y ){
  //initialize parameters
  img = data;  
  dim = 2;
  size[0]=n; size[1]=k; size[2]=-1;
  pix_size[0]=pixsize_x; pix_size[1]=pixsize_y; pix_size[2]=1;

  //allocate memory
  map = std::vector<double> (n*k, 0);
}


/**
 * Allocates necessary memory for the transformed map and initializes
 * class parameters. The data must be provided as follwed: with
 * increasing index, first the stacks are filled up, then the rows are
 * filled with stacks, and finally the columns are filled with the
 * layers made of rows of stacks.
 *
 * \param[in] data The array holding the image data
 * \param[in] n the width of the grid (=nr of cols)
 * \param[in] k the height of the grid (=nr of rows)
 * \param[in] l the depth of the grid (=height of "stacks")
 * \param[in] pixsize_x The size of a pixel in x(n) direction in length units
 * \param[in] pixsize_y The size of a pixel in y(k) direction in length units
 * \param[in] pixsize_z The size of a pixel in z(l) direction in length units
 */
template <class T>
distance_transform<T>::distance_transform( std::vector<T> data, int n, int k, int l, double pixsize_x, double pixsize_y, double pixsize_z ){
  //initialize parameters
  img = data;
  dim = 3;
  size[0]=n; size[1]=k; size[2]=l;
  pix_size[0]=pixsize_x; pix_size[1]=pixsize_y; pix_size[2]=pixsize_z;
  
  //allocate memory
  map = std::vector<double> (n*k*l, 0);
}


/**
 * Destructor, just frees memory
 */
template <class T>
distance_transform<T>::~distance_transform(){
}



/**
 * This function performs the actual 1D distance transform. grid is
 * the function evaluated on the grid. In case of the "standard" EDM
 * it is infinity for background pixels and 0 for object pixels.
 *
 * This algorithm looks for the lowest envelope of parabolas situated
 * at poins (i, grid[i])
 *
 * \param[in] grid The 1D array to transform
 * \param[in] pixsize The size of a pixel in length units
 */
template <class T>
std::vector<double> distance_transform<T>::do_distance_transform( std::vector<double> &grid, double pixsize ){

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
      v[k] = ii;
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
    
    //get the "real" distances in length units by multiplying squared
    //index distances by pixelsize
    transform[ii] = pixsize * pixsize *  ( (ii-v[k])*(ii-v[k]) ) + grid[ v[k] ];
  }

  return transform;
}


/**
 * Applies a cost function on the given image/grid. The cost function
 * is chosen so an Euclidean Distance Map (EDM) is computed.
 *
 * This is a weird construction with two templates. However, this is
 * needed, since this function will be called with a vector of type of
 * the class parameter (int, float, ...) and also with double (see
 * do_distance_trasnform ). Therefore we need instanciations of this
 * function with int and double, thus we can't use T here in the
 * argument, but need an additional type.
 */
template <class T> template <class U>
std::vector<double> distance_transform<T>::eval_grid_function( std::vector<U> &data ){

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
template <class T>
void distance_transform<T>::compute_distance_map(){

  if( dim == 1 ){
    //one dimensional case
    
    //get grid
    std::vector<double> grid = eval_grid_function( img );

    //do the distance transform
    map = do_distance_transform( grid, pix_size[0] );
  } //end dim=1 case

  else if ( dim == 2 ){
    //2d case

    // the object holding the 2d distance map col wise
    std::vector<double> cols_map ( map.size(), 0 );
    
    //do the distance transform for each column first
    for( unsigned int ii=0; ii < size[0]; ii++ ){
      
      //get the ii-th column
      std::vector<double> col (size[1], 0);
      for( unsigned int jj=0; jj<size[1]; jj++){
	col[jj] = img[ii + jj * size[0]];
      }
      
      //get the grid
      std::vector<double> col_grid = eval_grid_function( col );
      
      //get distance map of column
      std::vector<double> col_map = do_distance_transform( col_grid, pix_size[1] );
      
      //cpy preliminary map in to map vector
      memcpy( cols_map.data()+(ii*size[1]), col_map.data(), sizeof(double) * col_map.size() );
      
    } //end column wise distance transform

    
    /*
     * careful, map is now not row based, but column based. We need to
     * reverse that in the following step
     */

    //now iterate over all rows and compute the distance map row wise
    for( unsigned int ii=0; ii<size[1]; ii++){

      //get the iith row
      std::vector<double> row (size[0], 0);
      for( unsigned int jj=0; jj<size[0]; jj++){
	row[jj] = cols_map[ ii + jj*size[1] ];
      }

      //now get the transform
      std::vector<double> row_map = do_distance_transform( row, pix_size[0] );

      //cpy the final map data to map object
      memcpy( map.data()+(ii*size[0]), row_map.data(), sizeof(double) * row_map.size() );
    } // end row wise (and final) distance transform

    
  } //end dim=2 case
  
  else if ( dim == 3 ){

    //same story as with dim = 2

    // the distance map after the stack wise distance map computation
    std::vector<double> map_stacks ( map.size(), 0 );
    // the distance map after the column wise distance map computation
    std::vector<double> map_cols ( map.size(), 0 );    
    // the final transofrm, but in the wron order
    std::vector<double> map_unsorted ( map.size(), 0 );
    
    //start with distance transform of each stack
    for( unsigned int ii=0; ii<size[0]*size[1]; ii++){

      //get the stack
      std::vector<double> stack ( size[2], 0 );
      for(unsigned int jj=0; jj<size[2]; jj++){
	stack[jj] = img[jj + ii*size[2] ];
      }
      
      // apply edm cost function
      std::vector<double> stack_grid = eval_grid_function( stack );

      // compute distance map of stack
      std::vector<double> stack_map = do_distance_transform( stack_grid, pix_size[2] );

      //copy stack distance into map
      memcpy( map_stacks.data() + ii*size[2], stack_map.data(), stack_map.size() * sizeof(double) );
    } //end stack transform

    
    /*
     * map_stacks now should have the same index structure as our
     * original grid
     */
    
    // do the column wise transform
    for( unsigned int ii=0; ii<size[0]*size[2]; ii++){

      // get the columns
      std::vector<double> col (size[1], 0 );
      for( unsigned int jj=0; jj<size[1]; jj++){
	int ind = ii + jj*(size[2]*size[0]);
	col[jj] = map_stacks[ ind ];
      }

      //distance transform of columns
      std::vector<double> col_map = do_distance_transform( col, pix_size[1] );
      
      //copy data
      memcpy( map_cols.data()+ii*col_map.size(), col_map.data(), sizeof(double) * col_map.size() );
    }

    /*
     * there is a new order. To access the rows in the original array,
     * we now need to access the columns in the new order. Not that
     * the dimensions changed: n-->l ; k --> n ; k --> l
     */

    //do the distance transofrm on each row(old order)/col(current
    //order)
    for(unsigned int ii=0; ii<size[1]*size[2]; ii++){

      //get the row/col
      std::vector<double> row (size[0], 0 );
      for( unsigned jj=0; jj<size[0]; jj++ ){
	int ind = ii + jj*size[1]*size[2];
	row[jj] = map_cols[ind];
      }

      //do the distance transform
      std::vector<double> row_map = do_distance_transform( row, pix_size[0] );

      //store map in temp map
      memcpy( map_unsorted.data()+ii*row_map.size(), row_map.data(), sizeof(double) * row_map.size());
    } // end row distance map

    //we need some sorting to arrange the grid in the original order
    for(unsigned int ii=0; ii<size[0]*size[1]; ii++){

      //get the original stacks
      std::vector<double> stack (size[2], 0);
      for(unsigned int jj=0; jj<size[2]; jj++){
	stack[jj] = map_unsorted[ii + jj*size[0]*size[1]];
      }

      //store final distance map in the map array in correct order
      memcpy( map.data()+ii*stack.size(), stack.data(), sizeof(double)*stack.size() );
    } // end final sorting
    
  }

}


/*
 * getters
 */

/**
 * Returns a copy of the distance map. The distances are given in
 * squared distances of the same length units as the pixel
 * sizes. Dimensionless if pixelsizes=1
 */
template <class T>
std::vector<double> distance_transform<T>::get_distance_map() const {
  return map;
}
