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
 * Initializes an empty object. Needs a call of \ref set_parameters()
 * before being able to do stuff!
 */
template <class T, class M>
distance_transform<T, M>::distance_transform(){

  img = std::vector<T> (0, 0);
  dim = 1;
  size[0]=1; size[1]=1; size[2];
  pix_size[0]=1; pix_size[1]=1; pix_size[2]=1;

  map = std::vector<M> (0, 0);
}



/**
 * Allocates necessary memory for the transformed map and initializes
 * class parameters
 *
 * \param[in] data The data to transform
 * \param[in] n The length of the data - not necessary, but for consistency
 * \param[in] pixsize The size of a pixel in length units
 */
template <class T, class M>
distance_transform<T, M>::distance_transform( std::vector<T> data, int n, double pixsize ){
  //initialize parameters
  img = data;
  dim = 1;  
  size[0]=n; size[1]=1; size[2]=1;
  pix_size[0]=pixsize; pix_size[1]=1; pix_size[2]=1;
  
  //allocate memory
  map = std::vector<M> (n, 0);
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
template <class T, class M>
distance_transform<T, M>::distance_transform( std::vector<T> data, int n, int k, double pixsize_x, double pixsize_y ){
  //initialize parameters
  img = data;  
  dim = 2;
  size[0]=n; size[1]=k; size[2]=1;
  pix_size[0]=pixsize_x; pix_size[1]=pixsize_y; pix_size[2]=1;

  //allocate memory
  map = std::vector<M> (n*k, 0);
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
template <class T, class M>
distance_transform<T, M>::distance_transform( std::vector<T> data, int n, int k, int l, double pixsize_x, double pixsize_y, double pixsize_z ){
  //initialize parameters
  img = data;
  dim = 3;
  size[0]=n; size[1]=k; size[2]=l;
  pix_size[0]=pixsize_x; pix_size[1]=pixsize_y; pix_size[2]=pixsize_z;

  //allocate memory
  map = std::vector<M> (n*k*l, 0);
}



/**
 * Updates the object with a new grid, dimension and pixsize. The
 * function ref compute_distance_map() needs to be called again of
 * course to update the distance map returned by ref
 * get_distance_map();
 * \param[in] data The array holding the image data 
 * \param[in] dim The dimension of the image: [0] the width (=nr of
 * cols); [1] the height (=nr of rows); [3] the depth (=height of
 * "stacks"). The size of the vector determines the dimension of the
 * distance map
 * \param[in] pixsize The size of a pixel in its respective direction
 * in length units
 */
template<class T, class M>
void distance_transform<T, M>::set_parameters( std::vector<T> data,
					 std::vector<unsigned int> _dim,
					 std::vector<double> pixsize ){
  
  img = data;
  dim = _dim.size();

  if( _dim.size() == 1 ){

    size[0]=_dim[0]; size[1]=1; size[2]=1;
    pix_size[0]=pixsize[0]; pix_size[1]=1; pix_size[2]=1;

    map = std::vector<M> ( _dim[0], 0);
    
  } else if( _dim.size() == 2 ){

    size[0]=_dim[0]; size[1]=_dim[1]; size[2]=1;
    pix_size[0]=pixsize[0]; pix_size[1]=pixsize[1]; pix_size[2]=1;

    map = std::vector<M> ( _dim[0]*_dim[1], 0 );
    
  } else if( _dim.size() == 3 ){

    size[0]=_dim[0]; size[1]=_dim[1]; size[2]=_dim[2];
    pix_size[0]=pixsize[0]; pix_size[1]=pixsize[1]; pix_size[2]=pixsize[2];

    map = std::vector<M> ( _dim[0]*_dim[1]*_dim[2], 0 );
    
  }
  
}


/**
 * Destructor, just frees memory
 */
template <class T, class M>
distance_transform<T, M>::~distance_transform(){
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
template <class T, class M>
std::vector<M> distance_transform<T, M>::do_distance_transform( std::vector<M> &grid, double pixsize ){

  // The grid will be in length units. In it's initial call in
  // 1d/2d/3d it is just max or 0 (in both cases units don't
  // matter). After a second call in 2d/3d the units will be in length
  // units, so any grid[*] is in length units


  std::vector<M> transform ( grid.size(), 0 );
  
  //Index of rightmost parabola in lower envelope
  unsigned int k = 0;
  //indeces/locations of parabolas in lower envelope
  std::vector<int> v = {0}; 
  //Locations of boundaries between parabolas in length units
  std::vector<M> z = {-pixsize*std::numeric_limits<M>::max(), pixsize*std::numeric_limits<M>::max()};
  
  //find the envelope
  for( unsigned int ii=1; ii<grid.size(); ii++){
    //intersection of two parabolas

    //the position of the index in length units
    M ii_l = ii*pixsize;
    
    //this needs to be in length units already! so convert pixel
    //values to length units where neccessary
    M s = ( ( grid[ii] + (ii_l*ii_l) ) - (grid[ v[k] ] + (v[k]*v[k]*pixsize*pixsize)) ) / ( 2*( ii_l - (pixsize*v[k]) ) );

    //check where the intersection is. If it is left of the previous
    //parabola, the current one is the new lowest envelope
    while( s <= z[k] ){
      k--;
      s = ( ( grid[ii] + (ii_l*ii_l) ) - (grid[ v[k] ] + (v[k]*v[k]*pixsize*pixsize)) ) / (2*(ii_l-(v[k]*pixsize)));
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
      z[k+1] = std::numeric_limits<M>::max();
    } else {
      z.push_back( std::numeric_limits<M>::max() );
    }    
  }


  //fill in the transform by evaluating the parabola envelope at the
  //grid
  k=0;
  for(unsigned int ii=0; ii<grid.size(); ii++){
    while( z[k+1] < ii*pixsize ){
      k++;
    }
    
    //get the "real" distances in length units by multiplying squared
    //index distances by pixelsize
    transform[ii] = pixsize*pixsize*( (ii-v[k])*(ii-v[k]) ) + grid[ v[k] ];
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
 *
 *\param[in] data the grid on which to evaluate the function
 *\param[in] max The max value the cost function will take 
 */
template <class T, class M> template <class U>
std::vector<M> distance_transform<T, M>::eval_grid_function( std::vector<U> &data, double max ){

  //initialize memory and set it to zero
  std::vector<M> grid (data.size(), 0);
  
  for(unsigned int ii=0; ii<data.size(); ii++){
    if( data[ii] == 0 ){
      grid[ii] = max;
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
template <class T, class M>
void distance_transform<T, M>::compute_distance_map(){

  if( dim == 1 ){
    //one dimensional case
    
    //get grid

    //set the max value to a couple of lengths of the entire grid
    M max = size[0] * pix_size[0] * 10;
    std::vector<M> grid = eval_grid_function<T>( img, max );

    //do the distance transform
    map = do_distance_transform( grid, pix_size[0] );
  } //end dim=1 case

  else if ( dim == 2 ){
    //2d case

    // the object holding the 2d distance map col wise
    std::vector<M> cols_map ( map.size(), 0 );
    
    //do the distance transform for each column first
    for( unsigned int ii=0; ii < size[0]; ii++ ){
      
      //get the ii-th column
      std::vector<M> col (size[1], 0);
      for( unsigned int jj=0; jj<size[1]; jj++){
	col[jj] = img[ii + jj * size[0]];
      }
      
      //get the grid
      M max = std::max(size[0],size[1]) * std::max(pix_size[0],pix_size[1]) * 10;
      std::vector<M> col_grid = eval_grid_function<M>( col, max );
      
      //get distance map of column
      std::vector<M> col_map = do_distance_transform( col_grid, pix_size[1] );
      
      //cpy preliminary map in to map vector
      memcpy( cols_map.data()+(ii*size[1]), col_map.data(), sizeof(M) * col_map.size() );
      
    } //end column wise distance transform

    
    /*
     * careful, map is now not row based, but column based. We need to
     * reverse that in the following step
     */

    //now iterate over all rows and compute the distance map row wise
    for( unsigned int ii=0; ii<size[1]; ii++){

      //get the iith row
      std::vector<M> row (size[0], 0);
      for( unsigned int jj=0; jj<size[0]; jj++){
	row[jj] = cols_map[ ii + jj*size[1] ];
      }

      //now get the transform
      std::vector<M> row_map = do_distance_transform( row, pix_size[0] );

      //cpy the final map data to map object
      memcpy( map.data()+(ii*size[0]), row_map.data(), sizeof(M) * row_map.size() );
    } // end row wise (and final) distance transform

    
  } //end dim=2 case
  
  else if ( dim == 3 ){

    //same story as with dim = 2

    // the distance map after the stack wise distance map computation
    std::vector<M> map_stacks ( map.size(), 0 );
    // the distance map after the column wise distance map computation
    std::vector<M> map_cols ( map.size(), 0 );    
    // the final transofrm, but in the wron order
    std::vector<M> map_unsorted ( map.size(), 0 );
    
    //start with distance transform of each stack
    for( unsigned int ii=0; ii<size[0]*size[1]; ii++){

      //get the stack
      std::vector<M> stack ( size[2], 0 );
      for(unsigned int jj=0; jj<size[2]; jj++){
	stack[jj] = img[jj + ii*size[2] ];
      }
      
      // apply edm cost function
      M max = std::max( std::max( size[0], size[1] ), size[2]) *
	std::max( std::max( pix_size[0], pix_size[1] ), pix_size[2] ) * 10;
      std::vector<M> stack_grid = eval_grid_function<M>( stack, max );

      // compute distance map of stack
      std::vector<M> stack_map = do_distance_transform( stack_grid, pix_size[2] );

      //copy stack distance into map
      memcpy( map_stacks.data() + ii*size[2], stack_map.data(), stack_map.size() * sizeof(M) );
    } //end stack transform

    
    /*
     * map_stacks now should have the same index structure as our
     * original grid
     */
    
    // do the column wise transform
    for( unsigned int ii=0; ii<size[0]*size[2]; ii++){

      // get the columns
      std::vector<M> col (size[1], 0 );
      for( unsigned int jj=0; jj<size[1]; jj++){
	int ind = ii + jj*(size[2]*size[0]);
	col[jj] = map_stacks[ ind ];
      }

      //distance transform of columns
      std::vector<M> col_map = do_distance_transform( col, pix_size[1] );
      
      //copy data
      memcpy( map_cols.data()+ii*col_map.size(), col_map.data(), sizeof(M) * col_map.size() );
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
      std::vector<M> row (size[0], 0 );
      for( unsigned jj=0; jj<size[0]; jj++ ){
	int ind = ii + jj*size[1]*size[2];
	row[jj] = map_cols[ind];
      }

      //do the distance transform
      std::vector<M> row_map = do_distance_transform( row, pix_size[0] );

      //store map in temp map
      memcpy( map_unsorted.data()+ii*row_map.size(), row_map.data(), sizeof(M) * row_map.size());
    } // end row distance map
    
    //we need some sorting to arrange the grid in the original order
    for(unsigned int ii=0; ii<size[0]*size[1]; ii++){

      //get the original stacks
      std::vector<M> stack (size[2], 0);
      for(unsigned int jj=0; jj<size[2]; jj++){
	stack[jj] = map_unsorted[ii + jj*size[0]*size[1]];
      }

      //store final distance map in the map array in correct order
      memcpy( map.data()+ii*stack.size(), stack.data(), sizeof(M)*stack.size() );
    } // end final sorting
    
  }

}


/**
 * Outputs the distance map to the console
 */
template <class T, class M>
void distance_transform<T, M>::print_map() const {
  
  for( unsigned int ii=0; ii < size[2]; ii++){ //layers
    std::cout << "Layer " << ii << std::endl;
    for( unsigned int jj=0; jj < size[1]; jj++){ //columns
      for( unsigned int kk=0; kk< size[0]; kk++){ //rows

	int ind = ii + kk * size[2] + jj * size[2]*size[0];	
	
	std::cout << " " << std::setprecision(3) << std::setw(6) << sqrt(map[ind]) << " ";	  
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
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
template <class T, class M>
std::vector<M> distance_transform<T, M>::get_distance_map() const {
  return map;
}


// explicit template instantiations
template class distance_transform<short, float>;
template class distance_transform<int, float>;
template class distance_transform<float, float>;
template class distance_transform<double, float>;


