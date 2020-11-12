#include "percolation_analysis.hpp"


/**
 * Constructor
 *
 * \param[in] data The input grid to analyse
 * \param[in] distance_map The distance map of the input grid
 * \param[in] sx The dimension (nr of pixels) of input grid in x direction
 * \param[in] sy The dimension (nr of pixels) of input grid in y direction
 * \param[in] sz The dimension (nr of pixels) of input grid in z direction
 * \param[in] is_periodic Is the input grid periodic
 */
template <class T, class M>
percolation_analysis<T, M>::percolation_analysis(std::vector<T> data, std::vector<M> distance_map, unsigned int sx_, unsigned int sy_, unsigned int sz_, bool is_periodic) : sx(sx_), sy(sy_), sz(sz_), periodic( is_periodic), structure(data), dmap (distance_map ), it( 0, sx_, sy_, sz_ ) {

  cluster_labels.resize( structure.size(), -1 );
  cluster_sizes.clear();
}


/*
 * \param[in] id The index of the particle to get the neighbors for
 * \param[in] ch_id The id of the channel currently considered
 * \param[in] nbs The vector to write the id's of neighbors into
 */
template <class T, class M>
void percolation_analysis<T, M>::get_nbs( int id, int ch_id, std::vector<int> &nbs ){

  it.set(id);
  nbs.clear();

  if( periodic ){
    // get 6 neighbors, but only backwards directions where we already
    // labeled paritcles!
    if( cluster_labels[ it.d()() ] >= 0 ){
      nbs.push_back( it.d()() );
    }
    if( cluster_labels[ it.l()() ] >= 0 ){
      nbs.push_back( it.l()() );
    }
    if( cluster_labels[ it.b()() ] >= 0 ){
      nbs.push_back( it.b()() );
    }
    if( cluster_labels[ it.u()() ] >= 0 ){
      nbs.push_back( it.u()() );
    }
    if( cluster_labels[ it.r()() ] >= 0 ){
      nbs.push_back( it.r()() );
    }
    if( cluster_labels[ it.f()() ] >= 0 ){
      nbs.push_back( it.f()() );
    }
  } else {

    int x, y, z, c;
    z = id % sz;
    c = int( id / sz );
    y = int( c / sx );
    x = c % sx ;
    
    if( z != 0 ){ // if at lower border, don't add neighbor to avoid
		  // pbcs
      if( cluster_labels[ it.d()() ] >= 0 ){
	nbs.push_back( it.d()() );
      }
    }

    if( x != 0 ){
      if( cluster_labels[ it.l()() ] >= 0 ){
	nbs.push_back( it.l()() );
      }
    }

    if( y != 0 ){
      if( cluster_labels[ it.b()() ] >= 0 ){
	nbs.push_back( it.b()() );
      }
    }

  }
}


/**
 * Implementation of Hoshen-Koppelmann algorithm to find clusters in
 * the structure
 *
 * \param[in] ch_id The value of the 'occupied' sites, thus the sites to assume to belong to clusters
 *
 */
template <class T, class M>
void percolation_analysis<T, M>::find_clusters( int ch_id ){

  // reset the label map and sizes array
  memset( cluster_labels.data(), -1, sizeof(int) * cluster_labels.size() );
  cluster_sizes.clear();

  
  unsigned int id = 0;
  std::vector<int> nbs;

  // max label
  int k = -1; //start at -1, since increment is 1, so first cluster will have label 0


  // iterate over all voxels
  while( id < structure.size()  ){
    // only consider occupied sites of correct ch_id
    if( std::fabs( structure[id] ) == std::fabs( ch_id ) ){
	    
      // get neighbors
      it.set(id);
      get_nbs( id, ch_id, nbs );

      // no neighbor, create new cluster
      if( nbs.size() == 0 ){

	// increase label
	k++;
	// assigne label
	cluster_labels[id] = k;
	// initialise cluster size
	cluster_sizes.push_back( 1 );
	
      } else if ( nbs.size() == 1 ){

	// only one neighbor. Find the true label
	int true_label = cluster_labels[ nbs[0] ];

	while( cluster_sizes[ true_label ] <= 0 ){
	  true_label = -cluster_sizes[ true_label ];
	}

	// now we found the true label and can assign it
	cluster_labels[id] = true_label;
	cluster_sizes[true_label]++;
	
      } else {
	// found multiple neighbors need to join clusters
	
	std::vector<int> true_labels ( nbs.size(), 0 );
	// first get true cluster labels for each cluster
	for( unsigned int ii=0; ii<nbs.size(); ii++){

	  int true_label = cluster_labels[ nbs[ii] ];
	  while( cluster_sizes[ true_label ] <= 0 ){
	    true_label = -cluster_sizes[ true_label ];
	  }

	  true_labels[ii] = true_label;
	}

	// find minimum label, this is the new true label of merged cluster
	int min_label = *std::min_element( true_labels.begin(), true_labels.end() );

	// now we can assign the new label
	cluster_labels[id] = min_label;

	/*
	 * in case we are merging several clusters with identical true
	 * labels, make sure, we only adding the cluster mass once!
	 */
	std::vector<int> already_added (0, 0);
	
	// get the new cluster size
	for(unsigned int ii=0; ii<true_labels.size(); ii++){	  
	  if( true_labels[ii] != min_label ){

	    if( std::find( already_added.begin(), already_added.end(), true_labels[ii] ) == already_added.end() ){
	      cluster_sizes[ min_label ] += cluster_sizes[ true_labels[ii] ];
	      already_added.push_back( true_labels[ii] );
	    }

	    cluster_sizes[ true_labels[ii] ] = -min_label;
	  }
	}

	cluster_sizes[min_label]++;	
      }         
    }

    id++;
  }
}


/** 
 * Iterates over cluster map and assigns each site its true
 * label. This is more a cosmetic thing.
 */
template <class T, class M>
void percolation_analysis<T, M>::assign_true_labels(){

  for( unsigned int ii=0; ii < cluster_labels.size(); ii++){

    if( cluster_labels[ii] < 0 ){
      cluster_labels[ii] = 0;
    } else {
    
      // find true label
      int true_label = cluster_labels[ ii ];

      while( cluster_sizes[ true_label ] <= 0 ){
	true_label = -cluster_sizes[ true_label ];
    }
      
      cluster_labels[ii] = true_label + 1;
    }
  }
}
  

/**
 * Returns an array containing the masses of the clusters. Essentially
 * picks out the non-negative values from the \ref cluster_sizes array
 */
template <class T, class M>
std::vector<unsigned int> percolation_analysis<T, M>::get_cluster_sizes() const {
  std::vector<unsigned int> ret;
  for( auto it : cluster_sizes ){
    if( it >= 1 ){
      ret.push_back( it );
    }
  }
  return ret;
}

/**
 * Returns the number of clusters found, by counting non-negative
 * entries in \ref cluster_sizes
 */
template <class T, class M>
unsigned int percolation_analysis<T, M>::get_nr_clusters() const {

  unsigned int count = 0;
  for( auto it : cluster_sizes ){
    if( it >= 1 ){
      count++;
    }
  }
  return count;
}


/**
 * Prints label map to a file. Wrapper around \ref print function
 * \param[in] out The filename to write to
 */
template <class T, class M>
void percolation_analysis<T, M>::print( std::string out ){
  std::ofstream fout ( out );
  print( fout );
  fout.close();  
}


/*
 *
 */
template <class T, class M>
void percolation_analysis<T, M>::print( std::ostream &out) {

  int max = *std::max_element( cluster_labels.begin(), cluster_labels.end() );
  
  int nr_digits = int( log( float(max) ) / log (10.0) );
  
  for( unsigned int ii=0; ii < sz; ii++){ //layers
    for( unsigned int jj=0; jj < sy; jj++){ //columns
      for( unsigned int kk=0; kk< sx; kk++){ //rows

	int ind = ii + kk * sz + jj * sz*sx;	
	
	out << std::setw(nr_digits + 1) << cluster_labels[ind] << " ";	  
      }
      out << std::endl;
    }
    out << std::endl << std::endl;
  }
}


template class percolation_analysis<short, int>;
template class percolation_analysis<int, int>;
template class percolation_analysis<unsigned int, int>;
template class percolation_analysis<float, int>;
template class percolation_analysis<double, int>;

template class percolation_analysis<short, unsigned int>;
template class percolation_analysis<int, unsigned int>;
template class percolation_analysis<unsigned int, unsigned int>;
template class percolation_analysis<float, unsigned int>;
template class percolation_analysis<double, unsigned int>;

template class percolation_analysis<short, float>;
template class percolation_analysis<int, float>;
template class percolation_analysis<unsigned int, float>;
template class percolation_analysis<float, float>;
template class percolation_analysis<double, float>;

template class percolation_analysis<short, double>;
template class percolation_analysis<int, double>;
template class percolation_analysis<unsigned int, double>;
template class percolation_analysis<float, double>;
template class percolation_analysis<double, double>;
