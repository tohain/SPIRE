#include "percolation_analysis.hpp"

template <class T>
percolation_analysis<T>::percolation_analysis(std::vector<T> data, std::vector<T> distance_map, unsigned int sx_, unsigned int sy_, unsigned int sz_) : sx(sx_), sy(sy_), sz(sz_), structure(data), dmap (distance_map ), it( 0, sx_, sy_, sz_ ) {

  cluster_labels.resize( structure.size(), -1 );
  cluster_sizes.clear();
}


/*
 * \param[in] id The index of the particle to get the neighbors for
 * \param[in] ch_id The id of the channel currently considered
 * \param[in] nbs The vector to write the id's of neighbors into
 */
template <class T>
void percolation_analysis<T>::get_nbs( int id, int ch_id, std::vector<int> &nbs ){

  it.set(id);
  nbs.clear();
  
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
}


/**
 * Implementation of Hoshen-Koppelmann algorithm to find clusters in
 * the structure
 */
template <class T>
void percolation_analysis<T>::find_clusters(){

  int ch_id = 1;

  unsigned int id = 0;
  std::vector<int> nbs;

  // max label
  int k = -1; //start at -1, since increment is 1, so first cluster will have label 0
  
  while( id < structure.size()  ){


    if( std::fabs( structure[id] ) == std::fabs( ch_id ) ){
	    
      // get neighbors
      it.set(id);
      get_nbs( id, ch_id, nbs );

      if( nbs.size() == 0 ){

	// increase label
	k++;
	cluster_labels[id] = k;
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

	// find minimum label
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

	    // in case, we are 
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



template <class T>
void percolation_analysis<T>::assign_true_labels(){

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
  


template <class T>
std::vector<unsigned int> percolation_analysis<T>::get_cluster_sizes() const {
  std::vector<unsigned int> ret;
  for( auto it : cluster_sizes ){
    if( it >= 1 ){
      ret.push_back( it );
    }
  }
  return ret;
}

template <class T>
unsigned int percolation_analysis<T>::get_nr_clusters() const {

  unsigned int count = 0;
  for( auto it : cluster_sizes ){
    if( it >= 1 ){
      count++;
    }
  }
  return count;
}


template <class T>
void percolation_analysis<T>::print( std::string out ){
  std::ofstream fout ( out );
  print( fout );
  fout.close();  
}


template <class T>
void percolation_analysis<T>::print( std::ostream &out) {
  
  for( unsigned int ii=0; ii < sz; ii++){ //layers
    //out << "Layer " << ii << std::endl;
    for( unsigned int jj=0; jj < sy; jj++){ //columns
      for( unsigned int kk=0; kk< sx; kk++){ //rows

	int ind = ii + kk * sz + jj * sz*sx;	
	
	out << " " << std::setw(4) << cluster_labels[ind] << " ";	  
      }
      out << std::endl;
    }
    out << std::endl << std::endl;
  }
}

template class percolation_analysis<int>;
template class percolation_analysis<unsigned int>;
template class percolation_analysis<float>;
template class percolation_analysis<double>;
