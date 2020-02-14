#include "homotopic_thinning.hpp"


/*
 * Constructor
 * \param[in] nx size of the image in x direction
 * \param[in] ny size of the image in y direction
 * \param[in] nz size of the image in z direction
 */
homotopic_thinning::homotopic_thinning ( int nx, int ny, int nz, std::vector<int> image_, std::vector<double> dmap_ )  : n_points_x (nx), n_points_y (ny), n_points_z (nz), image( image_ ), dmap( dmap_ ) {
  
}

/*
 * Just a small helper function to compute the intersection of two
 * sets. Unfortunatly O(N^2)
 */
std::unordered_set<int> intersection( std::unordered_set<int> a, std::unordered_set<int> b ){

  std::unordered_set<int> intersection_set;

  for( auto it_a : a ){
    for( auto it_b : b ){
      if( it_a == it_b){
	intersection_set.insert(it_a);
	break;
      }
    }
  }
  
  return intersection_set;
}




/*
 * Returns a set of all exterior points of the channel with
 * channel_id. Exterior points are defined by being m_adjacent to the
 * channel as well as the background ( being any voxel with a channel
 * id different from channel_id
 *
 * param[in] channel_id Id of the channel to get the exterior points for
 * param[in] n background adjacency
 */
std::unordered_set<int> homotopic_thinning::get_exterior_points( int channel_id, int n ){

  std::unordered_set<int> exterior_points;

  for( unsigned int ii=0; ii< image.size(); ii++){

    if( std::abs(image[ii]) == std::abs( channel_id ) ){
    
      iterable_voxel point ( ii, n_points_x, n_points_y, n_points_z );
      
      std::unordered_set<int> nbs;
      if( n == 6 ){
	nbs = point.get_6_neighbors();
      } else if( n == 18 ){
	nbs = point.get_18_neighbors();
      } else if ( n == 26 ){
	nbs = point.get_26_neighbors();
      } else {
	throw std::string("no adjancy model found for n=" + std::to_string(n) );
      }

      bool found_background = false;
      for( auto it : nbs ){
	if( std::abs( image[it] ) != std::abs( channel_id ) ){
	  found_background = true;
	  break;
	}
      }
            
      if( found_background ){
	exterior_points.insert( point() );
      }

    }
  }

  return exterior_points;
  
}



/*
 * Checks, if a point is deletable. So far, only the Bertrand and
 * Malandain definition is implemented, maybe the Pudney definition
 * will be added later on
 *
 * \param[in] id The id of the point to check
 * \param[in] m The connectivity of the image
 * \param[in] n The connectivity of the background
 */
bool homotopic_thinning::point_deletable( int id, int m, int n ){

  auto components_image = get_connected_components( id, m, true );
  auto components_background = get_connected_components( id, n, false );  

  //now we need to check connectivity with the components. Lets start
  //of with the image components. We can only be connected to one of
  //them

  //it's probably easiest to compute the neighborhood of id and see,
  //if that has overlap with the component

  int nr_img_connected_components = 0;
  for( unsigned int ii=0; ii<components_image.size(); ii++){

    iterable_voxel point ( id, n_points_x, n_points_y, n_points_z);
    auto nbs = point.get_26_neighbors();

    auto common_points = intersection( nbs, components_image[ii] );

    //they are connected, if size >0
    if( common_points.size() > 0 ){
      nr_img_connected_components++;
    }
    
  }
  
  int nr_background_connected_components = 0;
  for( unsigned int ii=0; ii<components_background.size(); ii++){

    iterable_voxel point ( id, n_points_x, n_points_y, n_points_z);
    auto nbs = point.get_6_neighbors();

    auto common_points = intersection( nbs, components_background[ii] );

    //they are connected, if size >0
    if( common_points.size() > 0 ){
      nr_background_connected_components++;
    }
    
  }
 
  if( nr_img_connected_components == 1 &&
      nr_background_connected_components == 1){
    return true;
  } else {
    return false;
  }

}



/*
 * This function computes the homotopic, medial skeleton of the 3D
 * object, which is the channel with the id ch_id. Implementation of
 * the pseudo-code in Pudney (1998) (I think it is inspired by Vincent
 * (1991) )
 *
 * \param[in] ch_id The id of the channel to compute the skeleton of
 *
 * \return A set of points making up the skeleton
 */
std::unordered_set<int> homotopic_thinning::find_channel_skeleton( int ch_id ){



  //adjacency: (m, n) = (foreground, background)connectivity = (26, 6)
  int m = 26, n = 6;

  //get exterior points. They are 
  auto ext_points = get_exterior_points( ch_id, n );
  
  
  std::unordered_set<int> graph;


  //label if point is already queued for processing
  std::vector<int> queued (image.size(), 0);

  // add all skeleton points to the queue
  std::multimap<double, int> to_process;
  for( auto it : ext_points ){
    to_process.emplace( dmap[it], it );
    queued[it] = 1;
  }


  while( !to_process.empty() ){

    // pop the closest point to the membrane 
    iterable_voxel point ( to_process.begin()->second, n_points_x, n_points_y, n_points_z );
    to_process.erase( to_process.begin() );

    //point is no longer queued, since it's processed now
    queued[point()] = 0;

    //check if the point is safe to delete
    if( point_deletable( point(), m, n ) ){

      //yes, dump it!!!
      image[point()] = 0;

      //now add all neighbors of the point
      
      //get the neighbors
      std::unordered_set<int> nbs = point.get_26_neighbors();
      
      for( auto it : nbs ){
	if( std::abs( image[it] ) == std::abs(ch_id)   &&
	    queued[it] == 0 ){
	  
	  to_process.emplace( dmap[it], it );
	  queued[it] = 1;
	}

      }
      
    } else {
      graph.insert( point() );
    }
    
  }


  return graph;  
   
  
}


/*
 * This function counts the number of connected components around
 * center with the same id as the center point. The image is defined
 * as all points with the same \ref channel id as the center point,
 * background are all other channel ids
 *
 * \param[in] center The id of the particle to look for components for
 * \param[in] m The adjacency (6, 18, 26)
 * \param[in] img If true, count image components, false count background components
 *
 * \return A vector of unordered sets containing the points of the connected components
 */
std::vector< std::unordered_set<int> >  homotopic_thinning::get_connected_components( int center, int m, bool img ){


  // at first, lets find the suitable neighborhood of center
  // labeled N*(p)


  //depending on if we are working with (m,n)=(26,6) or (m,n)=(6,26)
  // the neighborhoods are defined differently
  iterable_voxel midpoint( center, n_points_x, n_points_y, n_points_z );
  std::unordered_set<int> nbs;
  if( m == 6 ){
    nbs = midpoint.get_18_neighbors();
  } else if ( m == 26 ) {
    nbs = midpoint.get_26_neighbors();
  } else {
    throw ("error");
  }


  // Now we need to compute N*(p) intersected with F (the image)
  // Or in case we are computing components of the background
  // N*(p) - F
  
  // delete all image/background points
  auto del_it = nbs.begin();
  while( del_it != nbs.end() ){

    if( img ){
      // delete all background points
      if( image[*del_it] != image[ midpoint() ] ){
	del_it = nbs.erase( del_it );
      } else {
	del_it++;
      }
    } else {
      // delete all image points
      if( image[*del_it] == image[ midpoint() ] ){
	del_it = nbs.erase( del_it );
      } else {
	del_it++;
      }
    }
    
  }

  //so now for img = ture, nbs is N*(p) intersect F. in case img =
  //false nbs = N*(p) -F
  
  // Now lets find the components
  
  std::unordered_set<int> to_visit;
  std::unordered_set<int> visited;
  std::vector< std::unordered_set<int> > components;

  int component_count = 0;
  
  while( visited.size() < nbs.size() ){

    to_visit.clear();
    components.push_back( std::unordered_set<int>() );
    
    //start somewhere
    auto find_next_point = nbs.begin();
    while( visited.find( *find_next_point ) != visited.end() ){
      find_next_point++;
    }

    
    to_visit.insert( *find_next_point );
    components[component_count].insert( *find_next_point );
    
    //now get all points connected to that point
    while( !to_visit.empty() ){
      
      // start somewhere
      iterable_voxel it ( *to_visit.begin(), n_points_x, n_points_y, n_points_z );
      to_visit.erase( it() );
      visited.insert( it() );
      
      // get all image points around it within the shell around the
      // midpoint
      
      // depending on if we are computing components for image or
      // background we ahve different neighbourhood definitions
      std::unordered_set<int> local_nbs;
      if( img ){
	local_nbs = it.get_26_neighbors();
      } else {
	local_nbs = it.get_6_neighbors();
      } 
      
      // now we have *all* points around nbs, we only want the
      // neighbouring points in the shell around midpoint. Compute the
      // intersection therefore
      local_nbs = intersection( nbs, local_nbs );
      
      // add all points to the component set and to the to_visit set
      for( auto el : local_nbs ){
	components[component_count].insert( el );
	if( visited.find( el ) == visited.end() ){
	  to_visit.insert( el );
	}
      }
      
    }

    //finish component and initialize next one
    component_count++;

  }

  return components;  
}
