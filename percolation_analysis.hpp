#ifndef PERCOLATION_ANALYSIS_H
#define PERCOLATION_ANALYSIS_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstring>

#include "iterable_voxel.hpp"

/**
 * This class implements some percolation analysis to find the largest
 * sphere that can travel through the structures
 */
template <class T, class M>
class percolation_analysis {

public:

  /// Constructor. Copies and initalises data and arrays
  percolation_analysis( std::vector<T> data, std::vector<M> distance_map, unsigned int sx_, unsigned int sy_, unsigned int sz_, bool periodic);

  /// Using Hoshen-Kopelmann to find clusters
  void find_clusters( int ch_id );

  /// makes the cluster map a bit nicer
  void assign_true_labels();
  
  /// prints the cluster map to a file
  void print( std::ostream &out );
  void print( std::string out );

  /// retrieves all the cluster sizes from \ref cluster_sizes
  std::vector<unsigned int> get_cluster_sizes() const;

  /// returns the number of clusters found
  unsigned int get_nr_clusters() const;

  
  
private:

  /// returns the neighbors of voxel id which belong to the channel id ch_id
  void get_nbs( int id, int ch_id, std::vector<int> &nbs );

  /// the grid to analyse
  std::vector<T> structure;
  /// the distance map of the grid provided. Needs to dilute surfaces
  std::vector<M> dmap;

  /// the map assigning each voxel its cluster label
  std::vector<int> cluster_labels;
  /// an array keeping track of true labels as well as cluster sizes
  std::vector<int> cluster_sizes;
  


  /// size of the simbox
  unsigned int sx, sy, sz;
  /// is it periodic
  bool periodic;

  
  iterable_voxel it;  
};





#endif
