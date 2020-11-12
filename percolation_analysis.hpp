#ifndef PERCOLATION_ANALYSIS_H
#define PERCOLATION_ANALYSIS_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "iterable_voxel.hpp"

/**
 * This class implements some percolation analysis to find the largest
 * sphere that can travel through the structures
 */
template <class T>
class percolation_analysis {

public:

  percolation_analysis( std::vector<T> data, std::vector<T> distance_map, unsigned int sx_, unsigned int sy_, unsigned int sz_, bool periodic);

  /// Using Hoshen-Kopelmann to find clusters
  void find_clusters();

  /// makes the cluster map a bit nicer
  void assign_true_labels();
  
  /// prints the cluster map to a file
  void print( std::ostream &out );
  void print( std::string out );

  ///
  std::vector<unsigned int> get_cluster_sizes() const;

  unsigned int get_nr_clusters() const;

  
  
private:

  void get_nbs( int id, int ch_id, std::vector<int> &nbs );
  
  std::vector<T> structure;
  std::vector<T> dmap; // have to think about the types here!  

  std::vector<int> cluster_labels;
  std::vector<int> cluster_sizes;
  


  /// size of the simbox
  unsigned int sx, sy, sz;
  /// is it periodic
  bool periodic;
  
  iterable_voxel it;
  
};





#endif
