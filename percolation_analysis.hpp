/* Projection tool - compute planar projection of triply periodic
 * minimal surfaces 
 * Copyright (C) 2020 Tobias Hain
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see https://www.gnu.org/licenses.
 */



#ifndef PERCOLATION_ANALYSIS_H
#define PERCOLATION_ANALYSIS_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstring>  //memset
#include <unordered_set>
#include <set>

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
  void find_clusters( int ch_id_ );

  /// makes the cluster map a bit nicer
  void assign_true_labels();
  
  /// Checks, if clusters are percolating
  std::unordered_set<int> get_percolating_clusters( bool in_x, bool in_y, bool in_z);

  /// Get the percolation threshold
  M get_percolation_threshold( int ch_id );
  
  /// retrieves all the cluster sizes from \ref cluster_sizes
  std::vector<unsigned int> get_cluster_sizes() const;

  /// returns the number of clusters found
  unsigned int get_nr_clusters() const;

  /// prints the cluster map to a file
  void print( std::ostream &out );
  void print( std::string out );

  /// more for debugging
  void print_points( std::ostream &out );
  void print_points( std::string out );

  
  
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


  /// dilutes the surface on pixel
  void dilute_surface( M threshold, T marker = 1 );
  
  /// The ch_id of the latest cluster analysis
  int ch_id;

  /// size of the simbox
  unsigned int sx, sy, sz;
  /// is it periodic
  bool periodic;

  
  iterable_voxel it;  
};





#endif
