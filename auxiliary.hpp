/* Projection tool - compute planar projections of triply periodic
 * minimal surfaces 
 * Copyright (C) 2021 Tobias Hain
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

#ifndef AUX_HPP
#define AUX_HPP

#include <algorithm>
#include <cmath>
#include <vector>
#include <ostream>
#include <fstream>
#include <boost/random.hpp>


///General class to provide some useful methods
class my_utility {

public:

  /** Function to split a string at given characters
   * \param[in] in The string to splie
   * \param[in] del The character to split string at
   * \return A vector of strings with parts of the input string
   */
  static std::vector<std::string> str_split ( std::string in, char del ){

  std::string current_word;

  std::vector<std::string> out;
  
  for( auto it = in.begin(); it != in.end(); it++){
    if(*it == del){
      if( current_word != "" ){
	out.push_back( current_word );
	current_word = "";
      }
    } else {
      current_word += *it;
    }
  }
  //add last word
  if( current_word != "" )
    out.push_back( current_word );
  
  return out;
  
  }


  static void strip ( std::string &in, char c = ' '){

    while( *in.begin() == c ){
      in.erase( in.begin() );
    }

    while( *(--in.end()) == c ){
      in.erase( (--in.end()) );
    }    

  }
  


  /** euclid's algorithm to find the greatest common divisor of two
   * integers implemented as presented in
   *   
   * https://www.khanacademy.org/computing/computer-science/cryptography/modarithmetic/a/the-euclidean-algorithm
   */
  template <class T>
  static T gcd_euclid( T a, T b ){    
    T mod;
    while(1){
      if( a == 0 ){
	return b;
      }
      if( b == 0 ){
	return a;
      }
      mod = a % b;
      a = b;
      b = mod;
    }
  }  


  /**
   * a wrapper around above function to compute the gcd of a set of
   * integer of arbitrary cardinality
   */
  template <class T>
  static T gcd_euclid( std::vector<T> set ){

    if( set.size() == 0 ){
      return 0;
    }
    if( set.size() == 1 ){
      return set.at(0);
    }

    int pos = 0;
    T gcd = set.at(pos);
    while( pos < set.size() - 1 ){
      gcd = gcd_euclid( gcd, set.at(pos+1) );
      pos++;
    }

    return gcd;    
  }

  
};



/** \brief A tiny wrapper class around the boost random generator
 */
class boost_rand {

public:

  ///Constructor. Initialize boost objects to current seed
  boost_rand( int _seed ) : seed(_seed) {
    rng = boost::random::mt19937 ( seed );
    pdf = boost::random::uniform_real_distribution<> (0, 1);    
  }

  ///Returns uniformly distributed doubles between 0 and 1
  double rand_uniform(){
    return pdf(rng);
  }

  /** \brief Returns uniformly distributed doubles between lo and hi
   * \param[in] lo Lower boundary of distribution
   * \param[in] hi Upper boundary of distribution
   */  
  double rand_uniform(double lo, double hi ){

    double diff = hi - lo;
    return (rand_uniform()*diff)+lo;

  }

private:  
  boost::random::mt19937 rng;
  boost::random::uniform_real_distribution<> pdf;
  std::size_t seed;
};




class geometry {

public:
  
  /*
   * Computes the area of a polygon
   * \param[in] x the x values of the vertices
   * \param[in] y the y values of the vertices
   */
  static double compute_polygon_area( std::vector<double> x, std::vector<double> y){

    if( x.size() != y.size() ){
      throw std::string( "vertices arrays have different lengths" );
    }

    double area = 0;

    for(unsigned int ii=0; ii<x.size()-1; ii++){
      area += (x.at(ii+1) + x.at(ii))*(y.at(ii+1) - y.at(ii));
    }
    //do the last one manually due to "periodic boundary conditions;
    area += ( x.at( 0 ) + x.at( x.size() - 1) ) * ( y.at( 0 ) - y.at( y.size() - 1 ));
           
    return 0.5*area;    
  }

  /* computes the perimter (cirumference) of a polygon
   * \param[in] x the x values of the vertices
   * \param[in] y the y values of the vertices
   */
  static double compute_polygon_perimeter( std::vector<double> x, std::vector<double> y){

    if( x.size() != y.size() ){
      throw std::string( "vertices arrays have different lengths" );
    }

    double perimeter = 0;

    for(unsigned int ii=0; ii<x.size()-1; ii++){
      perimeter += sqrt( (x[ii+1] - x[ii])*(x[ii+1] - x[ii]) + ((y[ii+1] - y[ii])*(y[ii+1] - y[ii])) );
    }
    //do the last one manually due to "periodic boundary conditions;
    perimeter += sqrt( (( x.at(0) - x.at( x.size() - 1) )*( x.at(0) - x.at( x.size() - 1) ) ) + (( y.at(0) - y.at( y.size() - 1) )*( y.at(0) - y.at( y.size() - 1) ) ) );
           
    return perimeter;
  }
  

};


#endif
