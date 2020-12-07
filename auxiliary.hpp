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



#ifndef AUX_HPP
#define AUX_HPP

#include <vector>
#include <ostream>
#include <fstream>


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


};



#endif
