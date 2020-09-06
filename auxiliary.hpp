
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
