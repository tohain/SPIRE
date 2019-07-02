#include "surface_projection.hpp"

#include <iostream>
#include <string>
#include <vector>

void print_menue(){

  std::vector<std::string> options = {"SET_NTUCS", "SET_A", "SET_TYPE", "SET_LEVEL_SET", "EXIT"};

  for(unsigned int ii=0; ii<options.size(); ii++){
    std::cout << ii << ":   " << options.at(ii) << std::endl;;
  }
    

}


int main( int argc, char* argv[] ){

  surface_projection sp;

  int option;

  do {
    print_menue();
    std::cin >> option;

    double d_val;
    int i_val;
    
    switch( option ){

    case 0:
	std::cout << "\nntucs=";	
	std::cin >> i_val;
	try {
	  sp.set_ntucs( i_val );
	} catch (invalid_parameter_exception e ){
	  std::cout << e.what() << ": " << e.details() << std::endl;
	}
	break;
      
	
    case 1:
	std::cout << "\na=";
	std::cin >> d_val;
	try {
	  sp.set_a( d_val );
	} catch (invalid_parameter_exception e ){
	  std::cout << e.what() << ": " << e.details() << std::endl;
	}
	break;
	
    case 2:
	std::cout << "\ntype=";
	std::cin >> i_val;
	try {
	  sp.set_type( i_val );
	} catch (invalid_parameter_exception e ){
	  std::cout << e.what() << ": " << e.details() << std::endl;
	}
	break;

    case 3:
	std::cout << "\ntype=";
	std::cin >> d_val;
	try {
	  sp.set_mem_width( d_val );
	} catch (invalid_parameter_exception e ){
	  std::cout << e.what() << ": " << e.details() << std::endl;
	}
	break;
	
    case 4:
	option = -1;
	break;
      
    }
  } while( option != -1  );

  
  return EXIT_SUCCESS;
}
