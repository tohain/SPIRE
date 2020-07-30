#include "surface_projection.hpp"
#include "distance_transform.hpp"
#include "img_out.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

void print_menue(){

  std::vector<std::string> options = {"SET_NTUCS", "SET_A", "SET_TYPE", "SET_LEVEL_SET", "EXIT"};

  for(unsigned int ii=0; ii<options.size(); ii++){
    std::cout << ii << ":   " << options.at(ii) << std::endl;;
  }
    

}


int main( int argc, char* argv[] ){



  SURFACE_TABLES tt;

  double level = tt.get_level( "Gyroid", 0.223 );
  double prop = tt.get_prop( "Gyroid", level );

  std::cout << level << " " << prop << std::endl;
  
  /*
  
  std::vector<double> fd;
  std::ifstream in ( argv[1] );

  double x,y,z;
  int i;
  while( !in.eof() ){
    in >> x >> y >> z >> i;
    fd.push_back( i );
  }
    

  int n = 25;
  int k = 25;
  int l = 1;
  std::vector<double> data ( n*k*l, 0 );


  //data[ (int(k/2)+1)*(n) + int(n/2) ] = 1;
  //data[ (int(k/2))*(n) + int(n/2) ] = 1;  

  //data[1]=1;
  //data[6]=1;
  //data[12]=1;
  
  data = fd;
  //data[382] = 0;
  
  for(unsigned int ii=0; ii<k; ii++){
    for(unsigned int jj=0; jj<n; jj++){
      std::cout << std::setw(3) << data[ii*n + jj] << " ";
      //std::cout << std::setw(3) << ii*n+jj << "|" << data[ii*n + jj] << " ";      
    }
    std::cout << std::endl;
  }
  
  
  distance_transform<double> dt ( data, n, k, l, 0.04f,  0.04f, 0.1f );
  
  dt.compute_distance_map();
  auto dmap = dt.get_distance_map();
  dt.print_map();


  
  write_image( "img.pgm", data, n, k );
  write_image( "edm.pgm", dmap, n, k );
  
  /*  
  int rad = 9;
  int width = 5;
  
  std::ofstream out ( "grid.dat" );
  //  n   k   l
  int q=0,w=0,e=0;
  for(int i=0; i<data.size(); i++){

    out << " " << q << " " << w << " " << e << " ";
    if( (rad*rad) - width < dmap[i] && dmap[i] < (rad*rad)+width  ){
      out << 1;
    } else {
      out << 0;
    }

    e++;
    if( e == l ) e = 0;

    if( i % l == (l-1) && i>0 ) {q++;}
    if( q == n ) q = 0;

    if( i % (n*l) == (n*l)-1 && i>0 ) {w++;}
    if( w == k ) w = 0;
			      
    out << std::endl;
  }
  
  */
  
  /*
  Do {
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
1	  std::cout << e.what() << ": " << e.details() << std::endl;
	}
	break;
	
    case 4:
	option = -1;
	break;
      
    }
  } while( option != -1  );

  */

  
  return EXIT_SUCCESS;
}
