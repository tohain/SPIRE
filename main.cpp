#include "img_out.hpp"

#include <sstream>
#include <iomanip>
#include <cmath>

double dot_prod ( std::vector<double> v, std::vector<double> w ){
  return v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
}


struct Matrix {
  std::vector<double> v;
  std::vector<double> w;
  std::vector<double> z;  
};



std::vector<double> dot_prod( Matrix m, std::vector<double> v ){
  std::vector<double> r (3, 0);
  r[0] = dot_prod( m.v, v);
  r[1] = dot_prod( m.w, v);
  r[2] = dot_prod( m.z, v);
  return r;
}



void set_up_points( double slice_width, double slice_height, double L, double res, double theta, double phi,  std::vector<std::vector<double> > &points ){
  
  //vectors of slice
  std::vector<double> nx = {1, 0, 0};
  std::vector<double> ny = {0, 1, 0};
  std::vector<double> nz = {0, 0, 1};  

  //rotation matrices
  Matrix Rx, Rz;
  Rx.v = { 1, 0,           0          };
  Rx.w = { 0, cos(theta), -sin(theta) };
  Rx.z = { 0, sin(theta),  cos(theta) };

  Rz.v = { cos(phi), -sin(phi), 0 };
  Rz.w = { sin(phi),  cos(phi), 0 };
  Rz.z = { 0,         0,        1 };

  //rotate
  nx = dot_prod( Rz, dot_prod( Rx, nx));
  ny = dot_prod( Rz, dot_prod( Rx, ny));
  nz = dot_prod( Rz, dot_prod( Rx, nz));  
  
  //get number of points
  unsigned int n_points = int(L/res);
  unsigned int n_points_z = int(slice_width/res);

  for(unsigned int ii=0; ii<n_points; ii++){
    for(unsigned int jj=0; jj<n_points; jj++){
      for(unsigned int kk=0; kk<n_points_z; kk++){
        
	int ind=ii*(n_points*n_points_z) + jj*n_points_z + kk;
	
	double x=((ii*res)-L/2.)*nx[0] + ((jj*res)-L/2.)*ny[0] + ((kk*res)-L/2)*nz[0] + (slice_height+L/2.)*nz[0] - 0.5*slice_width;
	double y=((ii*res)-L/2.)*nx[1] + ((jj*res)-L/2.)*ny[1] + ((kk*res)-L/2)*nz[1] + (slice_height+L/2.)*nz[1] - 0.5*slice_width;
	double z=((ii*res)-L/2.)*nx[2] + ((jj*res)-L/2.)*ny[2] + ((kk*res)-L/2)*nz[2] + (slice_height+L/2.)*nz[2] - 0.5*slice_width;
	  
	points[ind][0]=x;
	points[ind][1]=y;
	points[ind][2]=z;
      }
    }
  }
  
}



double level_set_gyroid( double x, double y, double z, double a) {
  double s = 1./(2*M_PI*a);
  return sin(s*x)*cos(s*y) + sin(s*y)*cos(s*z) + cos(s*x)*sin(s*z);
}

double level_set_diamond( double x, double y, double z, double a) {
  double s = 1./(2*M_PI*a);
  return cos(s*x)*cos(s*y)*cos(s*z) - sin(s*x)*sin(s*y)*sin(s*z);
}


void set_grid ( std::string type, double d, double a, std::vector<std::vector<double> > &points, std::vector<int> &grid ){
  
  for( unsigned int ii=0; ii<points.size(); ii++ ){

    double level = 0;
    if( type == "gyroid" ){
      level = level_set_gyroid( points[ii][0], points[ii][1], points[ii][2], a );
    } else if( type == "diamond" ){
      level = level_set_diamond( points[ii][0], points[ii][1], points[ii][2], a );
    }


    if( level < d && level > -d ){
      grid[ii] = 1;
    }

  }

}

void get_projection (unsigned int n_points, unsigned int n_points_z, std::vector<int> &grid, std::vector< std::vector<double> > &projection ){
   
  for(unsigned int ii=0; ii<n_points; ii++){
    for(unsigned int jj=0; jj<n_points; jj++){
      for(unsigned int kk=0; kk<n_points_z; kk++){
        
	int ind=ii*(n_points*n_points_z) + jj*n_points_z + kk;

	if(grid[ind] == 1){
	  projection[ii][jj] += 1;
	}
	
      }
    }
  }

}

template <class T>
void set_to_zero( std::vector< std::vector<T> > &data ){
  for(unsigned int i=0; i<data.size(); i++){
    for(unsigned int j=0; j<data[i].size(); j++){
      data[i][j]=0;
    }
  }
}


template <class T>
void set_to_zero( std::vector<T> &data ){
  for(unsigned int j=0; j<data.size(); j++){
    data[j]=0;    
  }
}

int main( int argc, char* argv[] ){


  //set up parameters
  double L = 10;
  double res = .02;

  //slice dimension
  double slice_width = 1;
  double slice_height = 1;

  //membrane thickness
  double mem_width = 0.5;
  //gyroid unit cell length
  double a = 0.1;
  
  //slice orientation
  double theta = 0.25*M_PI;
  double phi = 0*M_PI;


  //get number of points
  unsigned int n_points = int(L/res);
  unsigned int n_points_z = int(slice_width/res);

  //array holding coordinates of the points
  std::vector<std::vector<double> > points ( n_points*n_points*n_points_z, std::vector<double> (3, 0) );

  std::vector<int> grid ( points.size(), 0 );

  //2d projection
  std::vector< std::vector<double> > projection (n_points, std::vector<double> (n_points, 0) );
  
  int slices = 50;
  double min=0,max=L;

  double step = (max-min)/(double)slices;
  
  for(int ii=0; ii<slices; ii++){

    std::cout << 100*ii/(float)slices << "\%\r" << std::flush;

    slice_height = step*ii;
    
    //get the points in the slice    
    set_up_points( slice_width, slice_height, L, res, theta, phi, points );

    //reset grid
    set_to_zero(grid);
    //get grid
    set_grid( "diamond", mem_width, a, points, grid );

    //get projection
    set_to_zero(projection);
    get_projection( n_points, n_points_z, grid, projection );
  
    //write image
    std::stringstream fn ("img_");
    fn << "img_" << std::setfill('0') << std::setw(3) << ii << ".pgm";
    write_image( fn.str(), projection );

  }
  
  return EXIT_SUCCESS;
}
