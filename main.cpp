#include "img_out.hpp"

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



std::vector<std::vector<double> > set_up_points( double slice_width, double slice_height, double L, double res, double theta, double phi ){
  
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

  std::cout << nx[0] << " " << nx[1] << " " << nx[2] << std::endl;
  std::cout << ny[0] << " " << ny[1] << " " << ny[2] << std::endl;
  std::cout << nz[0] << " " << nz[1] << " " << nz[2] << std::endl;  
  
  //get number of points
  unsigned int n_points = int(L/res)+1;
  unsigned int n_points_z = int(slice_width/res);

    std::vector<std::vector<double> > points ( n_points*n_points*n_points_z, std::vector<double> (3, 0) );

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

  return points;
  
}



double level_set_gyroid( double x, double y, double z, double a) {
  double s = 1./(2*M_PI*a);
  return sin(s*x)*cos(s*y) + sin(s*y)*cos(s*z) + cos(s*x)*sin(s*z);
}


std::vector<int> set_grid ( std::string type, double d, double a, std::vector<std::vector<double> > &points ){

  std::vector<int> grid ( points.size(), 0 );
  
  for( unsigned int ii=0; ii<points.size(); ii++ ){
    
    double level = level_set_gyroid( points[ii][0], points[ii][1], points[ii][2], a );
    if( level < d && level > -d ){
      grid[ii] = 1;
    }

  }

  return grid;

}

std::vector< std::vector<double> > get_projection (unsigned int n_points, unsigned int n_points_z, std::vector<int> &grid ){

  std::vector< std::vector<double> > projection (n_points, std::vector<double> (n_points, 0) );
  
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

  return projection;
}


int main( int argc, char* argv[] ){


  //set up parameters
  double L = 10;
  double res = .05;

  //slice dimension
  double slice_width = 3;
  double slice_height = 1;
  
  //slice orientation
  double theta = 0;
  double phi = 0;


  //get number of points
  unsigned int n_points = int(L/res)+1;
  unsigned int n_points_z = int(slice_width/res);
  
  //get the points in the slice
  std::vector<std::vector<double> > points = set_up_points( slice_width, slice_height, L, res, theta, phi );

  //get grid
  auto grid = set_grid( "test", 0.3, .1, points );


  //get projection
  auto projection = get_projection( n_points, n_points_z, grid );
  
  //write image
  write_image( "img.pgm", projection );
  
  return EXIT_SUCCESS;
}
