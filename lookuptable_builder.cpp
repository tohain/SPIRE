#include <iostream>
#include <fstream>
#include <algorithm>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/tuple.h>
#include <CGAL/Vector_3.h>

#include "surface_projection.hpp"
#include "img_out.hpp"

//set double precision as the precision for the code
typedef CGAL::Simple_cartesian<double> K;

// type def to define a point
typedef K::Point_3  Point_3;
typedef K::Vector_3 Vector_3;

//A facet is an array of 3 objects of type size_t (some integer like data type)
typedef std::array<std::size_t,3> Facet;

// overload the output operator for facets 
namespace std {
  std::ostream&
  operator<<(std::ostream& os, const Facet& f)
  {
    os << "3 " << f[0] << " " << f[1] << " " << f[2];
    return os;
  }
}

// define a struct call perimeter
struct Perimeter {
  double bound;
  
  Perimeter(double bound): bound(bound){
  }


  // overload the () operator of this class making it a Functor
  template <typename AdvancingFront, typename Cell_handle>
  double operator() (const AdvancingFront& adv, Cell_handle& c,
                     const int& index) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return adv.smallest_radius_delaunay_sphere (c, index);
    }
    
    // If perimeter > bound, return infinity so that facet is not used
    double d  = 0;
    d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                              c->vertex((index+2)%4)->point()));
    if(d>bound) return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();
    // Otherwise, return usual priority value: smallest radius of
    // delaunay sphere
    return adv.smallest_radius_delaunay_sphere (c, index);
  }
};



int main(int argc, char* argv[])
{


  //get the surface projection done
  double progress;// = new double;
  std::string status;// = new char[20];

  
  surface_projection sp (progress, status);
  sp.set_slice_width( 1 );
  sp.set_slice_height( 0 );
  sp.set_n_points_x( std::stoi( argv[1] ) );
  sp.set_n_points_y( std::stoi( argv[1] ) );
  sp.set_n_points_z( std::stoi( argv[2] ) );
  sp.set_type( 0 );
  sp.set_surface_level( 0.0 ) ;
  //sp.set_ntucs(1);

  sp.update_containers();
  sp.update_geometry();

  std::cerr << "compute projection" << std::endl;

  sp.compute_projection();

  //std::cerr << "writing image" << std::endl;
  //write_image( "projection.ppm", sp.get_projection(), sp.get_width(), sp.get_height(), true );

  //sp.print_grid( "grid.dat");
  

  std::cerr << "finding relevant surface points\n";
  
  //write grid
  std::vector<double> membrane_points;
  //auto grid = sp.get_grid();
  auto grid_pos = sp.get_points( );
  auto interface_points = sp.get_surface_points( std::stoi( argv[3] ), 6 );
  std::cerr << "found " << interface_points.size() << " points\n";

  std::ofstream ip ("interface_points.dat");
  
  double per = std::stod( argv[4] );

  //copy points
  std::vector<Point_3> points;
  for(int ii=0; ii<interface_points.size(); ii++){
    double x = grid_pos.at( 3 * interface_points.at(ii) );
    double y = grid_pos.at( 3 * interface_points.at(ii) + 1 );
    double z = grid_pos.at( 3 * interface_points.at(ii) + 2 );
    points.push_back( Point_3( x, y, z) );
    ip << x << " " << y << " " << z << std::endl;
  }


  std::cerr << "computing surfce\n";

  std::vector<Facet> facets;
  Perimeter perimeter(per);
  CGAL::advancing_front_surface_reconstruction(points.begin(),
                                               points.end(),
                                               std::back_inserter(facets),
                                               perimeter);




  //compute surface area
  double total_area = 0;
  
  for( auto it : facets ){
    Vector_3 ab = points.at( it[0] ) - points.at( it[1] );
    Vector_3 ac = points.at( it[0] ) - points.at( it[2] );    
    Vector_3 cross = CGAL::cross_product( ab, ac );
    double area = 0.5 * sqrt( cross.squared_length() );
    total_area+=area;
  }

  std::cerr << "total_area=" << total_area << std::endl;
  sp.compute_volume();
  auto vols = sp.get_channel_volumes();
  for( auto it : vols ){
    std::cerr << "vol=" << it << std::endl;
  }
  std::cerr << std::endl;
    

  std::cout << "OFF\n" << points.size() << " " << facets.size() << " 0\n";
  std::copy(points.begin(),
            points.end(),
            std::ostream_iterator<Point_3>(std::cout, "\n"));
  std::copy(facets.begin(),
            facets.end(),
            std::ostream_iterator<Facet>(std::cout, "\n"));
  
  return 0;
}
