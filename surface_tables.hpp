#ifndef SURFACE_TABLES_H
#define SURFACE_TABLES_H


#include <map>

/* \brief This class holds ready to use values for some parameters the
 * built in surfaces.
 *
 *
 * This class is primarily there to store some basic properties of the
 * basic surfaces built in the system. So far it is a lookup table
 * surface_level_set <--> channel volume proportions
 *
 */
class SURFACE_TABLES {

public:
  SURFACE_TABLES(){}

  const double get_level( std::string surface, double val ) const {

    if( val > 1 || val < 0 ){
      throw invalid_parameter_exception( "Invalid volume proportions, must be in [0,1]" );
    }


    std::map<double, double>::const_iterator it;

    if( surface == "Primitive" ){
      it = P_SURFACE_VOL.begin();
    } else if ( surface == "Diamond" ){
      it = D_SURFACE_VOL.begin();
    } else if ( surface == "Gyroid" ){
      it = G_SURFACE_VOL.begin();
    } else {
      throw invalid_parameter_exception( "Unkown surface type" );
    }
   
    
    while( it->first < val && it != P_SURFACE_VOL.end() )
      it++;

    //leaves the iterator at the upper border

    if( it == P_SURFACE_VOL.begin() ){
      return it->second;
    } else if ( it == (--P_SURFACE_VOL.end()) ){
      return (P_SURFACE_VOL.rbegin())->second;
    } else {

      double hi_v = it->second, hi_p = it->first;
      it--;
      double lo_v = it->second, lo_p = it->first;

      double fraction = (val - lo_p) / (hi_p - lo_p);
      return lo_v + fraction * (hi_v - lo_v);
    }
   
  }




  const double get_prop( std::string surface, double val ) const {

    std::map<double, double>::const_iterator it;

    if( surface == "Primitive" ){
      it = P_SURFACE_VOL_INV.begin();
    } else if ( surface == "Diamond" ){
      it = D_SURFACE_VOL_INV.begin();
    } else if ( surface == "Gyroid" ){
      it = G_SURFACE_VOL_INV.begin();
    } else {
      throw invalid_parameter_exception( "Unkown surface type" );
    }
   
    
    while( it->first < val && it != P_SURFACE_VOL.end() )
      it++;

    //leaves the iterator at the upper border

    if( it == P_SURFACE_VOL.begin() ){
      return it->second;
    } else if ( it == (--P_SURFACE_VOL.end()) ){
      return (P_SURFACE_VOL.rbegin())->second;
    } else {

      double hi_v = it->second, hi_p = it->first;
      it--;
      double lo_v = it->second, lo_p = it->first;

      double fraction = (val - lo_p) / (hi_p - lo_p);
      return lo_v + fraction * (hi_v - lo_v);
    }
   
  }  
  
private:
  
  const std::map<double, double> P_SURFACE_VOL_INV = {{-3.10000e+00, 1.00000e+00}, {-3.00000e+00, 9.99728e-01}, {-2.90000e+00, 9.97256e-01}, {-2.80000e+00, 9.93584e-01}, {-2.70000e+00, 9.88966e-01}, {-2.60000e+00, 9.83501e-01}, {-2.50000e+00, 9.77357e-01}, {-2.40000e+00, 9.70413e-01}, {-2.30000e+00, 9.62674e-01}, {-2.20000e+00, 9.54176e-01}, {-2.10000e+00, 9.44801e-01}, {-2.00000e+00, 9.34763e-01}, {-1.90000e+00, 9.23607e-01}, {-1.80000e+00, 9.11672e-01}, {-1.70000e+00, 8.98637e-01}, {-1.60000e+00, 8.84633e-01}, {-1.50000e+00, 8.69219e-01}, {-1.40000e+00, 8.52717e-01}, {-1.30000e+00, 8.34405e-01}, {-1.20000e+00, 8.14189e-01}, {-1.10000e+00, 7.91519e-01}, {-1.00000e+00, 7.64265e-01}, {-9.00000e-01, 7.34592e-01}, {-8.00000e-01, 7.04732e-01}, {-7.00000e-01, 6.75816e-01}, {-6.00000e-01, 6.46495e-01}, {-5.00000e-01, 6.17794e-01}, {-4.00000e-01, 5.88631e-01}, {-3.00000e-01, 5.60364e-01}, {-2.00000e-01, 5.31192e-01}, {-1.00000e-01, 5.02752e-01}, {1.52656e-15, 4.74859e-01}, {1.00000e-01, 4.45855e-01}, {2.00000e-01, 4.17641e-01}, {3.00000e-01, 3.88896e-01}, {4.00000e-01, 3.60851e-01}, {5.00000e-01, 3.32371e-01}, {6.00000e-01, 3.04092e-01}, {7.00000e-01, 2.75808e-01}, {8.00000e-01, 2.47735e-01}, {9.00000e-01, 2.19345e-01}, {1.00000e+00, 1.92730e-01}, {1.10000e+00, 1.70489e-01}, {1.20000e+00, 1.51406e-01}, {1.30000e+00, 1.34150e-01}, {1.40000e+00, 1.18433e-01}, {1.50000e+00, 1.04335e-01}, {1.60000e+00, 9.11997e-02}, {1.70000e+00, 7.93502e-02}, {1.80000e+00, 6.82172e-02}, {1.90000e+00, 5.82166e-02}, {2.00000e+00, 4.89804e-02}, {2.10000e+00, 4.06279e-02}, {2.20000e+00, 3.29893e-02}, {2.30000e+00, 2.60661e-02}, {2.40000e+00, 1.99517e-02}, {2.50000e+00, 1.44969e-02}, {2.60000e+00, 9.83615e-03}, {2.70000e+00, 5.93630e-03}, {2.80000e+00, 2.82400e-03}, {2.90000e+00, 7.17333e-04}, {3.00000e+00, 0.00000e+00}};

  const std::map<double, double> P_SURFACE_VOL = {{1.00000e+00, -3.10000e+00}, {9.99728e-01, -3.00000e+00}, {9.97256e-01, -2.90000e+00}, {9.93584e-01, -2.80000e+00}, {9.88966e-01, -2.70000e+00}, {9.83501e-01, -2.60000e+00}, {9.77357e-01, -2.50000e+00}, {9.70413e-01, -2.40000e+00}, {9.62674e-01, -2.30000e+00}, {9.54176e-01, -2.20000e+00}, {9.44801e-01, -2.10000e+00}, {9.34763e-01, -2.00000e+00}, {9.23607e-01, -1.90000e+00}, {9.11672e-01, -1.80000e+00}, {8.98637e-01, -1.70000e+00}, {8.84633e-01, -1.60000e+00}, {8.69219e-01, -1.50000e+00}, {8.52717e-01, -1.40000e+00}, {8.34405e-01, -1.30000e+00}, {8.14189e-01, -1.20000e+00}, {7.91519e-01, -1.10000e+00}, {7.64265e-01, -1.00000e+00}, {7.34592e-01, -9.00000e-01}, {7.04732e-01, -8.00000e-01}, {6.75816e-01, -7.00000e-01}, {6.46495e-01, -6.00000e-01}, {6.17794e-01, -5.00000e-01}, {5.88631e-01, -4.00000e-01}, {5.60364e-01, -3.00000e-01}, {5.31192e-01, -2.00000e-01}, {5.02752e-01, -1.00000e-01}, {4.74859e-01, 1.52656e-15}, {4.45855e-01, 1.00000e-01}, {4.17641e-01, 2.00000e-01}, {3.88896e-01, 3.00000e-01}, {3.60851e-01, 4.00000e-01}, {3.32371e-01, 5.00000e-01}, {3.04092e-01, 6.00000e-01}, {2.75808e-01, 7.00000e-01}, {2.47735e-01, 8.00000e-01}, {2.19345e-01, 9.00000e-01}, {1.92730e-01, 1.00000e+00}, {1.70489e-01, 1.10000e+00}, {1.51406e-01, 1.20000e+00}, {1.34150e-01, 1.30000e+00}, {1.18433e-01, 1.40000e+00}, {1.04335e-01, 1.50000e+00}, {9.11997e-02, 1.60000e+00}, {7.93502e-02, 1.70000e+00}, {6.82172e-02, 1.80000e+00}, {5.82166e-02, 1.90000e+00}, {4.89804e-02, 2.00000e+00}, {4.06279e-02, 2.10000e+00}, {3.29893e-02, 2.20000e+00}, {2.60661e-02, 2.30000e+00}, {1.99517e-02, 2.40000e+00}, {1.44969e-02, 2.50000e+00}, {9.83615e-03, 2.60000e+00}, {5.93630e-03, 2.70000e+00}, {2.82400e-03, 2.80000e+00}, {7.17333e-04, 2.90000e+00}, {0.00000e+00, 3.00000e+00}};

const std::map<double, double> D_SURFACE_VOL = {{1.00000e+00, -1.10000e+00}, {9.97783e-01, -1.00000e+00}, {9.76542e-01, -9.00000e-01}, {9.40908e-01, -8.00000e-01}, {8.83144e-01, -7.00000e-01}, {8.16670e-01, -6.00000e-01}, {7.53381e-01, -5.00000e-01}, {6.91841e-01, -4.00000e-01}, {6.31743e-01, -3.00000e-01}, {5.72569e-01, -2.00000e-01}, {5.14204e-01, -1.00000e-01}, {4.56372e-01, 1.52656e-15}, {3.98769e-01, 1.00000e-01}, {3.41266e-01, 2.00000e-01}, {2.84252e-01, 3.00000e-01}, {2.27486e-01, 4.00000e-01}, {1.70099e-01, 5.00000e-01}, {1.13509e-01, 6.00000e-01}, {5.93529e-02, 7.00000e-01}, {2.55656e-02, 8.00000e-01}, {5.85126e-03, 9.00000e-01}, {0.00000e+00, 1.00000e+00}};
const std::map<double, double> D_SURFACE_VOL_INV = {{-1.10000e+00, 1.00000e+00}, {-1.00000e+00, 9.97783e-01}, {-9.00000e-01, 9.76542e-01}, {-8.00000e-01, 9.40908e-01}, {-7.00000e-01, 8.83144e-01}, {-6.00000e-01, 8.16670e-01}, {-5.00000e-01, 7.53381e-01}, {-4.00000e-01, 6.91841e-01}, {-3.00000e-01, 6.31743e-01}, {-2.00000e-01, 5.72569e-01}, {-1.00000e-01, 5.14204e-01}, {1.52656e-15, 4.56372e-01}, {1.00000e-01, 3.98769e-01}, {2.00000e-01, 3.41266e-01}, {3.00000e-01, 2.84252e-01}, {4.00000e-01, 2.27486e-01}, {5.00000e-01, 1.70099e-01}, {6.00000e-01, 1.13509e-01}, {7.00000e-01, 5.93529e-02}, {8.00000e-01, 2.55656e-02}, {9.00000e-01, 5.85126e-03}, {1.00000e+00, 0.00000e+00}};

  const std::map<double, double> G_SURFACE_VOL = {{1.00000e+00, -1.50000e+00}, {9.77529e-01, -1.40000e+00}, {9.38475e-01, -1.30000e+00}, {9.01056e-01, -1.20000e+00}, {8.64996e-01, -1.10000e+00}, {8.29777e-01, -1.00000e+00}, {7.95304e-01, -9.00000e-01}, {7.61223e-01, -8.00000e-01}, {7.27761e-01, -7.00000e-01}, {6.94590e-01, -6.00000e-01}, {6.61833e-01, -5.00000e-01}, {6.29167e-01, -4.00000e-01}, {5.96772e-01, -3.00000e-01}, {5.64520e-01, -2.00000e-01}, {5.32018e-01, -1.00000e-01}, {5.00000e-01, 1.52656e-15}, {4.67982e-01, 1.00000e-01}, {4.35480e-01, 2.00000e-01}, {4.03228e-01, 3.00000e-01}, {3.70833e-01, 4.00000e-01}, {3.38167e-01, 5.00000e-01}, {3.05410e-01, 6.00000e-01}, {2.72239e-01, 7.00000e-01}, {2.38777e-01, 8.00000e-01}, {2.04696e-01, 9.00000e-01}, {1.70223e-01, 1.00000e+00}, {1.35004e-01, 1.10000e+00}, {9.89440e-02, 1.20000e+00}, {6.15253e-02, 1.30000e+00}, {2.24711e-02, 1.40000e+00}, {0.00000e+00, 1.50000e+00}};

  const std::map<double, double> G_SURFACE_VOL_INV = {{-1.50000e+00, 1.00000e+00}, {-1.40000e+00, 9.77529e-01}, {-1.30000e+00, 9.38475e-01}, {-1.20000e+00, 9.01056e-01}, {-1.10000e+00, 8.64996e-01}, {-1.00000e+00, 8.29777e-01}, {-9.00000e-01, 7.95304e-01}, {-8.00000e-01, 7.61223e-01}, {-7.00000e-01, 7.27761e-01}, {-6.00000e-01, 6.94590e-01}, {-5.00000e-01, 6.61833e-01}, {-4.00000e-01, 6.29167e-01}, {-3.00000e-01, 5.96772e-01}, {-2.00000e-01, 5.64520e-01}, {-1.00000e-01, 5.32018e-01}, {1.52656e-15, 5.00000e-01}, {1.00000e-01, 4.67982e-01}, {2.00000e-01, 4.35480e-01}, {3.00000e-01, 4.03228e-01}, {4.00000e-01, 3.70833e-01}, {5.00000e-01, 3.38167e-01}, {6.00000e-01, 3.05410e-01}, {7.00000e-01, 2.72239e-01}, {8.00000e-01, 2.38777e-01}, {9.00000e-01, 2.04696e-01}, {1.00000e+00, 1.70223e-01}, {1.10000e+00, 1.35004e-01}, {1.20000e+00, 9.89440e-02}, {1.30000e+00, 6.15253e-02}, {1.40000e+00, 2.24711e-02}, {1.50000e+00, 0.00000e+00}};
  

};





#endif
