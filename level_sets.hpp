#ifndef LEVEL_SETS
#define LEVEL_SETS

#include <cmath>
#include <vector>

namespace level_set {
  
  /**
   * Nodal representation of a Gyroid Ia\bar(3)d surface. From
   *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
   * (1991) 83: 407. https://doi.org/10.1007/BF01313411
   */
  double gyroid( double x, double y, double z, std::vector<double> inv_a);

  /**
   * Nodal representation of a Diamond Pn\bar(3)m surface. From
   *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
   * (1991) 83: 407. https://doi.org/10.1007/BF01313411
   */
  double diamond( double x, double y, double z, std::vector<double> inv_a);
    
  /**
   * Nodal representation of a Primitive Im\bar(3)m surface. From
   *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
   * (1991) 83: 407. https://doi.org/10.1007/BF01313411
   */
  double primitive( double x, double y, double z, std::vector<double> inv_a);
  
  /**
   * 
   */
  double lonsdaleite( double x, double y, double z, std::vector<double> inv_a);
  
  /**
   * This one is the formula Matthias Saba derived analytically.
   */
  double lonsdaleite_topo( double x, double y, double z, std::vector<double> inv_a);
      
  /**
   * level-set surface derived of the beta-m rod packing from Myfanwy Evans.
   */
  double beta_m( double x, double y, double z, std::vector<double> inv_a);

  /**
   * level-set surface derived of the sigma_plus rod packing from Myfanwy Evans.
   */
  double sigma_p( double x, double y, double z, std::vector<double> inv_a);

}


#endif
