/* SPIRE - Structure Projection Image Recognition Environment
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
