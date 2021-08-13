
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


#include "level_sets.hpp"

/**
 * Nodal representation of a Gyroid Ia\bar(3)d surface. From
 *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
 * (1991) 83: 407. https://doi.org/10.1007/BF01313411
 */
double level_set::gyroid( double x, double y, double z, std::vector<double> inv_a) {
  return sin(inv_a[0]*x)*cos(inv_a[1]*y) + sin(inv_a[1]*y)*cos(inv_a[2]*z) + cos(inv_a[0]*x)*sin(inv_a[2]*z);
}


/**
 * Nodal representation of a Diamond Pn\bar(3)m surface. From
 *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
 * (1991) 83: 407. https://doi.org/10.1007/BF01313411
 */
double level_set::diamond( double x, double y, double z, std::vector<double> inv_a) {
  return cos(inv_a[0]*x)*cos(inv_a[1]*y)*cos(inv_a[2]*z) - sin(inv_a[0]*x)*sin(inv_a[1]*y)*sin(inv_a[2]*z);
}


/**
 * Nodal representation of a Primitive Im\bar(3)m surface. From
 *  Schnering, H.G. & Nesper, R. Z. Physik B - Condensed Matter
 * (1991) 83: 407. https://doi.org/10.1007/BF01313411
 */
double level_set::primitive( double x, double y, double z, std::vector<double> inv_a) {
  return cos(inv_a[0]*x)+cos(inv_a[1]*y)+cos(inv_a[2]*z);
}

/**
 * P.J.F. Gandy et al. / Chemical Physics Letters 336 (2001) 187-195
 */
double level_set::iwf( double x, double y, double z, std::vector<double> inv_a){

  return 2 * (cos(inv_a[0]*x)*cos(inv_a[1]*y)
	      + cos(inv_a[2]*z)*cos(inv_a[0]*x)
	      + cos(inv_a[1]*y)*cos(inv_a[2]*z) )
    - ( cos(2*inv_a[0]*x)
	+ cos(2*inv_a[1]*y)
	+ cos(2*inv_a[2]*z) );
}


/**
 *
 *
 */
double level_set::lonsdaleite( double x, double y, double z, std::vector<double> inv_a) {
  
  double freq_x = inv_a[0];
  double freq_y = inv_a[1];
  double freq_z = inv_a[2];      
  
    return + 0.0010124 - 0.360068*cos(2*+freq_z*z) + 0.221998*cos(2*+freq_y*y) + 0.225139*sin(2*+freq_y*y)*sin(+freq_z*z) + 0.445276*cos(freq_x*x)*cos(+freq_y*y) - 0.451854*cos(freq_x*x)*sin(+freq_y*y)*sin(+freq_z*z) + 0.211026*cos(freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z);
    
  }
  
  /*
   * This one is the formula Matthias Saba derived analytically.
   */
  double level_set::lonsdaleite_topo( double x, double y, double z, std::vector<double> inv_a){
    
    
    /* we need to adapt our x,y,z coordinates, since this series is
     * derived using hexgonal (not cubical) coordinates, so we need to
     * transform the into the canonical hex base
     */
    double x_ = inv_a[0] * (x + y/sqrt(3));
    double y_ = inv_a[1] * (2.0/sqrt(3)) * y;
    double z_ = inv_a[2] * z;  
    
    return  -cos(2*z_) + cos(x_) + cos(y_) + cos(x_ - y_) + sin(z_)*(-sin(x_) + sin(y_) + sin(x_ - y_) );
    
  }
  
  /**
   * level-set surface derived of the beta-m rod packing from Myfanwy Evans.
   */
  double level_set::beta_m( double x, double y, double z, std::vector<double> inv_a){
    
    double freq_x = inv_a[0];
    double freq_y = inv_a[1];
    double freq_z = inv_a[2];
    
    //
    // created this representation using a rod diameter of 0.05 in a
    // cubic box of edge length a=1. Surface voxelized with a resolution
    // of 64x64x64 voxels and included all terms above a threshold of
    // 0.02.
    //
    
    return  + 0.453125 + 0.0501726*cos(4*+freq_z*z) - 0.0227673*cos(8*+freq_z*z) + 0.0303563*cos(+freq_y*y)*cos(+freq_z*z) + 0.0275134*cos(+freq_y*y)*sin(+freq_z*z) - 0.0334931*sin(+freq_y*y)*cos(+freq_z*z) - 0.0303563*sin(+freq_y*y)*sin(+freq_z*z) - 0.0221057*cos(+freq_y*y)*cos(3*+freq_z*z) + 0.029806*cos(+freq_y*y)*sin(3*+freq_z*z) + 0.0243899*sin(+freq_y*y)*cos(3*+freq_z*z) - 0.0328859*sin(+freq_y*y)*sin(3*+freq_z*z) - 0.0257894*cos(+freq_y*y)*cos(5*+freq_z*z) + 0.0284542*sin(+freq_y*y)*cos(5*+freq_z*z) + 0.0210055*sin(+freq_y*y)*sin(7*+freq_z*z) + 0.0561277*sin(2*+freq_y*y)*sin(2*+freq_z*z) - 0.0474348*sin(2*+freq_y*y)*cos(4*+freq_z*z) - 0.0349778*sin(2*+freq_y*y)*sin(6*+freq_z*z) + 0.0212318*sin(2*+freq_y*y)*cos(8*+freq_z*z) - 0.0328859*cos(3*+freq_y*y)*cos(+freq_z*z) - 0.029806*cos(3*+freq_y*y)*sin(+freq_z*z) - 0.0243899*sin(3*+freq_y*y)*cos(+freq_z*z) - 0.0221057*sin(3*+freq_y*y)*sin(+freq_z*z) + 0.0238845*cos(3*+freq_y*y)*cos(3*+freq_z*z) - 0.0322045*cos(3*+freq_y*y)*sin(3*+freq_z*z) - 0.0238845*sin(3*+freq_y*y)*sin(3*+freq_z*z) + 0.0276904*cos(3*+freq_y*y)*cos(5*+freq_z*z) + 0.0205366*sin(3*+freq_y*y)*cos(5*+freq_z*z) + 0.0201749*cos(3*+freq_y*y)*sin(7*+freq_z*z) + 0.0501726*cos(4*+freq_y*y) + 0.0474348*cos(4*+freq_y*y)*sin(2*+freq_z*z) - 0.0398148*cos(4*+freq_y*y)*cos(4*+freq_z*z) - 0.0289311*cos(4*+freq_y*y)*sin(6*+freq_z*z) + 0.0284542*sin(5*+freq_y*y)*cos(+freq_z*z) + 0.0257894*sin(5*+freq_y*y)*sin(+freq_z*z) - 0.0205366*sin(5*+freq_y*y)*cos(3*+freq_z*z) + 0.0276904*sin(5*+freq_y*y)*sin(3*+freq_z*z) - 0.023452*sin(5*+freq_y*y)*cos(5*+freq_z*z) - 0.0349778*sin(6*+freq_y*y)*sin(2*+freq_z*z) + 0.0289311*sin(6*+freq_y*y)*cos(4*+freq_z*z) + 0.0203486*sin(6*+freq_y*y)*sin(6*+freq_z*z) + 0.0210055*cos(7*+freq_y*y)*cos(+freq_z*z) + 0.0201749*cos(7*+freq_y*y)*sin(3*+freq_z*z) - 0.0227673*cos(8*+freq_y*y) - 0.0212318*cos(8*+freq_y*y)*sin(2*+freq_z*z) + 0.0303563*cos(freq_x*x)*cos(+freq_z*z) - 0.0334931*cos(freq_x*x)*sin(+freq_z*z) + 0.0275134*sin(freq_x*x)*cos(+freq_z*z) - 0.0303563*sin(freq_x*x)*sin(+freq_z*z) - 0.0328859*cos(freq_x*x)*cos(3*+freq_z*z) - 0.0243899*cos(freq_x*x)*sin(3*+freq_z*z) - 0.029806*sin(freq_x*x)*cos(3*+freq_z*z) - 0.0221057*sin(freq_x*x)*sin(3*+freq_z*z) + 0.0284542*cos(freq_x*x)*sin(5*+freq_z*z) + 0.0257894*sin(freq_x*x)*sin(5*+freq_z*z) + 0.0210055*cos(freq_x*x)*cos(7*+freq_z*z) + 0.0303563*cos(freq_x*x)*cos(+freq_y*y) + 0.0275134*cos(freq_x*x)*sin(+freq_y*y) - 0.0334931*sin(freq_x*x)*cos(+freq_y*y) - 0.0303563*sin(freq_x*x)*sin(+freq_y*y) - 0.0221057*cos(freq_x*x)*cos(3*+freq_y*y) + 0.029806*cos(freq_x*x)*sin(3*+freq_y*y) + 0.0243899*sin(freq_x*x)*cos(3*+freq_y*y) - 0.0328859*sin(freq_x*x)*sin(3*+freq_y*y) - 0.0257894*cos(freq_x*x)*cos(5*+freq_y*y) + 0.0284542*sin(freq_x*x)*cos(5*+freq_y*y) + 0.0210055*sin(freq_x*x)*sin(7*+freq_y*y) + 0.0561277*sin(2*freq_x*x)*sin(2*+freq_z*z) + 0.0474348*sin(2*freq_x*x)*cos(4*+freq_z*z) - 0.0349778*sin(2*freq_x*x)*sin(6*+freq_z*z) - 0.0212318*sin(2*freq_x*x)*cos(8*+freq_z*z) + 0.0561277*sin(2*freq_x*x)*sin(2*+freq_y*y) - 0.0474348*sin(2*freq_x*x)*cos(4*+freq_y*y) - 0.0349778*sin(2*freq_x*x)*sin(6*+freq_y*y) + 0.0212318*sin(2*freq_x*x)*cos(8*+freq_y*y) - 0.0221057*cos(3*freq_x*x)*cos(+freq_z*z) + 0.0243899*cos(3*freq_x*x)*sin(+freq_z*z) + 0.029806*sin(3*freq_x*x)*cos(+freq_z*z) - 0.0328859*sin(3*freq_x*x)*sin(+freq_z*z) + 0.0238845*cos(3*freq_x*x)*cos(3*+freq_z*z) - 0.0322045*sin(3*freq_x*x)*cos(3*+freq_z*z) - 0.0238845*sin(3*freq_x*x)*sin(3*+freq_z*z) - 0.0205366*cos(3*freq_x*x)*sin(5*+freq_z*z) + 0.0276904*sin(3*freq_x*x)*sin(5*+freq_z*z) + 0.0201749*sin(3*freq_x*x)*cos(7*+freq_z*z) - 0.0328859*cos(3*freq_x*x)*cos(+freq_y*y) - 0.029806*cos(3*freq_x*x)*sin(+freq_y*y) - 0.0243899*sin(3*freq_x*x)*cos(+freq_y*y) - 0.0221057*sin(3*freq_x*x)*sin(+freq_y*y) + 0.0238845*cos(3*freq_x*x)*cos(3*+freq_y*y) - 0.0322045*cos(3*freq_x*x)*sin(3*+freq_y*y) - 0.0238845*sin(3*freq_x*x)*sin(3*+freq_y*y) + 0.0276904*cos(3*freq_x*x)*cos(5*+freq_y*y) + 0.0205366*sin(3*freq_x*x)*cos(5*+freq_y*y) + 0.0201749*cos(3*freq_x*x)*sin(7*+freq_y*y) + 0.0501726*cos(4*freq_x*x) - 0.0474348*cos(4*freq_x*x)*sin(2*+freq_z*z) - 0.0398148*cos(4*freq_x*x)*cos(4*+freq_z*z) + 0.0289311*cos(4*freq_x*x)*sin(6*+freq_z*z) + 0.0474348*cos(4*freq_x*x)*sin(2*+freq_y*y) - 0.0398148*cos(4*freq_x*x)*cos(4*+freq_y*y) - 0.0289311*cos(4*freq_x*x)*sin(6*+freq_y*y) - 0.0257894*cos(5*freq_x*x)*cos(+freq_z*z) + 0.0284542*cos(5*freq_x*x)*sin(+freq_z*z) + 0.0276904*cos(5*freq_x*x)*cos(3*+freq_z*z) + 0.0205366*cos(5*freq_x*x)*sin(3*+freq_z*z) - 0.023452*cos(5*freq_x*x)*sin(5*+freq_z*z) + 0.0284542*sin(5*freq_x*x)*cos(+freq_y*y) + 0.0257894*sin(5*freq_x*x)*sin(+freq_y*y) - 0.0205366*sin(5*freq_x*x)*cos(3*+freq_y*y) + 0.0276904*sin(5*freq_x*x)*sin(3*+freq_y*y) - 0.023452*sin(5*freq_x*x)*cos(5*+freq_y*y) - 0.0349778*sin(6*freq_x*x)*sin(2*+freq_z*z) - 0.0289311*sin(6*freq_x*x)*cos(4*+freq_z*z) + 0.0203486*sin(6*freq_x*x)*sin(6*+freq_z*z) - 0.0349778*sin(6*freq_x*x)*sin(2*+freq_y*y) + 0.0289311*sin(6*freq_x*x)*cos(4*+freq_y*y) + 0.0203486*sin(6*freq_x*x)*sin(6*+freq_y*y) + 0.0210055*sin(7*freq_x*x)*sin(+freq_z*z) + 0.0201749*sin(7*freq_x*x)*cos(3*+freq_z*z) + 0.0210055*cos(7*freq_x*x)*cos(+freq_y*y) + 0.0201749*cos(7*freq_x*x)*sin(3*+freq_y*y) - 0.0227673*cos(8*freq_x*x) + 0.0212318*cos(8*freq_x*x)*sin(2*+freq_z*z) - 0.0212318*cos(8*freq_x*x)*sin(2*+freq_y*y);
      
}

  /**
   * level-set surface derived of the sigma_plus rod packing from Myfanwy Evans.
   */
  double level_set::sigma_p( double x, double y, double z, std::vector<double> inv_a){
    
    double freq_x = inv_a[0];
    double freq_y = inv_a[1];
    double freq_z = inv_a[2];
    
    return  - 0.434654 - 0.053616*cos(+freq_y*y)*cos(+freq_z*z) - 0.0501469*cos(+freq_y*y)*sin(+freq_z*z) + 0.0605525*sin(+freq_y*y)*cos(+freq_z*z) + 0.054602*sin(+freq_y*y)*sin(+freq_z*z) - 0.0568908*sin(2*+freq_y*y)*sin(2*+freq_z*z) - 0.0368302*cos(4*+freq_y*y)*cos(4*+freq_z*z) - 0.035254*sin(5*+freq_y*y)*cos(5*+freq_z*z) - 0.0546025*cos(freq_x*x)*cos(+freq_z*z) + 0.0624045*cos(freq_x*x)*sin(+freq_z*z) - 0.0485864*sin(freq_x*x)*cos(+freq_z*z) + 0.0556542*sin(freq_x*x)*sin(+freq_z*z) - 0.0557864*cos(freq_x*x)*cos(+freq_y*y) - 0.0478919*cos(freq_x*x)*sin(+freq_y*y) + 0.0586734*sin(freq_x*x)*cos(+freq_y*y) + 0.0545809*sin(freq_x*x)*sin(+freq_y*y) + 0.0583669*cos(freq_x*x)*cos(+freq_y*y)*cos(2*+freq_z*z) - 0.0676233*cos(freq_x*x)*sin(+freq_y*y)*cos(2*+freq_z*z) + 0.0539172*sin(freq_x*x)*cos(+freq_y*y)*cos(2*+freq_z*z) - 0.0569372*sin(freq_x*x)*sin(+freq_y*y)*cos(2*+freq_z*z) + 0.0584311*cos(freq_x*x)*cos(2*+freq_y*y)*cos(+freq_z*z) + 0.0553288*cos(freq_x*x)*cos(2*+freq_y*y)*sin(+freq_z*z) - 0.0641567*sin(freq_x*x)*cos(2*+freq_y*y)*cos(+freq_z*z) - 0.0583956*sin(freq_x*x)*cos(2*+freq_y*y)*sin(+freq_z*z) + 0.0338027*cos(freq_x*x)*cos(2*+freq_y*y)*sin(3*+freq_z*z) + 0.05659*cos(freq_x*x)*sin(2*+freq_y*y)*cos(3*+freq_z*z) + 0.0450123*sin(freq_x*x)*sin(2*+freq_y*y)*cos(3*+freq_z*z) - 0.0303097*sin(freq_x*x)*cos(2*+freq_y*y)*sin(3*+freq_z*z) + 0.0369109*cos(freq_x*x)*sin(2*+freq_y*y)*sin(3*+freq_z*z) + 0.0382261*sin(freq_x*x)*sin(2*+freq_y*y)*sin(3*+freq_z*z) - 0.0324769*cos(freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z) - 0.033186*cos(freq_x*x)*cos(3*+freq_y*y)*sin(2*+freq_z*z) + 0.0420654*sin(freq_x*x)*cos(3*+freq_y*y)*sin(2*+freq_z*z) + 0.0531061*cos(freq_x*x)*sin(3*+freq_y*y)*sin(2*+freq_z*z) - 0.0513322*sin(freq_x*x)*sin(3*+freq_y*y)*sin(2*+freq_z*z) - 0.0315818*cos(freq_x*x)*cos(3*+freq_y*y)*sin(4*+freq_z*z) - 0.0421543*cos(freq_x*x)*sin(3*+freq_y*y)*cos(4*+freq_z*z) - 0.0372372*sin(freq_x*x)*cos(3*+freq_y*y)*cos(4*+freq_z*z) + 0.042761*sin(freq_x*x)*sin(3*+freq_y*y)*cos(4*+freq_z*z) + 0.0453281*cos(freq_x*x)*cos(4*+freq_y*y)*cos(3*+freq_z*z) + 0.0359237*sin(freq_x*x)*cos(4*+freq_y*y)*cos(3*+freq_z*z) + 0.0351133*sin(freq_x*x)*cos(4*+freq_y*y)*sin(3*+freq_z*z) + 0.0349508*cos(freq_x*x)*sin(4*+freq_y*y)*cos(5*+freq_z*z) - 0.0422584*sin(freq_x*x)*sin(4*+freq_y*y)*cos(5*+freq_z*z) + 0.0385901*cos(freq_x*x)*sin(5*+freq_y*y)*sin(4*+freq_z*z) + 0.0381063*sin(freq_x*x)*sin(5*+freq_y*y)*sin(4*+freq_z*z) - 0.0604244*sin(2*freq_x*x)*sin(2*+freq_z*z) + 0.0578792*cos(2*freq_x*x)*cos(+freq_y*y)*cos(+freq_z*z) - 0.0669029*cos(2*freq_x*x)*cos(+freq_y*y)*sin(+freq_z*z) + 0.0518448*cos(2*freq_x*x)*sin(+freq_y*y)*cos(+freq_z*z) - 0.0598102*cos(2*freq_x*x)*sin(+freq_y*y)*sin(+freq_z*z) - 0.0339156*cos(2*freq_x*x)*cos(+freq_y*y)*cos(3*+freq_z*z) - 0.0319576*sin(2*freq_x*x)*cos(+freq_y*y)*cos(3*+freq_z*z) + 0.0420931*sin(2*freq_x*x)*sin(+freq_y*y)*cos(3*+freq_z*z) + 0.0538362*sin(2*freq_x*x)*cos(+freq_y*y)*sin(3*+freq_z*z) - 0.054398*sin(2*freq_x*x)*sin(+freq_y*y)*sin(3*+freq_z*z) - 0.0543821*sin(2*freq_x*x)*sin(2*+freq_y*y) - 0.0864182*sin(2*freq_x*x)*sin(2*+freq_y*y)*cos(4*+freq_z*z) + 0.0345527*cos(2*freq_x*x)*sin(3*+freq_y*y)*cos(+freq_z*z) + 0.0547613*sin(2*freq_x*x)*cos(3*+freq_y*y)*cos(+freq_z*z) + 0.0350496*sin(2*freq_x*x)*sin(3*+freq_y*y)*cos(+freq_z*z) + 0.0471496*sin(2*freq_x*x)*cos(3*+freq_y*y)*sin(+freq_z*z) + 0.0382912*sin(2*freq_x*x)*sin(3*+freq_y*y)*sin(+freq_z*z) + 0.0433837*sin(2*freq_x*x)*cos(3*+freq_y*y)*cos(5*+freq_z*z) - 0.0855712*sin(2*freq_x*x)*cos(4*+freq_y*y)*sin(2*+freq_z*z) - 0.0390934*cos(2*freq_x*x)*sin(4*+freq_y*y)*cos(6*+freq_z*z) + 0.0331577*sin(2*freq_x*x)*sin(5*+freq_y*y)*cos(3*+freq_z*z) - 0.0402394*sin(2*freq_x*x)*sin(5*+freq_y*y)*sin(3*+freq_z*z) + 0.0371987*cos(2*freq_x*x)*cos(6*+freq_y*y)*sin(4*+freq_z*z) + 0.0584293*cos(3*freq_x*x)*cos(+freq_y*y)*sin(2*+freq_z*z) + 0.0324406*sin(3*freq_x*x)*cos(+freq_y*y)*cos(2*+freq_z*z) + 0.0376846*sin(3*freq_x*x)*cos(+freq_y*y)*sin(2*+freq_z*z) + 0.0468598*cos(3*freq_x*x)*sin(+freq_y*y)*sin(2*+freq_z*z) + 0.0401109*sin(3*freq_x*x)*sin(+freq_y*y)*sin(2*+freq_z*z) + 0.0478707*cos(3*freq_x*x)*cos(+freq_y*y)*cos(4*+freq_z*z) + 0.0343564*cos(3*freq_x*x)*sin(+freq_y*y)*cos(4*+freq_z*z) + 0.0340071*sin(3*freq_x*x)*sin(+freq_y*y)*cos(4*+freq_z*z) - 0.0328556*sin(3*freq_x*x)*cos(+freq_y*y)*sin(4*+freq_z*z) - 0.0384777*cos(3*freq_x*x)*cos(2*+freq_y*y)*cos(+freq_z*z) - 0.0302012*cos(3*freq_x*x)*sin(2*+freq_y*y)*cos(+freq_z*z) + 0.0490262*sin(3*freq_x*x)*sin(2*+freq_y*y)*cos(+freq_z*z) + 0.0433543*cos(3*freq_x*x)*sin(2*+freq_y*y)*sin(+freq_z*z) - 0.0516353*sin(3*freq_x*x)*sin(2*+freq_y*y)*sin(+freq_z*z) + 0.0326691*cos(3*freq_x*x)*sin(2*+freq_y*y)*sin(5*+freq_z*z) - 0.0445752*sin(3*freq_x*x)*sin(2*+freq_y*y)*sin(5*+freq_z*z) - 0.0331402*sin(3*freq_x*x)*cos(3*+freq_y*y)*cos(6*+freq_z*z) - 0.0360347*cos(3*freq_x*x)*cos(4*+freq_y*y)*sin(+freq_z*z) - 0.0456296*sin(3*freq_x*x)*cos(4*+freq_y*y)*cos(+freq_z*z) + 0.0376895*sin(3*freq_x*x)*cos(4*+freq_y*y)*sin(+freq_z*z) + 0.039363*cos(3*freq_x*x)*cos(5*+freq_y*y)*sin(2*+freq_z*z) + 0.0306598*sin(3*freq_x*x)*cos(5*+freq_y*y)*sin(2*+freq_z*z) - 0.0324927*cos(4*freq_x*x)*cos(4*+freq_z*z) - 0.0439643*cos(4*freq_x*x)*cos(+freq_y*y)*sin(3*+freq_z*z) - 0.0348071*cos(4*freq_x*x)*sin(+freq_y*y)*cos(3*+freq_z*z) - 0.0356451*sin(4*freq_x*x)*cos(+freq_y*y)*cos(3*+freq_z*z) + 0.0404075*cos(4*freq_x*x)*sin(+freq_y*y)*sin(3*+freq_z*z) + 0.0424661*sin(4*freq_x*x)*cos(+freq_y*y)*sin(5*+freq_z*z) + 0.0378461*sin(4*freq_x*x)*sin(+freq_y*y)*sin(5*+freq_z*z) - 0.0864813*cos(4*freq_x*x)*sin(2*+freq_y*y)*sin(2*+freq_z*z) + 0.0373625*sin(4*freq_x*x)*cos(2*+freq_y*y)*cos(6*+freq_z*z) + 0.0506671*cos(4*freq_x*x)*cos(3*+freq_y*y)*cos(+freq_z*z) + 0.0349921*cos(4*freq_x*x)*cos(3*+freq_y*y)*sin(+freq_z*z) + 0.0305694*cos(4*freq_x*x)*sin(3*+freq_y*y)*cos(+freq_z*z) + 0.0318146*cos(4*freq_x*x)*sin(3*+freq_y*y)*sin(+freq_z*z) - 0.0400586*cos(4*freq_x*x)*cos(4*+freq_y*y) + 0.0380193*sin(4*freq_x*x)*cos(5*+freq_y*y)*cos(+freq_z*z) - 0.0383586*sin(4*freq_x*x)*cos(5*+freq_y*y)*sin(+freq_z*z) - 0.0387042*sin(4*freq_x*x)*cos(6*+freq_y*y)*cos(2*+freq_z*z) - 0.0359808*cos(5*freq_x*x)*sin(5*+freq_z*z) + 0.0374416*cos(5*freq_x*x)*cos(+freq_y*y)*sin(4*+freq_z*z) - 0.0409826*cos(5*freq_x*x)*sin(+freq_y*y)*sin(4*+freq_z*z) + 0.0467104*cos(5*freq_x*x)*sin(2*+freq_y*y)*cos(3*+freq_z*z) - 0.0301981*cos(5*freq_x*x)*cos(3*+freq_y*y)*cos(2*+freq_z*z) + 0.0316691*sin(5*freq_x*x)*cos(3*+freq_y*y)*sin(2*+freq_z*z) - 0.0381817*sin(5*freq_x*x)*sin(3*+freq_y*y)*sin(2*+freq_z*z) + 0.0421826*sin(5*freq_x*x)*sin(4*+freq_y*y)*cos(+freq_z*z) + 0.0336355*sin(5*freq_x*x)*sin(4*+freq_y*y)*sin(+freq_z*z) - 0.0362188*sin(5*freq_x*x)*cos(5*+freq_y*y) + 0.0304612*sin(6*freq_x*x)*sin(6*+freq_z*z) - 0.0363384*cos(6*freq_x*x)*cos(2*+freq_y*y)*sin(4*+freq_z*z) - 0.033385*cos(6*freq_x*x)*sin(3*+freq_y*y)*cos(3*+freq_z*z) + 0.0411412*cos(6*freq_x*x)*sin(4*+freq_y*y)*cos(2*+freq_z*z);

  }
