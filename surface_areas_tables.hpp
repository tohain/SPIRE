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

#ifndef SURFACE_TABLES_H
#define SURFACE_TABLES_H

#include <map>

/** \brief This class holds ready to use values for some parameters the
 * built in surfaces.
 *
 *
 * This class is primarily there to store some basic properties of the
 * basic surfaces built in the system. So far it is a lookup table
 * surface_level_set <--> channel volume proportions
 *
 */
class SURFACE_AREAS_TABLES {

public:
  SURFACE_AREAS_TABLES(){}

  /** 
   * \brief returns the level set value for a given channel proportion
   * 
   * Computes the level set value for a given surface and channel
   * proportion. Interpolates linear bewtween two samples and returns
   * the minimum/maximum values if the input is outside of the sampled
   * area
   *
   * \param[in] surface The surface type
   * \param[in] val The channel proportion value
   * \param[out] The level set value
   */
  const double get_level( std::string surface, double val ) const {


    const std::map<double, double>* map;


    if( surface == "Primitive" ){
      map = &P_SURFACE_VOL;
    } else if ( surface == "Diamond" ){
      map = &D_SURFACE_VOL;
    } else if ( surface == "Gyroid" ){
      map = &G_SURFACE_VOL;
    } else if ( surface.substr(0,11) == "Lonsdaleite" ){
      map = &L_SURFACE_VOL;
    } else {
      map = &GENERAL_INV;
      std::cerr << "surface_tables.hpp (line: "
		<< std::to_string(__LINE__)
		<< "): Unkown surface type" << std::endl;
    }


    auto it = map->upper_bound( val );
    

    if( it == map->begin() ){
      return it->second;
    } else if ( it == map->end() ){
      return map->rbegin()->second;
    } else {

      double hi_v = it->second, hi_p = it->first;
      it--;
      double lo_v = it->second, lo_p = it->first;

      double fraction = (val - lo_p) / (hi_p - lo_p);
      return lo_v + fraction * (hi_v - lo_v);
    }
   
  }




  /** 
   * \brief returns the volume proportion value for a given level set
   * 
   * Computes the volume proportion for a given surface and level
   * set. Interpolates linear bewtween two samples and returns the
   * minimum/maximum values if the input is outside of the sampled
   * area
   *
   * \param[in] surface The surface type
   * \param[in] val the level set value
   * \param[out] The channel proportion value
   */
  const double get_prop( std::string surface, double val ) const {

    const std::map<double, double>* map;

    if( surface == "Primitive" ){
      map = &P_SURFACE_VOL_INV;
    } else if ( surface == "Diamond" ){
      map = &D_SURFACE_VOL_INV;
    } else if ( surface == "Gyroid" ){
      map = &G_SURFACE_VOL_INV;
    } else if ( surface.substr(0, 11)  == "Lonsdaleite" ){
      map = &L_SURFACE_VOL_INV;      
    } else {
      map = &GENERAL_INV;
      std::cerr << "surface_tables.hpp (line: "
		<< std::to_string(__LINE__)
		<< "): Unkown surface type" << std::endl;
    }
   
    
    auto it = map->upper_bound( val );

    if( it == map->begin() ){
      return it->second;
    } else if ( it == map->end() ){
      return map->rbegin()->second;
    } else {
      double hi_v = it->second, hi_p = it->first;
      it--;
      double lo_v = it->second, lo_p = it->first;

      double fraction = (val - lo_p) / (hi_p - lo_p);
      return lo_v + fraction * (hi_v - lo_v);
    }   
  }  
  
private:

  


  const std::map<double, double> GENERAL = {{-1, -1}, {0, 0}, {1, 1}};
  const std::map<double, double> GENERAL_INV = {{-1, -1}, {0, 0}, {1, 1}};  
  

  const std::map<double, double> D_SURFACE_VOL = {{1, -3.5}, {1, -3.4}, {1, -3.3}, {1, -3.2}, {1, -3.1}, {1, -3}, {1, -2.9}, {1, -2.8}, {1, -2.7}, {1, -2.6}, {1, -2.5}, {1, -2.4}, {1, -2.3}, {1, -2.2}, {1, -2.1}, {1, -2}, {1, -1.9}, {1, -1.8}, {1, -1.7}, {1, -1.6}, {1, -1.5}, {1, -1.4}, {1, -1.3}, {1, -1.2}, {1, -1.1}, {0.998869, -1}, {0.985201, -0.9}, {0.957571, -0.8}, {0.91154, -0.7}, {0.851888, -0.6}, {0.791613, -0.5}, {0.732402, -0.4}, {0.673542, -0.3}, {0.615889, -0.2}, {0.557682, -0.1}, {0.5, 1.97065e-15}, {0.442318, 0.1}, {0.384112, 0.2}, {0.326458, 0.3}, {0.267598, 0.4}, {0.208387, 0.5}, {0.148112, 0.6}, {0.08846, 0.7}, {0.0424287, 0.8}, {0.0147988, 0.9}, {0.001131, 1}, {0, 1.1}, {0, 1.2}, {0, 1.3}, {0, 1.4}, {0, 1.5}, {0, 1.6}, {0, 1.7}, {0, 1.8}, {0, 1.9}, {0, 2}, {0, 2.1}, {0, 2.2}, {0, 2.3}, {0, 2.4}, {0, 2.5}, {0, 2.6}, {0, 2.7}, {0, 2.8}, {0, 2.9}, {0, 3}, {0, 3.1}, {0, 3.2}, {0, 3.3}, {0, 3.4}, {0, 3.5}};


  const std::map<double, double> D_SURFACE_VOL_INV = {{-3.5, 1}, {-3.4, 1}, {-3.3, 1}, {-3.2, 1}, {-3.1, 1}, {-3, 1}, {-2.9, 1}, {-2.8, 1}, {-2.7, 1}, {-2.6, 1}, {-2.5, 1}, {-2.4, 1}, {-2.3, 1}, {-2.2, 1}, {-2.1, 1}, {-2, 1}, {-1.9, 1}, {-1.8, 1}, {-1.7, 1}, {-1.6, 1}, {-1.5, 1}, {-1.4, 1}, {-1.3, 1}, {-1.2, 1}, {-1.1, 1}, {-1, 0.998869}, {-0.9, 0.985201}, {-0.8, 0.957571}, {-0.7, 0.91154}, {-0.6, 0.851888}, {-0.5, 0.791613}, {-0.4, 0.732402}, {-0.3, 0.673542}, {-0.2, 0.615889}, {-0.1, 0.557682}, {1.97065e-15, 0.5}, {0.1, 0.442318}, {0.2, 0.384112}, {0.3, 0.326458}, {0.4, 0.267598}, {0.5, 0.208387}, {0.6, 0.148112}, {0.7, 0.08846}, {0.8, 0.0424287}, {0.9, 0.0147988}, {1, 0.001131}, {1.1, 0}, {1.2, 0}, {1.3, 0}, {1.4, 0}, {1.5, 0}, {1.6, 0}, {1.7, 0}, {1.8, 0}, {1.9, 0}, {2, 0}, {2.1, 0}, {2.2, 0}, {2.3, 0}, {2.4, 0}, {2.5, 0}, {2.6, 0}, {2.7, 0}, {2.8, 0}, {2.9, 0}, {3, 0}, {3.1, 0}, {3.2, 0}, {3.3, 0}, {3.4, 0}, {3.5, 0}};


const std::map<double, double> P_SURFACE_VOL = {{1, -3.5}, {1, -3.4}, {1, -3.3}, {1, -3.2}, {1, -3.1}, {0.999868, -3}, {0.998268, -2.9}, {0.995374, -2.8}, {0.991538, -2.7}, {0.986879, -2.6}, {0.981446, -2.5}, {0.975221, -2.4}, {0.968308, -2.3}, {0.960604, -2.2}, {0.952153, -2.1}, {0.942825, -2}, {0.932753, -1.9}, {0.921687, -1.8}, {0.909752, -1.7}, {0.896707, -1.6}, {0.88256, -1.5}, {0.867096, -1.4}, {0.850201, -1.3}, {0.831506, -1.2}, {0.810543, -1.1}, {0.785793, -1}, {0.757404, -0.9}, {0.728574, -0.8}, {0.699913, -0.7}, {0.671256, -0.6}, {0.642623, -0.5}, {0.614062, -0.4}, {0.585641, -0.3}, {0.55707, -0.2}, {0.5284, -0.1}, {0.5, 1.97065e-15}, {0.4716, 0.1}, {0.442931, 0.2}, {0.414359, 0.3}, {0.385939, 0.4}, {0.357377, 0.5}, {0.328744, 0.6}, {0.300087, 0.7}, {0.271426, 0.8}, {0.242596, 0.9}, {0.214207, 1}, {0.189457, 1.1}, {0.168494, 1.2}, {0.149799, 1.3}, {0.132904, 1.4}, {0.11744, 1.5}, {0.103293, 1.6}, {0.0902477, 1.7}, {0.0783134, 1.8}, {0.0672467, 1.9}, {0.0571755, 2}, {0.047847, 2.1}, {0.0393964, 2.2}, {0.0316923, 2.3}, {0.0247788, 2.4}, {0.0185541, 2.5}, {0.0131214, 2.6}, {0.00846163, 2.7}, {0.00462563, 2.8}, {0.001732, 2.9}, {0.000131812, 3}, {0, 3.1}, {0, 3.2}, {0, 3.3}, {0, 3.4}, {0, 3.5}};


const std::map<double, double> P_SURFACE_VOL_INV = {{-3.5, 1}, {-3.4, 1}, {-3.3, 1}, {-3.2, 1}, {-3.1, 1}, {-3, 0.999868}, {-2.9, 0.998268}, {-2.8, 0.995374}, {-2.7, 0.991538}, {-2.6, 0.986879}, {-2.5, 0.981446}, {-2.4, 0.975221}, {-2.3, 0.968308}, {-2.2, 0.960604}, {-2.1, 0.952153}, {-2, 0.942825}, {-1.9, 0.932753}, {-1.8, 0.921687}, {-1.7, 0.909752}, {-1.6, 0.896707}, {-1.5, 0.88256}, {-1.4, 0.867096}, {-1.3, 0.850201}, {-1.2, 0.831506}, {-1.1, 0.810543}, {-1, 0.785793}, {-0.9, 0.757404}, {-0.8, 0.728574}, {-0.7, 0.699913}, {-0.6, 0.671256}, {-0.5, 0.642623}, {-0.4, 0.614062}, {-0.3, 0.585641}, {-0.2, 0.55707}, {-0.1, 0.5284}, {1.97065e-15, 0.5}, {0.1, 0.4716}, {0.2, 0.442931}, {0.3, 0.414359}, {0.4, 0.385939}, {0.5, 0.357377}, {0.6, 0.328744}, {0.7, 0.300087}, {0.8, 0.271426}, {0.9, 0.242596}, {1, 0.214207}, {1.1, 0.189457}, {1.2, 0.168494}, {1.3, 0.149799}, {1.4, 0.132904}, {1.5, 0.11744}, {1.6, 0.103293}, {1.7, 0.0902477}, {1.8, 0.0783134}, {1.9, 0.0672467}, {2, 0.0571755}, {2.1, 0.047847}, {2.2, 0.0393964}, {2.3, 0.0316923}, {2.4, 0.0247788}, {2.5, 0.0185541}, {2.6, 0.0131214}, {2.7, 0.00846163}, {2.8, 0.00462563}, {2.9, 0.001732}, {3, 0.000131812}, {3.1, 0}, {3.2, 0}, {3.3, 0}, {3.4, 0}, {3.5, 0}};



const std::map<double, double> G_SURFACE_VOL = {{1, -3.5}, {1, -3.4}, {1, -3.3}, {1, -3.2}, {1, -3.1}, {1, -3}, {1, -2.9}, {1, -2.8}, {1, -2.7}, {1, -2.6}, {1, -2.5}, {1, -2.4}, {1, -2.3}, {1, -2.2}, {1, -2.1}, {1, -2}, {1, -1.9}, {1, -1.8}, {1, -1.7}, {1, -1.6}, {0.998461, -1.5}, {0.974875, -1.4}, {0.937104, -1.3}, {0.900082, -1.2}, {0.864131, -1.1}, {0.828969, -1}, {0.794599, -0.9}, {0.760732, -0.8}, {0.727404, -0.7}, {0.694383, -0.6}, {0.661546, -0.5}, {0.628882, -0.4}, {0.596636, -0.3}, {0.564302, -0.2}, {0.532064, -0.1}, {0.5, 1.97065e-15}, {0.467936, 0.1}, {0.435698, 0.2}, {0.403364, 0.3}, {0.371119, 0.4}, {0.338453, 0.5}, {0.305617, 0.6}, {0.272596, 0.7}, {0.239268, 0.8}, {0.205401, 0.9}, {0.171031, 1}, {0.135869, 1.1}, {0.0999181, 1.2}, {0.0628956, 1.3}, {0.0251254, 1.4}, {0.0015395, 1.5}, {0, 1.6}, {0, 1.7}, {0, 1.8}, {0, 1.9}, {0, 2}, {0, 2.1}, {0, 2.2}, {0, 2.3}, {0, 2.4}, {0, 2.5}, {0, 2.6}, {0, 2.7}, {0, 2.8}, {0, 2.9}, {0, 3}, {0, 3.1}, {0, 3.2}, {0, 3.3}, {0, 3.4}, {0, 3.5}};


  const std::map<double, double> G_SURFACE_VOL_INV = {{-3.5, 1}, {-3.4, 1}, {-3.3, 1}, {-3.2, 1}, {-3.1, 1}, {-3, 1}, {-2.9, 1}, {-2.8, 1}, {-2.7, 1}, {-2.6, 1}, {-2.5, 1}, {-2.4, 1}, {-2.3, 1}, {-2.2, 1}, {-2.1, 1}, {-2, 1}, {-1.9, 1}, {-1.8, 1}, {-1.7, 1}, {-1.6, 1}, {-1.5, 0.998461}, {-1.4, 0.974875}, {-1.3, 0.937104}, {-1.2, 0.900082}, {-1.1, 0.864131}, {-1, 0.828969}, {-0.9, 0.794599}, {-0.8, 0.760732}, {-0.7, 0.727404}, {-0.6, 0.694383}, {-0.5, 0.661546}, {-0.4, 0.628882}, {-0.3, 0.596636}, {-0.2, 0.564302}, {-0.1, 0.532064}, {1.97065e-15, 0.5}, {0.1, 0.467936}, {0.2, 0.435698}, {0.3, 0.403364}, {0.4, 0.371119}, {0.5, 0.338453}, {0.6, 0.305617}, {0.7, 0.272596}, {0.8, 0.239268}, {0.9, 0.205401}, {1, 0.171031}, {1.1, 0.135869}, {1.2, 0.0999181}, {1.3, 0.0628956}, {1.4, 0.0251254}, {1.5, 0.0015395}, {1.6, 0}, {1.7, 0}, {1.8, 0}, {1.9, 0}, {2, 0}, {2.1, 0}, {2.2, 0}, {2.3, 0}, {2.4, 0}, {2.5, 0}, {2.6, 0}, {2.7, 0}, {2.8, 0}, {2.9, 0}, {3, 0}, {3.1, 0}, {3.2, 0}, {3.3, 0}, {3.4, 0}, {3.5, 0}};



  const std::map<double, double> L_SURFACE_VOL = {{1, -3.5}, {1, -3.4}, {1, -3.3}, {1, -3.2}, {1, -3.1}, {1, -3}, {1, -2.9}, {1, -2.8}, {1, -2.7}, {1, -2.6}, {1, -2.5}, {1, -2.4}, {1, -2.3}, {1, -2.2}, {1, -2.1}, {1, -2}, {1, -1.9}, {1, -1.8}, {1, -1.7}, {1, -1.6}, {1, -1.5}, {1, -1.4}, {1, -1.3}, {1, -1.2}, {1, -1.1}, {1, -1}, {1, -0.9}, {0.999936, -0.8}, {0.994677, -0.7}, {0.941418, -0.6}, {0.786456, -0.5}, {0.662825, -0.4}, {0.613586, -0.3}, {0.578215, -0.2}, {0.546349, -0.1}, {0.515848, 1.97065e-15}, {0.485498, 0.1}, {0.454102, 0.2}, {0.415819, 0.3}, {0.350657, 0.4}, {0.224036, 0.5}, {0.0670281, 0.6}, {0.009234, 0.7}, {0, 0.8}, {0, 0.9}, {0, 1}, {0, 1.1}, {0, 1.2}, {0, 1.3}, {0, 1.4}, {0, 1.5}, {0, 1.6}, {0, 1.7}, {0, 1.8}, {0, 1.9}, {0, 2}, {0, 2.1}, {0, 2.2}, {0, 2.3}, {0, 2.4}, {0, 2.5}, {0, 2.6}, {0, 2.7}, {0, 2.8}, {0, 2.9}, {0, 3}, {0, 3.1}, {0, 3.2}, {0, 3.3}, {0, 3.4}, {0, 3.5}};


  const std::map<double, double> L_SURFACE_VOL_INV = {{-3.5, 1}, {-3.4, 1}, {-3.3, 1}, {-3.2, 1}, {-3.1, 1}, {-3, 1}, {-2.9, 1}, {-2.8, 1}, {-2.7, 1}, {-2.6, 1}, {-2.5, 1}, {-2.4, 1}, {-2.3, 1}, {-2.2, 1}, {-2.1, 1}, {-2, 1}, {-1.9, 1}, {-1.8, 1}, {-1.7, 1}, {-1.6, 1}, {-1.5, 1}, {-1.4, 1}, {-1.3, 1}, {-1.2, 1}, {-1.1, 1}, {-1, 1}, {-0.9, 1}, {-0.8, 0.999936}, {-0.7, 0.994677}, {-0.6, 0.941418}, {-0.5, 0.786456}, {-0.4, 0.662825}, {-0.3, 0.613586}, {-0.2, 0.578215}, {-0.1, 0.546349}, {1.97065e-15, 0.515848}, {0.1, 0.485498}, {0.2, 0.454102}, {0.3, 0.415819}, {0.4, 0.350657}, {0.5, 0.224036}, {0.6, 0.0670281}, {0.7, 0.009234}, {0.8, 0}, {0.9, 0}, {1, 0}, {1.1, 0}, {1.2, 0}, {1.3, 0}, {1.4, 0}, {1.5, 0}, {1.6, 0}, {1.7, 0}, {1.8, 0}, {1.9, 0}, {2, 0}, {2.1, 0}, {2.2, 0}, {2.3, 0}, {2.4, 0}, {2.5, 0}, {2.6, 0}, {2.7, 0}, {2.8, 0}, {2.9, 0}, {3, 0}, {3.1, 0}, {3.2, 0}, {3.3, 0}, {3.4, 0}, {3.5, 0}};

  

};





#endif
