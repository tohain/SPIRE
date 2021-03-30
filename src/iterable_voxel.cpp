/* Projection tool - compute planar projection of triply periodic
 * minimal surfaces 
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



#include "iterable_voxel.hpp"

/// +/-x = f/b, +/-y = l/r,
/// +/-z=u/d
std::vector< std::pair<int, std::string> > iterable_voxel::get_6_neighbors_dir(){
  
  std::vector< std::pair<int, std::string> > nbs (6, std::pair<int, std::string> (0, "") );
  
  nbs[0].first = u()();
  nbs[0].second = "001";
  nbs[1].first = d()();
  nbs[1].second = "00-1";
  
  nbs[2].first = l()();
  nbs[2].second = "010";
  nbs[3].first = r()();
  nbs[3].second = "0-10";
  
  nbs[4].first = f()();
  nbs[4].second = "100";
  nbs[5].first = b()();
  nbs[5].second = "-100";
  
  return nbs;
}

/// +/-x = f/b, +/-y = l/r,
/// +/-z=u/d
std::vector< std::pair<int, std::string> > iterable_voxel::get_18_neighbors_dir(){
  
  std::vector< std::pair<int, std::string> > nbs (18, std::pair<int, std::string> (0, "") );
  
  // get the six neighbors
  auto six_nbs = get_6_neighbors_dir();
  std::copy( six_nbs.begin(), six_nbs.end(), nbs.begin() );
  
  nbs[6].first = f().l()();
  nbs[6].second = "110";
  nbs[7].first = f().r()();
  nbs[7].second = "1-10";
  nbs[8].first = b().l()();
  nbs[8].second = "-110";
  nbs[9].first = b().r()();
  nbs[9].second = "-1-10";
  
  nbs[10].first = l().u()();
  nbs[10].second = "011";
  nbs[11].first = l().d()();
  nbs[11].second = "01-1";
  
  nbs[12].first = r().u()();
  nbs[12].second = "0-11";
  nbs[13].first = r().d()();
  nbs[13].second = "0-1-1";
  
  nbs[14].first = f().u()();
  nbs[14].second = "101";
  nbs[15].first = f().d()();
  nbs[15].second = "10-1";
  
  nbs[16].first = b().u()();
  nbs[16].second = "-101";      
  nbs[17].first = b().d()();
  nbs[17].second = "-10-1";
  
  return nbs;
}    

/// Returns an array containin all adjacent voxels (including
/// diagonals )

/// +/-x = f/b, +/-y = l/r,
/// +/-z=u/d
std::vector< std::pair<int, std::string> > iterable_voxel::get_26_neighbors_dir(){
  
  std::vector< std::pair<int, std::string> > nbs (26, std::pair<int, std::string> (0, "") );
  
  auto tmp_nbs = get_18_neighbors_dir();
  std::copy( tmp_nbs.begin(), tmp_nbs.end(), nbs.begin() );
  
  iterable_voxel up = u();
  nbs[18].first = up.f().l()();
  nbs[18].second = "111";
  nbs[19].first = up.f().r()();
  nbs[19].second = "1-11";      
  nbs[20].first = up.b().l()();
  nbs[20].second = "-111";
  nbs[21].first = up.b().r()();
  nbs[21].second = "-1-11";
  
  iterable_voxel down = d();
  nbs[22].first  = down.f().l()();
  nbs[22].second = "11-1";
  nbs[23].first = down.f().r()();
  nbs[23].second = "1-1-1";
  nbs[24].first = down.b().l()();
  nbs[24].second = "-11-1";
  nbs[25].first = down.b().r()();
  nbs[25].second = "-1-1-1";
  
  return nbs;
}


/// retruns an array containing all voxels in a circle around the
/// current voxel dir = (000, 001, 011, 111, ..., xyz ) give the
/// vector perpendicular to the circle, where +/-x = f/b, +/-y = l/r,
/// +/-z=u/d
std::vector<int> iterable_voxel::get_circle( std::string dir ){
  
  std::vector<int> ret_voxels (8, 0);
  
  if( dir == "100" || dir == "-100" ){
    ret_voxels[0] = u()();
    ret_voxels[1] = u().r()();
    ret_voxels[2] = r()();
    ret_voxels[3] = d().r()();
    ret_voxels[4] = d()();
    ret_voxels[5] = d().l()();
	ret_voxels[6] = l()();
	ret_voxels[7] = u().l()();	
  }    
  else if( dir == "010" || dir == "0-10" ){
    ret_voxels[0] = u()();	
    ret_voxels[1] = f().u()();
    ret_voxels[2] = f()();
    ret_voxels[3] = f().d()();
    ret_voxels[4] = d()();	
    ret_voxels[5] = b().d()();
    ret_voxels[6] = b()();
    ret_voxels[7] = b().u()();
  }
  else if( dir == "001" || dir == "00-1" ){
    ret_voxels[0] = f()();
    ret_voxels[1] = f().r()();
    ret_voxels[2] = r()();
    ret_voxels[3] = b().r()();
    ret_voxels[4] = b()();
    ret_voxels[5] = b().l()();
    ret_voxels[6] = l()();
    ret_voxels[7] = f().l()();
  }
  
  
  else if( dir == "1-10" || dir == "-110" ){
    ret_voxels[0] = u()();
    ret_voxels[1] = u().r().b()();
    ret_voxels[2] = r().b()();
    ret_voxels[3] = r().b().d()();
    ret_voxels[4] = d()();	
    ret_voxels[5] = f().l().d()();
    ret_voxels[6] = f().l().u()();
    ret_voxels[7] = f().l()();	
  }
  
  else if( dir == "110" || dir == "-1-10" ){
    ret_voxels[0] = u()();	
    ret_voxels[1] = u().f().r()();
    ret_voxels[2] = f().r()();
    ret_voxels[3] = d().f().r()();
    ret_voxels[4] = d()();	
    ret_voxels[5] = d().b().l()();
    ret_voxels[6] = b().l()();
    ret_voxels[7] = u().b().l()();
  }
  
  else if( dir == "101" || dir == "-10-1" ){
    ret_voxels[0] = f().d()();
    ret_voxels[1] = f().d().r()();
    ret_voxels[2] = r()();
    ret_voxels[3] = u().r().b()();
    ret_voxels[4] = u().b()();
    ret_voxels[5] = u().b().l()();
    ret_voxels[6] = l()();
    ret_voxels[7] = l().f().d()();
  }
  
  else if( dir == "-101" || dir == "10-1" ){
    ret_voxels[0] = u().f()();
    ret_voxels[1] = u().f().r()();
    ret_voxels[2] = r()();
    ret_voxels[3] = r().b().d()();
    ret_voxels[4] = b().d()();
    ret_voxels[5] = b().d().l()();
    ret_voxels[6] = l()();
    ret_voxels[7] = l().f().u()();
  }
  
  
  else if( dir == "0-11" || dir == "01-1" ){
    ret_voxels[0] = f()();
    ret_voxels[1] = f().r().d()();
    ret_voxels[2] = r().d()();
    ret_voxels[3] = r().d().b()();
    ret_voxels[4] = b()();
    ret_voxels[5] = b().l().u()();
    ret_voxels[6] = l().u()();
    ret_voxels[7] = l().f().u()();
  }      
  
  
  else if( dir == "011" || dir == "0-1-1" ){
    ret_voxels[0] = f()();
    ret_voxels[1] = f().r().u()();
    ret_voxels[2] = r().u()();
    ret_voxels[3] = r().u().b()();
    ret_voxels[4] = b()();
    ret_voxels[5] = b().l().d()();
    ret_voxels[6] = l().d()();
    ret_voxels[7] = f().l().d()();
  }      
  
  else if( dir == "1-11" || dir == "-11-1" ){
    ret_voxels.resize( 6, 0 );
    ret_voxels[0] = r().d()();
    ret_voxels[1] = r().b()();
    ret_voxels[2] = b().u()();
    ret_voxels[3] = l().u()();
    ret_voxels[4] = f().l()();
    ret_voxels[5] = f().d()();
  }      
  
  
  else if( dir == "111" || dir == "-1-1-1" ){
    ret_voxels.resize( 6, 0 );
    ret_voxels[0] = f().d()();
    ret_voxels[1] = f().r()();
    ret_voxels[2] = u().r()();
    ret_voxels[3] = u().b()();
    ret_voxels[4] = l().b()();
    ret_voxels[5] = l().d()();
  }
  
  
  else if( dir == "-1-11" || dir == "11-1" ){
    ret_voxels.resize( 6, 0 );
    ret_voxels[0] = r().d()();
    ret_voxels[1] = b().d()();
    ret_voxels[2] = b().l()();
    ret_voxels[3] = l().u()();
    ret_voxels[4] = f().u()();
    ret_voxels[5] = f().r()();
  }
  
  else if( dir == "-111" || dir == "1-1-1" ){
    ret_voxels.resize( 6, 0 );
    ret_voxels[0] = f().u()();
    ret_voxels[1] = r().u()();
    ret_voxels[2] = d().b()();
    ret_voxels[3] = d().l()();
    ret_voxels[4] = f().l()();
    ret_voxels[5] = b().r()();
  } else {
    throw std::string ("get_circle: Direction unkown");
  } 
  
  return ret_voxels;
  
}
