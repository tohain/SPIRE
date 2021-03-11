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



#ifndef SLICE_ORIENTATION_VISUALISATION_H
#define  SLICE_ORIENTATION_VISUALISATION_H

#include <iostream>
#include <vector>
#include <string>

#include "drawme.hpp"
#include "global_settings.hpp"

/// creates an svg in a single string visualizing the slice orientation
/// in the simulation box
std::string draw_slice_visualization( std::vector<double> base,
				      std::vector<double> L,
				      global_settings &gs  );

#endif
