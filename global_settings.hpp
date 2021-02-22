/* Projection tool - compute planar projections of triply periodic
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


/*
 * this file just holds some globally used constants
 */

#ifndef GLOBAL_SETTINGS_H
#define GLOBAL_SETTINGS_H

#include <string>

const std::string g_color_b1 = "#5aa02c";
const std::string g_color_b2 = "#ffcc00";
const std::string g_color_n  = "#c83737";

//
// the resolution used to measure the volume of the channels.
// -1 means the current resolution from the GUI is used
//
const int res_measure_vol = 76;


//
// the resolution used to measure the volume of the channels.
// -1 means the current resolution from the GUI is used
//
const int res_measure_area = 76;

//
// the resolution used to measure the percolation thresdhold
// -1 means the current resolution from the GUI is used
//
const int res_measure_perc = 77;


/* 
 * add some code that tries to read these parameters from a config
 * file, if it does not succeed, just let them be as is
 */ 

#endif
