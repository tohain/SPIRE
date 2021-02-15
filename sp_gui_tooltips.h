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

#include <string>

class tooltips {
 public:
  tooltips(){}
  ~tooltips(){}

  const std::string type_tooltip = "The type of surface to project.";
  const std::string aa_tooltip = "The scaling factor of the unit cell in x and y direction, given in length units.";
  const std::string ac_tooltip = "The scaling factor of the unit cell in z direction, given in length units. (only active for structures with hexagonal symmetry)";
  const std::string level_set_tooltip = "Surface control parameter; interpreted as chosen next to it";
  const std::string level_par_tooltip = "Chooses if the current surface control parameter is interpreted as a volume fill fraction, or the mathematical level set value used in the nodal representation.";
  

  const std::string slicewidth_tooltip = "The thickness of the slice (sample) in length units";
  const std::string sliceheight_tooltip = "The height (y) of the slice (sample) in absolute length units";
  const std::string slicelength_tooltip = "The width (x) of the slice (sample) in absolute length units";
  const std::string sliceposition_tooltip = "The position of the center of the slice (sample) along the normal vector in absolute length units";
  const std::string set_to_uc_dim_tooltip = "Sets the dimensions of the slice such as it fits one unit cell plus some margin";
  const std::string draw_uc_tooltip = "Draws an outline of a single unit cell";    
  const std::string hkl_tooltip = "The Miller indeces indicating the orientation of the projection box in the periodic surface structure. These can be seen as the components of the normal vector of the xy-plane of the projection box";
  
  const std::string membranes_settings_tooltip = "Tune the width and distance of additional membranes.";
  const std::string membranes_add_tooltip = "Add an additional membrane.";
  const std::string membranes_remove_tooltip = "Remove the currently selected membrane";
  const std::string channel_fill_tooltip = "Controls wether a channel or membrane is filled solid or not";


  const std::string nr_points_xy_tooltip = "The number of samples to put in the projection box in x direction. Corresponds to the final resolution of the projection. IMPORTANT: large (>500) values significantly increase memory and ressource usage and may crash your maschine";
  const std::string nr_points_z_tooltip = "Number of samples to put in the projection box in z direction. Increasing this value may help if artifacts and \"steps\" in the image occur. IMPORTANT: large (>250) values significantly increase memory and ressource usage and may crash your maschine";
   
  const std::string image_scaling_tooltip = "Chose between a linear or logarithmic scale between pixel brightness and value";
  const std::string invert_tooltip = "Inverts the colors in the picture. If inverted the pixel with the highest value is black.";
  const std::string autoupdate_tooltip = "Automatically update thr projection as soon as a parameter changes, can result in significant delays and lags in the usage of the user interface!";

  const std::string button_write_pars = "Writes the current parameters in an ASCI file";
  const std::string button_read_pars = "Reads parameters from an ASCI file";
  const std::string button_save = "Save the currently displayed image to a PNG file";

  const std::string button_measure_va = "Measures volume and areas of channels and membranes. This may take up to several minutes. Wait for the status indicator to turn green.";
  const std::string button_measure_percthres = "Computes the percolation threshold of all channels. That is the diameter of the largest sphere which barely can travel unrestricted through the structure. This may take up to several minutes. Wait for the status indicator to turn green.";  
  const std::string button_render = "Computes and displays the projection for the current parameters. This may take up to several minutes. Wait for the status indicator to turn green";  

  const std::string button_save_grid = "Saves all xyz positions, labels and color of all points in the slice to an ASCI files. This may take up several minutes and hundreds of MB disk space, depending on the resolution.";
  const std::string button_save_membranes = "Saves the xyz positions of all the points in a membranes to an ASCI file. Each membrane is saved in to a seperate file";

  const std::string save_prefix_tooltip = "Enter or choose a path and a prefix which is prepended to each of the files saved below. If empty the current directory is chosen";
  const std::string choose_prefix_path = "Open dialog to choose the location and a filename prefix for the files to save";

  const std::string status_uco_tooltip = "The size of the unit cell in the given orientation";
  const std::string status_or_tooltip = "The angles of rotation of the primitive unit cell around the y (theta) and z (phi) axis to reach the current orientation";
  const std::string status_pixs_tooltip = "Gives the nr (and size of the voxels in length units) in all spatial dimensions, dx/dy is the size of a single pixel in the projection in length units";
  const std::string status_uc_tooltip = "The size of the primitive unit cell, i.e. in (001) orientation";    
    
};

