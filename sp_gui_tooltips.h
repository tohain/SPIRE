
#include <string>


class tooltips {
 public:
  tooltips(){}
  ~tooltips(){}
  
  //const std::string ntucs_tooltip = "The number of translational unit cells to put in the projection box.";

  const std::string aa_tooltip = "The unitcell size in x and y direction in absolute length units.";
  const std::string ac_tooltip = "The unitcell size in z direction in absolute length units (only active for structures with hexagonal symmetry)";

  const std::string level_set_tooltip = "Surface control parameter; interpreted as chosen next to it";
  const std::string level_par_tooltip = "Chooses if the current surface control parameter is interpreted as a volume fill fraction, or the mathematical level set value used in the nodal representation.";
  
  const std::string type_tooltip = "The type of surface to project.";

  const std::string theta_tooltip = "The theta angle of the normal vector denoting the orientation of the slice.";
  const std::string phi_tooltip = "The phi angle of the normal vector denoting the orientation of the slice.";
  
  const std::string hkl_tooltip = "The Miller indeces indicating the orientation of the projection box in the periodic surface structure. These can be seen as the components of the normal vector of the xy-plane of the projection box";

  const std::string set_to_uc_dim_tooltip = "Sets the dimensions of the slice so it is periodic in all directions. Leaves values unchanged if no periodicity is found in a direction";
  
  const std::string slicewidth_tooltip = "The width (thickness of slice) of the projection box in absolute length units";
  const std::string sliceheight_tooltip = "The height (y) of the projection box in absolute length units";
  const std::string slicelength_tooltip = "The length (x)  of the projection box in absolute length units";
  const std::string sliceposition_tooltip = "The position of the center of the projection along the normal vector of the slice in absolute length units";
  

  
  const std::string nr_points_xy_tooltip = "The number of samples to put in the projection box in x direction. Corresponds to the final resolution of the projection. IMPORTANT: large (>500) values significantly increase memory and ressource usage and may crash your maschine";
  const std::string nr_points_z_tooltip = "Number of samples to put in the projection box in z direction. Increasing this value may help if artifacts and \"steps\" in the image occur. IMPORTANT: large (>250) values significantly increase memory and ressource usage and may crash your maschine";
   
  const std::string invert_tooltip = "Inverts the colors in the picture. If inverted the pixel with the highest value is black.";
  const std::string autoupdate_tooltip = "Automatically update thr projection as soon as a parameter changes, can result in significant delays and lags in the usage of the user interface!";
  const std::string image_scaling_tooltip = "Chose between a linear or logarithmic scale between pixel brightness and value";  

  const std::string membranes_settings_tooltip = "Tune the width and distance of additional membranes.";
  const std::string membranes_add_tooltip = "Add an additional membrane.";
  const std::string membranes_remove_tooltip = "Remove the currently selected membrane";

  const std::string channel_fill_tooltip = "Controls wether a channel or membrane is filled solid or not";
  
  const std::string button_save = "Save the currentyl displayed image to a PNG file";
  const std::string button_quit = "Closes this window and quits the application";
  const std::string button_measure_va = "Measures volume and areas of channels and membranes. This may take up to several minutes. Wait for the status indicator to turn green.";
  const std::string button_measure_percthres = "Computes the percolation threshold of all channels. That is the diameter of the largest sphere which parely can travel unrestricted through the structure. This may take up to several minutes. Wait for the status indicator to turn green.";  
  const std::string button_render = "Computes and displays the projection for the current parameters. This may take up to several minutes. Wait for the status indicator to turn green";

  const std::string button_write_pars = "Writes the current parameters in an ASCI file";
  const std::string button_read_pars = "Reads parameters from an ASCI file";

  const std::string button_save_grid = "Saves all xyz positions, labels and color of all points in the slice to an ASCI files. This may take up several minutes and hundreds of MB disk space, depending on the resolution.";
  const std::string button_save_membranes = "Saves the xyz positions of all the points in a membranes to an ASCI file. Each membrane is saved in to a seperate file";
  const std::string button_save_network = "Saves the xyz positions of the minimal topological network of the two main channels (inner and outermost channel) into an ASCI file.";

  const std::string save_prefix_tooltip = "Enter or choose a path and a prefix which is prepended to each of the files saved below. If empty the current directory is chosen";
  
  const std::string update_prev_tooltip = "Controls up to how many points the live preview should be computed. Over a certain number of points (depends on your machine) it will be too slow to compute live";



  const std::string choose_prefix_path = "Open dialog to choose the location and a filename prefix for the files to save";

  
  const std::string manual = "How it works: at first a set of points is generated, these points are located within a rectangular box, which will be called *slice* herein. The size and orientation as well as the number of points allocated in this slice essentially controls the projected image you will compute in the end: the size and orientation determines the region of the 3D structure to be projected, whereas the number of points controls the resolution and hence the quality of the projection: every column of points will collaps into a single pixel in the final projection. The available parameters to tune this slice can be found further down in this document.\n\n"
    "To compute the projection, each point in this slice begins as non-colored, or non-marked. In this picture, these points appear red. The slice with the points is then positioned in the 3D structure, as can be seen in the image (ref here!). Then every point which intersects the 3D structure (here a gyroid membrane) is assigned a different color, here black. To compute the final projection, the number of black points is summed up for each column in the slice. This sum then represents the brightness of the pixel in the 2D projection.\n\n"
    "So far, the 3D structure must be given as an isosurface, i.e. as a implicit formula: f(x,y,z)=c. A point is then marked, if c - T/2 < f(px,py,pz) < c + T/2, where T is the thickness of the membrane. Note, that the thickness T here is given as a unitless parameter, which is difficult to grasp in physical terms. An easier parameter would be a thickness in length units. This can be achieved using a distance map. A (euclidean) distance map assigns each points in the slice a scalar value, in this case it is the distance of each point to the nearest marked (black) points. For all marked points this is of course 0. Since these values are in length units, we now have a tool to control the thickness of the membrane in length units.\n\n"
    "To begin this, however, a first membrane is needed, herein calle \"main membrane\", or \"innermost membrane\", from which the distance map can be computed, we empirically determined a parameter of T=0.03, which yielded a membrane with a minimal thickness, however, without holes either. From this, the thickness can be controlled as mentioned above. Further membranes at distance d with width w (both in length units) can be realised by just marking all pixels, which values in the distance map fulfil d - w/2 < distance_map_value < d + w/2.\n\n"
    "Currently, only TPMS are implemented in the software, in this case, the parameter c can be interpreted as a more natural parameter: since the membrane separates the space into two channels, the position of the membrane, and hence c, controls the proportions of the volumes of these two channels. We use a lookup table, so the user can tune the channel volume proportion, but internally the correct value for the parameter c is computed.";
  
};

