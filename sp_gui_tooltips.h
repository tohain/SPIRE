#include <string>


class tooltips {
 public:
  tooltips(){}

  const std::string ntucs_tooltip = "This controls the edge length of the quadratic slice. The length is the number of unit cells times the unit cell length in {100},{010},{001}.";

  const std::string ntucs_p_tooltip = "This controls the edge length of the quadratic slice. The length is the number of unit cells times the periodicity length in the current orientatin.";

  const std::string a_tooltip = "The unit cell length in {100},{010},{001} orientation in arbitrary length units";

  const std::string type_tooltip = "The type of surface to project.";

  const std::string memwidth_tooltip = "The width of the membrane, therefore the distance between the two surfaces enclosing the memebrane";

  const std::string theta_tooltip = "The theta angle of the normal vector denoting the orientation of the slice.";
  const std::string phi_tooltip = "The phi angle of the normal vector denoting the orientation of the slice.";
  
  const std::string hkl_tooltip = "Miller indeces describing the orientation of the slice";

  const std::string slicewidth_tooltip = "The width of the slice in length units.";
  
  const std::string slicewidth_p_tooltip = "The width of the slice as a fraction of the periodicitiy length in the current orientation. 1 means the slice width covers an entire unit cell.";
  
  const std::string sliceheight_tooltip = "The position of the center of the slice along the normal vector of the slice in absolute length units. 0 means the slice cuts the origin.";
  
  const std::string sliceheight_p_tooltip = "The position of the center of the slice along the normal vector of the slice as a fraction of the periodicity length. 0 and 1 are therefore equivalent.";
  
  const std::string nr_points_xy_tooltip = "This parameter controls the resolution and is the number of samples along a single edge of the slice in the x and y direction (i.e. perpendicular to the normal vector). High number of points means better resolution but slower computation. The distance between the single samples (should equal picel size) is given below (dx,dy).";
  
  const std::string nr_points_z_tooltip = "Number of points (resolution) in the direction of the normal vector. Distance between the samples is given below";
  
  
  const std::string invert_tooltip = "Inverts the colors in the picture. If inverted the brightest pixel is black.";

  const std::string update_prev_tooltip = "Controls up to how many points the live preview should be computed. Over a certain number of points (depends on your machine) it will be too slow to compute live";

  const std::string filename_tooltip = "The path and filename to save the currently displayed image. Use the browse button to open a file dialog. The extensions determines the format the image is saved to. Currently supported: *.png and *.tiff";

};

