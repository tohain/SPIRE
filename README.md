# surface_projection
A tool to compute projections of 3D surfaces


## Install

### Linux/Unix (and also Mac)

This software uses CMake as build system. To build from the sources, following packages are needed:

- CMake (tested with 3.16.3)
- C++ compiler with C++ 2011 support (tested with g++ vers. 8+)
- Qt min vers. 5.12.5 ( and its dependencies, including pixel image libraries, e.g. libpng )
- CGAL (tested with 4.14) (and its dependencies, including boost, gmp and mpfr)

To build, navigate in the source directory (containing qt_main.cpp) and execute the commands:

mkdir build && cd build
cmake -CMAKE_BUILD_TYPE=Release ../
make

If successfull, the compiled binary ( projection_tool_qt ) should be in the directory <sources>/build/projection_tool_qt.

There are some options to be used with the cmake command to be considered:

-BUILD_QT_GUI=ON(OFF): uses QT as GUI library
-BUILD_WX_GUI=ON(OFF): uses wxWidgets as GUI library (deprecated and will be removed in future versions)
-USE_CGAL=ON(OFF): use CGAL libraries

If no GUI option is chosen, a minimial CLI interface will be used. However, for now, it is just a dummy code and does nothing.

Although the code will build and have most of its functionality without CGAL, some measurements (surface area of membranes) can **NOT** be computed.


### Windows

A pre-compiled binary for Windows users including all necessary libraries is provided. Just download the archive, extract it and run surface_projection.exe.

To compile from source, use your favourite build system. Tested with mingw32-64.



## Usage

Though we do hope most of the tool is self-explanatory, here are some instructions:

The main window is split in two: the right panel shows the projected image of the 3D surface, the left panel controls the parameter of the surface and the computations. The bottom of the window shows a status bar. The right half of the status bar is used to display messages to the user, the left part shows some measurements of the structures, as well as the current status of the tool.

### General concepts

The basic functionality of the software is as followed:

At first a set of points is generated, these points are located within a rectangular region, which will be called *slice* herein. The size and orientation as well as the number of points allocated in this slice essentially controls the projected image you will compute in the end: the size and orientation determines the region of the 3D structure to be projected, whereas the number of points controls the resolution and hence the quality of the projection: every column of points will collaps into a single pixel in the final projection. The available parameters to tune this slice can be found further down in this document.

To compute the projection, each point in this slice begins as non-colored, or non-marked. In this picture, these points appear red. The slice with the points is then positioned in the 3D structure, as can be seen in the image (ref here!). Then every point which intersects the 3D structure (here a gyroid membrane) is assigned a different color, here black. To compute the final projection, the number of black points is summed up for each column in the slice. This sum then represents the brightness of the pixel in the 2D projection.

So far, the 3D structure must be given as an isosurface, i.e. as a implicit formula: f(x,y,z)=c. A point is then marked, if c - T/2 < f(px,py,pz) < c + T/2, where T is the thickness of the membrane. Note, that the thickness T here is given as a unitless parameter, which is difficult to grasp in physical terms. An easier parameter would be a thickness in length units. This can be achieved using a distance map. A (euclidean) distance map assigns each points in the slice a scalar value, in this case it is the distance of each point to the nearest marked (black) points. For all marked points this is of course 0. Since these values are in length units, we now have a tool to control the thickness of the membrane in length units.

To begin this, however, a first membrane is needed, herein calle "main membrane", or "innermost membrane", from which the distance map can be computed, we empirically determined a parameter of T=0.03, which yielded a membrane with a minimal thickness, however, without holes either. From this, the thickness can be controlled as mentioned above. Further membranes at distance d with width w (both in length units) can be realised by just marking all pixels, which values in the distance map fulfil d - w/2 < distance_map_value < d + w/2.

Currently, only TPMS are implemented in the software, in this case, the parameter c can be interpreted as a more natural parameter: since the membrane separates the space into two channels, the position of the membrane, and hence c, controls the proportions of the volumes of these two channels. We use a lookup table, so the user can tune the channel volume proportion, but internally the correct value for the parameter c is computed.


### GUI

#### Status Bar

##### Status indicators

The right part of the status bar shows some measurements as well as the status. Since the computations of the projection and the measurements can take up to several minutes, the status indicator is implemented to show, if the software is still working. A read box and the label "Busy" means, the software is still working and the view in the projection image panel is **NOT** up to date. The same goes for the Measurement box: a red box indicates, the software is still working and the shown  are not up to date.

##### Measurements

The "box size XY" value shows the size of the entire box (all of what is seen in the projection image panel) in the same unit as the unit cell size (see extra section). I.e. if unit cell size is given in nm, so is the box size.

The same goes for the box height, it shows the height of the box in same units as the unit cell size.

The left side of the status bar shows important messages, for example the status of the software or if some choices of parameters are invalid.

#### Parameter control

Here we want to give a short overview of the parameters and their usage

##### Unit cells

This parameter controls the size of the slice compared to the 3D structure. It represents how many unitcell length the x and y dimensions of the slice will have.

##### Unit cell size

This parameter sets the lengthscale of the entire 3D structure and the slice. It does not really have any influence on the projection per se, but scales all values in units of lengths.

##### Volume proportions

See above for detailed explanations. Sets the proportion of the two channels separated by the main mebrane (i.e. the one directly evaluated with the isosurface formula, see above). Note that further membranes do not have any influence on this.

##### Surface type

Choses the 3D structre to be projected

##### Slice width

This is the width of the slize, meaning its length in z direction. This value is always a fraction of the current **periodicity length**. The periodicity length is the distance on must go in the currently set up direction (see hkl), to reach a mathematically identical point (i.e. advance one unit cell in the current orientation). For a direction of (001) and a unit cell size of 1nm, this would be 1nm, for a direction of (111) this would be sqrt(3)*1nm.

The **real** width of the slice in length units can be found in the status bar, in the bottom right corner of the main window.

##### Slice height

The slice height controls the location of the slice in z direction. It is also a fraction of the current **periodicity length**. For a slice width value of 1, this parameter has no influence, since always an entire unitcell is covered.

##### Membranes

The position and distance of all membranes can directly be entered in to the table. The two buttons add or remove the currently selected membrane. Please note that due to technical reasons, the innermost membrane **can not** be deleted, **must** be at a distance of 0, and **must have** a minimal width of 0.02. Further details are given above.

##### Resolution

These values control the number of points in the slice and the pixels of the projection. While XY resolution directly translates into the resolution of the projected pixel image, the Z resolution reduced artifacts in the projection. The scaling gives you different options for scaling, e.g. for a highly dynamical image with strongly varying pixel values, a logarithmic scaling reveals more details.

Please note that these values also determine the performance of this tool. Smaller resolutions of around 200x200x150 pixels only take a couple seconds, whereas large projections comprising 750*750*500 pixels can take a really long time. **More important though: the amount of memory needed increases also, there is no check to avoid large runs, so this could crash your maschine!!!**

The invert colors box should be self-explanatory

##### Autoupdate

This box enables autoupdate. This means, as soon as a parameter is changed, the projection will be updated. Depending on the resolution this may take up to several minutes, so ticking this box is really only recommended for small resolutions.