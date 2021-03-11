# surface_projection
A tool to compute projections of 3D surfaces


## Install

### Linux/Unix (and also Mac)

This software uses CMake as build system. To build from the sources, following packages are needed:

- CMake (tested with 3.16.3)
- C++ compiler with C++ 2011 support (tested with g++ vers. 8+)
- Qt min vers. 5.12.5 ( and its dependencies ). Components required: Widgets, Gui, Svg
- libpng and libz
- libgmp
- libgmpfr
- Integer Matrix Library (included)
- any BLAS system (openblas, intel mkl, ...)

Optional:

- CGAL (tested with 4.14) (and its dependencies)

The code will build and have most of its functionality without CGAL, some measurements (surface area of membranes) can **NOT** be computed correctly.

To build, navigate in the source directory (containing qt_main.cpp) and execute the commands:

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_QT_GUI=On ../
make
```

If successfull, the compiled binary ( projection_tool_qt ) should be in the directory <sources>/build/projection_tool_qt.

There are some options to be used with the cmake command to be considered:

```-DBUILD_QT_GUI=ON(OFF)```: build the QT GUI tool

```-DBUILD_BATCH_TOOL=ON(OFF)``` build the batch tool for command line




### Windows

A pre-compiled binary for Windows users including all necessary libraries is provided. Just download the archive, extract it and run surface_projection.exe.

To compile from source, use your favourite build system. Tested with mingw32-64.



## Usage

The main window is split in two: the right panel shows the projected
image of the 3D surface as well as some information about the
structure and the projection, the left panel contains several tabs
controlling the tool: tune parameters, run measurements, bulk creation
of projections and export and saving of data and images.

### General concepts

For details about the concepts of compute projections we refer to our
publication XXX. Here only the most important concepts are provided.

#### Units

The only relevant unit in the software is that of a length. This is
not fixed to a specific value, however, all values representing a
length are provided consistently. That means if one length value is
interpreted as being in nanometers, each other length value has to be
read as nanometers as well. For example, if the unit cell scale factor
is set to 5nm, then the slice dimension parameters has to be read to
be in nanometers as well.

### Status Bar

The status bar provides important information about the currently
displayed structure
- The size of the primitive unit cell of the current structure in length units
- The orientation of the normal vector onto the slice as rotation angles in radians
- The resolution of the voxelized slice, as well as the size of a single voxel in length units. Note that the resolution of the projection is always identical to the X and Y resolution of the voxelized slice
- Unit cell size in the current orientation
- Two status indicators (projection, measurement): if red, the software is either busy computing a projection or measuring

While the software is still responsive and allows changes while
computing projections or measuring, additional computations will only
be queued. We strongly recommend not to start another compuatation
until the software finsihed its current task.

### Control Tabs

#### Parameters

This is the main view of the software. Here each parameter can be
individually tuned. Here only information relevant to the usage of the
GUI will be provided about the parameters. For a more comprehensive
explanation of how they are incorporated into the projection see XXX.

##### Unit Cell Scale Factor

This factor determines the size of the primitive unit cell. Each
surface type in the software has a specific primitive unit cell with a
fixed proportion of edge lengths in each dimension. Structures with a
cubic symmtery must have a cubic unit cell. This parameter allows to
scale that intrinsic primitive unit cell size to any desired
value. Note, that for surfaces with cubic symmetries the unit cell
scale factor in z direction must be identical to the on in xy
direction to keep a cubic unit cell. It is thus deactivated for
structures with cubic symmetries.

##### Surface control parameter

The position of the inital membrane can be tuned using the
mathematical level-set value, or the volume proportion of the two
channels. The value here is interpreted as either, depending on the
choice made here. Note that the selection here will carry over to the
batch creation tab, i.e. if the projection is tuned using the volume
proportion value, so will be the projections in the batch creation
tab.


##### Slice dimensions

These parameters control the dimensions in the slice and are given in
length units.

##### Orientation

The hkl value provide the orientation of the slice. The small sketch
next to it visualises the orientation of the slice inside a single
unit cell of the structure. The red line is the normal vector, the
green and yello lines represent the inplane lattice vectors. (see XXX
for details)

##### Resolution

WARNING: chosing a high resolution may result in a extremely high
deman in memory an may crash your computer. So please be aware and
reasonable to choose the resolution


##### Autoupdate

If this box is ticked, the software will automatically recompute the
projection as soon as a parameter is changed. Note, that we strongly
discourage to use this feature for high resolutions, since it wil
result in a slow and poor user experience.

##### Render Parameters

If checked, all parameters used to create this projection will be
written into a margin in the saved pixel image of the projection

#### Measurements

This tabs contains a list of measurement about the current
structure. The buttons on the bottom allow to compute some
quantities. Depending on the choice of the dropdown menue on the
bottom left it computes the latte for the entire slice or just for a
single, primitive unit cell.

Note, that the measurements can only be accurate down to a couple of
voxel sizes. As a help to judge this, the dimension of the voxels used
to compute the quantities are listed below the tables containing the
measurement quantities.

The software does *NOT* use the resolution from the parameters tab,
but rather a predefined resolution for which the algorithm have the
best cost-yield efficancy. They may be tuned in the config file found
in the same directoy as the binary.

