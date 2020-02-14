cmake_minimum_required (VERSION 2.6.4)
project (surface_projection)
set (CMAKE_CXX_STANDARD 11) # needed for Mac OS X



#add an option to build the gui or only terminal appliaction
option(BUILD_GUI "Builds the GUI instead of terminal interface" ON)

#add an option to build the gui or only terminal appliaction
option(BUILD_TABLEBUILDER "Builds a tool to generate surface/volume lookup tables" OFF)

#source files independent of interface
SET(SRC_BASE img_out.cpp surface_projection.cpp distance_transform.cpp homotopic_thinning.cpp iterable_voxel.cpp)

#the source files for the gui
SET(SRC_GUI sp_app.cpp sp_gui.cpp sp_gui_imp.cpp)

#source files for terminal interface
SET(SRC_TER minimal_terminal.cpp)


if(BUILD_GUI)  
  #build the gui, bring in all necessary libraries

  set(wxWidgets_CONFIGURATION mswu)
  find_package(wxWidgets COMPONENTS core base REQUIRED)
  include(${wxWidgets_USE_FILE})
  
  add_executable(projection_tool ${SRC_BASE} ${SRC_GUI})
  target_link_libraries(projection_tool ${wxWidgets_LIBRARIES})
else(BUILD_GUI)
  add_executable(projection_tool ${SRC_BASE} ${SRC_TER})
endif(BUILD_GUI)
  

if(BUILD_TABLEBUILDER)
  cmake_policy(SET CMP0028 NEW)
  
  find_package(CGAL)
  if( CGAL_FOUND )
    add_executable(lookuptable_builder ${SRC_BASE} lookuptable_builder.cpp)  
    target_link_libraries(lookuptable_builder CGAL::CGAL)
  else()
    message("CGAL not found, can't build the executable!")
  endif()
endif(BUILD_TABLEBUILDER)

      

