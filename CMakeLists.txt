cmake_minimum_required (VERSION 2.6.4)
project (surface_projection)
set (CMAKE_CXX_STANDARD 11) # needed for Mac OS X


#add an option to build the gui or only terminal appliaction
option(BUILD_WX_GUI "Builds the GUI (wxWidgets) instead of terminal interface" OFF)
option(BUILD_QT_GUI "Builds the GUI (QT) instead of terminal interface" OFF)

#add an option to build the gui or only terminal appliaction
option(USE_CGAL "Includes CGAL to compute the surface area of the membranes" OFF)

#add an option to build the gui or only terminal appliaction
option(BUILD_COMPUTE_CHANNEL_PROPORTIONS "Builds a tool to compute volume channel proportions" OFF)

#source files independent of interface
SET(SRC_BASE img_out.cpp surface_projection.cpp distance_transform.cpp homotopic_thinning.cpp iterable_voxel.cpp)

#the source files for the wxWidgets gui
SET(SRC_WX_GUI sp_app.cpp sp_gui.cpp sp_gui_imp.cpp)

#the source files for the QT gui
SET(SRC_QT_GUI qt_gui.cpp qt_main.cpp qt_sp.cpp)

#source files for terminal interface
SET(SRC_TER minimal_terminal.cpp)


if(BUILD_WX_GUI)
  #build the gui, bring in all necessary libraries

  set(wxWidgets_CONFIGURATION mswu)
  find_package(wxWidgets COMPONENTS core base REQUIRED)
  include(${wxWidgets_USE_FILE})
  
  add_executable(projection_tool_wx ${SRC_BASE} ${SRC_WX_GUI})

  if( USE_CGAL )
    find_package(CGAL)
    if( CGAL_FOUND )
      add_compile_definitions(USE_CGAL)
      target_link_libraries(projection_tool_wx CGAL::CGAL)
    else()
      message("CGAL not found, can't build the executable!")
    endif()
  endif( USE_CGAL )

  target_link_libraries(projection_tool_wx ${wxWidgets_LIBRARIES})
endif(BUILD_WX_GUI)


#  add_executable(projection_tool ${SRC_BASE} ${SRC_TER})

if(BUILD_QT_GUI)
  #build the gui, bring in all necessary libraries

  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTORCC ON)
  set(CMAKE_AUTOUIC ON)

  find_package(Qt5 COMPONENTS Widgets Gui REQUIRED)
  
  add_executable(projection_tool_qt ${SRC_BASE} ${SRC_QT_GUI})

  if( USE_CGAL )
    find_package(CGAL)
    if( CGAL_FOUND )
      add_compile_definitions(USE_CGAL)
      target_link_libraries(projection_tool_qt CGAL::CGAL)
    else()
      message("CGAL not found, can't build the executable!")
    endif()
  endif( USE_CGAL )

  target_link_libraries(projection_tool_qt Qt5::Widgets Qt5::Gui)
endif(BUILD_QT_GUI)
  


if(BUILD_COMPUTE_CHANNEL_PROPORTIONS)
    add_executable(get_channel_proportions ${SRC_BASE} get_channel_proportions.cpp)  
endif(BUILD_COMPUTE_CHANNEL_PROPORTIONS)

      

