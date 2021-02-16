cmake_minimum_required (VERSION 2.6.4)
project (surface_projection)
set (CMAKE_CXX_STANDARD 11) # needed for Mac OS X


#add an option to build the gui or only terminal appliaction
option(BUILD_QT_GUI "Builds the GUI (QT) instead of terminal interface" ON)

#add an option to build the gui or only terminal appliaction
option(BUILD_COMPUTE_CHANNEL_PROPORTIONS "Builds a tool to compute volume channel proportions" OFF)

#add an option to build the gui or only terminal appliaction
option(BUILD_BATCH_TOOL "Builds a tool to generate batches of files" OFF)

#source files independent of interface
SET(SRC_BASE img_out.cpp surface_projection.cpp distance_transform.cpp homotopic_thinning.cpp iterable_voxel.cpp percolation_analysis.cpp batch_lib.cpp)

#the source files for the QT gui
SET(SRC_QT_GUI qt_gui.cpp qt_main.cpp qt_sp.cpp slice_orientation_visualisation.cpp drawme.cpp)

#
# external libraries
#
SET(IML_DIR ./ext_libs/iml-1.0.5/src/)
SET(SRC_IML ${IML_DIR}/nullspace.c ${IML_DIR}/basisop.c /
  ${IML_DIR}/certsolve.c ${IML_DIR}/error.c /
  ${IML_DIR}/latreduce.c ${IML_DIR}/memalloc.c /
  ${IML_DIR}/mtrans.c ${IML_DIR}/nonsysolve.c /
  ${IML_DIR}/nullspace.c ${IML_DIR}/padiclift.c /
  ${IML_DIR}/reconsolu.c ${IML_DIR}/RNSop.c )
include_directories(./ext_libs/iml-1.0.5/src)

if(BUILD_QT_GUI)
  #build the gui, bring in all necessary libraries

  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTORCC ON)
  set(CMAKE_AUTOUIC ON)

  find_package(Qt5 COMPONENTS Widgets Gui Svg REQUIRED)
  find_library(gmp gmp)
  find_package(BLAS)  
  
  add_executable(projection_tool_qt ${SRC_BASE} ${SRC_QT_GUI} ${SRC_IML} resources.qrc)


  #check if we have CGAL
  find_package(CGAL REQUIRED)    
  if( CGAL_FOUND )
    add_compile_definitions(USE_CGAL)
    target_link_libraries(projection_tool_qt CGAL::CGAL)
  endif()

  find_package(PNG REQUIRED)
  if( PNG_FOUND )
    add_compile_definitions(HAVE_PNG)
    target_link_libraries(projection_tool_qt PNG::PNG)
    target_include_directories(projection_tool_qt PRIVATE PNG_INCLUDE_DIRS)
  endif()
  
  target_link_libraries(projection_tool_qt BLAS::BLAS gmp)
  
  target_link_libraries(projection_tool_qt Qt5::Widgets Qt5::Gui Qt5::Svg)
endif(BUILD_QT_GUI)


if(BUILD_COMPUTE_CHANNEL_PROPORTIONS)
  add_executable(get_channel_proportions ${SRC_BASE} get_channel_proportions.cpp)  
endif(BUILD_COMPUTE_CHANNEL_PROPORTIONS)


if( BUILD_BATCH_TOOL )

  add_executable(batch_tool ${SRC_BASE} ${SRC_IML} batch_tool.cpp resources.qrc)

  find_library(gmp gmp)
  find_package(BLAS)  
  find_package(CGAL REQUIRED)

  find_package(PNG REQUIRED)

  if( PNG_FOUND )
    add_compile_definitions(HAVE_PNG)
    target_link_libraries(batch_tool PNG::PNG)
    target_include_directories(batch_tool PRIVATE PNG_INCLUDE_DIRS)
  endif()
  
  target_link_libraries(batch_tool BLAS::BLAS gmp)

  if( CGAL_FOUND )
    add_compile_definitions(USE_CGAL)
    target_link_libraries(batch_tool CGAL::CGAL)
  endif()

endif( BUILD_BATCH_TOOL )
