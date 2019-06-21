cmake_minimum_required (VERSION 2.6.4)
project (surface_projection)
set (CMAKE_CXX_STANDARD 11) # needed for Mac OS X

#Collects the names of all the source files in the specified directory
#and stores the list in the <variable> provided
aux_source_directory(${PROJECT_SOURCE_DIR} SRC_LIST)



set(wxWidgets_CONFIGURATION mswu)
find_package(wxWidgets COMPONENTS core base REQUIRED)
include(${wxWidgets_USE_FILE})

#add_executable(${PROJECT_NAME} main.cpp)

add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${wxWidgets_LIBRARIES})
      

