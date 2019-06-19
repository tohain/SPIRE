cmake_minimum_required (VERSION 2.6.4)
project (surface_projection)

#Collects the names of all the source files in the specified directory
#and stores the list in the <variable> provided
aux_source_directory(${PROJECT_SOURCE_DIR} SRC_LIST)

#add the executable
add_executable(${PROJECT_NAME} ${SRC_LIST})

      

