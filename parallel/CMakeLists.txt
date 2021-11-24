cmake_minimum_required(VERSION 3.10)

# set the project name
project(SCH VERSION 1.0)

# set the main cmake compilation variables
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(SCH_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

# include files
include_directories(${SCH_SOURCE_DIR})

# glob source files into one variable
file(GLOB SCH_SOURCE_FILES ${SCH_SOURCE_DIR}/*.cpp)

# add compiler options
add_compile_options(-Wall -std=c++11 -DVERBOSITY=1)

# MPI support
find_package( MPI COMPONENTS CXX )
include_directories( ${MPI_INCLUDE_PATH} )

# add the executable
add_executable(sch ${SCH_SOURCE_FILES})

target_link_libraries(sch ${MPI_CXX_LIBRARIES})
