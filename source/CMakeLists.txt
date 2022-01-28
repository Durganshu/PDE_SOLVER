CMAKE_MINIMUM_REQUIRED (VERSION 3.4)

# name of the project.
PROJECT (PDE_SOLVER)

#LINK_DIRECTORIES(usr/include/lib/)

#SET(CMAKE_CXX_FLAGS " -fopenmp -O3")

SET( PDE_CMAKE_DIR "${PDE_SOURCE_DIR}/cmake")
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PDE_CMAKE_DIR}")

# the sources that need to be compiled for the required executable.
SET (sources
    iterative_schemes.cpp
    pde_solver.cpp
    write_plot.cpp
    main.cpp
    )

#INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH}  Eigen/)
#LINK_DIRECTORIES ( ${MPI_LIBRARY_DIR} Eigen/)

# setting the compiler to use.
SET( CMAKE_CXX_COMPILER "g++" )

# general flags for the compiler.
set(CMAKE_CXX_STANDARD 17)

# flags for the compiler in debug mode.
#SET( CMAKE_CXX_FLAGS_DEBUG " -pg -g -Weffc++" )
SET( CMAKE_CXX_FLAGS_DEBUG " -pg -g " )

# flags for the compiler in release mode
SET( CMAKE_CXX_FLAGS_RELEASE " -O3 -s -DNDEBUG" )


#set(EIGEN_DIR "/home/durganshu/pde-solver_-group-dh/source/Eigen/")
#INCLUDE_DIRECTORIES( ${EIGEN_DIR} )
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(pybind11 REQUIRED)

#LINK_DIRECTORIES ( ~/install/lib )
# set the include paths

#target_link_libraries(vec -fsanitize=address)

# Setting the default build type.
IF( NOT CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE Debug )
ENDIF()

###############################################################################
# Setting the target name depending on the build mode.
###############################################################################
IF( CMAKE_BUILD_TYPE MATCHES Release )
  SET( target_name pdesolver )
  SET( compiler_flags ${CMAKE_CXX_FLAGS_RELEASE}  ${CMAKE_CXX_FLAGS} )
  # Add the executable to be generated.
  ADD_EXECUTABLE ( ${target_name} ${sources} )
  TARGET_LINK_LIBRARIES (${target_name} jsoncpp pybind11::embed)
ELSEIF( CMAKE_BUILD_TYPE MATCHES Debug )
  SET( target_name pdesolver_debug )
  SET( compiler_flags ${CMAKE_CXX_FLAGS_DEBUG}  ${CMAKE_CXX_FLAGS} )
  # Add the executable to be generated.
  ADD_EXECUTABLE ( ${target_name} ${sources} )
  TARGET_LINK_LIBRARIES (${target_name} jsoncpp pybind11::embed)
ELSEIF( CMAKE_BUILD_TYPE MATCHES Performance )
  SET( target_name pdesolver_performance )
  SET( compiler_flags ${CMAKE_CXX_FLAGS_RELEASE}  ${CMAKE_CXX_FLAGS} )
  # Add the executable to be generated.
  ADD_EXECUTABLE ( ${target_name} ${sources} )
  TARGET_LINK_LIBRARIES (${target_name} jsoncpp pybind11::embed -pg)
ELSE()
  MESSAGE( "Build type is" ${CMAKE_BUILD_TYPE} )
  MESSAGE( FATAL_ERROR "Build type can be Debug or Release" )
ENDIF()