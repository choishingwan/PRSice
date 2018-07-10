cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(PRSice)
include(ExternalProject)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_CXX_FLAGS "-Wall -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
# Eigen from http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
set(EIGEN3_INCLUDE_DIR ${LIBRARY_OUTPUT_PATH}/)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# Header files
include_directories(inc)

file(GLOB SOURCES "src/*.cpp")
find_package (Threads)

add_executable(PRSice ${SOURCES})
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( PRSice ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
target_link_libraries (PRSice ${CMAKE_THREAD_LIBS_INIT})
target_compile_features(PRSice PRIVATE cxx_range_for)

