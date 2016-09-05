cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(PRSice)
include(ExternalProject)
set(CMAKE_BUILD_TYPE Release)

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
#Eigen from http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
set(EIGEN3_INCLUDE_DIR ${LIBRARY_OUTPUT_PATH}/)
include_directories(${EIGEN3_INCLUDE_DIR})

# Header files
include_directories(inc)

file(GLOB SOURCES "src/*.cpp")

add_executable(PRSice ${SOURCES})
target_compile_features(PRSice PRIVATE cxx_range_for)

