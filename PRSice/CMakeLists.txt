cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(PRSice)
include(ExternalProject)
set(CMAKE_BUILD_TYPE Release)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# To compile a static library on Rosalind, we can use 
# g++ --std=c++11 -I ../inc/ -isystem ../lib/ -L /usr/lib/x86_64-redhat-linux5E/lib64 -DNDEBUG -Wl,--whole-archive -lpthread -Wl,--no-whole-archive  -static-libstdc++  -static-libgcc -static ../src/*.cpp -o PRSice

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
# Eigen from http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
set(EIGEN3_INCLUDE_DIR ${LIBRARY_OUTPUT_PATH}/)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# Header files
include_directories(inc)

file(GLOB SOURCES "src/*.cpp")
find_package (Threads)

add_executable(PRSice ${SOURCES})
target_link_libraries (PRSice ${CMAKE_THREAD_LIBS_INIT})
target_compile_features(PRSice PRIVATE cxx_range_for)

