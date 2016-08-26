cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(PRSice)
set(CMAKE_BUILD_TYPE Release)

# Header files
include_directories(inc)

file(GLOB SOURCES "src/*.cpp")

add_executable(PRSice ${SOURCES})
target_compile_features(PRSice PRIVATE cxx_range_for)

