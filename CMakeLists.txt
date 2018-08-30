cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(PRSice)
include(ExternalProject)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds
# the following targets: gtest, gtest_main, gmock
# and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build")

# The gtest/gmock targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
                        "${gmock_SOURCE_DIR}/include")
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

file(GLOB HEADERS "test/inc/*.hpp")
set(TEST_INPUT "test/data/")
include_directories( ${PROJECT_SOURCE_DIR}/test/inc )
# We don't want to include the PRSice main to our test program
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/main.cpp")
add_executable(PRSice_Test ${PROJECT_SOURCE_DIR}/test/main.cpp
    ${HEADERS}
    ${SOURCES})
target_include_directories(PRSice_Test PRIVATE ${PROJECT_SOURCE_DIR}/inc)
target_link_libraries(PRSice_Test gtest gtest_main)
target_link_libraries( PRSice_Test ${ZLIB_LIBRARIES} )
target_link_libraries( PRSice_Test  ${CMAKE_THREAD_LIBS_INIT} )
target_compile_features(PRSice_Test PRIVATE cxx_range_for)
enable_testing()
add_test(PRSice_Test COMMAND PRSice_Test ${TEST_INPUT})
