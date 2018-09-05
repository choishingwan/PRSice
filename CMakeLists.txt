cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(PROJECT_NAME PRSice)
project(${PROJECT_NAME})
set(CMAKE_CXX_FLAGS "-g -Wall")

option(march "Use --march." OFF)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_CXX_FLAGS "-Wall -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
if(march)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()


################################
#          Add EIGEN
################################
# Eigen from http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
set(EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/lib)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

################################
# Normal Libraries & Executables
################################

include_directories(${CMAKE_SOURCE_DIR}/inc)
file(GLOB SOURCES "src/*.cpp")
LIST(REMOVE_ITEM SOURCES
    ${CMAKE_SOURCE_DIR}/src/main.cpp)
add_library(prsice_lib ${SOURCES})
add_executable(PRSice ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_compile_features(PRSice PRIVATE cxx_range_for)
target_compile_features(prsice_lib PRIVATE cxx_range_for)
################################
#           Add zlib
################################
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( PRSice ${ZLIB_LIBRARIES} )
    target_link_libraries( prsice_lib ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
################################
#          Add pthread
################################
find_package (Threads)
target_link_libraries (PRSice ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries (prsice_lib ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(PRSice prsice_lib)




################################
#  Googletest configuration
################################
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



# Header files

file(GLOB SOURCES "src/*.cpp")
include_directories( ${PROJECT_SOURCE_DIR}/test/inc )
# We don't want to include the PRSice main to our test program
add_executable(runUnitTests ${PROJECT_SOURCE_DIR}/test/main.cpp)
target_include_directories(runUnitTests PRIVATE ${PROJECT_SOURCE_DIR}/inc)
target_link_libraries(runUnitTests prsice_lib)
target_link_libraries(runUnitTests gtest gtest_main)
target_link_libraries( runUnitTests ${ZLIB_LIBRARIES} )
target_link_libraries( runUnitTests  ${CMAKE_THREAD_LIBS_INIT} )
target_compile_features(runUnitTests PRIVATE cxx_range_for)
enable_testing()
add_test(runUnitTests COMMAND PRSice_Test ${TEST_INPUT})
