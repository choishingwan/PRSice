cmake_minimum_required (VERSION 3.1)
################################
#          Add EIGEN
################################
# Eigen from http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
set(EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/lib)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})



add_library(bgen bgen_lib.cpp)
add_library(binarygen binarygen.cpp)
add_library(binaryplink binaryplink.cpp)
add_library(commander commander.cpp)
add_library(genotype genotype.cpp)
add_library(gzstream gzstream.cpp)
add_library(misc misc.cpp)
add_library(plink plink_common.cpp)
add_library(prsice prsice.cpp)
add_library(region region.cpp)
add_library(regression regression.cpp)
add_library(reporter reporter.cpp)
add_library(snp snp.cpp)
include_directories(${CMAKE_SOURCE_DIR}/inc)


add_executable(PRSice main.cpp)
target_link_libraries(PRSice PRIVATE
    bgen
    binarygen
    binaryplink
    commander
    genotype
    gzstream
    misc
    plink
    prsice
    region
    regression
    reporter
    snp
    coverage_config)
################################
#           Add zlib
################################
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( PRSice PUBLIC ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
################################
#          Add pthread
################################
find_package (Threads)
target_link_libraries (PRSice PUBLIC ${CMAKE_THREAD_LIBS_INIT})

