cmake_minimum_required (VERSION 3.1)
################################
#          Add EIGEN
################################
# Eigen from http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/inc)
SET(prsice_header
    binarygen.hpp
    binaryplink.hpp
    commander.hpp
    enumerators.h
    dcdflib.h
    family.hpp
    fastlm.hpp
    genotype.hpp
    genotypefactory.hpp
    glm.hpp
    memoryread.hpp
    misc.hpp
    prsice.hpp
    region.hpp
    regression.hpp
    reporter.hpp
    snp.hpp
    storage.hpp
    thread_queue.hpp)

add_library(bgen ${CMAKE_SOURCE_DIR}/src/bgen_lib.cpp)
add_library(gzstream ${CMAKE_SOURCE_DIR}/src/gzstream.cpp)
add_library(plink
    ${CMAKE_SOURCE_DIR}/src/plink_common.cpp
    ${CMAKE_SOURCE_DIR}/src/dcdflib.cpp)
add_library(prsice_lib
    ${CMAKE_SOURCE_DIR}/src/binarygen.cpp
    ${CMAKE_SOURCE_DIR}/src/binaryplink.cpp
    ${CMAKE_SOURCE_DIR}/src/commander.cpp
    ${CMAKE_SOURCE_DIR}/src/fastlm.cpp
    ${CMAKE_SOURCE_DIR}/src/genotype.cpp
    ${CMAKE_SOURCE_DIR}/src/misc.cpp
    ${CMAKE_SOURCE_DIR}/src/prsice.cpp
    ${CMAKE_SOURCE_DIR}/src/region.cpp
    ${CMAKE_SOURCE_DIR}/src/regression.cpp
    ${CMAKE_SOURCE_DIR}/src/reporter.cpp
    ${CMAKE_SOURCE_DIR}/src/snp.cpp)
target_link_libraries( prsice_lib PRIVATE
    bgen
    gzstream
    plink
    coverage_config)

add_executable(PRSice ${CMAKE_SOURCE_DIR}/src/main.cpp)
################################
#           Add zlib
################################
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( PRSice PRIVATE ${ZLIB_LIBRARIES} )
    target_link_libraries( bgen PRIVATE ${ZLIB_LIBRARIES} )
    target_link_libraries( plink PRIVATE ${ZLIB_LIBRARIES} )
    target_link_libraries( gzstream PRIVATE ${ZLIB_LIBRARIES} )
    target_link_libraries( prsice_lib PRIVATE ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
################################
#          Add pthread
################################
find_package (Threads)
target_link_libraries (PRSice PRIVATE ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries (plink PRIVATE ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries( prsice_lib PRIVATE ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries(PRSice PRIVATE
    bgen
    gzstream
    plink
    prsice_lib
    coverage_config)

