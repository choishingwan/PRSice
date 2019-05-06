cmake_minimum_required (VERSION 3.1)
################################
#          Add EIGEN
################################
# Eigen from http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/inc)

add_library(bgen ${CMAKE_SOURCE_DIR}/src/bgen_lib.cpp)
target_link_libraries(bgen PRIVATE coverage_config)
add_library(binarygen ${CMAKE_SOURCE_DIR}/src/binarygen.cpp)
target_link_libraries(binarygen PRIVATE coverage_config)
target_link_libraries(bgen PRIVATE coverage_config)
add_library(binaryplink ${CMAKE_SOURCE_DIR}/src/binaryplink.cpp)
target_link_libraries(binaryplink PRIVATE coverage_config)
add_library(commander ${CMAKE_SOURCE_DIR}/src/commander.cpp)
target_link_libraries(commander PRIVATE coverage_config)
add_library(genotype ${CMAKE_SOURCE_DIR}/src/genotype.cpp)
target_link_libraries(genotype PRIVATE coverage_config)
add_library(gzstream ${CMAKE_SOURCE_DIR}/src/gzstream.cpp)
add_library(misc ${CMAKE_SOURCE_DIR}/src/misc.cpp)
target_link_libraries(misc PRIVATE coverage_config)
add_library(plink ${CMAKE_SOURCE_DIR}/src/plink_common.cpp)
add_library(prsice_class ${CMAKE_SOURCE_DIR}/src/prsice.cpp)
target_link_libraries(prsice_class PRIVATE coverage_config)
add_library(region ${CMAKE_SOURCE_DIR}/src/region.cpp)
target_link_libraries(region PRIVATE coverage_config)
add_library(regression ${CMAKE_SOURCE_DIR}/src/regression.cpp)
target_link_libraries(regression PRIVATE coverage_config)
add_library(reporter ${CMAKE_SOURCE_DIR}/src/reporter.cpp)
target_link_libraries(reporter PRIVATE coverage_config)
add_library(snp ${CMAKE_SOURCE_DIR}/src/snp.cpp)
target_link_libraries(snp PRIVATE coverage_config)

add_executable(PRSice ${CMAKE_SOURCE_DIR}/src/main.cpp)
################################
#           Add zlib
################################
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( PRSice PUBLIC ${ZLIB_LIBRARIES} )
    target_link_libraries( bgen PUBLIC ${ZLIB_LIBRARIES} )
    target_link_libraries( plink PUBLIC ${ZLIB_LIBRARIES} )
    target_link_libraries( genotype PRIVATE ${ZLIB_LIBRARIES} )
    target_link_libraries( commander PRIVATE ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
################################
#          Add pthread
################################
find_package (Threads)
target_link_libraries (PRSice PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(PRSice PRIVATE
    bgen
    binarygen
    binaryplink
    commander
    genotype
    gzstream
    misc
    plink
    prsice_class
    region
    regression
    reporter
    snp
    coverage_config)

