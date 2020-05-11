cmake_minimum_required (VERSION 3.1)
################################
#          Add EIGEN
################################
# Eigen from http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2

# bgen
add_library(bgen
    ${CMAKE_SOURCE_DIR}/src/bgen_lib.cpp)
target_include_directories(bgen SYSTEM PUBLIC
    ${ZLIB_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries(bgen ${ZLIB_LIBRARIES})
# gzstream
add_library(gzstream
    ${CMAKE_SOURCE_DIR}/src/gzstream.cpp)
target_include_directories(gzstream SYSTEM PUBLIC
    ${ZLIB_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(gzstream ${ZLIB_LIBRARIES})


# Useful helpers
add_library(utility
    ${CMAKE_SOURCE_DIR}/src/misc.cpp
    ${CMAKE_SOURCE_DIR}/src/commander.cpp
    ${CMAKE_SOURCE_DIR}/src/reporter.cpp)
target_include_directories(utility PUBLIC
    ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries(utility PUBLIC
    gzstream
    coverage_config)

# plink
add_library(plink
    ${CMAKE_SOURCE_DIR}/src/plink_common.cpp
    ${CMAKE_SOURCE_DIR}/src/dcdflib.cpp)
target_include_directories(plink SYSTEM PUBLIC
    ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(plink PUBLIC utility)


add_library(genotyping
    ${CMAKE_SOURCE_DIR}/src/binarygen.cpp
    ${CMAKE_SOURCE_DIR}/src/binaryplink.cpp
    ${CMAKE_SOURCE_DIR}/src/genotype.cpp
    ${CMAKE_SOURCE_DIR}/src/genotype_pool.cpp
    ${CMAKE_SOURCE_DIR}/src/snp.cpp)
target_include_directories(genotyping PUBLIC
    ${CMAKE_SOURCE_DIR}/inc)
target_include_directories(genotyping SYSTEM PUBLIC
    ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(genotyping PUBLIC
    plink
    utility
    bgen
    coverage_config)


add_library(regression
    ${CMAKE_SOURCE_DIR}/src/fastlm.cpp
    ${CMAKE_SOURCE_DIR}/src/regression.cpp)
target_include_directories(regression PUBLIC
    ${CMAKE_SOURCE_DIR}/inc)
target_include_directories(regression SYSTEM PUBLIC
    ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(regression PUBLIC plink)


add_library(prsice_lib
    ${CMAKE_SOURCE_DIR}/src/prset.cpp
    ${CMAKE_SOURCE_DIR}/src/prsice.cpp
    ${CMAKE_SOURCE_DIR}/src/region.cpp)
target_include_directories(prsice_lib INTERFACE
    ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries(prsice_lib PUBLIC
    genotyping
    regression
    utility
    ${CMAKE_THREAD_LIBS_INIT}
    coverage_config)


add_executable(PRSice
    ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_include_directories(PRSice PUBLIC
    ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries(PRSice PUBLIC
    bgen
    gzstream
    plink
    prsice_lib
    genotyping
    regression
    utility
    coverage_config)

