cmake_minimum_required(VERSION 3.1 FATAL_ERROR)


include_directories(${CMAKE_SOURCE_DIR}/test/inc)

add_executable(runUnitTests main.cpp)
target_link_libraries(runUnitTests PRIVATE
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
target_link_libraries(runUnitTests PRIVATE gtest gtest_main)
target_link_libraries( runUnitTests  ${CMAKE_THREAD_LIBS_INIT} )
target_compile_features(runUnitTests PRIVATE cxx_range_for)
add_test(NAME unitTest COMMAND runUnitTests "${CMAKE_CURRENT_LIST_DIR}/test/data/")
