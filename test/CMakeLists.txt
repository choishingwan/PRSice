cmake_minimum_required(VERSION 3.1 FATAL_ERROR)


include_directories(${CMAKE_SOURCE_DIR}/test/inc)
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/inc)

add_executable(runUnitTests
    main.cpp
    src/binplink_test.cpp
    src/genotype_test.cpp
    src/misc_test.cpp
    src/region_test.cpp
    src/snp_test.cpp
    src/commander_test.cpp
    src/prsice_test.cpp)
target_link_libraries(runUnitTests PRIVATE
    bgen
    gzstream
    plink
    prsice_lib
    coverage_config)
################################
#           Add zlib
################################
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( runUnitTests PUBLIC ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
################################
#          Add pthread
################################
find_package (Threads)
target_link_libraries (runUnitTests PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(runUnitTests PRIVATE gtest gtest_main)
target_link_libraries( runUnitTests  PUBLIC ${CMAKE_THREAD_LIBS_INIT} )
target_compile_features(runUnitTests PRIVATE cxx_range_for)
add_test(NAME unitTest COMMAND runUnitTests "${CMAKE_CURRENT_LIST_DIR}/test/data/")
