cmake_minimum_required(VERSION 3.1 FATAL_ERROR)


include_directories(${CMAKE_SOURCE_DIR}/test/)
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/inc)
set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/test)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# Make test executable
set(TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/catch-main.cpp)
add_executable(tests ${TEST_SOURCES}
    csrc/commander_test.cpp)
target_link_libraries(tests PRIVATE
    Catch
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
    target_link_libraries( tests PUBLIC ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
################################
#          Add pthread
################################
find_package (Threads)
target_link_libraries (tests PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries( tests  PUBLIC ${CMAKE_THREAD_LIBS_INIT} )
target_compile_features(tests PRIVATE cxx_range_for)
add_test(NAME unitTest COMMAND runUnitTests "${CMAKE_CURRENT_LIST_DIR}/test/data/")
