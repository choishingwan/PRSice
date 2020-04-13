cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(PROJECT_NAME PRSice)
project(${PROJECT_NAME} CXX)

add_compile_options(-g)
add_compile_options(-Wall)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

option(march "Use --march." OFF)
if(march)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

find_path(EIGEN_INCLUDE_DIR
    NAME EIGEN
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen-git-mirror/)
if((NOT ${EIGEN_INCLUDE_DIR}) OR (NOT EXISTS ${EIGEN_INCLUDE_DIR}))
    execute_process(COMMAND git submodule update --init -- lib/eigen/
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set(EIGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen/)
endif()
include_directories(${EIGEN_INCLUDE_DIR})
################################
#           Add zlib
################################
find_package( ZLIB REQUIRED )
# if found, will set ${ZLIB_INCLUDE_DIRS} which can be added
################################
#          Add pthread
################################
find_package (Threads REQUIRED)
# if found, will include ${CMAKE_THREAD_LIBS_INIT}

add_library(coverage_config INTERFACE)
option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER 3.13 OR CMAKE_VERSION VERSION_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif()

add_subdirectory(src)

option (BUILD_TESTING "Build the unit test." OFF)
# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
    enable_testing()
    add_subdirectory(test)
endif()
