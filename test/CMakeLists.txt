cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

add_library(Catch INTERFACE)
set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/test/inc)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})


set(TEST_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/test/inc)
set(TEST_SRC_DIR ${CMAKE_SOURCE_DIR}/test/csrc)

# Make test executable
set(TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/catch-main.cpp)
add_executable(tests ${TEST_SOURCES}
    ${TEST_SRC_DIR}/commander_test.cpp
    ${TEST_SRC_DIR}/command_loading.cpp
    ${TEST_SRC_DIR}/command_validation.cpp
    ${TEST_SRC_DIR}/misc_test.cpp
    ${TEST_SRC_DIR}/genotype_basic.cpp
    ${TEST_SRC_DIR}/genotype_read_base.cpp
    ${TEST_SRC_DIR}/genotype_read_sample.cpp
    ${TEST_SRC_DIR}/genotype_load_snp.cpp
    ${TEST_SRC_DIR}/snp_test.cpp
    ${TEST_SRC_DIR}/binaryplink_sample_load.cpp
    ${TEST_SRC_DIR}/binaryplink_snp_load.cpp
    ${TEST_SRC_DIR}/binaryplink_filtering.cpp
    ${TEST_SRC_DIR}/binarygen_sample_load.cpp
    ${TEST_SRC_DIR}/binarygen_snp_load.cpp
    ${TEST_SRC_DIR}/binarygen_filtering.cpp
    ${TEST_SRC_DIR}/region_basic.cpp
    ${TEST_SRC_DIR}/region_exclusion.cpp
    ${TEST_SRC_DIR}/region_process.cpp
    ${TEST_SRC_DIR}/prsice_pheno.cpp
    )
target_link_libraries(tests PUBLIC
    Catch
    genotyping
    prsice_lib
    plink
    utility
    coverage_config)

add_test(NAME unitTest COMMAND tests)

add_custom_command(
     TARGET tests
     COMMENT "Run tests"
     POST_BUILD
     COMMAND tests
)
