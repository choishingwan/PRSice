#ifndef BIN_GEN_TEST_H
#define BIN_GEN_TEST_H
#include "binarygen.hpp"
#include "global.hpp"
#include "gtest/gtest.h"
#include <string>

/*
TEST(SAMPLE_FILE_CHECK, CHECK_SAMPLE_FILE)
{
    ASSERT_TRUE(
        BinaryGen::check_is_sample_format(std::string(path + "valid.sample")));
    ASSERT_TRUE(BinaryGen::check_is_sample_format(
        std::string(path + "valid.sample.sex")));
    ASSERT_TRUE(BinaryGen::check_is_sample_format(
        std::string(path + "valid.sample.long")));
    // This contain only the phenotype column, which violates the bgen format
    ASSERT_FALSE(BinaryGen::check_is_sample_format(
        std::string(path + "invalid.sample")));
    // This is considered as invalid by PRSice because it doesn't have 3 column
    // of 0, but according to latest bgen specification, this might be a valid
    // format
    ASSERT_FALSE(BinaryGen::check_is_sample_format(
        std::string(path + "invalid.sample.special")));
    // Contain unknow type character
    ASSERT_FALSE(BinaryGen::check_is_sample_format(
        std::string(path + "invalid.sample.unknown")));
}*/
#endif
