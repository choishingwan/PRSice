#ifndef MISC_TEST_HPP
#define MISC_TEST_HPP
#include "reporter.hpp"
#include "global.hpp"
#include "gtest/gtest.h"
#include <vector>

TEST(REPORTER, FILE_NOT_WRITABLE){
    try {
        // shouldn't be able to write in the base file
        Reporter reporter(std::string("/LOG"));
        FAIL();
    } catch (...) {
        SUCCEED();
    }
}

TEST(REPORTER, CHANGE_WIDTH){
    try {
        // shouldn't be able to write in the base file
        Reporter reporter(std::string(path+"LOG"), 80);
        FAIL();
    } catch (...) {
        SUCCEED();
    }
}


TEST(REPORTER, INIT){
    try {
        // initialize with nothing
        Reporter reporter;
        FAIL();
    } catch (...) {
        SUCCEED();
    }
}
#endif //MISC_TEST_HPP
