#ifndef MISC_TEST_HPP
#define MISC_TEST_HPP
#include "global.hpp"
#include "reporter.hpp"
#include "gtest/gtest.h"
#include <vector>

TEST(REPORTER, CHANGE_WIDTH)
{
    try
    {
        // shouldn't be able to write in the base file
        Reporter reporter(std::string(path + "LOG"), 80);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(REPORTER, FILE_WRITE_FAILED)
{
    try
    {
        // shouldn't be able to write in the base file
        Reporter reporter("", 80);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}


TEST(REPORTER, EMPTY_INIT)
{
    try
    {
        // initialize with nothing
        Reporter reporter;
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(REPORTER, INITIALIZE)
{
    try
    {
        // initialize with nothing
        Reporter reporter;
        reporter.initiailize(std::string(path + "LOC"));
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(REPORTER, FAIL_INITIALIZE)
{
    try
    {
        // initialize with nothing
        Reporter reporter;
        reporter.initiailize("");
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}


TEST(REPORTER, REPORTING_MESSAGE)
{
    try
    {
        // initialize with nothing
        Reporter reporter(std::string(path + "LOG"));
        reporter.report("OUTPUT");
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(REPORTER, LIST_MESSAGE)
{
    try
    {
        // initialize with nothing
        Reporter reporter(std::string(path + "LOG"));
        reporter.report("1) Testing\n2)If this is ok\n");
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}

#endif // MISC_TEST_HPP
