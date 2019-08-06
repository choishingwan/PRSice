#ifndef MISC_TEST_HPP
#define MISC_TEST_HPP
#include "global.hpp"
#include "gzstream.h"
#include "misc.hpp"
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


TEST(REPORTER, NO_WRAP_MESSAGE)
{
    try
    {
        // initialize with nothing
        Reporter reporter(std::string(path + "LOG"), 10);
        reporter.report("Should not wrap this message and should allow this to "
                        "go on and on\n",
                        false);
        // need to figure out how to do mocking and actually test if the output
        // is correct. However, reporter is kinda a minor class
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}


TEST(UTILITY, IS_GZ_FILE)
{
    GZSTREAM_NAMESPACE::ogzstream file;
    file.open("DEBUG.gz");
    file.close();
    ASSERT_TRUE(misc::is_gz_file("DEBUG.gz"));
    std::remove("DEBUG.gz");
    file.open("DEBUG");
    file.close();
    ASSERT_TRUE(misc::is_gz_file("DEBUG"));
    std::remove("DEBUG");
    try
    {
        misc::is_gz_file("DEBUG");
        // not found
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    // should be robust against the suffix
    std::ofstream test;
    test.open("DEBUG.gz");
    test.close();
    ASSERT_FALSE(misc::is_gz_file("DEBUG.gz"));
    std::remove("DEBUG.gz");
}
#endif // MISC_TEST_HPP
