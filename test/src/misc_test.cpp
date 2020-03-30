#ifndef MISC_TEST_HPP
#define MISC_TEST_HPP
#include "global.hpp"
#include "gzstream.h"
#include "misc.hpp"
#include "reporter.hpp"
#include "gtest/gtest.h"
#include <string_view>
#include <vector>

TEST(HAS_ENDING, CHECK_VALIDITY)
{
    std::string full_text;
    std::string ending;
    try
    {
        misc::hasEnding(full_text, ending);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    full_text = "Hello";
    try
    {
        misc::hasEnding(full_text, ending);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    ending = "lo";
    ASSERT_TRUE(misc::hasEnding(full_text, ending));
    ending = "Lo";
    ASSERT_FALSE(misc::hasEnding(full_text, ending));
    ending = "Hello";
    ASSERT_TRUE(misc::hasEnding(full_text, ending));
    ending = "HEllO";
    ASSERT_FALSE(misc::hasEnding(full_text, ending));
}

TEST(SPLIT, STRING_VIEW_FORMAT)
{
    std::string input = "test,the,splitter,and,make,sure-it,works well";
    std::vector<std::string_view> token = misc::tokenize(input, ",");
    ASSERT_EQ(token.size(), 7);
    ASSERT_EQ(token[0], "test");
    ASSERT_EQ(token[1], "the");
    ASSERT_EQ(token[2], "splitter");
    ASSERT_EQ(token[3], "and");
    ASSERT_EQ(token[4], "make");
    ASSERT_EQ(token[5], "sure-it");
    ASSERT_EQ(token[6], "works well");
    input = "sure-it\tworks well ok";
    token = misc::tokenize(input);
    ASSERT_EQ(token.size(), 4);
    ASSERT_EQ(token[0], "sure-it");
    ASSERT_EQ(token[1], "works");
    ASSERT_EQ(token[2], "well");
    ASSERT_EQ(token[3], "ok");
}
TEST(OVERFLOW_CHECK, OVERFLOW_CHECK)
{
    ASSERT_FALSE(misc::overflow(10, 1));
    ASSERT_FALSE(misc::overflow(1000000000, 0));
}

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
        reporter.initialize(std::string(path + "LOC"));
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
        reporter.initialize("");
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
        Reporter reporter(std::string(path + "LOG"), 60, true);
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
        Reporter reporter(std::string(path + "LOG"), 60, true);
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

TEST(CONVERTOR, CONVERTOR)
{
    size_t res = misc::Convertor::convert<size_t>("10000");
    ASSERT_EQ(res, 10000);
    double res_db = misc::Convertor::convert<double>("0.05");
    ASSERT_DOUBLE_EQ(res_db, 0.05);
    try
    {
        res = misc::Convertor::convert<size_t>("-10000");
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    try
    {
        std::cerr << std::numeric_limits<double>::max() << "\t"
                  << std::numeric_limits<double>::min() << std::endl;
        res_db = misc::Convertor::convert<double>("1e-400");
        std::cerr << "Converted into: " << res_db << std::endl;
        std::cerr << "Check erange: " << ERANGE << std::endl;
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
}
TEST(TRIM, STRING_VIEW)
{
    std::string ref = " testing \n";
    std::string_view s = ref;
    std::string_view t = ref;
    s = misc::trimmed(s);
    ASSERT_EQ(s, "testing");
    ASSERT_STREQ(ref.c_str(), " testing \n");
    misc::trim(t);
    ASSERT_EQ(t, "testing");
}
#endif // MISC_TEST_HPP
