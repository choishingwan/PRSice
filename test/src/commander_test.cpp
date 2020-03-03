#ifndef COMMANDER_TEST_H
#define COMMANDER_TEST_H
#include "commander.hpp"
#include "global.hpp"
#include "reporter.hpp"
#include "storage.hpp"
#include "gtest/gtest.h"

TEST(COMMANDER_BASIC, INIT)
{
    Commander commander;
    ASSERT_FALSE(commander.all_scores());
    ASSERT_FALSE(commander.use_inter());
    ASSERT_STREQ(commander.delim().c_str(), " ");
    ASSERT_STREQ(commander.out().c_str(), "PRSice");
    ASSERT_TRUE(commander.exclusion_range().empty());
    ASSERT_TRUE(commander.exclude_file().empty());
    ASSERT_TRUE(commander.extract_file().empty());
    // we will use the parameter if memory is not provided, otherwise,
    // we will return the actual memory allowed
    ASSERT_DOUBLE_EQ(commander.max_memory(1.0), 1.0);
    ASSERT_DOUBLE_EQ(commander.max_memory(2.0), 2.0);
}

TEST(COMMANDER_BASIC, USAGE)
{
    Commander commander;
    Reporter reporter(std::string(path + "LOG"), 60, true);
    int argc = 2;
    char name[7], help[7];
    strcpy(name, "PRSice");
    strcpy(help, "--help");
    char* argv[2] = {name, help};
    try
    {
        ASSERT_FALSE(commander.init(argc, argv, reporter));
    }
    catch (...)
    {
        FAIL();
    }
}


TEST(COMMANDER_BASIC, NO_ARG)
{
    Commander commander;
    Reporter reporter(std::string(path + "LOG"), 60, true);
    int argc = 1;
    std::string name = "PRSice";
    char name_c[7];
    strcpy(name_c, name.c_str());
    char* argv[1] = {name_c};
    try
    {
        commander.init(argc, argv, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

class mockCommander : public Commander
{
public:
    static std::vector<std::string>
    transform_covariate(const std::string& cov_in)
    {
        return Commander::transform_covariate(cov_in);
    }
    bool check_parse_unit_value(const std::string& input, const std::string& c,
                                const size_t default_power, size_t& target,
                                bool memory = false)
    {
        return parse_unit_value(input, c, default_power, target, memory);
    }

    static bool find_first_end_wrapper(const std::string_view& cov,
                                       const size_t idx, size_t& res)
    {
        try
        {
            res = find_first_end(cov, idx);
            return true;
        }
        catch (const std::runtime_error&)
        {
            return false;
        }
    }
    static bool parse_range_wrapper(std::string_view cov,
                                    std::vector<size_t>& res)
    {
        try
        {
            res = parse_range(cov);
            return true;
        }
        catch (std::runtime_error&)
        {
            return false;
        }
    }
    static bool get_range_wrapper(std::string_view cov, size_t start,
                                  size_t end, std::vector<size_t>& res)
    {
        try
        {
            res = get_range(cov, start, end);
            return true;
        }
        catch (std::runtime_error&)
        {
            return false;
        }
    }
};


void invalid_cov_input(const std::string& cov_string)
{
    try
    {
        // invalid input
        std::vector<std::string> results =
            mockCommander::transform_covariate(cov_string);
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
}
TEST(COVARIATE_TRANSFORM, RANGE_CHECK)
{

    std::string cov = "PC[1-5]";
    size_t res;
    ASSERT_FALSE(mockCommander::find_first_end_wrapper(cov, 0, res));
    ASSERT_TRUE(mockCommander::find_first_end_wrapper(cov, 2, res));
    ASSERT_EQ(res, 6);
    cov = "PC[1-5[1-5]]";
    ASSERT_FALSE(mockCommander::find_first_end_wrapper(cov, 0, res));
    ASSERT_FALSE(mockCommander::find_first_end_wrapper(cov, 2, res));
}
TEST(COVARIATE_TRANSFORM, PARSE_RANGE)
{
    std::string cov = "[1-5]";
    std::vector<size_t> res;
    // we expect [] to be removed
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
    cov = "1-5";
    res.clear();
    ASSERT_TRUE(mockCommander::parse_range_wrapper(cov, res));
    ASSERT_EQ(res.size(), 5);
    for (size_t i = 0; i < res.size(); ++i) ASSERT_EQ(res[i], i + 1);
    cov = "10-50";
    res.clear();
    ASSERT_TRUE(mockCommander::parse_range_wrapper(cov, res));
    ASSERT_EQ(res.size(), 41);
    for (size_t i = 0; i < res.size(); ++i) ASSERT_EQ(res[i], i + 10);
    cov = "50-10";
    res.clear();
    ASSERT_TRUE(mockCommander::parse_range_wrapper(cov, res));
    ASSERT_EQ(res.size(), 41);
    for (size_t i = 0; i < res.size(); ++i) ASSERT_EQ(res[i], i + 10);
    cov = "10";
    res.clear();
    ASSERT_TRUE(mockCommander::parse_range_wrapper(cov, res));
    ASSERT_EQ(res.size(), 1);
    ASSERT_EQ(res[0], 10);
    // should not work for negative
    cov = "-1";
    res.clear();
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
    cov = "1--5";
    res.clear();
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
    cov = "-1--5";
    res.clear();
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
    res.clear();
    cov = "1,5";
    // we assume , is already dealt with
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
}

TEST(COVARIATE_TRANSFORM, GET_RANGE)
{
    std::string cov = "PC[1-5]";
    std::vector<size_t> res;
    // format cannot be converted
    ASSERT_FALSE(mockCommander::get_range_wrapper(cov, 0, 6, res));
    // still wrong format
    ASSERT_FALSE(mockCommander::get_range_wrapper(cov, 2, 4, res));
    // out of bound
    ASSERT_FALSE(mockCommander::get_range_wrapper(cov, 2, 7, res));
    res.clear();
    ASSERT_TRUE(mockCommander::get_range_wrapper(cov, 2, 6, res));
    ASSERT_EQ(res.size(), 5);
    for (size_t i = 0; i < res.size(); ++i) { ASSERT_EQ(res[i], i + 1); }
    // complex options
    cov = "PC[1-5,8,7-10]";
    res.clear();
    ASSERT_TRUE(mockCommander::get_range_wrapper(cov, 2, 13, res));
    // should be sorted and removed the duplicates (8)
    std::vector<size_t> expected = {1, 2, 3, 4, 5, 7, 8, 9, 10};
    for (size_t i = 0; i < res.size(); ++i) { ASSERT_EQ(res[i], expected[i]); }
    cov = "PC[1-5,-6]";
    // One fail, all fail
    res.clear();
    ASSERT_FALSE(mockCommander::get_range_wrapper(cov, 2, 9, res));
}
/*
TEST(COVARIATE_TRANSFORM, TRANSFORMATION)
{
    // should not do transformation when not start with @
    std::string cov_string = "PC1";
    std::string expected = cov_string;
    ASSERT_STREQ(
        expected.c_str(),
        mockCommander::transform_covariate(cov_string).front().c_str());
    // same for empty string
    cov_string = expected = "";
    ASSERT_STREQ(
        expected.c_str(),
        mockCommander::transform_covariate(cov_string).front().c_str());
    // should be fine if the @ is in middle of the string
    cov_string = expected = "PC1@Home";
    ASSERT_STREQ(
        expected.c_str(),
        mockCommander::transform_covariate(cov_string).front().c_str());
    // when start with @ but not with any [], we will just remove the @
    cov_string = "@PC1";
    expected = "PC1";
    ASSERT_STREQ(
        expected.c_str(),
        mockCommander::transform_covariate(cov_string).front().c_str());
    cov_string = "@PC[1-5]";
    // in this order
    std::vector<std::string> expected_outputs = {"PC1", "PC2", "PC3", "PC4",
                                                 "PC5"};
    auto results = mockCommander::transform_covariate(cov_string);
    EXPECT_EQ(results.size(), expected_outputs.size());
    for (size_t i = 0; i < results.size(); ++i)
    { EXPECT_STREQ(expected_outputs[i].c_str(), results[i].c_str()); }
    invalid_cov_input("@PC[[1-5]]");
    invalid_cov_input("@PC[1-5");
    invalid_cov_input("@PC1-5]");
    invalid_cov_input("@PC[[1-5]");
    invalid_cov_input("@PC[1-5]]");
    invalid_cov_input("@PC[1-5][");
    invalid_cov_input("@PC[1-5,]");
    invalid_cov_input("@PC[,1-5]");
}
*/
void quick_check_unit(const std::string& input_str, const size_t exp_output,
                      const size_t def_power = 0, const bool memory = false)
{
    mockCommander commander;
    size_t value;
    commander.check_parse_unit_value(input_str, "", def_power, value, memory);
    ASSERT_EQ(exp_output, value);
}

TEST(PARSE_UNIT, VALIDITY)
{
    mockCommander commander;
    size_t value = 0;
    // Check if valid
    ASSERT_FALSE(commander.check_parse_unit_value("m", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("b", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("mb", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("hi", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("TB", "--mem", 0, value));
}
TEST(PARSE_UNIT, OUT_BOUND)
{
    mockCommander commander;
    size_t value = 0;
    // out of bound
    ASSERT_FALSE(commander.check_parse_unit_value("1", "", 7, value));
    ASSERT_FALSE(
        commander.check_parse_unit_value("1000000000tb", "", 1, value));
}
TEST(PARSE_UNIT, NEGATIVES)
{
    // Check for negative values
    mockCommander commander;
    size_t value = 0;
    ASSERT_FALSE(commander.check_parse_unit_value("-1", "", 1, value));
    ASSERT_FALSE(commander.check_parse_unit_value("-1tb", "", 1, value));
}
TEST(PARSE_UNIT, WITH_UNIT)
{
    // default value should be ignored when user provide a unit
    quick_check_unit("1b", 1, 0);
    quick_check_unit("1b", 1, 1);
    quick_check_unit("1b", 1, 3);
    quick_check_unit("1b", 1, 4);
    quick_check_unit("1b", 1, 5);
    quick_check_unit("1b", 1, 6);
}
TEST(PARSE_UNIT, DIFFERENT_UNIT)
{
    // check unit works as expected
    quick_check_unit("1k", 1000, 0);
    quick_check_unit("1kb", 1000, 0);
    quick_check_unit("1m", 1000000, 0);
    quick_check_unit("1mb", 1000000, 0);
    quick_check_unit("1g", 1000000000, 0);
    quick_check_unit("1gb", 1000000000, 0);
    quick_check_unit("1t", 1000000000000, 0);
    quick_check_unit("1tb", 1000000000000, 0);
}
TEST(PARSE_UNIT, NON_INTEGER)
{
    mockCommander commander;
    size_t value = 0;
    // Check non-integer scenarios
    quick_check_unit("1.5k", 1500, 0);
    quick_check_unit("1.004k", 1004, 0);
    ASSERT_FALSE(commander.check_parse_unit_value("1.5b", "", 1, value));
}
TEST(PARSE_UNIT, DEFAULT_VALUE)
{
    // check default value works
    quick_check_unit("1", 1, 0);
    quick_check_unit("1", 1000, 1);
    quick_check_unit("1", 1000000000, 3);
    quick_check_unit("1", 1000000000000, 4);
    quick_check_unit("1", 1000000000000000, 5);
    quick_check_unit("1", 1000000000000000000, 6);
}
#endif // COMMANDER_TEST_H
