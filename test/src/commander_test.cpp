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
    Reporter reporter(std::string(path + "LOG"), true);
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
    Reporter reporter(std::string(path + "LOG"), true);
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
TEST(COVARITE_TRANSFORM, TRANSFORMATION)
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

void quick_check_unit(const std::string& input_str, const size_t exp_output,
                      const size_t def_power = 0, const bool memory = false)
{
    mockCommander commander;
    size_t value;
    commander.check_parse_unit_value(input_str, "", def_power, value, memory);
    ASSERT_EQ(exp_output, value);
}

TEST(PARSE_UNIT, CHECK_UNIT)
{
    mockCommander commander;
    size_t value = 0;
    // Check if valid
    ASSERT_FALSE(commander.check_parse_unit_value("m", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("b", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("mb", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("hi", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("TB", "--mem", 0, value));
    // check default value works
    quick_check_unit("1", 1, 0);
    quick_check_unit("1", 1000, 1);
    quick_check_unit("1", 1000000000, 3);
    quick_check_unit("1", 1000000000000, 4);
    quick_check_unit("1", 1000000000000000, 5);
    quick_check_unit("1", 1000000000000000000, 6);
    // default value should be ignored when user provide a unit
    quick_check_unit("1b", 1, 0);
    quick_check_unit("1b", 1, 1);
    quick_check_unit("1b", 1, 3);
    quick_check_unit("1b", 1, 4);
    quick_check_unit("1b", 1, 5);
    quick_check_unit("1b", 1, 6);
    // out of bound
    ASSERT_FALSE(commander.check_parse_unit_value("1", "", 7, value));
}
#endif // COMMANDER_TEST_H
