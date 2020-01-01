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
    Reporter reporter(std::string(path + "LOG"));
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
    Reporter reporter(std::string(path + "LOG"));
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

class CovariateTest : public Commander
{
public:
    static std::vector<std::string>
    transform_covariate(const std::string& cov_in)
    {
        return Commander::transform_covariate(cov_in);
    }
};

TEST(COVARITE_TRANSFORM, TRANSFORMATION)
{
    // should not do transformation when not start with @
    std::string cov_string = "PC1";
    std::string expected = cov_string;
    ASSERT_STREQ(
        expected.c_str(),
        CovariateTest::transform_covariate(cov_string).front().c_str());
    // same for empty string
    cov_string = expected = "";
    ASSERT_STREQ(
        expected.c_str(),
        CovariateTest::transform_covariate(cov_string).front().c_str());
    // should be fine if the @ is in middle of the string
    cov_string = expected = "PC1@Home";
    ASSERT_STREQ(
        expected.c_str(),
        CovariateTest::transform_covariate(cov_string).front().c_str());
    // when start with @ but not with any [], we will just remove the @
    cov_string = "@PC1";
    expected = "PC1";
    ASSERT_STREQ(
        expected.c_str(),
        CovariateTest::transform_covariate(cov_string).front().c_str());
    cov_string = "@PC[1-5]";
    // in this order
    std::vector<std::string> expected_outputs = {"PC1", "PC2", "PC3", "PC4",
                                                 "PC5"};
    std::vector<std::string> results =
        CovariateTest::transform_covariate(cov_string);
    EXPECT_EQ(results.size(), expected_outputs.size());
    for (size_t i = 0; i < results.size(); ++i)
    { EXPECT_STREQ(expected_outputs[i].c_str(), results[i].c_str()); }
}
#endif // COMMANDER_TEST_H
