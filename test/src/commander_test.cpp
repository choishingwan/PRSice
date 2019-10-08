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
#endif // COMMANDER_TEST_H
