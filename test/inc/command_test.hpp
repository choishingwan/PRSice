#ifndef COMMAND_TEST_HPP
#define COMMAND_TEST_HPP
#include "commander.hpp"
#include "reporter.hpp"
#include "gtest/gtest.h"
// Guess I will have to assume the data is always located at
// ../test/data
// This is currently causing me too much stress as for some reason, the command
// line parse in commander doesn't work as expected in the test case (but work
// as expected when the same input is provided in the command line). I am not
// sure what caused this difference
// i.e BASE_ONLY test case should in theory register the base file input and
// output not having the target file but that's not the case. Instead, we need
// to add a "dummy" string in front of "--base" for that to work as expected.
// But then when we try to include both --target and --base, the whole thing
// failed

TEST(COMMANDER, NOT_ENOUGH_ARGUEMENT)
{
    Reporter reporter;
    Commander commander;
    int argc = 1;
    char* argv[] = {"PRSice", NULL};
    try
    {
        commander.init(argc, argv, reporter);
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(COMMANDER, USAGE)
{
    Reporter reporter;
    Commander commander;
    int argc = 2;
    char* argv[] = {"PRSice", "-h", NULL};
    try
    {
        ASSERT_FALSE(commander.init(argc, argv, reporter));
    }
    catch (...)
    {
        FAIL();
    }
}


TEST(COMMANDER, BASE_ONLY)
{
    Reporter reporter;
    Commander commander;
    int argc = 3;
    char* argv[] = {"PRSice", "--base", "../test/data/TEST.assoc", NULL};
    try
    {
        ASSERT_FALSE(commander.init(argc, argv, reporter));
    }
    catch (...)
    {
        SUCCEED();
    }
}

#endif
