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
    std::vector<size_t> base_col_index = commander.index();
    std::vector<bool> base_has_col = commander.has_col();
    ASSERT_EQ(base_col_index.size(), base_has_col.size());
    ASSERT_EQ(base_col_index.size(), +BASE_INDEX::MAX + 1);
    // everything else should be empty or set to default
    ASSERT_TRUE(commander.base_name().empty());
    ASSERT_FALSE(commander.all_scores());
    ASSERT_FALSE(commander.fastscore());
    ASSERT_FALSE(commander.no_clump());
    ASSERT_FALSE(commander.no_full());
    ASSERT_FALSE(commander.no_regress());
    ASSERT_TRUE(commander.background().empty());
    ASSERT_FALSE(commander.beta());
}

#endif // COMMANDER_TEST_H
