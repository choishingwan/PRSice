#ifndef REGION_TEST_HPP
#define REGION_TEST_HPP
#include "region.hpp"
#include "reporter.hpp"
#include "gtest/gtest.h"
#include <fstream>
TEST(REGION, SINGLE_INIT)
{
    Reporter reporter;
    std::string range = "chr2:1234";
    try
    {
        Region region(range, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}
TEST(REGION, SINGLE_RANGE_INIT)
{
    Reporter reporter;
    std::string range = "chr2:1234-1357";
    try
    {
        Region region(range, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}
TEST(REGION, SINGLE_RANGE_WRONG)
{
    Reporter reporter;
    // start must be smaller than end
    std::string range = "chr2:12341-1357";
    try
    {
        Region region(range, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION, MULTI_RANGE_INIT)
{
    Reporter reporter;
    std::string range = "chr6:369-4321,chr2:1234-1357";
    try
    {
        Region region(range, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}
TEST(REGION, MULTI_MIX_INIT)
{
    Reporter reporter;
    std::string range = "chr6:369-4321,chr2:1234";
    try
    {
        Region region(range, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}
TEST(REGION, MULTI_MORE_MIX_INIT)
{
    Reporter reporter;
    std::string range = "chr6:369-4321,chr2:1234,chr1:312345-9437690";
    try
    {
        Region region(range, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}
TEST(REGION, WRONG_INPUT)
{

    Reporter reporter;
    try
    {
        Region region("chr1", reporter);
        // in this case, we will assume this is a bed file, but we can't read
        // it, so we will have 0 regions
        ASSERT_EQ(region.size(), 0);
    }
    catch (...)
    {
        FAIL();
    }
}
TEST(REGION, WRONG_RANGE_FORMAT)
{

    Reporter reporter;
    try
    {
        Region region("chr1:1-2-3", reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION, RANGE_PARSE_PROBLEM)
{

    Reporter reporter;
    try
    {
        Region region("chr1:1-,2", reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
class REGION_EX_STRING : public ::testing::Test
{

protected:
    Region region;
    void SetUp() override
    {
        std::string range = "chr2:1234";
        Reporter reporter;
        region = Region(range, reporter);
    }
};
TEST_F(REGION_EX_STRING, FOUND)
{
    ASSERT_TRUE(region.check_exclusion(2, 1234));
}
TEST_F(REGION_EX_STRING, WRONG_CHR)
{
    ASSERT_FALSE(region.check_exclusion(1, 1234));
}
TEST_F(REGION_EX_STRING, WRONG_BIGGER_CHR)
{
    ASSERT_FALSE(region.check_exclusion(3, 1234));
}
TEST_F(REGION_EX_STRING, BP_TOO_SMALL)
{
    ASSERT_FALSE(region.check_exclusion(2, 1233));
}
TEST_F(REGION_EX_STRING, BP_TOO_LARGE)
{
    // we expect the end to be non-inclusive, therefore chr1:11 should be
    // excluded
    ASSERT_FALSE(region.check_exclusion(2, 1235));
}
TEST_F(REGION_EX_STRING, SEQUENT_SEARCH)
{
    ASSERT_FALSE(region.check_exclusion(1, 1234));
    ASSERT_TRUE(region.check_exclusion(2, 1234));
    ASSERT_FALSE(region.check_exclusion(2, 1235));
}
TEST_F(REGION_EX_STRING, FINE_SEQUENT_SEARCH)
{
    ASSERT_FALSE(region.check_exclusion(1, 1234));
    ASSERT_FALSE(region.check_exclusion(2, 1233));
    ASSERT_TRUE(region.check_exclusion(2, 1234));
    ASSERT_FALSE(region.check_exclusion(2, 1235));
}
class REGION_EX_STRING_REGION : public ::testing::Test
{

protected:
    Region region;
    void SetUp() override
    {
        std::string range = "chr2:1234-1357";
        Reporter reporter;
        region = Region(range, reporter);
    }
};
TEST_F(REGION_EX_STRING_REGION, SEQ_STRING_WITHIN_RANGE)
{
    // the following should all work
    for (intptr_t i = 1234; i < 1357; ++i)
        ASSERT_TRUE(region.check_exclusion(2, i));
}
TEST_F(REGION_EX_STRING_REGION, SEQ_STRING_START)
{
    // the following should all work
    ASSERT_TRUE(region.check_exclusion(2, 1234));
}
TEST_F(REGION_EX_STRING_REGION, SEQ_STRING_END)
{
    // the following should all work
    ASSERT_TRUE(region.check_exclusion(2, 1356));
}
TEST_F(REGION_EX_STRING_REGION, SEQ_UP_BOUND)
{
    // the following should all work
    ASSERT_FALSE(region.check_exclusion(2, 1357));
}
class REGION_STRING_MIX : public ::testing::Test
{

protected:
    Region region;
    void SetUp() override
    {
        // we need to account for both range and single base input
        // also we want to check chromosome overrun (e.g.
        // Range: chr1:4601-5678,chr2:1357-2468, SNP input 1:5679, 2:134,2:1357)
        // and make sure the input is not in sorted order
        std::string range =
            "chr2:1234-1357,chr1:4601-5678,chr12:314,chr6:98741-102380";
        Reporter reporter;
        region = Region(range, reporter);
    }
};
TEST_F(REGION_STRING_MIX, SINGLE_CHECK_FIRST)
{
    ASSERT_TRUE(region.check_exclusion(2, 1234));
}
TEST_F(REGION_STRING_MIX, SINGLE_CHECK_SECOND)
{
    ASSERT_TRUE(region.check_exclusion(1, 4890));
}
TEST_F(REGION_STRING_MIX, SINGLE_CHECK_THIRD)
{
    ASSERT_TRUE(region.check_exclusion(12, 314));
}
TEST_F(REGION_STRING_MIX, SINGLE_CHECK_FORTH)
{
    ASSERT_TRUE(region.check_exclusion(6, 102379));
}
TEST_F(REGION_STRING_MIX, CHECK_SEQ)
{
    ASSERT_TRUE(region.check_exclusion(1, 4890));
    ASSERT_TRUE(region.check_exclusion(2, 1234));
    ASSERT_TRUE(region.check_exclusion(6, 102379));
    ASSERT_TRUE(region.check_exclusion(12, 314));
}
TEST_F(REGION_STRING_MIX, CHECK_SEQ_UNSORTED)
{
    ASSERT_TRUE(region.check_exclusion(2, 1234));
    // we don't expect the input of the exclusion list is sorted as the listed
    // input might have chromosome ordered in different sequence than we've
    // anticipated
    ASSERT_TRUE(region.check_exclusion(1, 4890));
    ASSERT_TRUE(region.check_exclusion(12, 314));
    ASSERT_TRUE(region.check_exclusion(6, 102379));
}
TEST_F(REGION_STRING_MIX, MIXED_FOUND)
{
    ASSERT_TRUE(region.check_exclusion(2, 1234));
    // we don't expect the input of the exclusion list is sorted as the listed
    // input might have chromosome ordered in different sequence than we've
    // anticipated
    ASSERT_FALSE(region.check_exclusion(1, 1));
    ASSERT_TRUE(region.check_exclusion(12, 314));
    ASSERT_FALSE(region.check_exclusion(6, 12432145));
}
TEST_F(REGION_STRING_MIX, RUN_OVER)
{
    // As check exclusion now uses binary search, running over shouldn't really
    // be a problem
    ASSERT_FALSE(region.check_exclusion(1, 5679));
    ASSERT_TRUE(region.check_exclusion(2, 1234));
}
class REGION_BED_MIN : public ::testing::Test
{

protected:
    Region region;
    void SetUp() override
    {
        // we need to account for both range and single base input
        // also we want to check chromosome overrun (e.g.
        // Range: chr1:4601-5678,chr2:1357-2468, SNP input 1:5679, 2:134,2:1357)
        // and make sure the input is not in sorted order
        std::string range =
            "chr2:1234-1357,chr1:4601-5678,chr12:314,chr6:98741-102380";
        Reporter reporter;
        region = Region(range, reporter);
    }
};

// bed file exclusion region performance (we don't allow multibed)
// bed file read
// multibed without name (should use the file name)
// multibed with name
// gtf read
// msigdb read
#endif // REGION_TEST_HPP
