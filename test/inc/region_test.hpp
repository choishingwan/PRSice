#ifndef REGION_TEST_HPP
#define REGION_TEST_HPP
#include "genotype.hpp"
#include "plink_common.hpp"
#include "region.hpp"
#include "reporter.hpp"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <string>

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
        // it, so we will have throw an error
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
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
class REGION_BED_MIN_TAB_NO_OVER : public ::testing::Test
{

protected:
    Region region;
    void SetUp() override
    {
        // we need to account for both range and single base input
        // also we want to check chromosome overrun (e.g.
        // Range: chr1:4601-5678,chr2:1357-2468, SNP input 1:5679, 2:134,2:1357)
        // and make sure the input is not in sorted order
        std::ofstream bed_file;
        std::string bed_name = "Test.bed";
        bed_file.open(bed_name.c_str());
        if (!bed_file.is_open()) {
            throw std::runtime_error("Error: Cannot open bed file");
        }
        //  now generate the output required
        bed_file << "2\t19182\t32729\n"
                 << "2\t94644\t98555\n"
                 << "3\t3209\t18821\n"
                 << "3\t29863\t38285\n"
                 << "4\t20139\t97433\n"
                 << "5\t13998\t35076\n"
                 << "5\t50433\t97855\n"
                 << "6\t34611\t45099\n"
                 << "7\t7080\t45054\n"
                 << "10\t54504\t62968\n"
                 << "11\t20844\t26475\n"
                 << "12\t38890\t50405\n"
                 << "13\t56146\t67102\n"
                 << "14\t1694\t47285\n"
                 << "15\t4706\t10214\n"
                 << "15\t26926\t85344\n"
                 << "16\t12143\t36596\n"
                 << "16\t43942\t85160\n"
                 << "19\t22463\t39329\n"
                 << "19\t46559\t49131\n"
                 << "20\t64037\t98171\n"
                 << "21\t9363\t49431\n";
        bed_file.close();
        Reporter reporter("LOG");
        region = Region(bed_name, reporter);
    }
};
TEST_F(REGION_BED_MIN_TAB_NO_OVER, CHECK_INPUT_PARSING)
{
    // for exclusion set, we will only have one set
    try
    {
        ASSERT_EQ(region.num_bound(0), 22);
    }
    catch (...)
    {
        FAIL();
    }
    try
    {
        // and we will through error if we are out of bound
        region.num_bound(1);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST_F(REGION_BED_MIN_TAB_NO_OVER, CHECK_INCLUSION)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input
    EXPECT_FALSE(region.check_exclusion(2, 19181 + 1));
    EXPECT_TRUE(region.check_exclusion(2, 19182 + 1));
    EXPECT_TRUE(region.check_exclusion(2, 19183 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(2, 32728 + 1));
    EXPECT_FALSE(region.check_exclusion(2, 32729 + 1));
    EXPECT_FALSE(region.check_exclusion(2, 94643 + 1));
    EXPECT_TRUE(region.check_exclusion(2, 94644 + 1));
    EXPECT_TRUE(region.check_exclusion(2, 94645 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(2, 98554 + 1));
    EXPECT_FALSE(region.check_exclusion(2, 98555 + 1));
    EXPECT_FALSE(region.check_exclusion(3, 3208 + 1));
    EXPECT_TRUE(region.check_exclusion(3, 3209 + 1));
    EXPECT_TRUE(region.check_exclusion(3, 3210 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(3, 18820 + 1));
    EXPECT_FALSE(region.check_exclusion(3, 18821 + 1));
    EXPECT_FALSE(region.check_exclusion(13, 56145 + 1));
    EXPECT_TRUE(region.check_exclusion(13, 56146 + 1));
    EXPECT_TRUE(region.check_exclusion(13, 56147 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(13, 67101 + 1));
    EXPECT_FALSE(region.check_exclusion(13, 67102 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 9362 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9363 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9364 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(21, 49430 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 49431 + 1));
}
TEST_F(REGION_BED_MIN_TAB_NO_OVER, RANDOM_ORDER)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input
    EXPECT_FALSE(region.check_exclusion(13, 56145 + 1));
    EXPECT_TRUE(region.check_exclusion(2, 19183 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 49431 + 1));
    EXPECT_FALSE(region.check_exclusion(3, 3208 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(2, 32728 + 1));
    EXPECT_FALSE(region.check_exclusion(2, 32729 + 1));
    EXPECT_TRUE(region.check_exclusion(2, 94645 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(2, 19182 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9363 + 1));
    EXPECT_TRUE(region.check_exclusion(2, 98554 + 1));
    EXPECT_FALSE(region.check_exclusion(2, 98555 + 1));
    EXPECT_TRUE(region.check_exclusion(3, 3209 + 1));
    EXPECT_TRUE(region.check_exclusion(3, 3210 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(3, 18820 + 1));
    EXPECT_FALSE(region.check_exclusion(3, 18821 + 1));
    EXPECT_TRUE(region.check_exclusion(13, 56147 + 1));
    EXPECT_FALSE(region.check_exclusion(2, 19181 + 1));
    EXPECT_FALSE(region.check_exclusion(2, 94643 + 1));
    EXPECT_TRUE(region.check_exclusion(2, 94644 + 1));
    EXPECT_TRUE(region.check_exclusion(13, 56146 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(13, 67101 + 1));
    EXPECT_FALSE(region.check_exclusion(13, 67102 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 9362 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9364 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(region.check_exclusion(21, 49430 + 1));
}
class REGION_BED_MIN_TAB : public ::testing::Test
{

protected:
    Region region;
    void SetUp() override
    {
        // we need to account for both range and single base input
        // also we want to check chromosome overrun (e.g.
        // Range: chr1:4601-5678,chr2:1357-2468, SNP input 1:5679, 2:134,2:1357)
        // and make sure the input is not in sorted order
        std::ofstream bed_file;
        std::string bed_name = "Test.bed";
        bed_file.open(bed_name.c_str());
        if (!bed_file.is_open()) {
            throw std::runtime_error("Error: Cannot open bed file");
        }
        //  now generate the output required
        bed_file << "2\t19182\t32729\n"
                 << "2\t94644\t98555\n"
                 << "3\t3209\t18821\n"
                 << "3\t29863\t38285\n"
                 << "4\t20139\t97433\n"
                 << "5\t13998\t35076\n"
                 << "5\t50433\t97855\n"
                 << "6\t34611\t45099\n"
                 << "6\t45503\t49751\n"
                 << "7\t7080\t45054\n"
                 << "7\t30305\t45723\n" // overlap
                 << "10\t54504\t62968\n"
                 << "11\t20844\t26475\n"
                 << "12\t38890\t50405\n"
                 << "13\t56146\t67102\n"
                 << "14\t1694\t47285\n"
                 << "14\t5225\t78548\n"  // overlap
                 << "14\t13102\t45658\n" // overlap
                 << "15\t4706\t10214\n"
                 << "15\t26926\t85344\n"
                 << "15\t78969\t98716\n" // overlap
                 << "16\t7139\t73747\n"
                 << "16\t12143\t36596\n" // overlap
                 << "16\t31326\t56532\n" // overlap
                 << "16\t43942\t85160\n" // overlap
                 << "19\t22463\t39329\n"
                 << "19\t46559\t49131\n"
                 << "20\t64037\t98171\n"
                 << "21\t9363\t49431\n"
                 << "21\t43440\t82120\n"; // overlap
        bed_file.close();
        Reporter reporter("LOG");
        region = Region(bed_name, reporter);
    }
};
TEST_F(REGION_BED_MIN_TAB, CHECK_INPUT_PARSING)
{
    // for exclusion set, we will only have one set
    try
    {
        ASSERT_EQ(region.num_bound(0), 22);
    }
    catch (...)
    {
        FAIL();
    }
    try
    {
        // and we will through error if we are out of bound
        region.num_bound(1);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST_F(REGION_BED_MIN_TAB, CHECK_INCLUSION_OVERLAPPED)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input

    EXPECT_FALSE(region.check_exclusion(7, 7079 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 7080 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 7081 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45053 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45054 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45055 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30303 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30305 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30306 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45722 + 1));
    EXPECT_FALSE(region.check_exclusion(7, 45723 + 1));
    EXPECT_FALSE(region.check_exclusion(7, 45724 + 1));

    EXPECT_FALSE(region.check_exclusion(14, 1693 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 1694 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 1695 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47284 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47285 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47286 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5224 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5225 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5226 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13101 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13102 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13103 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45657 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45658 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45659 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 78547 + 1));
    EXPECT_FALSE(region.check_exclusion(14, 78548 + 1));
    EXPECT_FALSE(region.check_exclusion(14, 78549 + 1));

    EXPECT_FALSE(region.check_exclusion(21, 9362 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9363 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9364 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49430 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49431 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49432 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43439 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43440 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43441 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 82119 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 82120 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 82121 + 1));
}
class REGION_BED_MIN_SPACE : public ::testing::Test
{

protected:
    Region region;
    void SetUp() override
    {
        // we need to account for both range and single base input
        // also we want to check chromosome overrun (e.g.
        // Range: chr1:4601-5678,chr2:1357-2468, SNP input 1:5679, 2:134,2:1357)
        // and make sure the input is not in sorted order
        std::ofstream bed_file;
        std::string bed_name = "Test.bed";
        bed_file.open(bed_name.c_str());
        if (!bed_file.is_open()) {
            throw std::runtime_error("Error: Cannot open bed file");
        }
        //  now generate the output required
        bed_file << "2 19182 32729\n"
                 << "2 94644 98555\n"
                 << "3 3209 18821\n"
                 << "3 29863 38285\n"
                 << "4 20139 97433\n"
                 << "5 13998 35076\n"
                 << "5 50433 97855\n"
                 << "6 34611 45099\n"
                 << "6 45503 49751\n"
                 << "7 7080 45054\n"
                 << "7 30305 45723\n" // overlap
                 << "10 54504 62968\n"
                 << "11 20844 26475\n"
                 << "12 38890 50405\n"
                 << "13 56146 67102\n"
                 << "14 1694 47285\n"
                 << "14 5225 78548\n"  // overlap
                 << "14 13102 45658\n" // overlap
                 << "15 4706 10214\n"
                 << "15 26926 85344\n"
                 << "15 78969 98716\n" // overlap
                 << "16 7139 73747\n"
                 << "16 12143 36596\n" // overlap
                 << "16 31326 56532\n" // overlap
                 << "16 43942 85160\n" // overlap
                 << "19 22463 39329\n"
                 << "19 46559 49131\n"
                 << "20 64037 98171\n"
                 << "21 9363 49431\n"
                 << "21 43440 82120\n"; // overlap
        bed_file.close();
        Reporter reporter("LOG");
        region = Region(bed_name, reporter);
    }
};
TEST_F(REGION_BED_MIN_SPACE, CHECK_INPUT_PARSING)
{
    // for exclusion set, we will only have one set
    try
    {
        ASSERT_EQ(region.num_bound(0), 22);
    }
    catch (...)
    {
        FAIL();
    }
    try
    {
        // and we will through error if we are out of bound
        region.num_bound(1);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST_F(REGION_BED_MIN_SPACE, CHECK_INCLUSION_OVERLAPPED)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input

    EXPECT_FALSE(region.check_exclusion(7, 7079 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 7080 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 7081 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45053 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45054 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45055 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30303 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30305 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30306 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45722 + 1));
    EXPECT_FALSE(region.check_exclusion(7, 45723 + 1));
    EXPECT_FALSE(region.check_exclusion(7, 45724 + 1));

    EXPECT_FALSE(region.check_exclusion(14, 1693 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 1694 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 1695 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47284 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47285 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47286 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5224 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5225 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5226 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13101 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13102 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13103 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45657 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45658 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45659 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 78547 + 1));
    EXPECT_FALSE(region.check_exclusion(14, 78548 + 1));
    EXPECT_FALSE(region.check_exclusion(14, 78549 + 1));

    EXPECT_FALSE(region.check_exclusion(21, 9362 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9363 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9364 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49430 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49431 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49432 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43439 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43440 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43441 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 82119 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 82120 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 82121 + 1));
}
class REGION_BED_5_COLUMN : public ::testing::Test
{

protected:
    Region region;
    void SetUp() override
    {
        // not enough for stand yet
        std::ofstream bed_file;
        std::string bed_name = "Test.bed";
        bed_file.open(bed_name.c_str());
        if (!bed_file.is_open()) {
            throw std::runtime_error("Error: Cannot open bed file");
        }
        //  now generate the output required
        bed_file << "2 19182 32729 . .\n"
                 << "2 94644 98555 . .\n"
                 << "3 3209 18821 . .\n"
                 << "3 29863 38285 . .\n"
                 << "4 20139 97433 . .\n"
                 << "5 13998 35076 . .\n"
                 << "5 50433 97855 . .\n"
                 << "6 34611 45099 . .\n"
                 << "6 45503 49751 . .\n"
                 << "7 7080 45054 . .\n"
                 << "7 30305 45723 . .\n" // overlap
                 << "10 54504 62968 . .\n"
                 << "11 20844 26475 . .\n"
                 << "12 38890 50405 . .\n"
                 << "13 56146 67102 . .\n"
                 << "14 1694 47285 . .\n"
                 << "14 5225 78548 . .\n"  // overlap
                 << "14 13102 45658 . .\n" // overlap
                 << "15 4706 10214 . .\n"
                 << "15 26926 85344 . .\n"
                 << "15 78969 98716 . .\n" // overlap
                 << "16 7139 73747 . .\n"
                 << "16 12143 36596 . .\n" // overlap
                 << "16 31326 56532 . .\n" // overlap
                 << "16 43942 85160 . .\n" // overlap
                 << "19 22463 39329 . .\n"
                 << "19 46559 49131 . .\n"
                 << "20 64037 98171 . .\n"
                 << "21 9363 49431 . .\n"
                 << "21 43440 82120 . .\n"; // overlap
        bed_file.close();
        Reporter reporter("LOG");
        region = Region(bed_name, reporter);
    }
};
TEST_F(REGION_BED_5_COLUMN, CHECK_INPUT_PARSING)
{
    // for exclusion set, we will only have one set
    try
    {
        ASSERT_EQ(region.num_bound(0), 22);
    }
    catch (...)
    {
        FAIL();
    }
    try
    {
        // and we will through error if we are out of bound
        region.num_bound(1);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST_F(REGION_BED_5_COLUMN, CHECK_INCLUSION_OVERLAPPED)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input

    EXPECT_FALSE(region.check_exclusion(7, 7079 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 7080 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 7081 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45053 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45054 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45055 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30303 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30305 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30306 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45722 + 1));
    EXPECT_FALSE(region.check_exclusion(7, 45723 + 1));
    EXPECT_FALSE(region.check_exclusion(7, 45724 + 1));

    EXPECT_FALSE(region.check_exclusion(14, 1693 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 1694 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 1695 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47284 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47285 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47286 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5224 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5225 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5226 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13101 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13102 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13103 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45657 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45658 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45659 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 78547 + 1));
    EXPECT_FALSE(region.check_exclusion(14, 78548 + 1));
    EXPECT_FALSE(region.check_exclusion(14, 78549 + 1));

    EXPECT_FALSE(region.check_exclusion(21, 9362 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9363 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9364 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49430 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49431 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49432 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43439 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43440 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43441 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 82119 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 82120 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 82121 + 1));
}
class REGION_BED_5_WITH_STRAND : public ::testing::Test
{
    // For exclusion, strand information should not alter result (window padding
    // should all be 0)
protected:
    Region region;
    void SetUp() override
    {
        std::ofstream bed_file;
        std::string bed_name = "Test.bed";
        bed_file.open(bed_name.c_str());
        if (!bed_file.is_open()) {
            throw std::runtime_error("Error: Cannot open bed file");
        }
        //  now generate the output required
        bed_file << "2 19182 32729 . . .\n"
                 << "2 94644 98555 . . .\n"
                 << "3 3209 18821 . . .\n"
                 << "3 29863 38285 . . .\n"
                 << "4 20139 97433 . . .\n"
                 << "5 13998 35076 . . .\n"
                 << "5 50433 97855 . . .\n"
                 << "6 34611 45099 . . .\n"
                 << "6 45503 49751 . . .\n"
                 << "7 7080 45054 . . .\n"
                 << "7 30305 45723 . . .\n" // overlap
                 << "10 54504 62968 . . .\n"
                 << "11 20844 26475 . . .\n"
                 << "12 38890 50405 . . .\n"
                 << "13 56146 67102 . . .\n"
                 << "14 1694 47285 . . .\n"
                 << "14 5225 78548 . . .\n"  // overlap
                 << "14 13102 45658 . . .\n" // overlap
                 << "15 4706 10214 . . .\n"
                 << "15 26926 85344 . . .\n"
                 << "15 78969 98716 . . .\n" // overlap
                 << "16 7139 73747 . . .\n"
                 << "16 12143 36596 . . .\n" // overlap
                 << "16 31326 56532 . . .\n" // overlap
                 << "16 43942 85160 . . .\n" // overlap
                 << "19 22463 39329 . . .\n"
                 << "19 46559 49131 . . .\n"
                 << "20 64037 98171 . . .\n"
                 << "21 9363 49431 . . .\n"
                 << "21 43440 82120 . . .\n"; // overlap
        bed_file.close();
        Reporter reporter("LOG");
        region = Region(bed_name, reporter);
    }
};
TEST_F(REGION_BED_5_WITH_STRAND, CHECK_INPUT_PARSING)
{
    // for exclusion set, we will only have one set
    try
    {
        ASSERT_EQ(region.num_bound(0), 22);
    }
    catch (...)
    {
        FAIL();
    }
    try
    {
        // and we will through error if we are out of bound
        region.num_bound(1);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST_F(REGION_BED_5_WITH_STRAND, CHECK_INCLUSION_OVERLAPPED)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input

    EXPECT_FALSE(region.check_exclusion(7, 7079 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 7080 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 7081 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45053 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45054 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45055 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30303 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30305 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 30306 + 1));
    EXPECT_TRUE(region.check_exclusion(7, 45722 + 1));
    EXPECT_FALSE(region.check_exclusion(7, 45723 + 1));
    EXPECT_FALSE(region.check_exclusion(7, 45724 + 1));

    EXPECT_FALSE(region.check_exclusion(14, 1693 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 1694 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 1695 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47284 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47285 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 47286 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5224 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5225 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 5226 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13101 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13102 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 13103 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45657 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45658 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 45659 + 1));
    EXPECT_TRUE(region.check_exclusion(14, 78547 + 1));
    EXPECT_FALSE(region.check_exclusion(14, 78548 + 1));
    EXPECT_FALSE(region.check_exclusion(14, 78549 + 1));

    EXPECT_FALSE(region.check_exclusion(21, 9362 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9363 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 9364 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49430 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49431 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 49432 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43439 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43440 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 43441 + 1));
    EXPECT_TRUE(region.check_exclusion(21, 82119 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 82120 + 1));
    EXPECT_FALSE(region.check_exclusion(21, 82121 + 1));
}
// now test different type of malform BED file
TEST(REGION_MALFORM_BED, NOT_ENOUGH_COLUMN)
{
    std::ofstream bed_file;
    std::string bed_name = "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "2 19182\n"
             << "2 94644 \n"
             << "3 3209\n"
             << "21 43440\n"; // overlap
    bed_file.close();
    Reporter reporter("LOG");
    try
    {
        // we want to penalize any form of malformed input
        Region region(bed_name, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_MALFORM_BED, INCONSISTEN_COLUMN_STRAND)
{
    std::ofstream bed_file;
    std::string bed_name = "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "2 19182 123141 . . +\n"
             << "2 94644 123567 .  \n"
             << "3 3209 123141 . . .\n"
             << "21 43440 123141 . . +\n"; // overlap
    bed_file.close();
    Reporter reporter("LOG");
    try
    {
        // we want to penalize any form of malformed input
        Region region(bed_name, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_MALFORM_BED, NOT_FOUND)
{
    std::string bed_name = "404.bed";
    Reporter reporter("LOG");
    try
    {
        // we want to penalize any form of malformed input
        Region region(bed_name, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_MALFORM_BED, UNSUPPORTED_STRAND)
{
    std::ofstream bed_file;
    std::string bed_name = "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "2 19182 123141 . . +\n"
             << "2 94644 123567 . . L\n"
             << "3 3209 123141 . . .\n"
             << "21 43440 123141 . . +\n"; // overlap
    bed_file.close();
    Reporter reporter("LOG");
    try
    {
        // we want to penalize any form of malformed input
        Region region(bed_name, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_MALFORM_BED, NEGATIVE_COORDINATE)
{
    std::ofstream bed_file;
    std::string bed_name = "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "2 19182 123141 . . +\n"
             << "2 -94644 123567 . . +\n"
             << "3 3209 123141 . . .\n"
             << "21 43440 123141 . . +\n"; // overlap
    bed_file.close();
    Reporter reporter("LOG");
    try
    {
        // we want to penalize any form of malformed input
        Region region(bed_name, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_STD_BED_INPUT, NO_RUN)
{
    std::vector<uintptr_t> not_found = {0};
    std::vector<uintptr_t> found = {0};
    std::ofstream bed_file;
    std::string bed_name = "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "2 19182 32729 . . .\n"
             << "2 94644 98555 . . .\n"
             << "3 3209 18821 . . .\n"
             << "3 29863 38285 . . .\n"
             << "4 20139 97433 . . .\n"
             << "5 13998 35076 . . .\n"
             << "5 50433 97855 . . .\n"
             << "6 34611 45099 . . .\n"
             << "6 45503 49751 . . .\n"
             << "7 7080 45054 . . .\n"
             << "7 30305 45723 . . .\n" // overlap
             << "10 54504 62968 . . .\n"
             << "11 20844 26475 . . .\n"
             << "12 38890 50405 . . .\n"
             << "13 56146 67102 . . .\n"
             << "14 1694 47285 . . .\n"
             << "14 5225 78548 . . .\n"  // overlap
             << "14 13102 45658 . . .\n" // overlap
             << "15 4706 10214 . . .\n"
             << "15 26926 85344 . . .\n"
             << "15 78969 98716 . . .\n" // overlap
             << "16 7139 73747 . . .\n"
             << "16 12143 36596 . . .\n" // overlap
             << "16 31326 56532 . . .\n" // overlap
             << "16 43942 85160 . . .\n" // overlap
             << "19 22463 39329 . . .\n"
             << "19 46559 49131 . . .\n"
             << "20 64037 98171 . . .\n"
             << "21 9363 49431 . . .\n"
             << "21 43440 82120 . . .\n"; // overlap
    bed_file.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    SET_BIT(0, not_found.data());
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    ASSERT_EQ(region.size(), 1);
    std::vector<uintptr_t> index = {0};
    // this is a SNP found in the bed file, but as we have not generated the
    // region (we haven't use the bed file), this will always be considered as
    // not found
    region.update_flag(5, "", 50533 + 1, index);
    ASSERT_EQ(index, not_found);
}
class REGION_STD_BED : public ::testing::Test
{
    // For exclusion, strand information should not alter result (window
    // padding should all be 0)
protected:
    Region region;
    std::vector<uintptr_t> not_found = {0};
    std::vector<uintptr_t> found = {0};
    void SetUp() override
    {
        std::ofstream bed_file;
        std::string bed_name = "Test.bed";
        bed_file.open(bed_name.c_str());
        if (!bed_file.is_open()) {
            throw std::runtime_error("Error: Cannot open bed file");
        }
        //  now generate the output required
        bed_file << "2 19182 32729 . . .\n"
                 << "2 94644 98555 . . .\n"
                 << "3 3209 18821 . . .\n"
                 << "3 29863 38285 . . .\n"
                 << "4 20139 97433 . . .\n"
                 << "5 13998 35076 . . .\n"
                 << "5 50433 97855 . . .\n"
                 << "6 34611 45099 . . .\n"
                 << "6 45503 49751 . . .\n"
                 << "7 7080 45054 . . .\n"
                 << "7 30305 45723 . . .\n" // overlap
                 << "10 54504 62968 . . .\n"
                 << "11 20844 26475 . . .\n"
                 << "12 38890 50405 . . .\n"
                 << "13 56146 67102 . . .\n"
                 << "14 1694 47285 . . .\n"
                 << "14 5225 78548 . . .\n"  // overlap
                 << "14 13102 45658 . . .\n" // overlap
                 << "15 4706 10214 . . .\n"
                 << "15 26926 85344 . . .\n"
                 << "15 78969 98716 . . .\n" // overlap
                 << "16 7139 73747 . . .\n"
                 << "16 12143 36596 . . .\n" // overlap
                 << "16 31326 56532 . . .\n" // overlap
                 << "16 43942 85160 . . .\n" // overlap
                 << "19 22463 39329 . . .\n"
                 << "19 46559 49131 . . .\n"
                 << "20 64037 98171 . . .\n"
                 << "21 9363 49431 . . .\n"
                 << "21 43440 82120 . . .\n"; // overlap
        bed_file.close();
        Reporter reporter("LOG");
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        region = Region(feature, 0, 0, false, false);
        std::vector<std::string> bed_names = {bed_name};
        Genotype dummy;
        region.generate_regions("", "", bed_names, "", "", "", dummy, reporter);
        SET_BIT(0, not_found.data());
        SET_BIT(0, found.data());
        SET_BIT(1, found.data());
    }
};
TEST_F(REGION_STD_BED, CHECK_INPUT_PARSING)
{
    // for exclusion set, we will only have one set
    try
    {
        ASSERT_EQ(region.num_bound(0), 1);
    }
    catch (...)
    {
        FAIL();
    }
    try
    {
        // and we will through error if we are out of bound
        ASSERT_EQ(region.num_bound(1), 22);
    }
    catch (...)
    {
        FAIL();
    }
    try
    {
        // and we will through error if we are out of bound
        region.num_bound(2);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST_F(REGION_STD_BED, CHECK_INCLUSION_OVERLAPPED)
{
    // with standard input, we can no longer use check_exclusion function as
    // that always uses the base region, which doesn't contain any boundary
    // instead, we must use the update flag function
    std::vector<uintptr_t> input = {0};
    input.front() = 0;
    region.update_flag(7, "", 7079 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 7080 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(7, "", 7081 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(7, "", 45053 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(7, "", 45054 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(7, "", 45055 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(7, "", 30303 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(7, "", 30305 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(7, "", 30306 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(7, "", 45722 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(7, "", 45723 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 45724 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(14, "", 1693 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(14, "", 1695 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 47284 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 47285 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 47286 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    // normally, unordered input will not work. But here, it work, because
    // we will not iterate to the next bound unless we have passed the
    // current bound. As the previous check and the current check falls
    // within the same bound, we should be able to get true for inclusion
    region.update_flag(14, "", 5224 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 5225 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 5226 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 13101 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 13102 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 13103 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 45657 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 45658 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 45659 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 78547 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 78548 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(14, "", 78549 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
}
TEST_F(REGION_STD_BED, UNORDERED_INCLUSION)
{
    // When the input isn't sorted. We will encounter false negative
    std::vector<uintptr_t> input = {0};
    input.front() = 0;
    region.update_flag(14, "", 1693 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(14, "", 1695 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 47284 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 47285 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 47286 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    // normally, unordered input will not work. But here, it work, because
    // we will not iterate to the next bound unless we have passed the
    // current bound. As the previous check and the current check falls
    // within the same bound, we should be able to get true for inclusion
    region.update_flag(14, "", 5224 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 5225 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 5226 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 13101 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 13102 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 13103 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 45657 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 45658 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 45659 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    // even though some of the below SNPs are within the bed file, they will
    // always return false as we have already moved onto chromosome 14
    input.front() = 0;
    region.update_flag(7, "", 7079 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 7080 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 7081 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 45053 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 45054 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 45055 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 30303 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 30305 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 30306 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 45722 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 45723 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(7, "", 45724 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    // but we should continue to be able to identify the chr14 findings as
    // we should not move on to the next bound when we encounter a chr
    // smaller than the current one
    input.front() = 0;
    region.update_flag(14, "", 78547 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), found.front());
    input.front() = 0;
    region.update_flag(14, "", 78548 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(14, "", 78549 + 1, input);
    // the flag should be set as 110 if not found
    EXPECT_EQ(input.front(), not_found.front());
}
TEST_F(REGION_STD_BED, RUN_OVER)
{
    // when a SNP is bigger than any region within the same chromosome, we
    // should move onto the first region on the next chromosome
    std::vector<uintptr_t> input = {0};
    input.front() = 0;
    region.update_flag(19, "", 49131 + 1, input);
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(20, "", 64037 + 1, input);
    EXPECT_EQ(input.front(), found.front());
}
TEST_F(REGION_STD_BED, MID_NOT_FOUND)
{
    std::vector<uintptr_t> input = {0};
    input.front() = 0;
    region.update_flag(19, "", 39329 + 1, input);
    EXPECT_EQ(input.front(), not_found.front());
    input.front() = 0;
    region.update_flag(20, "", 64037 + 1, input);
    EXPECT_EQ(input.front(), found.front());
}
class REGION_STD_BED_PAD : public ::testing::Test
{
    // For exclusion, strand information should not alter result (window
    // padding should all be 0)
protected:
    Region region;
    std::vector<uintptr_t> not_found = {0};
    std::vector<uintptr_t> found = {0};
    void SetUp() override
    {
        std::ofstream bed_file;
        std::string bed_name = "Test.bed";
        bed_file.open(bed_name.c_str());
        if (!bed_file.is_open()) {
            throw std::runtime_error("Error: Cannot open bed file");
        }
        //  now generate the output required
        bed_file << "2 19182 32729 . . .\n"
                 << "2 94644 98555 . . .\n"
                 << "3 3209 18821 . . .\n"
                 << "3 29863 38285 . . .\n"
                 << "4 20139 97433 . . +\n"
                 << "5 13998 35076 . . .\n"
                 << "5 50433 97855 . . .\n"
                 << "6 34611 45099 . . -\n"
                 << "6 45503 49751 . . .\n"
                 << "7 7080 45054 . . +\n"
                 << "7 30305 45723 . . .\n" // overlap
                 << "10 54504 62968 . . .\n"
                 << "11 20844 26475 . . .\n"
                 << "12 38890 50405 . . .\n"
                 << "13 56146 67102 . . .\n"
                 << "14 1694 47285 . . .\n"
                 << "14 5225 78548 . . .\n"  // overlap
                 << "14 13102 45658 . . -\n" // overlap
                 << "15 4706 10214 . . .\n"
                 << "15 26926 85344 . . .\n"
                 << "15 78969 98716 . . .\n" // overlap
                 << "16 7139 73747 . . .\n"
                 << "16 12143 36596 . . .\n" // overlap
                 << "16 31326 56532 . . .\n" // overlap
                 << "16 43942 85160 . . .\n" // overlap
                 << "19 22463 39329 . . .\n"
                 << "19 46559 49131 . . .\n"
                 << "20 64037 98171 . . .\n"
                 << "21 9363 49431 . . .\n"
                 << "21 43440 82120 . . .\n"; // overlap
        bed_file.close();
        Reporter reporter("LOG");
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        region = Region(feature, 10, 20, false, false);
        std::vector<std::string> bed_names = {bed_name};
        Genotype dummy;
        region.generate_regions("", "", bed_names, "", "", "", dummy, reporter);
        SET_BIT(0, not_found.data());
        SET_BIT(0, found.data());
        SET_BIT(1, found.data());
    }
};
TEST_F(REGION_STD_BED_PAD, CHECK_PAD)
{
    // normally, with standard input, we need to use the update_flag option
    // to check inclusion, but here we only have one set

    // We will see how the padding change the inclusion criteria
    std::vector<uintptr_t> index = {0};
    // this SNP doesn't contain the strand info, we should assume the start
    // is the 5' end
    index.front() = 0;
    region.update_flag(3, "", 29863 + 1 - 11, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(3, "", 29863 + 1 - 10, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(3, "", 29863 + 1, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(3, "", 38285 + 1, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(3, "", 38285 + 1 + 19, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(3, "", 38285 + 1 + 20, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), not_found.front());

    index.front() = 0;
    region.update_flag(4, "", 20139 + 1 - 11, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(4, "", 20139 + 1 - 10, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(4, "", 20139 + 1, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(4, "", 97433 + 1, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(4, "", 97433 + 1 + 19, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(4, "", 97433 + 1 + 20, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), not_found.front());

    // negative strand
    index.front() = 0;
    region.update_flag(6, "", 34611 + 1 - 21, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(6, "", 34611 + 1 - 20, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(6, "", 34611 + 1, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(6, "", 45099 + 1, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(6, "", 45099 + 1 + 9, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(6, "", 45099 + 1 + 10, index);
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(index.front(), not_found.front());
}
TEST(REGION_MULTI_BED, CHECK_NAME)
{
    Reporter reporter("LOG");
    std::ofstream bed_file;
    std::string bed_name = "Test.bed";
    std::string second_bed_name = "Test2.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "2 19182 32729 . . .\n"
             << "2 94644 98555 . . .\n";
    bed_file.close();
    bed_file.open(second_bed_name.c_str());
    bed_file << "2 19182 32729 . . .\n"
             << "2 94644 98555 . . .\n";
    bed_file.close();
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 10, 20, false, false);
    std::vector<std::string> bed_names = {std::string(bed_name + ":Name"),
                                          second_bed_name};
    Genotype dummy;
    region.generate_regions("", "", bed_names, "", "", "", dummy, reporter);
    ASSERT_STREQ(region.get_name(0).c_str(), "Base");
    ASSERT_STREQ(region.get_name(1).c_str(), "Name");
    ASSERT_STREQ(region.get_name(2).c_str(), "Test2.bed");
}
TEST(REGION_MULTI_BED, CHECK_NAME2)
{
    Reporter reporter("LOG");
    std::ofstream bed_file;
    std::string bed_name = "Test.bed";
    std::string second_bed_name = "Test2.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "2 19182 32729 . . .\n"
             << "2 94644 98555 . . .\n";
    bed_file.close();
    bed_file.open(second_bed_name.c_str());
    bed_file << "2 19182 32729 . . .\n"
             << "2 94644 98555 . . .\n";
    bed_file.close();
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 10, 20, false, false);
    std::vector<std::string> bed_names = {
        bed_name, std::string(second_bed_name + ":Name"),
    };
    Genotype dummy;
    region.generate_regions("", "", bed_names, "", "", "", dummy, reporter);
    ASSERT_STREQ(region.get_name(0).c_str(), "Base");
    ASSERT_STREQ(region.get_name(1).c_str(), "Test.bed");
    ASSERT_STREQ(region.get_name(2).c_str(), "Name");
}
// gtf read
// msigdb read
// problem is, with the current design of region class, we can't test gtf file
// and msigdb file separately and we need a mock to Genotype such that it can
// provide the fake max_chr
class GenotypeTest : public Genotype
{
public:
    GenotypeTest() { m_max_code = 22; }
};
// Any error in the GTF file will lead to throw
TEST(REGION_GTF_BASIC, NOT_EXIST)
{
    std::string gtf_name = "Test.gtf";
    remove(gtf_name.c_str());
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    Genotype dummy = GenotypeTest();
    try
    {
        region.generate_regions(gtf_name, "", bed_names, "", "", "", dummy,
                                reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, EMPTY)
{
    std::string gtf_name = "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    GenotypeTest dummy;
    try
    {
        region.generate_regions(gtf_name, "", bed_names, "", "", "", dummy,
                                reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, ALL_REGION_REMOVE)
{
    std::string gtf_name = "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "#!genome - version GRCh38\n"
           "#!genome - date 2013 - 12\n"
           "#!genome - build - accession NCBI : GCA_000001405 .22\n"
           "#!genebuild - last - updated 2016 - 06\n"
           "1\thavana\tstop_codon\t11869\t14409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"
           "1\thavana\tfive_prime_utr\t11869\t14409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n";
    gtf.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    GenotypeTest dummy;
    try
    {
        region.generate_regions(gtf_name, "", bed_names, "", "", "", dummy,
                                reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, MALFORMAT_SPACE)
{
    std::string gtf_name = "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "#!genome - version GRCh38\n"
           "#!genome - date 2013 - 12\n"
           "#!genome - build - accession NCBI : GCA_000001405 .22\n"
           "#!genebuild - last - updated 2016 - 06\n"
           "1 havana gene 11869 14409 . + . gene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n";
    gtf.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    Genotype dummy = GenotypeTest();
    try
    {
        region.generate_regions(gtf_name, "", bed_names, "", "", "", dummy,
                                reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, NEGATIVE_COORDINATE)
{
    std::string gtf_name = "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "#!genome - version GRCh38\n"
           "#!genome - date 2013 - 12\n"
           "#!genome - build - accession NCBI : GCA_000001405 .22\n"
           "#!genebuild - last - updated 2016 - 06\n"
           "1\thavana\tgene\t-11869\t14409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"
           "1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n";
    gtf.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    Genotype dummy = GenotypeTest();
    try
    {
        region.generate_regions(gtf_name, "", bed_names, "", "", "", dummy,
                                reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, BIGGER_START)
{
    std::string gtf_name = "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "#!genome - version GRCh38\n"
           "#!genome - date 2013 - 12\n"
           "#!genome - build - accession NCBI : GCA_000001405 .22\n"
           "#!genebuild - last - updated 2016 - 06\n"
           "1\thavana\tgene\t14409\t11869\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"
           "1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n";
    gtf.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    Genotype dummy = GenotypeTest();
    try
    {
        region.generate_regions(gtf_name, "", bed_names, "", "", "", dummy,
                                reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, UNDEFINED_STRAND)
{

    std::string gtf_name = "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "#!genome - version GRCh38\n"
           "#!genome - date 2013 - 12\n"
           "#!genome - build - accession NCBI : GCA_000001405 .22\n"
           "#!genebuild - last - updated 2016 - 06\n"
           "1\thavana\tgene\t11869\t14409\t.\t@\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"
           "1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n";
    gtf.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    Genotype dummy = GenotypeTest();
    try
    {
        region.generate_regions(gtf_name, "", bed_names, "", "", "", dummy,
                                reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, TAB_ATTRIBUTE)
{

    std::string gtf_name = "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "#!genome - version GRCh38\n"
           "#!genome - date 2013 - 12\n"
           "#!genome - build - accession NCBI : GCA_000001405 .22\n"
           "#!genebuild - last - updated 2016 - 06\n"
           "1\thavana\tgene\t11869\t14409\t.\t@\t.\tgene_id\t"
           "\"ENSG00000223972\";\t"
           "gene_version\t\"5\";\tgene_name\t\"DDX11L1\";\tgene_source\t"
           "\"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\";\thavana_gene_version \"2\";\n"
           "1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n";
    gtf.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    Genotype dummy = GenotypeTest();
    try
    {
        region.generate_regions(gtf_name, "", bed_names, "", "", "", dummy,
                                reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, NO_GENE_ID)
{

    std::string gtf_name = "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "#!genome - version GRCh38\n"
           "#!genome - date 2013 - 12\n"
           "#!genome - build - accession NCBI : GCA_000001405 .22\n"
           "#!genebuild - last - updated 2016 - 06\n"
           "1\thavana\tgene\t11869\t14409\t.\t@\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"
           "1\thavana\tgene\t11869\t14409\t.\t+\t.\t"
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n";
    gtf.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    Genotype dummy = GenotypeTest();
    try
    {
        region.generate_regions(gtf_name, "", bed_names, "", "", "", dummy,
                                reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
class REGION_GTF_FEATURE : public ::testing::Test
{
    // For exclusion, strand information should not alter result (window
    // padding should all be 0)
protected:
    Region region;
    std::vector<uintptr_t> not_found = {0};
    void SetUp() override
    {
        std::string gtf_name = "Test.gtf";
        std::string gmt_name = "Test.gmt";
        std::ofstream gtf, gmt;
        gtf.open(gtf_name.c_str());
        gmt.open(gmt_name.c_str());
        gtf << "#!genome-build GRCh38.p7\n"
               "#!genome - version GRCh38\n"
               "#!genome - date 2013 - 12\n"
               "#!genome - build - accession NCBI : GCA_000001405 .22\n"
               "#!genebuild - last - updated 2016 - 06\n"
               "1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "
               "\"ENSG00000223972\"; "
               "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source "
               "\"havana\"; "
               "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
               "havana_gene "
               "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

               "1\thavana\tfive_prime_utr\t15869\t16409\t.\t+\t.\tgene_id "
               "\"ENSG00000223973\"; "
               "gene_version \"5\"; gene_name \"DDX11L2\"; gene_source "
               "\"havana\"; "
               "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
               "havana_gene "
               "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

               "12\thavana\ttranscript\t11399381\t11486678\t.\t-\t.\tgene_id "
               "\"ENSG00000255790\"; gene_version \"5\"; transcript_id "
               "\"ENST00000538349\"; transcript_version \"5\"; gene_name "
               "\"RP11-711K1.7\"; gene_source \"havana\"; gene_biotype "
               "\"sense_intronic\"; havana_gene \"OTTHUMG00000169117\"; "
               "havana_gene_version \"2\"; transcript_name "
               "\"RP11-711K1.7-001\"; transcript_source \"havana\"; "
               "transcript_biotype \"sense_intronic\"; havana_transcript "
               "\"OTTHUMT00000402313\"; havana_transcript_version \"2\"; tag "
               "\"basic\"; transcript_support_level \"4\";\n"

               "12\tensembl_havana\tCDS\t119697659\t119697838\t.\t-\t1\tgene_"
               "id \"ENSG00000122966\"; gene_version \"14\"; transcript_id "
               "\"ENST00000261833\"; transcript_version \"11\"; exon_number "
               "\"45\"; gene_name \"CIT\"; gene_source \"ensembl_havana\"; "
               "gene_biotype \"protein_coding\"; havana_gene "
               "\"OTTHUMG00000134325\"; havana_gene_version \"8\"; "
               "transcript_name \"CIT-001\"; transcript_source "
               "\"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag "
               "\"CCDS\"; ccds_id \"CCDS9192\"; havana_transcript "
               "\"OTTHUMT00000259410\"; havana_transcript_version \"4\"; "
               "protein_id \"ENSP00000261833\"; protein_version \"7\"; tag "
               "\"basic\"; transcript_support_level \"1\";\n"

               "15\tensembl\tCDS\t55320275\t55320410\t.\t+\t2\tgene_id "
               "\"ENSG00000069943\"; gene_version \"9\"; transcript_id "
               "\"ENST00000539642\"; transcript_version \"5\"; exon_number "
               "\"2\"; gene_name \"PIGB\"; gene_source \"ensembl_havana\"; "
               "gene_biotype \"protein_coding\"; havana_gene "
               "\"OTTHUMG00000172654\"; havana_gene_version \"1\"; "
               "transcript_name \"PIGB-201\"; transcript_source \"ensembl\"; "
               "transcript_biotype \"protein_coding\"; protein_id "
               "\"ENSP00000438963\"; protein_version \"2\"; tag \"basic\"; "
               "transcript_support_level \"5\";\n";
        gtf.close();
        gmt << "SET1 ENSG00000223972" << std::endl;
        gmt << "SET2 ENSG00000223973" << std::endl; // should be filtered
        gmt << "SET3 ENSG00000255790" << std::endl; // should also be filtered
        gmt << "SET4 CIT" << std::endl;
        gmt << "SET5 www.google.com ENSG00000069943" << std::endl;
        gmt << "SET6 www.google.com ENSG00000223972 ENSG00000255790 "
               "ENSG00000223973 ENSG00000255790 ENSG00000122966 "
            << std::endl;
        gmt.close();
        Reporter reporter("LOG");
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        region = Region(feature, 0, 0, false, false);
        std::vector<std::string> bed_names = {};
        Genotype dummy = GenotypeTest();
        region.generate_regions(gtf_name, gmt_name, bed_names, "", "", "",
                                dummy, reporter);
        SET_BIT(0, not_found.data());
    }
};
TEST_F(REGION_GTF_FEATURE, FEATURE_FILTER)
{
    // we don't "remove" from this stage. Benefit = simplier, and also better
    // capturing of duplicated gene set name?
    ASSERT_EQ(region.size(), 7);
}
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET1)
{
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(6, found.data());
    // 1 havana gene 11869 14409
    region.update_flag(1, "", 11868, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 11869, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 11870, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 14408, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 14409, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 14410, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET2)
{
    // should all be failed
    std::vector<uintptr_t> index = {0};
    region.update_flag(1, "", 15868, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 15869, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 15870, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 16408, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 16409, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 16410, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET3)
{
    // should all be failed
    std::vector<uintptr_t> index = {0};
    region.update_flag(12, "", 11399380, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11399381, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11399382, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11486677, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11486678, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11486679, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET4)
{
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, found.data());
    SET_BIT(4, found.data());
    SET_BIT(6, found.data());
    region.update_flag(12, "", 119697658, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697659, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697660, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697837, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697838, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697839, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET5)
{
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, found.data());
    SET_BIT(5, found.data());
    region.update_flag(15, "", 55320274, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320275, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320276, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320409, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320410, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320411, index);
    ASSERT_EQ(index.front(), not_found.front());
}
class REGION_GTF_PAD : public ::testing::Test
{
    // For exclusion, strand information should not alter result (window
    // padding should all be 0)
protected:
    Region region;
    std::vector<uintptr_t> not_found = {0};
    void SetUp() override
    {
        std::string gtf_name = "Test.gtf";
        std::string gmt_name = "Test.gmt";
        std::ofstream gtf, gmt;
        gtf.open(gtf_name.c_str());
        gmt.open(gmt_name.c_str());
        gtf << "#!genome-build GRCh38.p7\n"
               "#!genome - version GRCh38\n"
               "#!genome - date 2013 - 12\n"
               "#!genome - build - accession NCBI : GCA_000001405 .22\n"
               "#!genebuild - last - updated 2016 - 06\n"
               "1\thavana\tgene\t11869\t14409\t.\t.\t.\tgene_id "
               "\"ENSG00000223972\"; "
               "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source "
               "\"havana\"; "
               "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
               "havana_gene "
               "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

               "1\thavana\tfive_prime_utr\t15869\t16409\t.\t+\t.\tgene_id "
               "\"ENSG00000223973\"; "
               "gene_version \"5\"; gene_name \"DDX11L2\"; gene_source "
               "\"havana\"; "
               "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
               "havana_gene "
               "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

               "12\thavana\ttranscript\t11399381\t11486678\t.\t-\t.\tgene_id "
               "\"ENSG00000255790\"; gene_version \"5\"; transcript_id "
               "\"ENST00000538349\"; transcript_version \"5\"; gene_name "
               "\"RP11-711K1.7\"; gene_source \"havana\"; gene_biotype "
               "\"sense_intronic\"; havana_gene \"OTTHUMG00000169117\"; "
               "havana_gene_version \"2\"; transcript_name "
               "\"RP11-711K1.7-001\"; transcript_source \"havana\"; "
               "transcript_biotype \"sense_intronic\"; havana_transcript "
               "\"OTTHUMT00000402313\"; havana_transcript_version \"2\"; tag "
               "\"basic\"; transcript_support_level \"4\";\n"

               "12\tensembl_havana\tCDS\t119697659\t119697838\t.\t-\t1\tgene_"
               "id \"ENSG00000122966\"; gene_version \"14\"; transcript_id "
               "\"ENST00000261833\"; transcript_version \"11\"; exon_number "
               "\"45\"; gene_name \"CIT\"; gene_source \"ensembl_havana\"; "
               "gene_biotype \"protein_coding\"; havana_gene "
               "\"OTTHUMG00000134325\"; havana_gene_version \"8\"; "
               "transcript_name \"CIT-001\"; transcript_source "
               "\"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag "
               "\"CCDS\"; ccds_id \"CCDS9192\"; havana_transcript "
               "\"OTTHUMT00000259410\"; havana_transcript_version \"4\"; "
               "protein_id \"ENSP00000261833\"; protein_version \"7\"; tag "
               "\"basic\"; transcript_support_level \"1\";\n"

               "15\tensembl\tCDS\t55320275\t55320410\t.\t+\t2\tgene_id "
               "\"ENSG00000069943\"; gene_version \"9\"; transcript_id "
               "\"ENST00000539642\"; transcript_version \"5\"; exon_number "
               "\"2\"; gene_name \"PIGB\"; gene_source \"ensembl_havana\"; "
               "gene_biotype \"protein_coding\"; havana_gene "
               "\"OTTHUMG00000172654\"; havana_gene_version \"1\"; "
               "transcript_name \"PIGB-201\"; transcript_source \"ensembl\"; "
               "transcript_biotype \"protein_coding\"; protein_id "
               "\"ENSP00000438963\"; protein_version \"2\"; tag \"basic\"; "
               "transcript_support_level \"5\";\n";
        gtf.close();
        gmt << "SET1 ENSG00000223972" << std::endl;
        gmt << "SET2 ENSG00000223973" << std::endl; // should be filtered
        gmt << "SET3 ENSG00000255790" << std::endl; // should also be filtered
        gmt << "SET4 CIT" << std::endl;
        gmt << "SET5 www.google.com ENSG00000069943" << std::endl;
        gmt << "SET6 www.google.com ENSG00000223972 ENSG00000255790 "
               "ENSG00000223973 ENSG00000255790 ENSG00000122966 "
            << std::endl;
        gmt.close();
        Reporter reporter("LOG");
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        region = Region(feature, 10, 20, false, false);
        std::vector<std::string> bed_names = {};
        Genotype dummy = GenotypeTest();
        region.generate_regions(gtf_name, gmt_name, bed_names, "", "", "",
                                dummy, reporter);
        SET_BIT(0, not_found.data());
    }
};
TEST_F(REGION_GTF_PAD, FEATURE_FILTER)
{
    // we don't "remove" from this stage. Benefit = simplier, and also better
    // capturing of duplicated gene set name?
    ASSERT_EQ(region.size(), 7);
}
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET1)
{
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(6, found.data());
    // 1 havana gene 11869 14409
    region.update_flag(1, "", 11858, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 11859, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 11860, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 14428, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 14429, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 14430, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET2)
{
    // should all be failed
    std::vector<uintptr_t> index = {0};
    region.update_flag(1, "", 15858, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 15859, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 15860, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 16428, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 16429, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 16430, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET3)
{
    // should all be failed
    std::vector<uintptr_t> index = {0};
    region.update_flag(12, "", 11399360, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11399361, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11399362, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11486687, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11486688, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11486689, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET4)
{
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, found.data());
    SET_BIT(4, found.data());
    SET_BIT(6, found.data());
    // 1 havana gene 11869 14409
    region.update_flag(12, "", 119697638, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697639, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697640, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697847, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697848, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697849, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET5)
{
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, found.data());
    SET_BIT(5, found.data());
    // 1 havana gene 11869 14409
    region.update_flag(15, "", 55320264, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320265, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320266, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320429, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320430, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(15, "", 55320431, index);
    ASSERT_EQ(index.front(), not_found.front());
}
// need class for following
class REGION_GTF_MULTI_EX : public ::testing::Test
{
    // For exclusion, strand information should not alter result (window
    // padding should all be 0)
protected:
    Region region;
    std::vector<uintptr_t> not_found = {0};
    void SetUp() override
    {
        std::string gtf_name = "Test.gtf";
        std::string gmt_name = "Test.gmt";
        std::ofstream gtf, gmt;
        gtf.open(gtf_name.c_str());
        gmt.open(gmt_name.c_str());
        gtf << "#!genome-build GRCh38.p7\n"
               "#!genome - version GRCh38\n"
               "#!genome - date 2013 - 12\n"
               "#!genome - build - accession NCBI : GCA_000001405 .22\n"
               "#!genebuild - last - updated 2016 - 06\n"
               "1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "
               "\"ENSG00000223972\"; "
               "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source "
               "\"havana\"; "
               "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
               "havana_gene "
               "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

               "1\thavana\tgene\t15869\t16409\t.\t+\t.\tgene_id "
               "\"ENSG00000223973\"; "
               "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source "
               "\"havana\"; "
               "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
               "havana_gene "
               "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

               "12\thavana\tCDS\t11399381\t11486678\t.\t-\t.\tgene_id "
               "\"ENSG00000122966\"; gene_version \"5\"; transcript_id "
               "\"ENST00000538349\"; transcript_version \"5\"; gene_name "
               "\"CIT\"; gene_source \"havana\"; gene_biotype "
               "\"sense_intronic\"; havana_gene \"OTTHUMG00000169117\"; "
               "havana_gene_version \"2\"; transcript_name "
               "\"RP11-711K1.7-001\"; transcript_source \"havana\"; "
               "transcript_biotype \"sense_intronic\"; havana_transcript "
               "\"OTTHUMT00000402313\"; havana_transcript_version \"2\"; tag "
               "\"basic\"; transcript_support_level \"4\";\n"

               "12\tensembl_havana\tCDS\t119697659\t119697838\t.\t-\t1\tgene_"
               "id \"ENSG00000122966\"; gene_version \"14\"; transcript_id "
               "\"ENST00000261833\"; transcript_version \"11\"; exon_number "
               "\"45\"; gene_name \"CIT\"; gene_source \"ensembl_havana\"; "
               "gene_biotype \"protein_coding\"; havana_gene "
               "\"OTTHUMG00000134325\"; havana_gene_version \"8\"; "
               "transcript_name \"CIT-001\"; transcript_source "
               "\"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag "
               "\"CCDS\"; ccds_id \"CCDS9192\"; havana_transcript "
               "\"OTTHUMT00000259410\"; havana_transcript_version \"4\"; "
               "protein_id \"ENSP00000261833\"; protein_version \"7\"; tag "
               "\"basic\"; transcript_support_level \"1\";\n";
        gtf.close();
        gmt << "SET1 DDX11L1" << std::endl;
        gmt << "SET2 CIT" << std::endl;
        gmt << "SET3 DDX11L1 CIT" << std::endl;
        gmt << "SET4 ENSG00000223972" << std::endl;
        gmt << "SET5 ENSG00000223973" << std::endl;
        gmt << "SET6 ENSG00000122966" << std::endl;
        gmt.close();
        Reporter reporter("LOG");
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        region = Region(feature, 0, 0, false, false);
        std::vector<std::string> bed_names = {};
        Genotype dummy = GenotypeTest();
        region.generate_regions(gtf_name, gmt_name, bed_names, "", "", "",
                                dummy, reporter);
        SET_BIT(0, not_found.data());
    }
};
TEST_F(REGION_GTF_MULTI_EX, MULTI_NAME_MAP)
{
    // same name mapped to multiple transcript
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(3, found.data());
    SET_BIT(4, found.data());
    // 1 11869 14409
    // 1 15869 16409
    region.update_flag(1, "", 11868, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 11869, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 11870, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 14408, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 14409, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 14410, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    found.front() = 0;
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(3, found.data());
    SET_BIT(5, found.data());
    region.update_flag(1, "", 15868, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(1, "", 15869, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 15870, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 16408, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 16409, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(1, "", 16410, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST_F(REGION_GTF_MULTI_EX, SIMPLE_MULTI)
{
    // 12 11399381 11486678
    // 12 119697659 119697838
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, found.data());
    SET_BIT(2, found.data());
    SET_BIT(3, found.data());
    SET_BIT(6, found.data());
    region.update_flag(12, "", 11399380, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 11399381, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 11399382, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 11486677, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 11486678, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 11486679, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697658, index);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697659, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697660, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697837, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697838, index);
    ASSERT_EQ(index.front(), found.front());
    index.front() = 0;
    region.update_flag(12, "", 119697839, index);
    ASSERT_EQ(index.front(), not_found.front());
}
TEST(REGION_MSIGDB, NAME_CROSS_CHR)
{
    // This should be ok? Transplicing or something like that?
    std::string gtf_name = "Test.gtf";
    std::string gmt_name = "Test.gmt";
    std::ofstream gtf, gmt;
    gtf.open(gtf_name.c_str());
    gmt.open(gmt_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "#!genome - version GRCh38\n"
           "#!genome - date 2013 - 12\n"
           "#!genome - build - accession NCBI : GCA_000001405 .22\n"
           "#!genebuild - last - updated 2016 - 06\n"
           "1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source "
           "\"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
           "havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

           "2\thavana\tgene\t15869\t16409\t.\t+\t.\tgene_id "
           "\"ENSG00000223973\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source "
           "\"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
           "havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

           "12\thavana\tCDS\t11399381\t11486678\t.\t-\t.\tgene_id "
           "\"ENSG00000122966\"; gene_version \"5\"; transcript_id "
           "\"ENST00000538349\"; transcript_version \"5\"; gene_name "
           "\"CIT\"; gene_source \"havana\"; gene_biotype "
           "\"sense_intronic\"; havana_gene \"OTTHUMG00000169117\"; "
           "havana_gene_version \"2\"; transcript_name "
           "\"RP11-711K1.7-001\"; transcript_source \"havana\"; "
           "transcript_biotype \"sense_intronic\"; havana_transcript "
           "\"OTTHUMT00000402313\"; havana_transcript_version \"2\"; tag "
           "\"basic\"; transcript_support_level \"4\";\n"

           "12\tensembl_havana\tCDS\t119697659\t119697838\t.\t-\t1\tgene_"
           "id \"ENSG00000122966\"; gene_version \"14\"; transcript_id "
           "\"ENST00000261833\"; transcript_version \"11\"; exon_number "
           "\"45\"; gene_name \"CIT\"; gene_source \"ensembl_havana\"; "
           "gene_biotype \"protein_coding\"; havana_gene "
           "\"OTTHUMG00000134325\"; havana_gene_version \"8\"; "
           "transcript_name \"CIT-001\"; transcript_source "
           "\"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag "
           "\"CCDS\"; ccds_id \"CCDS9192\"; havana_transcript "
           "\"OTTHUMT00000259410\"; havana_transcript_version \"4\"; "
           "protein_id \"ENSP00000261833\"; protein_version \"7\"; tag "
           "\"basic\"; transcript_support_level \"1\";\n";
    gtf.close();
    gmt << "SET1 DDX11L1" << std::endl;
    gmt << "SET2 CIT" << std::endl;
    gmt.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    Genotype dummy = GenotypeTest();
    try
    {
        region.generate_regions(gtf_name, gmt_name, bed_names, "", "", "",
                                dummy, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}
TEST(REGION_MSIGDB, ID_CROSS_CHR)
{
    // This should be ok? Transplicing or something like that?
    std::string gtf_name = "Test.gtf";
    std::string gmt_name = "Test.gmt";
    std::ofstream gtf, gmt;
    gtf.open(gtf_name.c_str());
    gmt.open(gmt_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "#!genome - version GRCh38\n"
           "#!genome - date 2013 - 12\n"
           "#!genome - build - accession NCBI : GCA_000001405 .22\n"
           "#!genebuild - last - updated 2016 - 06\n"
           "1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source "
           "\"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
           "havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

           "2\thavana\tgene\t15869\t16409\t.\t+\t.\tgene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source "
           "\"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
           "havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

           "12\thavana\tCDS\t11399381\t11486678\t.\t-\t.\tgene_id "
           "\"ENSG00000122966\"; gene_version \"5\"; transcript_id "
           "\"ENST00000538349\"; transcript_version \"5\"; gene_name "
           "\"CIT\"; gene_source \"havana\"; gene_biotype "
           "\"sense_intronic\"; havana_gene \"OTTHUMG00000169117\"; "
           "havana_gene_version \"2\"; transcript_name "
           "\"RP11-711K1.7-001\"; transcript_source \"havana\"; "
           "transcript_biotype \"sense_intronic\"; havana_transcript "
           "\"OTTHUMT00000402313\"; havana_transcript_version \"2\"; tag "
           "\"basic\"; transcript_support_level \"4\";\n"

           "12\tensembl_havana\tCDS\t119697659\t119697838\t.\t-\t1\tgene_"
           "id \"ENSG00000122966\"; gene_version \"14\"; transcript_id "
           "\"ENST00000261833\"; transcript_version \"11\"; exon_number "
           "\"45\"; gene_name \"CIT\"; gene_source \"ensembl_havana\"; "
           "gene_biotype \"protein_coding\"; havana_gene "
           "\"OTTHUMG00000134325\"; havana_gene_version \"8\"; "
           "transcript_name \"CIT-001\"; transcript_source "
           "\"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag "
           "\"CCDS\"; ccds_id \"CCDS9192\"; havana_transcript "
           "\"OTTHUMT00000259410\"; havana_transcript_version \"4\"; "
           "protein_id \"ENSP00000261833\"; protein_version \"7\"; tag "
           "\"basic\"; transcript_support_level \"1\";\n";
    gtf.close();
    gmt << "SET1 DDX11L1" << std::endl;
    gmt << "SET2 CIT" << std::endl;
    gmt.close();
    Reporter reporter("LOG");
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    Region region(feature, 0, 0, false, false);
    std::vector<std::string> bed_names = {};
    Genotype dummy = GenotypeTest();
    try
    {
        region.generate_regions(gtf_name, gmt_name, bed_names, "", "", "",
                                dummy, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}

//TODO  SNP Test (this will require information from Genotype)
#endif // REGION_TEST_HPP
