#ifndef REGION_TEST_HPP
#define REGION_TEST_HPP
#include "genotype.hpp"
#include "global.hpp"
#include "plink_common.hpp"
#include "region.hpp"
#include "reporter.hpp"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <string>

TEST(REGION, SINGLE_INIT)
{
    std::string range = "chr2:1234";
    try
    {
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}
class RegionTest : public Region
{
public:
    bool in_feature(std::string in, const std::vector<std::string>&feature){
        return Region::in_feature(in, feature);
    }
};

TEST(REGION, IN_FEATURE){
    RegionTest test;
    const std::vector<std::string> feature={"gene", "intron", "exon","5'-UTR"};
    ASSERT_TRUE(test.in_feature("gene", feature));
    // case sensitive
    ASSERT_FALSE(test.in_feature("Gene", feature));
    ASSERT_TRUE(test.in_feature("intron", feature));
    ASSERT_TRUE(test.in_feature("exon", feature));
    ASSERT_TRUE(test.in_feature("5'-UTR", feature));
    ASSERT_FALSE(test.in_feature("3'-UTR", feature));
    ASSERT_FALSE(test.in_feature("Genes", feature));
    ASSERT_FALSE(test.in_feature("transcript", feature));

}
TEST(REGION, INVALID_INPUT)
{
    std::string range = "chr2:1234:456";
    try
    {
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

TEST(REGION, SINGLE_RANGE_INIT)
{
    std::string range = "chr2:1234-1357";
    try
    {
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}

TEST(REGION, SINGLE_RANGE_WRONG)
{
    // start must be smaller than end
    std::string range = "chr2:12341-1357";
    try
    {
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION, MULTI_RANGE_INIT)
{
    std::string range = "chr6:369-4321,chr2:1234-1357";
    try
    {
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}
TEST(REGION, MULTI_MIX_INIT)
{
    std::string range = "chr6:369-4321,chr2:1234";
    try
    {
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}
TEST(REGION, MULTI_MORE_MIX_INIT)
{
    std::string range = "chr6:369-4321,chr2:1234,chr1:312345-9437690";
    try
    {
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
    }
    catch (...)
    {
        FAIL();
    }
    SUCCEED();
}
TEST(REGION, WRONG_INPUT)
{
    try
    {
        std::string range = "chr1";
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
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
    try
    {
        std::string range = "chr1:1-2-3";
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION, RANGE_PARSE_PROBLEM)
{
    try
    {
        std::string range = "chr1:1-,2";
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        exclusion_region.clear();
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

TEST(REGION, MULTI_CHROMOSOME)
{
    try
    {
        // we want to test if the chromosome parsing work as expected
        std::string range = "chr1:1,chr10:1,chr2:132,chr20:12,chrX:10";
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        ASSERT_EQ(exclusion_region.size(), CHROM_X + 1);
        exclusion_region.clear();
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(REGION, MULTI_CHROMOSOME_INVALID)
{
    try
    {
        // we want to test if the chromosome parsing work as expected
        std::string range = "chr1:1,chr10:1,chr2:132,chr20:12,chrA:10";
        std::vector<IITree<int, int>> exclusion_region;
        Region::generate_exclusion(exclusion_region, range);
        ASSERT_EQ(exclusion_region.size(), 21);
        exclusion_region.clear();
    }
    catch (...)
    {
        FAIL();
    }
}
class REGION_EX_STRING : public ::testing::Test
{

protected:
    std::vector<IITree<int, int>> exclusion_region;
    void SetUp() override
    {
        std::string range = "chr2:1234";
        Region::generate_exclusion(exclusion_region, range);
    }
    void TearDown() override { exclusion_region.clear(); }
    bool in_region(int chr, int loc)
    {
        return Genotype::within_region(exclusion_region, chr, loc);
    }
};

TEST_F(REGION_EX_STRING, FOUND) { ASSERT_TRUE(in_region(2, 1234)); }
TEST_F(REGION_EX_STRING, WRONG_CHR) { ASSERT_FALSE(in_region(1, 1234)); }
TEST_F(REGION_EX_STRING, WRONG_BIGGER_CHR) { ASSERT_FALSE(in_region(3, 1234)); }
TEST_F(REGION_EX_STRING, BP_TOO_SMALL) { ASSERT_FALSE(in_region(2, 1233)); }
TEST_F(REGION_EX_STRING, BP_TOO_LARGE) { ASSERT_FALSE(in_region(2, 1235)); }
TEST_F(REGION_EX_STRING, SEQUENT_SEARCH)
{
    ASSERT_FALSE(in_region(1, 1234));
    ASSERT_TRUE(in_region(2, 1234));
    ASSERT_FALSE(in_region(2, 1235));
}
TEST_F(REGION_EX_STRING, FINE_SEQUENT_SEARCH)
{
    // new implementation shouldn't be affected by sequence of check
    ASSERT_FALSE(in_region(1, 1234));
    ASSERT_FALSE(in_region(2, 1233));
    ASSERT_TRUE(in_region(2, 1234));
    ASSERT_FALSE(in_region(2, 1235));
    ASSERT_TRUE(in_region(2, 1234));
}

class REGION_EX_STRING_REGION : public ::testing::Test
{

protected:
    std::vector<IITree<int, int>> exclusion_region;
    void SetUp() override
    {
        // inclusion range
        std::string range = "chr2:1234-1357";
        Region::generate_exclusion(exclusion_region, range);
    }
    void TearDown() override { exclusion_region.clear(); }
    bool in_region(int chr, int loc)
    {
        return Genotype::within_region(exclusion_region, chr, loc);
    }
};
TEST_F(REGION_EX_STRING_REGION, SEQ_STRING_WITHIN_RANGE)
{
    // the following should all work
    for (int i = 1234; i < 1357; ++i) ASSERT_TRUE(in_region(2, i));
}
TEST_F(REGION_EX_STRING_REGION, SEQ_STRING_START)
{
    ASSERT_TRUE(in_region(2, 1234));
}
TEST_F(REGION_EX_STRING_REGION, SEQ_STRING_END)
{
    ASSERT_TRUE(in_region(2, 1356));
}
TEST_F(REGION_EX_STRING_REGION, SEQ_UP_BOUND)
{
    ASSERT_FALSE(in_region(2, 1358));
}
TEST_F(REGION_EX_STRING_REGION, SEQ_LOW_BOUND)
{
    ASSERT_FALSE(in_region(2, 1233));
}

class REGION_STRING_MIX : public ::testing::Test
{

protected:
    std::vector<IITree<int, int>> exclusion_region;
    void SetUp() override
    {
        // we need to account for both range and single base input
        // also we want to check chromosome overrun (e.g.
        // Range: chr1:4601-5678,chr2:1357-2468, SNP input 1:5679, 2:134,2:1357)
        // and make sure the input is not in sorted order
        std::string range =
            "chr2:1234-1357,chr1:4601-5678,chr12:314,chr6:98741-102380";

        Region::generate_exclusion(exclusion_region, range);
    }
    void TearDown() override { exclusion_region.clear(); }
    bool in_region(int chr, int loc)
    {
        return Genotype::within_region(exclusion_region, chr, loc);
    }
};
TEST_F(REGION_STRING_MIX, SINGLE_CHECK_FIRST)
{
    ASSERT_TRUE(in_region(2, 1234));
}
TEST_F(REGION_STRING_MIX, SINGLE_CHECK_SECOND)
{
    ASSERT_TRUE(in_region(1, 4890));
}
TEST_F(REGION_STRING_MIX, SINGLE_CHECK_THIRD)
{
    ASSERT_TRUE(in_region(12, 314));
}
TEST_F(REGION_STRING_MIX, SINGLE_CHECK_FORTH)
{
    ASSERT_TRUE(in_region(6, 102379));
}
TEST_F(REGION_STRING_MIX, CHECK_SEQ)
{
    ASSERT_TRUE(in_region(1, 4890));
    ASSERT_TRUE(in_region(2, 1234));
    ASSERT_TRUE(in_region(6, 102379));
    ASSERT_TRUE(in_region(12, 314));
}
TEST_F(REGION_STRING_MIX, CHECK_SEQ_UNSORTED)
{
    ASSERT_TRUE(in_region(2, 1234));
    ASSERT_TRUE(in_region(1, 4890));
    ASSERT_TRUE(in_region(12, 314));
    ASSERT_TRUE(in_region(6, 102379));
}
TEST_F(REGION_STRING_MIX, MIXED_FOUND)
{
    ASSERT_TRUE(in_region(2, 1234));
    // we don't expect the input of the exclusion list is sorted as the listed
    // input might have chromosome ordered in different sequence than we've
    // anticipated
    ASSERT_FALSE(in_region(1, 1));
    ASSERT_TRUE(in_region(12, 314));
    ASSERT_FALSE(in_region(6, 12432145));
}
TEST_F(REGION_STRING_MIX, RUN_OVER)
{
    ASSERT_FALSE(in_region(1, 5679));
    ASSERT_TRUE(in_region(2, 1234));
}

class REGION_BED_MIN_TAB_NO_OVER : public ::testing::Test
{

protected:
    std::vector<IITree<int, int>> exclusion_region;
    void SetUp() override
    {
        std::ofstream bed_file;
        std::string bed_name = path + "/Test.bed";
        bed_file.open(bed_name.c_str());
        if (!bed_file.is_open()) {
            throw std::runtime_error(std::string(
                "Error: Cannot open bed file to write " + bed_name));
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
        Region::generate_exclusion(exclusion_region, bed_name);
    }
    void TearDown() override { exclusion_region.clear(); }
    bool in_region(int chr, int loc)
    {
        return Genotype::within_region(exclusion_region, chr, loc);
    }
};
TEST_F(REGION_BED_MIN_TAB_NO_OVER, CHECK_INCLUSION)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input
    EXPECT_FALSE(in_region(2, 19181 + 1));
    EXPECT_TRUE(in_region(2, 19182 + 1));
    EXPECT_TRUE(in_region(2, 19183 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(2, 32728 + 1));
    EXPECT_FALSE(in_region(2, 32729 + 1));
    EXPECT_FALSE(in_region(2, 94643 + 1));
    EXPECT_TRUE(in_region(2, 94644 + 1));
    EXPECT_TRUE(in_region(2, 94645 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(2, 98554 + 1));
    EXPECT_FALSE(in_region(2, 98555 + 1));
    EXPECT_FALSE(in_region(3, 3208 + 1));
    EXPECT_TRUE(in_region(3, 3209 + 1));
    EXPECT_TRUE(in_region(3, 3210 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(3, 18820 + 1));
    EXPECT_FALSE(in_region(3, 18821 + 1));
    EXPECT_FALSE(in_region(13, 56145 + 1));
    EXPECT_TRUE(in_region(13, 56146 + 1));
    EXPECT_TRUE(in_region(13, 56147 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(13, 67101 + 1));
    EXPECT_FALSE(in_region(13, 67102 + 1));
    EXPECT_FALSE(in_region(21, 9362 + 1));
    EXPECT_TRUE(in_region(21, 9363 + 1));
    EXPECT_TRUE(in_region(21, 9364 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(21, 49430 + 1));
    EXPECT_FALSE(in_region(21, 49431 + 1));
}
TEST_F(REGION_BED_MIN_TAB_NO_OVER, RANDOM_ORDER)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input
    EXPECT_FALSE(in_region(13, 56145 + 1));
    EXPECT_TRUE(in_region(2, 19183 + 1));
    EXPECT_FALSE(in_region(21, 49431 + 1));
    EXPECT_FALSE(in_region(3, 3208 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(2, 32728 + 1));
    EXPECT_FALSE(in_region(2, 32729 + 1));
    EXPECT_TRUE(in_region(2, 94645 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(2, 19182 + 1));
    EXPECT_TRUE(in_region(21, 9363 + 1));
    EXPECT_TRUE(in_region(2, 98554 + 1));
    EXPECT_FALSE(in_region(2, 98555 + 1));
    EXPECT_TRUE(in_region(3, 3209 + 1));
    EXPECT_TRUE(in_region(3, 3210 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(3, 18820 + 1));
    EXPECT_FALSE(in_region(3, 18821 + 1));
    EXPECT_TRUE(in_region(13, 56147 + 1));
    EXPECT_FALSE(in_region(2, 19181 + 1));
    EXPECT_FALSE(in_region(2, 94643 + 1));
    EXPECT_TRUE(in_region(2, 94644 + 1));
    EXPECT_TRUE(in_region(13, 56146 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(13, 67101 + 1));
    EXPECT_FALSE(in_region(13, 67102 + 1));
    EXPECT_FALSE(in_region(21, 9362 + 1));
    EXPECT_TRUE(in_region(21, 9364 + 1));
    // end of bed is non-inclusive
    EXPECT_TRUE(in_region(21, 49430 + 1));
}

class REGION_BED_MIN_TAB : public ::testing::Test
{

protected:
    std::vector<IITree<int, int>> exclusion_region;
    void SetUp() override
    {
        std::ofstream bed_file;
        std::string bed_name = path + "Test.bed";
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
        Region::generate_exclusion(exclusion_region, bed_name);
    }
    void TearDown() override { exclusion_region.clear(); }
    bool in_region(int chr, int loc)
    {
        return Genotype::within_region(exclusion_region, chr, loc);
    }
};

TEST_F(REGION_BED_MIN_TAB, CHECK_INCLUSION_OVERLAPPED)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input

    EXPECT_FALSE(in_region(7, 7079 + 1));
    EXPECT_TRUE(in_region(7, 7080 + 1));
    EXPECT_TRUE(in_region(7, 7081 + 1));
    EXPECT_TRUE(in_region(7, 45053 + 1));
    EXPECT_TRUE(in_region(7, 45054 + 1));
    EXPECT_TRUE(in_region(7, 45055 + 1));
    EXPECT_TRUE(in_region(7, 30303 + 1));
    EXPECT_TRUE(in_region(7, 30305 + 1));
    EXPECT_TRUE(in_region(7, 30306 + 1));
    EXPECT_TRUE(in_region(7, 45722 + 1));
    EXPECT_FALSE(in_region(7, 45723 + 1));
    EXPECT_FALSE(in_region(7, 45724 + 1));

    EXPECT_FALSE(in_region(14, 1693 + 1));
    EXPECT_TRUE(in_region(14, 1694 + 1));
    EXPECT_TRUE(in_region(14, 1695 + 1));
    EXPECT_TRUE(in_region(14, 47284 + 1));
    EXPECT_TRUE(in_region(14, 47285 + 1));
    EXPECT_TRUE(in_region(14, 47286 + 1));
    EXPECT_TRUE(in_region(14, 5224 + 1));
    EXPECT_TRUE(in_region(14, 5225 + 1));
    EXPECT_TRUE(in_region(14, 5226 + 1));
    EXPECT_TRUE(in_region(14, 13101 + 1));
    EXPECT_TRUE(in_region(14, 13102 + 1));
    EXPECT_TRUE(in_region(14, 13103 + 1));
    EXPECT_TRUE(in_region(14, 45657 + 1));
    EXPECT_TRUE(in_region(14, 45658 + 1));
    EXPECT_TRUE(in_region(14, 45659 + 1));
    EXPECT_TRUE(in_region(14, 78547 + 1));
    EXPECT_FALSE(in_region(14, 78548 + 1));
    EXPECT_FALSE(in_region(14, 78549 + 1));

    EXPECT_FALSE(in_region(21, 9362 + 1));
    EXPECT_TRUE(in_region(21, 9363 + 1));
    EXPECT_TRUE(in_region(21, 9364 + 1));
    EXPECT_TRUE(in_region(21, 49430 + 1));
    EXPECT_TRUE(in_region(21, 49431 + 1));
    EXPECT_TRUE(in_region(21, 49432 + 1));
    EXPECT_TRUE(in_region(21, 43439 + 1));
    EXPECT_TRUE(in_region(21, 43440 + 1));
    EXPECT_TRUE(in_region(21, 43441 + 1));
    EXPECT_TRUE(in_region(21, 82119 + 1));
    EXPECT_FALSE(in_region(21, 82120 + 1));
    EXPECT_FALSE(in_region(21, 82121 + 1));
}

class REGION_BED_MIN_SPACE : public ::testing::Test
{

protected:
    std::vector<IITree<int, int>> exclusion_region;
    void SetUp() override
    {
        std::ofstream bed_file;
        std::string bed_name = path + "Test.bed";
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
        Region::generate_exclusion(exclusion_region, bed_name);
    }
    void TearDown() override { exclusion_region.clear(); }
    bool in_region(int chr, int loc)
    {
        return Genotype::within_region(exclusion_region, chr, loc);
    }
};

TEST_F(REGION_BED_MIN_SPACE, CHECK_INCLUSION_OVERLAPPED)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input

    EXPECT_FALSE(in_region(7, 7079 + 1));
    EXPECT_TRUE(in_region(7, 7080 + 1));
    EXPECT_TRUE(in_region(7, 7081 + 1));
    EXPECT_TRUE(in_region(7, 45053 + 1));
    EXPECT_TRUE(in_region(7, 45054 + 1));
    EXPECT_TRUE(in_region(7, 45055 + 1));
    EXPECT_TRUE(in_region(7, 30303 + 1));
    EXPECT_TRUE(in_region(7, 30305 + 1));
    EXPECT_TRUE(in_region(7, 30306 + 1));
    EXPECT_TRUE(in_region(7, 45722 + 1));
    EXPECT_FALSE(in_region(7, 45723 + 1));
    EXPECT_FALSE(in_region(7, 45724 + 1));

    EXPECT_FALSE(in_region(14, 1693 + 1));
    EXPECT_TRUE(in_region(14, 1694 + 1));
    EXPECT_TRUE(in_region(14, 1695 + 1));
    EXPECT_TRUE(in_region(14, 47284 + 1));
    EXPECT_TRUE(in_region(14, 47285 + 1));
    EXPECT_TRUE(in_region(14, 47286 + 1));
    EXPECT_TRUE(in_region(14, 5224 + 1));
    EXPECT_TRUE(in_region(14, 5225 + 1));
    EXPECT_TRUE(in_region(14, 5226 + 1));
    EXPECT_TRUE(in_region(14, 13101 + 1));
    EXPECT_TRUE(in_region(14, 13102 + 1));
    EXPECT_TRUE(in_region(14, 13103 + 1));
    EXPECT_TRUE(in_region(14, 45657 + 1));
    EXPECT_TRUE(in_region(14, 45658 + 1));
    EXPECT_TRUE(in_region(14, 45659 + 1));
    EXPECT_TRUE(in_region(14, 78547 + 1));
    EXPECT_FALSE(in_region(14, 78548 + 1));
    EXPECT_FALSE(in_region(14, 78549 + 1));

    EXPECT_FALSE(in_region(21, 9362 + 1));
    EXPECT_TRUE(in_region(21, 9363 + 1));
    EXPECT_TRUE(in_region(21, 9364 + 1));
    EXPECT_TRUE(in_region(21, 49430 + 1));
    EXPECT_TRUE(in_region(21, 49431 + 1));
    EXPECT_TRUE(in_region(21, 49432 + 1));
    EXPECT_TRUE(in_region(21, 43439 + 1));
    EXPECT_TRUE(in_region(21, 43440 + 1));
    EXPECT_TRUE(in_region(21, 43441 + 1));
    EXPECT_TRUE(in_region(21, 82119 + 1));
    EXPECT_FALSE(in_region(21, 82120 + 1));
    EXPECT_FALSE(in_region(21, 82121 + 1));
}

class REGION_BED_5_COLUMN : public ::testing::Test
{

protected:
    std::vector<IITree<int, int>> exclusion_region;
    void SetUp() override
    {
        // not enough for stand yet
        std::ofstream bed_file;
        std::string bed_name = path + "Test.bed";
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
        Region::generate_exclusion(exclusion_region, bed_name);
    }
    void TearDown() override { exclusion_region.clear(); }
    bool in_region(int chr, int loc)
    {
        return Genotype::within_region(exclusion_region, chr, loc);
    }
};

TEST_F(REGION_BED_5_COLUMN, CHECK_INCLUSION_OVERLAPPED)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input

    EXPECT_FALSE(in_region(7, 7079 + 1));
    EXPECT_TRUE(in_region(7, 7080 + 1));
    EXPECT_TRUE(in_region(7, 7081 + 1));
    EXPECT_TRUE(in_region(7, 45053 + 1));
    EXPECT_TRUE(in_region(7, 45054 + 1));
    EXPECT_TRUE(in_region(7, 45055 + 1));
    EXPECT_TRUE(in_region(7, 30303 + 1));
    EXPECT_TRUE(in_region(7, 30305 + 1));
    EXPECT_TRUE(in_region(7, 30306 + 1));
    EXPECT_TRUE(in_region(7, 45722 + 1));
    EXPECT_FALSE(in_region(7, 45723 + 1));
    EXPECT_FALSE(in_region(7, 45724 + 1));

    EXPECT_FALSE(in_region(14, 1693 + 1));
    EXPECT_TRUE(in_region(14, 1694 + 1));
    EXPECT_TRUE(in_region(14, 1695 + 1));
    EXPECT_TRUE(in_region(14, 47284 + 1));
    EXPECT_TRUE(in_region(14, 47285 + 1));
    EXPECT_TRUE(in_region(14, 47286 + 1));
    EXPECT_TRUE(in_region(14, 5224 + 1));
    EXPECT_TRUE(in_region(14, 5225 + 1));
    EXPECT_TRUE(in_region(14, 5226 + 1));
    EXPECT_TRUE(in_region(14, 13101 + 1));
    EXPECT_TRUE(in_region(14, 13102 + 1));
    EXPECT_TRUE(in_region(14, 13103 + 1));
    EXPECT_TRUE(in_region(14, 45657 + 1));
    EXPECT_TRUE(in_region(14, 45658 + 1));
    EXPECT_TRUE(in_region(14, 45659 + 1));
    EXPECT_TRUE(in_region(14, 78547 + 1));
    EXPECT_FALSE(in_region(14, 78548 + 1));
    EXPECT_FALSE(in_region(14, 78549 + 1));

    EXPECT_FALSE(in_region(21, 9362 + 1));
    EXPECT_TRUE(in_region(21, 9363 + 1));
    EXPECT_TRUE(in_region(21, 9364 + 1));
    EXPECT_TRUE(in_region(21, 49430 + 1));
    EXPECT_TRUE(in_region(21, 49431 + 1));
    EXPECT_TRUE(in_region(21, 49432 + 1));
    EXPECT_TRUE(in_region(21, 43439 + 1));
    EXPECT_TRUE(in_region(21, 43440 + 1));
    EXPECT_TRUE(in_region(21, 43441 + 1));
    EXPECT_TRUE(in_region(21, 82119 + 1));
    EXPECT_FALSE(in_region(21, 82120 + 1));
    EXPECT_FALSE(in_region(21, 82121 + 1));
}
class REGION_BED_WITH_STRAND : public ::testing::Test
{
    // For exclusion, strand information should not alter result (window padding
    // should all be 0)
protected:
    std::vector<IITree<int, int>> exclusion_region;
    void SetUp() override
    {
        std::ofstream bed_file;
        std::string bed_name = path + "Test.bed";
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
        Region::generate_exclusion(exclusion_region, bed_name);
    }
    void TearDown() override { exclusion_region.clear(); }
    bool in_region(int chr, int loc)
    {
        return Genotype::within_region(exclusion_region, chr, loc);
    }
};

TEST_F(REGION_BED_WITH_STRAND, CHECK_INCLUSION_OVERLAPPED)
{
    // NOTE: +1 here because the number is w.r.t. the bed file, which has a 0
    // base, but check_exclusion expect a 1 based input

    EXPECT_FALSE(in_region(7, 7079 + 1));
    EXPECT_TRUE(in_region(7, 7080 + 1));
    EXPECT_TRUE(in_region(7, 7081 + 1));
    EXPECT_TRUE(in_region(7, 45053 + 1));
    EXPECT_TRUE(in_region(7, 45054 + 1));
    EXPECT_TRUE(in_region(7, 45055 + 1));
    EXPECT_TRUE(in_region(7, 30303 + 1));
    EXPECT_TRUE(in_region(7, 30305 + 1));
    EXPECT_TRUE(in_region(7, 30306 + 1));
    EXPECT_TRUE(in_region(7, 45722 + 1));
    EXPECT_FALSE(in_region(7, 45723 + 1));
    EXPECT_FALSE(in_region(7, 45724 + 1));

    EXPECT_FALSE(in_region(14, 1693 + 1));
    EXPECT_TRUE(in_region(14, 1694 + 1));
    EXPECT_TRUE(in_region(14, 1695 + 1));
    EXPECT_TRUE(in_region(14, 47284 + 1));
    EXPECT_TRUE(in_region(14, 47285 + 1));
    EXPECT_TRUE(in_region(14, 47286 + 1));
    EXPECT_TRUE(in_region(14, 5224 + 1));
    EXPECT_TRUE(in_region(14, 5225 + 1));
    EXPECT_TRUE(in_region(14, 5226 + 1));
    EXPECT_TRUE(in_region(14, 13101 + 1));
    EXPECT_TRUE(in_region(14, 13102 + 1));
    EXPECT_TRUE(in_region(14, 13103 + 1));
    EXPECT_TRUE(in_region(14, 45657 + 1));
    EXPECT_TRUE(in_region(14, 45658 + 1));
    EXPECT_TRUE(in_region(14, 45659 + 1));
    EXPECT_TRUE(in_region(14, 78547 + 1));
    EXPECT_FALSE(in_region(14, 78548 + 1));
    EXPECT_FALSE(in_region(14, 78549 + 1));

    EXPECT_FALSE(in_region(21, 9362 + 1));
    EXPECT_TRUE(in_region(21, 9363 + 1));
    EXPECT_TRUE(in_region(21, 9364 + 1));
    EXPECT_TRUE(in_region(21, 49430 + 1));
    EXPECT_TRUE(in_region(21, 49431 + 1));
    EXPECT_TRUE(in_region(21, 49432 + 1));
    EXPECT_TRUE(in_region(21, 43439 + 1));
    EXPECT_TRUE(in_region(21, 43440 + 1));
    EXPECT_TRUE(in_region(21, 43441 + 1));
    EXPECT_TRUE(in_region(21, 82119 + 1));
    EXPECT_FALSE(in_region(21, 82120 + 1));
    EXPECT_FALSE(in_region(21, 82121 + 1));
}

// now test different type of malform BED file
TEST(REGION_MALFORM_BED, NOT_ENOUGH_COLUMN)
{
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
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
    std::vector<IITree<int, int>> exclusion_region;
    try
    {
        // we want to penalize any form of malformed input
        Region::generate_exclusion(exclusion_region, bed_name);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    exclusion_region.clear();
}

TEST(REGION_MALFORM_BED, INCONSISTEN_COLUMN_STRAND)
{
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
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
    std::vector<IITree<int, int>> exclusion_region;
    try
    {
        // we want to penalize any form of malformed input
        Region::generate_exclusion(exclusion_region, bed_name);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    exclusion_region.clear();
}
TEST(REGION_MALFORM_BED, NOT_FOUND)
{
    std::string bed_name = path + "404.bed";
    std::vector<IITree<int, int>> exclusion_region;
    try
    {
        // we want to penalize any form of malformed input
        Region::generate_exclusion(exclusion_region, bed_name);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    exclusion_region.clear();
}
TEST(REGION_MALFORM_BED, UNSUPPORTED_STRAND)
{
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
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
    std::vector<IITree<int, int>> exclusion_region;
    try
    {
        Region::generate_exclusion(exclusion_region, bed_name);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    exclusion_region.clear();
}
TEST(REGION_MALFORM_BED, NEGATIVE_COORDINATE)
{
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
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
    std::vector<IITree<int, int>> exclusion_region;
    try
    {
        Region::generate_exclusion(exclusion_region, bed_name);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    exclusion_region.clear();
}

TEST(REGION_STD_BED_INPUT, NO_RUN)
{
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
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
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    std::vector<std::string> region_names;
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string gtf = "";
    std::string msigdb = "";
    std::string snp_set = "";
    std::vector<std::string> bed;
    std::string background = "";
    Reporter reporter(std::string(path + "LOG"));
    std::vector<IITree<int, int>> gene_sets;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    size_t num_regions = Region::generate_regions(
        gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
        genome_wide_background, gtf, msigdb, bed, snp_set, background, 22,
        reporter);
    ASSERT_STREQ(region_names[0].c_str(), "Base");
    ASSERT_STREQ(region_names[1].c_str(), "Background");
    // there will always be 2 regions
    ASSERT_EQ(num_regions, 2);
    std::vector<uintptr_t> not_found = {0};
    std::vector<uintptr_t> found = {0};
    SET_BIT(0, not_found.data());
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    const size_t required_size = BITCT_TO_WORDCT(num_regions);
    std::vector<uintptr_t> index(required_size, 0);
    // this is a SNP found in the bed file, but as we have not generated the
    // region (we haven't use the bed file), this will always be considered as
    // not found
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             5, 50533 + 1, genome_wide_background);
    ASSERT_EQ(index, not_found);
}


TEST(REGION_STD_BED_INPUT, WITH_HEADER_TRACE)
{
    // test BED file with valid header
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "track name=pairedReads description=\"Clone Paired Reads\" "
                "useScore=1\n"
             << "2 19182 32729 . . .\n";
    bed_file.close();
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    std::vector<std::string> region_names;
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string gtf = "";
    std::string msigdb = "";
    std::string snp_set = "";
    std::vector<std::string> bed;
    std::string background = "";
    Reporter reporter(std::string(path + "LOG"));
    std::vector<IITree<int, int>> gene_sets;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf, msigdb, bed, snp_set, background, 22,
                                 reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(REGION_STD_BED_INPUT, WITH_HEADER_BROWSER)
{
    // test BED file with valid header
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "browser position chr7:127471196-127495720\n"
             << "browser hide all\n"
                "useScore=1\n"
             << "2 19182 32729 . . .\n";
    bed_file.close();
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    std::vector<std::string> region_names;
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string gtf = "";
    std::string msigdb = "";
    std::string snp_set = "";
    std::vector<std::string> bed;
    std::string background = "";
    Reporter reporter(std::string(path + "LOG"));
    std::vector<IITree<int, int>> gene_sets;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf, msigdb, bed, snp_set, background, 22,
                                 reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}


TEST(REGION_STD_BED_INPUT, EXCLUSION_WITH_HEADER_BROWSER)
{
    // test BED file with valid header
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "browser position chr7:127471196-127495720\n"
             << "browser hide all\n"
                "useScore=1\n"
             << "2 19182 32729 . . .\n";
    bed_file.close();
    std::vector<IITree<int, int>> exclusion_region;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    try
    {
        Region::generate_exclusion(exclusion_region, bed_name);
        SUCCEED();
    }
    catch (const std::runtime_error& re)
    {
        std::cerr << re.what() << std::endl;
        FAIL();
    }
    FAIL();
}
TEST(REGION_STD_BED_INPUT, EXCLUSION_WITH_HEADER_TRACK)
{
    // test BED file with valid header
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "track name=pairedReads description=\"Clone Paired Reads\" "
                "useScore=1\n"
             << "2 19182 32729 . . .\n";
    bed_file.close();
    std::vector<IITree<int, int>> exclusion_region;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    try
    {
        Region::generate_exclusion(exclusion_region, bed_name);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}

// now test different type of malform BED file
TEST(REGION_MALFORM_BED, INVALID_HEADER_FOR_EXCLUSION)
{
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "#CHR START END\n"
             << "2 19182 32729 . . .\n";
    bed_file.close();
    std::vector<IITree<int, int>> exclusion_region;
    try
    {
        // we want to penalize any form of malformed input
        Region::generate_exclusion(exclusion_region, bed_name);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

TEST(REGION_MALFORM_BED, INVALID_HEADER_FOR_SET_SELECT)
{
    // test BED file with invalid header
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
    bed_file.open(bed_name.c_str());
    if (!bed_file.is_open()) {
        throw std::runtime_error("Error: Cannot open bed file");
    }
    //  now generate the output required
    bed_file << "#CHR START END\n"
             << "2 19182 32729 . . .\n";
    bed_file.close();
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    std::vector<std::string> region_names;
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string gtf = "";
    std::string msigdb = "";
    std::string snp_set = "";
    std::vector<std::string> bed = {bed_name};
    std::string background = "";
    Reporter reporter(std::string(path + "LOG"));
    std::vector<IITree<int, int>> gene_sets;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    try
    {
        // malformed anything are considered as fatal
        Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf, msigdb, bed, snp_set, background, 22,
            reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

class REGION_STD_BED : public ::testing::Test
{
protected:
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<uintptr_t> not_found = {0};
    std::vector<uintptr_t> found = {0};
    std::vector<std::string> region_names;
    size_t num_regions;
    size_t required_size;
    bool genome_wide_background = false;
    void SetUp() override
    {
        std::ofstream bed_file;
        std::string bed_name = path + "Test.bed";
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
        Reporter reporter(std::string(path + "LOG"));
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        int window_5 = 0;
        int window_3 = 0;
        std::string gtf = "";
        std::string msigdb = "";
        std::string snp_set = "";
        std::string background = "";
        std::vector<std::string> bed_names = {bed_name};
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf, msigdb, bed_names, snp_set, background,
            22, reporter);
        SET_BIT(0, not_found.data());
        SET_BIT(0, found.data());
        // 2 because 1 is reserved for background and 2 is the first set
        SET_BIT(2, found.data());
        required_size = BITCT_TO_WORDCT(num_regions);
    }
    std::vector<uintptr_t> get_flag(const int chr, const int bp)
    {
        std::vector<uintptr_t> index(required_size, 0);
        Genotype::construct_flag("", gene_sets, snp_in_sets, index,
                                 required_size, chr, bp,
                                 genome_wide_background);
        return index;
    }
};
TEST_F(REGION_STD_BED, CHECK_INCLUSION_OVERLAPPED)
{
    // with standard input, we can no longer use check_exclusion function as
    // that always uses the base region, which doesn't contain any boundary
    // instead, we must use the update flag function

    EXPECT_EQ(get_flag(7, 7079 + 1).front(), not_found.front());
    EXPECT_EQ(get_flag(7, 7080 + 1).front(), found.front());
    EXPECT_EQ(get_flag(7, 7081 + 1).front(), found.front());
    EXPECT_EQ(get_flag(7, 45053 + 1).front(), found.front());
    EXPECT_EQ(get_flag(7, 45054 + 1).front(), found.front());
    EXPECT_EQ(get_flag(7, 45055 + 1).front(), found.front());
    EXPECT_EQ(get_flag(7, 30303 + 1).front(), found.front());
    EXPECT_EQ(get_flag(7, 30305 + 1).front(), found.front());
    EXPECT_EQ(get_flag(7, 30306 + 1).front(), found.front());
    EXPECT_EQ(get_flag(7, 45722 + 1).front(), found.front());
    EXPECT_EQ(get_flag(7, 45723 + 1).front(), not_found.front());
    EXPECT_EQ(get_flag(7, 45724 + 1).front(), not_found.front());
    EXPECT_EQ(get_flag(14, 1693 + 1).front(), not_found.front());
    EXPECT_EQ(get_flag(14, 1695 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 47284 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 47285 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 47286 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 5224 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 5225 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 5226 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 13101 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 13102 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 13103 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 45657 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 45658 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 45659 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 78547 + 1).front(), found.front());
    EXPECT_EQ(get_flag(14, 78548 + 1).front(), not_found.front());
    EXPECT_EQ(get_flag(14, 78549 + 1).front(), not_found.front());
}
TEST_F(REGION_STD_BED, MID_NOT_FOUND)
{
    EXPECT_EQ(get_flag(19, 39329 + 1).front(), not_found.front());
    EXPECT_EQ(get_flag(20, 64037 + 1).front(), found.front());
}
class REGION_STD_BED_PAD : public ::testing::Test
{
    // test window padding
protected:
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<uintptr_t> not_found = {0};
    std::vector<uintptr_t> found = {0};
    std::vector<std::string> region_names;
    size_t num_regions;
    size_t required_size;
    bool genome_wide_background = false;
    void SetUp() override
    {
        std::ofstream bed_file;
        std::string bed_name = path + "Test.bed";
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
        Reporter reporter(std::string(path + "LOG"));
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        int window_5 = 10;
        int window_3 = 20;
        std::string gtf = "";
        std::string msigdb = "";
        std::string snp_set = "";
        std::string background = "";
        std::vector<std::string> bed_names = {bed_name};
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf, msigdb, bed_names, snp_set, background,
            22, reporter);
        SET_BIT(0, not_found.data());
        SET_BIT(0, found.data());
        // 2 because 1 is reserved for background and 2 is the first set
        SET_BIT(2, found.data());
        required_size = BITCT_TO_WORDCT(num_regions);
    }
    std::vector<uintptr_t> get_flag(const int chr, const int bp)
    {
        std::vector<uintptr_t> index(required_size, 0);
        Genotype::construct_flag("", gene_sets, snp_in_sets, index,
                                 required_size, chr, bp,
                                 genome_wide_background);
        return index;
    }
};
TEST_F(REGION_STD_BED_PAD, CHECK_PAD)
{
    // We will see how the padding change the inclusion criteria
    // this SNP doesn't contain the strand info, we should assume the start
    // is the 5' end
    // we have pad 10 bp to the 5' and 20 to the 3'
    EXPECT_EQ(get_flag(3, 29863 + 1 - 11).front(), not_found.front());
    EXPECT_EQ(get_flag(3, 29863 + 1 - 10).front(), found.front());
    EXPECT_EQ(get_flag(3, 29863 + 1).front(), found.front());
    EXPECT_EQ(get_flag(3, 38285 + 1).front(), found.front());
    EXPECT_EQ(get_flag(3, 38285 + 1 + 19).front(), found.front());
    EXPECT_EQ(get_flag(3, 38285 + 1 + 20).front(), not_found.front());
    EXPECT_EQ(get_flag(4, 20139 + 1 - 11).front(), not_found.front());
    EXPECT_EQ(get_flag(4, 20139 + 1 - 10).front(), found.front());
    EXPECT_EQ(get_flag(4, 20139 + 1).front(), found.front());
    EXPECT_EQ(get_flag(4, 97433 + 1).front(), found.front());
    EXPECT_EQ(get_flag(4, 97433 + 1 + 19).front(), found.front());
    EXPECT_EQ(get_flag(4, 97433 + 1 + 20).front(), not_found.front());
    // 6 34611 45099 . . -
    EXPECT_EQ(get_flag(6, 34611 + 1 - 21).front(), not_found.front());
    EXPECT_EQ(get_flag(6, 34611 + 1 - 20).front(), found.front());
    EXPECT_EQ(get_flag(6, 34611 + 1).front(), found.front());
    EXPECT_EQ(get_flag(6, 45099 + 1).front(), found.front());
    EXPECT_EQ(get_flag(6, 45099 + 1 + 9).front(), found.front());
    EXPECT_EQ(get_flag(6, 45099 + 1 + 10).front(), not_found.front());
}

TEST(REGION_MULTI_BED, CHECK_NAME)
{
    Reporter reporter(std::string(path + "LOG"));
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
    std::string second_bed_name = path + "Test2.bed";
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
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string gtf = "";
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> bed_names = {std::string(bed_name + ":Name"),
                                          second_bed_name};
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    size_t num_regions = Region::generate_regions(
        gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
        genome_wide_background, gtf, msigdb, bed_names, snp_set, background, 22,
        reporter);
    ASSERT_EQ(num_regions, 4);
    ASSERT_STREQ(region_names[0].c_str(), "Base");
    ASSERT_STREQ(region_names[1].c_str(), "Background");
    ASSERT_STREQ(region_names[2].c_str(), "Name");
    ASSERT_STREQ(region_names[3].c_str(),
                 std::string(path + "Test2.bed").c_str());
}

TEST(REGION_MULTI_BED, CHECK_NAME2)
{
    Reporter reporter(std::string(path + "LOG"));
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
    std::string second_bed_name = path + "Test2.bed";
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

    std::vector<std::string> bed_names = {
        bed_name, std::string(second_bed_name + ":Name"),
    };
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string gtf = "";
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::vector<IITree<int, int>> gene_sets;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    size_t num_regions = Region::generate_regions(
        gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
        genome_wide_background, gtf, msigdb, bed_names, snp_set, background, 22,
        reporter);
    ASSERT_EQ(num_regions, 4);
    ASSERT_STREQ(region_names[0].c_str(), "Base");
    ASSERT_STREQ(region_names[1].c_str(), "Background");
    ASSERT_STREQ(region_names[2].c_str(),
                 std::string(path + "Test.bed").c_str());
    ASSERT_STREQ(region_names[3].c_str(), "Name");
}

// gtf read
// msigdb read

// Any error in the GTF file will lead to throw
TEST(REGION_GTF_BASIC, NOT_EXIST)
{
    std::string gtf_name = path + "Test.gtf";
    remove(gtf_name.c_str());
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

TEST(REGION_GTF_BASIC, EMPTY)
{
    std::string gtf_name = path + "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf.close();
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

TEST(REGION_GTF_BASIC, ALL_REGION_REMOVE)
{
    std::string gtf_name = path + "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

TEST(REGION_GTF_BASIC, MALFORMAT_SPACE)
{
    std::string gtf_name = path + "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
           "1 havana gene 11869 14409 . + . gene_id "
           "\"ENSG00000223972\"; "
           "gene_version \"5\"; gene_name \"DDX11L1\"; gene_source \"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n";
    gtf.close();
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, NEGATIVE_COORDINATE)
{
    std::string gtf_name = path + "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, BIGGER_START)
{
    std::string gtf_name = path + "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::vector<IITree<int, int>> gene_sets;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<std::string> bed_names = {};
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, UNDEFINED_STRAND)
{

    std::string gtf_name = path + "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, TAB_ATTRIBUTE)
{

    std::string gtf_name = path + "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}
TEST(REGION_GTF_BASIC, NO_GENE_ID)
{

    std::string gtf_name = path + "Test.gtf";
    std::ofstream gtf;
    gtf.open(gtf_name.c_str());
    gtf << "#!genome-build GRCh38.p7\n"
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
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
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<uintptr_t> not_found = {0};
    size_t num_regions;
    size_t required_size;
    // make sure genome_wide_background is true, or each set will
    // have a different not_found bit
    bool genome_wide_background = true;
    void SetUp() override
    {
        std::string gtf_name = path + "Test.gtf";
        std::string gmt_name = path + "Test.gmt";
        std::ofstream gtf, gmt;
        gtf.open(gtf_name.c_str());
        gmt.open(gmt_name.c_str());
        gtf << "#!genome-build GRCh38.p7\n"
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
        Reporter reporter(std::string(path + "LOG"));
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        int window_5 = 0;
        int window_3 = 0;
        std::string msigdb = "";
        std::string snp_set = "";
        std::string background = "";
        std::vector<std::string> region_names;
        std::vector<std::string> bed_names = {};
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set,
            background, 22, reporter);
        SET_BIT(0, not_found.data());
        // because we use genome_wide_background, which should have same bit set
        // as base
        SET_BIT(1, not_found.data());
        required_size = BITCT_TO_WORDCT(num_regions);
    }
    std::vector<uintptr_t> get_flag(const int chr, const int bp)
    {
        std::vector<uintptr_t> index(required_size, 0);
        Genotype::construct_flag("", gene_sets, snp_in_sets, index,
                                 required_size, chr, bp,
                                 genome_wide_background);
        return index;
    }
};
TEST_F(REGION_GTF_FEATURE, FEATURE_FILTER) { ASSERT_EQ(num_regions, 8); }
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET1)
{
    std::vector<uintptr_t> found = {0};
    // both base and background are set
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(2, found.data());
    SET_BIT(7, found.data());
    // 1 havana gene 11869 14409
    ASSERT_EQ(get_flag(1, 11868).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 11869).front(), found.front());
    ASSERT_EQ(get_flag(1, 11870).front(), found.front());
    ASSERT_EQ(get_flag(1, 14408).front(), found.front());
    ASSERT_EQ(get_flag(1, 14409).front(), found.front());
    ASSERT_EQ(get_flag(1, 14410).front(), not_found.front());
}
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET2)
{
    // should all be failed as filtered by feature
    ASSERT_EQ(get_flag(1, 15868).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 15869).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 15870).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 16408).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 16409).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 16410).front(), not_found.front());
}
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET3)
{
    // should all be failed as filtered by feature
    ASSERT_EQ(get_flag(12, 11399380).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11399381).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11399382).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11486677).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11486678).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11486679).front(), not_found.front());
}
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET4)
{
    std::vector<uintptr_t> found = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    // + 1 because 0 base, otherwise, we need to +2 to set number as base and
    // background
    SET_BIT(4 + 1, found.data());
    SET_BIT(6 + 1, found.data());
    ASSERT_EQ(get_flag(12, 119697658).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 119697659).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697660).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697837).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697838).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697839).front(), not_found.front());
}
TEST_F(REGION_GTF_FEATURE, FOUND_SNP_SET5)
{
    std::vector<uintptr_t> found = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(5 + 1, found.data());
    ASSERT_EQ(get_flag(15, 55320274).front(), not_found.front());
    ASSERT_EQ(get_flag(15, 55320275).front(), found.front());
    ASSERT_EQ(get_flag(15, 55320276).front(), found.front());
    ASSERT_EQ(get_flag(15, 55320409).front(), found.front());
    ASSERT_EQ(get_flag(15, 55320410).front(), found.front());
    ASSERT_EQ(get_flag(15, 55320411).front(), not_found.front());
}
class REGION_GTF_PAD : public ::testing::Test
{
protected:
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<uintptr_t> not_found = {0};
    size_t num_regions;
    size_t required_size;
    // make sure genome_wide_background is true, or each set will
    // have a different not_found bit
    bool genome_wide_background = true;
    void SetUp() override
    {
        std::string gtf_name = path + "Test.gtf";
        std::string gmt_name = path + "Test.gmt";
        std::ofstream gtf, gmt;
        gtf.open(gtf_name.c_str());
        gmt.open(gmt_name.c_str());
        gtf << "#!genome-build GRCh38.p7\n"
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
        Reporter reporter(std::string(path + "LOG"));
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        int window_5 = 10;
        int window_3 = 20;
        std::string msigdb = "";
        std::string snp_set = "";
        std::string background = "";
        std::vector<std::string> region_names;
        std::vector<std::string> bed_names = {};
        std::unordered_map<std::string, std::vector<int>> snp_in_sets;
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set,
            background, 22, reporter);
        SET_BIT(0, not_found.data());
        // because we use genome_wide_background, which should have same bit set
        // as base
        SET_BIT(1, not_found.data());
        required_size = BITCT_TO_WORDCT(num_regions);
    }
    std::vector<uintptr_t> get_flag(const int chr, const int bp)
    {
        std::vector<uintptr_t> index(required_size, 0);
        Genotype::construct_flag("", gene_sets, snp_in_sets, index,
                                 required_size, chr, bp,
                                 genome_wide_background);
        return index;
    }
};
TEST_F(REGION_GTF_PAD, FEATURE_FILTER) { ASSERT_EQ(num_regions, 8); }
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET1)
{
    std::vector<uintptr_t> found = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(1 + 1, found.data());
    SET_BIT(6 + 1, found.data());
    // 1 havana gene 11869 14409
    ASSERT_EQ(get_flag(1, 11858).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 11859).front(), found.front());
    ASSERT_EQ(get_flag(1, 11860).front(), found.front());
    ASSERT_EQ(get_flag(1, 14428).front(), found.front());
    ASSERT_EQ(get_flag(1, 14429).front(), found.front());
    ASSERT_EQ(get_flag(1, 14430).front(), not_found.front());
}
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET2)
{
    // should all be failed
    ASSERT_EQ(get_flag(1, 15858).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 15859).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 15860).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 16428).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 16429).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 16430).front(), not_found.front());
}
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET3)
{
    // should all be failed
    ASSERT_EQ(get_flag(12, 11399360).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11399361).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11399362).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11486687).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11486688).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11486689).front(), not_found.front());
}
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET4)
{
    std::vector<uintptr_t> found = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(4 + 1, found.data());
    SET_BIT(6 + 1, found.data());
    // 1 havana gene 11869 14409
    ASSERT_EQ(get_flag(12, 119697638).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 119697639).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697640).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697847).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697848).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697849).front(), not_found.front());
}
TEST_F(REGION_GTF_PAD, FOUND_SNP_SET5)
{
    std::vector<uintptr_t> found = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(5 + 1, found.data());
    // 1 havana gene 11869 14409
    ASSERT_EQ(get_flag(15, 55320264).front(), not_found.front());
    ASSERT_EQ(get_flag(15, 55320265).front(), found.front());
    ASSERT_EQ(get_flag(15, 55320266).front(), found.front());
    ASSERT_EQ(get_flag(15, 55320429).front(), found.front());
    ASSERT_EQ(get_flag(15, 55320430).front(), found.front());
    ASSERT_EQ(get_flag(15, 55320431).front(), not_found.front());
}
// need class for following
class REGION_GTF_MULTI_EX : public ::testing::Test
{
protected:
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<uintptr_t> not_found = {0};
    size_t num_regions;
    size_t required_size;
    // make sure genome_wide_background is true, or each set will
    // have a different not_found bit
    bool genome_wide_background = true;
    void SetUp() override
    {
        std::string gtf_name = path + "Test.gtf";
        std::string gmt_name = path + "Test.gmt";
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
        Reporter reporter(std::string(path + "LOG"));
        std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                            "CDS"};
        int window_5 = 0;
        int window_3 = 0;
        std::string msigdb = "";
        std::string snp_set = "";
        std::string background = "";
        std::vector<std::string> region_names;
        std::vector<std::string> bed_names = {};
        std::unordered_map<std::string, std::vector<int>> snp_in_sets;
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set,
            background, 22, reporter);
        SET_BIT(0, not_found.data());
        // because we use genome_wide_background, which should have same bit set
        // as base
        SET_BIT(1, not_found.data());
        required_size = BITCT_TO_WORDCT(num_regions);
    }
    std::vector<uintptr_t> get_flag(const int chr, const int bp)
    {
        std::vector<uintptr_t> index(required_size, 0);
        Genotype::construct_flag("", gene_sets, snp_in_sets, index,
                                 required_size, chr, bp,
                                 genome_wide_background);
        return index;
    }
};
TEST_F(REGION_GTF_MULTI_EX, MULTI_NAME_MAP)
{
    // same name mapped to multiple transcript
    std::vector<uintptr_t> found = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(1 + 1, found.data());
    SET_BIT(3 + 1, found.data());
    SET_BIT(4 + 1, found.data());
    // 1 11869 14409
    // 1 15869 16409
    ASSERT_EQ(get_flag(1, 11868).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 11869).front(), found.front());
    ASSERT_EQ(get_flag(1, 11870).front(), found.front());
    ASSERT_EQ(get_flag(1, 14408).front(), found.front());
    ASSERT_EQ(get_flag(1, 14409).front(), found.front());
    ASSERT_EQ(get_flag(1, 14410).front(), not_found.front());
    found.front() = 0;
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(1 + 1, found.data());
    SET_BIT(3 + 1, found.data());
    SET_BIT(5 + 1, found.data());
    ASSERT_EQ(get_flag(1, 15868).front(), not_found.front());
    ASSERT_EQ(get_flag(1, 15869).front(), found.front());
    ASSERT_EQ(get_flag(1, 15870).front(), found.front());
    ASSERT_EQ(get_flag(1, 16408).front(), found.front());
    ASSERT_EQ(get_flag(1, 16409).front(), found.front());
    ASSERT_EQ(get_flag(1, 16410).front(), not_found.front());
}
TEST_F(REGION_GTF_MULTI_EX, SIMPLE_MULTI)
{
    // 12 11399381 11486678
    // 12 119697659 119697838
    std::vector<uintptr_t> found = {0};
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(2 + 1, found.data());
    SET_BIT(3 + 1, found.data());
    SET_BIT(6 + 1, found.data());
    ASSERT_EQ(get_flag(12, 11399380).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 11399381).front(), found.front());
    ASSERT_EQ(get_flag(12, 11399382).front(), found.front());
    ASSERT_EQ(get_flag(12, 11486677).front(), found.front());
    ASSERT_EQ(get_flag(12, 11486678).front(), found.front());
    ASSERT_EQ(get_flag(12, 11486679).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 119697658).front(), not_found.front());
    ASSERT_EQ(get_flag(12, 119697659).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697660).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697837).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697838).front(), found.front());
    ASSERT_EQ(get_flag(12, 119697839).front(), not_found.front());
}
TEST(REGION_MSIGDB, NAME_CROSS_CHR)
{
    // This should be ok? Transplicing or something like that?
    std::string gtf_name = path + "Test.gtf";
    std::string gmt_name = path + "Test.gmt";
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
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
    std::string gtf_name = path + "Test.gtf";
    std::string gmt_name = path + "Test.gmt";
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string msigdb = "";
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    try
    {
        Region::generate_regions(gene_sets, region_names, snp_in_sets, feature,
                                 window_5, window_3, genome_wide_background,
                                 gtf_name, msigdb, bed_names, snp_set,
                                 background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}
TEST(REGION_BACKGROUND, GTF_BACKGROUND)
{
    std::string gtf_name = path + "Test.gtf";
    std::string gmt_name = path + "Test.gmt";
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
           "\"ENSG00000122967\"; gene_version \"5\"; transcript_id "
           "\"ENST00000538349\"; transcript_version \"5\"; gene_name "
           "\"CCTV\"; gene_source \"havana\"; gene_biotype "
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    size_t num_regions;
    try
    {
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set,
            background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    // we should have the 2 sets, the base and the background
    ASSERT_EQ(num_regions, 4);
    const size_t required_size = BITCT_TO_WORDCT(num_regions);
    std::vector<uintptr_t> found = {0}, not_found = {0}, index = {0};
    SET_BIT(0, not_found.data());
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    // 12 11399381 11486678
    // SNP not in any location, not even in the background
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 1139, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    index.front() = 0;
    // SNP not in any region, but in background
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 11399381, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
}

TEST(REGION_BACKGROUND, GENOME_BACKGROUND)
{
    std::string gtf_name = path + "Test.gtf";
    std::string gmt_name = path + "Test.gmt";
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
           "\"ENSG00000122967\"; gene_version \"5\"; transcript_id "
           "\"ENST00000538349\"; transcript_version \"5\"; gene_name "
           "\"CCTV\"; gene_source \"havana\"; gene_biotype "
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
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = true;
    std::string snp_set = "";
    std::string background = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::vector<std::string> bed_names = {};
    size_t num_regions;
    try
    {
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set,
            background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    // we should have the 2 sets, the base and the background
    ASSERT_EQ(num_regions, 4);
    const size_t required_size = BITCT_TO_WORDCT(num_regions);
    std::vector<uintptr_t> found = {0}, not_found = {0}, index = {0};
    SET_BIT(0, not_found.data());
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    // 12 11399381 11486678
    // SNP not in any location will still be included in the background (genome
    // wide background)
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 1139, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    index[0] = 0;
    // SNP not in any region, but in background
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 11399381, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
}

TEST(REGION_BACKGROUND, BED_BACKGROUND)
{
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
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
    std::string background = path + "background";
    bed_file.open(background);
    bed_file << "1 89557 96038\n"
                "4  3016 87782\n"
                "10 14013 68802\n"
                "13 53964 90572\n"
                "14 22104 47572\n";
    bed_file.close();
    background.append(":bed");
    std::vector<std::string> bed_names = {bed_name};
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string snp_set = "", gtf_name = "", gmt_name = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    size_t num_regions;
    try
    {
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set,
            background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    // we should have the 3 sets, the base, the bed and the background
    ASSERT_EQ(num_regions, 3);
    const size_t required_size = BITCT_TO_WORDCT(num_regions);
    std::vector<uintptr_t> not_found = {0};
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, not_found.data());
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    // 4  3016 87782
    //"4 20139 97433 . . .\n"
    // 4  87000 should be found in both the set and the background
    // 13 53970 should only be found in the background
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 3015 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 3016 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 3017 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 20138 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    // found in both background and the bed file
    SET_BIT(2, found.data());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 20139 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 20140 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 87781 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    found.front() = 0;
    // only found in bed but not background
    SET_BIT(0, found.data());
    SET_BIT(2, found.data());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 87782 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 87783 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    found.front() = 0;
    // found in all
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             13, 53970 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    // anything on chromosome 17 should only be found in the base
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             17, 53970 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
}

TEST(REGION_BACKGROUND, RANGE_BACKGROUND)
{
    std::ofstream bed_file;
    std::string bed_name = path + "Test.bed";
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
    std::string background = path + "background";
    bed_file.open(background);
    bed_file << "1 89557 96038\n"
                "4  3016 87782\n"
                "10 14013 68802\n"
                "13 53964 90572\n"
                "14 22104 47572\n";
    bed_file.close();
    background.append(":range");
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    std::vector<std::string> bed_names = {bed_name};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string snp_set = "", gtf_name = "", gmt_name = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    size_t num_regions;
    try
    {
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set,
            background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    // we should have the 3 sets, the base, the bed and the background
    ASSERT_EQ(num_regions, 3);
    const size_t required_size = BITCT_TO_WORDCT(num_regions);
    std::vector<uintptr_t> not_found = {0};
    std::vector<uintptr_t> found = {0}, index = {0};
    SET_BIT(0, not_found.data());
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    // 4  3016 87782
    //"4 20139 97433 . . .\n"
    // 4  87000 should be found in both the set and the background
    // 13 53970 should only be found in the background
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 3015, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 3016, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 3017, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 20138 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    // also found in set
    SET_BIT(2, found.data());
    // the non-background range should be bed format
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 20139 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 20140 + 1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 87781, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 87782, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    // now only found in set
    found.front() = 0;
    SET_BIT(0, found.data());
    SET_BIT(2, found.data());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             4, 87783, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    found.front() = 0;
    // found in background only
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             13, 53970, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    // anything on chromosome 17 should only be found in the base
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             17, 53970, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
}

TEST(REGION_BACKGROUND, GENE_NAME_BACKGROUND)
{
    // This should be ok? Transplicing or something like that?
    std::string gtf_name = path + "Test.gtf";
    std::string gmt_name = path + "Test.gmt";
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
           "\"ENSG00000223974\"; "
           "gene_version \"5\"; gene_name \"DDX11L3\"; gene_source "
           "\"havana\"; "
           "gene_biotype \"transcribed_unprocessed_pseudogene\"; "
           "havana_gene "
           "\"OTTHUMG00000000961\"; havana_gene_version \"2\";\n"

           "12\thavana\tCDS\t11399381\t11486678\t.\t-\t.\tgene_id "
           "\"ENSG00000122967\"; gene_version \"5\"; transcript_id "
           "\"ENST00000538349\"; transcript_version \"5\"; gene_name "
           "\"CCTV\"; gene_source \"havana\"; gene_biotype "
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
    std::string background = path + "background";
    gmt.open(background.c_str());
    // accept both single line and multi-line (maybe a little bit too
    // flexible?
    gmt << "DDX11L1" << std::endl;
    gmt << "CIT CCTV" << std::endl;
    gmt.close();
    background.append(":gene");
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    std::vector<std::string> bed_names = {};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::string snp_set = "";
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    size_t num_regions;
    try
    {
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set,
            background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    // we should have the 4 sets, the base, the 2 regions and the background
    ASSERT_EQ(num_regions, 4);
    const size_t required_size = BITCT_TO_WORDCT(num_regions);
    std::vector<uintptr_t> found = {0}, not_found = {0}, index = {0};
    // 1 11869 14409 1110

    SET_BIT(0, not_found.data());
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(2, found.data());
    index.front() = 0;
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             1, 11868, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             1, 11869, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             1, 11870, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             1, 14408, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             1, 14409, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             1, 14410, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    // 2 15869 16409 1000
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             2, 15868, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             2, 15869, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             2, 15870, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             2, 16408, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             2, 16409, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             2, 16410, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    // 12 11399381 11486678 1100
    // only found in background
    found.front() = 0;
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 11399380, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 11399381, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 11399382, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 11486677, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 11486678, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 11486679, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    // 12 119697659 119697838 1101
    //
    found.front() = 0;
    SET_BIT(0, found.data());
    SET_BIT(1, found.data());
    SET_BIT(3, found.data());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 119697658, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 119697659, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 119697660, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 119697837, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 119697838, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("", gene_sets, snp_in_sets, index, required_size,
                             12, 119697839, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
}

TEST(REGION_SNP_SET, VERTICAL_SNP_SET)
{
    // This should be ok? Transplicing or something like that?
    std::string snp_set_name = path + "snp_set";
    std::ofstream snp_set;
    snp_set.open(snp_set_name.c_str());
    snp_set << "SNP_1\nSNP_2\nSNP_4\nSNP_5\n";
    snp_set.close();
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    std::vector<std::string> bed_names = {};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::string gtf_name = "", gmt_name = "", background = "";
    size_t num_regions;
    try
    {
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set_name,
            background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    // we should have the 3 sets, the base, the SNP_SET and the background
    ASSERT_EQ(num_regions, 3);
    // for single set, we use the file name as the set name
    ASSERT_STREQ(region_names[2].c_str(), snp_set_name.c_str());
    // Or we allow user defined name
    snp_set_name.append(":SNP_SET");
    region_names.clear();
    num_regions = Region::generate_regions(
        gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
        genome_wide_background, gtf_name, gmt_name, bed_names, snp_set_name,
        background, 22, reporter);
    ASSERT_EQ(num_regions, 3);
    ASSERT_STREQ(region_names[2].c_str(), "SNP_SET");
    ASSERT_EQ(snp_in_sets.size(), 4);
    const size_t required_size = BITCT_TO_WORDCT(num_regions);
    // we can simply check if the target SNPs are located in snp_in_sets
    // we have 1245
    std::vector<uintptr_t> found(required_size, 0), not_found(required_size, 0),
        index(required_size, 0);
    // here, we don't provide anything for background construction,
    // and as we set genome_wide_background as false, we will never
    // set the bit for background
    SET_BIT(0, found.data());
    SET_BIT(2, found.data());
    SET_BIT(0, not_found.data());
    Genotype::construct_flag("SNP_1", gene_sets, snp_in_sets, index,
                             required_size, -1, -1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("SNP_2", gene_sets, snp_in_sets, index,
                             required_size, -1, -1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("SNP_3", gene_sets, snp_in_sets, index,
                             required_size, -1, -1, genome_wide_background);
    ASSERT_EQ(index.front(), not_found.front());
    Genotype::construct_flag("SNP_4", gene_sets, snp_in_sets, index,
                             required_size, -1, -1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    Genotype::construct_flag("SNP_5", gene_sets, snp_in_sets, index,
                             required_size, -1, -1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
}

TEST(REGION_SNP_SET, MULTI_SNP_SET)
{
    // This should be ok? Transplicing or something like that?
    std::string snp_set_name = path + "snp_set";
    std::ofstream snp_set;
    snp_set.open(snp_set_name.c_str());
    snp_set << "SET_1 SNP_1 SNP_2 SNP_4 SNP_5\n";
    snp_set << "SET_2 SNP_12 SNP_8974 SNP_82 SNP_98\n";
    snp_set << "SET_3 www.google.com SNP_32 SNP_2 SNP_137 SNP_824\n";
    snp_set << "SET_4 SNP_86 SNP_478 SNP_155 SNP_743\n";
    snp_set << "SET_5 SNP_97 SNP_912 SNP_132 SNP_53\n";
    snp_set.close();
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    std::vector<std::string> bed_names = {};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::string gtf_name = "", gmt_name = "", background = "";
    size_t num_regions;
    try
    {
        // we don't want multi-set
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names,
            snp_set_name + ":SNP_SET", background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    ASSERT_STREQ(region_names[2].c_str(), "SET_1");
    ASSERT_STREQ(region_names[3].c_str(), "SET_2");
    ASSERT_STREQ(region_names[4].c_str(), "SET_3");
    ASSERT_STREQ(region_names[5].c_str(), "SET_4");
    ASSERT_STREQ(region_names[6].c_str(), "SET_5");
    try
    {
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set_name,
            background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    // we should have the 7 sets, the base, the SNP_SET and the background
    ASSERT_EQ(num_regions, 7);
    // For multi-set, we always use the first row as their name
    ASSERT_STREQ(region_names[2].c_str(), "SET_1");
    ASSERT_STREQ(region_names[3].c_str(), "SET_2");
    ASSERT_STREQ(region_names[4].c_str(), "SET_3");
    ASSERT_STREQ(region_names[5].c_str(), "SET_4");
    ASSERT_STREQ(region_names[6].c_str(), "SET_5");

    // now check inclusion
    const size_t required_size = BITCT_TO_WORDCT(num_regions);
    // we can simply check if the target SNPs are located in snp_in_sets
    // we have 1245
    std::vector<uintptr_t> found(required_size, 0), not_found(required_size, 0),
        index(required_size, 0);
    SET_BIT(0, found.data());
    SET_BIT(2, found.data());
    SET_BIT(4, found.data());
    SET_BIT(0, not_found.data());
    // SNP_2 1,3
    Genotype::construct_flag("SNP_2", gene_sets, snp_in_sets, index,
                             required_size, -1, -1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    // SNP_32 3
    found.front() = 0;
    SET_BIT(0, found.data());
    SET_BIT(4, found.data());
    Genotype::construct_flag("SNP_32", gene_sets, snp_in_sets, index,
                             required_size, -1, -1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
    // SNP_912 5
    found.front() = 0;
    SET_BIT(0, found.data());
    SET_BIT(6, found.data());
    Genotype::construct_flag("SNP_912", gene_sets, snp_in_sets, index,
                             required_size, -1, -1, genome_wide_background);
    ASSERT_EQ(index.front(), found.front());
}


TEST(REGION, UNINIT_INDEX)
{
    // This should be ok? Transplicing or something like that?
    std::string snp_set_name = path + "snp_set";
    std::ofstream snp_set;
    snp_set.open(snp_set_name.c_str());
    snp_set << "SET_1 SNP_1 SNP_2 SNP_4 SNP_5\n";
    snp_set << "SET_2 SNP_12 SNP_8974 SNP_82 SNP_98\n";
    snp_set << "SET_3 www.google.com SNP_32 SNP_2 SNP_137 SNP_824\n";
    snp_set << "SET_4 SNP_86 SNP_478 SNP_155 SNP_743\n";
    snp_set << "SET_5 SNP_97 SNP_912 SNP_132 SNP_53\n";
    snp_set.close();
    Reporter reporter(std::string(path + "LOG"));
    std::vector<std::string> feature = {"exon", "gene", "protein_coding",
                                        "CDS"};
    std::vector<std::string> bed_names = {};
    int window_5 = 0;
    int window_3 = 0;
    bool genome_wide_background = false;
    std::vector<std::string> region_names;
    std::unordered_map<std::string, std::vector<int>> snp_in_sets;
    std::vector<IITree<int, int>> gene_sets;
    std::string gtf_name = "", gmt_name = "", background = "";
    size_t num_regions;
    try
    {
        // we don't want multi-set
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names,
            snp_set_name + ":SNP_SET", background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    ASSERT_STREQ(region_names[2].c_str(), "SET_1");
    ASSERT_STREQ(region_names[3].c_str(), "SET_2");
    ASSERT_STREQ(region_names[4].c_str(), "SET_3");
    ASSERT_STREQ(region_names[5].c_str(), "SET_4");
    ASSERT_STREQ(region_names[6].c_str(), "SET_5");
    try
    {
        num_regions = Region::generate_regions(
            gene_sets, region_names, snp_in_sets, feature, window_5, window_3,
            genome_wide_background, gtf_name, gmt_name, bed_names, snp_set_name,
            background, 22, reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
    // we should have the 7 sets, the base, the SNP_SET and the background
    ASSERT_EQ(num_regions, 7);
    // For multi-set, we always use the first row as their name
    ASSERT_STREQ(region_names[2].c_str(), "SET_1");
    ASSERT_STREQ(region_names[3].c_str(), "SET_2");
    ASSERT_STREQ(region_names[4].c_str(), "SET_3");
    ASSERT_STREQ(region_names[5].c_str(), "SET_4");
    ASSERT_STREQ(region_names[6].c_str(), "SET_5");

    // now check inclusion
    const size_t required_size = BITCT_TO_WORDCT(num_regions);
    // we can simply check if the target SNPs are located in snp_in_sets
    // we have 1245
    // should still work without initializing the index
    std::vector<uintptr_t> index;
    Genotype::construct_flag("SNP_2", gene_sets, snp_in_sets, index,
                             required_size, -1, -1, genome_wide_background);
}
#endif // REGION_TEST_HPP
