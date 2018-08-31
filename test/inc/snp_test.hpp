#ifndef SNP_TEST_HPP
#define SNP_TEST_HPP
#include "snp.hpp"
#include "gtest/gtest.h"
#include <vector>

TEST(SNP_TEST, INITIALIZE_NO_COUNT)
{
    // check if the initialization sets all the parameters correctly
    uint32_t homcom, het, homrar, missing;
    SNP snp("Test", 1, 1, "A", "C", "Input", 1);
    ASSERT_STREQ(snp.rs().c_str(), "Test");
    ASSERT_STREQ(snp.ref().c_str(), "A");
    ASSERT_STREQ(snp.alt().c_str(), "C");
    ASSERT_EQ(snp.chr(), 1);
    ASSERT_EQ(snp.loc(), 1);
    // We want the target and reference to be the same unless set_ref is used
    ASSERT_STREQ(snp.file_name().c_str(), "Input");
    ASSERT_EQ(snp.byte_pos(), 1);
    ASSERT_STREQ(snp.ref_file_name().c_str(), "Input");
    ASSERT_EQ(snp.ref_byte_pos(), 1);
    // When initialize without count, has_count (return value of get_counts)
    // should be false
    ASSERT_FALSE(snp.get_counts(homcom, het, homrar, missing));
    // default of clump should be false
    ASSERT_FALSE(snp.clumped());
    // default of flipped should be false
    ASSERT_FALSE(snp.is_flipped());
    // default statistic is 0.0
    ASSERT_DOUBLE_EQ(snp.stat(), 0.0);
    // default p-value is 2.0, therefore always the least significant
    ASSERT_DOUBLE_EQ(snp.p_value(), 2.0);
    // default threshold should be 0.0
    ASSERT_DOUBLE_EQ(snp.get_threshold(), 0.0);
    // default bounaries is always 0
    ASSERT_EQ(snp.low_bound(), 0);
    ASSERT_EQ(snp.up_bound(), 0);
    // default category is -1,
    ASSERT_EQ(snp.category(), -1);
    ASSERT_EQ(homcom, 0);
    ASSERT_EQ(het, 0);
    ASSERT_EQ(homrar, 0);
    ASSERT_EQ(missing, 0);
}
TEST(SNP_TEST, INITIALIZE_COUNT)
{
    // check if the initialization sets all the parameters correctly
    uint32_t homcom, het, homrar, missing;
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    ASSERT_STREQ(snp.rs().c_str(), "Test");
    ASSERT_STREQ(snp.ref().c_str(), "A");
    ASSERT_STREQ(snp.alt().c_str(), "C");
    ASSERT_EQ(snp.chr(), 1);
    ASSERT_EQ(snp.loc(), 1);
    // We want the target and reference to be the same unless set_ref is used
    ASSERT_STREQ(snp.file_name().c_str(), "Input");
    ASSERT_EQ(snp.byte_pos(), 1);
    ASSERT_STREQ(snp.ref_file_name().c_str(), "Input");
    ASSERT_EQ(snp.ref_byte_pos(), 1);
    // When initialize with count, has_count (return value of get_counts)
    // should be true
    ASSERT_TRUE(snp.get_counts(homcom, het, homrar, missing));
    // default of clump should be false
    ASSERT_FALSE(snp.clumped());
    // default of flipped should be false
    ASSERT_FALSE(snp.is_flipped());
    // default statistic is 0.0
    ASSERT_DOUBLE_EQ(snp.stat(), 0.0);
    // default p-value is 2.0, therefore always the least significant
    ASSERT_DOUBLE_EQ(snp.p_value(), 2.0);
    // default threshold should be 0.0
    ASSERT_DOUBLE_EQ(snp.get_threshold(), 0.0);
    // default category is -1,
    ASSERT_EQ(snp.category(), -1);
    // default bounaries is always 0
    ASSERT_EQ(snp.low_bound(), 0);
    ASSERT_EQ(snp.up_bound(), 0);
    ASSERT_EQ(homcom, 1);
    ASSERT_EQ(het, 2);
    ASSERT_EQ(homrar, 3);
    ASSERT_EQ(missing, 4);
}
TEST(SNP_TEST, SET_STATISTIC)
{
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    // default statistic is 0.0
    ASSERT_DOUBLE_EQ(snp.stat(), 0.0);
    // default p-value is 2.0, therefore always the least significant
    ASSERT_DOUBLE_EQ(snp.p_value(), 2.0);
    // default threshold should be 0.0
    ASSERT_DOUBLE_EQ(snp.get_threshold(), 0.0);
    // default category is -1,
    ASSERT_EQ(snp.category(), -1);
    // setting statistic
    snp.set_statistic(0.23498, 0.05, 1, 0.05);
    ASSERT_DOUBLE_EQ(snp.stat(), 0.23498);
    ASSERT_DOUBLE_EQ(snp.p_value(), 0.05);
    ASSERT_DOUBLE_EQ(snp.get_threshold(), 0.05);
    ASSERT_EQ(snp.category(), 1);
}
// Might want a test to test if set_statistic throw the correct assertion error
// when category == -1
TEST(SNP_TEST, INVALIDATE)
{
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    // default is valid
    ASSERT_TRUE(snp.valid());
    // we then invalidate it
    snp.invalidate();
    ASSERT_FALSE(snp.valid());
}
TEST(SNP_TEST, ADD_REF)
{
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    // default, reference and the target are the same unless set_ref is used
    ASSERT_STREQ(snp.file_name().c_str(), "Input");
    ASSERT_EQ(snp.byte_pos(), 1);
    ASSERT_STREQ(snp.ref_file_name().c_str(), "Input");
    ASSERT_EQ(snp.ref_byte_pos(), 1);
    snp.add_reference("Reference", 1);
    ASSERT_STREQ(snp.ref_file_name().c_str(), "Reference");
    ASSERT_EQ(snp.ref_byte_pos(), 1);
    snp.add_reference("Reference", 13789560123);
    ASSERT_EQ(snp.ref_byte_pos(), 13789560123);
}
TEST(SNP_MATCHING, FLIPPING_AC)
{
    // Flipping occurrs
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    bool flipped = false;
    std::string ref = "A", alt = "C";
    ASSERT_TRUE(snp.matching(1, 1, alt, ref, flipped));
    // the flipped boolean should change to true
    ASSERT_TRUE(flipped);
    flipped = false;
    // we should get the same result if we have complements
    ref = "T";
    alt = "G";
    ASSERT_TRUE(snp.matching(1, 1, alt, ref, flipped));
    // the flipped boolean should change to true
    ASSERT_TRUE(flipped);
}
TEST(SNP_MATCHING, FLIPPING_GT)
{
    // Flipping occurrs
    SNP snp("Test", 1, 1, "G", "T", "Input", 1, 1, 2, 3, 4);
    bool flipped = false;
    std::string ref = "G", alt = "T";
    ASSERT_TRUE(snp.matching(1, 1, alt, ref, flipped));
    // the flipped boolean should change to true
    ASSERT_TRUE(flipped);
    flipped = false;
    // we should get the same result if we have complements
    ref = "C";
    alt = "A";
    ASSERT_TRUE(snp.matching(1, 1, alt, ref, flipped));
    // the flipped boolean should change to true
    ASSERT_TRUE(flipped);
}
TEST(SNP_MATCHING, NO_ALT_AC)
{
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    bool flipped = false;
    std::string ref = "A", alt = "";
    // should still return true
    ASSERT_TRUE(snp.matching(1, 1, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    // We won't do flipping if we don't have the information of the
    // alternative allele and will assume this is a mismatch
    ref = "C";
    flipped = false;
    ASSERT_FALSE(snp.matching(1, 1, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    // Same for the complement
    ref = "G";
    flipped = false;
    ASSERT_FALSE(snp.matching(1, 1, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    // but we will assume true if this is a complement
    ref = "T";
    flipped = false;
    ASSERT_TRUE(snp.matching(1, 1, ref, alt, flipped));
    ASSERT_FALSE(flipped);
}
TEST(SNP_MATCHING, NO_ALT_GT)
{
    SNP snp("Test", 1, 1, "G", "T", "Input", 1, 1, 2, 3, 4);
    bool flipped = false;
    std::string ref = "G", alt = "";
    // should still return true
    ASSERT_TRUE(snp.matching(1, 1, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    // We won't do flipping if we don't have the information of the
    // alternative allele and will consider it as mismatch
    ref = "T";
    flipped = false;
    ASSERT_FALSE(snp.matching(1, 1, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    // Same for the complement
    ref = "A";
    flipped = false;
    ASSERT_FALSE(snp.matching(1, 1, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    // but we will assume true if this is a complement
    ref = "C";
    flipped = false;
    ASSERT_TRUE(snp.matching(1, 1, ref, alt, flipped));
    ASSERT_FALSE(flipped);
}
TEST(SNP_MATCHING, CHR_POS_MATCHING)
{
    // we assume we always got the chr and loc information when we
    // initialize our SNP object (as we are either reading from bim or from
    // bgen which contain those information as part of the file
    // specification Matching has an assumption that the alleles should be
    // in the same case
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    bool flipped = false;
    std::string ref = "A", alt = "C";
    // chromosome mismatch
    ASSERT_FALSE(snp.matching(2, 1, ref, alt, flipped));
    // the flipped boolean should remain the same
    ASSERT_FALSE(flipped);
    flipped = false;
    // base pair mismatch
    ASSERT_FALSE(snp.matching(1, 2, ref, alt, flipped));
    // the flipped boolean should remain the same
    ASSERT_FALSE(flipped);
    flipped = false;
    // No chromosome information mismatch
    ASSERT_TRUE(snp.matching(-1, 1, ref, alt, flipped));
    // the flipped boolean should remain the same
    ASSERT_FALSE(flipped);
    flipped = false;
    // No BP information mismatch
    ASSERT_TRUE(snp.matching(1, -1, ref, alt, flipped));
    // the flipped boolean should remain the same
    ASSERT_FALSE(flipped);
    flipped = false;
    // No BP & CHR information mismatch
    ASSERT_TRUE(snp.matching(-1, -1, ref, alt, flipped));
    // the flipped boolean should remain the same
    ASSERT_FALSE(flipped);
}
TEST(SNP_MATCHING, SET_FLIPPED)
{
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    ASSERT_FALSE(snp.is_flipped());
    snp.set_flipped();
    ASSERT_TRUE(snp.is_flipped());
}
TEST(SNP_CLUMP, SET_CLUMP)
{
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    // default of clump should be false
    ASSERT_FALSE(snp.clumped());
    snp.set_clumped();
    // set clumped should set clump to true
    ASSERT_TRUE(snp.clumped());
}
TEST(SNP_TEST, SET_LOW_BOUND)
{
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    ASSERT_EQ(snp.low_bound(), 0);
    ASSERT_EQ(snp.up_bound(), 0);
    snp.set_low_bound(10);
    ASSERT_EQ(snp.low_bound(), 10);
    ASSERT_EQ(snp.up_bound(), 0);
}
TEST(SNP_TEST, SET_UP_BOUND)
{
    SNP snp("Test", 1, 1, "A", "C", "Input", 1, 1, 2, 3, 4);
    ASSERT_EQ(snp.low_bound(), 0);
    ASSERT_EQ(snp.low_bound(), 0);
    snp.set_up_bound(10);
    ASSERT_EQ(snp.up_bound(), 10);
    ASSERT_EQ(snp.low_bound(), 0);
}
TEST(SNP_TEST, SORT_BY_P_CHR)
{
    std::vector<SNP> snps;
    // first generate SNPs for sorting
    // we need to test the following
    // 1. Same chromosome
    //    a. Different p-value    : SNP_A SNP_C
    //       i. Same everything (except name)
    //    b. Same p-value
    //       i. different location : SNP_B SNP_D
    //      ii. same location:  SNP_C SNP_E
    //          - Different name
    // 2. Different chromosome
    //    a. Same everything (except name) : SNP_A SNP_B

    // by definition, we don't allow multiple SNPs with the same name. Using SNP
    // name as the last comparison condition should allow us to avoid troubles


    snps.emplace_back(SNP("SNP_A", 1, 10, "A", "", "Input", 1));
    snps.back().set_statistic(1, 0.05, 1, 0.05);
    snps.emplace_back(SNP("SNP_B", 2, 10, "A", "", "Input", 1));
    snps.back().set_statistic(1, 0.05, 1, 0.05);
    snps.emplace_back(SNP("SNP_C", 1, 10, "A", "", "Input", 1));
    snps.back().set_statistic(1, 0.01, 1, 0.05);
    snps.emplace_back(SNP("SNP_D", 2, 11, "A", "", "Input", 1));
    snps.back().set_statistic(1, 0.05, 1, 0.05);
    snps.emplace_back(SNP("SNP_E", 1, 10, "A", "", "Input", 1));
    snps.back().set_statistic(1, 0.01, 1, 0.05);

    // our desired order for the above input should be
    // 2,4,0,1,3
    std::vector<size_t> result = SNP::sort_by_p_chr(snps);
    ASSERT_EQ(result.size(), snps.size());
    ASSERT_EQ(result[0], 2);
    ASSERT_EQ(result[1], 4);
    ASSERT_EQ(result[2], 0);
    ASSERT_EQ(result[3], 1);
    ASSERT_EQ(result[4], 3);
}
// clump will need a special test just for that (Need to take into account of
// both region and SNP)
#endif // SNP_TEST_HPP
