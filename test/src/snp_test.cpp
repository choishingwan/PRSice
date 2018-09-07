#ifndef SNP_TEST_HPP
#define SNP_TEST_HPP
#include "genotype.hpp"
#include "global.hpp"
#include "region.hpp"
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
    snp.set_statistic(0.23498, 0.05, 0.123, 0.06, 1, 0.05);
    ASSERT_DOUBLE_EQ(snp.stat(), 0.23498);
    ASSERT_DOUBLE_EQ(snp.p_value(), 0.05);
    ASSERT_DOUBLE_EQ(snp.get_threshold(), 0.05);
    ASSERT_DOUBLE_EQ(snp.get_se(), 0.123);
    ASSERT_DOUBLE_EQ(snp.get_maf(), 0.06);
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
    snps.back().set_statistic(1, 0.05, 0.123, 0.5, 1, 0.05);
    snps.emplace_back(SNP("SNP_B", 2, 10, "A", "", "Input", 1));
    snps.back().set_statistic(1, 0.05, 0.123, 0.5, 1, 0.05);
    snps.emplace_back(SNP("SNP_C", 1, 10, "A", "", "Input", 1));
    snps.back().set_statistic(1, 0.01, 0.123, 0.5, 1, 0.05);
    snps.emplace_back(SNP("SNP_D", 2, 11, "A", "", "Input", 1));
    snps.back().set_statistic(1, 0.05, 0.123, 0.5, 1, 0.05);
    snps.emplace_back(SNP("SNP_E", 1, 10, "A", "", "Input", 1));
    snps.back().set_statistic(1, 0.01, 0.123, 0.5, 1, 0.05);

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
class GenotypeSNPTest : public Genotype
{
public:
    GenotypeSNPTest() { m_max_code = 22; }
};
// Any error in the GTF file will lead to throw
class SNP_REGION : public ::testing::Test
{
    // For exclusion, strand information should not alter result (window
    // padding should all be 0)
protected:
    Region region;
    std::vector<uintptr_t> not_found = {0};
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
        region = Region(feature, 0, 0, false, false);
        std::vector<std::string> bed_names = {};
        GenotypeSNPTest dummy;
        region.generate_regions(gtf_name, gmt_name, bed_names, "", "", "",
                                dummy, reporter);
        SET_BIT(0, not_found.data());
    }
};

TEST_F(SNP_REGION, BASE_SET1_STANDARD)
{
    SNP base_snp("Base_SNP", 1, 1, "A", "C", "Test", 1);
    SNP set_snp("Set_SNP", 1, 11869, "A", "C", "Test", 1);
    base_snp.set_flag(region);
    set_snp.set_flag(region);
    ASSERT_FALSE(base_snp.clumped());
    ASSERT_FALSE(set_snp.clumped());

    ASSERT_TRUE(base_snp.in(0));
    for (size_t i = 1; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        ASSERT_FALSE(base_snp.in(i));
    }

    for (size_t i = 0; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        if (i == 0 || i == 1 || i == 6) {
            ASSERT_TRUE(set_snp.in(i));
        }
        else
        {
            ASSERT_FALSE(set_snp.in(i));
        }
    }
    base_snp.clump(set_snp, 0.1, false, 0.8);
    // will set clumped to true to protect itself
    ASSERT_TRUE(base_snp.clumped());
    ASSERT_FALSE(set_snp.clumped());
    ASSERT_TRUE(base_snp.in(0));
    for (size_t i = 1; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        ASSERT_FALSE(base_snp.in(i));
    }
    // the set SNP should've lost the base flag
    for (size_t i = 0; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        if (i == 1 || i == 6) {
            ASSERT_TRUE(set_snp.in(i));
        }
        else
        {
            ASSERT_FALSE(set_snp.in(i));
        }
    }
}
TEST_F(SNP_REGION, BASE_SET1_PROXY_NO_GO)
{
    SNP base_snp("Base_SNP", 1, 1, "A", "C", "Test", 1);
    SNP set_snp("Set_SNP", 1, 11869, "A", "C", "Test", 1);
    base_snp.set_flag(region);
    set_snp.set_flag(region);
    ASSERT_FALSE(base_snp.clumped());
    ASSERT_FALSE(set_snp.clumped());
    ASSERT_TRUE(base_snp.in(0));
    for (size_t i = 1; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        ASSERT_FALSE(base_snp.in(i));
    }

    for (size_t i = 0; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        if (i == 0 || i == 1 || i == 6) {
            ASSERT_TRUE(set_snp.in(i));
        }
        else
        {
            ASSERT_FALSE(set_snp.in(i));
        }
    }
    base_snp.clump(set_snp, 0.1, true, 0.8);
    // will set clumped to true to protect itself
    ASSERT_TRUE(base_snp.clumped());
    ASSERT_FALSE(set_snp.clumped());
    ASSERT_TRUE(base_snp.in(0));
    for (size_t i = 1; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        ASSERT_FALSE(base_snp.in(i));
    }
    // the set SNP should've lost the base flag
    for (size_t i = 0; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        if (i == 1 || i == 6) {
            ASSERT_TRUE(set_snp.in(i));
        }
        else
        {
            ASSERT_FALSE(set_snp.in(i));
        }
    }
}

TEST_F(SNP_REGION, BASE_SET1_PROXY_GO)
{
    SNP base_snp("Base_SNP", 1, 1, "A", "C", "Test", 1);
    SNP set_snp("Set_SNP", 1, 11869, "A", "C", "Test", 1);
    base_snp.set_flag(region);
    set_snp.set_flag(region);
    ASSERT_FALSE(base_snp.clumped());
    ASSERT_FALSE(set_snp.clumped());

    ASSERT_TRUE(base_snp.in(0));
    for (size_t i = 1; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        ASSERT_FALSE(base_snp.in(i));
    }

    for (size_t i = 0; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        if (i == 0 || i == 1 || i == 6) {
            ASSERT_TRUE(set_snp.in(i));
        }
        else
        {
            ASSERT_FALSE(set_snp.in(i));
        }
    }
    base_snp.clump(set_snp, 0.1, true, 0.01);
    // will set clumped to true to protect itself
    ASSERT_TRUE(base_snp.clumped());
    ASSERT_TRUE(set_snp.clumped());
    for (size_t i = 0; i < 7; ++i) {
        // we should have acquire the flag from set_snp
        if (i == 0 || i == 1 || i == 6) {
            ASSERT_TRUE(base_snp.in(i));
        }
        else
        {
            ASSERT_FALSE(base_snp.in(i));
        }
    }
    // we simply set the set_snp to clumped without setting its flags to zero
}
TEST_F(SNP_REGION, BASE_BASE_STANDARD)
{
    SNP base_snp("Base_SNP", 1, 1, "A", "C", "Test", 1);
    SNP set_snp("Set_SNP", 1, 2, "A", "C", "Test", 1);
    base_snp.set_flag(region);
    set_snp.set_flag(region);
    ASSERT_FALSE(base_snp.clumped());
    ASSERT_FALSE(set_snp.clumped());

    ASSERT_TRUE(base_snp.in(0));
    ASSERT_TRUE(set_snp.in(0));
    for (size_t i = 1; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        ASSERT_FALSE(base_snp.in(i));
        ASSERT_FALSE(set_snp.in(i));
    }

    base_snp.clump(set_snp, 0.1, false, 0.8);
    // will set clumped to true to protect itself
    ASSERT_TRUE(base_snp.clumped());
    ASSERT_TRUE(set_snp.clumped());
    ASSERT_TRUE(base_snp.in(0));
    for (size_t i = 1; i < 7; ++i) {
        // we should not acquire the flag from set_snp
        ASSERT_FALSE(base_snp.in(i));
    }
    // the set SNP should've lost the base flag
    for (size_t i = 0; i < 7; ++i) {
        ASSERT_FALSE(set_snp.in(i));
    }
}
#endif // SNP_TEST_HPP
