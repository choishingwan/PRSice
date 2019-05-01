#ifndef SNP_TEST_HPP
#define SNP_TEST_HPP
#include "genotype.hpp"
#include "global.hpp"
#include "region.hpp"
#include "snp.hpp"
#include "gtest/gtest.h"
#include <vector>
class SNP_INIT_TEST : public ::testing::Test
{
protected:
    SNP snp;
    std::string rs = "Test";
    std::string ref = "A";
    std::string alt = "C";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = 1, loc = 1;
    void SetUp() override
    {
        snp=SNP(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
    }
    void TearDown() override {}


};

TEST_F(SNP_INIT_TEST, INIT_TEST)
{
    // check if the initialization sets all the parameters correctly
    size_t homcom, het, homrar, missing;
    ASSERT_STREQ(snp.rs().c_str(), rs.c_str());
    ASSERT_STREQ(snp.ref().c_str(), ref.c_str());
    ASSERT_STREQ(snp.alt().c_str(), alt.c_str());
    ASSERT_EQ(snp.chr(), chr);
    ASSERT_EQ(snp.loc(), loc);
    // The file information should be missing thus far
    ASSERT_TRUE(snp.file_name().empty());
    ASSERT_EQ(snp.byte_pos(), 0);
    ASSERT_TRUE(snp.ref_file_name().empty());
    ASSERT_EQ(snp.ref_byte_pos(), 0);
    // When initialize without count, has_count (return value of get_counts)
    // should be false
    ASSERT_FALSE(snp.get_counts(homcom, het, homrar, missing, false));
    // should be false for both using reference maf and not using it
    ASSERT_FALSE(snp.get_counts(homcom, het, homrar, missing, true));
    // default of clump should be false
    ASSERT_FALSE(snp.clumped());
    // default of flipped should be false
    ASSERT_FALSE(snp.is_flipped());
    ASSERT_DOUBLE_EQ(snp.stat(), stat);
    ASSERT_DOUBLE_EQ(snp.p_value(), p);
    ASSERT_DOUBLE_EQ(snp.get_threshold(), p_threshold);
    ASSERT_EQ(snp.category(), category);
    // default bounaries is always 0
    ASSERT_EQ(snp.low_bound(), 0);
    ASSERT_EQ(snp.up_bound(), 0);
    ASSERT_EQ(homcom, 0);
    ASSERT_EQ(het, 0);
    ASSERT_EQ(homrar, 0);
    ASSERT_EQ(missing, 0);
}

TEST_F(SNP_INIT_TEST, ADD_REF)
{
    // default, reference are empty
    ASSERT_TRUE(snp.file_name().empty());
    // and the bytepos is 0
    ASSERT_EQ(snp.byte_pos(), 0);
    ASSERT_TRUE(snp.ref_file_name().empty());
    ASSERT_EQ(snp.ref_byte_pos(), 0);
    // we also need to know if we are flipping
    snp.add_reference("Reference", 1, true);
    ASSERT_STREQ(snp.ref_file_name().c_str(), "Reference");
    ASSERT_EQ(snp.ref_byte_pos(), 1);
    // should not touch target's flip flag
    ASSERT_FALSE(snp.is_flipped());
    ASSERT_TRUE(snp.is_ref_flipped());
    snp.add_reference("Reference", 1, false);
    ASSERT_STREQ(snp.ref_file_name().c_str(), "Reference");
    ASSERT_EQ(snp.ref_byte_pos(), 1);
    ASSERT_FALSE(snp.is_flipped());
    ASSERT_FALSE(snp.is_ref_flipped());
    snp.add_reference("Reference", 13789560123, true);
    ASSERT_EQ(snp.ref_byte_pos(), 13789560123);
}


TEST_F(SNP_INIT_TEST, ADD_TARGET)
{
    // default, target are empty
    ASSERT_TRUE(snp.file_name().empty());
    // and the bytepos is 0
    ASSERT_EQ(snp.byte_pos(), 0);
    ASSERT_TRUE(snp.ref_file_name().empty());
    ASSERT_EQ(snp.ref_byte_pos(), 0);
    // we also need to know if we are flipping
    int new_chr = 2;
    int new_loc = 3;
    std::string new_ref = "C";
    std::string new_alt = "G";
    std::string target_name = "Target";
    std::streampos new_pos = 1;
    snp.add_target(target_name, new_pos, new_chr, new_loc, new_ref, new_alt, true);
    // check if the names are updated correctly
    ASSERT_STREQ(snp.file_name().c_str(), target_name.c_str());
    ASSERT_EQ(snp.byte_pos(), new_pos);
    // Reference should follow target's name (always do target before ref)
    ASSERT_STREQ(snp.ref_file_name().c_str(), target_name.c_str());
    ASSERT_EQ(snp.ref_byte_pos(), new_pos);
    ASSERT_EQ(snp.chr(), new_chr);
    ASSERT_EQ(snp.loc(), new_loc);
    ASSERT_STREQ(snp.ref().c_str(), new_ref.c_str());
    ASSERT_STREQ(snp.alt().c_str(), new_alt.c_str());
    ASSERT_TRUE(snp.is_flipped());
    ASSERT_FALSE(snp.is_ref_flipped());
    // check none-flip
    snp.add_target(target_name, new_pos, new_chr, new_loc, new_ref, new_alt, false);
    ASSERT_FALSE(snp.is_flipped());
    ASSERT_FALSE(snp.is_ref_flipped());
    new_loc = 189560123;
    snp.add_target(target_name, new_pos, new_chr, new_loc, new_ref, new_alt, false);
    ASSERT_EQ(snp.byte_pos(), new_pos);
}

TEST_F(SNP_INIT_TEST, TARGET_AND_REF)
{
    // default, target are empty
    ASSERT_TRUE(snp.file_name().empty());
    // and the bytepos is 0
    ASSERT_EQ(snp.byte_pos(), 0);
    ASSERT_TRUE(snp.ref_file_name().empty());
    ASSERT_EQ(snp.ref_byte_pos(), 0);
    // we also need to know if we are flipping
    int new_chr = 2;
    int new_loc = 3;
    std::string new_ref = "C";
    std::string new_alt = "G";
    std::string target_name = "Target";
    std::streampos new_pos = 1;
    std::string ref_name = "reference";
    std::streampos new_ref_pos = 19;
    snp.add_target(target_name, new_pos, new_chr, new_loc, new_ref, new_alt, true);
    // check if the names are updated correctly
    ASSERT_STREQ(snp.file_name().c_str(), target_name.c_str());
    ASSERT_EQ(snp.byte_pos(), new_pos);
    // Reference should follow target's name (always do target before ref)
    ASSERT_STREQ(snp.ref_file_name().c_str(), target_name.c_str());
    ASSERT_EQ(snp.ref_byte_pos(), new_pos);
    ASSERT_EQ(snp.chr(), new_chr);
    ASSERT_EQ(snp.loc(), new_loc);
    ASSERT_STREQ(snp.ref().c_str(), new_ref.c_str());
    ASSERT_STREQ(snp.alt().c_str(), new_alt.c_str());
    ASSERT_TRUE(snp.is_flipped());
    ASSERT_FALSE(snp.is_ref_flipped());
    // check none-flip
    snp.add_reference(ref_name, new_ref_pos, false);
    ASSERT_TRUE(snp.is_flipped());
    ASSERT_FALSE(snp.is_ref_flipped());
    ASSERT_STREQ(snp.file_name().c_str(), target_name.c_str());
    ASSERT_EQ(snp.byte_pos(), new_pos);
    // Reference should follow target's name (always do target before ref)
    ASSERT_STREQ(snp.ref_file_name().c_str(), ref_name.c_str());
    ASSERT_EQ(snp.ref_byte_pos(), new_ref_pos);
}

TEST(SNP_MATCHING, FLIPPING_AC)
{
    std::string rs = "Test";
    std::string ref = "A";
    std::string alt = "C";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = 1, loc = 1;
    SNP snp(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
    // Flipping occurrs
    bool flipped = false;
    ASSERT_TRUE(snp.matching(chr, loc, alt, ref, flipped));
    // the flipped boolean should change to true
    ASSERT_TRUE(flipped);
    flipped = false;
    ASSERT_TRUE(snp.matching(chr, loc, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    flipped = false;
    // we should get the same result if we have complements
    ref = "T";
    alt = "G";
    ASSERT_TRUE(snp.matching(chr, loc, alt, ref, flipped));
    // the flipped boolean should change to true
    ASSERT_TRUE(flipped);
    flipped = false;
    ASSERT_TRUE(snp.matching(chr, loc, ref, alt, flipped));
    ASSERT_FALSE(flipped);
}
TEST(SNP_MATCHING, FLIPPING_GT)
{
    // Flipping occurrs
    bool flipped = false;
    std::string ref = "G", alt = "T";
    std::string rs = "Test";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = 1, loc = 1;
    SNP snp(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
    // change referenace and alt using add_target function
    ASSERT_TRUE(snp.matching(chr, loc, alt, ref, flipped));
    // the flipped boolean should change to true
    ASSERT_TRUE(flipped);
    flipped = false;
    ASSERT_TRUE(snp.matching(chr, loc, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    flipped = false;
    // we should get the same result if we have complements
    ref = "C";
    alt = "A";
    ASSERT_TRUE(snp.matching(chr, loc, alt, ref, flipped));
    // the flipped boolean should change to true
    ASSERT_TRUE(flipped);
    flipped = false;
    ASSERT_TRUE(snp.matching(chr, loc, ref, alt, flipped));
    ASSERT_FALSE(flipped);
}

TEST(SNP_MATCHING, NO_ALT_AC)
{
    std::string rs = "Test";
    std::string ref = "A";
    std::string alt = "";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = 1, loc = 1;
    SNP snp(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
    bool flipped = false;
    alt = "C";
    // should still return true
    ASSERT_TRUE(snp.matching(chr, loc, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    flipped = false;
    // as the original alt is "" and ref is A
    // and the input is ref C alt A we will consider this mismatch
    // as we are not comfortable in flipping it
    ASSERT_FALSE(snp.matching(chr, loc, alt, ref, flipped));
    ref = "T";
    // now we have A"" and TC and we will allow this to match
    flipped = false;
    ASSERT_TRUE(snp.matching(chr, loc, ref, alt, flipped));
    ASSERT_FALSE(flipped);
    // but again, we don't allow A"" and CT to match
    flipped = false;
    ASSERT_FALSE(snp.matching(chr, loc, alt, ref, flipped));
    ASSERT_FALSE(flipped);
    // similarly, we don't allow A"" and GT to match
    ref="G";
    alt="T";
    flipped = false;
    ASSERT_FALSE(snp.matching(chr, loc, ref, alt, flipped));
    ASSERT_FALSE(flipped);

}

TEST(SNP_MATCHING, CHR_POS_MATCHING)
{
    // we assume we always got the chr and loc information when we
    // initialize our SNP object (as we are either reading from bim or from
    // bgen which contain those information as part of the file
    // specification Matching has an assumption that the alleles should be
    // in the same case
    std::string rs = "Test";
    std::string ref = "G";
    std::string alt = "";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = 1, loc = 1;
    SNP snp(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
    bool flipped = false;
    // chromosome mismatch
    ASSERT_FALSE(snp.matching(chr+1, loc, ref, alt, flipped));
    // the flipped boolean should remain the same
    ASSERT_FALSE(flipped);
    flipped = false;
    // base pair mismatch
    ASSERT_FALSE(snp.matching(chr, loc+1, ref, alt, flipped));
    // the flipped boolean should remain the same
    ASSERT_FALSE(flipped);
    flipped = false;
}

TEST(SNP_MATCHING, NO_CHR_MATCHING)
{
    // we assume we always got the chr and loc information when we
    // initialize our SNP object (as we are either reading from bim or from
    // bgen which contain those information as part of the file
    // specification Matching has an assumption that the alleles should be
    // in the same case
    std::string rs = "Test";
    std::string ref = "G";
    std::string alt = "";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = -1, loc = 1;
    SNP snp(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
    bool flipped = false;
    // When chr is -1, we don't care if it is different as we assume it is missing
    ASSERT_TRUE(snp.matching(chr+10, loc, ref, alt, flipped));
    // the flipped boolean should remain the same
    ASSERT_FALSE(flipped);
    flipped = false;
}
TEST(SNP_MATCHING, NO_BP_MATCHING)
{
    // we assume we always got the chr and loc information when we
    // initialize our SNP object (as we are either reading from bim or from
    // bgen which contain those information as part of the file
    // specification Matching has an assumption that the alleles should be
    // in the same case
    std::string rs = "Test";
    std::string ref = "G";
    std::string alt = "";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = 1, loc = -1;
    SNP snp(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
    bool flipped = false;
    // When bp is -1, we don't care if it is different as we assume it is missing
    ASSERT_TRUE(snp.matching(chr, loc+10, ref, alt, flipped));
    // the flipped boolean should remain the same
    ASSERT_FALSE(flipped);
    flipped = false;
}

TEST(SNP_CLUMP, SET_CLUMP)
{
    std::string rs = "Test";
    std::string ref = "G";
    std::string alt = "";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = 1, loc = -1;
    SNP snp(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
    // default of clump should be false
    ASSERT_FALSE(snp.clumped());
    snp.set_clumped();
    // set clumped should set clump to true
    ASSERT_TRUE(snp.clumped());
}

TEST(SNP_BOUND, SET_LOW_BOUND)
{
    std::string rs = "Test";
    std::string ref = "G";
    std::string alt = "";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = 1, loc = -1;
    SNP snp(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
    ASSERT_EQ(snp.low_bound(), 0);
    ASSERT_EQ(snp.up_bound(), 0);
    snp.set_low_bound(10);
    ASSERT_EQ(snp.low_bound(), 10);
    ASSERT_EQ(snp.up_bound(), 0);
}
TEST(SNP_BOUND, SET_UP_BOUND)
{
    std::string rs = "Test";
    std::string ref = "G";
    std::string alt = "";
    double stat = 0.0;
    double p = 0.0;
    double p_threshold = 1;
    int category = 1;
    int chr = 1, loc = -1;
    SNP snp(rs, chr, loc, ref, alt, stat, p, category, p_threshold);
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
    //    b. Same p-value
    //       i. different location : SNP_B SNP_D
    //      ii. same location:  SNP_C SNP_E
    //          - Different name
    // 2. Different chromosome
    //    a. Same everything (except name) : SNP_A SNP_B

    // by definition, we don't allow multiple SNPs with the same name. Using SNP
    // name as the last comparison condition should allow us to avoid troubles
    snps.emplace_back(SNP("SNP_A", 1, 10, "A", "", 1, 0.05, 1, 0.05));
    snps.emplace_back(SNP("SNP_B", 2, 10, "A", "", 1, 0.05,  1, 0.05));
    snps.emplace_back(SNP("SNP_C", 1, 10, "A", "", 1, 0.01, 1, 0.05));
    snps.emplace_back(SNP("SNP_D", 2, 11, "A", "", 1, 0.05, 1, 0.05));
    snps.emplace_back(SNP("SNP_E", 1, 10, "A", "", 1, 0.01, 1, 0.05));

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
/*
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
*/
#endif // SNP_TEST_HPP
