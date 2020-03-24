#ifndef GENOTYPE_TEST_HPP
#define GENOTYPE_TEST_HPP
#include "genotype.hpp"
#include "gtest/gtest.h"
class GENOTYPE_BASIC : public Genotype, public ::testing::Test
{
public:
    void cleanup()
    {
        m_keep_file = "";
        m_sample_file = "";
        m_remove_file = "";
        m_genotype_file_names.clear();
        m_delim = "";
        m_ignore_fid = false;
    }
};
// This should be the simplest of the three genotype class. Test anything that
// doesn't require binaryplink and binarygen
// set_genotype_files

TEST_F(GENOTYPE_BASIC, CHR_CONVERT)
{
    ASSERT_FALSE(chr_prefix("1"));
    ASSERT_TRUE(chr_prefix("chr1"));
    auto chr = get_chrom_code("1");
    ASSERT_EQ(chr, 1);
    chr = get_chrom_code("chr10");
    ASSERT_EQ(chr, 10);
    chr = get_chrom_code("CHRX");
    ASSERT_EQ(chr, CHROM_X);
    chr = get_chrom_code("M");
    ASSERT_EQ(chr, CHROM_MT);
    chr = get_chrom_code("Y");
    ASSERT_EQ(chr, CHROM_Y);
    chr = get_chrom_code("XY");
    ASSERT_EQ(chr, CHROM_XY);
    chr = get_chrom_code("mt");
    ASSERT_EQ(chr, CHROM_MT);
    chr = get_chrom_code("CHR Happy");
    ASSERT_EQ(chr, -1);
    chr = get_chrom_code("0X");
    ASSERT_EQ(chr, CHROM_X);
    chr = get_chrom_code("0Y");
    ASSERT_EQ(chr, CHROM_Y);
    chr = get_chrom_code("0M");
    ASSERT_EQ(chr, CHROM_MT);
    // TODO: Add more crazy chromosome use cases
}

TEST_F(GENOTYPE_BASIC, INITIALIZE)
{
    GenoFile geno;
    geno.num_autosome = 22;
    Phenotype pheno;
    pheno.ignore_fid = false;
    std::string delim = " ";
    std::string type = "bed";
    Reporter reporter(std::string("LOG"), 60, true);
    // need to check the following
    // m_genotype_file_name is set
    // m_sample_file is set when external file is provided and empty if not
    // m_keep_file, m_remove_file and m_delim should be set
    // there should be no error check
    // Normal input
    geno.file_name = "Genotype";
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_STREQ(delim.c_str(), m_delim.c_str());
    ASSERT_EQ(m_genotype_file_names.size(), 1);
    ASSERT_STREQ(geno.file_name.c_str(), m_genotype_file_names.front().c_str());
    ASSERT_TRUE(m_keep_file.empty());
    ASSERT_TRUE(m_remove_file.empty());
    ASSERT_TRUE(m_sample_file.empty());
    ASSERT_FALSE(m_ignore_fid);
    cleanup();
    // now with external sample_file
    geno.file_name = "Genotype,External";
    delim = "\t";
    pheno.ignore_fid = true;
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_STREQ(delim.c_str(), m_delim.c_str());
    ASSERT_EQ(m_genotype_file_names.size(), 1);
    ASSERT_STREQ("Genotype", m_genotype_file_names.front().c_str());
    ASSERT_TRUE(m_keep_file.empty());
    ASSERT_TRUE(m_remove_file.empty());
    ASSERT_TRUE(m_ignore_fid);
    ASSERT_STREQ(m_sample_file.c_str(), "External");
    // now check replacement
    cleanup();
    geno.file_name = "chr#_geno";
    geno.keep = "Hi";
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_EQ(m_genotype_file_names.size(), geno.num_autosome);
    for (size_t i = 0; i < geno.num_autosome; ++i)
    {
        std::string cur_file = "chr" + std::to_string(i + 1) + "_geno";
        ASSERT_STREQ(m_genotype_file_names[i].c_str(), cur_file.c_str());
    }
    ASSERT_STREQ(m_keep_file.c_str(), "Hi");
    ASSERT_TRUE(m_remove_file.empty());
    ASSERT_TRUE(m_sample_file.empty());
    cleanup();
    // invalid format
    geno.file_name = "chr,x,y";
    try
    {
        initialize(geno, pheno, delim, type, &reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    cleanup();
    geno.file_name = "chr#_geno,Sample";
    geno.keep = "Hi";
    geno.num_autosome = 30;
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_EQ(m_genotype_file_names.size(), geno.num_autosome);
    for (size_t i = 0; i < geno.num_autosome; ++i)
    {
        std::string cur_file = "chr" + std::to_string(i + 1) + "_geno";
        ASSERT_STREQ(m_genotype_file_names[i].c_str(), cur_file.c_str());
    }
    ASSERT_STREQ(m_keep_file.c_str(), "Hi");
    ASSERT_TRUE(m_remove_file.empty());
    ASSERT_STREQ(m_sample_file.c_str(), "Sample");
    // finally, with list
    cleanup();
    // we need a dummy file
    std::ofstream dummy("DUMMY");
    std::vector<std::string> expected = {"A", "B", "C D", "E F G", "H#I"};
    for (auto& d : expected) { dummy << d << std::endl; }
    dummy.close();
    geno.file_list = "DUMMY";
    geno.file_name = "";
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_EQ(m_genotype_file_names.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i)
    { ASSERT_STREQ(m_genotype_file_names[i].c_str(), expected[i].c_str()); }
    ASSERT_TRUE(m_sample_file.empty());
    // now try with external sample
    cleanup();
    geno.file_list = "DUMMY,sample";
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_EQ(m_genotype_file_names.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i)
    { ASSERT_STREQ(m_genotype_file_names[i].c_str(), expected[i].c_str()); }
    ASSERT_FALSE(m_sample_file.empty());
    ASSERT_STREQ(m_sample_file.c_str(), "sample");

    // invalid file input
    cleanup();
    geno.file_list = "DUMMY,B,C";
    try
    {
        initialize(geno, pheno, delim, type, &reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    cleanup();
    std::remove("DUMMY");
    cleanup();
    // should fail here when the file is not found
    try
    {
        initialize(geno, pheno, delim, type, &reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

TEST_F(GENOTYPE_BASIC, SET_FILE_NAME_WITHOUT_HASH)
{
    std::string name = "Test";
    m_genotype_file_names = set_genotype_files(name);
    ASSERT_STREQ(m_genotype_file_names.front().c_str(), name.c_str());
    ASSERT_EQ(m_genotype_file_names.size(), 1);
}

TEST_F(GENOTYPE_BASIC, SET_FILE_NAME_WITH_HASH)
{
    std::string name = "chr#test";
    m_autosome_ct = 22;
    m_genotype_file_names = set_genotype_files(name);
    ASSERT_EQ(m_genotype_file_names.size(), m_autosome_ct);
    for (size_t i = 1; i <= m_autosome_ct; ++i)
    {
        std::string name = "chr" + std::to_string(i) + "test";
        ASSERT_STREQ(m_genotype_file_names[i - 1].c_str(), name.c_str());
    }
}

TEST_F(GENOTYPE_BASIC, SET_FILE_NAME_MULTI_HASH)
{
    // don't think this is a useful case but in theory all # should be
    // substituted
    std::string name = "chr#test#";
    m_autosome_ct = 22;
    m_genotype_file_names = set_genotype_files(name);
    ASSERT_EQ(m_genotype_file_names.size(), m_autosome_ct);
    for (size_t i = 1; i <= m_autosome_ct; ++i)
    {
        std::string name =
            "chr" + std::to_string(i) + "test" + std::to_string(i);
        ASSERT_STREQ(m_genotype_file_names[i - 1].c_str(), name.c_str());
    }
}


TEST_F(GENOTYPE_BASIC, GET_RS_COLUMN)
{
    Reporter reporter(std::string("LOG"), 60, true);
    m_reporter = &reporter;
    // simplest format
    std::string input = "rs1234";
    auto result = get_rs_column(input);
    ASSERT_EQ(result, 0);
    // might want to test the rest, RS.ID RS_ID, RSID, SNP.ID, SNP_ID,
    // Variant_ID, Variant.ID
    input = "A B CHR snp Weight P A1";
    // read from header
    result = get_rs_column(input);
    ASSERT_EQ(result, 3);
    input = "1 rs1234 123 0 A T";
    // this should be bim file
    result = get_rs_column(input);
    ASSERT_EQ(result, 1);
    input = "rs1234 A BC T A D E T";
    // no idea what this is, will take first column
    result = get_rs_column(input);
    ASSERT_EQ(result, 0);
}
TEST_F(GENOTYPE_BASIC, LOAD_SNP_LIST)
{
    Reporter reporter("LOG", 60, true);
    m_reporter = &reporter;
    std::vector<std::string> expected = {"rs1234", "rs76591", "rs139486"};
    std::ofstream dummy("DUMMY");
    for (auto& e : expected) { dummy << e << std::endl; }
    auto res = load_snp_list("DUMMY");
    ASSERT_EQ(res.size(), expected.size());
    for (auto& e : expected) { ASSERT_FALSE(res.find(e) == res.end()); }
    dummy.close();
    std::remove("DUMMY");
    // now different format, just in case the idx is correct (use bim)
    dummy.open("DUMMY");
    for (auto& e : expected) { dummy << "1 " << e << " 3 4 5 6" << std::endl; }
    dummy.close();
    res = load_snp_list("DUMMY");

    ASSERT_EQ(res.size(), expected.size());
    for (auto& e : expected) { ASSERT_FALSE(res.find(e) == res.end()); }
    std::remove("DUMMY");
    // should fail if file not found
    try
    {
        res = load_snp_list("DUMMY");
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

// init_chr
// chr_code_check
// load_snp_list
// load_ref

TEST_F(GENOTYPE_BASIC, AMBIGUOUS)
{
    ASSERT_TRUE(ambiguous("A", "A"));
    ASSERT_TRUE(ambiguous("G", "G"));
    ASSERT_TRUE(ambiguous("C", "C"));
    ASSERT_TRUE(ambiguous("T", "T"));
    ASSERT_TRUE(ambiguous("A", "T"));
    ASSERT_TRUE(ambiguous("G", "C"));
    ASSERT_TRUE(ambiguous("C", "G"));
    ASSERT_TRUE(ambiguous("T", "A"));
    ASSERT_FALSE(ambiguous("A", "G"));
    ASSERT_FALSE(ambiguous("G", "A"));
    ASSERT_FALSE(ambiguous("C", "T"));
    ASSERT_FALSE(ambiguous("T", "G"));
    ASSERT_FALSE(ambiguous("A", "C"));
    ASSERT_FALSE(ambiguous("G", "T"));
    ASSERT_FALSE(ambiguous("C", "A"));
    ASSERT_FALSE(ambiguous("T", "C"));
    ASSERT_FALSE(ambiguous("A", ""));
    ASSERT_FALSE(ambiguous("G", ""));
    ASSERT_FALSE(ambiguous("C", ""));
    ASSERT_FALSE(ambiguous("T", ""));
    ASSERT_FALSE(ambiguous("A", "TA"));
    ASSERT_FALSE(ambiguous("G", "CG"));
    ASSERT_FALSE(ambiguous("C", "GC"));
    ASSERT_FALSE(ambiguous("T", "AT"));
}

TEST_F(GENOTYPE_BASIC, CATEGORY)
{
    unsigned long long category;
    double pthres;
    // anything less than lowest threshold is consider 0
    PThresholding mock;
    mock.lower = 0.05;
    mock.inter = 0.01;
    mock.upper = 0.5;
    mock.no_full = true;
    double p_value = 0;
    category = calculate_category(p_value, pthres, mock);
    ASSERT_EQ(category, 0);
    ASSERT_DOUBLE_EQ(pthres, 0.05);
    // anything above the upper threshold is considered as 1.0 and with category
    // higher than the biggest category
    // though in theory, this should never happen as we will always filter SNPs
    // that are bigger than the last required threshold
    // the pthres is meaningless in this situation
    p_value = 0.6;
    category = calculate_category(p_value, pthres, mock);
    ASSERT_GT(category, mock.upper / mock.inter - mock.lower / mock.inter);
    // if we want full model, we will still do the same
    mock.no_full = false;
    category = calculate_category(p_value, pthres, mock);
    ASSERT_GT(category, mock.upper / mock.inter - mock.lower / mock.inter);
    ASSERT_DOUBLE_EQ(pthres, 1);
    mock.no_full = true;
    p_value = 0.05;
    category = calculate_category(p_value, pthres, mock);
    // this should be the first threshold
    ASSERT_EQ(category, 0);
    ASSERT_DOUBLE_EQ(pthres, 0.05);
    p_value = 0.055;
    category = calculate_category(p_value, pthres, mock);
    // this should be the first threshold
    ASSERT_EQ(category, 1);
    ASSERT_DOUBLE_EQ(pthres, 0.06);
}

void category_threshold(unsigned long long category,
                        unsigned long long expected_category, double pthres,
                        double expected_pthres)
{
    ASSERT_EQ(category, expected_category);
    ASSERT_DOUBLE_EQ(pthres, expected_pthres);
}
TEST_F(GENOTYPE_BASIC, BAR_LEVELS)
{
    unsigned long long category;
    double pthres;
    std::vector<double> barlevels = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    // anything less than lowest threshold is consider 0
    category = calculate_category(0.0, barlevels, pthres);
    category_threshold(category, 0, pthres, 0.001);
    category = calculate_category(0.001, barlevels, pthres);
    category_threshold(category, 0, pthres, 0.001);
    category = calculate_category(0.5, barlevels, pthres);
    category_threshold(category, 6, pthres, 0.5);
    category = calculate_category(0.7, barlevels, pthres);
    category_threshold(category, 7, pthres, 1);
}


#endif // GENOTYPE_TEST_HPP
