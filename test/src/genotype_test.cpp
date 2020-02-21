#ifndef GENOTYPE_TEST_HPP
#define GENOTYPE_TEST_HPP
#include "genotype.hpp"
#include "gtest/gtest.h"
class GENOTYPE_BASIC : public Genotype, public ::testing::Test
{
};
// This should be the simplest of the three genotype class. Test anything that
// doesn't require binaryplink and binarygen
// set_genotype_files
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
