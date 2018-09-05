#ifndef PRSICE_TEST_HPP
#define PRSICE_TEST_HPP
#include "prsice.hpp"
#include "gtest/gtest.h"
#include <fstream>
class PRSICE_PHENO_PARSE : public PRSice, public ::testing::Test
{
protected:
    Reporter reporter;
    std::string pheno_name;
    void SetUp() override
    {
        pheno_name = "pheno";
        std::ofstream pheno_file;
        pheno_file.open(pheno_name.c_str());
        pheno_file << "FID IID Phenotype A B C D" << std::endl;
        pheno_file.close();
    }
};

TEST_F(PRSICE_PHENO_PARSE, ALL_PHENO_OCCUR)
{
    std::vector<std::string> col_names = {"Phenotype", "A", "B", "C", "D"};
    std::vector<bool> is_binary = {true, false, true, true, false};
    try
    {
        pheno_check(pheno_name, col_names, is_binary, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    ASSERT_EQ(num_phenotype(), col_names.size());
}

TEST_F(PRSICE_PHENO_PARSE, SOME_PHENO_MISS)
{
    std::vector<std::string> col_names = {"Phenotype", "E", "C"};
    std::vector<bool> is_binary = {true, false, true, true, false};
    try
    {
        pheno_check(pheno_name, col_names, is_binary, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    ASSERT_EQ(num_phenotype(), 2);
}

TEST_F(PRSICE_PHENO_PARSE, NO_PHENO_FILE)
{
    std::vector<std::string> col_names = {""};
    std::vector<bool> is_binary = {true};
    try
    {
        pheno_check("", col_names, is_binary, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    ASSERT_EQ(num_phenotype(), 1);
}

#endif // PRSICE_TEST_HPP
