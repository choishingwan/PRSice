#ifndef PRSICE_TEST_HPP
#define PRSICE_TEST_HPP
#include "global.hpp"
#include "gtest/gtest.h"
#include <prsice.hpp>
#include <reporter.hpp>
#include <storage.hpp>
#include <string>

TEST(PRSICE, CONSTRUCT)
{
    // should do nothing but assignment
    Reporter reporter(std::string(path + "LOG"));
    std::string output;
    CalculatePRS prs_info;
    PThresholding p_info;
    Phenotype pheno;
    Permutations perm;
    try
    {
        PRSice prsice(prs_info, p_info, pheno, perm, output, &reporter);
        SUCCEED();
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(PRSICE, PHENO_CHECK)
{
    // Test if the phenotype checking function works
    Reporter reporter(std::string(path + "LOG"));
    std::string output;
    CalculatePRS prs_info;
    // will still run if no_regress is use as we don't actually read the
    // phenotype file
    prs_info.no_regress = true;
    PThresholding p_info;
    Phenotype pheno;
    Permutations perm;
    PRSice no_regress(prs_info, p_info, pheno, perm, output, &reporter);
    // no phenotype file but no-regress requested, should set num pheno to 1 so
    // that we will run the loop in main
    no_regress.pheno_check();
    ASSERT_EQ(no_regress.num_phenotype(), 1);
    prs_info.no_regress = false;
    PRSice no_binary_info(prs_info, p_info, pheno, perm, output, &reporter);
    try
    {
        // we assume there's always at least one binary trait info
        no_binary_info.pheno_check();
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    pheno.binary.push_back(true);
    PRSice valid_no_pheno(prs_info, p_info, pheno, perm, output, &reporter);
    valid_no_pheno.pheno_check();
    ASSERT_EQ(valid_no_pheno.num_phenotype(), 1);

    pheno.pheno_file = std::string(path + "404");
    PRSice pheno_not_found(prs_info, p_info, pheno, perm, output, &reporter);
    try
    {
        pheno_not_found.pheno_check();
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
}

#endif
