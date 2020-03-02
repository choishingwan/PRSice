#ifndef PRSICE_TEST_HPP
#define PRSICE_TEST_HPP
#include "global.hpp"
#include "gtest/gtest.h"
#include <prsice.hpp>
#include <reporter.hpp>
#include <storage.hpp>
#include <string>
#include <vector>

TEST(PRSICE, CONSTRUCT)
{
    // should do nothing but assignment
    Reporter reporter(std::string(path + "LOG"), 60, true);
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
    Reporter reporter(std::string(path + "LOG"), 60, true);
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
    pheno.pheno_file = std::string(path + "Phenotype");
    // Phenotype file, but without phenotype column and not ignoring FID
    PRSice pheno_single_no_col(prs_info, p_info, pheno, perm, output,
                               &reporter);
    try
    {
        pheno_single_no_col.pheno_check();
        ASSERT_EQ(pheno_single_no_col.num_phenotype(), 1);
        // will use header as the phenotype name
        ASSERT_STREQ(pheno_single_no_col.pheno_name(0).c_str(), "CC_Valid");
    }
    catch (const std::runtime_error&)
    {
        FAIL();
    }
    pheno.pheno_col.clear();
    pheno.pheno_col.push_back("404");
    PRSice pheno_single_not_found(prs_info, p_info, pheno, perm, output,
                                  &reporter);
    // if no header, and pheno_col is provided, then it will be similar to
    // situation where none of the items were found
    try
    {
        pheno_single_not_found.pheno_check();
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    pheno.pheno_col.clear();
    pheno.binary.clear();
    // multiple check
    pheno.pheno_col = {"CC_Valid", "QT", "QT", "404", "QT_miss", "CC_miss"};
    pheno.binary = {true, false, false, true, false, true};
    pheno.prevalence = {0.1, 0.4, 0.3};
    PRSice pheno_multi_check(prs_info, p_info, pheno, perm, output, &reporter);
    pheno_multi_check.pheno_check();
    // 4 because we should remove 1 duplicate and 1 not found
    ASSERT_EQ(pheno_multi_check.num_phenotype(), 4);
    ASSERT_STREQ(pheno_multi_check.pheno_name(0).c_str(), "CC_Valid");
    ASSERT_STREQ(pheno_multi_check.pheno_name(1).c_str(), "QT");
    ASSERT_STREQ(pheno_multi_check.pheno_name(2).c_str(), "QT_miss");
    ASSERT_STREQ(pheno_multi_check.pheno_name(3).c_str(), "CC_miss");
    std::vector<double> prevalence = pheno_multi_check.get_prevalence();
    ASSERT_EQ(prevalence.size(), 2);
    ASSERT_DOUBLE_EQ(prevalence[0], 0.1);
    ASSERT_DOUBLE_EQ(prevalence[1], 0.3);

    pheno.pheno_col.clear();
    pheno.binary.clear();
    pheno.binary = {true};
    pheno.pheno_file = std::string(path + "Pheno.misshead");
    PRSice pheno_miss_header(prs_info, p_info, pheno, perm, output, &reporter);
    try
    {
        pheno_miss_header.pheno_check();
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    pheno.pheno_file = std::string(path + "Pheno.noheader");
    // Phenotype file doesn't have a header, should use Phenotype as default
    // name
    PRSice pheno_no_header(prs_info, p_info, pheno, perm, output, &reporter);
    pheno_no_header.pheno_check();
    ASSERT_STREQ(pheno_no_header.pheno_name(0).c_str(), "Phenotype");

    pheno.pheno_file = std::string(path + "invalid.pheno");
    // this file only contain 1 column, which can't be a phenotype file
    PRSice invalid_pheno(prs_info, p_info, pheno, perm, output, &reporter);
    try
    {
        invalid_pheno.pheno_check();
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    pheno.pheno_file = std::string(path + "Pheno.dup");
    pheno.pheno_col.clear();
    pheno.pheno_col.push_back("CC_Valid");
    PRSice dup_pheno_invalid_col(prs_info, p_info, pheno, perm, output,
                                 &reporter);
    try
    {
        // don't allow duplicated phenotype in file that we need
        dup_pheno_invalid_col.pheno_check();
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    pheno.pheno_col.clear();
    pheno.pheno_col.push_back("CC_miss");
    PRSice dup_pheno_valid_col(prs_info, p_info, pheno, perm, output,
                               &reporter);
    try
    {
        // allow duplicated column as long as we don't need it
        dup_pheno_valid_col.pheno_check();
        ASSERT_EQ(dup_pheno_valid_col.num_phenotype(), 1);
    }
    catch (const std::runtime_error&)
    {
        FAIL();
    }
}


#endif
