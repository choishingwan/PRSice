#ifndef BIN_PLINK_TEST_H
#define BIN_PLINK_TEST_H
// this'd be the second most difficult but also one of the most important test
// unit we need to do, to test the behaviour of the binaryplink class. But if I
// am able to manage it, it will truely mark the time where we can finally
// officially release PRSice
#include "binaryplink.hpp"
#include "genotype.hpp"
#include "global.hpp"
#include "reporter.hpp"
#include "gtest/gtest.h"
#include <fstream>
#include <math.h>
#include <storage.hpp>

class BPLINK_GEN_SAMPLE_TARGET : public ::testing::Test
{
protected:
    BinaryPlink* plink = nullptr;
    Reporter* reporter = nullptr;
    void SetUp() override
    {
        std::string file_list;
        std::string file = std::string(path + "TOY_TARGET_DATA");
        GenoFile geno;
        geno.file_list = file_list;
        geno.file_name = file;
        geno.is_ref = false;
        Phenotype pheno;
        pheno.ignore_fid = false;
        reporter = new Reporter(std::string(path + "LOG"));
        plink = new BinaryPlink(geno, pheno, " ", reporter);
    }
    void TearDown() override
    {
        delete plink;
        delete reporter;
    }
};

TEST_F(BPLINK_GEN_SAMPLE_TARGET, NO_SELECTION)
{
    std::string keep_file = "";
    std::string remove_file = "";
    bool verbose = true;
    std::string delim = " ";
    plink->load_samples(verbose);
    ASSERT_EQ(plink->num_sample(), 2000);
}
/*
TEST_F(BPLINK_GEN_SAMPLE_TARGET, KEEP_SAMPLE)
{
    std::string keep_file = path + "keep";
    std::string delim = " ";
    std::ofstream keep;
    keep.open(keep_file.c_str());
    // We are using toy sample, so the sample naming convention should be
    // CAS_XXX CONT_XXX where XXX range from 1 to 1000
    for (size_t i = 123; i < 234; ++i)
    {
        keep << "CAS_" << i << "\t"
             << "CAS_" << i << std::endl;
        keep << "CONT_" << i << "\t"
             << "CONT_" << i << std::endl;
    }

    keep.close();
    std::string remove_file = "";
    bool verbose = true;
    plink->load_samples(verbose);
    ASSERT_EQ(plink->num_sample(), 222);
}
TEST_F(BPLINK_GEN_SAMPLE_TARGET, REMOVE_SAMPLE)
{
    std::string keep_file = "";
    std::string remove_file = path + "remove";
    std::string delim = " ";
    std::ofstream remove;
    remove.open(remove_file.c_str());
    // We are using toy sample, so the sample naming convention should be
    // CAS_XXX CONT_XXX where XXX range from 1 to 1000
    for (size_t i = 321; i < 432; ++i)
    {
        remove << "CAS_" << i << "\t"
               << "CAS_" << i << std::endl;
        remove << "CONT_" << i << "\t"
               << "CONT_" << i << std::endl;
    }

    remove.close();
    bool verbose = true;
    plink->load_samples(keep_file, remove_file, delim, verbose);
    ASSERT_EQ(plink->num_sample(), 1778);
}
TEST_F(BPLINK_GEN_SAMPLE_TARGET, KEEP_SAMPLE_IID)
{
    std::string keep_file = path + "keep";
    std::string delim = " ";
    std::ofstream keep;
    keep.open(keep_file.c_str());
    // We are using toy sample, so the sample naming convention should be
    // CAS_XXX CONT_XXX where XXX range from 1 to 1000
    for (size_t i = 123; i < 234; ++i)
    {
        keep << "CAS_" << i << std::endl;
        keep << "CONT_" << i << std::endl;
    }

    keep.close();
    std::string remove_file = "";
    bool verbose = true;
    // we have set ignore_fid to false, therefore when only IID is provided,
    // PRSice should error out
    try
    {
        plink->load_samples(keep_file, remove_file, delim, verbose);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

TEST(BPLINK_EXTERNAL, EXTERNAL_SAMPLE)
{
    std::string file_list;
    std::string delim = " ";
    std::string file =
        std::string(path + "TOY_TARGET_DATA," + path + "TOY_TARGET_DATA.fam");
    uint32_t thread = 1;
    bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
         is_ref = false;
    Reporter reporter(std::string(path + "LOG"));
    BinaryPlink plinkBinary(file_list, file, thread, ignore_fid,
                            keep_nonfounder, keep_ambig, is_ref, &reporter);
    plinkBinary.load_samples("", "", delim, true);
    ASSERT_EQ(plinkBinary.num_sample(), 2000);
}
TEST(BPLINK_SAMPLE_CHECK, DUPLICATE_SAMPLE)
{
    std::string file_list;
    std::string file = std::string(path + "TOY_TARGET_DATA," + path
                                   + "TOY_TARGET_DATA.dup.fam");
    std::string delim = " ";
    size_t thread = 1;
    bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
         is_ref = false;
    Reporter reporter(std::string(path + "LOG"));
    BinaryPlink plinkBinary(file_list, file, thread, ignore_fid,
                            keep_nonfounder, keep_ambig, is_ref, &reporter);
    try
    {
        plinkBinary.load_samples("", "", delim, true);
        FAIL();
    }
    catch (...)
    {
        // with duplicate samples, we should error out, though we won't generate
        // a valid sample file
        SUCCEED();
    }
}
// check founder samples
TEST(BPLINK_FOUNDER, FOUNDER_REMOVE)
{
    std::string file_list;
    std::string delim = " ";
    std::string file = std::string(path + "TOY_TARGET_DATA," + path
                                   + "TOY_TARGET_DATA.founder.fam");
    uint32_t thread = 1;
    bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
         is_ref = false;
    Reporter reporter(std::string(path + "LOG"));
    BinaryPlink plinkBinary(file_list, file, thread, ignore_fid,
                            keep_nonfounder, keep_ambig, is_ref, &reporter);
    plinkBinary.load_samples("", "", delim, true);
    // we still keep the non-founders, just not using them for regression
    ASSERT_EQ(plinkBinary.num_sample(), 2000);
    int sum_founder = 0;
    for (size_t i = 0; i < 2000; ++i)
    { sum_founder += plinkBinary.is_founder(i); }
    // NOTE: This was intended to be 859. But our dummy has two of the
    // founders set to have a different FID, the script correctly account for
    // samples in these family to be a fonuder, because while the IID of
    // paternal and maternal is here, that sample is from a different family
    ASSERT_EQ(sum_founder, 862);
}
*/
/*
class BPLINK_GEN_SNP_TARGET : public ::testing::Test
{
protected:
    BinaryPlink* plink = nullptr;
    void SetUp() override
    {
        std::string file_list;
        std::string file = std::string(path + "TOY_TARGET_MISS");
        uint32_t thread = 1;
        bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
             is_ref = false;
        Reporter reporter(std::string(path + "LOG"));
        plink = new BinaryPlink(file_list, file, thread, ignore_fid,
                                keep_nonfounder, keep_ambig, is_ref, reporter);
        plink->load_samples("", "", true, reporter);
    }
    void TearDown() override { delete plink; }
};
// Test maf filtering
TEST_F(BPLINK_GEN_SNP_TARGET, SIMPLE_READ)
{
    std::string out = path + "test";
    Reporter reporter(std::string(path + "LOG"));
    Region exclusion_region("", reporter);

    double maf = 0.0;
    double geno = 0.0;
    double info = 0.0;
    double hard_threshold = 0.0;
    bool maf_filter = false;
    bool geno_filter = false;
    bool hard_coded = false;
    bool info_filter = false;
    bool verbose = true;
    try
    {
        plink->load_snps(out, "", "", maf, geno, info, hard_threshold,
                         maf_filter, geno_filter, hard_coded, info_filter,
                         exclusion_region, verbose, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    // remember, we need to remove chr 22+
    ASSERT_EQ(plink->num_snps(), 88836);
}
TEST_F(BPLINK_GEN_SNP_TARGET, MAF_FILTERING_1)
{
    std::string out = path + "test";
    Reporter reporter(std::string(path + "LOG"));
    Region exclusion_region("", reporter);

    double maf = 0.2;
    double geno = 0.0;
    double info = 0.0;
    double hard_threshold = 0.0;
    bool maf_filter = true;
    bool geno_filter = false;
    bool hard_coded = false;
    bool info_filter = false;
    bool verbose = true;
    try
    {
        plink->load_snps(out, "", "", maf, geno, info, hard_threshold,
                         maf_filter, geno_filter, hard_coded, info_filter,
                         exclusion_region, verbose, reporter);
    }
    catch (...)
    {
        FAIL();
    }

    ASSERT_EQ(plink->num_snps(), 29465);
}
TEST_F(BPLINK_GEN_SNP_TARGET, MAF_FILTERING_2)
{
    std::string out = path + "test";
    Reporter reporter(std::string(path + "LOG"));
    Region exclusion_region("", reporter);

    double maf = 0.2068;
    double geno = 0.0;
    double info = 0.0;
    double hard_threshold = 0.0;
    bool maf_filter = true;
    bool geno_filter = false;
    bool hard_coded = false;
    bool info_filter = false;
    bool verbose = true;
    try
    {
        plink->load_snps(out, "", "", maf, geno, info, hard_threshold,
                         maf_filter, geno_filter, hard_coded, info_filter,
                         exclusion_region, verbose, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    ASSERT_EQ(plink->num_snps(), 7414);
}
TEST_F(BPLINK_GEN_SNP_TARGET, GENO_FILTERING)
{
    std::string out = path + "test";
    Reporter reporter(std::string(path + "LOG"));
    Region exclusion_region("", reporter);

    double maf = 0.2068;
    double geno = 0.05;
    double info = 0.0;
    double hard_threshold = 0.0;
    bool maf_filter = false;
    bool geno_filter = true;
    bool hard_coded = false;
    bool info_filter = false;
    bool verbose = true;
    try
    {
        plink->load_snps(out, "", "", maf, geno, info, hard_threshold,
                         maf_filter, geno_filter, hard_coded, info_filter,
                         exclusion_region, verbose, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    // We follow convension of PLINK
    ASSERT_EQ(plink->num_snps(), 88772);
}
TEST_F(BPLINK_GEN_SNP_TARGET, MAF_GENO_FILTERING)
{
    std::string out = path + "test";
    Reporter reporter(std::string(path + "LOG"));
    Region exclusion_region("", reporter);

    double maf = 0.2068;
    double geno = 0.05;
    double info = 0.0;
    double hard_threshold = 0.0;
    bool maf_filter = true;
    bool geno_filter = true;
    bool hard_coded = false;
    bool info_filter = false;
    bool verbose = true;
    try
    {
        plink->load_snps(out, "", "", maf, geno, info, hard_threshold,
                         maf_filter, geno_filter, hard_coded, info_filter,
                         exclusion_region, verbose, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    // We follow convension of PLINK
    ASSERT_EQ(plink->num_snps(), 7407);
}
TEST_F(BPLINK_GEN_SNP_TARGET, TEST_EXCLUSION_FUNCTION)
{
    // we will only test it once to ensure BinaryPlink does call the exclusion
    // function. We don't need to test the exclusion performance as that should
    // be tested in region's unit testing
    std::string out = path + "test";
    Reporter reporter(std::string(path + "LOG"));
    // NOTE: The end range boundary is exclusion. Any SNP with exact match on
    // that number will not be excluded. E.g. SNP on chr1 2842568
    Region exclusion_region("1:2832179-2842568", reporter);

    double maf = 0.0;
    double geno = 0.0;
    double info = 0.0;
    double hard_threshold = 0.0;
    bool maf_filter = false;
    bool geno_filter = false;
    bool hard_coded = false;
    bool info_filter = false;
    bool verbose = true;
    try
    {
        plink->load_snps(out, "", "", maf, geno, info, hard_threshold,
                         maf_filter, geno_filter, hard_coded, info_filter,
                         exclusion_region, verbose, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    // remember, we need to remove chr 22+
    ASSERT_EQ(plink->num_snps(), 88834);
}
TEST(BPLINK_GEN_SNP, DUP_SNP)
{
    std::string file_list;
    std::string file = std::string(path + "TOY_TARGET_DUP");
    uint32_t thread = 1;
    bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
         is_ref = false;
    Reporter reporter(std::string(path + "LOG"));
    BinaryPlink plink(file_list, file, thread, ignore_fid, keep_nonfounder,
                      keep_ambig, is_ref, reporter);
    plink.load_samples("", "", true, reporter);
    std::string out = path + "test";
    Region exclusion_region("", reporter);

    double maf = 0.2068;
    double geno = 0.05;
    double info = 0.0;
    double hard_threshold = 0.0;
    bool maf_filter = false;
    bool geno_filter = false;
    bool hard_coded = false;
    bool info_filter = false;
    bool verbose = true;
    try
    {
        plink.load_snps(out, "", "", maf, geno, info, hard_threshold,
                        maf_filter, geno_filter, hard_coded, info_filter,
                        exclusion_region, verbose, reporter);
        // there are duplicated SNPs, therefore it should fail immediately
        FAIL();
    }
    catch (std::runtime_error& e)
    {
        reporter.report(e.what());
        std::ifstream valid(std::string(out + ".valid").c_str());
        if (!valid.is_open()) {
            // we should have the valid file
            FAIL();
        }
        int num_retain = 0;
        std::string line;
        while (std::getline(valid, line)) {
            num_retain++;
        }
        ASSERT_EQ(num_retain, 88640);
    }
    catch (...)
    {
        // should only throw runtime_error
        FAIL();
    }
}
TEST(BPLINK_GEN_SNP, SEQ_INPUT)
{
    std::string file_list;
    std::string file = std::string(path + "chr#");
    uint32_t thread = 1;
    bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
         is_ref = false;
    Reporter reporter(std::string(path + "LOG"));
    BinaryPlink plink(file_list, file, thread, ignore_fid, keep_nonfounder,
                      keep_ambig, is_ref, reporter);
    plink.load_samples("", "", true, reporter);
    std::string out = path + "test";
    Region exclusion_region("", reporter);

    double maf = 0.2068;
    double geno = 0.05;
    double info = 0.0;
    double hard_threshold = 0.0;
    bool maf_filter = false;
    bool geno_filter = false;
    bool hard_coded = false;
    bool info_filter = false;
    bool verbose = true;
    try
    {
        plink.load_snps(out, "", "", maf, geno, info, hard_threshold,
                        maf_filter, geno_filter, hard_coded, info_filter,
                        exclusion_region, verbose, reporter);
    }
    catch (...)
    {
        FAIL();
    }
    ASSERT_EQ(plink.num_snps(), 88836);
}
// Now test performance of reading base
class BPLINK_BASE_READ : public ::testing::Test
{
protected:
    BinaryPlink* plink = nullptr;
    void SetUp() override
    {
        std::string file_list;
        std::string file = std::string(path + "TOY_TARGET_MISS");
        uint32_t thread = 1;
        bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
             is_ref = false;
        Reporter reporter(std::string(path + "LOG"));
        plink = new BinaryPlink(file_list, file, thread, ignore_fid,
                                keep_nonfounder, keep_ambig, is_ref, reporter);
        plink->load_samples("", "", true, reporter);
        std::string out = path + "test";
        Region exclusion_region("", reporter);
        double maf = 0.0;
        double geno = 0.0;
        double info = 0.0;
        double hard_threshold = 0.0;
        bool maf_filter = false;
        bool geno_filter = false;
        bool hard_coded = false;
        bool info_filter = false;
        bool verbose = true;
        plink->load_snps(out, "", "", maf, geno, info, hard_threshold,
                         maf_filter, geno_filter, hard_coded, info_filter,
                         exclusion_region, verbose, reporter);
    }
    void TearDown() override { delete plink; }
};
// we know the mismatch checking and ambiguous algorithm works (from snp_test).
// Therefore we only need to check if read_base involve those command.
TEST_F(BPLINK_BASE_READ, SIMPLE_BETA)
{
    // TOY file should still be OR, but should be ok
    std::string base = path + "TOY_BASE_GWAS.short.assoc";
    std::string out = path + "test";
    std::vector<int> index = {1, 3, 4, 6, 0, 2, -1, 5, -1, -1, -1, 6};
    std::vector<double> barlevels = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    double lower = 5e-5;
    double inter = 0.0001;
    double upper = 0.5;
    double maf_control = 0;
    double maf_case = 0;
    double info_score = 0;
    bool maf_control_filter = false;
    bool maf_case_filter = false;
    bool info_filter = false;
    bool fastscore = false;
    bool no_full = false;
    bool beta = true;
    bool is_index = false;
    bool perform_shrinkage = false;
    std::vector<std::string> feature;
    Region region(feature, 0, 0, false, false);
    Reporter reporter(std::string(path + "LOG"));
    try
    {
        plink->read_base(base, out, index, barlevels, lower, inter, upper,
                         maf_control, maf_case, info_score, maf_control_filter,
                         maf_case_filter, info_filter, fastscore, no_full, beta,
                         is_index, perform_shrinkage, region, reporter);
    }
    catch (const std::runtime_error& e)
    {
        reporter.report(e.what());
        FAIL();
    }
    ASSERT_EQ(plink->num_snps(), 1947);
    // NOTE: the order of SNPs is based on the TOY_TARGET_DATA, not based on the
    // base. Therefore the expected statistic should be 0.08
    auto&& snp = plink->get_snp(0);
    ASSERT_DOUBLE_EQ(snp.stat(), 9.323);
    ASSERT_DOUBLE_EQ(snp.p_value(), 0.2348);
    SUCCEED();
}
TEST_F(BPLINK_BASE_READ, SIMPLE_OR)
{
    std::string base = path + "TOY_BASE_GWAS.short.assoc";
    std::string out = path + "test";
    std::vector<int> index = {1, 3, 4, 6, 0, 2, -1, 5, -1, -1, -1, 6};
    std::vector<double> barlevels = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    double lower = 5e-5;
    double inter = 0.0001;
    double upper = 0.5;
    double maf_control = 0;
    double maf_case = 0;
    double info_score = 0;
    bool maf_control_filter = false;
    bool maf_case_filter = false;
    bool info_filter = false;
    bool fastscore = false;
    bool no_full = false;
    bool beta = false;
    bool is_index = false;
    bool perform_shrinkage = false;
    std::vector<std::string> feature;
    Region region(feature, 0, 0, false, false);
    Reporter reporter(std::string(path + "LOG"));
    try
    {
        plink->read_base(base, out, index, barlevels, lower, inter, upper,
                         maf_control, maf_case, info_score, maf_control_filter,
                         maf_case_filter, info_filter, fastscore, no_full, beta,
                         is_index, perform_shrinkage, region, reporter);
    }
    catch (const std::runtime_error& e)
    {
        reporter.report(e.what());
        FAIL();
    }
    ASSERT_EQ(plink->num_snps(), 1947);
    auto&& snp = plink->get_snp(0);
    ASSERT_DOUBLE_EQ(snp.stat(), std::log(9.323));
    SUCCEED();
}
TEST_F(BPLINK_BASE_READ, PROBLEM_OR)
{
    std::string base = path + "Wrong.OR.assoc";
    std::string out = path + "test";
    std::vector<int> index = {1, 3, 4, 6, 0, 2, -1, 5, -1, -1, -1, 6};
    std::vector<double> barlevels = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    double lower = 5e-5;
    double inter = 0.0001;
    double upper = 0.5;
    double maf_control = 0;
    double maf_case = 0;
    double info_score = 0;
    bool maf_control_filter = false;
    bool maf_case_filter = false;
    bool info_filter = false;
    bool fastscore = false;
    bool no_full = false;
    bool beta = false;
    bool is_index = false;
    bool perform_shrinkage = false;
    std::vector<std::string> feature;
    Region region(feature, 0, 0, false, false);
    Reporter reporter(std::string(path + "LOG"));
    try
    {
        plink->read_base(base, out, index, barlevels, lower, inter, upper,
                         maf_control, maf_case, info_score, maf_control_filter,
                         maf_case_filter, info_filter, fastscore, no_full, beta,
                         is_index, perform_shrinkage, region, reporter);
    }
    catch (const std::runtime_error& e)
    {
        // this should still run, but all SNP with negative OR or an OR of 0
        // will be removed, but should give a warning
        // Or do we want to make it more stringent?
        // TODO: Ask users (though it seems like people ain't too active)
        reporter.report(e.what());
        FAIL();
    }
    ASSERT_EQ(plink->num_snps(), 977);
    auto&& snp = plink->get_snp(0);
    ASSERT_DOUBLE_EQ(snp.stat(), std::log(0.95912597246771));
    SUCCEED();
}
TEST_F(BPLINK_BASE_READ, WRONG_P)
{
    // TOY file should still be OR, but should be ok
    std::string base = path + "Wrong.P.assoc";
    std::string out = path + "test";
    std::vector<int> index = {1, 3, 4, 6, 0, 2, -1, 5, -1, -1, -1, 6};
    std::vector<double> barlevels = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    double lower = 5e-5;
    double inter = 0.0001;
    double upper = 0.5;
    double maf_control = 0;
    double maf_case = 0;
    double info_score = 0;
    bool maf_control_filter = false;
    bool maf_case_filter = false;
    bool info_filter = false;
    bool fastscore = false;
    bool no_full = false;
    bool beta = true;
    bool is_index = false;
    bool perform_shrinkage = false;
    std::vector<std::string> feature;
    Region region(feature, 0, 0, false, false);
    Reporter reporter(std::string(path + "LOG"));
    try
    {
        plink->read_base(base, out, index, barlevels, lower, inter, upper,
                         maf_control, maf_case, info_score, maf_control_filter,
                         maf_case_filter, info_filter, fastscore, no_full, beta,
                         is_index, perform_shrinkage, region, reporter);
    }
    catch (const std::runtime_error& e)
    {
        reporter.report(e.what());
        FAIL();
    }
    ASSERT_EQ(plink->num_snps(), 1756);
    // NOTE: the order of SNPs is based on the TOY_TARGET_DATA, not based on the
    // base. Therefore the expected statistic should be 0.08
    auto&& snp = plink->get_snp(0);
    ASSERT_DOUBLE_EQ(snp.stat(), 8.044);
    ASSERT_DOUBLE_EQ(snp.p_value(), 0.548476538760588);
    SUCCEED();
}
// negative coordinates
TEST_F(BPLINK_BASE_READ, WRONG_COORDINATE)
{
    // TOY file should still be OR, but should be ok
    std::string base = path + "Wrong.BP.assoc";
    std::string out = path + "test";
    std::vector<int> index = {1, 3, 4, 6, 0, 2, -1, 5, -1, -1, -1, 6};
    std::vector<double> barlevels = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    double lower = 5e-5;
    double inter = 0.0001;
    double upper = 0.5;
    double maf_control = 0;
    double maf_case = 0;
    double info_score = 0;
    bool maf_control_filter = false;
    bool maf_case_filter = false;
    bool info_filter = false;
    bool fastscore = false;
    bool no_full = false;
    bool beta = true;
    bool is_index = false;
    bool perform_shrinkage = false;
    std::vector<std::string> feature;
    Region region(feature, 0, 0, false, false);
    Reporter reporter(std::string(path + "LOG"));
    try
    {
        plink->read_base(base, out, index, barlevels, lower, inter, upper,
                         maf_control, maf_case, info_score, maf_control_filter,
                         maf_case_filter, info_filter, fastscore, no_full, beta,
                         is_index, perform_shrinkage, region, reporter);
        FAIL();
    }
    catch (const std::runtime_error& e)
    {
        // negative coordinate is very bad as it is almost certain something
        // went very wrong
        reporter.report(e.what());
        SUCCEED();
    }
}
// Test LD Reference sample read and SNP read.
*/
#endif // BIN_PLINK_TEST_H
