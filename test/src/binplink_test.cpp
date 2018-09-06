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

class BPLINK_GEN_SAMPLE_TARGET : public ::testing::Test
{
protected:
    BinaryPlink* plink = nullptr;
    void SetUp() override
    {
        std::string file_list;
        std::string file = std::string(path + "TOY_TARGET_DATA");
        uint32_t thread = 1;
        bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
             is_ref = false;
        Reporter reporter(std::string(path + "LOG"));
        plink = new BinaryPlink(file_list, file, thread, ignore_fid,
                                keep_nonfounder, keep_ambig, is_ref, reporter);
    }
    void TearDown() override { delete plink; }
};
TEST_F(BPLINK_GEN_SAMPLE_TARGET, NO_SELECTION)
{
    std::string keep_file = "";
    std::string remove_file = "";
    bool verbose = true;

    Reporter reporter(std::string(path + "LOG"));
    plink->load_samples(keep_file, remove_file, verbose, reporter);
    ASSERT_EQ(plink->num_sample(), 2000);
}
TEST_F(BPLINK_GEN_SAMPLE_TARGET, KEEP_SAMPLE)
{
    std::string keep_file = path + "keep";
    std::ofstream keep;
    keep.open(keep_file.c_str());
    // We are using toy sample, so the sample naming convention should be
    // CAS_XXX CONT_XXX where XXX range from 1 to 1000
    for (size_t i = 123; i < 234; ++i) {
        keep << "CAS_" << i << "\t"
             << "CAS_" << i << std::endl;
        keep << "CONT_" << i << "\t"
             << "CONT_" << i << std::endl;
    }

    keep.close();
    std::string remove_file = "";
    bool verbose = true;

    Reporter reporter(std::string(path + "LOG"));
    plink->load_samples(keep_file, remove_file, verbose, reporter);
    ASSERT_EQ(plink->num_sample(), 222);
}
TEST_F(BPLINK_GEN_SAMPLE_TARGET, REMOVE_SAMPLE)
{
    std::string keep_file = "";
    std::string remove_file = path + "remove";
    std::ofstream remove;
    remove.open(remove_file.c_str());
    // We are using toy sample, so the sample naming convention should be
    // CAS_XXX CONT_XXX where XXX range from 1 to 1000
    for (size_t i = 321; i < 432; ++i) {
        remove << "CAS_" << i << "\t"
               << "CAS_" << i << std::endl;
        remove << "CONT_" << i << "\t"
               << "CONT_" << i << std::endl;
    }

    remove.close();
    bool verbose = true;

    Reporter reporter(std::string(path + "LOG"));
    plink->load_samples(keep_file, remove_file, verbose, reporter);
    ASSERT_EQ(plink->num_sample(), 1778);
}
TEST_F(BPLINK_GEN_SAMPLE_TARGET, KEEP_SAMPLE_IID)
{
    std::string keep_file = path + "keep";
    std::ofstream keep;
    keep.open(keep_file.c_str());
    // We are using toy sample, so the sample naming convention should be
    // CAS_XXX CONT_XXX where XXX range from 1 to 1000
    for (size_t i = 123; i < 234; ++i) {
        keep << "CAS_" << i << std::endl;
        keep << "CONT_" << i << std::endl;
    }

    keep.close();
    std::string remove_file = "";
    bool verbose = true;
    // we have set ignore_fid to false, therefore when only IID is provided,
    // PRSice should error out
    Reporter reporter(std::string(path + "LOG"));
    try
    {
        plink->load_samples(keep_file, remove_file, verbose, reporter);
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
    std::string file =
        std::string(path + "TOY_TARGET_DATA," + path + "TOY_TARGET_DATA.fam");
    uint32_t thread = 1;
    bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
         is_ref = false;
    Reporter reporter(std::string(path + "LOG"));
    BinaryPlink plinkBinary(file_list, file, thread, ignore_fid,
                            keep_nonfounder, keep_ambig, is_ref, reporter);
    plinkBinary.load_samples("", "", true, reporter);
    ASSERT_EQ(plinkBinary.num_sample(), 2000);
}
TEST(BPLINK_SAMPLE_CHECK, DUPLICATE_SAMPLE)
{
    std::string file_list;
    std::string file = std::string(path + "TOY_TARGET_DATA," + path
                                   + "TOY_TARGET_DATA.dup.fam");
    uint32_t thread = 1;
    bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
         is_ref = false;
    Reporter reporter(std::string(path + "LOG"));
    BinaryPlink plinkBinary(file_list, file, thread, ignore_fid,
                            keep_nonfounder, keep_ambig, is_ref, reporter);
    try
    {
        plinkBinary.load_samples("", "", true, reporter);
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
    std::string file = std::string(path + "TOY_TARGET_DATA," + path
                                   + "TOY_TARGET_DATA.founder.fam");
    uint32_t thread = 1;
    bool ignore_fid = false, keep_ambig = false, keep_nonfounder = false,
         is_ref = false;
    Reporter reporter(std::string(path + "LOG"));
    BinaryPlink plinkBinary(file_list, file, thread, ignore_fid,
                            keep_nonfounder, keep_ambig, is_ref, reporter);
    plinkBinary.load_samples("", "", true, reporter);
    // we still keep the non-founders, just not using them for regression
    ASSERT_EQ(plinkBinary.num_sample(), 2000);
    int sum_founder = 0;
    for (size_t i = 0; i < 2000; ++i) {
        sum_founder += plinkBinary.is_founder(i);
    }
    // NOTE: This was intended to be 859. But our dummy has two of the
    // founders set to have a different FID, the script correctly account for
    // samples in these family to be a fonuder, because while the IID of
    // paternal and maternal is here, that sample is from a different family
    ASSERT_EQ(sum_founder, 862);
}
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
TEST(BPLINK_GEN_SNP, SEQ_INPUT) {}
#endif // BIN_PLINK_TEST_H
