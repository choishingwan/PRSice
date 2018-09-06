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

class BPLINK_LOAD_SAMPLE_TARGET : public ::testing::Test
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

TEST_F(BPLINK_LOAD_SAMPLE_TARGET, NO_SELECTION)
{
    std::string keep_file = "";
    std::string remove_file = "";
    bool verbose = true;

    Reporter reporter(std::string(path + "LOG"));
    plink->load_samples(keep_file, remove_file, verbose, reporter);
    ASSERT_EQ(plink->num_sample(), 2000);
}
#endif // BIN_PLINK_TEST_H
