#ifndef BIN_GEN_TEST_H
#define BIN_GEN_TEST_H
#include "bgen_lib.hpp"
#include "binarygen.hpp"
#include "global.hpp"
#include "gtest/gtest.h"
#include <string>

// BinaryGen file can be rather fragile with sample read
// if we failed to detect the sample file, but the file is in pseudo sample file
// format, we will get an error as the number of read samples and number of
// sample indicated in the bgen file differs

// we assume the bgen lib is valid and have done the required unit tests

void check_sample_file(const std::vector<std::string>& input, bool is_sample)
{
    std::ofstream dummy("DUMMY");
    for (auto&& i : input) { dummy << i << std::endl; }
    dummy.close();
    ASSERT_EQ(BinaryGen::check_is_sample_format("DUMMY"), is_sample);
    std::remove("DUMMY");
}
TEST(BINARY_GEN, CHECK_SAMPLE_FILE)
{
    // sample format
    const bool is_sample = true;
    // minimum sample file
    check_sample_file(
        std::vector<std::string> {"ID1 ID2 missing", "0 0 0", "1 1 1"},
        is_sample);
    // with sex info
    check_sample_file(
        std::vector<std::string> {"ID1 ID2 missing Sex", "0 0 0 D", "1 1 1 F"},
        is_sample);
    // malform format
    check_sample_file(
        std::vector<std::string> {"ID1 ID2 missing ", "0 0", "1 1 1"},
        !is_sample);
    // forgot missing
    check_sample_file(std::vector<std::string> {"ID1 ID2", "0 0", "1 1"},
                      !is_sample);
    // unknown format
    check_sample_file(std::vector<std::string> {"ID1 ID2 missing sex Pheno",
                                                "0 0 0 D N", "S1 S1 0 F -1"},
                      !is_sample);
    // missing replaced by sex, caught by looking at the second row
    check_sample_file(
        std::vector<std::string> {"ID1 ID2 sex", "0 0 D", "S1 S1 F"},
        !is_sample);
    // rather long format, but is valid
    check_sample_file(
        std::vector<std::string> {
            "ID_1 ID_2 missing sex category binary positive disgrete",
            "0 0 0 D C B P D", "1 1 1 F 10.0 1 2 Hi"},
        is_sample);
    // this is a fam file
    check_sample_file(
        std::vector<std::string> {"ID1 ID2 0 0 1 Pheno", "ID2 ID3 0 0 2 Pheno"},
        !is_sample);
    // one sample fam file?
    check_sample_file(std::vector<std::string> {"ID2 ID3 0 0 2 Pheno"},
                      !is_sample);
    // phenotype file
    check_sample_file(
        std::vector<std::string> {"ID1 ID2 Pheno", "ID2 ID3 Pheno"},
        !is_sample);
}

// test sample file reading

class BGEN_TEST : public ::BinaryGen, public ::testing::Test
{
public:
    std::ofstream dummy;
    BGEN_TEST()
        : BinaryGen(GenoFile("DUMMY,dummy.sample"), Phenotype(), " ",
                    new Reporter("LOG", 60, true))
    {
        dummy.open("DUMMY.bgen", std::ofstream::binary);
    }
    ~BGEN_TEST()
    {
        if (m_reporter != nullptr) { delete m_reporter; }
        std::remove("DUMMY.bgen");
        std::remove("dummy.sample");
    }
    void gen_bgen(uint32_t number_of_snp_blocks, uint32_t number_of_samples,
                  std::string free_data, uint32_t flags)
    {
        genfile::bgen::Context context;
        context.number_of_variants = number_of_snp_blocks;
        context.number_of_samples = number_of_samples;
        context.free_data = free_data;
        context.flags = flags;
        genfile::bgen::write_header_block(dummy, context);
        dummy.close();
    }

    void check_sample(const std::vector<Sample_ID>& result,
                      const std::string& fid, const std::string& iid,
                      const bool keep, size_t& idx, size_t& bit_idx)
    {
        ASSERT_TRUE(result.size() > idx);
        if (keep)
        {
            ASSERT_STREQ(result[idx].FID.c_str(), fid.c_str());
            ASSERT_STREQ(result[idx].IID.c_str(), iid.c_str());
            ASSERT_TRUE(result[idx].in_regression);
            ++idx;
        }
        ASSERT_EQ(IS_SET(m_sample_for_ld.data(), bit_idx), keep);
        ASSERT_EQ(IS_SET(m_calculate_prs.data(), bit_idx), keep);
        ++bit_idx;
    }
};
TEST_F(BGEN_TEST, LOAD_SAMPLE_SFILE)
{
    const size_t num_snp = 10;
    // now generate sample file
    std::string header = "ID_1 ID_2 missing\n0 0 0";
    // not as complicated as we don't have non-founders with bgen format
    std::vector<std::string> expected = {"ID1 ID1 1", // normal
                                         "FAM1 DAD1 1", "REMOVE1 REMOVE1 1",
                                         "FAM2 GIRL1"};
    m_sample_selection_list.insert("REMOVE1 REMOVE1");
    std::ofstream sample("dummy.sample");
    sample << header << std::endl;
    for (auto e : expected) { sample << e << std::endl; }
    sample.close();
    // now generate the bgen header file
    gen_bgen(num_snp, static_cast<uint32_t>(expected.size()), "Testing data",
             4294967295u);
    get_context(0);
    auto result = gen_sample_vector();
    const bool keep = true;
    size_t idx = 0, bit_idx = 0;
    check_sample(result, "ID1", "ID1", keep, idx, bit_idx);
    check_sample(result, "FAM1", "DAD1", keep, idx, bit_idx);
    check_sample(result, "REMOVE1", "REMOVE1", !keep, idx, bit_idx);
    check_sample(result, "FAM2", "GIRL1", keep, idx, bit_idx);
}
TEST_F(BGEN_TEST, LOAD_SAMPLE_SFILE_WITH_SEX)
{
    const size_t num_snp = 10;
    // now generate sample file
    std::string header = "ID_1 ID_2 missing Pheno Sex\n0 0 0 C D";
    // not as complicated as we don't have non-founders with bgen format
    std::vector<std::string> expected = {
        "ID1 ID1 1 0.1 F", // normal
        "FAM1 DAD1 1 0.5 M", "REMOVE1 REMOVE1 1 1.23 M", "FAM2 GIRL1 1 0.19 F"};
    m_sample_selection_list.insert("REMOVE1 REMOVE1");
    std::ofstream sample("dummy.sample");
    sample << header << std::endl;
    for (auto e : expected) { sample << e << std::endl; }
    sample.close();
    // now generate the bgen header file
    gen_bgen(num_snp, static_cast<uint32_t>(expected.size()), "Testing data",
             4294967295u);
    get_context(0);
    auto result = gen_sample_vector();
    const bool keep = true;
    size_t idx = 0, bit_idx = 0;
    check_sample(result, "ID1", "ID1", keep, idx, bit_idx);
    check_sample(result, "FAM1", "DAD1", keep, idx, bit_idx);
    check_sample(result, "REMOVE1", "REMOVE1", !keep, idx, bit_idx);
    check_sample(result, "FAM2", "GIRL1", keep, idx, bit_idx);
}
TEST_F(BGEN_TEST, LOAD_SAMPLE_SFILE_IGNORE_FID)
{
    m_ignore_fid = true;
    const size_t num_snp = 10;
    // now generate sample file
    std::string header = "ID_1 ID_2 missing Pheno Sex\n0 0 0 C D";
    // not as complicated as we don't have non-founders with bgen format
    std::vector<std::string> expected = {
        "ID1 ID1 1 0.1 F", // normal
        "FAM1 DAD1 1 0.5 M", "REMOVE1 REMOVE1 1 1.23 M", "FAM2 GIRL1 1 0.19 F"};
    m_sample_selection_list.insert("REMOVE1");
    std::ofstream sample("dummy.sample");
    sample << header << std::endl;
    for (auto e : expected) { sample << e << std::endl; }
    sample.close();
    // now generate the bgen header file
    gen_bgen(num_snp, static_cast<uint32_t>(expected.size()), "Testing data",
             4294967295u);
    get_context(0);
    auto result = gen_sample_vector();
    const bool keep = true;
    size_t idx = 0, bit_idx = 0;
    check_sample(result, "", "ID1", keep, idx, bit_idx);
    check_sample(result, "", "DAD1", keep, idx, bit_idx);
    check_sample(result, "", "REMOVE1", !keep, idx, bit_idx);
    check_sample(result, "", "GIRL1", keep, idx, bit_idx);
}
TEST_F(BGEN_TEST, LOAD_SAMPLE_PHENO_FILE_NO_HEADER)
{
    const size_t num_snp = 10;
    // now generate sample file
    // not as complicated as we don't have non-founders with bgen format
    std::vector<std::string> expected = {
        "ID1 ID1 0.1", // normal
        "FAM1 DAD1 0.5", "REMOVE1 REMOVE1 1.23", "FAM2 GIRL1 0.19"};
    m_sample_selection_list.insert("REMOVE1 REMOVE1");
    std::ofstream sample("dummy.sample");
    for (auto e : expected) { sample << e << std::endl; }
    sample.close();
    // now generate the bgen header file
    gen_bgen(num_snp, static_cast<uint32_t>(expected.size()), "Testing data",
             4294967295u);
    get_context(0);
    auto result = gen_sample_vector();
    const bool keep = true;
    size_t idx = 0, bit_idx = 0;
    check_sample(result, "ID1", "ID1", keep, idx, bit_idx);
    check_sample(result, "FAM1", "DAD1", keep, idx, bit_idx);
    check_sample(result, "REMOVE1", "REMOVE1", !keep, idx, bit_idx);
    check_sample(result, "FAM2", "GIRL1", keep, idx, bit_idx);
}
TEST_F(BGEN_TEST, LOAD_SAMPLE_PHENO_FILE_HEADER)
{
    const size_t num_snp = 10;
    // now generate sample file
    // not as complicated as we don't have non-founders with bgen format
    std::string header = "FID IID Pheno";
    std::vector<std::string> expected = {
        "ID1 ID1 0.1", // normal
        "FAM1 DAD1 0.5", "REMOVE1 REMOVE1 1.23", "FAM2 GIRL1 0.19"};
    m_sample_selection_list.insert("REMOVE1 REMOVE1");
    std::ofstream sample("dummy.sample");
    sample << header << std::endl;
    for (auto e : expected) { sample << e << std::endl; }
    sample.close();
    // now generate the bgen header file
    gen_bgen(num_snp, static_cast<uint32_t>(expected.size()), "Testing data",
             4294967295u);
    get_context(0);
    auto result = gen_sample_vector();
    const bool keep = true;
    size_t idx = 0, bit_idx = 0;
    check_sample(result, "ID1", "ID1", keep, idx, bit_idx);
    check_sample(result, "FAM1", "DAD1", keep, idx, bit_idx);
    check_sample(result, "REMOVE1", "REMOVE1", !keep, idx, bit_idx);
    check_sample(result, "FAM2", "GIRL1", keep, idx, bit_idx);
}
TEST_F(BGEN_TEST, LOAD_SAMPLE_PHENO_FILE_IGNORE_FID)
{
    const size_t num_snp = 10;
    // now generate sample file
    // not as complicated as we don't have non-founders with bgen format
    m_ignore_fid = true;
    // with ignore fid and phenotype file, we will use the first column as the
    // ID
    std::vector<std::string> expected = {
        "ID1 ID1 0.1", // normal
        "FAM1 DAD1 0.5", "REMOVE1 REMOVE1 1.23", "FAM2 GIRL1 0.19"};
    m_sample_selection_list.insert("REMOVE1");
    std::ofstream sample("dummy.sample");
    for (auto e : expected) { sample << e << std::endl; }
    sample.close();
    // now generate the bgen header file
    gen_bgen(num_snp, static_cast<uint32_t>(expected.size()), "Testing data",
             4294967295u);
    get_context(0);
    auto result = gen_sample_vector();
    const bool keep = true;
    size_t idx = 0, bit_idx = 0;
    check_sample(result, "", "ID1", keep, idx, bit_idx);
    check_sample(result, "", "FAM1", keep, idx, bit_idx);
    check_sample(result, "", "REMOVE1", !keep, idx, bit_idx);
    check_sample(result, "", "FAM2", keep, idx, bit_idx);
}
#endif
