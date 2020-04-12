#ifndef GENOTYPE_TEST_HPP
#define GENOTYPE_TEST_HPP
#include "genotype.hpp"
#include "region.hpp"
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

    void
    sample_generation_check(const std::unordered_set<std::string>& founder_info,
                            const bool in_regress, const bool cal_prs,
                            const bool cal_ld, const bool keep,
                            const bool duplicated,
                            const std::vector<std::string>& input,
                            std::unordered_set<std::string>& processed_sample,
                            std::vector<Sample_ID>& result,
                            std::vector<std::string>& duplicated_sample,
                            size_t& cur_idx, bool remove_sample = true,
                            bool ignore_fid = false, size_t sex_col = +FAM::SEX)
    {
        m_ignore_fid = ignore_fid;
        m_remove_sample = remove_sample;
        const size_t before = duplicated_sample.size();
        std::vector<std::string> token = input;
        gen_sample(+FAM::FID, +FAM::IID, sex_col, +FAM::FATHER, +FAM::MOTHER,
                   cur_idx, founder_info, token[+FAM::PHENOTYPE], token, result,
                   processed_sample, duplicated_sample);
        if (!keep)
        {
            ASSERT_EQ(result.size(), cur_idx);
            // shouldn't update the duplicated sample list
            ASSERT_EQ(before, duplicated_sample.size());
        }
        else if (duplicated)
        {
            ASSERT_EQ(duplicated_sample.size(), before + 1);
        }
        else
        {
            ASSERT_EQ(result.size(), cur_idx + 1);
            if (!m_ignore_fid)
            {
                ASSERT_NE(processed_sample.find(token[+FAM::FID] + m_delim
                                                + token[+FAM::IID]),
                          processed_sample.end());
            }
            else
            {
                ASSERT_NE(processed_sample.find(token[+FAM::IID]),
                          processed_sample.end());
            }
            ASSERT_STREQ(result.back().pheno.c_str(),
                         token[+FAM::PHENOTYPE].c_str());
            ASSERT_EQ(result.back().in_regression, in_regress);
            if (!m_ignore_fid)
                ASSERT_STREQ(result.back().FID.c_str(),
                             token[+FAM::FID].c_str());
            ASSERT_STREQ(result.back().IID.c_str(), token[+FAM::IID].c_str());
            ASSERT_EQ(IS_SET(m_calculate_prs.data(), cur_idx), cal_prs);
            ASSERT_EQ(IS_SET(m_sample_for_ld.data(), cur_idx), cal_ld);
            ASSERT_EQ(before, duplicated_sample.size());
            ++cur_idx;
        }
    }
};
// This should be the simplest of the three genotype class. Test anything that
// doesn't require binaryplink and binarygen
// set_genotype_files


TEST_F(GENOTYPE_BASIC, READ_BASE_FULL)
{
    Reporter reporter(std::string("LOG"), 60, true);
    init_chr();
    m_reporter = &reporter;
    m_snp_selection_list.insert("exclude");
    // wrong loc need to test separately as that is a throw condition
    const double maf_case = 0.04, maf_control = 0.05, info = 0.8,
                 max_thres = 0.5;
    std::vector<std::string> base = {
        "CHR BP RS A1 A2 P STAT MAF INFO MAF_CASE",
        "chr1 1234 exclude A C 0.05 1.96 0.1 0.9 0.05",
        "chr1 1234 normal A C 0.05 1.96 0.1 0.9 0.05",
        "chr1 1234 dup A C 0.05 1.96 0.1 0.9 0.05",
        "chr1 1234 dup A C 0.05 1.96 0.1 0.9 0.05",
        "chrX 1234 sex_chr A C 0.07 1.98 0.1 0.9 0.05",
        "chromosome 1234 wrong_chr A C 0.07 1.98 0.1 0.9 0.05",
        "chr1 1234 filter_control A C 0.05 1.96 0.01 0.9 0.05",
        "chr1 1234 filter_case A C 0.05 1.96 0.1 0.9 0.03",
        "chr1 1234 filter_info A C 0.05 1.96 0.1 0.02 0.05",
        "chr1 1234 p_not_convert A C NA 1.96 0.1 0.9 0.05",
        "chr1 1234 p_exclude A C 0.6 1.96 0.1 0.9 0.05",
        "chr1 1234 stat_not_convert A C 0.05 NA 0.1 0.9 0.05",
        "chr1 1234 negative_stat A C 0.05 -1.96 0.1 0.9 0.05",
        "chr6 1234 region_exclude A C 0.05 -1.96 0.1 0.9 0.05",
        "chr1 1234 ambiguous A T 0.05 1.96 0.1 0.9 0.05"};
    std::vector<size_t> expected(+FILTER_COUNT::MAX, 1);
    // skip header, as we are not using is_index
    expected[+FILTER_COUNT::NUM_LINE] = base.size() - 1;
    // both P and stat shares the same count
    expected[+FILTER_COUNT::NOT_CONVERT] = 2;
    // same for maf case and maf control
    expected[+FILTER_COUNT::MAF] = 2;
    std::ofstream dummy("DUMMY");
    for (auto e : base) { dummy << e << std::endl; }
    dummy.close();
    BaseFile base_file;
    base_file.file_name = "DUMMY";
    base_file.has_column[+BASE_INDEX::CHR] = true;
    base_file.has_column[+BASE_INDEX::BP] = true;
    base_file.has_column[+BASE_INDEX::EFFECT] = true;
    base_file.has_column[+BASE_INDEX::INFO] = true;
    base_file.has_column[+BASE_INDEX::MAF] = true;
    base_file.has_column[+BASE_INDEX::MAF_CASE] = true;
    base_file.has_column[+BASE_INDEX::NONEFFECT] = true;
    base_file.has_column[+BASE_INDEX::P] = true;
    base_file.has_column[+BASE_INDEX::RS] = true;
    base_file.has_column[+BASE_INDEX::STAT] = true;
    base_file.column_index[+BASE_INDEX::CHR] = 0;
    base_file.column_index[+BASE_INDEX::BP] = 1;
    base_file.column_index[+BASE_INDEX::RS] = 2;
    base_file.column_index[+BASE_INDEX::EFFECT] = 3;
    base_file.column_index[+BASE_INDEX::NONEFFECT] = 4;
    base_file.column_index[+BASE_INDEX::P] = 5;
    base_file.column_index[+BASE_INDEX::STAT] = 6;
    base_file.column_index[+BASE_INDEX::MAF] = 7;
    base_file.column_index[+BASE_INDEX::INFO] = 8;
    base_file.column_index[+BASE_INDEX::MAF_CASE] = 9;
    base_file.column_index[+BASE_INDEX::MAX] = 9;
    base_file.is_or = true;
    QCFiltering base_qc;
    base_qc.info_score = info;
    base_qc.maf = maf_control;
    base_qc.maf_case = maf_case;
    PThresholding threshold_info;
    threshold_info.no_full = true;
    threshold_info.fastscore = true;
    threshold_info.bar_levels = {max_thres};
    std::vector<IITree<size_t, size_t>> exclusion_regions;
    Region::generate_exclusion(exclusion_regions, "chr6:1-2000");
    auto [filter_count, dup_idx] =
        read_base(base_file, base_qc, threshold_info, exclusion_regions);
    std::remove("DUMMY");
    // 2 because we have both normal and dup in it
    ASSERT_EQ(m_existed_snps.size(), 2);
    // normal must be here
    ASSERT_TRUE(m_existed_snps_index.find("normal")
                != m_existed_snps_index.end());
    ASSERT_EQ(dup_idx.size(), 1);
    ASSERT_TRUE(dup_idx.find("dup") != dup_idx.end());
    ASSERT_EQ(expected.size(), filter_count.size());
    for (size_t i = 0; i < expected.size(); ++i)
    { ASSERT_EQ(filter_count[i], expected[i]); }
    dummy.open("DUMMY");
    dummy << "CHR BP RS A1 A2 P STAT MAF INFO MAF_CASE" << std::endl;
    dummy << "chr1 -1234 negative_loc A C 0.05 1.96 0.1 0.9 0.05" << std::endl;
    dummy.close();
    // should fail due to invalid coordinates
    try
    {
        read_base(base_file, base_qc, threshold_info, exclusion_regions);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    // invalid p-value
    dummy.open("DUMMY");
    dummy << "CHR BP RS A1 A2 P STAT MAF INFO MAF_CASE" << std::endl;
    dummy << "chr1 1234 invalid_p A C -0.05 1.96 0.1 0.9 0.05" << std::endl;
    dummy.close();
    // should fail due to invalid coordinates
    try
    {
        read_base(base_file, base_qc, threshold_info, exclusion_regions);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    std::remove("LOG");
    std::remove("DUMMY");
}
// init_chr
// chr_code_check
// load_snp_list
// load_ref


TEST_F(GENOTYPE_BASIC, INIT_SAMPLE_VECTOR)
{
    // should fail as m_unfiltered_sample_ct is 0
    try
    {
        init_sample_vectors();
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    // should sucess
    m_unfiltered_sample_ct = 1025;
    init_sample_vectors();
    ASSERT_EQ(m_num_male, 0);
    ASSERT_EQ(m_num_female, 0);
    ASSERT_EQ(m_num_ambig, 0);
    ASSERT_EQ(m_num_non_founder, 0);
    ASSERT_EQ(m_sample_ct, 0);
    ASSERT_EQ(m_founder_ct, 0);
    ASSERT_TRUE(!m_sample_for_ld.empty());
    ASSERT_TRUE(!m_calculate_prs.empty());
}


TEST_F(GENOTYPE_BASIC, SAMPLE_GENERATION)
{
    /* Check the following:
     * 1. Founder
     * 2. Non-founder
     * 3. Non-founder but keep
     * 4. Remove sample
     * 5. Keep sample
     * 6. Sex isn't too important here. Will just make sure we don't fail when
     * sex is not provided
     */
    // we will first assume fam format for easier handling
    // for bgen, the founder info structure should be empty
    //{"FID", "IID", "Dad",
    //"Mum", "Sex", "Pheno"}
    m_delim = " ";
    std::unordered_set<std::string> founder_info = {"FAM1 DAD1", "FAM1 MUM1",
                                                    "FAM2 MUM2", "FAM3 DAD1"};
    m_sample_selection_list.clear();
    m_sample_selection_list.insert("REMOVE1 REMOVE1");
    m_sample_selection_list.insert("REMOVE2 REMOVE2");
    m_remove_sample = true;
    std::vector<Sample_ID> result;
    std::unordered_set<std::string> processed_sample;
    std::vector<std::string> duplicated_sample;
    size_t cur_idx = 0;
    // must first initialize the vector
    m_unfiltered_sample_ct = 65;
    init_sample_vectors();
    const bool in_regress = true, cal_prs = true, cal_ld = true, keep = true,
               duplicated = true;
    // now start testing, this is a founder;
    sample_generation_check(
        founder_info, in_regress, cal_prs, cal_ld, keep, !duplicated,
        std::vector<std::string> {"ID1", "ID1", "0", "0", "1", "1"},
        processed_sample, result, duplicated_sample, cur_idx);
    // now non-founder, but not using it for regression
    sample_generation_check(
        founder_info, !in_regress, cal_prs, !cal_ld, keep, !duplicated,
        std::vector<std::string> {"FAM1", "BOY1", "DAD1", "MUM1", "1", "0"},
        processed_sample, result, duplicated_sample, cur_idx);
    // now if we change it to keep non_founder
    m_keep_nonfounder = true;
    // will use it for everything except LD
    sample_generation_check(
        founder_info, in_regress, cal_prs, !cal_ld, keep, !duplicated,
        std::vector<std::string> {"FAM2", "GIRL1", "DAD2", "MUM2", "2", "NA"},
        processed_sample, result, duplicated_sample, cur_idx);
    // another non-founder, but both parents are not found in the system
    // so should treat as a normal founders
    sample_generation_check(
        founder_info, in_regress, cal_prs, cal_ld, keep, !duplicated,
        std::vector<std::string> {"FAM4", "GIRL2", "DAD4", "MUM4", "2", "1"},
        processed_sample, result, duplicated_sample, cur_idx);
    // should be ok if we are in the same family as a founder but not the same
    // parent (e.g. skipped a generation?)
    sample_generation_check(
        founder_info, in_regress, cal_prs, cal_ld, keep, !duplicated,
        std::vector<std::string> {"FAM3", "BOY2", "0", "MUM2", "1", "1"},
        processed_sample, result, duplicated_sample, cur_idx);
    // and if we have any parent founder, we are defined as non-founder
    sample_generation_check(
        founder_info, in_regress, cal_prs, !cal_ld, keep, !duplicated,
        std::vector<std::string> {"FAM1", "BOY3", "0", "MUM1", "1", "1"},
        processed_sample, result, duplicated_sample, cur_idx);
    // sample removal
    sample_generation_check(
        founder_info, in_regress, cal_prs, !cal_ld, !keep, !duplicated,
        std::vector<std::string> {"REMOVE1", "REMOVE1", "0", "MUM1", "1", "1"},
        processed_sample, result, duplicated_sample, cur_idx);
    // must match both IID and FID to remove it
    sample_generation_check(
        founder_info, in_regress, cal_prs, cal_ld, keep, !duplicated,
        std::vector<std::string> {"REMOVE1", "BOY4", "0", "MUM1", "1", "1"},
        processed_sample, result, duplicated_sample, cur_idx);
    // now switch to keep
    sample_generation_check(
        founder_info, in_regress, cal_prs, cal_ld, keep, !duplicated,
        std::vector<std::string> {"REMOVE2", "REMOVE2", "0", "MUM1", "1", "1"},
        processed_sample, result, duplicated_sample, cur_idx, false);
    // test ignore FID
    sample_generation_check(
        founder_info, in_regress, cal_prs, cal_ld, keep, !duplicated,
        std::vector<std::string> {"ABC", "ID2", "0", "0", "1", "1"},
        processed_sample, result, duplicated_sample, cur_idx, true, true);
    // now check duplicated sample
    sample_generation_check(
        founder_info, in_regress, cal_prs, cal_ld, keep, duplicated,
        std::vector<std::string> {"ID1", "ID1", "0", "0", "1", "1"},
        processed_sample, result, duplicated_sample, cur_idx);
    // check when we have no sex information
    sample_generation_check(
        founder_info, in_regress, cal_prs, cal_ld, keep, !duplicated,
        std::vector<std::string> {"NO_SEX", "Ambig", "0", "0", "1", "3"},
        processed_sample, result, duplicated_sample, cur_idx, true, false,
        ~size_t(0));
}

#endif // GENOTYPE_TEST_HPP
