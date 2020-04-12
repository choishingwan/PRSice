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
