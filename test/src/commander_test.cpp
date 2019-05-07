#ifndef COMMANDER_TEST_H
#define COMMANDER_TEST_H
#include "commander.hpp"
#include "global.hpp"
#include "reporter.hpp"
#include "storage.hpp"
#include "gtest/gtest.h"

TEST(COMMANDER_BASIC, INIT)
{
    Commander commander;
    std::vector<size_t> base_col_index = commander.index();
    std::vector<bool> base_has_col = commander.has_col();
    ASSERT_EQ(base_col_index.size(), base_has_col.size());
    ASSERT_EQ(base_col_index.size(), +BASE_INDEX::MAX + 1);
    // everything else should be empty or set to default
    ASSERT_TRUE(commander.base_name().empty());
    ASSERT_FALSE(commander.all_scores());
    ASSERT_FALSE(commander.use_inter());
    ASSERT_FALSE(commander.fastscore());
    ASSERT_FALSE(commander.no_clump());
    ASSERT_FALSE(commander.no_full());
    ASSERT_FALSE(commander.no_regress());
    ASSERT_TRUE(commander.background().empty());
    ASSERT_FALSE(commander.genome_wide_background());
    ASSERT_FALSE(commander.ignore_fid());
    ASSERT_FALSE(commander.nonfounders());
    ASSERT_FALSE(commander.is_index());
    ASSERT_FALSE(commander.keep_ambig());
    ASSERT_FALSE(commander.logit_perm());
    ASSERT_FALSE(commander.non_cumulate());
    ASSERT_FALSE(commander.pearson());
    ASSERT_FALSE(commander.print_snp());
    ASSERT_FALSE(commander.beta());
    ASSERT_FALSE(commander.hard_coded());
    ASSERT_FALSE(commander.use_ref_maf());
    ASSERT_FALSE(commander.use_ref());
    ASSERT_FALSE(commander.beta());
    ASSERT_TRUE(commander.get_cov_file().empty());
    ASSERT_TRUE(commander.get_cov_name().empty());
    ASSERT_TRUE(commander.get_cov_index().empty());
    ASSERT_TRUE(commander.get_factor_cov_index().empty());
    ASSERT_TRUE(commander.bed().empty());
    ASSERT_TRUE(commander.feature().empty());
    ASSERT_TRUE(commander.pheno_col().empty());
    ASSERT_TRUE(commander.pheno_file().empty());
    ASSERT_TRUE(commander.bar_levels().empty());
    ASSERT_TRUE(commander.prevalence().empty());
    ASSERT_TRUE(commander.is_binary().empty());
    ASSERT_TRUE(commander.target_name().empty());
    ASSERT_TRUE(commander.target_list().empty());
    ASSERT_TRUE(commander.keep_sample_file().empty());
    ASSERT_TRUE(commander.remove_sample_file().empty());
    ASSERT_STREQ(commander.target_type().c_str(), "bed");
    ASSERT_TRUE(commander.gtf().empty());
    ASSERT_TRUE(commander.msigdb().empty());
    ASSERT_TRUE(commander.background().empty());
    ASSERT_TRUE(commander.snp_set().empty());
    ASSERT_TRUE(commander.base_name().empty());
    ASSERT_STREQ(commander.delim().c_str(), " ");
    ASSERT_STREQ(commander.out().c_str(), "PRSice");
    ASSERT_TRUE(commander.exclusion_range().empty());
    ASSERT_TRUE(commander.exclude_file().empty());
    ASSERT_TRUE(commander.extract_file().empty());
    ASSERT_TRUE(commander.ref_name().empty());
    ASSERT_TRUE(commander.ref_list().empty());
    ASSERT_STREQ(commander.ref_type().c_str(), "bed");
    ASSERT_TRUE(commander.ref_keep_file().empty());
    ASSERT_TRUE(commander.ref_remove_file().empty());
    double dummy=-1;
    ASSERT_FALSE(commander.target_maf(dummy));
    ASSERT_FALSE(commander.target_geno(dummy));
    ASSERT_FALSE(commander.target_info(dummy));
    ASSERT_FALSE(commander.target_hard_threshold(dummy));
    ASSERT_FALSE(commander.ref_maf(dummy));
    ASSERT_FALSE(commander.ref_geno(dummy));
    ASSERT_FALSE(commander.ref_info(dummy));
    ASSERT_FALSE(commander.ref_hard_threshold(dummy));
    ASSERT_FALSE(commander.proxy(dummy));
    ASSERT_DOUBLE_EQ(commander.base_info_score(),0.0);
    ASSERT_DOUBLE_EQ(commander.clump_p(), 1.0);
    ASSERT_DOUBLE_EQ(commander.clump_r2(), 0.1);
    ASSERT_DOUBLE_EQ(commander.lower(), 5e-8);
    ASSERT_DOUBLE_EQ(commander.inter(), 0.00005);
    ASSERT_DOUBLE_EQ(commander.upper(), 0.5);
    // we will use the parameter if memory is not provided, otherwise,
    // we will return the actual memory allowed
    ASSERT_DOUBLE_EQ(commander.max_memory(1.0),1.0);
    ASSERT_DOUBLE_EQ(commander.max_memory(2.0),2.0);
    size_t dummy_int;
    ASSERT_FALSE(commander.num_perm(dummy_int));
    ASSERT_FALSE(commander.set_perm(dummy_int));
    ASSERT_EQ(commander.get_missing_score(), MISSING_SCORE::MEAN_IMPUTE);
    ASSERT_EQ(commander.get_score(), SCORING::AVERAGE);
    ASSERT_EQ(commander.model(), MODEL::ADDITIVE);
    ASSERT_EQ(commander.clump_dist(), 1000000);
    ASSERT_EQ(commander.thread(), 1);
    ASSERT_EQ(commander.window_3(), 0);
    ASSERT_EQ(commander.window_5(), 0);
}

TEST(COMMANDER_BASIC, USAGE)
{
    Commander commander;
    Reporter reporter(std::string(path + "LOG"));
    int argc=2;
    char name[7], help[7];
    strcpy(name, "PRSice");
    strcpy(help, "--help");
    char* argv[2] ={name, help};
    try {
        ASSERT_FALSE(commander.init(argc, argv, reporter));
    } catch (...) {
        FAIL();
    }
}


TEST(COMMANDER_BASIC, NO_ARG)
{
    Commander commander;
    Reporter reporter(std::string(path + "LOG"));
    int argc=1;
    std::string name = "PRSice";
    char name_c[7];
    strcpy(name_c, name.c_str());
    char* argv[1] ={name_c};
    try {
        commander.init(argc, argv, reporter);
        FAIL();
    } catch (...) {
        SUCCEED();
    }
}
#endif // COMMANDER_TEST_H
