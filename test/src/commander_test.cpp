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
    ASSERT_FALSE(commander.all_scores());
    ASSERT_FALSE(commander.use_inter());
    ASSERT_STREQ(commander.delim().c_str(), " ");
    ASSERT_STREQ(commander.out().c_str(), "PRSice");
    ASSERT_TRUE(commander.exclusion_range().empty());
    ASSERT_TRUE(commander.exclude_file().empty());
    ASSERT_TRUE(commander.extract_file().empty());
    // we will use the parameter if memory is not provided, otherwise,
    // we will return the actual memory allowed
    ASSERT_DOUBLE_EQ(commander.max_memory(1.0), 1.0);
    ASSERT_DOUBLE_EQ(commander.max_memory(2.0), 2.0);
}

class mockCommander : public Commander
{
public:
    static std::vector<std::string>
    transform_covariate(const std::string& cov_in)
    {
        return Commander::transform_covariate(cov_in);
    }
    bool check_parse_unit_value(const std::string& input, const std::string& c,
                                const size_t default_power, size_t& target,
                                bool memory = false)
    {
        return parse_unit_value(input, c, default_power, target, memory);
    }

    static bool find_first_end_wrapper(const std::string_view& cov,
                                       const size_t idx, size_t& res)
    {
        try
        {
            res = find_first_end(cov, idx);
            return true;
        }
        catch (const std::runtime_error&)
        {
            return false;
        }
    }
    static bool parse_range_wrapper(std::string_view cov,
                                    std::vector<size_t>& res)
    {
        try
        {
            res = parse_range(cov);
            return true;
        }
        catch (std::runtime_error&)
        {
            return false;
        }
    }
    static bool get_range_wrapper(std::string_view cov, size_t start,
                                  size_t end, std::vector<size_t>& res)
    {
        try
        {
            res = get_range(cov, start, end);
            return true;
        }
        catch (std::runtime_error&)
        {
            return false;
        }
    }
    static bool
    update_covariate_ranges_wrapper(std::vector<std::string>& result,
                                    std::vector<size_t> ranges)
    {
        try
        {
            update_covariate_range(ranges, result);
            return true;
        }
        catch (...)
        {
            return false;
        }
    }
    static bool transform_wrapper(const std::string& str,
                                  std::vector<std::string>& result)
    {
        try
        {
            result = transform_covariate(str);
            return true;
        }
        catch (...)
        {
            return false;
        }
    }
    bool parse_command_wrapper(const std::string& command)
    {
        bool early_terminate = false;
        return parse_command_wrapper(command, early_terminate);
    }
    bool parse_command_wrapper(const std::string& command,
                               bool& early_terminate)
    {
        Reporter reporter(std::string("LOG"), 60, true);
        std::vector<std::string> argv_str = misc::split("PRSice " + command);
        std::vector<char*> cstrings;
        cstrings.reserve(argv_str.size());
        for (size_t i = 0; i < argv_str.size(); ++i)
        { cstrings.push_back(const_cast<char*>(argv_str[i].c_str())); }
        int argc = static_cast<int>(argv_str.size());
        try
        {
            early_terminate = false;
            // return false if error
            return !init(argc, &cstrings[0], early_terminate, reporter);
        }
        catch (...)
        {
            // error = false
            return false;
        }
    }

    bool no_default() const { return m_user_no_default; }
    int32_t max_thread() { return maximum_thread(); }
};

TEST(COMMAND_PARSING, USAGE)
{
    mockCommander commander;
    bool early_terminate = false;
    ASSERT_FALSE(commander.parse_command_wrapper("--help", early_terminate));
    ASSERT_TRUE(early_terminate);
    ASSERT_FALSE(commander.parse_command_wrapper("-h", early_terminate));
    ASSERT_TRUE(early_terminate);
    // this is a throw error
    ASSERT_FALSE(commander.parse_command_wrapper("", early_terminate));
    ASSERT_FALSE(early_terminate);
    // this should fail, as ? is reserved for invalid operators
    ASSERT_FALSE(commander.parse_command_wrapper("-?", early_terminate));
    ASSERT_FALSE(early_terminate);
    // version check should be similar to --help
    ASSERT_FALSE(commander.parse_command_wrapper("-v", early_terminate));
    ASSERT_TRUE(early_terminate);
    ASSERT_FALSE(commander.parse_command_wrapper("--version", early_terminate));
    ASSERT_TRUE(early_terminate);
}

void check_bar_threshold(const std::string& command,
                         const std::vector<double>& expected,
                         const bool expect_fail)
{
    mockCommander commander;
    if (expect_fail) { ASSERT_FALSE(commander.parse_command_wrapper(command)); }
    else
    {
        ASSERT_TRUE(commander.parse_command_wrapper(command));
        auto res = commander.get_p_threshold();
        ASSERT_EQ(expected.size(), res.bar_levels.size());
        for (size_t i = 0; i < res.bar_levels.size(); ++i)
        { ASSERT_DOUBLE_EQ(res.bar_levels[i], expected[i]); }
    }
}
TEST(COMMAND_PARSING, BAR_LEVELS_VALID)
{
    // valid
    check_bar_threshold("--bar-levels 0.1,0.2,0.3,0.4,0.5",
                        std::vector<double> {0.1, 0.2, 0.3, 0.4, 0.5}, false);
    // we have not deal with duplicates yet
    check_bar_threshold("--bar-levels 0.1,0.2,0.3,0.3,0.4,0.5",
                        std::vector<double> {0.1, 0.2, 0.3, 0.3, 0.4, 0.5},
                        false);
    // Have not sorted either
    check_bar_threshold("--bar-levels 0.5,0.2,0.3,0.3,0.4,0.1",
                        std::vector<double> {0.5, 0.2, 0.3, 0.3, 0.4, 0.1},
                        false);
    // supposed to fail but init didn't check for these either
    // negative number is no no
    check_bar_threshold("--bar-levels 0.1,-0.2,0.3,0.4,0.5",
                        std::vector<double> {0.1, -0.2, 0.3, 0.4, 0.5}, false);
    // No zero surely?
    check_bar_threshold("--bar-levels 0,0.2,0.3,0.3,0.4,0.5",
                        std::vector<double> {0, 0.2, 0.3, 0.3, 0.4, 0.5},
                        false);
    // Number that is too big is also prohibited
    check_bar_threshold("--bar-levels 0.5,0.2,0.3,0.3,0.4,0.1,2",
                        std::vector<double> {0.5, 0.2, 0.3, 0.3, 0.4, 0.1, 2},
                        false);
}
TEST(COMMAND_PARSING, BAR_LEVELS_INVALID)
{
    // the only situation where this will fail in init is if there are
    // non-numeric inputs
    check_bar_threshold("--bar-levels 0.1,0.2,0.3,a,0.4,0.5",
                        std::vector<double> {}, true);
    // if the value overflow,it should also error out
    check_bar_threshold("--bar-levels 0.1,0.2,0.3,0.4,1.79769e+309",
                        std::vector<double> {}, true);
}
TEST(COMMAND_PARSING, FASTSCORE)
{
    mockCommander set;
    std::string command = "--fastscore";
    // default is not using --fastscore
    ASSERT_FALSE(set.get_p_threshold().fastscore);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_p_threshold().fastscore);
}
TEST(COMMAND_PARSING, NO_FULL)
{
    mockCommander set;
    std::string command = "--no-full";
    ASSERT_FALSE(set.get_p_threshold().no_full);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_p_threshold().no_full);
}
TEST(COMMAND_PARSING, NO_CLUMP)
{
    mockCommander set;
    std::string command = "--no-clump";
    ASSERT_FALSE(set.get_clump_info().no_clump);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_clump_info().no_clump);
}
TEST(COMMAND_PARSING, HARD_CODED)
{
    mockCommander set;
    std::string command = "--hard";
    ASSERT_FALSE(set.get_target().hard_coded);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_target().hard_coded);
}
TEST(COMMAND_PARSING, ALLOW_INTER)
{
    mockCommander set;
    std::string command = "--allow-inter";
    ASSERT_FALSE(set.use_inter());
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.use_inter());
}
TEST(COMMAND_PARSING, NON_FOUNDERS)
{
    mockCommander set;
    std::string command = "--nonfounders";
    ASSERT_FALSE(set.nonfounders());
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.nonfounders());
}
TEST(COMMAND_PARSING, BETA)
{
    mockCommander set;
    std::string command = "--beta";
    ASSERT_FALSE(set.get_base().is_beta);
    ASSERT_FALSE(set.get_base().is_or);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_base().is_beta);
    ASSERT_FALSE(set.get_base().is_or);
}
TEST(COMMAND_PARSING, OR)
{
    mockCommander set;
    std::string command = "--or";
    ASSERT_FALSE(set.get_base().is_or);
    ASSERT_FALSE(set.get_base().is_beta);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_base().is_or);
    ASSERT_FALSE(set.get_base().is_beta);
}
TEST(COMMAND_PARSING, INDEX)
{
    mockCommander set;
    std::string command = "--index";
    ASSERT_FALSE(set.get_base().is_index);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_base().is_index);
}
TEST(COMMAND_PARSING, ALLSCORE)
{
    mockCommander set;
    std::string command = "--all-score";
    // default is not using --allscore
    ASSERT_FALSE(set.all_scores());
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.all_scores());
}
TEST(COMMAND_PARSING, IGNORE_FID)
{
    mockCommander set;
    std::string command = "--ignore-fid";
    ASSERT_FALSE(set.get_pheno().ignore_fid);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_pheno().ignore_fid);
}
TEST(COMMAND_PARSING, KEEP_AMBIG)
{
    mockCommander set;
    std::string command = "--keep-ambig";
    ASSERT_FALSE(set.keep_ambig());
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.keep_ambig());
}
TEST(COMMAND_PARSING, KEEP_AMBIG_AS_IS)
{
    mockCommander set;
    std::string command = "--keep-ambig-as-is";
    ASSERT_FALSE(set.ambig_no_flip());
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.ambig_no_flip());
}
TEST(COMMAND_PARSING, NON_CUMULATE)
{
    mockCommander set;
    std::string command = "--non-cumulate";
    ASSERT_FALSE(set.get_prs_instruction().non_cumulate);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_prs_instruction().non_cumulate);
}
TEST(COMMAND_PARSING, NO_REGRESS)
{
    mockCommander set;
    std::string command = "--no-regress";
    ASSERT_FALSE(set.get_prs_instruction().no_regress);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_prs_instruction().no_regress);
}
TEST(COMMAND_PARSING, PRINT_SNP)
{
    mockCommander set;
    std::string command = "--print-snp";
    ASSERT_FALSE(set.print_snp());
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.print_snp());
}
TEST(COMMAND_PARSING, USE_REF_MAF)
{
    mockCommander set;
    std::string command = "--use-ref-maf";
    ASSERT_FALSE(set.get_prs_instruction().use_ref_maf);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_prs_instruction().use_ref_maf);
}
TEST(COMMAND_PARSING, LOGIT_PERM)
{
    mockCommander set;
    std::string command = "--logit-perm";
    ASSERT_FALSE(set.get_perm().logit_perm);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_perm().logit_perm);
}
TEST(COMMAND_PARSING, FULL_BACK)
{
    mockCommander set;
    std::string command = "--full-back";
    ASSERT_FALSE(set.get_set().full_as_background);
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.get_set().full_as_background);
}
TEST(COMMAND_PARSING, NO_DEFAULT)
{
    mockCommander set;
    std::string command = "--no-default";
    ASSERT_FALSE(set.no_default());
    ASSERT_TRUE(set.parse_command_wrapper(command));
    ASSERT_TRUE(set.no_default());
}
std::string get_base_name(const mockCommander& commander, size_t idx)
{
    return commander.get_base().column_name[idx];
}
bool get_has_base(const mockCommander& commander, size_t idx)
{
    return commander.get_base().has_column[idx];
}
void check_set_base_flag(const std::string& command,
                         const std::string& expected,
                         const std::string& default_str, size_t idx)
{
    mockCommander commander;
    ASSERT_STREQ(get_base_name(commander, idx).c_str(), default_str.c_str());
    ASSERT_FALSE(get_has_base(commander, idx));
    if (default_str != expected)
        ASSERT_STRNE(get_base_name(commander, idx).c_str(), expected.c_str());
    ASSERT_TRUE(commander.parse_command_wrapper(command + " " + expected));
    ASSERT_STREQ(get_base_name(commander, idx).c_str(), expected.c_str());
    ASSERT_TRUE(get_has_base(commander, idx));
}

TEST(COMMAND_PARSING, SET_BASE)
{
    check_set_base_flag("--A1", "a", "A1", +BASE_INDEX::EFFECT);
    check_set_base_flag("--a1", "b", "A1", +BASE_INDEX::EFFECT);
    check_set_base_flag("--A2", "c", "A2", +BASE_INDEX::NONEFFECT);
    check_set_base_flag("--a2", "d", "A2", +BASE_INDEX::NONEFFECT);
    check_set_base_flag("--stat", "statistic", "", +BASE_INDEX::STAT);
    check_set_base_flag("--pvalue", "insignificant", "P", +BASE_INDEX::P);
    check_set_base_flag("-p", "postdoc", "P", +BASE_INDEX::P);
    check_set_base_flag("--chr", "chromosome", "CHR", +BASE_INDEX::CHR);
    check_set_base_flag("--bp", "location", "BP", +BASE_INDEX::BP);
    check_set_base_flag("--snp", "cnv", "SNP", +BASE_INDEX::RS);
    mockCommander commander;
    ASSERT_TRUE(commander.get_base().file_name.empty());
    ASSERT_TRUE(commander.parse_command_wrapper("--base BaseInfo"));
    ASSERT_STREQ(commander.get_base().file_name.c_str(), "BaseInfo");
    ASSERT_TRUE(commander.parse_command_wrapper("-b Basic"));
    ASSERT_STREQ(commander.get_base().file_name.c_str(), "Basic");
    check_set_base_flag("--base-info", "INFO_FILTER", "INFO,0.9",
                        +BASE_INDEX::INFO);
    check_set_base_flag("--base-maf", "MAF_FILTER", "", +BASE_INDEX::MAF);
}

void check_binary_target(const std::string& command,
                         const std::vector<bool> expected, bool expect_fail)
{
    mockCommander commander;
    // no default at the beginning
    ASSERT_TRUE(commander.get_pheno().binary.empty());
    bool success =
        commander.parse_command_wrapper("--binary-target " + command);
    if (expect_fail) { ASSERT_FALSE(success); }
    else
    {
        ASSERT_TRUE(success);
        ASSERT_EQ(commander.get_pheno().binary.size(), expected.size());
        for (size_t i = 0; i < expected.size(); ++i)
        { ASSERT_EQ(commander.get_pheno().binary[i], expected[i]); }
    }
}
TEST(COMMAND_PARSING, BINARY_TARGET_INVALID)
{
    // we no longer allow numeric representation of T/F as we need those for
    // parsing
    check_binary_target("1", std::vector<bool> {true}, true);
    check_binary_target("0", std::vector<bool> {false}, true);
    // way too much
    check_binary_target("1e200T", std::vector<bool> {false}, true);
    // Wrong spelling
    check_binary_target("Tru", std::vector<bool> {false}, true);
    // Non numeric multiplier
    check_binary_target("aT", std::vector<bool> {false}, true);
    // negative multiplier
    check_binary_target("-1T", std::vector<bool> {false}, true);
    check_binary_target("F,-10T", std::vector<bool> {false}, true);
    // this in theory is correct, but as the second argument starts with -, and
    // PRSice doesn't have a -1 parameter, it will cause an error
    check_binary_target("F, -10T", std::vector<bool> {false}, true);
    // this is "valid" but wrong in the sense that 3F will not be processed and
    // PRSice can in theory continue to run until we reach check
    check_binary_target("F,2T, 3F", std::vector<bool> {false, true, true},
                        false);
}
TEST(COMMAND_PARSING, BINARY_TARGET_VALID)
{
    // try differnent form of binary target input
    // valid
    check_binary_target("T", std::vector<bool> {true}, false);
    check_binary_target("True", std::vector<bool> {true}, false);
    check_binary_target("true", std::vector<bool> {true}, false);
    check_binary_target("1true", std::vector<bool> {true}, false);
    check_binary_target("1T", std::vector<bool> {true}, false);
    check_binary_target("F", std::vector<bool> {false}, false);
    check_binary_target("False", std::vector<bool> {false}, false);
    check_binary_target("false", std::vector<bool> {false}, false);
    check_binary_target("1false", std::vector<bool> {false}, false);
    check_binary_target("1F", std::vector<bool> {false}, false);
    // more complex
    check_binary_target("4T", std::vector<bool> {true, true, true, true},
                        false);
    check_binary_target(
        "6F", std::vector<bool> {false, false, false, false, false, false},
        false);
    check_binary_target("True,3F",
                        std::vector<bool> {true, false, false, false}, false);
    check_binary_target("True,3F",
                        std::vector<bool> {true, false, false, false}, false);
    // check if it append properly
    mockCommander commander;
    // no default at the beginning
    ASSERT_TRUE(commander.get_pheno().binary.empty());
    std::vector<bool> expected = {true, true, true, false};
    ASSERT_TRUE(commander.parse_command_wrapper("--binary-target 3T,F"));
    ASSERT_EQ(commander.get_pheno().binary.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i)
    { ASSERT_EQ(commander.get_pheno().binary[i], expected[i]); }
    // now second invoke of --binary-target
    ASSERT_TRUE(commander.parse_command_wrapper("--binary-target 2T"));
    expected.push_back(true);
    expected.push_back(true);
    for (size_t i = 0; i < expected.size(); ++i)
    { ASSERT_EQ(commander.get_pheno().binary[i], expected[i]); }
}
TEST(COMMAND_PARSING, DOSAGE)
{
    mockCommander commander;
    ASSERT_DOUBLE_EQ(commander.get_target_qc().dose_threshold, 0.0);
    ASSERT_DOUBLE_EQ(commander.get_target_qc().hard_threshold, 0.1);
    ASSERT_TRUE(commander.parse_command_wrapper("--dose-thres 1.0"));
    ASSERT_DOUBLE_EQ(commander.get_target_qc().dose_threshold, 1.0);
    ASSERT_TRUE(commander.parse_command_wrapper("--hard-thres -0.1"));
    ASSERT_DOUBLE_EQ(commander.get_target_qc().hard_threshold, -0.1);
    // out bound check
    ASSERT_FALSE(commander.parse_command_wrapper("--hard-thres 1e400"));
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().dose_threshold, 0.0);
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().hard_threshold, 0.1);
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-dose-thres 1.0"));
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().dose_threshold, 1.0);
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-hard-thres -0.1"));
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().hard_threshold, -0.1);
    // out bound check
    ASSERT_FALSE(commander.parse_command_wrapper("--ld-hard-thres 1e400"));
}
TEST(COMMAND_PARSING, TARGET_DEFAULT)
{
    mockCommander commander;
    // check default values
    ASSERT_FALSE(commander.get_target().is_ref);
    ASSERT_DOUBLE_EQ(commander.get_target_qc().geno, 1.0);
    ASSERT_DOUBLE_EQ(commander.get_target_qc().info_score, 0.0);
    ASSERT_DOUBLE_EQ(commander.get_target_qc().maf, 0.0);
    ASSERT_TRUE(commander.get_pheno().pheno_file.empty());
    ASSERT_TRUE(commander.get_pheno().pheno_col.empty());
    ASSERT_TRUE(commander.get_pheno().prevalence.empty());
    ASSERT_TRUE(commander.get_target().remove.empty());
    ASSERT_TRUE(commander.get_target().keep.empty());
    ASSERT_TRUE(commander.get_target().file_name.empty());
    ASSERT_TRUE(commander.get_target().file_list.empty());
    ASSERT_STREQ(commander.get_target().type.c_str(), "bed");
}

TEST(COMMAND_PARSING, CLUMP_DEFAULT)
{
    mockCommander commander;
    // check default values
    ASSERT_TRUE(commander.get_reference().is_ref);
    ASSERT_DOUBLE_EQ(commander.get_target_qc().geno, 1.0);
    ASSERT_DOUBLE_EQ(commander.get_target_qc().info_score, 0.0);
    ASSERT_DOUBLE_EQ(commander.get_target_qc().maf, 0.0);
    ASSERT_TRUE(commander.get_reference().remove.empty());
    ASSERT_TRUE(commander.get_reference().keep.empty());
    ASSERT_TRUE(commander.get_reference().file_name.empty());
    ASSERT_TRUE(commander.get_reference().file_list.empty());
    ASSERT_STREQ(commander.get_reference().type.c_str(), "bed");
    ASSERT_DOUBLE_EQ(commander.get_clump_info().r2, 0.1);
    ASSERT_DOUBLE_EQ(commander.get_clump_info().proxy, 0.0);
    ASSERT_EQ(commander.get_clump_info().distance, 250000);
    ASSERT_FALSE(commander.get_clump_info().provided_distance);
    ASSERT_DOUBLE_EQ(commander.get_clump_info().pvalue, 1.0);
}
TEST(COMMAND_PARSING, CLUMP_SETTINGS)
{
    mockCommander commander;
    ASSERT_TRUE(commander.parse_command_wrapper("--clump-p 0.1"));
    ASSERT_DOUBLE_EQ(commander.get_clump_info().pvalue, 0.1);
    ASSERT_TRUE(commander.parse_command_wrapper("--clump-r2 0.5"));
    ASSERT_DOUBLE_EQ(commander.get_clump_info().r2, 0.5);
    ASSERT_TRUE(commander.parse_command_wrapper("--clump-kb 100"));
    ASSERT_DOUBLE_EQ(commander.get_clump_info().distance, 100000);
    ASSERT_TRUE(commander.get_clump_info().provided_distance);
    ASSERT_TRUE(commander.parse_command_wrapper("--clump-kb 100kb"));
    ASSERT_DOUBLE_EQ(commander.get_clump_info().distance, 100000);
    ASSERT_TRUE(commander.parse_command_wrapper("--clump-kb 100b"));
    ASSERT_DOUBLE_EQ(commander.get_clump_info().distance, 100);
    ASSERT_TRUE(commander.parse_command_wrapper("--clump-kb 200mb"));
    ASSERT_DOUBLE_EQ(commander.get_clump_info().distance, 200000000);
}
TEST(COMMAND_PARSING, PRSET)
{
    // wind-3 and --wind-5 use the same function as clump
    mockCommander commander;
    ASSERT_EQ(commander.get_set().wind_3, 0);
    ASSERT_EQ(commander.get_set().wind_5, 0);
    // default is bp
    ASSERT_TRUE(commander.parse_command_wrapper("--wind-5 10"));
    ASSERT_EQ(commander.get_set().wind_5, 10);
    ASSERT_TRUE(commander.parse_command_wrapper("--wind-3 20k"));
    ASSERT_EQ(commander.get_set().wind_3, 20000);
    // now check the background stuff
    ASSERT_TRUE(commander.get_set().background.empty());
    ASSERT_TRUE(commander.get_set().msigdb.empty());
    ASSERT_TRUE(commander.get_set().bed.empty());
    ASSERT_TRUE(commander.get_set().snp.empty());
    ASSERT_TRUE(commander.get_set().feature.empty());
    ASSERT_TRUE(commander.get_set().gtf.empty());
    ASSERT_FALSE(commander.get_set().run);
    ASSERT_TRUE(commander.exclusion_range().empty());
    // now check if they are loaded correctly (doesn't have to be in correct
    // format at the moment)
    ASSERT_TRUE(commander.parse_command_wrapper("--background Name:0"));
    ASSERT_STREQ(commander.get_set().background.c_str(), "Name:0");
    ASSERT_TRUE(commander.parse_command_wrapper("--msigdb kegg"));
    ASSERT_EQ(commander.get_set().msigdb.size(), 1);
    ASSERT_STREQ(commander.get_set().msigdb[0].c_str(), "kegg");
    ASSERT_TRUE(commander.parse_command_wrapper("-m Reactome,MP"));
    // it append
    ASSERT_EQ(commander.get_set().msigdb.size(), 3);
    ASSERT_STREQ(commander.get_set().msigdb[0].c_str(), "kegg");
    ASSERT_STREQ(commander.get_set().msigdb[1].c_str(), "Reactome");
    ASSERT_STREQ(commander.get_set().msigdb[2].c_str(), "MP");
    // GTF
    ASSERT_TRUE(commander.parse_command_wrapper("--gtf Homo"));
    ASSERT_STREQ(commander.get_set().gtf.c_str(), "Homo");
    ASSERT_TRUE(commander.parse_command_wrapper("-g Misc"));
    ASSERT_STREQ(commander.get_set().gtf.c_str(), "Misc");
    // bed B
    ASSERT_TRUE(commander.parse_command_wrapper("--bed File:Name"));
    ASSERT_EQ(commander.get_set().bed.size(), 1);
    ASSERT_STREQ(commander.get_set().bed[0].c_str(), "File:Name");
    ASSERT_TRUE(commander.parse_command_wrapper("-B Something,oK"));
    ASSERT_STREQ(commander.get_set().bed[0].c_str(), "File:Name");
    ASSERT_STREQ(commander.get_set().bed[1].c_str(), "Something");
    ASSERT_STREQ(commander.get_set().bed[2].c_str(), "oK");
    // snp-set
    ASSERT_TRUE(commander.parse_command_wrapper("--snp-set list,of,snp"));
    ASSERT_EQ(commander.get_set().snp.size(), 3);
    ASSERT_STREQ(commander.get_set().snp[0].c_str(), "list");
    ASSERT_STREQ(commander.get_set().snp[1].c_str(), "of");
    ASSERT_STREQ(commander.get_set().snp[2].c_str(), "snp");
    // feature
    ASSERT_TRUE(commander.parse_command_wrapper("--feature gene"));
    ASSERT_EQ(commander.get_set().feature.size(), 1);
    ASSERT_STREQ(commander.get_set().feature[0].c_str(), "gene");
    // no duplicate check
    ASSERT_TRUE(commander.parse_command_wrapper("--feature protein,gene"));
    ASSERT_EQ(commander.get_set().feature.size(), 3);
    ASSERT_STREQ(commander.get_set().feature[0].c_str(), "gene");
    ASSERT_STREQ(commander.get_set().feature[1].c_str(), "protein");
    ASSERT_STREQ(commander.get_set().feature[2].c_str(), "gene");
    // Exclusion range is a direct loading
    ASSERT_TRUE(commander.parse_command_wrapper("--x-range chr6:1-10"));
    ASSERT_STREQ(commander.exclusion_range().c_str(), "chr6:1-10");
    // we don't even tokenize it
    ASSERT_TRUE(
        commander.parse_command_wrapper("--x-range chr6:1-10,chr22:133:288"));
    ASSERT_STREQ(commander.exclusion_range().c_str(),
                 "chr6:1-10,chr22:133:288");
}

TEST(COMMAND_PARSING, MISC)
{
    mockCommander commander;
    // check defaults
    ASSERT_STREQ(commander.out().c_str(), "PRSice");
    ASSERT_EQ(commander.get_prs_instruction().thread, 1);
    ASSERT_EQ(commander.memory(), 1e10);
    ASSERT_STREQ(commander.delim().c_str(), " ");
    ASSERT_TRUE(commander.exclude_file().empty());
    ASSERT_TRUE(commander.extract_file().empty());
    ASSERT_TRUE(commander.parse_command_wrapper("--out PRSet"));
    ASSERT_STREQ(commander.out().c_str(), "PRSet");
    int32_t max_thread = commander.max_thread();
    if (max_thread > 2)
    {
        ASSERT_TRUE(commander.parse_command_wrapper("--thread 2"));
        ASSERT_EQ(commander.get_prs_instruction().thread, 2);
    }
    ASSERT_TRUE(commander.parse_command_wrapper("--thread "
                                                + std::to_string(max_thread)));
    ASSERT_EQ(commander.get_prs_instruction().thread, max_thread);
    // reset it first
    ASSERT_TRUE(commander.parse_command_wrapper("--thread 1"));
    ASSERT_EQ(commander.get_prs_instruction().thread, 1);
    ASSERT_TRUE(commander.parse_command_wrapper("--thread max"));
    ASSERT_EQ(commander.get_prs_instruction().thread, max_thread);
    ASSERT_TRUE(commander.parse_command_wrapper("--extract Love"));
    ASSERT_STREQ(commander.extract_file().c_str(), "Love");
    ASSERT_TRUE(commander.parse_command_wrapper("--exclude Hate"));
    ASSERT_STREQ(commander.exclude_file().c_str(), "Hate");
    /*
     not sure how to do proper escape here, will only do very simple cases
    ASSERT_TRUE(commander.parse_command_wrapper("--id-delim \"-\""));
    ASSERT_STREQ(commander.delim().c_str(), "-");
    ASSERT_TRUE(commander.parse_command_wrapper("--id-delim \"  \""));
    ASSERT_STREQ(commander.delim().c_str(), "  ");
    */
    ASSERT_TRUE(commander.parse_command_wrapper("--id-delim -"));
    ASSERT_STREQ(commander.delim().c_str(), "-");
    ASSERT_TRUE(commander.parse_command_wrapper("--memory 1k"));
    ASSERT_EQ(commander.memory(), 1024);
    // default is mb
    ASSERT_TRUE(commander.parse_command_wrapper("--memory 10"));
    ASSERT_EQ(commander.memory(), 10485760);
    ASSERT_TRUE(commander.parse_command_wrapper("--memory 1gb"));
    ASSERT_EQ(commander.memory(), 1073741824);
    ASSERT_TRUE(commander.parse_command_wrapper("--memory 30tb"));
    ASSERT_EQ(commander.memory(), 32985348833280);
    // the default of seed is random, which is difficult to test. So we will
    // just check if we set the seed correctly
    ASSERT_TRUE(commander.parse_command_wrapper("--seed 123"));
    ASSERT_EQ(commander.get_perm().seed, 123);
    // check permutation default
    ASSERT_EQ(commander.get_perm().num_permutation, 0);
    ASSERT_FALSE(commander.get_perm().run_perm);
    ASSERT_FALSE(commander.get_perm().run_set_perm);
    ASSERT_TRUE(commander.parse_command_wrapper("--perm 100"));
    ASSERT_EQ(commander.get_perm().num_permutation, 100);
    ASSERT_TRUE(commander.get_perm().run_perm);
    ASSERT_FALSE(commander.get_perm().run_set_perm);
    ASSERT_TRUE(commander.parse_command_wrapper("--set-perm 1026"));
    ASSERT_EQ(commander.get_perm().num_permutation, 1026);
    ASSERT_TRUE(commander.get_perm().run_set_perm);
    // we won't change the other
    ASSERT_TRUE(commander.get_perm().run_perm);
    // now check for overflow
    ASSERT_FALSE(commander.parse_command_wrapper("--set-perm 1e200"));
}

TEST(COMMAND_PARSE, PRS_MODEL_THRESHOLD)
{
    mockCommander commander;
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().inter, 0.00005);
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().lower, 5e-8);
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().upper, 0.5);
    ASSERT_TRUE(commander.parse_command_wrapper("--inter 1e-10"));
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().inter, 1e-10);
    ASSERT_TRUE(commander.parse_command_wrapper("-i 110"));
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().inter, 110);
    ASSERT_TRUE(commander.parse_command_wrapper("--lower 1e-20"));
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().lower, 1e-20);
    ASSERT_TRUE(commander.parse_command_wrapper("-l 123"));
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().lower, 123);
    ASSERT_TRUE(commander.parse_command_wrapper("--upper 5e-70"));
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().upper, 5e-70);
    ASSERT_TRUE(commander.parse_command_wrapper("-u 5e10"));
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().upper, 5e10);
    //  extreme value
    ASSERT_TRUE(commander.parse_command_wrapper("-u 5e300"));
    ASSERT_DOUBLE_EQ(commander.get_p_threshold().upper, 5e300);
    // fail
    ASSERT_FALSE(commander.parse_command_wrapper("-u 5e400"));
    ASSERT_FALSE(commander.parse_command_wrapper("-l -5e400"));
    ASSERT_FALSE(commander.parse_command_wrapper("-i hi"));
    // MODEL and SCORES
    ASSERT_EQ(commander.get_prs_instruction().scoring_method, SCORING::AVERAGE);
    ASSERT_TRUE(commander.parse_command_wrapper("--score Sum"));
    ASSERT_EQ(commander.get_prs_instruction().scoring_method, SCORING::SUM);
    ASSERT_TRUE(commander.parse_command_wrapper("--score std"));
    ASSERT_EQ(commander.get_prs_instruction().scoring_method,
              SCORING::STANDARDIZE);
    ASSERT_TRUE(commander.parse_command_wrapper("--score con-std"));
    ASSERT_EQ(commander.get_prs_instruction().scoring_method,
              SCORING::CONTROL_STD);
    ASSERT_TRUE(commander.parse_command_wrapper("--score avg"));
    ASSERT_EQ(commander.get_prs_instruction().scoring_method, SCORING::AVERAGE);
    // we do exact match
    ASSERT_FALSE(commander.parse_command_wrapper("--score averaging"));

    ASSERT_EQ(commander.get_prs_instruction().missing_score,
              MISSING_SCORE::MEAN_IMPUTE);
    ASSERT_TRUE(commander.parse_command_wrapper("--missing SET_Zero"));
    ASSERT_EQ(commander.get_prs_instruction().missing_score,
              MISSING_SCORE::SET_ZERO);
    ASSERT_TRUE(commander.parse_command_wrapper("--missing Center"));
    ASSERT_EQ(commander.get_prs_instruction().missing_score,
              MISSING_SCORE::CENTER);
    ASSERT_TRUE(commander.parse_command_wrapper("--missing mean_impute"));
    ASSERT_EQ(commander.get_prs_instruction().missing_score,
              MISSING_SCORE::MEAN_IMPUTE);
    // Allowed, but don't think I have implemented this yet
    ASSERT_TRUE(commander.parse_command_wrapper("--missing IMPUTE_CONTROL"));
    ASSERT_EQ(commander.get_prs_instruction().missing_score,
              MISSING_SCORE::IMPUTE_CONTROL);
    // We only matched the first character
    ASSERT_TRUE(commander.parse_command_wrapper("--missing cat"));
    ASSERT_EQ(commander.get_prs_instruction().missing_score,
              MISSING_SCORE::CENTER);
    // but should fail if we have something that starts with different character
    ASSERT_FALSE(commander.parse_command_wrapper("--missing beatrice"));
    ASSERT_TRUE(commander.parse_command_wrapper("--model dom"));
    ASSERT_EQ(commander.get_prs_instruction().genetic_model, MODEL::DOMINANT);
    ASSERT_TRUE(commander.parse_command_wrapper("--model het"));
    ASSERT_EQ(commander.get_prs_instruction().genetic_model,
              MODEL::HETEROZYGOUS);
    ASSERT_TRUE(commander.parse_command_wrapper("--model rec"));
    ASSERT_EQ(commander.get_prs_instruction().genetic_model, MODEL::RECESSIVE);
    ASSERT_TRUE(commander.parse_command_wrapper("--model ADD"));
    ASSERT_EQ(commander.get_prs_instruction().genetic_model, MODEL::ADDITIVE);
    // similar to missing
    ASSERT_TRUE(commander.parse_command_wrapper("--model darwin"));
    ASSERT_EQ(commander.get_prs_instruction().genetic_model, MODEL::DOMINANT);
    ASSERT_FALSE(commander.parse_command_wrapper("--model mendel"));
}
void check_cov_loading(const std::string& command,
                       const std::vector<std::string>& expected,
                       const bool expect_fail, const bool factor = false)
{
    mockCommander commander;
    bool success = commander.parse_command_wrapper(command);
    if (expect_fail) { ASSERT_FALSE(success); }
    else if (!factor)
    {
        ASSERT_TRUE(success);
        ASSERT_EQ(expected.size(), commander.get_pheno().cov_colname.size());
        for (size_t i = 0; i < expected.size(); ++i)
        {
            ASSERT_STREQ(expected[i].c_str(),
                         commander.get_pheno().cov_colname[i].c_str());
        }
    }
    else
    {
        ASSERT_TRUE(success);
        ASSERT_EQ(expected.size(), commander.get_pheno().factor_cov.size());
        for (size_t i = 0; i < expected.size(); ++i)
        {
            ASSERT_STREQ(expected[i].c_str(),
                         commander.get_pheno().factor_cov[i].c_str());
        }
    }
}
TEST(COMMAND_PARSING, COVARIATE)
{
    mockCommander commander;
    ASSERT_TRUE(commander.get_pheno().cov_file.empty());
    ASSERT_TRUE(commander.get_pheno().cov_colname.empty());
    ASSERT_TRUE(commander.get_pheno().factor_cov.empty());
    ASSERT_TRUE(commander.parse_command_wrapper("--cov Covar"));
    ASSERT_STREQ(commander.get_pheno().cov_file.c_str(), "Covar");
    check_cov_loading("--cov-col Testing,@PC[1-55]",
                      std::vector<std::string> {"Testing", "@PC[1-55]"}, false);
    // this should be allowed
    check_cov_loading("--cov-col Testing,@PC[1.3.5]",
                      std::vector<std::string> {"Testing", "@PC[1.3.5]"},
                      false);
    // this will be stored but shouldn't pass the check

    check_cov_loading("--cov-col Testing,@PC[1,3,5]",
                      std::vector<std::string> {"Testing", "@PC[1", "3", "5]"},
                      false);
    // cov-factor uses the same function as cov-col, so will only test if it is
    // set properly
    check_cov_loading("--cov-factor Sex", std::vector<std::string> {"Sex"},
                      false, true);
    // can append
    ASSERT_TRUE(commander.parse_command_wrapper("--cov-col Testing,@PC[1-55]"));
    std::vector<std::string> expected {"Testing", "@PC[1-55]"};
    ASSERT_EQ(expected.size(), commander.get_pheno().cov_colname.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        ASSERT_STREQ(expected[i].c_str(),
                     commander.get_pheno().cov_colname[i].c_str());
    }
    ASSERT_TRUE(commander.parse_command_wrapper("--cov-col More,Covariate"));
    expected.push_back("More");
    expected.push_back("Covariate");
    ASSERT_EQ(expected.size(), commander.get_pheno().cov_colname.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        ASSERT_STREQ(expected[i].c_str(),
                     commander.get_pheno().cov_colname[i].c_str());
    }
}
TEST(COMMAND_PARSING, REFERENCE_FILE)
{
    mockCommander commander;
    ASSERT_TRUE(commander.get_reference().is_ref);
    ASSERT_TRUE(commander.parse_command_wrapper("--ld genotype"));
    ASSERT_STREQ(commander.get_reference().file_name.c_str(), "genotype");
    ASSERT_TRUE(commander.parse_command_wrapper("-L plink"));
    ASSERT_STREQ(commander.get_reference().file_name.c_str(), "plink");
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-list testing"));
    ASSERT_STREQ(commander.get_reference().file_list.c_str(), "testing");
    // default is bed, currently don't do any check
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-type bgen"));
    ASSERT_STREQ(commander.get_reference().type.c_str(), "bgen");
    // so in theory, we can set whatever string we like
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-type beatrice"));
    ASSERT_STREQ(commander.get_reference().type.c_str(), "beatrice");
    // check keep and remove is correct
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-keep fun"));
    ASSERT_STREQ(commander.get_reference().keep.c_str(), "fun");
    ASSERT_TRUE(commander.get_reference().remove.empty());
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-remove depression"));
    ASSERT_STREQ(commander.get_reference().remove.c_str(), "depression");
    // we should not change keep when we set keep
    ASSERT_STREQ(commander.get_reference().keep.c_str(), "fun");
    ASSERT_TRUE(commander.get_reference().is_ref);
}
TEST(COMMAND_PARSING, TARGET_FILE)
{
    mockCommander commander;
    ASSERT_FALSE(commander.get_target().is_ref);
    ASSERT_TRUE(commander.parse_command_wrapper("--target genotype"));
    ASSERT_STREQ(commander.get_target().file_name.c_str(), "genotype");
    ASSERT_TRUE(commander.parse_command_wrapper("-t plink"));
    ASSERT_STREQ(commander.get_target().file_name.c_str(), "plink");
    ASSERT_TRUE(commander.parse_command_wrapper("--target-list testing"));
    ASSERT_STREQ(commander.get_target().file_list.c_str(), "testing");
    // default is bed, currently don't do any check
    ASSERT_TRUE(commander.parse_command_wrapper("--type bgen"));
    ASSERT_STREQ(commander.get_target().type.c_str(), "bgen");
    // so in theory, we can set whatever string we like
    ASSERT_TRUE(commander.parse_command_wrapper("--type beatrice"));
    ASSERT_STREQ(commander.get_target().type.c_str(), "beatrice");
    // check keep and remove is correct
    ASSERT_TRUE(commander.parse_command_wrapper("--keep fun"));
    ASSERT_STREQ(commander.get_target().keep.c_str(), "fun");
    ASSERT_TRUE(commander.get_target().remove.empty());
    ASSERT_TRUE(commander.parse_command_wrapper("--remove depression"));
    ASSERT_STREQ(commander.get_target().remove.c_str(), "depression");
    // we should not change keep when we set keep
    ASSERT_STREQ(commander.get_target().keep.c_str(), "fun");
    ASSERT_FALSE(commander.get_target().is_ref);
}
TEST(COMMAND_PARSING, PHENO_SET)
{
    mockCommander commander;
    ASSERT_TRUE(commander.get_pheno().pheno_file.empty());
    ASSERT_TRUE(commander.get_pheno().pheno_col.empty());
    ASSERT_TRUE(commander.get_pheno().pheno_col_idx.empty());
    ASSERT_TRUE(commander.parse_command_wrapper("--pheno Phenotype"));
    ASSERT_STREQ(commander.get_pheno().pheno_file.c_str(), "Phenotype");
    ASSERT_TRUE(commander.get_pheno().pheno_col.empty());
    ASSERT_TRUE(commander.get_pheno().pheno_col_idx.empty());
    ASSERT_TRUE(commander.parse_command_wrapper("--pheno-col A1,B1"));
    std::vector<std::string> expected = {"A1", "B1"};
    ASSERT_EQ(commander.get_pheno().pheno_col.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        ASSERT_STREQ(commander.get_pheno().pheno_col[i].c_str(),
                     expected[i].c_str());
    }
    ASSERT_TRUE(commander.get_pheno().pheno_col_idx.empty());
    // We do allow multiple use of --pheno-col, though not sure if that is a
    // good idea or not
    ASSERT_TRUE(commander.parse_command_wrapper("--pheno-col C2,D2"));
    expected = {"A1", "B1", "C2", "D2"};
    ASSERT_EQ(commander.get_pheno().pheno_col.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        ASSERT_STREQ(commander.get_pheno().pheno_col[i].c_str(),
                     expected[i].c_str());
    }
    // now check prevalence
    ASSERT_TRUE(commander.get_pheno().prevalence.empty());
    ASSERT_TRUE(commander.parse_command_wrapper("-k 0.1,0.3,1,3"));
    // there is no bound check yet
    std::vector<double> expected_prev = {0.1, 0.3, 1, 3};
    ASSERT_EQ(commander.get_pheno().prevalence.size(), expected_prev.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        ASSERT_DOUBLE_EQ(commander.get_pheno().prevalence[i], expected_prev[i]);
    }
    // also check long flag
    // use new mockCommander, as prevalence should stack
    mockCommander second_command;
    ASSERT_TRUE(
        second_command.parse_command_wrapper("--prevalence -0.1,0.44,1e-5"));
    // there is no bound check yet
    expected_prev.clear();
    expected_prev = {-0.1, 0.44, 1e-5};
    ASSERT_EQ(second_command.get_pheno().prevalence.size(),
              expected_prev.size());
    for (size_t i = 0; i < expected_prev.size(); ++i)
    {
        ASSERT_DOUBLE_EQ(second_command.get_pheno().prevalence[i],
                         expected_prev[i]);
    }
    ASSERT_TRUE(
        second_command.parse_command_wrapper("--prevalence 0.1,0.3,0.5"));
    // check stacking
    expected_prev.push_back(0.1);
    expected_prev.push_back(0.3);
    expected_prev.push_back(0.5);
    ASSERT_EQ(second_command.get_pheno().prevalence.size(),
              expected_prev.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        ASSERT_DOUBLE_EQ(second_command.get_pheno().prevalence[i],
                         expected_prev[i]);
    }
    // check out of bound
    ASSERT_FALSE(second_command.parse_command_wrapper("-k 1e-400"));
    // non-numeric
    ASSERT_FALSE(second_command.parse_command_wrapper("-k common_disease"));
}
TEST(COMMAND_PARSING, TARGET_FILTER_CHECK)
{
    // now check the get set combo works
    mockCommander commander;
    // first check valid inputs
    ASSERT_TRUE(commander.parse_command_wrapper("--geno 0.4"));
    ASSERT_DOUBLE_EQ(commander.get_target_qc().geno, 0.4);
    ASSERT_TRUE(commander.parse_command_wrapper("--info 0.2"));
    ASSERT_DOUBLE_EQ(commander.get_target_qc().info_score, 0.2);
    ASSERT_TRUE(commander.parse_command_wrapper("--maf 0.01"));
    ASSERT_DOUBLE_EQ(commander.get_target_qc().maf, 0.01);
    // out of bound input (check later, so should still be valid as of now
    ASSERT_TRUE(commander.parse_command_wrapper("--geno -0.4"));
    ASSERT_DOUBLE_EQ(commander.get_target_qc().geno, -0.4);
    ASSERT_TRUE(commander.parse_command_wrapper("--info 20"));
    ASSERT_DOUBLE_EQ(commander.get_target_qc().info_score, 20);
    ASSERT_TRUE(commander.parse_command_wrapper("--maf -10.01"));
    ASSERT_DOUBLE_EQ(commander.get_target_qc().maf, -10.01);
    // the invalid input e.g non-numeric
    ASSERT_FALSE(commander.parse_command_wrapper("--geno --0.4"));
    ASSERT_FALSE(commander.parse_command_wrapper("--geno geno"));
    ASSERT_FALSE(commander.parse_command_wrapper("--info --0.2"));
    ASSERT_FALSE(commander.parse_command_wrapper("--info test_yourself"));
    ASSERT_FALSE(commander.parse_command_wrapper("--maf -+0.01"));
    ASSERT_FALSE(commander.parse_command_wrapper("--maf rare"));
}

TEST(COMMAND_PARSING, REFERENCE_FILTER_CHECK)
{
    // now check the get set combo works
    mockCommander commander;
    // first check valid inputs
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-geno 0.4"));
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().geno, 0.4);
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-info 0.2"));
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().info_score, 0.2);
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-maf 0.01"));
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().maf, 0.01);
    // out of bound input (check later, so should still be valid as of now
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-geno -0.4"));
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().geno, -0.4);
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-info 20"));
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().info_score, 20);
    ASSERT_TRUE(commander.parse_command_wrapper("--ld-maf -10.01"));
    ASSERT_DOUBLE_EQ(commander.get_ref_qc().maf, -10.01);
    // the invalid input e.g non-numeric
    ASSERT_FALSE(commander.parse_command_wrapper("--ld-geno --0.4"));
    ASSERT_FALSE(commander.parse_command_wrapper("--ld-geno geno"));
    ASSERT_FALSE(commander.parse_command_wrapper("--ld-info --0.2"));
    ASSERT_FALSE(commander.parse_command_wrapper("--ld-info test_yourself"));
    ASSERT_FALSE(commander.parse_command_wrapper("--ld-maf -+0.01"));
    ASSERT_FALSE(commander.parse_command_wrapper("--ld-maf rare"));
}
void invalid_cov_input(const std::string& cov_string)
{
    try
    {
        // invalid input
        std::vector<std::string> results =
            mockCommander::transform_covariate(cov_string);
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
}
TEST(COVARIATE_TRANSFORM, RANGE_CHECK)
{

    std::string cov = "PC[1-5]";
    size_t res;
    ASSERT_FALSE(mockCommander::find_first_end_wrapper(cov, 0, res));
    ASSERT_TRUE(mockCommander::find_first_end_wrapper(cov, 2, res));
    ASSERT_EQ(res, 6);
    cov = "PC[1-5[1-5]]";
    ASSERT_FALSE(mockCommander::find_first_end_wrapper(cov, 0, res));
    ASSERT_FALSE(mockCommander::find_first_end_wrapper(cov, 2, res));
}
TEST(COVARIATE_TRANSFORM, PARSE_RANGE)
{
    std::string cov = "[1-5]";
    std::vector<size_t> res;
    // we expect [] to be removed
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
    cov = "1-5";
    res.clear();
    ASSERT_TRUE(mockCommander::parse_range_wrapper(cov, res));
    ASSERT_EQ(res.size(), 5);
    for (size_t i = 0; i < res.size(); ++i) ASSERT_EQ(res[i], i + 1);
    cov = "10-50";
    res.clear();
    ASSERT_TRUE(mockCommander::parse_range_wrapper(cov, res));
    ASSERT_EQ(res.size(), 41);
    for (size_t i = 0; i < res.size(); ++i) ASSERT_EQ(res[i], i + 10);
    cov = "50-10";
    res.clear();
    ASSERT_TRUE(mockCommander::parse_range_wrapper(cov, res));
    ASSERT_EQ(res.size(), 41);
    for (size_t i = 0; i < res.size(); ++i) ASSERT_EQ(res[i], i + 10);
    cov = "10";
    res.clear();
    ASSERT_TRUE(mockCommander::parse_range_wrapper(cov, res));
    ASSERT_EQ(res.size(), 1);
    ASSERT_EQ(res[0], 10);
    // should not work for negative
    cov = "-1";
    res.clear();
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
    cov = "1--5";
    res.clear();
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
    cov = "-1--5";
    res.clear();
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
    res.clear();
    cov = "1,5";
    // we assume , is already dealt with
    ASSERT_FALSE(mockCommander::parse_range_wrapper(cov, res));
}

TEST(COVARIATE_TRANSFORM, GET_RANGE)
{
    std::string cov = "PC[1-5]";
    std::vector<size_t> res;
    // format cannot be converted
    ASSERT_FALSE(mockCommander::get_range_wrapper(cov, 0, 6, res));
    // still wrong format
    ASSERT_FALSE(mockCommander::get_range_wrapper(cov, 2, 4, res));
    // out of bound
    ASSERT_FALSE(mockCommander::get_range_wrapper(cov, 2, 7, res));
    res.clear();
    ASSERT_TRUE(mockCommander::get_range_wrapper(cov, 2, 6, res));
    ASSERT_EQ(res.size(), 5);
    for (size_t i = 0; i < res.size(); ++i) { ASSERT_EQ(res[i], i + 1); }
    // complex options
    cov = "PC[1-5.8.7-10]";
    res.clear();
    ASSERT_TRUE(mockCommander::get_range_wrapper(cov, 2, 13, res));
    // should be sorted and removed the duplicates (8)
    std::vector<size_t> expected = {1, 2, 3, 4, 5, 7, 8, 9, 10};
    ASSERT_EQ(res.size(), expected.size());
    for (size_t i = 0; i < res.size(); ++i) { ASSERT_EQ(res[i], expected[i]); }
    cov = "PC[1-5.-6]";
    // One fail, all fail
    res.clear();
    ASSERT_FALSE(mockCommander::get_range_wrapper(cov, 2, 9, res));
    // as we use . to separate each input, it means double value will be parsed
    // to something else. User will have to read the log to check if the parsing
    // is correct
    cov = "PC[1-5.6.0.005]";
    ASSERT_TRUE(mockCommander::get_range_wrapper(cov, 2, 14, res));
    expected.clear();
    expected = {0, 1, 2, 3, 4, 5, 6};
    ASSERT_EQ(res.size(), expected.size());
    for (size_t i = 0; i < res.size(); ++i) { ASSERT_EQ(res[i], expected[i]); }
}
TEST(COVARIATE_TRANSFORM, UPDATE_COVARIATE_WITH_RANGE)
{
    std::vector<std::string> result;
    std::vector<size_t> range;
    ASSERT_FALSE(mockCommander::update_covariate_ranges_wrapper(result, range));
    range = {1, 3, 5, 7, 9};
    result.clear();
    std::vector<std::string> expected = {"1", "3", "5", "7", "9"};
    ASSERT_TRUE(mockCommander::update_covariate_ranges_wrapper(result, range));
    ASSERT_EQ(result.size(), range.size());
    for (size_t i = 0; i < result.size(); ++i)
    { ASSERT_STREQ(result[i].c_str(), expected[i].c_str()); }
    result.clear();
    result = {"PC1AB", "PC2AB"};
    // sequence = iterate result then range
    expected = {"PC1AB1", "PC1AB3", "PC1AB5", "PC1AB7", "PC1AB9",
                "PC2AB1", "PC2AB3", "PC2AB5", "PC2AB7", "PC2AB9"};
    ASSERT_TRUE(mockCommander::update_covariate_ranges_wrapper(result, range));
    ASSERT_EQ(result.size(), expected.size());
    for (size_t i = 0; i < result.size(); ++i)
    { ASSERT_STREQ(result[i].c_str(), expected[i].c_str()); }
}
void transform_test(const std::string& cov,
                    const std::vector<std::string>& expected,
                    bool expect_success)
{
    std::vector<std::string> result;
    if (!expect_success)
    { ASSERT_FALSE(mockCommander::transform_wrapper(cov, result)); }
    else
    {
        ASSERT_TRUE(mockCommander::transform_wrapper(cov, result));
        ASSERT_EQ(result.size(), expected.size());
        for (size_t i = 0; i < result.size(); ++i)
        { ASSERT_STREQ(result[i].c_str(), expected[i].c_str()); }
    }
}
TEST(COVARIATE_TRANSFORM, TRANSFORMATION)
{
    transform_test("@PC", std::vector<std::string> {"PC"}, true);
    transform_test("@@PC", std::vector<std::string> {"@PC"}, true);
    transform_test("PC1", std::vector<std::string> {"PC1"}, true);
    transform_test("PC1-5", std::vector<std::string> {"PC1-5"}, true);
    transform_test("PC1@5", std::vector<std::string> {"PC1@5"}, true);
    transform_test("@PC1-5", std::vector<std::string> {"PC1-5"}, true);
    transform_test("@PC[1-5]",
                   std::vector<std::string> {"PC1", "PC2", "PC3", "PC4", "PC5"},
                   true);
    transform_test("@PC[1-2.5]", std::vector<std::string> {"PC1", "PC2", "PC5"},
                   true);
    transform_test(
        "@PC[1-2.4.3-6]",
        std::vector<std::string> {"PC1", "PC2", "PC3", "PC4", "PC5", "PC6"},
        true);
    transform_test("@PC[1-2]A", std::vector<std::string> {"PC1A", "PC2A"},
                   true);
    transform_test(
        "@PC[1-2]A[1-2]",
        std::vector<std::string> {"PC1A1", "PC1A2", "PC2A1", "PC2A2"}, true);
}

void quick_check_unit(const std::string& input_str, const size_t exp_output,
                      const size_t def_power = 0, const bool memory = false)
{
    mockCommander commander;
    size_t value;
    commander.check_parse_unit_value(input_str, "", def_power, value, memory);
    ASSERT_EQ(exp_output, value);
}

TEST(PARSE_UNIT, VALIDITY)
{
    mockCommander commander;
    size_t value = 0;
    // Check if valid
    ASSERT_FALSE(commander.check_parse_unit_value("m", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("b", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("mb", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("hi", "--mem", 0, value));
    ASSERT_FALSE(commander.check_parse_unit_value("TB", "--mem", 0, value));
}
TEST(PARSE_UNIT, OUT_BOUND)
{
    mockCommander commander;
    size_t value = 0;
    // out of bound
    ASSERT_FALSE(commander.check_parse_unit_value("1", "", 7, value));
    ASSERT_FALSE(
        commander.check_parse_unit_value("1000000000tb", "", 1, value));
}
TEST(PARSE_UNIT, NEGATIVES)
{
    // Check for negative values
    mockCommander commander;
    size_t value = 0;
    ASSERT_FALSE(commander.check_parse_unit_value("-1", "", 1, value));
    ASSERT_FALSE(commander.check_parse_unit_value("-1tb", "", 1, value));
}
TEST(PARSE_UNIT, WITH_UNIT)
{
    // default value should be ignored when user provide a unit
    quick_check_unit("1b", 1, 0);
    quick_check_unit("1b", 1, 1);
    quick_check_unit("1b", 1, 3);
    quick_check_unit("1b", 1, 4);
    quick_check_unit("1b", 1, 5);
    quick_check_unit("1b", 1, 6);
}
TEST(PARSE_UNIT, DIFFERENT_UNIT)
{
    // check unit works as expected
    quick_check_unit("1k", 1000, 0);
    quick_check_unit("1kb", 1000, 0);
    quick_check_unit("1m", 1000000, 0);
    quick_check_unit("1mb", 1000000, 0);
    quick_check_unit("1g", 1000000000, 0);
    quick_check_unit("1gb", 1000000000, 0);
    quick_check_unit("1t", 1000000000000, 0);
    quick_check_unit("1tb", 1000000000000, 0);
}
TEST(PARSE_UNIT, NON_INTEGER)
{
    mockCommander commander;
    size_t value = 0;
    // Check non-integer scenarios
    quick_check_unit("1.5k", 1500, 0);
    quick_check_unit("1.004k", 1004, 0);
    ASSERT_FALSE(commander.check_parse_unit_value("1.5b", "", 1, value));
}
TEST(PARSE_UNIT, DEFAULT_VALUE)
{
    // check default value works
    quick_check_unit("1", 1, 0);
    quick_check_unit("1", 1000, 1);
    quick_check_unit("1", 1000000000, 3);
    quick_check_unit("1", 1000000000000, 4);
    quick_check_unit("1", 1000000000000000, 5);
    quick_check_unit("1", 1000000000000000000, 6);
}
#endif // COMMANDER_TEST_H
