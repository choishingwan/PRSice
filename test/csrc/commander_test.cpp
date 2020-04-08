#include "catch.hpp"
#include "commander.hpp"

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
    bool target_check_wrapper() { return target_check(); }
    bool prsice_check_wrapper() { return prsice_check(); }
    bool clump_check_wrapper() { return clump_check(); }
    bool ref_check_wrapper() { return ref_check(); }
    bool misc_check_wrapper() { return misc_check(); }
    bool filter_check_wrapper() { return filter_check(); }
    bool prset_check_wrapper() { return prset_check(); }
    bool base_check_wrapper()
    {
        try
        {
            return base_check();
        }
        catch (const std::runtime_error&)
        {
            return false;
        }
    }
    bool base_column_check_wrapper(std::vector<std::string>& column_names)
    {
        return base_column_check(column_names);
    }
    bool pheno_check_wrapper(bool is_beta)
    {
        m_ran_base_check = true;
        m_base_info.is_beta = is_beta;
        m_base_info.is_or = !is_beta;
        return pheno_check();
    }
    std::string get_error() const { return m_error_message; }
    int32_t max_thread() { return maximum_thread(); }
    auto get_cov_names_wrap() { return get_cov_names(); }
    size_t
    find_cov_idx_wrap(const std::unordered_set<std::string>& included,
                      const std::unordered_map<std::string, size_t>& ref_index,
                      std::string& missing)
    {
        return find_cov_idx(included, ref_index, missing);
    }
    void reorganize_cov_name_wrap(const std::vector<std::string>& cov_header)
    {
        reorganize_cov_name(cov_header);
    }
    bool process_factor_cov_wrap(
        const std::unordered_set<std::string>& included,
        const std::unordered_map<std::string, size_t>& ref_index,
        const std::unordered_set<std::string>& ori_input)
    {
        return process_factor_cov(included, ref_index, ori_input);
    }
};


std::string get_base_name(const mockCommander& commander, size_t idx)
{
    return commander.get_base().column_name[idx];
}
bool get_has_base(const mockCommander& commander, size_t idx)
{
    return commander.get_base().has_column[idx];
}
void check_set_base_flag(mockCommander& commander, const std::string& command,
                         const std::string& expected,
                         const std::string& default_str, size_t idx)
{
    REQUIRE(get_base_name(commander, idx) == default_str);
    REQUIRE_FALSE(get_has_base(commander, idx));
    if (default_str != expected)
    { REQUIRE_FALSE(get_base_name(commander, idx) == expected); }
    REQUIRE(commander.parse_command_wrapper(command + " " + expected));
    REQUIRE(get_base_name(commander, idx) == expected);
    REQUIRE(get_has_base(commander, idx));
}

TEST_CASE("Commander")
{
    Commander base_commander;
    SECTION("Initialize")
    {
        SECTION("Test Default")
        {
            REQUIRE_FALSE(base_commander.all_scores());
            REQUIRE_FALSE(base_commander.use_inter());
            REQUIRE(base_commander.delim() == " ");
            REQUIRE_THAT(base_commander.out(),
                         Catch::Matchers::Equals("PRSice"));
            REQUIRE(base_commander.exclusion_range().empty());
            REQUIRE(base_commander.exclude_file().empty());
            REQUIRE(base_commander.extract_file().empty());
            REQUIRE(base_commander.max_memory(1.0) == Approx(1.0));
            REQUIRE(base_commander.max_memory(2.0) == Approx(2.0));
        }
    }

    mockCommander commander;
    SECTION("No argument should throw error")
    {
        bool early_terminate;
        REQUIRE_FALSE(commander.parse_command_wrapper("", early_terminate));
        REQUIRE_FALSE(early_terminate);
    }
    SECTION("Print Usage with --help")
    {
        bool early_terminate;
        REQUIRE_FALSE(
            commander.parse_command_wrapper("--help", early_terminate));
        REQUIRE(early_terminate);
    }
    SECTION("Print Usage with -h")
    {
        bool early_terminate;
        REQUIRE_FALSE(commander.parse_command_wrapper("-h", early_terminate));
        REQUIRE(early_terminate);
    }
    SECTION("Invalid operation -?")
    {
        bool early_terminate;
        REQUIRE_FALSE(commander.parse_command_wrapper("-?", early_terminate));
        REQUIRE_FALSE(early_terminate);
    }
    SECTION("Print version with --version")
    {
        bool early_terminate;
        REQUIRE_FALSE(
            commander.parse_command_wrapper("--version", early_terminate));
        REQUIRE(early_terminate);
    }
    SECTION("Print version with -v")
    {
        bool early_terminate;
        REQUIRE_FALSE(commander.parse_command_wrapper("-v", early_terminate));
        REQUIRE(early_terminate);
    }
    SECTION("Check barlevel")
    {
        SECTION("Valid bar level")
        {
            REQUIRE(commander.parse_command_wrapper(
                "--bar-levels 0.1,0.2,0.3,0.4,0.5"));
            auto res = commander.get_p_threshold().bar_levels;
            REQUIRE_THAT(res, Catch::Equals<double>({0.1, 0.2, 0.3, 0.4, 0.5}));
        }
        SECTION("Duplicated bar levels")
        {
            REQUIRE(commander.parse_command_wrapper(
                "--bar-levels 0.1,0.2,0.3,0.3,0.4,0.5"));
            auto res = commander.get_p_threshold().bar_levels;
            REQUIRE_THAT(res,
                         Catch::Equals<double>({0.1, 0.2, 0.3, 0.3, 0.4, 0.5}));
        }
        SECTION("Unsorted")
        {
            REQUIRE(commander.parse_command_wrapper(
                "--bar-levels 0.5,0.2,0.3,0.3,0.4,0.1"));
            auto res = commander.get_p_threshold().bar_levels;
            REQUIRE_THAT(res,
                         Catch::Equals<double>({0.5, 0.2, 0.3, 0.3, 0.4, 0.1}));
        }
        SECTION("Negative")
        {
            // supposed to fail but init didn't check for these either
            REQUIRE(commander.parse_command_wrapper(
                "--bar-levels 0.1,-0.2,0.3,0.4,0.5"));
            auto res = commander.get_p_threshold().bar_levels;
            REQUIRE_THAT(res,
                         Catch::Equals<double>({0.1, -0.2, 0.3, 0.4, 0.5}));
        }

        SECTION("Contain Zero")
        {
            // supposed to fail but init didn't check for these either
            REQUIRE(commander.parse_command_wrapper(
                "--bar-levels 0,0.2,0.3,0.3,0.4,0.5"));
            auto res = commander.get_p_threshold().bar_levels;
            REQUIRE_THAT(res,
                         Catch::Equals<double>({0, 0.2, 0.3, 0.3, 0.4, 0.5}));
        }
        SECTION("Out of p-value range")
        {
            // supposed to fail but init didn't check for these either
            REQUIRE(commander.parse_command_wrapper(
                "--bar-levels 0.5,0.2,0.3,0.3,0.4,0.1,2"));
            auto res = commander.get_p_threshold().bar_levels;
            REQUIRE_THAT(
                res, Catch::Equals<double>({0.5, 0.2, 0.3, 0.3, 0.4, 0.1, 2}));
        }
        SECTION("Invalid character")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper(
                "--bar-levels 0.1,0.2,0.3,a,0.4,0.5"));
        }
        SECTION("Overflow")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper(
                "--bar-levels 0.1,0.2,0.3,0.4,1.79769e+309"));
        }
    }
    SECTION("Set flags")
    {
        SECTION("fastscore")
        {
            REQUIRE_FALSE(commander.get_p_threshold().fastscore);
            REQUIRE(commander.parse_command_wrapper("--fastscore"));
            REQUIRE(commander.get_p_threshold().fastscore);
        }
        SECTION("no-full")
        {
            REQUIRE_FALSE(commander.get_p_threshold().no_full);
            REQUIRE(commander.parse_command_wrapper("--no-full"));
            REQUIRE(commander.get_p_threshold().no_full);
        }
        SECTION("no-clump")
        {
            REQUIRE_FALSE(commander.get_clump_info().no_clump);
            REQUIRE(commander.parse_command_wrapper("--no-clump"));
            REQUIRE(commander.get_clump_info().no_clump);
        }
        SECTION("hard-coded")
        {
            REQUIRE_FALSE(commander.get_target().hard_coded);
            REQUIRE(commander.parse_command_wrapper("--hard"));
            REQUIRE(commander.get_target().hard_coded);
        }
        SECTION("allow-inter")
        {
            REQUIRE_FALSE(commander.use_inter());
            REQUIRE(commander.parse_command_wrapper("--allow-inter"));
            REQUIRE(commander.use_inter());
        }
        SECTION("nonfounders")
        {
            REQUIRE_FALSE(commander.nonfounders());
            REQUIRE(commander.parse_command_wrapper("--nonfounders"));
            REQUIRE(commander.nonfounders());
        }
        SECTION("beta")
        {
            REQUIRE_FALSE(commander.get_base().is_beta);
            REQUIRE_FALSE(commander.get_base().is_or);
            REQUIRE(commander.parse_command_wrapper("--beta"));
            REQUIRE(commander.get_base().is_beta);
            REQUIRE_FALSE(commander.get_base().is_or);
        }
        SECTION("or")
        {
            REQUIRE_FALSE(commander.get_base().is_beta);
            REQUIRE_FALSE(commander.get_base().is_or);
            REQUIRE(commander.parse_command_wrapper("--or"));
            REQUIRE_FALSE(commander.get_base().is_beta);
            REQUIRE(commander.get_base().is_or);
        }
        SECTION("index")
        {
            REQUIRE_FALSE(commander.get_base().is_index);
            REQUIRE(commander.parse_command_wrapper("--index"));
            REQUIRE(commander.get_base().is_index);
        }
        SECTION("ignore-fid")
        {
            REQUIRE_FALSE(commander.get_pheno().ignore_fid);
            REQUIRE(commander.parse_command_wrapper("--ignore-fid"));
            REQUIRE(commander.get_pheno().ignore_fid);
        }

        SECTION("keep-ambig")
        {
            REQUIRE_FALSE(commander.keep_ambig());
            REQUIRE(commander.parse_command_wrapper("--keep-ambig"));
            REQUIRE(commander.keep_ambig());
        }

        SECTION("non-cumulate")
        {
            REQUIRE_FALSE(commander.get_prs_instruction().non_cumulate);
            REQUIRE(commander.parse_command_wrapper("--non-cumulate"));
            REQUIRE(commander.get_prs_instruction().non_cumulate);
        }
        SECTION("no-regress")
        {
            REQUIRE_FALSE(commander.get_prs_instruction().no_regress);
            REQUIRE(commander.parse_command_wrapper("--no-regress"));
            REQUIRE(commander.get_prs_instruction().no_regress);
        }
        SECTION("print-snp")
        {
            REQUIRE_FALSE(commander.print_snp());
            REQUIRE(commander.parse_command_wrapper("--print-snp"));
            REQUIRE(commander.print_snp());
        }
        SECTION("use-ref-maf")
        {
            REQUIRE_FALSE(commander.get_prs_instruction().use_ref_maf);
            REQUIRE(commander.parse_command_wrapper("--use-ref-maf"));
            REQUIRE(commander.get_prs_instruction().use_ref_maf);
        }
        SECTION("logit-perm")
        {
            REQUIRE_FALSE(commander.get_perm().logit_perm);
            REQUIRE(commander.parse_command_wrapper("--logit-perm"));
            REQUIRE(commander.get_perm().logit_perm);
        }
        SECTION("full-back")
        {
            REQUIRE_FALSE(commander.get_set().full_as_background);
            REQUIRE(commander.parse_command_wrapper("--full-back"));
            REQUIRE(commander.get_set().full_as_background);
        }
        SECTION("no-default")
        {
            REQUIRE_FALSE(commander.no_default());
            REQUIRE(commander.parse_command_wrapper("--no-default"));
            REQUIRE(commander.no_default());
        }
    }
    SECTION("Set base flags")
    {
        SECTION("A1")
        {
            check_set_base_flag(commander, "--A1", "a", "A1",
                                +BASE_INDEX::EFFECT);
        }
        SECTION("a1")
        {
            check_set_base_flag(commander, "--a1", "b", "A1",
                                +BASE_INDEX::EFFECT);
        }
        SECTION("A2")
        {
            check_set_base_flag(commander, "--A2", "c", "A2",
                                +BASE_INDEX::NONEFFECT);
        }
        SECTION("a2")
        {
            check_set_base_flag(commander, "--a2", "d", "A2",
                                +BASE_INDEX::NONEFFECT);
        }
        SECTION("stat")
        {
            check_set_base_flag(commander, "--stat", "statistic", "",
                                +BASE_INDEX::STAT);
        }
        SECTION("pvalue")
        {
            check_set_base_flag(commander, "--pvalue", "insignificant", "P",
                                +BASE_INDEX::P);
        }
        SECTION("p")
        {
            check_set_base_flag(commander, "-p", "postdoc", "P",
                                +BASE_INDEX::P);
        }
        SECTION("chr")
        {
            check_set_base_flag(commander, "--chr", "chromosome", "CHR",
                                +BASE_INDEX::CHR);
        }
        SECTION("bp")
        {
            check_set_base_flag(commander, "--bp", "location", "BP",
                                +BASE_INDEX::BP);
        }
        SECTION("snp")
        {
            check_set_base_flag(commander, "--snp", "cnv", "SNP",
                                +BASE_INDEX::RS);
        }
        SECTION("base")
        {
            REQUIRE(commander.get_base().file_name.empty());
            REQUIRE(commander.parse_command_wrapper(("--base BaseInfo")));
            REQUIRE(commander.get_base().file_name == "BaseInfo");
        }
        SECTION("b")
        {
            REQUIRE(commander.get_base().file_name.empty());
            REQUIRE(commander.parse_command_wrapper(("-b Basic")));
            REQUIRE(commander.get_base().file_name == "Basic");
        }
        SECTION("base-info")
        {
            check_set_base_flag(commander, "--base-info", "INFO_FILTER",
                                "INFO:0.9", +BASE_INDEX::INFO);
        }
        SECTION("base-maf")
        {
            check_set_base_flag(commander, "--base-maf", "MAF_FILTER", "",
                                +BASE_INDEX::MAF);
        }
    }
    SECTION("Parse binary target")
    {
        SECTION("Invalid")
        {

            SECTION("1 as T")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--binary-target 1"));
            }
            SECTION("0 as F")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--binary-target 0"));
            }
            SECTION("Out bound")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--binary-target 1e200T"));
            }
            SECTION("Mis-spell")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--binary-target Tru"));
            }
            SECTION("non-numeric")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--binary-target aT"));
            }
            SECTION("negative (-1)")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--binary-target -1T"));
            }
            SECTION("mixed negative (-10)")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--binary-target F,-10T"));
            }
            SECTION("negative and space")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--binary-target F, -10T"));
            }
        }
        SECTION("Grey area - space")
        {
            REQUIRE(commander.get_pheno().binary.empty());
            // 3F will be ignored
            REQUIRE(
                commander.parse_command_wrapper("--binary-target F,2T, 3F"));
            REQUIRE_THAT(commander.get_pheno().binary,
                         Catch::Equals<bool>({false, true, true}));
        }
        SECTION("Valid")
        {
            SECTION("T")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(commander.parse_command_wrapper("--binary-target T"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({true}));
            }
            SECTION("True")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(
                    commander.parse_command_wrapper("--binary-target True"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({true}));
            }
            SECTION("true")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(
                    commander.parse_command_wrapper("--binary-target true"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({true}));
            }
            SECTION("1true")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(
                    commander.parse_command_wrapper("--binary-target 1true"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({true}));
            }
            SECTION("1T")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(commander.parse_command_wrapper("--binary-target 1T"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({true}));
            }
            SECTION("F")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(commander.parse_command_wrapper("--binary-target F"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({false}));
            }
            SECTION("False")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(
                    commander.parse_command_wrapper("--binary-target False"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({false}));
            }
            SECTION("1false")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(
                    commander.parse_command_wrapper("--binary-target 1false"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({false}));
            }

            SECTION("4T")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(commander.parse_command_wrapper("--binary-target 4T"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({true, true, true, true}));
            }
            SECTION("6F")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                // 3F will be ignored
                REQUIRE(commander.parse_command_wrapper("--binary-target 6F"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>(
                                 {false, false, false, false, false, false}));
            }
            SECTION("True,3F")
            {
                REQUIRE(commander.get_pheno().binary.empty());
                REQUIRE(
                    commander.parse_command_wrapper("--binary-target True,3F"));
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({true, false, false, false}));
            }
        }
        SECTION("append check")
        {
            // we should append the vector

            REQUIRE(commander.get_pheno().binary.empty());
            REQUIRE(commander.parse_command_wrapper("--binary-target T,F"));
            REQUIRE_THAT(commander.get_pheno().binary,
                         Catch::Equals<bool>({true, false}));
            REQUIRE(commander.parse_command_wrapper("--binary-target F,T"));
            REQUIRE_THAT(commander.get_pheno().binary,
                         Catch::Equals<bool>({true, false, false, true}));
        }
    }
    SECTION("Set dose-threshold")
    {
        // make sure we specify the correct threshold
        REQUIRE(commander.get_target_qc().dose_threshold == 0.0);
        REQUIRE(commander.get_target_qc().hard_threshold == Approx(0.1));
        REQUIRE(commander.parse_command_wrapper("--dose-thres 1.0"));
        REQUIRE(commander.get_target_qc().dose_threshold == Approx(1.0));
        REQUIRE(commander.get_target_qc().hard_threshold == Approx(0.1));
    }
    SECTION("Set hard-thres")
    {
        REQUIRE(commander.get_target_qc().dose_threshold == 0.0);
        REQUIRE(commander.get_target_qc().hard_threshold == Approx(0.1));
        REQUIRE(commander.parse_command_wrapper("--hard-thres -0.1"));
        REQUIRE(commander.get_target_qc().dose_threshold == 0.0);
        REQUIRE(commander.get_target_qc().hard_threshold == Approx(-0.1));
    }
    SECTION("hard-thres out of bound")
    {
        REQUIRE_FALSE(commander.parse_command_wrapper("--hard-thres 1e400"));
    }
    SECTION("Reference dose-thres set")
    {
        REQUIRE(commander.get_ref_qc().dose_threshold == 0.0);
        REQUIRE(commander.get_ref_qc().hard_threshold == Approx(0.1));
        REQUIRE(commander.parse_command_wrapper("--ld-dose-thres 1.0"));
        REQUIRE(commander.get_ref_qc().dose_threshold == Approx(1.0));
        REQUIRE(commander.get_ref_qc().hard_threshold == Approx(0.1));
    }
    SECTION("Reference hard-thres set")
    {
        REQUIRE(commander.get_ref_qc().dose_threshold == 0.0);
        REQUIRE(commander.get_ref_qc().hard_threshold == Approx(0.1));
        REQUIRE(commander.parse_command_wrapper("--ld-hard-thres -0.1"));
        REQUIRE(commander.get_ref_qc().dose_threshold == 0.0);
        REQUIRE(commander.get_ref_qc().hard_threshold == Approx(-0.1));
    }
    SECTION("Target default")
    {
        SECTION("is_ref") { REQUIRE_FALSE(commander.get_target().is_ref); }
        SECTION("geno")
        {
            REQUIRE(commander.get_target_qc().geno == Approx(1.0));
        }
        SECTION("info")
        {
            REQUIRE(commander.get_target_qc().info_score == 0.0);
        }
        SECTION("maf") { REQUIRE(commander.get_target_qc().maf == 0.0); }
        SECTION("name") { REQUIRE(commander.get_target().file_name.empty()); }
        SECTION("name-list")
        {
            REQUIRE(commander.get_target().file_list.empty());
        }
        SECTION("keep") { REQUIRE(commander.get_target().keep.empty()); }
        SECTION("remove") { REQUIRE(commander.get_target().remove.empty()); }
        SECTION("type") { REQUIRE(commander.get_target().type == "bed"); }
    }
    SECTION("Pheno default")
    {
        SECTION("name") { REQUIRE(commander.get_pheno().pheno_file.empty()); }
        SECTION("name column")
        {
            REQUIRE(commander.get_pheno().pheno_col.empty());
        }
        SECTION("prevalence")
        {
            REQUIRE(commander.get_pheno().prevalence.empty());
        }
    }
    SECTION("Clump default")
    {
        SECTION("is_ref") { REQUIRE(commander.get_reference().is_ref); }
        SECTION("geno") { REQUIRE(commander.get_ref_qc().geno == Approx(1.0)); }
        SECTION("info") { REQUIRE(commander.get_ref_qc().info_score == 0.0); }
        SECTION("maf") { REQUIRE(commander.get_ref_qc().maf == 0.0); }
        SECTION("remove") { REQUIRE(commander.get_reference().remove.empty()); }
        SECTION("keep") { REQUIRE(commander.get_reference().keep.empty()); }
        SECTION("name")
        {
            REQUIRE(commander.get_reference().file_name.empty());
        }
        SECTION("list")
        {
            REQUIRE(commander.get_reference().file_list.empty());
        }
        SECTION("type") { REQUIRE(commander.get_reference().type == "bed"); }
        SECTION("r2") { REQUIRE(commander.get_clump_info().r2 == Approx(0.1)); }
        SECTION("proxy") { REQUIRE(commander.get_clump_info().proxy == 0.0); }
        SECTION("distance")
        {
            REQUIRE(commander.get_clump_info().distance == 250000);
        }
        SECTION("provided_distance")
        {
            REQUIRE_FALSE(commander.get_clump_info().provided_distance);
        }
        SECTION("pvalue")
        {
            REQUIRE(commander.get_clump_info().pvalue == Approx(1.0));
        }
    }
    SECTION("Set clump")
    {
        SECTION("clump-p")
        {
            REQUIRE(commander.parse_command_wrapper("--clump-p 0.1"));
            REQUIRE(commander.get_clump_info().pvalue == Approx(0.1));
        }
        SECTION("clump-r2")
        {
            REQUIRE(commander.parse_command_wrapper("--clump-r2 0.5"));
            REQUIRE(commander.get_clump_info().r2 == Approx(0.5));
        }
        SECTION("clump-kb without suffix")
        {
            REQUIRE(commander.parse_command_wrapper("--clump-kb 100"));
            REQUIRE(commander.get_clump_info().distance == 100000);
            REQUIRE(commander.get_clump_info().provided_distance);
        }
        SECTION("clump-kb with suffix (b)")
        {
            REQUIRE(commander.parse_command_wrapper("--clump-kb 100b"));
            REQUIRE(commander.get_clump_info().distance == 100);
            REQUIRE(commander.get_clump_info().provided_distance);
        }
        SECTION("clump-kb with suffix (kb)")
        {
            REQUIRE(commander.parse_command_wrapper("--clump-kb 100kb"));
            REQUIRE(commander.get_clump_info().distance == 100000);
            REQUIRE(commander.get_clump_info().provided_distance);
        }
        SECTION("clump-kb with suffix (mb)")
        {
            REQUIRE(commander.parse_command_wrapper("--clump-kb 200mb"));
            REQUIRE(commander.get_clump_info().distance == 200000000);
            REQUIRE(commander.get_clump_info().provided_distance);
        }
        SECTION("Negative clump-kb")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper("--clump-kb -100kb"));
        }
    }
    SECTION("PRSet defaults")
    {
        SECTION("wind-3") { REQUIRE(commander.get_set().wind_3 == 0); }
        SECTION("wind-5") { REQUIRE(commander.get_set().wind_5 == 0); }
        SECTION("background")
        {
            REQUIRE(commander.get_set().background.empty());
        }
        SECTION("msigdb") { REQUIRE(commander.get_set().msigdb.empty()); }
        SECTION("bed") { REQUIRE(commander.get_set().bed.empty()); }
        SECTION("snp") { REQUIRE(commander.get_set().snp.empty()); }
        SECTION("feature") { REQUIRE(commander.get_set().feature.empty()); }
        SECTION("gtf") { REQUIRE(commander.get_set().gtf.empty()); }
        SECTION("x-range") { REQUIRE(commander.exclusion_range().empty()); }
        SECTION("run") { REQUIRE_FALSE(commander.get_set().run); }
    }
    SECTION("Set PRSet parameters")
    {
        // default for window is bp
        SECTION("wind-3")
        {
            REQUIRE(commander.parse_command_wrapper("--wind-3 20"));
            REQUIRE(commander.get_set().wind_3 == 20);
        }
        SECTION("wind-5")
        {
            REQUIRE(commander.parse_command_wrapper("--wind-5 10"));
            REQUIRE(commander.get_set().wind_5 == 10);
        }
        SECTION("x-range single")
        {
            REQUIRE(commander.parse_command_wrapper("--x-range chr6:1-10"));
            REQUIRE(commander.exclusion_range() == "chr6:1-10");
        }
        SECTION("x-range multi")
        {
            // we don't tokenize nor do we check validity
            REQUIRE(commander.parse_command_wrapper(
                "--x-range chr6:1-10,chr22:133:288"));
            REQUIRE(commander.exclusion_range() == "chr6:1-10,chr22:133:288");
        }
        SECTION("background")
        {
            REQUIRE(commander.parse_command_wrapper("--background Name:0"));
            REQUIRE(commander.get_set().background == "Name:0");
        }
        SECTION("msigdb")
        {
            REQUIRE(commander.parse_command_wrapper("--msigdb kegg"));
            REQUIRE_THAT(commander.get_set().msigdb,
                         Catch::Equals<std::string>({"kegg"}));
            // it appends
            REQUIRE(commander.parse_command_wrapper("-m Reactome,MP"));
            REQUIRE_THAT(
                commander.get_set().msigdb,
                Catch::Equals<std::string>({"kegg", "Reactome", "MP"}));
        }
        SECTION("gtf")
        {
            REQUIRE(commander.parse_command_wrapper("--gtf Homo"));
            REQUIRE(commander.get_set().gtf == "Homo");
        }
        SECTION("g")
        {
            REQUIRE(commander.parse_command_wrapper("-g Misc"));
            REQUIRE(commander.get_set().gtf == "Misc");
        }
        SECTION("bed")
        {
            REQUIRE(commander.parse_command_wrapper("--bed File:Misc"));
            REQUIRE_THAT(commander.get_set().bed,
                         Catch::Equals<std::string>({"File:Misc"}));
            // also append
            REQUIRE(commander.parse_command_wrapper("-B Something,ok"));
            REQUIRE_THAT(
                commander.get_set().bed,
                Catch::Equals<std::string>({"File:Misc", "Something", "ok"}));
        }
        SECTION("snp-set")
        {
            REQUIRE(commander.parse_command_wrapper("--snp-set list,of,snp"));
            REQUIRE_THAT(commander.get_set().snp,
                         Catch::Equals<std::string>({"list", "of", "snp"}));
        }
        SECTION("feature")
        {
            REQUIRE(commander.parse_command_wrapper("--feature gene"));
            REQUIRE_THAT(commander.get_set().feature,
                         Catch::Equals<std::string>({"gene"}));
            // we don't check for duplicates here=
            REQUIRE(commander.parse_command_wrapper("--feature protein,gene"));
            REQUIRE_THAT(
                commander.get_set().feature,
                Catch::Equals<std::string>({"gene", "protein", "gene"}));
        }
    }
}
