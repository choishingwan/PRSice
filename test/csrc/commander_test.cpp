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
    void prepare_header_cov_check_wrapper(
        const std::vector<std::string>& cov_header,
        std::unordered_map<std::string, size_t>& ref_index,
        std::unordered_set<std::string>& included)
    {
        prepare_header_cov_check(cov_header, ref_index, included);
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
    SECTION("Target")
    {
        SECTION("default")
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
            SECTION("name")
            {
                REQUIRE(commander.get_target().file_name.empty());
            }
            SECTION("name-list")
            {
                REQUIRE(commander.get_target().file_list.empty());
            }
            SECTION("keep") { REQUIRE(commander.get_target().keep.empty()); }
            SECTION("remove")
            {
                REQUIRE(commander.get_target().remove.empty());
            }
            SECTION("type") { REQUIRE(commander.get_target().type == "bed"); }
        }
        SECTION("Set")
        {
            SECTION("name")
            {
                REQUIRE(commander.parse_command_wrapper("--target genotype"));
                REQUIRE(commander.get_target().file_name == "genotype");
            }
            SECTION("list")
            {
                REQUIRE(
                    commander.parse_command_wrapper("--target-list testing"));
                REQUIRE(commander.get_target().file_list == "testing");
            }
            SECTION("type")
            {
                REQUIRE(commander.parse_command_wrapper("--type bgen"));
                REQUIRE(commander.get_target().type == "bgen");
            }
            SECTION("keep")
            {
                REQUIRE(commander.parse_command_wrapper("--keep fun"));
                REQUIRE(commander.get_target().keep == "fun");
            }
            SECTION("remove")
            {
                REQUIRE(commander.parse_command_wrapper("--remove depression"));
                REQUIRE(commander.get_target().remove == "depression");
            }
            SECTION("dose-thres")
            {
                REQUIRE(commander.get_target_qc().dose_threshold == 0.0);
                REQUIRE(commander.get_target_qc().hard_threshold
                        == Approx(0.1));
                REQUIRE(commander.parse_command_wrapper("--dose-thres 1.0"));
                REQUIRE(commander.get_target_qc().dose_threshold
                        == Approx(1.0));
                REQUIRE(commander.get_target_qc().hard_threshold
                        == Approx(0.1));
            }

            SECTION("hard-thres")
            {
                REQUIRE(commander.get_target_qc().dose_threshold == 0.0);
                REQUIRE(commander.get_target_qc().hard_threshold
                        == Approx(0.1));
                REQUIRE(commander.parse_command_wrapper("--hard-thres -0.1"));
                REQUIRE(commander.get_target_qc().dose_threshold == 0.0);
                REQUIRE(commander.get_target_qc().hard_threshold
                        == Approx(-0.1));
            }
            SECTION("hard-thres out of bound")
            {
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--hard-thres 1e400"));
            }
        }
        SECTION("Filtering")
        {
            SECTION("geno")
            {
                REQUIRE(commander.parse_command_wrapper("--geno 0.4"));
                REQUIRE(commander.get_target_qc().geno == Approx(0.4));
            }
            SECTION("info")
            {
                REQUIRE(commander.parse_command_wrapper("--info 0.2"));
                REQUIRE(commander.get_target_qc().info_score == Approx(0.2));
            }
            SECTION("maf")
            {
                REQUIRE(commander.parse_command_wrapper("--maf 0.01"));
                REQUIRE(commander.get_target_qc().maf == Approx(0.01));
            }
        }
    }
    SECTION("Phenotype")
    {
        SECTION("default")
        {
            SECTION("name")
            {
                REQUIRE(commander.get_pheno().pheno_file.empty());
            }
            SECTION("name column")
            {
                REQUIRE(commander.get_pheno().pheno_col.empty());
            }
            SECTION("prevalence")
            {
                REQUIRE(commander.get_pheno().prevalence.empty());
            }
        }
        SECTION("set parameter")
        {
            SECTION("name")
            {
                REQUIRE(commander.parse_command_wrapper("--pheno Phenotype"));
                REQUIRE(commander.get_pheno().pheno_file == "Phenotype");
            }
            SECTION("pheno col")
            {
                REQUIRE(commander.parse_command_wrapper("--pheno-col A1,B1"));
                REQUIRE_THAT(commander.get_pheno().pheno_col,
                             Catch::Equals<std::string>({"A1", "B1"}));
                // append
                REQUIRE(commander.parse_command_wrapper("--pheno-col C1,D2"));
                REQUIRE_THAT(
                    commander.get_pheno().pheno_col,
                    Catch::Equals<std::string>({"A1", "B1", "C1", "D2"}));
            }
            SECTION("prevalence")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--prevalence 0.1,0.4,1,-0.5"));
                REQUIRE_THAT(commander.get_pheno().prevalence,
                             Catch::Equals<double>({0.1, 0.4, 1, -0.5}));
            }
        }
    }
    SECTION("Clump")
    {
        SECTION("default")
        {
            SECTION("is_ref") { REQUIRE(commander.get_reference().is_ref); }
            SECTION("geno")
            {
                REQUIRE(commander.get_ref_qc().geno == Approx(1.0));
            }
            SECTION("info")
            {
                REQUIRE(commander.get_ref_qc().info_score == 0.0);
            }
            SECTION("maf") { REQUIRE(commander.get_ref_qc().maf == 0.0); }
            SECTION("remove")
            {
                REQUIRE(commander.get_reference().remove.empty());
            }
            SECTION("keep") { REQUIRE(commander.get_reference().keep.empty()); }
            SECTION("name")
            {
                REQUIRE(commander.get_reference().file_name.empty());
            }
            SECTION("list")
            {
                REQUIRE(commander.get_reference().file_list.empty());
            }
            SECTION("type")
            {
                REQUIRE(commander.get_reference().type == "bed");
            }
            SECTION("r2")
            {
                REQUIRE(commander.get_clump_info().r2 == Approx(0.1));
            }
            SECTION("proxy")
            {
                REQUIRE(commander.get_clump_info().proxy == 0.0);
            }
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
        SECTION("Set parameter")
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
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--clump-kb -100kb"));
            }
        }
    }
    SECTION("PRSet")
    {
        SECTION("default")
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
        SECTION("Set parameter")
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
                REQUIRE(commander.exclusion_range()
                        == "chr6:1-10,chr22:133:288");
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
                REQUIRE_THAT(commander.get_set().bed,
                             Catch::Equals<std::string>(
                                 {"File:Misc", "Something", "ok"}));
            }
            SECTION("snp-set")
            {
                REQUIRE(
                    commander.parse_command_wrapper("--snp-set list,of,snp"));
                REQUIRE_THAT(commander.get_set().snp,
                             Catch::Equals<std::string>({"list", "of", "snp"}));
            }
            SECTION("feature")
            {
                REQUIRE(commander.parse_command_wrapper("--feature gene"));
                REQUIRE_THAT(commander.get_set().feature,
                             Catch::Equals<std::string>({"gene"}));
                // we don't check for duplicates here=
                REQUIRE(
                    commander.parse_command_wrapper("--feature protein,gene"));
                REQUIRE_THAT(
                    commander.get_set().feature,
                    Catch::Equals<std::string>({"gene", "protein", "gene"}));
            }
        }
    }
    SECTION("Set misc parameters")
    {
        SECTION("out")
        {
            REQUIRE(commander.out() == "PRSice");
            REQUIRE(commander.parse_command_wrapper("--out PRSet"));
            REQUIRE(commander.out() == "PRSet");
        }
        SECTION("thread")
        {
            REQUIRE(commander.get_prs_instruction().thread == 1);
            auto max_thread = commander.max_thread();
            SECTION("simple")
            {
                REQUIRE(commander.parse_command_wrapper("--thread 2"));
                REQUIRE(commander.get_prs_instruction().thread == 2);
            }
            SECTION("max")
            {
                REQUIRE(commander.parse_command_wrapper("--thread max"));
                REQUIRE(commander.get_prs_instruction().thread == max_thread);
            }
            SECTION("Too large")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--thread " + std::to_string(2 * max_thread)));
                REQUIRE(commander.get_prs_instruction().thread == max_thread);
            }
        }
        SECTION("extract")
        {
            REQUIRE(commander.extract_file().empty());
            REQUIRE(commander.parse_command_wrapper("--extract Love"));
            REQUIRE(commander.extract_file() == "Love");
        }
        SECTION("exclude")
        {
            REQUIRE(commander.exclude_file().empty());
            REQUIRE(commander.parse_command_wrapper("--exclude hate"));
            REQUIRE(commander.exclude_file() == "hate");
        }
        SECTION("delim")
        {
            REQUIRE(commander.delim() == " ");
            REQUIRE(commander.parse_command_wrapper("--id-delim -"));
            REQUIRE(commander.delim() == "-");
        }
        SECTION("memory")
        {
            REQUIRE(commander.memory() == 1e10);
            SECTION("suffix k")
            {
                REQUIRE(commander.parse_command_wrapper("--memory 1k"));
                REQUIRE(commander.memory() == 1024);
            }
            SECTION("suffix gb")
            {
                REQUIRE(commander.parse_command_wrapper("--memory 1gb"));
                REQUIRE(commander.memory() == 1073741824);
            }
            SECTION("suffix tb")
            {
                REQUIRE(commander.parse_command_wrapper("--memory 30tb"));
                REQUIRE(commander.memory() == 32985348833280);
            }
            SECTION("no suffix")
            {
                REQUIRE(commander.parse_command_wrapper("--memory 10"));
                REQUIRE(commander.memory() == 10485760);
            }
        }
        SECTION("seed")
        {
            // Default is random, so we won't check that
            REQUIRE(commander.parse_command_wrapper("--seed 123"));
            REQUIRE(commander.get_perm().seed == 123);
        }
        SECTION("num-autosome")
        {
            REQUIRE(commander.get_target().num_autosome == 22);
            REQUIRE(commander.get_reference().num_autosome == 22);
            SECTION("valid")
            {
                REQUIRE(commander.parse_command_wrapper("--num-auto 1"));
                REQUIRE(commander.get_target().num_autosome == 1);
                REQUIRE(commander.get_reference().num_autosome == 1);
            }
            SECTION("negative")
            {
                REQUIRE(commander.parse_command_wrapper("--num-auto -100"));
                REQUIRE(commander.get_target().num_autosome == -100);
                REQUIRE(commander.get_reference().num_autosome == -100);
            }
            SECTION("overflow")
            {
                REQUIRE_FALSE(
                    commander.parse_command_wrapper("--num-auto 1e100"));
            }
        }
    }
    SECTION("Set permutation")
    {
        REQUIRE(commander.get_perm().num_permutation == 0);
        REQUIRE_FALSE(commander.get_perm().run_perm);
        REQUIRE_FALSE(commander.get_perm().run_set_perm);
        SECTION("perm")
        {
            REQUIRE(commander.parse_command_wrapper("--perm 100"));
            REQUIRE(commander.get_perm().num_permutation == 100);
            REQUIRE(commander.get_perm().run_perm);
            REQUIRE_FALSE(commander.get_perm().run_set_perm);
        }
        SECTION("set perm")
        {
            REQUIRE(commander.parse_command_wrapper("--set-perm 1026"));
            REQUIRE(commander.get_perm().num_permutation == 1026);
            REQUIRE(commander.get_perm().run_set_perm);
            REQUIRE_FALSE(commander.get_perm().run_perm);
        }
        SECTION("overflow")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper("--set-perm 1e200"));
        }
    }
    SECTION("Set p-value thresholds")
    {
        REQUIRE(commander.get_p_threshold().inter == Approx(0.00005));
        REQUIRE(commander.get_p_threshold().lower == Approx(5e-8));
        REQUIRE(commander.get_p_threshold().upper == Approx(0.5));
        SECTION("inter")
        {
            REQUIRE(commander.parse_command_wrapper("--inter 1e-10"));
            REQUIRE(commander.get_p_threshold().inter == Approx(1e-10));
        }
        SECTION("i")
        {
            REQUIRE(commander.parse_command_wrapper("-i 110"));
            REQUIRE(commander.get_p_threshold().inter == Approx(110));
        }
        SECTION("lower")
        {
            REQUIRE(commander.parse_command_wrapper("--lower 1e-20"));
            REQUIRE(commander.get_p_threshold().lower == Approx(1e-20));
        }
        SECTION("l")
        {
            REQUIRE(commander.parse_command_wrapper("-l 123"));
            REQUIRE(commander.get_p_threshold().lower == Approx(123));
        }
        SECTION("upper")
        {
            REQUIRE(commander.parse_command_wrapper("--upper 5e10"));
            REQUIRE(commander.get_p_threshold().upper == Approx(5e10));
        }
        SECTION("u")
        {
            REQUIRE(commander.parse_command_wrapper("--upper 5e-70"));
            REQUIRE(commander.get_p_threshold().upper == Approx(5e-70));
        }
        SECTION("Very large but valid")
        {
            REQUIRE(commander.parse_command_wrapper("--upper 5e300"));
            REQUIRE(commander.get_p_threshold().upper == Approx(5e300));
        }
        SECTION("Very close to zero but valid")
        {
            REQUIRE(commander.parse_command_wrapper("--lower 5e-300"));
            REQUIRE(commander.get_p_threshold().lower == Approx(5e-300));
        }
        SECTION("Too large")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper("--lower 5e400"));
        }
        SECTION("Too closet to zero")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper("--lower 5e-400"));
        }
        SECTION("invalid input lower")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper("--lower lower"));
        }
        SECTION("invalid input upper")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper("--upper upper"));
        }
        SECTION("invalid input inter")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper("--inter inter"));
        }
    }
    SECTION("Score method")
    {
        REQUIRE(commander.get_prs_instruction().scoring_method
                == SCORING::AVERAGE);
        SECTION("set sum")
        {
            // test to capital works
            SECTION("Sum")
            {
                REQUIRE(commander.parse_command_wrapper("--score Sum"));
                REQUIRE(commander.get_prs_instruction().scoring_method
                        == SCORING::SUM);
            }
            SECTION("sUm")
            {
                REQUIRE(commander.parse_command_wrapper("--score sUm"));
                REQUIRE(commander.get_prs_instruction().scoring_method
                        == SCORING::SUM);
            }
        }
        SECTION("std")
        {
            REQUIRE(commander.parse_command_wrapper("--score std"));
            REQUIRE(commander.get_prs_instruction().scoring_method
                    == SCORING::STANDARDIZE);
        }
        SECTION("con-std")
        {
            REQUIRE(commander.parse_command_wrapper("--score con-std"));
            REQUIRE(commander.get_prs_instruction().scoring_method
                    == SCORING::CONTROL_STD);
        }
        SECTION("average")
        {
            REQUIRE(commander.parse_command_wrapper("--score avg"));
            REQUIRE(commander.get_prs_instruction().scoring_method
                    == SCORING::AVERAGE);
        }
        SECTION("require exact match")
        {
            REQUIRE_FALSE(commander.parse_command_wrapper("--score avgerage"));
        }
    }
    SECTION("missing")
    {
        REQUIRE(commander.get_prs_instruction().missing_score
                == MISSING_SCORE::MEAN_IMPUTE);
        SECTION("set_zero")
        {
            REQUIRE(commander.parse_command_wrapper("--missing Set_Zero"));
            REQUIRE(commander.get_prs_instruction().missing_score
                    == MISSING_SCORE::SET_ZERO);
        }
        SECTION("mean")
        {
            REQUIRE(commander.parse_command_wrapper("--missing mean_impute"));
            REQUIRE(commander.get_prs_instruction().missing_score
                    == MISSING_SCORE::MEAN_IMPUTE);
        }
        SECTION("Center")
        {
            REQUIRE(commander.parse_command_wrapper("--missing centrer"));
            REQUIRE(commander.get_prs_instruction().missing_score
                    == MISSING_SCORE::CENTER);
        }
        SECTION("impute control")
        {
            // allow parsing, but implementation not done yet
            REQUIRE(
                commander.parse_command_wrapper("--missing IMPUTE_CONTROL"));
            REQUIRE(commander.get_prs_instruction().missing_score
                    == MISSING_SCORE::IMPUTE_CONTROL);
        }
    }
    SECTION("model")
    {
        REQUIRE(commander.get_prs_instruction().genetic_model
                == MODEL::ADDITIVE);
        SECTION("dom")
        {
            REQUIRE(commander.parse_command_wrapper("--model dom"));
            REQUIRE(commander.get_prs_instruction().genetic_model
                    == MODEL::DOMINANT);
        }
        SECTION("add")
        {
            REQUIRE(commander.parse_command_wrapper("--model ADD"));
            REQUIRE(commander.get_prs_instruction().genetic_model
                    == MODEL::ADDITIVE);
        }
        SECTION("rec")
        {
            REQUIRE(commander.parse_command_wrapper("--model ReC"));
            REQUIRE(commander.get_prs_instruction().genetic_model
                    == MODEL::RECESSIVE);
        }
        SECTION("het")
        {
            REQUIRE(commander.parse_command_wrapper("--model het"));
            REQUIRE(commander.get_prs_instruction().genetic_model
                    == MODEL::HETEROZYGOUS);
        }
    }
    SECTION("Covariate")
    {
        REQUIRE(commander.get_pheno().cov_file.empty());
        REQUIRE(commander.get_pheno().cov_colname.empty());
        REQUIRE(commander.get_pheno().factor_cov.empty());
        SECTION("simple load")
        {
            REQUIRE(
                commander.parse_command_wrapper("--cov-col Testing,@PC[1-55]"));
            REQUIRE_THAT(commander.get_pheno().cov_colname,
                         Catch::Equals<std::string>({"Testing", "@PC[1-55]"}));
            SECTION("append")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--cov-col More,Covariate"));
                REQUIRE_THAT(commander.get_pheno().cov_colname,
                             Catch::Equals<std::string>({"Testing", "@PC[1-55]",
                                                         "More", "Covariate"}));
            }
        }
        SECTION("discontinue load")
        {
            REQUIRE(commander.parse_command_wrapper(
                "--cov-col Testing,@PC[1.3.5]"));
            REQUIRE_THAT(commander.get_pheno().cov_colname,
                         Catch::Equals<std::string>({"Testing", "@PC[1.3.5]"}));
        }
        SECTION("factor cov")
        {
            REQUIRE(commander.parse_command_wrapper("--cov-factor Sex"));
            REQUIRE_THAT(commander.get_pheno().factor_cov,
                         Catch::Equals<std::string>({"Sex"}));
        }
    }
    SECTION("Reference")
    {
        SECTION("default")
        {
            SECTION("is_ref") { REQUIRE(commander.get_reference().is_ref); }
            SECTION("name")
            {
                REQUIRE(commander.get_reference().file_name.empty());
            }
            SECTION("list")
            {
                REQUIRE(commander.get_reference().file_list.empty());
            }
            SECTION("keep") { REQUIRE(commander.get_reference().keep.empty()); }
            SECTION("remove")
            {
                REQUIRE(commander.get_reference().remove.empty());
            }
            SECTION("type")
            {
                REQUIRE(commander.get_reference().type == "bed");
            }
            SECTION("dose-thres")
            {
                REQUIRE(commander.get_ref_qc().dose_threshold == 0.0);
                REQUIRE(commander.get_ref_qc().hard_threshold == Approx(0.1));
                REQUIRE(commander.parse_command_wrapper("--ld-dose-thres 1.0"));
                REQUIRE(commander.get_ref_qc().dose_threshold == Approx(1.0));
                REQUIRE(commander.get_ref_qc().hard_threshold == Approx(0.1));
            }
            SECTION("hard-thres")
            {
                REQUIRE(commander.get_ref_qc().dose_threshold == 0.0);
                REQUIRE(commander.get_ref_qc().hard_threshold == Approx(0.1));
                REQUIRE(
                    commander.parse_command_wrapper("--ld-hard-thres -0.1"));
                REQUIRE(commander.get_ref_qc().dose_threshold == 0.0);
                REQUIRE(commander.get_ref_qc().hard_threshold == Approx(-0.1));
            }
        }
        SECTION("set parameter")
        {
            SECTION("name")
            {
                REQUIRE(commander.parse_command_wrapper("--ld geno"));
                REQUIRE(commander.get_reference().file_name == "geno");
            }
            SECTION("list")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-list test"));
                REQUIRE(commander.get_reference().file_list == "test");
            }
            SECTION("type")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-type bgen"));
                REQUIRE(commander.get_reference().type == "bgen");
            }
            SECTION("keep")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-keep fun"));
                REQUIRE(commander.get_reference().keep == "fun");
            }
            SECTION("remove")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-remove depress"));
                REQUIRE(commander.get_reference().remove == "depress");
            }
        }
        SECTION("filtering")
        {
            SECTION("geno")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-geno 0.4"));
                REQUIRE(commander.get_ref_qc().geno == Approx(0.4));
            }
            SECTION("info")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-info 0.2"));
                REQUIRE(commander.get_ref_qc().info_score == Approx(0.2));
            }
            SECTION("maf")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-maf 0.01"));
                REQUIRE(commander.get_ref_qc().maf == Approx(0.01));
            }
        }
    }
}


void transform_test(const std::string& cov,
                    const std::vector<std::string>& expected,
                    bool expect_success)
{
    std::vector<std::string> result;
    if (!expect_success)
    { REQUIRE_FALSE(mockCommander::transform_wrapper(cov, result)); }
    else
    {
        REQUIRE(mockCommander::transform_wrapper(cov, result));
        REQUIRE_THAT(result, Catch::UnorderedEquals<std::string>(expected));
    }
}
TEST_CASE("Covariate Transformation")
{
    SECTION("find_first_end")
    {
        std::string cov = "PC[1-5]";
        size_t res;
        SECTION("wrong start index")
        {
            REQUIRE_FALSE(mockCommander::find_first_end_wrapper(cov, 0, res));
        }
        SECTION("correct start")
        {
            REQUIRE(mockCommander::find_first_end_wrapper(cov, 2, res));
        }
        SECTION("nested list not allowed")
        {
            cov = "PC[1-5[1-5]]";
            SECTION("wrong start index")
            {
                REQUIRE_FALSE(
                    mockCommander::find_first_end_wrapper(cov, 0, res));
            }
            SECTION("correct start")
            {
                REQUIRE_FALSE(
                    mockCommander::find_first_end_wrapper(cov, 2, res));
            }
        }
    }
    SECTION("parse_range")
    {
        std::vector<size_t> res;
        SECTION("Expect [] removed")
        {
            REQUIRE_FALSE(mockCommander::parse_range_wrapper("[1-5]", res));
        }
        SECTION("valid ranges")
        {
            int start = GENERATE(take(10, random(1, 1000)));
            int end = GENERATE(take(10, random(1, 1000)));
            REQUIRE(mockCommander::parse_range_wrapper(
                std::to_string(start) + "-" + std::to_string(end), res));
            if (start > end) { std::swap(start, end); }
            std::vector<size_t> expected(end - start + 1, start);
            std::iota(expected.begin(), expected.end(), start);
            REQUIRE_THAT(res, Catch::Equals<size_t>(expected));
        }
        SECTION("start == end")
        {
            REQUIRE(mockCommander::parse_range_wrapper("1-1", res));
            REQUIRE_THAT(res, Catch::Equals<size_t>({1}));
        }
        SECTION("single value")
        {
            auto i = GENERATE(take(100, random(-100, 100)));
            if (i < 0)
            {
                REQUIRE_FALSE(
                    mockCommander::parse_range_wrapper(std::to_string(i), res));
            }
            else
            {
                REQUIRE(
                    mockCommander::parse_range_wrapper(std::to_string(i), res));
                REQUIRE_THAT(res,
                             Catch::Equals<size_t>({static_cast<size_t>(i)}));
            }
        }
        SECTION("not allow negative range")
        {
            REQUIRE_FALSE(mockCommander::parse_range_wrapper("1--5", res));
        }
        SECTION("not allow full negative")
        {
            REQUIRE_FALSE(mockCommander::parse_range_wrapper("-1--5", res));
        }
        SECTION("not allow comma")
        {
            REQUIRE_FALSE(mockCommander::parse_range_wrapper("1,5", res));
        }
    }
    SECTION("get_range")
    {
        std::vector<size_t> res;
        SECTION("invalid start")
        {
            REQUIRE_FALSE(
                mockCommander::get_range_wrapper("PC[1-5]", 0, 6, res));
        }
        SECTION("invalid end")
        {
            REQUIRE_FALSE(
                mockCommander::get_range_wrapper("PC[1-5]", 2, 4, res));
        }
        SECTION("out of bound end")
        {
            REQUIRE_FALSE(
                mockCommander::get_range_wrapper("PC[1-5]", 2, 7, res));
        }
        SECTION("Correct parse")
        {
            REQUIRE(mockCommander::get_range_wrapper("PC[1-5]", 2, 6, res));
            REQUIRE_THAT(res, Catch::Equals<size_t>({1, 2, 3, 4, 5}));
        }
        SECTION("mixed parse")
        {
            REQUIRE(
                mockCommander::get_range_wrapper("PC[1-5.8.7-10]", 2, 13, res));
            REQUIRE_THAT(res, Catch::UnorderedEquals<size_t>(
                                  {1, 2, 3, 4, 5, 7, 8, 9, 10}));
        }
        SECTION("double parse")
        {
            // should be disallowed, but can't detect, need to check log
            REQUIRE(mockCommander::get_range_wrapper("PC[1-5.6.0.005]", 2, 14,
                                                     res));
            REQUIRE_THAT(res,
                         Catch::UnorderedEquals<size_t>({0, 1, 2, 3, 4, 5, 6}));
        }
    }
    SECTION("update_covariate_ranges")
    {
        std::vector<std::string> result;
        std::vector<size_t> range;
        SECTION("empty vectors")
        {
            REQUIRE_FALSE(
                mockCommander::update_covariate_ranges_wrapper(result, range));
        }
        range = {1, 3, 5, 7, 9};
        SECTION("valid range")
        {
            REQUIRE(
                mockCommander::update_covariate_ranges_wrapper(result, range));
            std::vector<std::string> expected = {"1", "3", "5", "7", "9"};
            REQUIRE_THAT(result, Catch::UnorderedEquals<std::string>(expected));
        }
        SECTION("append range")
        {
            result = {"PC1AB", "PC2CD"};
            REQUIRE(
                mockCommander::update_covariate_ranges_wrapper(result, range));
            std::vector<std::string> expected = {
                "PC1AB1", "PC1AB3", "PC1AB5", "PC1AB7", "PC1AB9",
                "PC2CD1", "PC2CD3", "PC2CD5", "PC2CD7", "PC2CD9"};
            REQUIRE_THAT(result, Catch::UnorderedEquals<std::string>(expected));
        }
    }
    SECTION("transformation")
    {
        SECTION("without parsing")
        {
            transform_test("@PC", std::vector<std::string> {"PC"}, true);
        }
        SECTION("double @")
        {
            transform_test("@@PC", std::vector<std::string> {"@PC"}, true);
        }
        SECTION("without @")
        {
            // this is ok, as we'd removed the @ when we start transformation
            transform_test("PC1", std::vector<std::string> {"PC1"}, true);
        }
        SECTION("Range")
        {
            transform_test("PC1-5", std::vector<std::string> {"PC1-5"}, true);
        }
        SECTION("Internal @")
        {
            transform_test("PC1@5", std::vector<std::string> {"PC1@5"}, true);
        }
        SECTION("Without []")
        {
            // we will always remove first @
            transform_test("@PC1-5", std::vector<std::string> {"PC1-5"}, true);
        }
        SECTION("Normal transformation")
        {
            transform_test(
                "@PC[1-5]",
                std::vector<std::string> {"PC1", "PC2", "PC3", "PC4", "PC5"},
                true);
        }
        SECTION("Mixed tranform")
        {
            transform_test("@PC[1-2.5]",
                           std::vector<std::string> {"PC1", "PC2", "PC5"},
                           true);
        }
        SECTION("Complex transform")
        {
            transform_test("@PC[1-2.4.3-6]",
                           std::vector<std::string> {"PC1", "PC2", "PC3", "PC4",
                                                     "PC5", "PC6"},
                           true);
        }
        SECTION("Parse internal")
        {
            transform_test("@PC[1-2]A",
                           std::vector<std::string> {"PC1A", "PC2A"}, true);
        }
        SECTION("Multiple range")
        {
            transform_test(
                "@PC[1-2]A[1-2]",
                std::vector<std::string> {"PC1A1", "PC1A2", "PC2A1", "PC2A2"},
                true);
        }
    }
}

TEST_CASE("Unit parsing")
{
    mockCommander commander;
    size_t value = 0;
    SECTION("memory without number")
    {
        auto i = GENERATE("m", "b", "mb", "hi", "TB");
        REQUIRE_FALSE(commander.check_parse_unit_value(i, "--mem", 0, value));
    }
    SECTION("out of bound")
    {
        SECTION("high power")
        {
            REQUIRE_FALSE(commander.check_parse_unit_value("1", "", 7, value));
        }
        SECTION("negative input")
        {
            REQUIRE_FALSE(
                commander.check_parse_unit_value("1000000000tb", "", 1, value));
        }
    }
    SECTION("negative input")
    {
        SECTION("without unit")
        {
            REQUIRE_FALSE(commander.check_parse_unit_value("-1", "", 1, value));
        }
        SECTION("with unit")
        {
            REQUIRE_FALSE(
                commander.check_parse_unit_value("-1tb", "", 1, value));
        }
    }
    SECTION("unit provided")
    {
        SECTION("ignore default")
        {
            auto i = GENERATE(range(0, 7));
            REQUIRE(commander.check_parse_unit_value(
                "1b", "", static_cast<size_t>(i), value));
            REQUIRE(value == 1);
        }
        SECTION("different units")
        {
            using record = std::tuple<std::string, size_t>;
            auto expected = GENERATE(table<std::string, size_t>(
                {record {"1k", 1000}, record {"1kb", 1000},
                 record {"1m", 1000000}, record {"1mb", 1000000},
                 record {"1g", 1000000000}, record {"1gb", 1000000000},
                 record {"1t", 1000000000000}, record {"1tb", 1000000000000}}));
            auto input = std::get<0>(expected);
            auto res = std::get<1>(expected);
            REQUIRE(commander.check_parse_unit_value(input, "", 1, value));
            REQUIRE(value == res);
        }
        SECTION("non-integer")
        {
            using record = std::tuple<std::string, size_t>;
            auto expected = GENERATE(table<std::string, size_t>(
                {record {"1.5k", 1500}, record {"1.004kb", 1004}}));
            auto input = std::get<0>(expected);
            auto res = std::get<1>(expected);
            REQUIRE(commander.check_parse_unit_value(input, "", 1, value));
            REQUIRE(value == res);
            SECTION("result can't be double")
            {
                REQUIRE_FALSE(
                    commander.check_parse_unit_value("1.5b", "", 1, value));
            }
        }
    }
    SECTION("default value")
    {
        auto power = GENERATE(range(0, 6));
        auto exp = std::pow(1000, power);
        REQUIRE(commander.check_parse_unit_value(
            "1", "", static_cast<size_t>(power), value));
        REQUIRE(value == static_cast<size_t>(exp));
    }
}

TEST_CASE("Validation")
{
    mockCommander commander;
    SECTION("empty parameter")
    {
        REQUIRE_FALSE(commander.parse_command_wrapper(""));
    }
    SECTION("Target")
    {
        SECTION("name")
        {
            REQUIRE(commander.parse_command_wrapper("--target test"));
            REQUIRE(commander.target_check_wrapper());
        }
        SECTION("list")
        {
            REQUIRE(commander.parse_command_wrapper("--target-list test"));
            REQUIRE(commander.target_check_wrapper());
        }
        SECTION("both list and name")
        {
            REQUIRE(commander.parse_command_wrapper(
                "--target test --target-list list"));
            REQUIRE_FALSE(commander.target_check_wrapper());
        }
        SECTION("other target")
        {
            REQUIRE(commander.parse_command_wrapper("--target test"));
            SECTION("keep")
            {
                REQUIRE(commander.parse_command_wrapper("--keep more"));
                REQUIRE(commander.target_check_wrapper());
            }
            SECTION("remove")
            {
                REQUIRE(commander.parse_command_wrapper("--remove no-more"));
                REQUIRE(commander.target_check_wrapper());
            }
            SECTION("both keep and remove")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--keep more --remove no-more"));
                REQUIRE_FALSE(commander.target_check_wrapper());
            }
            SECTION("Valid type")
            {
                auto i = GENERATE("bgen", "bed", "ped");
                REQUIRE(commander.parse_command_wrapper("--type "
                                                        + std::string(i)));
                REQUIRE(commander.target_check_wrapper());
            }
            SECTION("universial bound check")
            {
                auto para = GENERATE("--geno ", "--maf ");
                auto value = GENERATE(take(100, random(-1.1, 1.1)));
                REQUIRE(commander.parse_command_wrapper(
                    std::string(para) + std::to_string(value)));
                if (value > 1 || value < 0)
                { REQUIRE_FALSE(commander.filter_check_wrapper()); }
                else
                {
                    REQUIRE(commander.filter_check_wrapper());
                }
            }
            SECTION("bgen bound check")
            {
                auto para =
                    GENERATE("--hard-thres ", "--dose-thres ", "--info ");
                auto value = GENERATE(take(100, random(-1.1, 1.1)));
                REQUIRE(commander.parse_command_wrapper(
                    "--type bgen " + std::string(para)
                    + std::to_string(value)));
                if (value > 1 || value < 0)
                { REQUIRE_FALSE(commander.filter_check_wrapper()); }
                else
                {
                    REQUIRE(commander.filter_check_wrapper());
                }
            }
            SECTION("invalid type")
            {
                REQUIRE(commander.parse_command_wrapper("--type vcf"));
                REQUIRE_FALSE(commander.target_check_wrapper());
            }
            SECTION("autosome")
            {
                auto i = GENERATE(take(100, random(-100, 100)));
                REQUIRE(commander.parse_command_wrapper("--num-auto "
                                                        + std::to_string(i)));
                if (i < 0) { REQUIRE_FALSE(commander.target_check_wrapper()); }
                else
                {
                    REQUIRE(commander.target_check_wrapper());
                }
            }
        }
    }
    SECTION("Thresholding")
    {
        SECTION("PRSice default")
        {
            SECTION("without no-full")
            {
                REQUIRE(commander.prsice_check_wrapper());
                REQUIRE_THAT(commander.get_p_threshold().bar_levels,
                             Catch::Approx<double>(
                                 {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1}));
            }
            SECTION("with no-full")
            {
                REQUIRE(commander.parse_command_wrapper("--no-full"));
                REQUIRE(commander.prsice_check_wrapper());
                REQUIRE_THAT(commander.get_p_threshold().bar_levels,
                             Catch::Approx<double>(
                                 {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5}));
            }
        }
        SECTION("PRSet default")
        {
            REQUIRE(commander.parse_command_wrapper("--msigdb PRSet"));
            SECTION("without thresholding")
            {
                REQUIRE(commander.prsice_check_wrapper());
                REQUIRE_THAT(commander.get_p_threshold().bar_levels,
                             Catch::Approx<double>({1}));
            }
            SECTION("with other parameters")
            {
                auto i = GENERATE("--fastscore", "--lower 1e-5", "--inter 1e-4",
                                  "--upper 0.6");
                REQUIRE(commander.parse_command_wrapper(i));
                REQUIRE(commander.prsice_check_wrapper());
                REQUIRE_THAT(commander.get_p_threshold().bar_levels,
                             Catch::Approx<double>(
                                 {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1}));
            }
            SECTION("no full")
            {
                // no effect
                REQUIRE(commander.prsice_check_wrapper());
                REQUIRE_THAT(commander.get_p_threshold().bar_levels,
                             Catch::Approx<double>({1}));
            }
        }
        SECTION("bar level")
        {
            auto i = GENERATE(take(100, random(-1.1, 1.1)));
            REQUIRE(commander.parse_command_wrapper("--bar-levels "
                                                    + std::to_string(i)));
            if (i <= 1.0 && i >= 0.0)
            { REQUIRE(commander.prsice_check_wrapper()); }
            else
            {
                REQUIRE_FALSE(commander.prsice_check_wrapper());
            }
        }
        SECTION("high res check")
        {
            auto para = GENERATE("--lower ", "--upper ", "--inter ");
            auto value = GENERATE(take(100, random(-1.1, 1.1)));
            REQUIRE(commander.parse_command_wrapper(std::string(para)
                                                    + std::to_string(value)));
            if (value > 1 || value < 0
                || commander.get_p_threshold().lower
                       > commander.get_p_threshold().upper
                || misc::logically_equal(commander.get_p_threshold().inter,
                                         0.0))
            { REQUIRE_FALSE(commander.prsice_check_wrapper()); }
            else
            {
                REQUIRE(commander.prsice_check_wrapper());
            }
        }
    }
    SECTION("Clump")
    {
        SECTION("PRSice default")
        {
            REQUIRE(commander.clump_check_wrapper());
            REQUIRE(commander.get_clump_info().distance == 250000);
        }
        SECTION("PRSet default")
        {
            // for prset, default is set to 1mb
            REQUIRE(commander.parse_command_wrapper("--msigdb RunPRSet"));
            REQUIRE(commander.clump_check_wrapper());
            REQUIRE(commander.get_clump_info().distance == 1000000);
        }
        SECTION("p, r2 and proxy")
        {
            // we just check the validation, set and get checked before
            auto para = GENERATE("--clump-p ", "--clump-r2 ", "--proxy ");
            auto value = GENERATE(take(100, random(-1.1, 1.1)));
            REQUIRE(commander.parse_command_wrapper(std::string(para)
                                                    + std::to_string(value)));
            if (value > 1 || value < 0)
            { REQUIRE_FALSE(commander.clump_check_wrapper()); }
            else
            {
                REQUIRE(commander.clump_check_wrapper());
                if (std::string(para) == "--proxy ")
                {
                    REQUIRE(commander.get_clump_info().proxy
                            == Approx(misc::Convertor::convert<double>(
                                std::to_string(value))));
                }
            }
        }
    }
    SECTION("Reference")
    {
        SECTION("name")
        {
            REQUIRE(commander.parse_command_wrapper("--ld test"));
            REQUIRE(commander.ref_check_wrapper());
        }
        SECTION("list")
        {
            REQUIRE(commander.parse_command_wrapper("--ld-list test"));
            REQUIRE(commander.ref_check_wrapper());
        }
        SECTION("both list and name")
        {
            REQUIRE(
                commander.parse_command_wrapper("--ld test --ld-list list"));
            REQUIRE_FALSE(commander.ref_check_wrapper());
        }
        SECTION("other target")
        {
            REQUIRE(commander.parse_command_wrapper("--ld test"));
            SECTION("keep")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-keep more"));
                REQUIRE(commander.ref_check_wrapper());
            }
            SECTION("remove")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-remove no-more"));
                REQUIRE(commander.ref_check_wrapper());
            }
            SECTION("both keep and remove")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--ld-keep more --ld-remove no-more"));
                REQUIRE_FALSE(commander.ref_check_wrapper());
            }
            SECTION("Valid type")
            {
                auto i = GENERATE("bgen", "bed", "ped");
                REQUIRE(commander.parse_command_wrapper("--ld-type "
                                                        + std::string(i)));
                REQUIRE(commander.ref_check_wrapper());
            }
            SECTION("universial bound check")
            {
                auto para = GENERATE("--ld-geno ", "--ld-maf ");
                auto value = GENERATE(take(100, random(-1.1, 1.1)));
                REQUIRE(commander.parse_command_wrapper(
                    std::string(para) + std::to_string(value)));
                if (value > 1 || value < 0)
                { REQUIRE_FALSE(commander.ref_check_wrapper()); }
                else
                {
                    REQUIRE(commander.ref_check_wrapper());
                }
            }
            SECTION("bgen bound check")
            {
                auto para = GENERATE("--ld-hard-thres ", "--ld-dose-thres ",
                                     "--ld-info ");
                auto value = GENERATE(take(100, random(-1.1, 1.1)));
                REQUIRE(commander.parse_command_wrapper(
                    "--ld-type bgen " + std::string(para)
                    + std::to_string(value)));
                if (value > 1 || value < 0)
                { REQUIRE_FALSE(commander.ref_check_wrapper()); }
                else
                {
                    REQUIRE(commander.ref_check_wrapper());
                }
            }
            SECTION("invalid type")
            {
                REQUIRE(commander.parse_command_wrapper("--ld-type vcf"));
                REQUIRE_FALSE(commander.ref_check_wrapper());
            }
        }
    }
    SECTION("Misc")
    {
        SECTION("Negative thread")
        {
            REQUIRE(commander.parse_command_wrapper("--thread -10"));
            REQUIRE_FALSE(commander.misc_check_wrapper());
        }
        SECTION("logit-perm")
        {
            REQUIRE(commander.parse_command_wrapper("--logit-perm"));
            // doesn't matter if we don't use perm, only add warning
            REQUIRE(commander.misc_check_wrapper());
        }
        SECTION("use-ref-maf")
        {
            REQUIRE(commander.parse_command_wrapper("--use-ref-maf"));
            SECTION("without ld")
            {
                REQUIRE_FALSE(commander.misc_check_wrapper());
            }
            SECTION("with ld")
            {
                REQUIRE(commander.parse_command_wrapper("--ld ref"));
                REQUIRE(commander.misc_check_wrapper());
            }
        }
        SECTION("Intermediate check")
        {
            REQUIRE(commander.parse_command_wrapper("--allow-inter"));
            SECTION("no bgen")
            {
                SECTION("no ld")
                {
                    REQUIRE(commander.misc_check_wrapper());
                    REQUIRE_FALSE(commander.use_inter());
                }
                SECTION("bed ld")
                {
                    REQUIRE(commander.parse_command_wrapper("--ld ref"));
                    REQUIRE(commander.misc_check_wrapper());
                    REQUIRE_FALSE(commander.use_inter());
                }
                SECTION("bgen ld")
                {
                    REQUIRE(commander.parse_command_wrapper(
                        "--ld-type bgen --ld ref"));
                    REQUIRE(commander.misc_check_wrapper());
                    REQUIRE(commander.use_inter());
                }
            }
            SECTION("bgen")
            {
                REQUIRE(commander.parse_command_wrapper("--type bgen"));
                SECTION("with bed reference")
                {
                    REQUIRE(commander.parse_command_wrapper(
                        "--ld-type bed --ld ref"));
                    SECTION("dosage score")
                    {
                        REQUIRE(commander.misc_check_wrapper());
                        REQUIRE_FALSE(commander.use_inter());
                    }
                    SECTION("hardcoded score")
                    {
                        REQUIRE(commander.parse_command_wrapper("--hard"));
                        REQUIRE(commander.misc_check_wrapper());
                        REQUIRE(commander.use_inter());
                    }
                }
                SECTION("only target")
                {
                    REQUIRE(commander.misc_check_wrapper());
                    REQUIRE(commander.use_inter());
                }
            }
        }
        SECTION("ultra")
        {
            REQUIRE(commander.parse_command_wrapper("--ultra"));
            SECTION("with bed")
            {
                REQUIRE(commander.parse_command_wrapper("--type bed"));
                REQUIRE(commander.misc_check_wrapper());
                REQUIRE(commander.ultra_aggressive());
            }
            SECTION("with bgen")
            {
                REQUIRE(commander.parse_command_wrapper("--type bgen"));
                REQUIRE(commander.misc_check_wrapper());
                REQUIRE_FALSE(commander.ultra_aggressive());
            }
        }
    }
    SECTION("snp selection")
    {
        SECTION("extract")
        {
            REQUIRE(commander.parse_command_wrapper("--extract keep"));
            REQUIRE(commander.misc_check_wrapper());
        }
        SECTION("exclude")
        {
            REQUIRE(commander.parse_command_wrapper("--exclude remove"));
            REQUIRE(commander.misc_check_wrapper());
        }
        SECTION("both")
        {
            REQUIRE(commander.parse_command_wrapper(
                "--extract keep --exclude remove"));
            REQUIRE_FALSE(commander.misc_check_wrapper());
        }
    }
    SECTION("PRSet")
    {
        SECTION("perm check")
        {
            // not allow both
            REQUIRE(commander.parse_command_wrapper(
                "--perm 100 --set-perm 100 --snp-set snps"));
            REQUIRE_FALSE(commander.prset_check_wrapper());
        }
        SECTION("bed")
        {
            REQUIRE(
                commander.parse_command_wrapper("--bed test --set-perm 10"));
            // only need background when we perform set based permutation
            SECTION("no gtf use full")
            {
                REQUIRE(commander.prset_check_wrapper());
                REQUIRE(commander.get_set().full_as_background);
            }
            SECTION("with gtf")
            {
                REQUIRE(commander.parse_command_wrapper("--gtf Homo"));
                REQUIRE(commander.prset_check_wrapper());
                REQUIRE_FALSE(commander.get_set().full_as_background);
            }
        }
        SECTION("msigdb")
        {
            REQUIRE(commander.parse_command_wrapper("--msigdb kegg"));
            SECTION("no gtf")
            {
                REQUIRE_FALSE(commander.prset_check_wrapper());
            }
            SECTION("with gtf")
            {
                REQUIRE(commander.parse_command_wrapper("--gtf Homo"));
                REQUIRE(commander.prset_check_wrapper());
            }
            SECTION("feature")
            {
                REQUIRE(commander.parse_command_wrapper("--gtf Homo"));
                SECTION("default")
                {
                    REQUIRE(commander.prset_check_wrapper());
                    REQUIRE_THAT(
                        commander.get_set().feature,
                        Catch::UnorderedEquals<std::string>(
                            {"exon", "gene", "protein_coding", "CDS"}));
                }
                SECTION("not default")
                {
                    REQUIRE(commander.parse_command_wrapper(
                        "--feature gene,intron"));
                    REQUIRE(commander.prset_check_wrapper());
                    REQUIRE_THAT(commander.get_set().feature,
                                 Catch::UnorderedEquals<std::string>(
                                     {"gene", "intron"}));
                }
            }
        }
    }
    SECTION("Covariate")
    {
        REQUIRE(commander.parse_command_wrapper("--cov test"));
        std::unordered_map<std::string, size_t> ref_index;
        SECTION("With valid cov-col")
        {
            REQUIRE(commander.parse_command_wrapper(
                "--cov-col Sex,@PC[1-5.6],@PC[6.8.9],@Hi"));
            auto included = commander.get_cov_names_wrap();
            auto ori_input = included;
            SECTION("Check include content")
            {
                std::vector<std::string> res;
                res.insert(res.end(), included.begin(), included.end());
                REQUIRE_THAT(res, Catch::UnorderedEquals<std::string>(
                                      {"Sex", "PC1", "PC2", "PC3", "PC4", "PC5",
                                       "PC6", "PC8", "PC9", "Hi"}));
            }
            SECTION("Mimic file read")
            {
                std::vector<std::string> cov_name = {"FID", "IID", "PC3", "PC2",
                                                     "PC1"};
                commander.prepare_header_cov_check_wrapper(cov_name, ref_index,
                                                           included);
                std::string missing = "";
                size_t valid_cov =
                    commander.find_cov_idx_wrap(included, ref_index, missing);
                REQUIRE_FALSE(valid_cov == 0);
                auto cov = commander.get_pheno().cov_colname;
                REQUIRE_THAT(cov,
                             Catch::UnorderedEquals<std::string>(
                                 {"Sex", "@PC[1-5.6]", "@PC[6.8.9]", "@Hi"}));
                auto cov_idx = commander.get_pheno().col_index_of_cov;
                REQUIRE_THAT(cov_idx,
                             Catch::UnorderedEquals<size_t>({4, 3, 2}));
                // now update cov_name
                SECTION("update cov_name")
                {
                    commander.reorganize_cov_name_wrap(cov_name);
                    cov = commander.get_pheno().cov_colname;
                    auto cov_idx = commander.get_pheno().col_index_of_cov;
                    REQUIRE_THAT(
                        cov, Catch::Equals<std::string>({"PC3", "PC2", "PC1"}));
                    REQUIRE_THAT(cov_idx, Catch::Equals<size_t>({2, 3, 4}));
                    SECTION("with invalid factor (not found)")
                    {
                        REQUIRE(commander.parse_command_wrapper(
                            "--cov-factor @PC[10-20]"));
                        REQUIRE_FALSE(commander.process_factor_cov_wrap(
                            included, ref_index, ori_input));
                    }
                }
            }
        }
        SECTION("with cov-col but not in file")
        {
            REQUIRE(commander.parse_command_wrapper("--cov-col Sex,Age"));
            std::unordered_set<std::string> included =
                commander.get_cov_names_wrap();
            std::vector<std::string> cov_name = {"FID", "IID", "PC3", "PC2",
                                                 "PC1"};
            commander.prepare_header_cov_check_wrapper(cov_name, ref_index,
                                                       included);
            std::string missing = "";
            size_t valid_cov =
                commander.find_cov_idx_wrap(included, ref_index, missing);
            REQUIRE(valid_cov == 0);
        }
        SECTION("with FID as cov")
        {
            REQUIRE(commander.parse_command_wrapper("--cov-col FID"));
            std::unordered_set<std::string> included =
                commander.get_cov_names_wrap();
            std::vector<std::string> cov_name = {"FID", "IID", "PC3", "PC2",
                                                 "PC1"};
            commander.prepare_header_cov_check_wrapper(cov_name, ref_index,
                                                       included);
            std::string missing = "";
            size_t valid_cov =
                commander.find_cov_idx_wrap(included, ref_index, missing);
            REQUIRE(valid_cov == 1);
        }
        SECTION("without cov-col")
        {
            std::vector<std::string> cov_name = {"FID", "IID", "PC3", "PC2",
                                                 "PC1"};
            std::unordered_set<std::string> included, ori_input;
            included = commander.get_cov_names_wrap();
            ori_input = included;
            REQUIRE(included.empty());
            commander.prepare_header_cov_check_wrapper(cov_name, ref_index,
                                                       included);
            std::vector<std::string> res;
            res.insert(res.end(), included.begin(), included.end());
            REQUIRE_THAT(res, Catch::UnorderedEquals<std::string>(
                                  {"PC3", "PC2", "PC1"}));

            SECTION("with factor cov")
            {
                REQUIRE(commander.parse_command_wrapper("--cov-factor PC3"));
                CHECK(commander.process_factor_cov_wrap(included, ref_index,
                                                        ori_input));
                REQUIRE_THAT(commander.get_pheno().col_index_of_factor_cov,
                             Catch::Equals<size_t>({2}));
            }
        }
    }
}
