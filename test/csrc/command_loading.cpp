#include "catch.hpp"
#include "commander.hpp"
#include "mock_commander.hpp"


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

TEST_CASE("Commander Initialization")
{
    Commander base_commander;
    SECTION("default")
    {
        REQUIRE_FALSE(base_commander.all_scores());
        REQUIRE_FALSE(base_commander.use_inter());
        REQUIRE(base_commander.delim() == " ");
        REQUIRE_THAT(base_commander.out(), Catch::Matchers::Equals("PRSice"));
        REQUIRE(base_commander.exclusion_range().empty());
        REQUIRE(base_commander.exclude_file().empty());
        REQUIRE(base_commander.extract_file().empty());
        REQUIRE(base_commander.max_memory(1.0) == Approx(1.0));
        REQUIRE(base_commander.max_memory(2.0) == Approx(2.0));
    }
}

TEST_CASE("Commander set flag")
{
    mockCommander commander;
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

TEST_CASE("usage and version")
{
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
}

TEST_CASE("binary target")
{
    mockCommander commander;
    SECTION("Invalid")
    {

        SECTION("1 as T")
        {
            REQUIRE(commander.get_pheno().binary.empty());
            REQUIRE_FALSE(commander.parse_command_wrapper("--binary-target 1"));
        }
        SECTION("0 as F")
        {
            REQUIRE(commander.get_pheno().binary.empty());
            REQUIRE_FALSE(commander.parse_command_wrapper("--binary-target 0"));
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
        REQUIRE(commander.parse_command_wrapper("--binary-target F,2T, 3F"));
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
            REQUIRE(commander.parse_command_wrapper("--binary-target True"));
            REQUIRE_THAT(commander.get_pheno().binary,
                         Catch::Equals<bool>({true}));
        }
        SECTION("true")
        {
            REQUIRE(commander.get_pheno().binary.empty());
            // 3F will be ignored
            REQUIRE(commander.parse_command_wrapper("--binary-target true"));
            REQUIRE_THAT(commander.get_pheno().binary,
                         Catch::Equals<bool>({true}));
        }
        SECTION("1true")
        {
            REQUIRE(commander.get_pheno().binary.empty());
            // 3F will be ignored
            REQUIRE(commander.parse_command_wrapper("--binary-target 1true"));
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
            REQUIRE(commander.parse_command_wrapper("--binary-target False"));
            REQUIRE_THAT(commander.get_pheno().binary,
                         Catch::Equals<bool>({false}));
        }
        SECTION("1false")
        {
            REQUIRE(commander.get_pheno().binary.empty());
            // 3F will be ignored
            REQUIRE(commander.parse_command_wrapper("--binary-target 1false"));
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
            REQUIRE(commander.parse_command_wrapper("--binary-target True,3F"));
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

TEST_CASE("bar level loading")
{
    mockCommander commander;
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
        REQUIRE_THAT(res, Catch::Equals<double>({0.1, -0.2, 0.3, 0.4, 0.5}));
    }

    SECTION("Contain Zero")
    {
        // supposed to fail but init didn't check for these either
        REQUIRE(commander.parse_command_wrapper(
            "--bar-levels 0,0.2,0.3,0.3,0.4,0.5"));
        auto res = commander.get_p_threshold().bar_levels;
        REQUIRE_THAT(res, Catch::Equals<double>({0, 0.2, 0.3, 0.3, 0.4, 0.5}));
    }
    SECTION("Out of p-value range")
    {
        // supposed to fail but init didn't check for these either
        REQUIRE(commander.parse_command_wrapper(
            "--bar-levels 0.5,0.2,0.3,0.3,0.4,0.1,2"));
        auto res = commander.get_p_threshold().bar_levels;
        REQUIRE_THAT(res,
                     Catch::Equals<double>({0.5, 0.2, 0.3, 0.3, 0.4, 0.1, 2}));
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

TEST_CASE("Base parameter loading")
{
    mockCommander commander;
    SECTION("A1")
    {
        check_set_base_flag(commander, "--A1", "a", "A1", +BASE_INDEX::EFFECT);
    }
    SECTION("a1")
    {
        check_set_base_flag(commander, "--a1", "b", "A1", +BASE_INDEX::EFFECT);
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
        check_set_base_flag(commander, "-p", "postdoc", "P", +BASE_INDEX::P);
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
        check_set_base_flag(commander, "--snp", "cnv", "SNP", +BASE_INDEX::RS);
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
        check_set_base_flag(commander, "--base-info", "INFO_FILTER", "INFO:0.9",
                            +BASE_INDEX::INFO);
    }
    SECTION("base-maf")
    {
        check_set_base_flag(commander, "--base-maf", "MAF_FILTER", "",
                            +BASE_INDEX::MAF);
    }
}

TEST_CASE("Target parameter loading")
{
    mockCommander commander;
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
        SECTION("name") { REQUIRE(commander.get_target().file_name.empty()); }
        SECTION("name-list")
        {
            REQUIRE(commander.get_target().file_list.empty());
        }
        SECTION("keep") { REQUIRE(commander.get_target().keep.empty()); }
        SECTION("remove") { REQUIRE(commander.get_target().remove.empty()); }
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
            REQUIRE(commander.parse_command_wrapper("--target-list testing"));
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
            REQUIRE(commander.get_target_qc().hard_threshold == Approx(0.1));
            REQUIRE(commander.parse_command_wrapper("--dose-thres 1.0"));
            REQUIRE(commander.get_target_qc().dose_threshold == Approx(1.0));
            REQUIRE(commander.get_target_qc().hard_threshold == Approx(0.1));
        }

        SECTION("hard-thres")
        {
            REQUIRE(commander.get_target_qc().dose_threshold == 0.0);
            REQUIRE(commander.get_target_qc().hard_threshold == Approx(0.1));
            REQUIRE(commander.parse_command_wrapper("--hard-thres -0.1"));
            REQUIRE(commander.get_target_qc().dose_threshold == 0.0);
            REQUIRE(commander.get_target_qc().hard_threshold == Approx(-0.1));
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

TEST_CASE("Phenotype parameter loading")
{
    mockCommander commander;
    SECTION("default")
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
            REQUIRE_THAT(commander.get_pheno().pheno_col,
                         Catch::Equals<std::string>({"A1", "B1", "C1", "D2"}));
        }
        SECTION("prevalence")
        {
            REQUIRE(
                commander.parse_command_wrapper("--prevalence 0.1,0.4,1,-0.5"));
            REQUIRE_THAT(commander.get_pheno().prevalence,
                         Catch::Equals<double>({0.1, 0.4, 1, -0.5}));
        }
    }
}

TEST_CASE("Clumping parameter loading")
{
    mockCommander commander;
    SECTION("default")
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
            REQUIRE_FALSE(commander.parse_command_wrapper("--clump-kb -100kb"));
        }
    }
}

TEST_CASE("PRSet parameter loading")
{
    mockCommander commander;
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

TEST_CASE("misc parameter loading")
{
    mockCommander commander;
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
            REQUIRE_FALSE(commander.parse_command_wrapper("--num-auto 1e100"));
        }
    }
}

TEST_CASE("Permutation parameter loading")
{
    mockCommander commander;
    SECTION("default")
    {
        REQUIRE(commander.get_perm().num_permutation == 0);
        REQUIRE_FALSE(commander.get_perm().run_perm);
        REQUIRE_FALSE(commander.get_perm().run_set_perm);
    }
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

TEST_CASE("P-value thresholding parameter loading")
{
    mockCommander commander;
    SECTION("default")
    {
        REQUIRE(commander.get_p_threshold().inter == Approx(0.00005));
        REQUIRE(commander.get_p_threshold().lower == Approx(5e-8));
        REQUIRE(commander.get_p_threshold().upper == Approx(0.5));
    }
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

TEST_CASE("Scoring parameter loading")
{
    mockCommander commander;
    SECTION("default")
    {
        REQUIRE(commander.get_prs_instruction().scoring_method
                == SCORING::AVERAGE);
    }
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

TEST_CASE("missing parameter loading")
{
    mockCommander commander;
    SECTION("default")
    {
        REQUIRE(commander.get_prs_instruction().missing_score
                == MISSING_SCORE::MEAN_IMPUTE);
    }
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
        REQUIRE(commander.parse_command_wrapper("--missing IMPUTE_CONTROL"));
        REQUIRE(commander.get_prs_instruction().missing_score
                == MISSING_SCORE::IMPUTE_CONTROL);
    }
}

TEST_CASE("model parameter loading")
{
    mockCommander commander;
    SECTION("default")
    {
        REQUIRE(commander.get_prs_instruction().genetic_model
                == MODEL::ADDITIVE);
    }
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

TEST_CASE("Covariate parameter loading")
{
    mockCommander commander;
    SECTION("default")
    {
        REQUIRE(commander.get_pheno().cov_file.empty());
        REQUIRE(commander.get_pheno().cov_colname.empty());
        REQUIRE(commander.get_pheno().factor_cov.empty());
    }
    SECTION("simple load")
    {
        REQUIRE(commander.parse_command_wrapper("--cov-col Testing,@PC[1-55]"));
        REQUIRE_THAT(commander.get_pheno().cov_colname,
                     Catch::Equals<std::string>({"Testing", "@PC[1-55]"}));
        SECTION("append")
        {
            REQUIRE(
                commander.parse_command_wrapper("--cov-col More,Covariate"));
            REQUIRE_THAT(commander.get_pheno().cov_colname,
                         Catch::Equals<std::string>(
                             {"Testing", "@PC[1-55]", "More", "Covariate"}));
        }
    }
    SECTION("discontinue load")
    {
        REQUIRE(
            commander.parse_command_wrapper("--cov-col Testing,@PC[1.3.5]"));
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

TEST_CASE("Reference parameter loading")
{
    mockCommander commander;
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
        SECTION("remove") { REQUIRE(commander.get_reference().remove.empty()); }
        SECTION("type") { REQUIRE(commander.get_reference().type == "bed"); }
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
            REQUIRE(commander.parse_command_wrapper("--ld-hard-thres -0.1"));
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
