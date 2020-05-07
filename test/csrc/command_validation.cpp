#include "catch.hpp"
#include "commander.hpp"
#include "mock_commander.hpp"

TEST_CASE("Target validation")
{
    mockCommander commander;
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
            REQUIRE(
                commander.parse_command_wrapper("--type " + std::string(i)));
            REQUIRE(commander.target_check_wrapper());
        }
        SECTION("universial bound check")
        {
            auto para = GENERATE("--geno ", "--maf ");
            auto value = GENERATE(take(10, random(-1.1, 1.1)));
            REQUIRE(commander.parse_command_wrapper(std::string(para)
                                                    + std::to_string(value)));
            if (value > 1 || value < 0)
            { REQUIRE_FALSE(commander.filter_check_wrapper()); }
            else
            {
                REQUIRE(commander.filter_check_wrapper());
            }
        }
        SECTION("bgen bound check")
        {
            auto para = GENERATE("--hard-thres ", "--dose-thres ", "--info ");
            auto value = GENERATE(take(10, random(-1.1, 1.1)));
            REQUIRE(commander.parse_command_wrapper(
                "--type bgen " + std::string(para) + std::to_string(value)));
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
            auto i = GENERATE(take(10, random(-100, 100)));
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

TEST_CASE("Thresholding validation")
{
    mockCommander commander;
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
            REQUIRE_THAT(
                commander.get_p_threshold().bar_levels,
                Catch::Approx<double>({0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5}));
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
        auto i = GENERATE(take(10, random(-1.1, 1.1)));
        REQUIRE(commander.parse_command_wrapper("--bar-levels "
                                                + std::to_string(i)));
        if (i <= 1.0 && i >= 0.0) { REQUIRE(commander.prsice_check_wrapper()); }
        else
        {
            REQUIRE_FALSE(commander.prsice_check_wrapper());
        }
    }
    SECTION("high res check")
    {
        auto para = GENERATE("--lower ", "--upper ", "--inter ");
        auto value = GENERATE(take(10, random(-1.1, 1.1)));
        REQUIRE(commander.parse_command_wrapper(std::string(para)
                                                + std::to_string(value)));
        if (value > 1 || value < 0
            || commander.get_p_threshold().lower
                   > commander.get_p_threshold().upper
            || misc::logically_equal(commander.get_p_threshold().inter, 0.0))
        { REQUIRE_FALSE(commander.prsice_check_wrapper()); }
        else
        {
            REQUIRE(commander.prsice_check_wrapper());
        }
    }
}

TEST_CASE("Clump")
{
    mockCommander commander;
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
        auto value = GENERATE(take(10, random(-1.1, 1.1)));
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

TEST_CASE("Test automatic reference copy over")
{
    mockCommander commander;
    auto bgen = GENERATE("--type bed", "--type bgen");
    REQUIRE(
        commander.parse_command_wrapper("--target name --maf 0.05 --geno 0.01 "
                                        "--hard-thres 0.2 --dose-thres 0.9 "
                                        "--keep A --remove B --info "
                                        "0.8 --info-type mach "
                                        + std::string(bgen)));
    // command, copy over
    SECTION("Test need target as reference function")
    {
        using record = std::tuple<std::string, bool, bool>;
        auto input = GENERATE(table<std::string, bool, bool>(
            {record {"--ld-maf 0.2", true, false},
             record {"--ld-maf 0.05", false, false},
             record {"--ld-geno 0.2", true, false},
             record {"--ld-geno 0.01", false, false},
             record {"--ld-hard-thres 0.3", true, true},
             record {"--ld-hard-thres 0.2", false, true},
             record {"--ld-dose-thres 0.8", true, true},
             record {"--ld-dose-thres 0.9", false, true},
             record {"--ld-info 0.2", true, true},
             record {"--ld-info 0.8", false, true},
             record {"--ld-keep C", true, false},
             record {"--ld-keep A", false, false},
             record {"--ld-remove C", true, false},
             record {"--ld-remove B", false, false}}));
        REQUIRE(commander.parse_command_wrapper(std::get<0>(input)));
        if (std::get<2>(input))
        {
            // bgen only
            if (std::string(bgen) == "--type bgen")
            {
                REQUIRE(commander.test_need_target_as_reference()
                        == std::get<1>(input));
            }
            else
            {
                // don't need the bgen related parameter, so won't copy over
                REQUIRE_FALSE(commander.test_need_target_as_reference());
            }
        }
        else
        {
            REQUIRE(commander.test_need_target_as_reference()
                    == std::get<1>(input));
        }
        SECTION("Combine with reference validation")
        {
            REQUIRE(commander.ref_check_wrapper());
            if (std::get<2>(input))
            {
                // bgen only
                if (std::string(bgen) == "--type bgen")
                {
                    REQUIRE(commander.has_parameter("ld")
                            == std::get<1>(input));
                }
                else
                {
                    // don't need the bgen related parameter, so won't copy over
                    REQUIRE_FALSE(commander.has_parameter("ld"));
                }
            }
            else
            {
                REQUIRE(commander.has_parameter("ld") == std::get<1>(input));
            }
        }
    }
}

TEST_CASE("Reference validation")
{
    mockCommander commander;
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
        REQUIRE(commander.parse_command_wrapper("--ld test --ld-list list"));
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
            REQUIRE(
                commander.parse_command_wrapper("--ld-type " + std::string(i)));
            REQUIRE(commander.ref_check_wrapper());
        }
        SECTION("universial bound check")
        {
            auto para = GENERATE("--ld-geno ", "--ld-maf ");
            auto value = GENERATE(take(10, random(-1.1, 1.1)));
            REQUIRE(commander.parse_command_wrapper(std::string(para)
                                                    + std::to_string(value)));
            if (value > 1 || value < 0)
            { REQUIRE_FALSE(commander.ref_check_wrapper()); }
            else
            {
                REQUIRE(commander.ref_check_wrapper());
            }
        }
        SECTION("bgen bound check")
        {
            auto para =
                GENERATE("--ld-hard-thres ", "--ld-dose-thres ", "--ld-info ");
            auto value = GENERATE(take(10, random(-1.1, 1.1)));
            REQUIRE(commander.parse_command_wrapper(
                "--ld-type bgen " + std::string(para) + std::to_string(value)));
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

TEST_CASE("Misc validation")
{
    mockCommander commander;
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
        SECTION("without ld") { REQUIRE_FALSE(commander.misc_check_wrapper()); }
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
                REQUIRE(
                    commander.parse_command_wrapper("--ld-type bgen --ld ref"));
                REQUIRE(commander.misc_check_wrapper());
                REQUIRE(commander.use_inter());
            }
        }
        SECTION("bgen")
        {
            REQUIRE(commander.parse_command_wrapper("--type bgen"));
            SECTION("with bed reference")
            {
                REQUIRE(
                    commander.parse_command_wrapper("--ld-type bed --ld ref"));
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
}

TEST_CASE("PRSet validation")
{
    mockCommander commander;
    SECTION("perm check")
    {
        // not allow both
        REQUIRE(commander.parse_command_wrapper(
            "--perm 100 --set-perm 100 --snp-set snps"));
        REQUIRE_FALSE(commander.prset_check_wrapper());
    }
    SECTION("bed")
    {
        REQUIRE(commander.parse_command_wrapper("--bed test --set-perm 10"));
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
        SECTION("no gtf") { REQUIRE_FALSE(commander.prset_check_wrapper()); }
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
                REQUIRE_THAT(commander.get_set().feature,
                             Catch::UnorderedEquals<std::string>(
                                 {"exon", "gene", "protein_coding", "CDS"}));
            }
            SECTION("not default")
            {
                REQUIRE(
                    commander.parse_command_wrapper("--feature gene,intron"));
                REQUIRE(commander.prset_check_wrapper());
                REQUIRE_THAT(
                    commander.get_set().feature,
                    Catch::UnorderedEquals<std::string>({"gene", "intron"}));
            }
        }
    }
}

TEST_CASE("Covariate validation")
{
    mockCommander commander;
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
            REQUIRE_THAT(cov, Catch::UnorderedEquals<std::string>(
                                  {"Sex", "@PC[1-5.6]", "@PC[6.8.9]", "@Hi"}));
            auto cov_idx = commander.get_pheno().col_index_of_cov;
            REQUIRE_THAT(cov_idx, Catch::UnorderedEquals<size_t>({4, 3, 2}));
            // now update cov_name
            SECTION("update cov_name")
            {
                commander.reorganize_cov_name_wrap(cov_name);
                cov = commander.get_pheno().cov_colname;
                auto cov_idx = commander.get_pheno().col_index_of_cov;
                REQUIRE_THAT(cov,
                             Catch::Equals<std::string>({"PC3", "PC2", "PC1"}));
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
        std::vector<std::string> cov_name = {"FID", "IID", "PC3", "PC2", "PC1"};
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
        std::vector<std::string> cov_name = {"FID", "IID", "PC3", "PC2", "PC1"};
        commander.prepare_header_cov_check_wrapper(cov_name, ref_index,
                                                   included);
        std::string missing = "";
        size_t valid_cov =
            commander.find_cov_idx_wrap(included, ref_index, missing);
        REQUIRE(valid_cov == 1);
    }
    SECTION("without cov-col")
    {
        std::vector<std::string> cov_name = {"FID", "IID", "PC3", "PC2", "PC1"};
        std::unordered_set<std::string> included, ori_input;
        included = commander.get_cov_names_wrap();
        ori_input = included;
        REQUIRE(included.empty());
        commander.prepare_header_cov_check_wrapper(cov_name, ref_index,
                                                   included);
        std::vector<std::string> res;
        res.insert(res.end(), included.begin(), included.end());
        REQUIRE_THAT(
            res, Catch::UnorderedEquals<std::string>({"PC3", "PC2", "PC1"}));

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

TEST_CASE("Phenotype validation")
{
    mockCommander commander;
    SECTION("col without pheno")
    {
        REQUIRE(commander.parse_command_wrapper("--pheno-col A,B,C"));
        REQUIRE_FALSE(commander.pheno_check_wrapper());
    }
    REQUIRE(commander.parse_command_wrapper("--pheno Pheno"));
    SECTION("Check default binary-target parsing")
    {
        SECTION("with beta flag")
        {
            REQUIRE(commander.parse_command_wrapper("--beta"));
            SECTION("with single col")
            {
                REQUIRE(commander.pheno_check_wrapper());
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({false}));
            }
            SECTION("with multiple col")
            {
                REQUIRE(commander.parse_command_wrapper("--pheno-col A,B,C,D"));
                REQUIRE(commander.pheno_check_wrapper());
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>(std::vector<bool>(4, false)));
            }
        }
        SECTION("with or flag")
        {
            REQUIRE(commander.parse_command_wrapper("--or"));
            SECTION("with single col")
            {
                REQUIRE(commander.pheno_check_wrapper());
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>({true}));
            }
            SECTION("with multiple col")
            {
                REQUIRE(commander.parse_command_wrapper("--pheno-col A,B,C"));
                REQUIRE(commander.pheno_check_wrapper());
                REQUIRE_THAT(commander.get_pheno().binary,
                             Catch::Equals<bool>(std::vector<bool>(3, true)));
            }
        }
    }

    SECTION("Duplicated phenotyp ename")
    {
        REQUIRE(commander.parse_command_wrapper("--pheno-col A,B,A"));
        REQUIRE_FALSE(commander.pheno_check_wrapper());
    }
    SECTION("matching col,  binary target and prevalence")
    {
        REQUIRE(commander.parse_command_wrapper("--pheno-col A,B,C"));
        SECTION("mismatch length")
        {
            SECTION("simple binary target")
            {
                REQUIRE(commander.parse_command_wrapper("--binary-target T"));
                REQUIRE_FALSE(commander.pheno_check_wrapper());
            }
            SECTION("complex binary target")
            {
                REQUIRE(
                    commander.parse_command_wrapper("--binary-target 3T,4F"));
                REQUIRE_FALSE(commander.pheno_check_wrapper());
            }
        }
        SECTION("matching")
        {
            REQUIRE(commander.parse_command_wrapper("--binary-target 3T"));
            REQUIRE(commander.pheno_check_wrapper());
        }
        SECTION("match prevalence")
        {
            SECTION("no binary")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--binary-target 3F --prevalence 0.05"));
                REQUIRE_FALSE(commander.pheno_check_wrapper());
            }
            SECTION("match length")
            {
                REQUIRE(
                    commander.parse_command_wrapper("--binary-target F,2T -k "
                                                    "0.1,0.3"));
                REQUIRE(commander.pheno_check_wrapper());
            }
            SECTION("mismatch length")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--binary-target F,2T,F,3T -k "
                    "0.1"));
                REQUIRE_FALSE(commander.pheno_check_wrapper());
            }
            SECTION("prevalence bound")
            {
                auto i = GENERATE(take(10, random(1.1, -1.0)));
                REQUIRE(commander.parse_command_wrapper(
                    "--binary-target F,2T -k 0.1," + std::to_string(i)));
                if (i <= 0 || i >= 1)
                { REQUIRE_FALSE(commander.pheno_check_wrapper()); }
                else
                {
                    REQUIRE(commander.pheno_check_wrapper());
                }
            }
        }
    }
}

TEST_CASE("Base validaation")
{
    mockCommander commander;
    SECTION("base file not provided")
    {
        REQUIRE_FALSE(commander.base_check_wrapper());
    }
    SECTION("get base name")
    {
        REQUIRE(commander.parse_command_wrapper(
            "--base /home/bin/Base.gz.summary"));
        REQUIRE(commander.get_base_name() == "Base.gz");
    }
    REQUIRE(commander.parse_command_wrapper("--base Base"));
    SECTION("check auto loading works")
    {
        std::vector<std::string> base_file_header = {"P",  "BETA", "CHR", "BP",
                                                     "A1", "A2",   "SNP"};
        std::vector<size_t> expected = {2, 5, 3, 0, 0, 0, 0, 4, 6, 0, 1, 6};
        std::vector<int> expected_has = {true,  true,  true,  false,
                                         false, false, false, true,
                                         true,  true,  true,  false};
        SECTION("both beta and or")
        {
            REQUIRE(commander.parse_command_wrapper("--beta --or"));
            REQUIRE_FALSE(
                commander.base_column_check_wrapper(base_file_header));
        }
        SECTION("check default")
        {
            REQUIRE(commander.base_column_check_wrapper(base_file_header));
            REQUIRE_THAT(commander.get_base().column_index,
                         Catch::Equals<size_t>(expected));
            REQUIRE_THAT(commander.get_base().has_column,
                         Catch::Equals<int>(expected_has));
        }
        SECTION("not found warning")
        {
            using record = std::tuple<std::string, size_t>;
            auto gen_record = GENERATE(table<std::string, size_t>(
                {record {"--chr", +BASE_INDEX::CHR},
                 record {"--bp", +BASE_INDEX::BP},
                 record {"--a2", +BASE_INDEX::NONEFFECT}}));
            auto cmd = std::get<0>(gen_record);
            auto idx = std::get<1>(gen_record);
            expected[idx] = 0;
            expected_has[idx] = false;
            REQUIRE(
                commander.parse_command_wrapper(std::string(cmd) + " XXXX"));
            REQUIRE(commander.base_column_check_wrapper(base_file_header));
            REQUIRE_THAT(commander.get_base().column_index,
                         Catch::Equals<size_t>(expected));
            REQUIRE_THAT(commander.get_base().has_column,
                         Catch::Equals<int>(expected_has));
        }
        SECTION("not found fatal")
        {
            using record = std::tuple<std::string, size_t>;
            auto gen_record = GENERATE(table<std::string, size_t>(
                {record {"--stat", +BASE_INDEX::STAT},
                 record {"--snp", +BASE_INDEX::RS},
                 record {"--a1", +BASE_INDEX::EFFECT},
                 record {"--pvalue", +BASE_INDEX::P}}));
            auto cmd = std::get<0>(gen_record);
            auto idx = std::get<1>(gen_record);
            expected[idx] = 0;
            expected_has[idx] = false;
            REQUIRE(
                commander.parse_command_wrapper(std::string(cmd) + " XXXX"));
            REQUIRE_FALSE(
                commander.base_column_check_wrapper(base_file_header));
        }
        SECTION("no default")
        {
            REQUIRE(commander.parse_command_wrapper("--no-default"));
            SECTION("without vital")
            {
                REQUIRE_FALSE(
                    commander.base_column_check_wrapper(base_file_header));
            }
            SECTION("with vital")
            {
                std::vector<size_t> expected = {0, 0, 0, 0, 0, 0,
                                                0, 4, 6, 0, 1, 6};
                std::vector<int> expected_has = {false, false, false, false,
                                                 false, false, false, true,
                                                 true,  true,  true,  false};
                REQUIRE(commander.parse_command_wrapper(
                    "--snp SNP --pvalue P --stat BETA --a1 A1"));
                REQUIRE(commander.base_column_check_wrapper(base_file_header));
                REQUIRE_THAT(commander.get_base().column_index,
                             Catch::Equals<size_t>(expected));
                REQUIRE_THAT(commander.get_base().has_column,
                             Catch::Equals<int>(expected_has));
            }
        }
        SECTION("index")
        {
            REQUIRE(commander.parse_command_wrapper("--index"));
            SECTION("without vital")
            {
                // similar to no-default
                REQUIRE_FALSE(
                    commander.base_column_check_wrapper(base_file_header));
            }
            SECTION("with colname input")
            {
                // need numeric input
                REQUIRE(commander.parse_command_wrapper(
                    "--snp SNP --pvalue P --stat BETA --a1 A1"));
                REQUIRE_FALSE(
                    commander.base_column_check_wrapper(base_file_header));
            }
            SECTION("with negative index")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--snp 6 --pvalue 0 --stat -1 --a1 4"));
                REQUIRE_FALSE(
                    commander.base_column_check_wrapper(base_file_header));
            }
            SECTION("with index input")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--snp 6 --pvalue 0 --stat 1 --a1 4"));
                std::vector<size_t> expected = {0, 0, 0, 0, 0, 0,
                                                0, 4, 6, 0, 1, 6};
                std::vector<int> expected_has = {false, false, false, false,
                                                 false, false, false, true,
                                                 true,  true,  true,  false};
                SECTION("with beta / or")
                {
                    auto i = GENERATE("--beta", "--or");
                    REQUIRE(commander.parse_command_wrapper(i));
                    REQUIRE(
                        commander.base_column_check_wrapper(base_file_header));
                    REQUIRE_THAT(commander.get_base().column_index,
                                 Catch::Equals<size_t>(expected));
                    REQUIRE_THAT(commander.get_base().has_column,
                                 Catch::Equals<int>(expected_has));
                }
                SECTION("without beta / or")
                {
                    // can't guess without header
                    REQUIRE_FALSE(
                        commander.base_column_check_wrapper(base_file_header));
                }
            }
        }
    }
    SECTION("check non-default input")
    {
        std::vector<std::string> base_file_header = {
            "pvalue", "z-score", "coordinate", "effect", "non-effect",
            "empty",  "pad",     "check",      "rsid",   "chrom"};
        std::vector<size_t> expected = {9, 4, 2, 0, 0, 0, 0, 3, 8, 0, 1, 9};
        std::vector<int> expected_has = {true,  true,  true,  false,
                                         false, false, false, true,
                                         true,  true,  true,  false};
        REQUIRE(commander.parse_command_wrapper(
            "--base Base --pvalue pvalue --a1 effect --a2 "
            "non-effect --snp rsid --chr chrom --bp coordinate"));
        SECTION("without beta / or")
        {
            REQUIRE(commander.parse_command_wrapper("--stat z-score"));
            REQUIRE_FALSE(
                commander.base_column_check_wrapper(base_file_header));
        }
        SECTION("with beta / or")
        {
            auto i = GENERATE("--beta", "--or");
            REQUIRE(commander.parse_command_wrapper(i));
            SECTION("stat provided")
            {
                commander.parse_command_wrapper("--stat z-score");
                REQUIRE(commander.base_column_check_wrapper(base_file_header));
                REQUIRE_THAT(commander.get_base().column_index,
                             Catch::Equals<size_t>(expected));
                REQUIRE_THAT(commander.get_base().has_column,
                             Catch::Equals<int>(expected_has));
            }
            SECTION("stat not provided")
            {
                // we look for beta or or in header, and non-were found
                REQUIRE_FALSE(
                    commander.base_column_check_wrapper(base_file_header));
            }
        }
    }
    SECTION("guess beta / or")
    {
        // header, command, beta flag, or flag, expect ok
        using record = std::tuple<std::string, std::string, bool, bool, bool>;
        auto tests = GENERATE(table<std::string, std::string, bool, bool, bool>(
            {record {"or", "", false, true, true},
             record {"Or", "", false, true, true},
             record {"beta", "", true, false, true},
             record {"BeTa", "", true, false, true},
             record {"or", "--stat or", false, true, true},
             record {"Or", "--stat oR", false, true, false}, // need exact match
             record {"Or", "--stat beta", false, true,
                     false}, // need exact match
             record {"beta", "--stat BETA", true, false,
                     false}, // need exact match
             record {"beta", "--stat or", true, false,
                     false}, // need exact match
             record {"BeTa", "--stat BeTa", true, false, true}}));
        std::vector<std::string> base_file_header = {
            "P", std::get<0>(tests), "CHR", "BP", "A1", "A2", "SNP"};
        std::vector<size_t> expected = {2, 5, 3, 0, 0, 0, 0, 4, 6, 0, 1, 6};
        std::vector<int> expected_has = {true,  true,  true,  false,
                                         false, false, false, true,
                                         true,  true,  true,  false};
        if (!std::get<1>(tests).empty())
        { REQUIRE(commander.parse_command_wrapper(std::get<1>(tests))); }
        if (std::get<4>(tests))
        {
            REQUIRE(commander.base_column_check_wrapper(base_file_header));
            REQUIRE(commander.get_base().is_beta == std::get<2>(tests));
            REQUIRE(commander.get_base().is_or == std::get<3>(tests));
        }
        else
        {
            REQUIRE_FALSE(
                commander.base_column_check_wrapper(base_file_header));
        }
    }
    SECTION("file contain both beta and or")
    {
        std::vector<std::string> base_file_header = {"P",  "OR", "CHR", "BP",
                                                     "A1", "A2", "SNP", "BETA"};
        SECTION("without beta or or flag")
        {
            REQUIRE_FALSE(
                commander.base_column_check_wrapper(base_file_header));
        }
        SECTION("with beta or or flag")
        {
            using record = std::tuple<std::string, size_t>;
            auto i = GENERATE(table<std::string, size_t>(
                {record {"--or", 1}, record {"--beta", 7}}));
            REQUIRE(commander.parse_command_wrapper(std::get<0>(i)));
            REQUIRE(commander.base_column_check_wrapper(base_file_header));
            REQUIRE(commander.get_base().column_index[+BASE_INDEX::STAT]
                    == std::get<1>(i));
        }
    }
}

TEST_CASE("Base filtering validation")
{
    mockCommander commander;
    SECTION("info score")
    {
        SECTION("info column has ,")
        {
            std::vector<std::string> col = {"P",   "CHR",  "BETA", "LOC",
                                            "A1",  "A2",   "SNP",  "INFO,Score",
                                            "MAF", "Cases"};
            REQUIRE(
                commander.parse_command_wrapper("--base-info INFO,Score:0.1"));
            REQUIRE(commander.base_column_check_wrapper(col));
        }
        std::vector<std::string> col = {"P",  "CHR", "BETA", "LOC", "A1",
                                        "A2", "SNP", "INFO", "MAF", "Cases"};
        REQUIRE(commander.parse_command_wrapper("--base Base"));
        SECTION("info default")
        {
            // we have tested parsing of other column, don't bother to do
            // extensive check
            SECTION("with user value")
            {
                auto i = GENERATE(take(10, random(-1.1, 1.1)));
                std::string input = std::to_string(i);
                REQUIRE(commander.parse_command_wrapper("--base-info INFO:"
                                                        + input));
                if (i > 1 || i < 0)
                { REQUIRE_FALSE(commander.base_column_check_wrapper(col)); }
                else
                {
                    REQUIRE(commander.base_column_check_wrapper(col));
                    REQUIRE(commander.get_base().column_index[+BASE_INDEX::INFO]
                            == 7);
                    REQUIRE(commander.get_base().has_column[+BASE_INDEX::INFO]);
                    REQUIRE(commander.get_base_qc().info_score
                            == Approx(misc::Convertor::convert<double>(input)));
                }
            }
            SECTION("without user input")
            {
                REQUIRE(commander.base_column_check_wrapper(col));
                REQUIRE(commander.get_base().column_index[+BASE_INDEX::INFO]
                        == 7);
                REQUIRE(commander.get_base().has_column[+BASE_INDEX::INFO]);
                REQUIRE(commander.get_base_qc().info_score == Approx(0.9));
            }
        }
        SECTION("invalid info format")
        {
            auto i = GENERATE("INFO,0.1", "INFO:0.1,hi");
            REQUIRE(commander.parse_command_wrapper("--base-info "
                                                    + std::string(i)));
            REQUIRE_FALSE(commander.base_column_check_wrapper(col));
        }
    }
    SECTION("maf, maf case or info column not found")
    {
        std::vector<std::string> col = {"P",  "CHR", "BETA", "LOC", "A1",
                                        "A2", "SNP", "INFO", "MAF", "Cases"};
        // info is warning
        using record = std::tuple<std::string, size_t>;
        auto i = GENERATE(table<std::string, size_t>(
            {record {"--base-info INFO_SCORE:0.1", +BASE_INDEX::INFO},
             record {"--base-maf maf_con:0.5", +BASE_INDEX::MAF},
             record {"--base-maf MAF:0.1,case:0.05", +BASE_INDEX::MAF_CASE}}));
        REQUIRE(commander.parse_command_wrapper(std::get<0>(i)));
        CHECK(commander.base_column_check_wrapper(col));
        REQUIRE_FALSE(commander.get_base().has_column[std::get<1>(i)]);
    }
    SECTION("maf")
    {
        std::vector<std::string> col = {"P",  "CHR", "BETA", "LOC", "A1",
                                        "A2", "SNP", "INFO", "MAF", "Cases"};
        SECTION("only control")
        {
            auto i = GENERATE(take(10, random(-1.1, 1.1)));
            std::string input = std::to_string(i);
            REQUIRE(commander.parse_command_wrapper("--base-maf MAF:" + input));
            // maf ain't essential

            if (i > 1.0 || i < 0)
            {
                REQUIRE_FALSE(commander.base_column_check_wrapper(col));
                REQUIRE_FALSE(
                    commander.get_base().has_column[+BASE_INDEX::MAF]);
            }
            else
            {
                REQUIRE(commander.base_column_check_wrapper(col));
                REQUIRE(commander.get_base().has_column[+BASE_INDEX::MAF]);
                REQUIRE(commander.get_base().column_index[+BASE_INDEX::MAF]
                        == 8);
                REQUIRE(commander.get_base_qc().maf
                        == Approx(misc::Convertor::convert<double>(input)));
            }
            REQUIRE_FALSE(
                commander.get_base().has_column[+BASE_INDEX::MAF_CASE]);
        }
        SECTION("with cases")
        {
            SECTION("invalid control")
            {
                REQUIRE(commander.parse_command_wrapper(
                    "--base-maf mAf:0.1,Cases:0.05"));
                REQUIRE(commander.base_column_check_wrapper(col));
                REQUIRE_FALSE(
                    commander.get_base().has_column[+BASE_INDEX::MAF]);
                REQUIRE_FALSE(
                    commander.get_base().has_column[+BASE_INDEX::MAF_CASE]);
            }
            SECTION("Valid control")
            {
                auto i = GENERATE(take(10, random(-1.1, 1.1)));
                std::string input = std::to_string(i);
                REQUIRE(commander.parse_command_wrapper(
                    "--base-maf MAF:0.1,Cases:" + input));
                if (i > 1.0 || i < 0)
                { REQUIRE_FALSE(commander.base_column_check_wrapper(col)); }
                else
                {
                    REQUIRE(commander.base_column_check_wrapper(col));
                    REQUIRE(commander.get_base().has_column[+BASE_INDEX::MAF]);
                    REQUIRE(
                        commander.get_base().has_column[+BASE_INDEX::MAF_CASE]);
                    REQUIRE(commander.get_base().column_index[+BASE_INDEX::MAF]
                            == 8);
                    REQUIRE(
                        commander.get_base().column_index[+BASE_INDEX::MAF_CASE]
                        == 9);
                    REQUIRE(commander.get_base_qc().maf == Approx(0.1));

                    REQUIRE(commander.get_base_qc().maf_case
                            == Approx(misc::Convertor::convert<double>(input)));
                }
            }
        }
    }
}
