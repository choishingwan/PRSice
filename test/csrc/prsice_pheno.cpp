#include "catch.hpp"
#include "mock_prsice.hpp"
#include "prsice.hpp"
TEST_CASE("Simple Phenotype check")
{
    Reporter reporter("log", 60, true);
    SECTION("Get Phenotype index")
    {
        mock_prsice prsice(&reporter);
        std::string input = "FID IID SCZ MDD Bipolar SCZ";
        std::vector<std::string_view> col_sv = misc::tokenize(input);

        SECTION("With duplicated target phenotype column")
        {
            REQUIRE_THROWS(prsice.test_get_pheno_idx(col_sv, "SCZ"));
        }
        SECTION("Phenotyep not found in header")
        {
            auto [idx, found] = prsice.test_get_pheno_idx(col_sv, "ANX");
            REQUIRE(idx == 0);
            REQUIRE_FALSE(found);
        }
        SECTION("Phenotype found")
        {
            auto [idx, found] = prsice.test_get_pheno_idx(col_sv, "MDD");
            REQUIRE(found);
            REQUIRE(idx == 3);
        }
    }
    SECTION("Parse phenotype file header")
    {
        PThresholding pthres;
        auto output = "PRSice";
        Permutations perm;
        CalculatePRS prs_info;
        Phenotype pheno;
        pheno.ignore_fid = GENERATE(true, false);
        prs_info.no_regress = false;
        SECTION("Invalid phenotype file")
        {
            auto invalid = GENERATE("\nFID IID Ok\n", "Invalid");
            auto input = std::make_unique<std::istringstream>(invalid);
            mock_prsice prsice(prs_info, pthres, pheno, perm, output,
                               &reporter);
            REQUIRE_THROWS(prsice.test_parse_pheno_header(std::move(input)));
        }
        SECTION("User did not provide a phenotype header")
        {
            using record = std::tuple<std::string, std::string>;
            auto headers = GENERATE(table<std::string, std::string>(
                {record {"FID IID SCZ MDD Bipolar SCZ", "SCZ"},
                 record {"FID IID 1 2 3", "Phenotype"}}));
            auto input =
                std::make_unique<std::istringstream>(std::get<0>(headers));
            mock_prsice prsice(prs_info, pthres, pheno, perm, output,
                               &reporter);
            REQUIRE_NOTHROW(prsice.test_parse_pheno_header(std::move(input)));
            auto res = prsice.get_pheno_info();
            if (pheno.ignore_fid)
            {
                REQUIRE_THAT(res.pheno_col,
                             Catch::Equals<std::string>({"IID"}));
            }
            else
            {
                REQUIRE_THAT(res.pheno_col, Catch::Equals<std::string>(
                                                {std::get<1>(headers)}));
            }
            size_t exp_idx = !pheno.ignore_fid + 1;
            REQUIRE_THAT(res.pheno_col_idx, Catch::Equals<size_t>({exp_idx}));
        }
        SECTION("User provided multiple header")
        {
            auto input = std::make_unique<std::istringstream>(
                "FID IID SCZ MDD Bipolar SCZ");
            SECTION("Non of the column were found")
            {
                pheno.pheno_col = {"ALZ", "ADD", "ANX"};
                mock_prsice prsice(prs_info, pthres, pheno, perm, output,
                                   &reporter);
                REQUIRE_THROWS(
                    prsice.test_parse_pheno_header(std::move(input)));
            }
            SECTION("Some column were found")
            {
                pheno.pheno_col = {"ALZ", "MDD"};
                mock_prsice prsice(prs_info, pthres, pheno, perm, output,
                                   &reporter);
                REQUIRE_NOTHROW(
                    prsice.test_parse_pheno_header(std::move(input)));
                auto res = prsice.get_pheno_info();
                REQUIRE_THAT(res.pheno_col,
                             Catch::Equals<std::string>({"ALZ", "MDD"}));
                REQUIRE_THAT(res.skip_pheno,
                             Catch::Equals<bool>({true, false}));
                REQUIRE_THAT(res.pheno_col_idx, Catch::Equals<size_t>({0, 3}));
            }
        }
    }
}
