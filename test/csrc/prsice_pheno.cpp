#include "catch.hpp"
#include "mock_prsice.hpp"
#include "prsice.hpp"
TEST_CASE("Simple Phenotype check")
{
    Reporter reporter("log", 60, true);
    Phenotype pheno;
    SECTION("Get Phenotype index")
    {
        mock_prsice prsice(&reporter);
        std::string input = "FID IID SCZ MDD Bipolar SCZ";
        std::vector<std::string_view> col_sv = misc::tokenize(input);
        SECTION("With duplicated target phenotype column")
        {
            REQUIRE_THROWS(
                mock_prsice::test_get_pheno_idx(col_sv, pheno, "SCZ"));
        }
        SECTION("Phenotyep not found in header")
        {
            auto [idx, found] =
                mock_prsice::test_get_pheno_idx(col_sv, pheno, "ANX");
            REQUIRE(idx == 0);
            REQUIRE_FALSE(found);
        }
        SECTION("Phenotype found")
        {
            auto [idx, found] =
                mock_prsice::test_get_pheno_idx(col_sv, pheno, "MDD");
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
            REQUIRE_THROWS(mock_prsice::test_parse_pheno_header(
                std::move(input), pheno, reporter));
        }
        SECTION("User did not provide a phenotype header")
        {
            using record = std::tuple<std::string, std::string>;
            auto headers = GENERATE(table<std::string, std::string>(
                {record {"FID IID SCZ MDD Bipolar SCZ", "SCZ"},
                 record {"FID IID 1 2 3", "Phenotype"}}));
            auto input =
                std::make_unique<std::istringstream>(std::get<0>(headers));
            Phenotype original = pheno;
            REQUIRE_NOTHROW(mock_prsice::test_parse_pheno_header(
                std::move(input), pheno, reporter));
            if (original.ignore_fid)
            {
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"IID"}));
            }
            else
            {
                REQUIRE_THAT(pheno.pheno_col, Catch::Equals<std::string>(
                                                  {std::get<1>(headers)}));
            }
            size_t exp_idx = !pheno.ignore_fid + 1;
            REQUIRE_THAT(pheno.pheno_col_idx, Catch::Equals<size_t>({exp_idx}));
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
                REQUIRE_THROWS(mock_prsice::test_parse_pheno_header(
                    std::move(input), pheno, reporter));
            }
            SECTION("Some column were found")
            {
                pheno.pheno_col = {"ALZ", "MDD"};
                mock_prsice prsice(prs_info, pthres, pheno, perm, output,
                                   &reporter);
                REQUIRE_NOTHROW(mock_prsice::test_parse_pheno_header(
                    std::move(input), pheno, reporter));
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"ALZ", "MDD"}));
                REQUIRE_THAT(pheno.skip_pheno,
                             Catch::Equals<bool>({true, false}));
                REQUIRE_THAT(pheno.pheno_col_idx,
                             Catch::Equals<size_t>({0, 3}));
            }
        }
    }
    SECTION("Full pheno check")
    {
        // ignore fid don't have any effect on  when phenotype file isn't
        // provided
        pheno.ignore_fid = GENERATE(true, false);
        pheno.binary = {true};
        SECTION("Did not provide phenotype file")
        {
            auto no_regress = GENERATE(true, false);
            REQUIRE_NOTHROW(PRSice::pheno_check(no_regress, pheno, reporter));
            if (no_regress)
            {
                // with no regress, user are not required to provide binary
                // flag, thus we will just add one
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"PlaceHolder"}));
                REQUIRE_THAT(pheno.pheno_col_idx,
                             Catch::Equals<size_t>({~size_t(0)}));
                REQUIRE(pheno.binary.size() == pheno.pheno_col_idx.size());
                REQUIRE_THAT(pheno.skip_pheno, Catch::Equals<bool>({false}));
            }
            else
            {
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"Phenotype"}));
                REQUIRE_THAT(pheno.pheno_col_idx,
                             Catch::Equals<size_t>({~size_t(0)}));
                REQUIRE_THAT(pheno.skip_pheno, Catch::Equals<bool>({false}));
            }
        }
        SECTION("Provided a phenotype file")
        {
            std::ofstream pheno_file("pheno_check.prs");
            pheno_file << "FID IID SCZ MDD" << std::endl;
            pheno_file << "ID1 ID1 0 1" << std::endl;
            pheno_file.close();
            pheno.pheno_file = "pheno_check.prs";
            SECTION("Did not provide a header")
            {
                PRSice::pheno_check(false, pheno, reporter);
                REQUIRE_THAT(pheno.pheno_col_idx,
                             Catch::Equals<size_t>({1ul + !pheno.ignore_fid}));
                // because header isn't digit, we will use the column header as
                // the phenotype name
                if (pheno.ignore_fid)
                {
                    REQUIRE_THAT(pheno.pheno_col,
                                 Catch::Equals<std::string>({"IID"}));
                }
                else
                {
                    REQUIRE_THAT(pheno.pheno_col,
                                 Catch::Equals<std::string>({"SCZ"}));
                }
                REQUIRE_THAT(pheno.skip_pheno, Catch::Equals<bool>({false}));
            }
            SECTION("provided a header")
            {
                pheno.pheno_col = {"MDD"};
                PRSice::pheno_check(false, pheno, reporter);
                REQUIRE_THAT(pheno.pheno_col_idx, Catch::Equals<size_t>({3}));
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"MDD"}));
                REQUIRE_THAT(pheno.skip_pheno, Catch::Equals<bool>({false}));
            }
        }
    }
}
