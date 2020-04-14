#include "binaryplink.hpp"
#include "catch.hpp"
#include "mock_binaryplink.hpp"
#include "reporter.hpp"

TEST_CASE("binary plink load sample")
{
    GenoFile geno;
    geno.file_name = "test";
    Phenotype pheno;
    pheno.ignore_fid = GENERATE(true, false);
    Reporter reporter("log", 60, true);
    auto keep_nonfounder = GENERATE(true, false);
    auto is_ref = GENERATE(true, false);
    SECTION("Normal file")
    {
        // generate a small test file
        std::ofstream mock_fam("test.fam");
        mock_fam << "CAS_1 CAS_1 0 0 2 2" << std::endl;
        mock_fam << "CAS_1 CAS_2 0 0 1 2" << std::endl;
        mock_fam << "CAS_1 CAS_3 CAS_1 CAS_2 1 2" << std::endl;
        mock_binaryplink bplink(geno, pheno, " ", &reporter);
        if (is_ref) { bplink.reference(); }
        bplink.keep_nonfounder(keep_nonfounder);
        // this is ok as we know we haven't check for the bed and bim file yet
        bplink.load_samples(false);
        auto res = bplink.sample_id();
        auto sample_ld = bplink.sample_for_ld();
        auto in_prs = bplink.calculate_prs();
        if (is_ref) { REQUIRE(res.empty()); }
        else
        {
            REQUIRE(res.size() == 3);
            for (size_t i = 0; i < 2; ++i) REQUIRE(res[i].in_regression);
            if (pheno.ignore_fid || keep_nonfounder)
            { REQUIRE(res[2].in_regression); }
            else
            {
                REQUIRE_FALSE(res[2].in_regression);
            }
        }
        REQUIRE(IS_SET(sample_ld.data(), 0));
        REQUIRE(IS_SET(sample_ld.data(), 1));
        if (pheno.ignore_fid) { REQUIRE(IS_SET(sample_ld.data(), 2)); }
        else
        {
            REQUIRE_FALSE(IS_SET(sample_ld.data(), 2));
        }
        for (size_t i = 0; i < 3; ++i) REQUIRE(IS_SET(in_prs.data(), i));
    }
    SECTION("malformed file")
    {
        std::ofstream mock_fam("test.fam");
        mock_fam << "CAS_1 CAS_1 0 0 2 2 3" << std::endl;
        mock_fam << "CAS_1 CAS_2" << std::endl;
        mock_fam << "CAS_1 CAS_3 CAS_1 CAS_2 1 2 3 4 5 6 7" << std::endl;
        mock_binaryplink bplink(geno, pheno, " ", &reporter);
        if (is_ref) { bplink.reference(); }
        bplink.keep_nonfounder(keep_nonfounder);
        // this is ok as we know we haven't check for the bed and bim file yet
        REQUIRE_THROWS(bplink.load_samples(false));
    }
}
