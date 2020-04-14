#include "binaryplink.hpp"
#include "catch.hpp"
#include "mock_binaryplink.hpp"


TEST_CASE("check bed format")
{
    mock_binaryplink plink;
    SECTION("modern format")
    {
        bool mordan_snp = GENERATE(true, false);
        bool sample_major = GENERATE(true, false);
        std::string name = "header_check.bed";
        size_t num_snp = GENERATE(take(5, random(1ul, 100ul)));
        size_t num_sample = GENERATE(take(5, random(1ul, 100ul)));
        plink.gen_bed_head(name, num_sample, num_snp, mordan_snp, sample_major);
        plink.set_sample(num_sample);
        uintptr_t offset;
        if (!sample_major)
        {
            REQUIRE_NOTHROW(plink.test_check_bed(name, num_snp, offset));
            if (mordan_snp) { REQUIRE(offset == 3); }
            else
            {
                REQUIRE(offset == 1);
            }
        }
        else
        {
            REQUIRE_THROWS(plink.test_check_bed(name, num_snp, offset));
        }
    }
}
