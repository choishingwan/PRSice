#include "catch.hpp"
#include "genotype.hpp"
#include "mock_genotype.hpp"

TEST_CASE("process snp")
{
    mockGenotype geno;
    Reporter reporter("log", 60, true);
    geno.set_reporter(&reporter);
    geno.test_init_chr();
    SECTION("check chromosome")
    {
        bool chr_error, sex_error;
        std::string prev = "chr1";
        size_t chr_num = 1;
        SECTION("valid input")
        {
            SECTION("same chromosome as before")
            {
                size_t prev_num = chr_num;
                std::string new_chr = prev;
                REQUIRE(geno.test_check_chr(new_chr, prev, chr_num, chr_error,
                                            sex_error));
                REQUIRE(prev_num == chr_num);
                REQUIRE(prev == new_chr);
            }
            SECTION("new chromosome")
            {
                std::string new_chr = "chr2";
                REQUIRE(geno.test_check_chr(new_chr, prev, chr_num, chr_error,
                                            sex_error));
                REQUIRE(chr_num == 2);
                REQUIRE(prev == new_chr);
            }
        }
        SECTION("invalid chromosomes")
        {
            auto chr = GENERATE(std::string("chrX"), "chromosome");
            auto chr_error = GENERATE(true, false);
            auto sex_error = GENERATE(true, false);
            bool ori_sex = sex_error;
            bool ori_chr = chr_error;
            REQUIRE_FALSE(
                geno.test_check_chr(chr, prev, chr_num, chr_error, sex_error));
            if (chr == "chrX")
            {
                REQUIRE(sex_error);
                REQUIRE(chr_error == ori_chr);
            }
            else
            {
                REQUIRE(chr_error);
                REQUIRE(sex_error == ori_sex);
            }
        }
    }
    SECTION("check rs")
    {
        std::unordered_set<std::string> processed_snps;
        std::unordered_set<std::string> duplicated_snps;
        std::string snp_id, rs_id;
        SECTION("both snp and rs id are . ")
        {
            rs_id = ".";
            snp_id = ".";
            REQUIRE_FALSE(geno.test_check_rs(rs_id, snp_id, processed_snps,
                                             duplicated_snps, &geno));
        }
    }
}
