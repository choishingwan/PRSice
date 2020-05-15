#include "catch.hpp"
#include "mock_genotype.hpp"

TEST_CASE("Prepare PRSice")
{
    // If there are no SNPs, false
    mockGenotype geno;
    Reporter reporter("log", 60, true);
    geno.set_reporter(&reporter);
    SECTION("No SNP") { REQUIRE_FALSE(geno.prepare_prsice()); }
    const std::string ref_allele = "A", alt_allele = "G";
    const double stat = 1.96;
    // p -> file -> byte
    // category -> file -> byte
    // byte within same file should never be identical
    std::vector<SNP> input = {
        SNP("rs1", 1, 123, ref_allele, alt_allele, 0, 1, stat, 0.05, 0, 0.01),
        SNP("rs2", 1, 123, ref_allele, alt_allele, 0, 2, stat, 0.05, 0, 0.02),
        SNP("rs3", 1, 123, ref_allele, alt_allele, 1, 1, stat, 0.05, 0, 0.03),
        SNP("rs4", 1, 123, ref_allele, alt_allele, 0, 3, stat, 0.07, 0, 0.04),
        SNP("rs5", 1, 123, ref_allele, alt_allele, 0, 1, stat, 0.07, 1, 0.05)};
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(input.begin(), input.end(), g);
    for (auto&& snp : input) { geno.load_snp(snp); }
    SECTION("Very tiny thresholds")
    {
        // we will sort by p-values
        geno.set_very_small_thresholds();
        geno.prepare_prsice();
        auto&& snp = geno.existed_snps();
        std::vector<std::string> observed;
        for (auto&& s : snp)
        {
            observed.push_back(s.rs());
            REQUIRE(s.p_value() == s.get_threshold());
        }
        REQUIRE_THAT(observed, Catch::Equals<std::string>(
                                   {"rs1", "rs2", "rs3", "rs5", "rs4"}));
    }
    SECTION("Normal preparation")
    {
        geno.prepare_prsice();
        auto&& snp = geno.existed_snps();
        std::vector<std::string> observed;
        for (auto&& s : snp) { observed.push_back(s.rs()); }
        REQUIRE_THAT(observed, Catch::Equals<std::string>(
                                   {"rs1", "rs2", "rs4", "rs3", "rs5"}));
    }
}
