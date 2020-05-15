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

TEST_CASE("Build membership matrix")
{
    // we need maybe 10 SNPs with different regions
    mockGenotype geno;
    Reporter reporter("log", 60, true);
    geno.set_reporter(&reporter);
    auto num_sets = 65ul;
    std::random_device rd;
    std::mt19937 mersenne_engine(rd());
    std::uniform_real_distribution<double> statistic(0.0, 1.0);
    std::uniform_int_distribution<size_t> set_idx(0, num_sets - 1);
    auto stats = [&statistic, &mersenne_engine]() {
        return statistic(mersenne_engine);
    };
    std::vector<std::vector<size_t>> expected_idx(num_sets);
    std::vector<std::string> fake_names(num_sets);
    for (size_t i = 0; i < num_sets; ++i)
    { fake_names[i] = "set" + std::to_string(i); }
    const size_t required_size = BITCT_TO_WORDCT(num_sets);
    std::ostringstream expected_output;
    expected_output << "CHR\tSNP\tBP\tP";
    for (size_t i = 0; i < fake_names.size(); ++i)
    { expected_output << "\t" << fake_names[i]; }
    expected_output << "\n";
    for (size_t i_snp = 0; i_snp < 10; ++i_snp)
    {
        auto cur_snp = SNP("rs" + std::to_string(i_snp), 1, 123, "A", "T",
                           stats(), stats());
        std::unordered_set<size_t> used;
        auto&& flag = cur_snp.get_flag();
        flag.resize(required_size);
        for (size_t set = 0; set < 30; ++set)
        {
            auto idx = set_idx(mersenne_engine);
            if (used.find(idx) == used.end())
            {
                SET_BIT(idx, flag.data());
                expected_idx[idx].push_back(i_snp);
                used.insert(idx);
            }
        }
        expected_output << cur_snp.chr() << "\t" << cur_snp.rs() << "\t"
                        << cur_snp.loc() << "\t" << cur_snp.p_value();
        for (size_t i_set = 0; i_set < num_sets; ++i_set)
        { expected_output << "\t" << IS_SET(flag.data(), i_set); }
        expected_output << "\n";
        geno.load_snp(cur_snp);
    }
    std::ostringstream fake_out;
    SECTION("Wrong set number")
    {
        REQUIRE_THROWS(geno.build_membership_matrix(num_sets + 1, fake_names,
                                                    true, fake_out));
    }
    SECTION("Correct input")
    {
        auto res =
            geno.build_membership_matrix(num_sets, fake_names, true, fake_out);
        REQUIRE(res.size() == expected_idx.size());
        for (size_t i = 0; i < res.size(); ++i)
        {
            REQUIRE_THAT(res[i],
                         Catch::UnorderedEquals<size_t>({expected_idx[i]}));
        }
        REQUIRE(fake_out.str() == expected_output.str());
    }
}
