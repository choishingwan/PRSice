#include "binaryplink.hpp"
#include "catch.hpp"
#include "mock_binaryplink.hpp"

// difficulties is to:
// 1. Generate samples with known amount of missingness & MAF
// 2. And to account for founder / non-founder status
TEST_CASE("Geno and MAF filtering")
{
    auto n_sample = GENERATE(repeat(5, range(124ul, 126ul)));
    // use random generators to generate our founder and genotypes
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()}; // Generates random integers
    std::uniform_int_distribution<size_t> dist {0, 3}, rand_founder {1, 10};
    auto gen = [&dist, &mersenne_engine]() { return dist(mersenne_engine); };
    // use code 0, 1, 2, 3 to represent dosage of genotype and missing (3)
    std::vector<size_t> sample_genotype(n_sample, 0);
    std::generate(begin(sample_genotype), end(sample_genotype), gen);
    // now generate random founder status
    // set ~70% of samples being founder
    auto founder_gen = [&rand_founder, &mersenne_engine]() {
        return rand_founder(mersenne_engine) > 3;
    };
    std::vector<bool> founder(n_sample, true);
    std::generate(begin(founder), end(founder), founder_gen);
    auto n_missing =
        std::count(sample_genotype.begin(), sample_genotype.end(), 3);
    double exp_geno =
        static_cast<double>(n_missing) / static_cast<double>(n_sample);
    size_t ref_founder_ct = 0, het_founder_ct = 0, alt_founder_ct = 0,
           miss_founder_ct = 0, ref_exp_ct = 0, het_exp_ct = 0, alt_exp_ct = 0,
           miss_exp_ct = 0, n_founder = 0;
    for (size_t i = 0; i < n_sample; ++i)
    {
        switch (sample_genotype[i])
        {
        case 0:
            ref_founder_ct += founder[i];
            ++ref_exp_ct;
            break;
        case 1:
            ++het_exp_ct;
            het_founder_ct += founder[i];
            break;
        case 2:
            ++alt_exp_ct;
            alt_founder_ct += founder[i];
            break;
        case 3:
            ++miss_exp_ct;
            miss_founder_ct += founder[i];
            break;
        }
        n_founder += founder[i];
    }
    double exp_maf =
        (2.0 * alt_founder_ct + het_founder_ct)
        / (2.0 * (ref_founder_ct + het_founder_ct + alt_founder_ct));
    exp_maf = (exp_maf > 0.5) ? 1 - exp_maf : exp_maf;
    // need to initialize the gneotype vectors
    mock_binaryplink plink;
    plink.set_sample(n_sample);
    plink.test_init_sample_vectors();
    plink.set_founder_vector(founder);
    plink.set_sample_vector(n_sample);
    plink.test_post_sample_read_init();
    Reporter reporter("log", 60, true);
    plink.set_reporter(&reporter);

    SECTION("Directly test algorithm from plink")
    {
        uint32_t ref_ct = 0, het_ct = 0, alt_ct = 0, miss = 0, ref_get_ct = 0,
                 het_get_ct = 0, alt_get_ct = 0, miss_get_ct = 0;
        plink.test_single_marker_freqs_and_hwe(sample_genotype, &ref_ct,
                                               &het_ct, &alt_ct, &ref_get_ct,
                                               &het_get_ct, &alt_get_ct);
        miss = static_cast<uint32_t>(n_sample) - ref_ct - het_ct - alt_ct;
        miss_get_ct = static_cast<uint32_t>(n_founder) - ref_get_ct - het_get_ct
                      - alt_get_ct;
        REQUIRE(ref_ct == ref_exp_ct);
        REQUIRE(het_ct == het_exp_ct);
        REQUIRE(alt_ct == alt_exp_ct);
        REQUIRE(miss == miss_exp_ct);
        REQUIRE(ref_get_ct == ref_founder_ct);
        REQUIRE(het_get_ct == het_founder_ct);
        REQUIRE(alt_get_ct == alt_founder_ct);
        REQUIRE(miss_get_ct == miss_founder_ct);
    }
    SECTION("Test with file")
    {
        plink.gen_fake_bed(sample_genotype, "filter_sim");
        QCFiltering filter_info;
        filter_info.maf = GENERATE(0.05, 0.5, 0.0);
        filter_info.geno = GENERATE(1.0, 0.5, 0.0);
        bool filtered =
            filter_info.maf > exp_maf || filter_info.geno < exp_geno;
        REQUIRE(plink.test_calc_freq_gen_inter(filter_info));
        auto res = plink.existed_snps();
        if (filtered) { REQUIRE(res.empty()); }
        else
        {
            REQUIRE(res.size() == 1);
            size_t ref_ct = 0, het_ct = 0, alt_ct = 0, miss = 0;
            res.front().get_counts(ref_ct, het_ct, alt_ct, miss, false);
            REQUIRE(ref_ct == ref_founder_ct);
            REQUIRE(het_ct == het_founder_ct);
            REQUIRE(alt_ct == alt_founder_ct);
            REQUIRE(miss == miss_founder_ct);
        }
    }
}
