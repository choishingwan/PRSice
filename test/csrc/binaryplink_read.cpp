#include "catch.hpp"
#include "mock_binaryplink.hpp"
#include "plink_common.hpp"

TEST_CASE("read plink genotype for LD")
{
    // first generate two object, the target and the reference
    Reporter reporter("log", 60, true);
    mock_binaryplink target, reference;
    target.set_reporter(&reporter);
    reference.set_reporter(&reporter);
    // use more sample in ref so that if we read the wrong vector, we will
    // fatally crushed the test
    auto n_target_sample = 126ul;
    auto n_ref_sample = 1025ul;
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()}; // Generates random integers
    std::uniform_int_distribution<size_t> dist {0, 3}, rand_founder {1, 10};
    auto gen = [&dist, &mersenne_engine]() { return dist(mersenne_engine); };
    std::vector<size_t> taget_genotype(n_target_sample, 0);
    std::generate(begin(taget_genotype), end(taget_genotype), gen);
    std::vector<size_t> ref_genotype(n_ref_sample, 0);
    std::generate(begin(ref_genotype), end(ref_genotype), gen);
    // two SNPs, to make sure we do the skipping right
    std::vector<size_t> ref_second_genotype(n_ref_sample, 0);
    std::generate(begin(ref_second_genotype), end(ref_second_genotype), gen);
    auto founder_gen = [&rand_founder, &mersenne_engine]() {
        return rand_founder(mersenne_engine) > 3;
    };
    // simulate different founder status
    std::vector<bool> target_founder(n_target_sample, true);
    std::generate(begin(target_founder), end(target_founder), founder_gen);
    std::vector<bool> ref_founder(n_ref_sample, true);
    std::generate(begin(ref_founder), end(ref_founder), founder_gen);
    // expected vector:
    const uintptr_t ref_ctl = BITCT_TO_WORDCT(n_ref_sample);
    const uintptr_t ref_ctv2 = 2 * ref_ctl;
    std::vector<uintptr_t> expected_geno(ref_ctv2, 0);
    size_t geno_idx = 0;
    for (size_t i = 0; i < n_ref_sample; ++i)
    {
        if (ref_founder[i])
        {
            switch (ref_second_genotype[i])
            {
            case 0: break;
            case 1: SET_BIT(geno_idx + 1, expected_geno.data()); break;
            case 2:
                SET_BIT(geno_idx + 1, expected_geno.data());
                SET_BIT(geno_idx, expected_geno.data());
                break;
            case 3: SET_BIT(geno_idx, expected_geno.data()); break;
            }
            geno_idx += 2;
        }
    }


    // now init the data

    target.set_sample(n_target_sample);
    target.test_init_sample_vectors();
    target.set_founder_vector(target_founder);
    target.set_sample_vector(n_target_sample);
    target.test_post_sample_read_init();
    target.gen_fake_bed(taget_genotype, "read_sim_target");


    reference.set_sample(n_ref_sample);
    reference.test_init_sample_vectors();
    reference.set_founder_vector(ref_founder);
    reference.set_sample_vector(n_ref_sample);
    reference.test_post_sample_read_init();
    reference.gen_fake_bed(
        std::vector<std::vector<size_t>> {ref_genotype, ref_second_genotype},
        "read_sim_ref");
    // now set the target SNP data correctly

    auto ref_sample_ct4 = (n_ref_sample + 3) / 4;
    auto ref_byte = static_cast<std::streampos>(3 + (ref_sample_ct4));
    auto target_byte = static_cast<std::streampos>(3);
    SNP cur_snp("rs123", 1, 1, "A", "C", 0, target_byte);
    cur_snp.update_file(0, ref_byte, true);
    target.manual_load_snp(cur_snp);
    std::vector<uintptr_t> observed_geno(ref_ctv2, 0);
    // read genotype should always only read from reference, therefore should
    // fail for target
    REQUIRE_THROWS(target.test_read_genotype(observed_geno.data(), cur_snp));
    REQUIRE_NOTHROW(
        reference.test_read_genotype(observed_geno.data(), cur_snp));
    REQUIRE_THAT(observed_geno, Catch::Equals<uintptr_t>(expected_geno));
}
