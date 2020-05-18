#include "catch.hpp"
#include "mock_binaryplink.hpp"
#include "plink_common.hpp"
TEST_CASE("plink read_genotype")
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
    const uintptr_t target_ctl = BITCT_TO_WORDCT(n_target_sample);
    const uintptr_t ref_ctv2 = 2 * ref_ctl;
    const uintptr_t target_ctv2 = 2 * target_ctl;
    std::vector<uintptr_t> expected_geno(ref_ctv2, 0);
    std::vector<uintptr_t> expected_target_geno(target_ctv2, 0),
        expected_memory(target_ctv2, 0);
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
    geno_idx = 0;
    size_t mem_geno_idx = 0;
    for (size_t i = 0; i < n_target_sample; ++i)
    {
        if (target_founder[i])
        {
            switch (taget_genotype[i])
            {
            case 0: break;
            case 1: SET_BIT(geno_idx + 1, expected_target_geno.data()); break;
            case 2:
                SET_BIT(geno_idx + 1, expected_target_geno.data());
                SET_BIT(geno_idx, expected_target_geno.data());
                break;
            case 3: SET_BIT(geno_idx, expected_target_geno.data()); break;
            }
            geno_idx += 2;
        }
        switch (taget_genotype[i])
        {
        case 0: break;
        case 1: SET_BIT(mem_geno_idx + 1, expected_memory.data()); break;
        case 2:
            SET_BIT(mem_geno_idx + 1, expected_memory.data());
            SET_BIT(mem_geno_idx, expected_memory.data());
            break;
        case 3: SET_BIT(mem_geno_idx, expected_memory.data()); break;
        }
        mem_geno_idx += 2;
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
    SECTION("read genotype")
    {
        std::vector<uintptr_t> observed_geno(ref_ctv2, 0);
        std::vector<uintptr_t> observed_target_geno(target_ctv2, 0);
        // read genotype should always only read from reference, therefore
        // should fail for target
        REQUIRE_NOTHROW(target.test_read_genotype(
            cur_snp, observed_target_geno.data(), false));
        REQUIRE_THAT(observed_target_geno,
                     Catch::Equals<uintptr_t>(expected_target_geno));
        REQUIRE_NOTHROW(
            reference.test_read_genotype(cur_snp, observed_geno.data(), true));
        REQUIRE_THAT(observed_geno, Catch::Equals<uintptr_t>(expected_geno));
    }
    SECTION("load genotype to memory")
    {
        target.load_genotype_to_memory();
        auto&& snps = target.existed_snps();
        // SNP current geno should now point to memory with the genotype
        REQUIRE_FALSE(snps.front().current_genotype() == nullptr);
        // we know size of cur_snp, which is target_ctv2
        auto&& snp_memory = snps.front().current_genotype();
        std::vector<uintptr_t> observed(snp_memory, snp_memory + target_ctv2);
        REQUIRE_THAT(observed, Catch::Equals<uintptr_t>(expected_memory));
    }
}
/*
void generate_expected_prs(const std::vector<size_t>& genotype,
                           const std::vector<bool>& selected,
                           const double ref_maf, const MODEL& genetic_model,
                           const MISSING_SCORE& missing_score,
                           const SCORING& scoring, const bool flipped,
                           const bool use_ref_maf, Genotype& geno,
                           std::vector<PRS>& expected_prs,
                           misc::RunningStat& rs)
{
}

TEST_CASE("Test PRS Read score")
{
    // need to prepare the followings:
    // 1: Generate samples
    //    - Need to do selection
    // 2: Fill the SNPs with proper byte pos
    // 3. Generate samples with known counts
    // 4. Preload genotype vs read on the fly
    // 5. Flip or no flip
    // 6. Different missing scores
    // 7. Use ref_maf
    auto genetic_model = GENERATE(MODEL::ADDITIVE, MODEL::DOMINANT,
                                  MODEL::RECESSIVE, MODEL::HETEROZYGOUS);
    auto missing_score =
        GENERATE(MISSING_SCORE::CENTER, MISSING_SCORE::SET_ZERO,
                 MISSING_SCORE::MEAN_IMPUTE, MISSING_SCORE::IMPUTE_CONTROL);
    auto use_ref_maf = GENERATE(true, false);
    size_t n_sample = 1000ul;
    Reporter reporter("log", 60, true);
    mock_binaryplink plink;
    plink.set_reporter(&reporter);
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()}; // Generates random integers
    std::uniform_int_distribution<size_t> dist {0, 3}, rand_selected {1, 10};
    auto gen = [&dist, &mersenne_engine]() { return dist(mersenne_engine); };
    std::vector<size_t> simulated_genotype(n_sample, 0);
    std::generate(begin(simulated_genotype), end(simulated_genotype), gen);
    auto selection_gen = [&rand_selected, &mersenne_engine]() {
        return rand_selected(mersenne_engine) > 3;
    };
    // simulate different founder status
    std::vector<bool> selected(n_sample, true);
    std::generate(begin(selected), end(selected), selection_gen);
    // expected vector:
    const uintptr_t sample_ctl = BITCT_TO_WORDCT(n_sample);
    const uintptr_t sample_ctv2 = 2 * sample_ctl;
    std::vector<uintptr_t> genotype_data(sample_ctv2, 0);
    size_t geno_idx = 0;
    for (size_t i = 0; i < n_sample; ++i)
    {
        switch (simulated_genotype[i])
        {
        case 0: break;
        case 1: SET_BIT(geno_idx + 1, genotype_data.data()); break;
        case 2:
            SET_BIT(geno_idx + 1, genotype_data.data());
            SET_BIT(geno_idx, genotype_data.data());
            break;
        case 3: SET_BIT(geno_idx, genotype_data.data()); break;
        }
        geno_idx += 2;
    }
    plink.set_sample_vector(selected);
    plink.set_founder_vector(selected);
    plink.test_post_sample_read_init();
    uint32_t ref_homcom = GENERATE(take(1, random(0, 1000)));
    uint32_t ref_het = GENERATE(take(1, random(0, 1000)));
    uint32_t ref_homrar = GENERATE(take(1, random(0, 1000)));
    uint32_t ref_missing = GENERATE(take(1, random(0, 10)));
    double ref_maf =
        static_cast<double>(ref_het + 2 * ref_homrar)
        / static_cast<double>(2.0 * (ref_homcom + ref_het + ref_homrar));
    // generate case control sample and do the exclude_std thingy
}
*/
