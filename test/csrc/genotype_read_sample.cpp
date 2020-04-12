#include "catch.hpp"
#include "genotype.hpp"
#include "mock_genotype.hpp"

TEST_CASE("initialize sample vectors")
{
    mockGenotype geno;
    SECTION("without any samples")
    {
        REQUIRE_THROWS(geno.init_sample_vectors());
    }
    SECTION("with samples")
    {
        auto n_sample = GENERATE(range(1023, 1025), range(127, 129));
        geno.set_sample(static_cast<uintptr_t>(n_sample));
        const uintptr_t expected_size =
            BITCT_TO_WORDCT(static_cast<uintptr_t>(n_sample));
        REQUIRE_NOTHROW(geno.init_sample_vectors());
        REQUIRE(geno.num_ambig() == 0);
        REQUIRE(geno.num_male() == 0);
        REQUIRE(geno.num_female() == 0);
        REQUIRE(geno.num_founder() == 0);
        REQUIRE(geno.num_nonfounder() == 0);
        REQUIRE(geno.num_sample() == 0);
        REQUIRE(geno.sample_for_ld().size() == expected_size);
        REQUIRE(geno.calculate_prs().size() == expected_size);
    }
}

TEST_CASE("Sample object generation")
{
    mockGenotype geno;
    std::unordered_set<std::string> founder_info;
    std::unordered_set<std::string> processed_samples;
    std::vector<std::string> duplicated_sample_names;
    std::vector<Sample_ID> sample_storage;
    // initialize the vector for later use
    auto n_sample = GENERATE(range(63, 65));
    // auto n_sample = GENERATE(range(63, 65),range(1023, 1025), range(127,
    // 129));
    geno.set_sample(static_cast<uintptr_t>(n_sample));
    geno.test_init_sample_vectors();
    // test handling of ignore fid
    auto ignore_fid = GENERATE(true, false);
    auto mock_file = std::make_unique<std::istringstream>(std::ios_base::ate
                                                          | std::ios_base::app);
    auto delim = GENERATE(" ", "-", "\t");
    geno.set_delim(delim);
    geno.set_ignore_fid(ignore_fid);
    SECTION("without family")
    {
        // select
        // extract
        // duplicates
        SECTION("all included")
        {
            // test different file format
            auto n_col = GENERATE(range(3, 10));
            // test handling of sex chromosome
            auto has_sex = GENERATE(true, false);
            std::random_device rd;
            std::mt19937 gen(rd());
            // make fake id
            std::uniform_int_distribution<> dis(1, static_cast<int>(n_sample)
                                                       * 1000);
            std::uniform_int_distribution<> sex_rand(1, 3);
            std::unordered_set<int> dup_gen;
            std::vector<std::string> included_samples;
            std::string file_input;
            size_t exp_male = 0, exp_female = 0, exp_ambig = 0;
            const size_t fid_idx = 0;
            const size_t iid_idx = ignore_fid ? 0 : 1;
            size_t sex_idx = 0;
            for (size_t i = 0; i < static_cast<uintptr_t>(n_sample); ++i)
            {
                auto b_in_prs = geno.calculate_prs();
                auto b_in_ld = geno.sample_for_ld();
                REQUIRE_FALSE(IS_SET(b_in_prs.data(), i));
                REQUIRE_FALSE(IS_SET(b_in_ld.data(), i));
                std::vector<std::string> mock_input(
                    static_cast<uintptr_t>(n_col), std::string("@"));
                int temp = dis(gen);
                // generate non-duplicated samples
                while (dup_gen.find(temp) != dup_gen.end()) temp = dis(gen);
                dup_gen.insert(temp);
                mock_input[0] = "F" + std::to_string(temp);
                std::string exp = mock_input[0];
                if (!ignore_fid) { mock_input[1] = exp; }
                included_samples.push_back(exp);
                if (has_sex)
                {
                    auto sex = sex_rand(gen);
                    switch (sex)
                    {
                    case 1:
                        mock_input[2] = "1";
                        ++exp_male;
                        break;
                    case 2:
                        mock_input[2] = "2";
                        ++exp_female;
                        break;
                    case 3:
                        mock_input[2] = "NA";
                        ++exp_ambig;
                        break;
                    }
                    sex_idx = 2;
                }
                else
                {
                    sex_idx = ~size_t(0);
                }
                REQUIRE_NOTHROW(geno.test_gen_sample(
                    fid_idx, iid_idx, sex_idx, ~size_t(0), ~size_t(0), i,
                    founder_info, "", mock_input, sample_storage,
                    processed_samples, duplicated_sample_names));
                REQUIRE(duplicated_sample_names.empty());
                auto in_prs = geno.calculate_prs();
                auto in_ld = geno.sample_for_ld();
                REQUIRE(IS_SET(in_prs.data(), i));
                REQUIRE(IS_SET(in_ld.data(), i));
                REQUIRE(geno.num_male() == exp_male);
                REQUIRE(geno.num_female() == exp_female);
                REQUIRE(geno.num_ambig_sex() == exp_ambig);
                // always have IID
                REQUIRE(sample_storage.back().IID == exp);
                REQUIRE(sample_storage.back().in_regression);
                if (!ignore_fid) { REQUIRE(sample_storage.back().FID == exp); }
                else
                {
                    REQUIRE(sample_storage.back().FID.empty());
                }
            }
        }
    }
    SECTION("with family") {}
}
