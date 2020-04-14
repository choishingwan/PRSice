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
        auto n_sample = GENERATE(range(1023, 1025));
        geno.set_sample(static_cast<uintptr_t>(n_sample));
        const uintptr_t expected_size =
            BITCT_TO_WORDCT(static_cast<uintptr_t>(n_sample));
        REQUIRE_NOTHROW(geno.init_sample_vectors());
        REQUIRE(geno.num_ambig_sex() == 0);
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
    auto n_sample = GENERATE(range(127, 129));
    geno.set_sample(static_cast<uintptr_t>(n_sample));
    geno.test_init_sample_vectors();
    // test handling of ignore fid
    auto ignore_fid = GENERATE(true, false);
    geno.set_ignore_fid(ignore_fid);
    auto keep_nonfounder = GENERATE(true, false);
    geno.set_keep_nonfounder(keep_nonfounder);
    // use funny symbol to make sure it is set and used
    auto delim = "@";
    geno.set_delim(delim);
    // select
    // extract
    // duplicates
    // test different file format
    auto n_col = GENERATE(take(2, range(5, 20)));
    // test handling of sex chromosome
    auto has_sex = GENERATE(true, false);
    auto has_father_id = GENERATE(true, false);
    auto has_mother_id = GENERATE(true, false);
    auto remove_samples = GENERATE(true, false);
    std::random_device rd;
    std::mt19937 gen(rd());
    // make fake id
    std::uniform_int_distribution<> dis(1, static_cast<int>(n_sample) * 1000);
    std::uniform_int_distribution<> sex_rand(1, 3);
    std::uniform_int_distribution<> remove_rand(0, 1);
    std::uniform_int_distribution<> has_father(0, 1);
    std::uniform_int_distribution<> has_mother(0, 1);
    std::uniform_int_distribution<> is_duplicate(0, 100);
    std::unordered_set<int> dup_gen;
    std::vector<std::string> included_samples;
    std::string file_input;
    size_t exp_male = 0, exp_female = 0, exp_ambig = 0;
    const size_t fid_idx = 0;
    const size_t iid_idx = ignore_fid ? 0 : 1;
    size_t sex_idx = 0;
    geno.change_sample_selection(remove_samples);
    for (size_t i = 0; i < static_cast<uintptr_t>(n_sample); ++i)
    {
        auto b_in_prs = geno.calculate_prs();
        auto b_in_ld = geno.sample_for_ld();
        REQUIRE_FALSE(IS_SET(b_in_prs.data(), i));
        REQUIRE_FALSE(IS_SET(b_in_ld.data(), i));
        std::vector<std::string> mock_input(static_cast<uintptr_t>(n_col),
                                            std::string("@"));
        int temp = dis(gen);
        // generate non-duplicated samples
        while (dup_gen.find(temp) != dup_gen.end()) temp = dis(gen);
        dup_gen.insert(temp);
        // around 2% duplication?
        auto duplicated = is_duplicate(gen) < 2;

        mock_input[0] = "F" + std::to_string(temp);
        std::string exp = mock_input[0];
        if (!ignore_fid) { mock_input[1] = exp; }
        if (duplicated)
        { processed_samples.insert((ignore_fid ? "" : exp + delim) + exp); }
        included_samples.push_back(exp);
        // always include this info
        mock_input[3] = "Dad" + std::to_string(dis(gen));
        mock_input[4] = "Mum" + std::to_string(dis(gen));
        if (has_sex)
        {
            auto sex = sex_rand(gen);
            switch (sex)
            {
            case 1:
                mock_input[2] = "1";
                if (!duplicated) ++exp_male;
                break;
            case 2:
                mock_input[2] = "2";
                if (!duplicated) ++exp_female;
                break;
            default:
                mock_input[2] = "NA";
                if (!duplicated) ++exp_ambig;
                break;
            }
            sex_idx = 2;
        }
        else
        {
            // when there's no sex information, we assume ambiguous
            sex_idx = ~size_t(0);
            if (!duplicated) ++exp_ambig;
        }
        auto has_dad = has_father(gen) == 1;
        auto has_mum = has_mother(gen) == 1;
        if (has_dad)
        {
            founder_info.insert((ignore_fid ? "" : exp + delim)
                                + mock_input[3]);
        }
        if (has_mum)
        {
            founder_info.insert((ignore_fid ? "" : exp + delim)
                                + mock_input[4]);
        }
        // sample selection
        auto remove = remove_rand(gen);
        if (remove == 1)
        {
            auto id = ((ignore_fid) ? "" : exp + delim) + exp;
            geno.add_select_sample(id);
        }
        auto dad_idx = has_father_id ? 3 : ~size_t(0);
        auto mum_idx = has_mother_id ? 4 : ~size_t(0);
        auto dup_size = duplicated_sample_names.size();
        REQUIRE_NOTHROW(
            geno.test_gen_sample(fid_idx, iid_idx, sex_idx, dad_idx, mum_idx, i,
                                 founder_info, "", mock_input, sample_storage,
                                 processed_samples, duplicated_sample_names));
        REQUIRE(duplicated_sample_names.size() == dup_size + duplicated);
        auto in_prs = geno.calculate_prs();
        auto in_ld = geno.sample_for_ld();
        if (remove_samples ^ (remove == 1) && !duplicated)
        {
            REQUIRE(IS_SET(in_prs.data(), i));
            // never use non-founder for LD calculation
            if (!ignore_fid
                && ((has_dad && has_father_id) || (has_mum && has_mother_id)))
            { REQUIRE_FALSE(IS_SET(in_ld.data(), i)); }
            else
            {
                REQUIRE(IS_SET(in_ld.data(), i));
            }
            // always have IID
            REQUIRE(sample_storage.back().IID == exp);
            // only use non-founder for regression if asked to
            if (!ignore_fid && !keep_nonfounder
                && ((has_dad && has_father_id) || (has_mum && has_mother_id)))
            { REQUIRE_FALSE(sample_storage.back().in_regression); }
            else
            {
                REQUIRE(sample_storage.back().in_regression);
            }
            if (!ignore_fid) { REQUIRE(sample_storage.back().FID == exp); }
            else
            {
                REQUIRE(sample_storage.back().FID.empty());
            }
        }
        else
        {
            REQUIRE_FALSE(IS_SET(in_prs.data(), i));
            REQUIRE_FALSE(IS_SET(in_ld.data(), i));
        }
        REQUIRE(geno.num_male() == exp_male);
        REQUIRE(geno.num_female() == exp_female);
        REQUIRE(geno.num_ambig_sex() == exp_ambig);
    }
}
