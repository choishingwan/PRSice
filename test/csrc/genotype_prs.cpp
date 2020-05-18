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


TEST_CASE("Read_PRS")
{
    Reporter reporter("log", 60, true);
    mockGenotype geno;
    geno.set_reporter(&reporter);
    // need to set the weight
    // need to initialize prs_list
    // need m_sample_ct
    /*
    auto model = GENERATE(MODEL::ADDITIVE, MODEL::DOMINANT, MODEL::RECESSIVE,
                          MODEL::HETEROZYGOUS);
*/
    auto model = MODEL::ADDITIVE;
    geno.set_weight(model);
    const size_t num_sample = 2000;
    std::vector<bool> selected_sample(num_sample);
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()};
    std::uniform_int_distribution<size_t> rand_selected {1, 10},
        geno_dist {0, 3};
    auto selected = [&rand_selected, &mersenne_engine]() {
        return rand_selected(mersenne_engine) > 3;
    };
    auto gen_geno = [&geno_dist, &mersenne_engine]() {
        return geno_dist(mersenne_engine);
    };
    std::generate(begin(selected_sample), end(selected_sample), selected);
    geno.set_sample_vector(selected_sample);
    const size_t num_selected = geno.num_sample();

    const uintptr_t sample_ctl = BITCT_TO_WORDCT(num_sample);
    const uintptr_t sample_ctv2 = 2 * sample_ctl;
    std::vector<uintptr_t> genotype_data(sample_ctv2, 0);
    auto miss_score = GENERATE(take(1, random(-1.0, 1.0)));
    auto stat = GENERATE(take(1, random(-1.0, 1.0)));
    const size_t ploidy = 2;
    auto adj_score = GENERATE(take(1, random(-1.0, 1.0)), 0.0);
    size_t homcom_weight = 0, het_weight = 1, homrar_weight = 2;
    switch (model)
    {
    case MODEL::HETEROZYGOUS:
        homcom_weight = 0;
        het_weight = 1;
        homrar_weight = 0;
        break;
    case MODEL::DOMINANT:
        homcom_weight = 0;
        het_weight = 1;
        homrar_weight = 1;
        break;
    case MODEL::RECESSIVE:
        homcom_weight = 0;
        het_weight = 0;
        homrar_weight = 1;
        break;
    default:
        homcom_weight = 0;
        het_weight = 1;
        homrar_weight = 2;
        break;
    }
    auto miss_count = GENERATE(0ul, 2ul);
    std::vector<double> expected_weighting = {
        homcom_weight * stat - adj_score, het_weight * stat - adj_score,
        miss_score, homrar_weight * stat - adj_score};
    std::vector<size_t> counts = {ploidy, ploidy, miss_count, ploidy};
    size_t geno_idx = 0;
    std::vector<double> expected_prs(num_selected, 0.0);
    std::vector<size_t> expected_num(num_selected, 0);
    for (size_t i = 0; i < num_selected; ++i)
    {
        auto geno = gen_geno();
        expected_prs[i] += expected_weighting[geno];
        expected_num[i] += counts[geno];

        // when parse, will do bitwise not
        switch (geno)
        {
        case 0:
            SET_BIT(geno_idx + 1, genotype_data.data());
            SET_BIT(geno_idx, genotype_data.data());
            break;
        case 1: SET_BIT(geno_idx + 1, genotype_data.data()); break;
        case 2: SET_BIT(geno_idx, genotype_data.data()); break;
        case 3: break;
        }
        geno_idx += 2;
    }
    std::vector<double> observed_prs(num_selected);
    std::vector<size_t> observed_num(num_selected);
    std::vector<PRS> observed(num_selected, PRS());
    const bool not_first = true;

    geno.test_read_prs(genotype_data.data(), observed, ploidy, stat, adj_score,
                       miss_score, miss_count, homcom_weight, het_weight,
                       homrar_weight, !not_first);
    for (size_t i = 0; i < num_selected; ++i)
    {
        observed_prs[i] = observed[i].prs;
        observed_num[i] = observed[i].num_snp;
    }

    REQUIRE_THAT(observed_prs, Catch::Equals<double>(expected_prs));
    REQUIRE_THAT(observed_num, Catch::Equals<size_t>(expected_num));
    SECTION("Reset value")
    {
        geno.test_read_prs(genotype_data.data(), observed, ploidy, stat,
                           adj_score, miss_score, miss_count, homcom_weight,
                           het_weight, homrar_weight, !not_first);
        for (size_t i = 0; i < num_selected; ++i)
        {
            observed_prs[i] = observed[i].prs;
            observed_num[i] = observed[i].num_snp;
        }
        REQUIRE_THAT(observed_prs, Catch::Equals<double>(expected_prs));
        REQUIRE_THAT(observed_num, Catch::Equals<size_t>(expected_num));
    }
    SECTION("Don't reset value")
    {
        // we expect value to double
        geno.test_read_prs(genotype_data.data(), observed, ploidy, stat,
                           adj_score, miss_score, miss_count, homcom_weight,
                           het_weight, homrar_weight, not_first);
        for (size_t i = 0; i < num_selected; ++i)
        {
            observed_prs[i] = observed[i].prs;
            observed_num[i] = observed[i].num_snp;
            expected_prs[i] = 2 * expected_prs[i];
            expected_num[i] = 2 * expected_num[i];
        }
        REQUIRE_THAT(observed_prs, Catch::Equals<double>(expected_prs));
        REQUIRE_THAT(observed_num, Catch::Equals<size_t>(expected_num));
    }
}
