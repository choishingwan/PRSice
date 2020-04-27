#include "binarygen.hpp"
#include "binarygen_setters.hpp"
#include "catch.hpp"
#include "mock_binarygen.hpp"

void generate_samples(const size_t n_entries, const size_t n_sample,
                      const QCFiltering qc,
                      std::vector<uintptr_t>& plink_genotype,
                      std::vector<double>& data_prob, double& exp_mach,
                      double& exp_impute, uint32_t& ref_ct, uint32_t& het_ct,
                      uint32_t& alt_ct, uint32_t& miss_ct)
{
    miss_ct = 0;
    ref_ct = 0;
    het_ct = 0;
    alt_ct = 0;
    misc::RunningStat statistic;
    std::vector<size_t> expected_geno(n_sample, 0);
    std::vector<double> cur_geno_prob(3, 0.0);
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()}; // Generates random integers
    std::uniform_real_distribution<double> dist {0, 1};
    std::uniform_int_distribution<size_t> rand_miss {1, 10};
    // round to 4 decimal places so that BGEN's imprecision would not cause our
    // test to break
    auto gen_prob = [&dist, &mersenne_engine]() {
        return dist(mersenne_engine);
        // double tmp = dist(mersenne_engine);
        // return std::round(tmp * 100000.0) / 100000.0;
    };
    auto missing_prob = [&rand_miss, &mersenne_engine]() {
        return rand_miss(mersenne_engine) <= 2;
    };
    bool missing = false;
    exp_mach = 0.0;
    exp_impute = 0.0;
    double sum_prob = 0.0, impute_sum = 0.0;
    for (size_t i = 0; i < n_sample; ++i)
    {
        double exp = 0.0, tmp = 0.0;
        expected_geno[i] = 3;
        if (missing_prob())
        {
            missing = true;
            ++miss_ct;
            SET_BIT(i * 2, plink_genotype.data());
        }
        else
        {
            missing = false;
            double prob_a1 = gen_prob();
            if (n_entries == 3)
            {
                data_prob[i * n_entries] = prob_a1 * prob_a1;
                data_prob[i * n_entries + 1] = prob_a1 * (1.0 - prob_a1) * 2;
                data_prob[i * n_entries + 2] =
                    (1.0 - prob_a1) * (1.0 - prob_a1);
                cur_geno_prob[0] = prob_a1 * prob_a1;
                cur_geno_prob[1] = 2 * prob_a1 * (1.0 - prob_a1);
                cur_geno_prob[2] = (1.0 - prob_a1) * (1.0 - prob_a1);
            }
            else
            {
                double prob_a2 = gen_prob();
                data_prob[i * n_entries] = prob_a1;
                data_prob[i * n_entries + 1] = 1.0 - prob_a1;
                data_prob[i * n_entries + 2] = prob_a2;
                data_prob[i * n_entries + 3] = 1.0 - prob_a2;
                cur_geno_prob[0] = prob_a1 * prob_a2;
                cur_geno_prob[1] =
                    prob_a1 * (1.0 - prob_a2) + (1.0 - prob_a1) * prob_a2;
                cur_geno_prob[2] = (1.0 - prob_a1) * (1.0 - prob_a2);
            }
            exp = cur_geno_prob[1] + cur_geno_prob[2] * 2.0;
            tmp = cur_geno_prob[1] + cur_geno_prob[2] * 4.0;
            for (auto&& g : cur_geno_prob) { sum_prob += g; }
            const double prob1 = cur_geno_prob[0] * 2 + cur_geno_prob[1];
            const double prob2 = cur_geno_prob[2] * 2 + cur_geno_prob[1];
            const double hard_score = (std::fabs(prob1 - std::round(prob1))
                                       + std::fabs(prob2 - std::round(prob2)))
                                      * 0.5;
            if (hard_score <= qc.hard_threshold)
            {
                double hard_prob = 0;
                for (size_t geno = 0; geno < 3; ++geno)
                {
                    if (cur_geno_prob[geno] > hard_prob
                        && cur_geno_prob[geno] >= qc.dose_threshold)
                    {
                        // +1 because geno ==1 represents missing
                        expected_geno[i] = geno;
                        hard_prob = cur_geno_prob[geno];
                    }
                }
            }


            if (!missing)
            {
                statistic.push(exp);
                impute_sum -= tmp - exp * exp;
            }
            switch (expected_geno[i])
            {
            case 0: ++ref_ct; break;
            case 1:
                ++het_ct;
                SET_BIT(i * 2 + 1, plink_genotype.data());
                break;
            case 2:
                ++alt_ct;
                SET_BIT(i * 2, plink_genotype.data());
                SET_BIT(i * 2 + 1, plink_genotype.data());
                break;
            case 3:
                ++miss_ct;
                SET_BIT(i * 2, plink_genotype.data());
                break;
            }
        }
    }
    double p = statistic.mean() / 2.0;
    double p_all = 2.0 * p * (1.0 - p);
    exp_mach = statistic.var() / p_all;
    exp_impute = 1.0 + ((impute_sum / p_all) / sum_prob);
}
TEST_CASE("BGEN Target filtering")
{
    auto n_sample = GENERATE(range(124u, 127u), range(124u, 127u));
    using record = std::tuple<genfile::bgen::Layout, genfile::OrderType>;
    auto settings = GENERATE(table<genfile::bgen::Layout, genfile::OrderType>(
        {record {genfile::bgen::e_Layout1, genfile::ePerUnorderedGenotype},
         record {genfile::bgen::e_Layout2,
                 genfile::ePerPhasedHaplotypePerAllele},
         record {genfile::bgen::e_Layout2, genfile::ePerUnorderedGenotype}}));

    auto compressed = GENERATE(genfile::bgen::e_ZlibCompression,
                               genfile::bgen::e_NoCompression);
    auto n_entries =
        std::get<1>(settings) == genfile::ePerUnorderedGenotype ? 3ul : 4ul;
    std::vector<double> data_prob(n_sample * n_entries);
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()}; // Generates random integers
    std::uniform_real_distribution<double> dist {0.7, 1};
    QCFiltering qc;
    qc.hard_threshold = dist(mersenne_engine);
    qc.dose_threshold = dist(mersenne_engine);
    const uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(n_sample);
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    std::vector<uintptr_t> plink_genotype(unfiltered_sample_ctv2, 0);
    uint32_t ref_ct = 0, het_ct = 0, alt_ct = 0, miss_ct = 0;
    double exp_impute = 0.0, exp_mach = 0.0;
    generate_samples(n_entries, n_sample, qc, plink_genotype, data_prob,
                     exp_mach, exp_impute, ref_ct, het_ct, alt_ct, miss_ct);
    Reporter reporter("log", 60, true);
    GenoFile geno;
    geno.num_autosome = 2;
    geno.file_name = "filtering,sample";
    Phenotype pheno;
    mock_binarygen bgen(geno, pheno, " ", &reporter);
    bgen.test_init_chr();
    size_t idx = 0;
    std::vector<SNP> input = {SNP("SNP_1", 1, 12382, "A", "C", idx, 1)};
    bgen.set_thresholds(qc);
    bgen.update_sample(n_sample);
    SECTION("Directly testing the setter")
    {
        auto str = bgen.gen_mock_snp(
            std::vector<std::vector<double>> {data_prob}, input, n_sample,
            std::get<1>(settings), std::get<0>(settings), compressed);
        auto in_file = std::make_unique<std::istringstream>(str);
        bgen.load_context(*in_file);
        auto [ct, impute, mach] = bgen.test_plink_generator(
            std::move(in_file), input.front().get_byte_pos());
        REQUIRE(ct.homcom == ref_ct);
        REQUIRE(ct.het == het_ct);
        REQUIRE(ct.homrar == alt_ct);
        REQUIRE(ct.missing == miss_ct);
        // BGEN only has accuracy up to four decimal places.
        REQUIRE(impute == Approx(exp_impute).epsilon(1e-4));
        REQUIRE(mach == Approx(exp_mach).epsilon(1e-4));
    }

    SECTION("Test filter from file")
    {
        auto str = bgen.gen_mock_snp(
            std::vector<std::vector<double>> {data_prob}, input, n_sample,
            std::get<1>(settings), std::get<0>(settings), compressed);
        auto in_file = std::make_unique<std::istringstream>(str);
        bgen.manual_load_snp(input.front());
        bgen.load_context(*in_file);
        std::ofstream file("filtering.bgen", std::ios::binary);
        file << str << std::endl;
        file.close();
        qc.info_type = GENERATE(INFO::IMPUTE2, INFO::MACH);
        qc.info_score = GENERATE(take(2, random(0.0, 1.0)));
        std::uniform_real_distribution<double> threshold {0.0, 1};
        qc.geno = threshold(mersenne_engine);
        qc.maf = threshold(mersenne_engine);
        SECTION("Check Intermediate file generation")
        {
            bgen.intermediate(true);
            REQUIRE(bgen.test_calc_freq_gen_inter(qc, "filter"));
            // we expect a file called filter.inter is generated
            std::ifstream intermediate("filter.inter", std::ios::binary);
            REQUIRE(intermediate.is_open());
            auto names = bgen.genotype_file_names();
            REQUIRE(names.size() == 2);
            REQUIRE(names.back() == "filter.inter");
            double exp_maf =
                (2.0 * alt_ct + het_ct) / (2.0 * (ref_ct + het_ct + alt_ct));
            if (exp_maf > 0.5) exp_maf = 1 - exp_maf;
            if (miss_ct == n_sample || miss_ct / n_sample > qc.geno
                || exp_maf < qc.maf
                || (qc.info_type == INFO::IMPUTE2 && qc.info_score > exp_impute)
                || (qc.info_type == INFO::MACH && qc.info_score > exp_mach))
            {}
            else
            {
                const uintptr_t unfiltered_sample_ct4 = (n_sample + 3) / 4;
                std::vector<uintptr_t> observed(unfiltered_sample_ctv2, 0);
                REQUIRE_FALSE(
                    !intermediate.read(reinterpret_cast<char*>(observed.data()),
                                       unfiltered_sample_ct4));
                for (size_t i = 0; i < unfiltered_sample_ctv2; ++i)
                { REQUIRE(observed[i] == plink_genotype[i]); }
            }
        }
        SECTION("Without intermediate")
        {
            REQUIRE(bgen.test_calc_freq_gen_inter(qc, "filter"));
            double exp_maf =
                (2.0 * alt_ct + het_ct) / (2.0 * (ref_ct + het_ct + alt_ct));
            if (exp_maf > 0.5) exp_maf = 1 - exp_maf;
            if (miss_ct == n_sample)
            {
                REQUIRE(bgen.num_miss_filter() == 1);
                REQUIRE(bgen.num_geno_filter() == 0);
                REQUIRE(bgen.num_maf_filter() == 0);
                REQUIRE(bgen.num_info_filter() == 0);
            }
            else if (miss_ct / n_sample > qc.geno)
            {
                REQUIRE(bgen.num_miss_filter() == 0);
                REQUIRE(bgen.num_maf_filter() == 0);
                REQUIRE(bgen.num_info_filter() == 0);
                REQUIRE(bgen.num_geno_filter() == 1);
            }
            else if (exp_maf < qc.maf)
            {
                REQUIRE(bgen.num_miss_filter() == 0);
                REQUIRE(bgen.num_info_filter() == 0);
                REQUIRE(bgen.num_geno_filter() == 0);
                REQUIRE(bgen.num_maf_filter() == 1);
            }
            else if ((qc.info_type == INFO::IMPUTE2
                      && qc.info_score > exp_impute)
                     || (qc.info_type == INFO::MACH
                         && qc.info_score > exp_mach))
            {
                REQUIRE(bgen.num_miss_filter() == 0);
                REQUIRE(bgen.num_geno_filter() == 0);
                REQUIRE(bgen.num_maf_filter() == 0);
                REQUIRE(bgen.num_info_filter() == 1);
            }
            else
            {
                REQUIRE(bgen.num_miss_filter() == 0);
                REQUIRE(bgen.num_geno_filter() == 0);
                REQUIRE(bgen.num_maf_filter() == 0);
                REQUIRE(bgen.num_info_filter() == 0);
            }
        }
    }
}
