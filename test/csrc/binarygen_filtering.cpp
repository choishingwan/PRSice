#include "binarygen.hpp"
#include "binarygen_setters.hpp"
#include "catch.hpp"
#include "mock_binarygen.hpp"

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
    std::vector<bool> founder(n_sample, true);
    mock_binarygen::generate_samples(
        n_entries, n_sample, qc, plink_genotype, data_prob, founder,
        std::get<0>(settings), std::get<1>(settings), exp_mach, exp_impute,
        ref_ct, het_ct, alt_ct, miss_ct);
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
        REQUIRE(impute == Approx(exp_impute));
        REQUIRE(mach == Approx(exp_mach));
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
        SECTION("Test using genotype function")
        {
            bgen.calc_freqs_and_intermediate(qc, "test", false);
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
        SECTION("Using bgen function Without intermediate")
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
