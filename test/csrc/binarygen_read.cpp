#include "catch.hpp"
#include "mock_binarygen.hpp"


TEST_CASE("bgen read_genotype")
{
    // different founder status (usually we don't have founder in bgen, but
    // whatever)
    using record = std::tuple<genfile::bgen::Layout, genfile::OrderType>;
    // different layout
    auto settings = GENERATE(table<genfile::bgen::Layout, genfile::OrderType>(
        {record {genfile::bgen::e_Layout1, genfile::ePerUnorderedGenotype},
         record {genfile::bgen::e_Layout2,
                 genfile::ePerPhasedHaplotypePerAllele},
         record {genfile::bgen::e_Layout2, genfile::ePerUnorderedGenotype}}));
    // compressed or not
    auto compressed = GENERATE(genfile::bgen::e_ZlibCompression,
                               genfile::bgen::e_NoCompression);
    auto n_entries =
        std::get<1>(settings) == genfile::ePerUnorderedGenotype ? 3ul : 4ul;
    // just use default
    QCFiltering qc;
    // different sample size
    auto n_target = 125ul, n_ref = 257ul;
    std::vector<uintptr_t> expected_geno, expected_target, tmp2, expected_mem,
        tmp_mem;
    std::vector<bool> ref_founder, target_founder;
    target_founder.resize(n_target, true);
    ref_founder.resize(n_ref, true);
    std::vector<double> target_prob(n_target * n_entries);
    mock_binarygen::generate_samples(
        n_entries, n_target, qc, expected_target, expected_mem, target_prob,
        target_founder, std::get<0>(settings), std::get<1>(settings));
    std::vector<double> ref_prob(n_ref * n_entries);
    std::vector<double> second_prob(n_ref * n_entries);
    mock_binarygen::generate_samples(
        n_entries, n_ref, qc, tmp2, tmp_mem, ref_prob, ref_founder,
        std::get<0>(settings), std::get<1>(settings));
    // we want read_genotype to read the second SNP just in case
    mock_binarygen::generate_samples(
        n_entries, n_ref, qc, expected_geno, tmp_mem, second_prob, ref_founder,
        std::get<0>(settings), std::get<1>(settings));
    Reporter reporter("log", 60, true);
    Phenotype pheno;
    GenoFile target_geno_info, ref_geno_info;
    target_geno_info.num_autosome = 22;
    target_geno_info.file_name = "target_read,sample";
    ref_geno_info.num_autosome = 22;
    ref_geno_info.file_name = "ref_read,sample";
    mock_binarygen target_bgen(target_geno_info, pheno, " ", &reporter),
        ref_bgen(ref_geno_info, pheno, " ", &reporter);
    ref_bgen.reference();
    target_bgen.expect_reference();
    target_bgen.test_init_chr();
    ref_bgen.test_init_chr();
    target_bgen.set_thresholds(qc);
    ref_bgen.set_thresholds(qc);
    target_bgen.update_sample(n_target);
    ref_bgen.update_sample(n_ref);
    auto allow_inter = GENERATE(true, false);
    target_bgen.intermediate(allow_inter);
    ref_bgen.intermediate(allow_inter);
    size_t idx = 0;
    std::vector<SNP> input = {SNP("SNP_1", 1, 12382, "A", "C", idx, 1)};
    std::vector<SNP> ref_input = {SNP("SNP_2", 1, 123, "C", "G", idx, 1),
                                  input.front()};
    auto target_str = target_bgen.gen_mock_snp(
        std::vector<std::vector<double>> {target_prob}, input, n_target,
        std::get<1>(settings), std::get<0>(settings), compressed);
    auto ref_str = ref_bgen.gen_mock_snp(
        std::vector<std::vector<double>> {ref_prob, second_prob}, ref_input,
        n_ref, std::get<1>(settings), std::get<0>(settings), compressed);
    // ref_input will have the new byte pos set to the target location and we
    // want to change that
    assert(input.size() == 1);
    input.front().update_file(ref_input.back().get_file_idx(),
                              ref_input.back().get_byte_pos(), true);
    auto target_file = std::make_unique<std::istringstream>(target_str);
    auto ref_file = std::make_unique<std::istringstream>(ref_str);
    target_bgen.manual_load_snp(input.front());
    target_bgen.load_context(*target_file);
    ref_bgen.load_context(*ref_file);

    std::ofstream t_out("target_read.bgen", std::ios::binary),
        r_out("ref_read.bgen", std::ios::binary);
    t_out << target_str << std::endl;
    t_out.close();
    r_out << ref_str << std::endl;
    r_out.close();
    auto hard_coded = GENERATE(true, false);
    target_bgen.set_hard_code(hard_coded);
    ref_bgen.set_hard_code(hard_coded);
    target_bgen.calc_freqs_and_intermediate(qc, "bgen_read", true);
    ref_bgen.calc_freqs_and_intermediate(qc, "bgen_read", true, &target_bgen);
    SECTION("Reading without storing")
    {
        if (allow_inter)
        {
            REQUIRE(ref_bgen.genotype_file_names().back() == "bgen_read.inter");
        }
        const uintptr_t ref_sample_ctl = BITCT_TO_WORDCT(n_ref);
        const uintptr_t ref_sample_ctv2 = 2 * ref_sample_ctl;
        std::vector<uintptr_t> observed_geno(ref_sample_ctv2, 0);
        const uintptr_t target_sample_ctl = BITCT_TO_WORDCT(n_target);
        const uintptr_t target_sample_ctv2 = 2 * target_sample_ctl;
        std::vector<uintptr_t> observed_target_geno(target_sample_ctv2, 0);
        auto target_snps = target_bgen.existed_snps();
        ref_bgen.test_read_genotype(target_snps.back(), n_ref,
                                    observed_geno.data(), true);
        target_bgen.test_read_genotype(target_snps.back(), n_target,
                                       observed_target_geno.data(), false);
        REQUIRE_THAT(observed_geno, Catch::Equals<uintptr_t>(expected_geno));
        REQUIRE_THAT(observed_target_geno,
                     Catch::Equals<uintptr_t>(expected_target));
    }

    SECTION("Read into memory")
    {
        target_bgen.load_genotype_to_memory();
        auto&& snp = target_bgen.existed_snps();
        if (!hard_coded) { REQUIRE(snp.front().current_genotype() == nullptr); }
        else
        {
            const uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(n_target);
            const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
            auto&& snp_geno = snp.front().current_genotype();
            REQUIRE(snp_geno != nullptr);
            std::vector<uintptr_t> observed(snp_geno,
                                            snp_geno + unfiltered_sample_ctv2);
            REQUIRE_THAT(observed, Catch::Equals<uintptr_t>(expected_mem));
        }
    }
}
