#include "binarygen.hpp"
#include "catch.hpp"
#include "mock_binarygen.hpp"

TEST_CASE("SNP loading check")
{
    // different version and different compression
    SECTION("Single file check")
    {
        // we only do single file here first, as it is slightly easier
        auto compressed = GENERATE(genfile::bgen::e_ZlibCompression,
                                   genfile::bgen::e_NoCompression);
        auto number_of_individuals = 3u;
        // now try to make the binarygen object
        Reporter reporter("log", 60, true);
        GenoFile geno;
        geno.num_autosome = 2;
        geno.file_name = "load_snp,sample";
        Phenotype pheno;
        mock_binarygen bgen(geno, pheno, " ", &reporter);
        bgen.test_init_chr();
        size_t idx = 0;
        std::vector<SNP> input;
        input.push_back(SNP("SNP_1", 1, 742429, "A", "C", idx, 1));
        input.push_back(SNP("SNP_2", 1, 933331, "C", "T", idx, 1));
        input.push_back(SNP("SNP_3", 1, 933331, "C", "T", idx, 1));
        input.push_back(SNP("SNP_4", 1, 1008567, "G", "A", idx, 1));
        input.push_back(SNP("SNP_5", 1, 742429, "A", "C", idx, 1));
        for (auto snp : input) { bgen.manual_load_snp(snp); }

        // load SNP
        std::vector<IITree<size_t, size_t>> exclusion_region;
        std::string mismatch_name = "mismatch";
        std::unordered_set<std::string> duplicated_snps;
        std::unordered_set<std::string> processed_snps;
        std::vector<bool> retain_snp(5, false);
        bool chr_error = false;
        bool sex_error = false;

        SECTION("fail chromosome")
        {
            input[3].add_snp_info(SNP("SNP_4",
                                      std::numeric_limits<size_t>::max() - 1,
                                      1008567, "G", "A", 0, 0),
                                  false, false);
            SECTION("version 1.1")
            {
                auto version = genfile::bgen::e_Layout1;
                auto phased = genfile::ePerUnorderedGenotype;
                auto str = bgen.gen_mock_snp(input, number_of_individuals,
                                             phased, version, compressed);
                auto in_file = std::make_unique<std::istringstream>(str);
                bgen.load_context(*in_file);
                // 4 because one of the SNP is removed due to chrX
                REQUIRE(bgen.test_transverse_bgen_for_snp(
                            exclusion_region, "load_bgen_snp", 0,
                            std::move(in_file), duplicated_snps, processed_snps,
                            retain_snp, chr_error, sex_error, &bgen)
                        == 4);
                auto res = bgen.existed_snps();
                // we shouldn't have reshaped our vector here
                REQUIRE(res.size() == 5);
                REQUIRE(sex_error);
            }
            SECTION("version 1.2")
            {
                /*
                auto version = genfile::bgen::e_Layout2;
                auto phased = GENERATE(genfile::ePerUnorderedGenotype,
                                       genfile::ePerPhasedHaplotypePerAllele);
                auto str = bgen.gen_mock_snp(input, number_of_individuals,
                                             phased, version, compressed);
                auto in_file = std::make_unique<std::istringstream>(str);
                bgen.load_context(*in_file, 0);
                REQUIRE_THROWS(bgen.test_transverse_bgen_for_snp(
                    exclusion_region, "load_bgen_snp", 0, std::move(in_file),
                    duplicated_snps, processed_snps, retain_snp, chr_error,
                    sex_error, &bgen));
                auto res = bgen.existed_snps();
                // we shouldn't have reshaped our vector here
                REQUIRE(res.size() == 5);
                REQUIRE(sex_error);*/
            }
            // auto in_file = std::make_unique<std::istringstream>(str);
            // bgen.load_context(*in_file, 0);
            /*
                        REQUIRE_THROWS(bgen.test_transverse_bgen_for_snp(
                            exclusion_region, "load_bgen_snp", 0,
               std::move(in_file), duplicated_snps, processed_snps, retain_snp,
               chr_error, sex_error, &bgen)); auto res = bgen.existed_snps();
                        // we shouldn't have reshaped our vector here
                        REQUIRE(res.size() == 5);
                        REQUIRE(sex_error);*/
        }
        SECTION("mismatch SNPs") {}
        SECTION("valid input") {}
        SECTION("duplicated SNP") {}
        SECTION("full run") {}
        SECTION("full run with duplicated") {}
    }
}
