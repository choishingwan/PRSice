#include "binarygen.hpp"
#include "catch.hpp"
#include "mock_binarygen.hpp"

TEST_CASE("SNP loading check")
{
    // different version and different compression
    using record = std::tuple<genfile::bgen::Layout, genfile::OrderType>;
    auto settings = GENERATE(table<genfile::bgen::Layout, genfile::OrderType>(
        {record {genfile::bgen::e_Layout1, genfile::ePerUnorderedGenotype},
         record {genfile::bgen::e_Layout2, genfile::ePerUnorderedGenotype},
         record {genfile::bgen::e_Layout2,
                 genfile::ePerPhasedHaplotypePerAllele}}));
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
        SECTION("filtered by chromosome")
        {
            input[3].add_snp_info(SNP("SNP_4",
                                      std::numeric_limits<size_t>::max() - 1,
                                      1008567, "G", "A", 0, 0),
                                  false, false);
            auto str = bgen.gen_mock_snp(input, number_of_individuals,
                                         std::get<1>(settings),
                                         std::get<0>(settings), compressed);
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
        SECTION("mismatch SNPs")
        {
            input[3].add_snp_info(SNP("SNP_4", 1, 1008568, "G", "A", 0, 0),
                                  false, false);
            auto str = bgen.gen_mock_snp(input, number_of_individuals,
                                         std::get<1>(settings),
                                         std::get<0>(settings), compressed);
            auto in_file = std::make_unique<std::istringstream>(str);
            bgen.load_context(*in_file);
            // 4 because one of the SNP is removed due to chrX
            REQUIRE(bgen.test_transverse_bgen_for_snp(
                        exclusion_region, "load_bgen_snp", 0,
                        std::move(in_file), duplicated_snps, processed_snps,
                        retain_snp, chr_error, sex_error, &bgen)
                    == 4);
            std::ifstream mismatch("load_bgen_snp");
            REQUIRE(mismatch.is_open());
            std::string line;
            size_t num_mismatch = 0;
            std::vector<std::string> token;
            while (std::getline(mismatch, line))
            {
                token = misc::split(line, "\t");
                ++num_mismatch;
            }
            // because of header
            REQUIRE(num_mismatch == 2);
            REQUIRE(token.size() == 10);
            REQUIRE(token[1] == "SNP_4");
        }
        SECTION("duplicated SNP")
        {
            input.push_back(SNP("SNP_4", 1, 1008567, "G", "A", idx, 1));
            auto str = bgen.gen_mock_snp(input, number_of_individuals,
                                         std::get<1>(settings),
                                         std::get<0>(settings), compressed);
            auto in_file = std::make_unique<std::istringstream>(str);
            bgen.load_context(*in_file);
            // 4 because one of the SNP is removed due to chrX
            REQUIRE(bgen.test_transverse_bgen_for_snp(
                        exclusion_region, "load_bgen_snp", 0,
                        std::move(in_file), duplicated_snps, processed_snps,
                        retain_snp, chr_error, sex_error, &bgen)
                    == 5);
            REQUIRE(duplicated_snps.size() == 1);
            REQUIRE(duplicated_snps.find("SNP_4") != duplicated_snps.end());
        }
        SECTION("valid input")
        {
            input.erase(input.begin() + 2);
            auto str = bgen.gen_mock_snp(input, number_of_individuals,
                                         std::get<1>(settings),
                                         std::get<0>(settings), compressed);
            auto in_file = std::make_unique<std::istringstream>(str);
            bgen.load_context(*in_file);
            // 4 as one of the SNPs were not found in target
            REQUIRE(bgen.test_transverse_bgen_for_snp(
                        exclusion_region, "load_bgen_snp", 0,
                        std::move(in_file), duplicated_snps, processed_snps,
                        retain_snp, chr_error, sex_error, &bgen)
                    == 4);
            auto res = bgen.existed_snps();
            REQUIRE(res.size() == 5);
            REQUIRE(res[0].rs() == "SNP_1");
            REQUIRE(res[0].get_byte_pos() != 0);
            REQUIRE(res[1].rs() == "SNP_2");
            REQUIRE(res[1].get_byte_pos() != 0);
            REQUIRE(res[2].rs() == "SNP_3");
            REQUIRE(res[2].get_byte_pos() != 0);
            REQUIRE(res[3].rs() == "SNP_4");
            REQUIRE(res[3].get_byte_pos() != 0);
            REQUIRE(res[4].rs() == "SNP_5");
            REQUIRE(res[4].get_byte_pos() != 0);
        }
        SECTION("sample identifier check") {}
    }
    SECTION("full run")
    {
        auto compressed = GENERATE(genfile::bgen::e_ZlibCompression,
                                   genfile::bgen::e_NoCompression);
        auto number_of_individuals = 10u;
        // now try to make the binarygen object
        Reporter reporter("log", 60, true);
        GenoFile geno;
        geno.num_autosome = 2;
        geno.file_name = "load_bgen#,sample";
        Phenotype pheno;
        mock_binarygen bgen(geno, pheno, " ", &reporter);
        bgen.test_init_chr();
        size_t idx = 0;
        std::vector<SNP> input;
        input.push_back(SNP("SNP_1", 1, 742429, "A", "C", idx, 1));
        input.push_back(SNP("SNP_2", 1, 933331, "C", "T", idx, 1));
        input.push_back(SNP("SNP_3", 1, 933331, "C", "T", idx, 1));
        auto first = bgen.gen_mock_snp(input, number_of_individuals,
                                       std::get<1>(settings),
                                       std::get<0>(settings), compressed);
        std::ofstream first_file("load_bgen1.bgen");
        first_file << first << std::endl;
        first_file.close();
        auto first_context =
            misc::load_stream("load_bgen1.bgen", std::ios_base::binary);
        bgen.load_context(*first_context, 0);
        first_context.reset();
        for (auto snp : input) { bgen.manual_load_snp(snp); }
        // now load the other two sample without printing out the file yet
        input.clear();
        input.push_back(SNP("SNP_4", 1, 1008567, "G", "A", idx, 1));
        input.push_back(SNP("SNP_5", 1, 742429, "A", "C", idx, 1));
        // load SNP
        std::vector<IITree<size_t, size_t>> exclusion_region;
        // make two file
        SECTION("without duplicate")
        {
            input.push_back(SNP("SNP_6", 1, 1008512, "G", "A", idx, 1));
            input.push_back(SNP("SNP_7", 1, 1008876, "C", "A", idx, 1));
            // ignore SNP 4 and 6
            for (auto snp : input) { bgen.manual_load_snp(snp); }
            input.erase(input.begin() + 2);
            input.erase(input.begin());
            auto second = bgen.gen_mock_snp(input, number_of_individuals,
                                            std::get<1>(settings),
                                            std::get<0>(settings), compressed);
            std::ofstream sec_file("load_bgen2.bgen");
            sec_file << second << std::endl;
            sec_file.close();
            auto sec_context =
                misc::load_stream("load_bgen2.bgen", std::ios_base::binary);
            bgen.load_context(*sec_context, 1);
            sec_context.reset();
            REQUIRE_NOTHROW(
                bgen.load_snps("load_bgen", exclusion_region, false));
            auto res = bgen.existed_snps();
            REQUIRE(res.size() == 5);
            std::vector<std::string> expected = {"SNP_1", "SNP_2", "SNP_3",
                                                 "SNP_5", "SNP_7"};
            for (size_t i = 0; i < res.size(); ++i)
            { REQUIRE(res[i].rs() == expected[i]); }
        }
        SECTION("with duplicate")
        {
            input.push_back(SNP("SNP_6", 1, 1008512, "G", "A", idx, 1));
            input.push_back(SNP("SNP_7", 1, 1008876, "C", "A", idx, 1));
            // ignore SNP 4 and 6
            for (auto snp : input) { bgen.manual_load_snp(snp); }
            input.erase(input.begin() + 2);
            input.erase(input.begin());
            // duplicate 7
            input.push_back(SNP("SNP_7", 1, 1008876, "C", "A", idx, 1));
            auto second = bgen.gen_mock_snp(input, number_of_individuals,
                                            std::get<1>(settings),
                                            std::get<0>(settings), compressed);
            std::ofstream sec_file("load_bgen2.bgen");
            sec_file << second << std::endl;
            sec_file.close();
            auto sec_context =
                misc::load_stream("load_bgen2.bgen", std::ios_base::binary);
            bgen.load_context(*sec_context, 1);
            sec_context.reset();
            REQUIRE_THROWS(
                bgen.load_snps("load_bgen", exclusion_region, false));
            std::ifstream dup("load_bgen.valid");
            std::string line;
            size_t num_valid = 0;
            std::vector<std::string> token;
            while (std::getline(dup, line))
            {
                token = misc::split(line, "\t");
                ++num_valid;
            }
            REQUIRE(num_valid == 4);
            REQUIRE(token.size() == 5);
            REQUIRE(token[0] == "SNP_5");
        }
    }
}
