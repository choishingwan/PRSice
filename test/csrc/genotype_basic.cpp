#include "catch.hpp"
#include "genotype.hpp"
#include "mock_genotype.hpp"
#include "plink_common.hpp"
#include "reporter.hpp"
#include <algorithm>
#include <random>

TEST_CASE("CHR_CONVERTION")
{
    mockGenotype geno;
    SECTION("CHECK has CHR as prefix")
    {
        SECTION("false") { REQUIRE_FALSE(geno.chr_prefix("1")); }
        SECTION("true") { REQUIRE(geno.chr_prefix("chr1")); }
    }
    SECTION("Check chr code parsing")
    {
        SECTION("normal")
        {
            auto i = GENERATE(range(1, 22));
            auto chr = GENERATE("chr", "ChR", "CHR", "chR", "Chr");
            std::string code = std::string(chr) + std::to_string(i);
            REQUIRE(geno.get_chrom_code(code) == i);
        }
        SECTION("Special")
        {
            using record = std::tuple<std::string, size_t>;
            auto i = GENERATE(table<std::string, int32_t>(
                {record {"X", CHROM_X}, record {"Y", CHROM_Y},
                 record {"XY", CHROM_XY}, record {"M", CHROM_MT},
                 record {"mt", CHROM_MT}, record {"0X", CHROM_X},
                 record {"0Y", CHROM_Y}, record {"0M", CHROM_MT}}));
            auto chr = GENERATE("", "chr");
            auto code = std::string(chr) + std::string(std::get<0>(i));
            REQUIRE(geno.get_chrom_code(code) == std::get<1>(i));
        }
        SECTION("fail")
        {
            auto i = GENERATE("CHR", "Happy", "chromosome");
            REQUIRE(geno.get_chrom_code(i) == -1);
        }
    }
}

TEST_CASE("Initialize genotype")
{
    mockGenotype genotype;
    Reporter reporter("log", 60, true);
    Phenotype pheno;
    GenoFile geno;
    auto delim = " ";
    SECTION("non-list input")
    {
        auto type = GENERATE("bed", "bgen");
        auto ignore_fid = GENERATE(true, false);
        auto keep = GENERATE("", "keep");
        auto remove = GENERATE("", "remove");
        geno.keep = keep;
        geno.remove = remove;
        pheno.ignore_fid = ignore_fid;
        SECTION("basic initialization")
        {
            geno.file_name = "Genotype";
            genotype.test_initialize(geno, pheno, delim, type, &reporter);
            REQUIRE(delim == genotype.delim());
            REQUIRE(genotype.ignore_fid() == ignore_fid);
            REQUIRE(genotype.keep_file() == keep);
            REQUIRE(genotype.remove_file() == remove);
            REQUIRE(genotype.sample_file().empty());
            REQUIRE_THAT(genotype.genotype_file_names(),
                         Catch::Equals<std::string>({"Genotype"}));
        }
        SECTION("intialize with external sample")
        {
            geno.file_name = "Genotype,external";
            genotype.test_initialize(geno, pheno, delim, type, &reporter);
            REQUIRE(delim == genotype.delim());
            REQUIRE(genotype.ignore_fid() == ignore_fid);
            REQUIRE(genotype.keep_file() == keep);
            REQUIRE(genotype.remove_file() == remove);
            REQUIRE(genotype.sample_file() == "external");
            REQUIRE_THAT(genotype.genotype_file_names(),
                         Catch::Equals<std::string>({"Genotype"}));
        }
        SECTION("chromosome number replacement")
        {
            auto num_auto = GENERATE(take(10, random(1, 40)));
            geno.num_autosome = num_auto;
            geno.file_name = "chr#_geno";
            std::vector<std::string> expected;
            for (int i = 0; i < num_auto; ++i)
            { expected.push_back("chr" + std::to_string(i + 1) + "_geno"); }
            genotype.test_initialize(geno, pheno, delim, type, &reporter);
            REQUIRE(delim == genotype.delim());
            REQUIRE(genotype.ignore_fid() == ignore_fid);
            REQUIRE(genotype.keep_file() == keep);
            REQUIRE(genotype.remove_file() == remove);
            REQUIRE(genotype.sample_file().empty());
            REQUIRE_THAT(genotype.genotype_file_names(),
                         Catch::Equals<std::string>(expected));
        }
        SECTION("chromosome number multiple replacement")
        {
            auto num_auto = GENERATE(take(10, random(1, 40)));
            geno.num_autosome = num_auto;
            geno.file_name = "chr#_geno#";
            std::vector<std::string> expected;
            for (int i = 0; i < num_auto; ++i)
            {
                expected.push_back("chr" + std::to_string(i + 1) + "_geno"
                                   + std::to_string(i + 1));
            }
            genotype.test_initialize(geno, pheno, delim, type, &reporter);
            REQUIRE(delim == genotype.delim());
            REQUIRE(genotype.ignore_fid() == ignore_fid);
            REQUIRE(genotype.keep_file() == keep);
            REQUIRE(genotype.remove_file() == remove);
            REQUIRE(genotype.sample_file().empty());
            REQUIRE_THAT(genotype.genotype_file_names(),
                         Catch::Equals<std::string>(expected));
        }
        SECTION("invalid formats")
        {
            auto code = GENERATE("a,b,c", ",");
            geno.file_name = code;
            REQUIRE_THROWS(
                genotype.test_initialize(geno, pheno, delim, type, &reporter));
        }
    }
    SECTION("Genotype list provided")
    {
        SECTION("Test load genotype prefix for list input")
        {
            std::unique_ptr<std::istringstream> input =
                std::make_unique<std::istringstream>();
            input->str("File_1\nFile_2");
            auto res = genotype.test_load_genotype_prefix(std::move(input));
            REQUIRE_THAT(res, Catch::Equals<std::string>({"File_1", "File_2"}));
        }
        SECTION("Check file name is properly extracted with external sample")
        {
            auto type = GENERATE("bgen", "bed");
            geno.file_list = "List_input,External";
            REQUIRE_THROWS_WITH(
                genotype.test_initialize(geno, pheno, delim, type, &reporter),
                "Error: Cannot open file: List_input");
        }
    }
}

TEST_CASE("Check ambiguous")
{
    using record = std::tuple<std::string, std::string, bool>;
    auto settings = GENERATE(table<std::string, std::string, bool>(
        {record {"A", "A", true}, record {"C", "C", true},
         record {"G", "G", true}, record {"T", "T", true},
         record {"A", "T", true}, record {"C", "G", true},
         record {"A", "", false}, record {"C", "", false},
         record {"G", "", false}, record {"T", "", false},
         record {"A", "C", false}, record {"A", "G", false},
         record {"T", "C", false}, record {"T", "G", false},
         record {"A", "AT", false}, record {"A", "AA", false}}));
    mockGenotype geno;
    REQUIRE(geno.test_ambiguous(std::get<0>(settings), std::get<1>(settings))
            == std::get<2>(settings));
    REQUIRE(geno.test_ambiguous(std::get<1>(settings), std::get<0>(settings))
            == std::get<2>(settings));
}

TEST_CASE("load snp selection")
{
    mockGenotype geno;
    Reporter reporter("log", 60, true);
    geno.set_reporter(&reporter);
    SECTION("Check rs column location")
    {
        SECTION("single column file or unidentified file")
        {
            auto i = GENERATE("rs1234", "rs1234 A B C D T E A");
            REQUIRE(geno.test_get_rs_column(i) == 0);
        }
        SECTION("Guess from header")
        {
            std::vector<std::string> input(6, "A");
            auto idx = GENERATE(take(10, random(0, 5)));
            auto header = GENERATE("RS.ID", "RS_ID", "RSID", "SNP.ID", "SNP_ID",
                                   "Variant_ID", "Variant.ID", "SNP");

            input[static_cast<size_t>(idx)] = header;
            std::string in;
            for (auto&& i : input) in.append(i + " ");
            REQUIRE(geno.test_get_rs_column(in) == static_cast<size_t>(idx));
        }
        SECTION("bim format")
        {
            std::vector<std::string> input(6, "A");
            std::string in;
            for (auto&& i : input) in.append(i + " ");
            REQUIRE(geno.test_get_rs_column(in) == 1);
        }
    }
    SECTION("load_snp_list")
    {
        std::unique_ptr<std::istringstream> input =
            std::make_unique<std::istringstream>();
        auto size = GENERATE(take(5, random(1, 20)));
        std::vector<std::string> expected;
        size_t idx = (size == 1 || size != 6) ? 0 : 1;
        std::string mock_file;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(1, 20);
        for (size_t nrow = 0; nrow < 3; ++nrow)
        {
            std::vector<std::string> row;
            for (int i = 0; i < size; ++i)
            { row.push_back(std::string(1, 'A' + i + dis(gen))); }
            std::string input_string = "";
            for (auto i : row) { input_string.append(i + "\t"); }
            mock_file.append(input_string + "\n");
            expected.push_back(row[idx]);
        }
        std::sort(expected.begin(), expected.end());
        expected.erase(std::unique(expected.begin(), expected.end()),
                       expected.end());
        input->str(mock_file);
        auto res = geno.test_load_snp_list(std::move(input));
        REQUIRE_THAT(res, Catch::UnorderedEquals<std::string>(expected));
    }
}
