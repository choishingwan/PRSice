#include "binaryplink.hpp"
#include "catch.hpp"
#include "mock_binaryplink.hpp"


TEST_CASE("check bed format")
{
    mock_binaryplink plink;
    SECTION("modern format")
    {
        bool mordan_snp = GENERATE(true, false);
        bool sample_major = GENERATE(true, false);
        std::string name = "header_check.bed";
        size_t num_snp = GENERATE(take(5, random(1ul, 100ul)));
        size_t num_sample = GENERATE(take(5, random(1ul, 100ul)));
        plink.gen_bed_head(name, num_sample, num_snp, mordan_snp, sample_major);
        plink.set_sample(num_sample);
        uintptr_t offset;
        if (!sample_major)
        {
            REQUIRE_NOTHROW(plink.test_check_bed(name, num_snp, offset));
            if (mordan_snp) { REQUIRE(offset == 3); }
            else
            {
                REQUIRE(offset == 1);
            }
        }
        else
        {
            REQUIRE_THROWS(plink.test_check_bed(name, num_snp, offset));
        }
    }
    SECTION("UCSC bed")
    {
        std::ofstream bed("ucsc.bed");
        bed << "# chr" << std::endl;
        bed.close();
        uintptr_t offset;
        REQUIRE_THROWS(plink.test_check_bed("ucsc.bed", 10, offset));
    }
}

TEST_CASE("generate snp vector")
{
    Reporter reporter("log", 60, true);
    GenoFile geno;
    geno.num_autosome = 2;
    geno.file_name = "load_snp#";
    Phenotype pheno;
    mock_binaryplink bplink(geno, pheno, " ", &reporter);
    size_t num_sample = 10;
    bplink.set_sample(num_sample);

    size_t idx = 0;
    bplink.gen_bed_head("load_snp1.bed", num_sample, 3, true, false);
    bplink.manual_load_snp(SNP("SNP_1", 1, 742429, "A", "C", idx, 1));
    bplink.manual_load_snp(SNP("SNP_2", 1, 933331, "C", "T", idx, 1));
    bplink.manual_load_snp(SNP("SNP_3", 1, 933331, "C", "T", idx, 1));
    bplink.manual_load_snp(SNP("SNP_4", 1, 1008567, "G", "A", idx, 1));
    bplink.manual_load_snp(SNP("SNP_5", 1, 742429, "A", "C", idx, 1));
    // load SNP
    std::vector<IITree<size_t, size_t>> exclusion_region;
    std::string mismatch_name = "mismatch";
    uintptr_t unfiltered_sample_ct4 = 0;
    uintptr_t bed_offset = 4;
    std::unordered_set<std::string> duplicated_snps;
    std::unordered_set<std::string> processed_snps;
    std::vector<bool> retain_snp(5, false);
    bool chr_error = false;
    bool sex_error = false;
    SECTION("invalid SNPs")
    {

        auto str =
            GENERATE("chr1	SNP_5	0	-742429	A	C", /* Negative coord*/
                     "chr1	SNP_5	0	742429	A" /* malformed bed*/);
        std::unique_ptr<std::istream> input =
            std::make_unique<std::istringstream>(str);
        REQUIRE_THROWS(bplink.test_transverse_bed_for_snp(
            exclusion_region, mismatch_name, idx, unfiltered_sample_ct4,
            bed_offset, std::move(input), duplicated_snps, processed_snps,
            retain_snp, chr_error, sex_error, &bplink));
    }
    SECTION("filtered out except mismatch")
    {
        std::unique_ptr<std::istream> input =
            std::make_unique<std::istringstream>(
                "chrX	SNP_5	0	742429	A	C");
        // == 0 as no SNP loaded
        REQUIRE(bplink.test_transverse_bed_for_snp(
                    exclusion_region, mismatch_name, idx, unfiltered_sample_ct4,
                    bed_offset, std::move(input), duplicated_snps,
                    processed_snps, retain_snp, chr_error, sex_error, &bplink)
                == 0);
        auto res = bplink.existed_snps();
        // we shouldn't have reshaped our vector here
        REQUIRE(res.size() == 5);
        REQUIRE(sex_error);
    }
    SECTION("mismatch")
    {
        auto str = GENERATE("chr1	SNP_5	0	742430	A	C",
                            "chr1	SNP_5	0	742429	A	G");
        std::unique_ptr<std::istream> input =
            std::make_unique<std::istringstream>(str);
        REQUIRE(bplink.test_transverse_bed_for_snp(
                    exclusion_region, mismatch_name, idx, unfiltered_sample_ct4,
                    bed_offset, std::move(input), duplicated_snps,
                    processed_snps, retain_snp, chr_error, sex_error, &bplink)
                == 0);
        std::ifstream mismatch(mismatch_name.c_str());
        REQUIRE(mismatch.is_open());
        std::string line;
        size_t num_dup = 0;
        std::vector<std::string> token;
        while (std::getline(mismatch, line))
        {
            token = misc::split(line, "\t");
            ++num_dup;
        }
        // because of header
        REQUIRE(num_dup == 2);
        REQUIRE(token.size() == 10);
        REQUIRE(token[1] == "SNP_5");
    }
    SECTION("valid")
    {
        std::unique_ptr<std::istream> input =
            std::make_unique<std::istringstream>(
                "chr1	SNP_5	0	742429	A	C");
        REQUIRE(bplink.test_transverse_bed_for_snp(
                    exclusion_region, mismatch_name, idx + 1,
                    unfiltered_sample_ct4, bed_offset, std::move(input),
                    duplicated_snps, processed_snps, retain_snp, chr_error,
                    sex_error, &bplink)
                == 1);
        auto res = bplink.existed_snps();
        REQUIRE(res.size() == 5);
        REQUIRE(res[0].rs() == "SNP_1");
        REQUIRE(res[1].rs() == "SNP_2");
        REQUIRE(res[2].rs() == "SNP_3");
        REQUIRE(res[3].rs() == "SNP_4");
        REQUIRE(res[4].rs() == "SNP_5");
        REQUIRE(res[4].get_file_idx() == idx + 1);
    }
    SECTION("duplicated SNP")
    {
        std::unique_ptr<std::istream> input =
            std::make_unique<std::istringstream>("1	SNP_5	0	742429	A	C\n"
                                                 "1	SNP_5	0	742429	A	C");
        REQUIRE(bplink.test_transverse_bed_for_snp(
                    exclusion_region, mismatch_name, idx + 1,
                    unfiltered_sample_ct4, bed_offset, std::move(input),
                    duplicated_snps, processed_snps, retain_snp, chr_error,
                    sex_error, &bplink)
                == 1);
        REQUIRE(duplicated_snps.size() == 1);
        REQUIRE(duplicated_snps.find("SNP_5") != duplicated_snps.end());
    }
    SECTION("Full test")
    {
        // this will test the gen_snp
        std::ofstream bim("load_snp1.bim");
        // first chr only good SNP * 3
        bim << "1	SNP_1	0	742429	A	C" << std::endl;
        bim << "1	SNP_2	0	933331	C	T" << std::endl;
        bim << "1	SNP_4	0	1008567	G	A" << std::endl;
        bim.close();
        bim.open("load_snp2.bim");
        bim << "chr1	SNP_5	0	742429	A	C" << std::endl;
        bim.close();
        bplink.gen_bed_head("load_snp2.bed", num_sample, 1, true, false);
        bplink.load_snps("load_snp", std::vector<IITree<size_t, size_t>> {},
                         false);
        auto res = bplink.existed_snps();
        REQUIRE(res.size() == 4);
        REQUIRE(res[0].rs() == "SNP_1");
        REQUIRE(res[1].rs() == "SNP_2");
        REQUIRE(res[2].rs() == "SNP_4");
        REQUIRE(res[3].rs() == "SNP_5");
    }

    // one invalid loc
    // one invalid chr
    // one invalid anything
    // one valid
}
