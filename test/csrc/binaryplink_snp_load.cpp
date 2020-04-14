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

    std::ofstream bim("load_snp1.bim");
    // first chr only good SNP * 3
    bim << "1	SNP_1	0	742429	A	C" << std::endl;
    bim << "1	SNP_2	0	933331	C	T" << std::endl;
    bim << "1	SNP_4	0	1008567	G	A" << std::endl;
    bim.close();
    bplink.gen_bed_head("load_snp1.bed", num_sample, 3, true, false);
    bplink.manual_load_snp(SNP("SNP_1", 1, 742429, "A", "C", 1, 1));
    bplink.manual_load_snp(SNP("SNP_2", 1, 933331, "C", "T", 1, 1));
    bplink.manual_load_snp(SNP("SNP_3", 1, 933331, "C", "T", 1, 1));
    bplink.manual_load_snp(SNP("SNP_4", 1, 1008567, "G", "A", 1, 1));
    bplink.manual_load_snp(SNP("SNP_5", 1, 742429, "A", "C", 1, 1));
    // load SNP

    SECTION("invalid SNPs")
    {
        bim.open("load_snp2.bim");
        auto input = GENERATE("chr1	SNP_5	0	-742429	A	C",
                              "chr1	SNP_5	0	742429	A",
                              "chr1	SNP_5	0	742429	A C G");
        bim << input << std::endl;
        bim.close();
        bplink.gen_bed_head("load_snp2.bed", num_sample, 1, true, false);
        REQUIRE_THROWS(bplink.load_snps(
            "load_snp", std::vector<IITree<size_t, size_t>> {}, false));
    }
    SECTION("filtered out except mismatch")
    {
        auto input = GENERATE("chrX	SNP_5	0	742429	A	C",
                              "chr1	SNP_5	0	742429	A	G");
        bim.open("load_snp2.bim");
        bim << input << std::endl;
        bim.close();
        bplink.gen_bed_head("load_snp2.bed", num_sample, 1, true, false);
        bplink.load_snps("load_snp", std::vector<IITree<size_t, size_t>> {},
                         false);
        auto res = bplink.existed_snps();
        REQUIRE(res.size() == 3);
        REQUIRE(res[0].rs() == "SNP_1");
        REQUIRE(res[1].rs() == "SNP_2");
        REQUIRE(res[2].rs() == "SNP_4");
    }
    SECTION("mismatch")
    {
        auto input = "chr1	SNP_5	0	742430	A	C";
        bim.open("load_snp2.bim");
        bim << input << std::endl;
        bim.close();
        bplink.gen_bed_head("load_snp2.bed", num_sample, 1, true, false);
        bplink.load_snps("load_snp", std::vector<IITree<size_t, size_t>> {},
                         false);
        auto res = bplink.existed_snps();
        REQUIRE(res.size() == 3);
        REQUIRE(res[0].rs() == "SNP_1");
        REQUIRE(res[1].rs() == "SNP_2");
        REQUIRE(res[2].rs() == "SNP_4");
        std::ifstream mismatch("load_snp.mismatch");
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
    SECTION("duplicated SNP")
    {
        bim.open("load_snp2.bim");
        bim << "1	SNP_5	0	742429	A	C" << std::endl;
        bim << "1	SNP_5	0	742429	A	C" << std::endl;
        bim.close();
        bplink.gen_bed_head("load_snp2.bed", num_sample, 2, true, false);
        std::cerr << "Start duplicate check" << std::endl;
        REQUIRE_THROWS(bplink.load_snps(
            "load_snp", std::vector<IITree<size_t, size_t>> {}, false));
        // we only output the valid SNPs
        std::ifstream dup("load_snp.valid");
        REQUIRE(dup.is_open());
        std::string line;
        size_t num_dup = 0;
        std::vector<std::string> token;
        while (std::getline(dup, line))
        {
            token = misc::split(line, "\t");
            ++num_dup;
        }
        // because of header
        REQUIRE(num_dup == 3);
        REQUIRE(token.size() == 5);
        REQUIRE(token[0] == "SNP_4");
    }

    // one invalid loc
    // one invalid chr
    // one invalid anything
    // one valid
}
