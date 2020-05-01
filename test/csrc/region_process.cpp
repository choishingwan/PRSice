#include "IITree.h"
#include "catch.hpp"
#include "gzstream.h"
#include "mock_region.hpp"
#include "region.hpp"
#include "snp.hpp"

TEST_CASE("Generate regions")
{
    // this is the super master test for region
    SECTION("No region")
    {
        Region region;
        region.generate_regions(std::unordered_map<std::string, size_t> {},
                                std::vector<SNP> {}, 22);
        REQUIRE_THAT(region.get_names(),
                     Catch::Equals<std::string>({"Base", "Background"}));
    }
}

TEST_CASE("Load background")
{
    mock_region region;
    Reporter reporter("log", 60, true);
    region.set_reporter(&reporter);
}
TEST_CASE("Load GTF")
{
    mock_region region;
    Reporter reporter("log", 60, true);
    region.set_reporter(&reporter);
    region.set_feature({"exon", "CDS", "gene"});
    std::unordered_map<std::string, std::vector<size_t>> msigdb_list;
    msigdb_list["DDX11L1"] = {GENERATE(take(1, random(1ul, 1025ul)))};
    msigdb_list["ENSG00000187634"] = {GENERATE(take(1, random(1ul, 1025ul)))};
    msigdb_list["ENSG00000175084"] = {GENERATE(take(1, random(1ul, 1025ul)))};
    msigdb_list["MT-TP"] = {GENERATE(take(1, random(1ul, 1025ul)))};
    size_t max_chr = 22;
    SECTION("malformed gtf")
    {
        auto invalid = GENERATE(
            "Space instead of tab though the size is "
            "correct",
            "Incorrect\tcolumn\tnumber",
            "1\tpseudogene\tgene\t11869\t"
            "14409\t.\t+\t.\tgene_name "
            "\"DDX11L1\"; gene_source \"havana\"; "
            "gene_biotype "
            "\"pseudogene\";" /*No gene id, therefore invalid*/,
            "1\tprocessed_transcript\texon\t-12613\t12721\t.\t+\t"
            ".\tgene_id \"ENSG00000223972\"; transcript_id "
            "\"ENST00000456328\"; "
            "exon_number \"2\"; gene_name \"DDX11L1\"; gene_source "
            "\"ensembl_havana\"; gene_biotype \"pseudogene\"; transcript_name "
            "\"DDX11L1-002\"; transcript_source \"havana\"; exon_id "
            "\"ENSE00003582793\";" /*Invalid start coordinate*/);
        auto input = std::make_unique<std::istringstream>(invalid);
        input->seekg(0, input->end);
        auto file_length = input->tellg();
        input->clear();
        input->seekg(0, input->beg);
        REQUIRE_THROWS(region.test_transverse_gtf(
            msigdb_list, file_length, max_chr, false, std::move(input)));
    }
    SECTION("transverse check")
    {
        auto has_background = GENERATE(true, false);
        auto BACKGROUND_IDX = 1ul;
        region.set_gwas_bk(has_background);
        // won't do gz here, can do it in the full file test where we can output
        // the gz file using gzstream
        auto input_str =
            "#!genome-build GRCh37.p13\n"
            "#!genome-version GRCh37\n"
            "#!genome-date 2009-02\n"
            "#!genome-build-accession NCBI:GCA_000001405.14\n"
            "#!genebuild-last-updated 2013-09\n"
            "1\tprocessed_transcript\texon\t12613\t12621\t.\t+\t.\tgene_id "
            "\"ENSG00000223972\"; transcript_id "
            "\"ENST00000456328\";exon_number \"2\"; gene_name \"DDX11L1\"; "
            "gene_source \"ensembl_havana\"; gene_biotype \"pseudogene\"; "
            "transcript_name \"DDX11L1-002\"; transcript_source "
            "\"havana\"; "
            "exon_id \"ENSE00003582793\";\n" /*chr1 12613 12621 DDX11L1*/
            "1\tprotein_coding\tUTR\t861302\t861321\t.\t+\t.\tgene_id "
            "\"ENSG00000187634\"; transcript_id \"ENST00000420190\"; "
            "gene_name "
            "\"SAMD11\"; gene_source \"ensembl_havana\"; gene_biotype "
            "\"protein_coding\"; transcript_name \"SAMD11-011\"; "
            "transcript_source \"havana\"; tag \"cds_end_NF\"; tag "
            "\"mRNA_end_NF\";\n" /* Feature excluded*/
            "2\tprotein_coding\tCDS\t220288499\t220288542\t.\t+\t1\tgene_"
            "id "
            "\"ENSG00000175084\"; transcript_id \"ENST00000373960\"; "
            "exon_number \"7\"; gene_name \"DES\"; gene_source "
            "\"ensembl_havana\"; gene_biotype \"protein_coding\"; "
            "transcript_name \"DES-001\"; transcript_source "
            "\"ensembl_havana\"; tag \"CCDS\"; ccds_id \"CCDS33383\"; "
            "protein_id \"ENSP00000363071\";\n" /*chr2 220288499 220288542
                                                   ENSG00000175084*/
            "2\tprotein_coding\tCDS\t179497253\t179497299\t.\t-\t2\tgene_"
            "id "
            "\"ENSG00000155657\"; transcript_id \"ENST00000359218\"; "
            "exon_number \"64\"; gene_name \"TTN\"; gene_source "
            "\"ensembl_havana\"; gene_biotype \"protein_coding\"; "
            "transcript_name \"TTN-202\"; transcript_source \"ensembl\"; "
            "tag "
            "\"CCDS\"; ccds_id \"CCDS54423\"; protein_id "
            "\"ENSP00000352154\";\n" /*Not in msigdb*/
            "MT\tMt_tRNA\tgene\t15956\t16023\t.\t-\t.\tgene_id "
            "\"ENSG00000210196\"; gene_name \"MT-TP\"; gene_source "
            "\"insdc\"; "
            "gene_biotype \"Mt_tRNA\";" /* chr filtered */;
        SECTION("direct call function")
        {
            auto input = std::make_unique<std::istringstream>(input_str);
            input->seekg(0, input->end);
            auto file_length = input->tellg();
            input->clear();
            input->seekg(0, input->beg);
            auto [num_line, exclude, chr_exclude] = region.test_transverse_gtf(
                msigdb_list, file_length, max_chr, false, std::move(input));
            REQUIRE(num_line == 5);
            REQUIRE(exclude == 1);
            REQUIRE(chr_exclude == 1);
            // check first gene is correctly set

            region.index();
            auto gene_set = region.get_gene_sets();
            REQUIRE(gene_set.size() == 3);
            auto idx = msigdb_list["DDX11L1"];
            if (!has_background) idx.push_back(BACKGROUND_IDX);
            std::vector<size_t> out;
            for (size_t i = 12613; i < 12621; ++i)
            {
                REQUIRE(gene_set[1].has_overlap(i, out));
                REQUIRE_THAT(out, Catch::UnorderedEquals<size_t>(idx));
            }
            REQUIRE_FALSE(gene_set[1].has_overlap(12612, out));
            REQUIRE_FALSE(gene_set[1].has_overlap(12622, out));
            idx = msigdb_list["ENSG00000175084"];
            if (!has_background) idx.push_back(BACKGROUND_IDX);
            for (size_t i = 220288499; i < 220288542; ++i)
            {
                REQUIRE(gene_set[2].has_overlap(i, out));
                REQUIRE_THAT(out, Catch::UnorderedEquals<size_t>(idx));
            }
            REQUIRE_FALSE(gene_set[2].has_overlap(220288498, out));
            REQUIRE_FALSE(gene_set[2].has_overlap(220288543, out));
            if (!has_background)
            {

                for (size_t i = 179497253; i < 179497299; ++i)
                {
                    REQUIRE(gene_set[2].has_overlap(i, out));
                    REQUIRE_THAT(out, Catch::UnorderedEquals<size_t>({1}));
                }
                REQUIRE_FALSE(gene_set[2].has_overlap(179497252, out));
                REQUIRE_FALSE(gene_set[2].has_overlap(179497300, out));
            }
        }
        SECTION("from load_gtf")
        {
            auto gz = GENERATE(true, false);
            if (!gz)
            {
                std::ofstream gtf("gtf_check");
                gtf << input_str << std::endl;
                gtf.close();
                region.set_gtf("gtf_check");
            }
            else
            {
                GZSTREAM_NAMESPACE::ogzstream gtf("gz_gtf_check");
                gtf << input_str << std::endl;
                gtf.close();
                region.set_gtf("gz_gtf_check");
            }
            region.test_load_gtf(msigdb_list, 22);
            region.index();
            auto gene_set = region.get_gene_sets();
            REQUIRE(gene_set.size() == 3);
            auto idx = msigdb_list["DDX11L1"];
            if (!has_background) idx.push_back(BACKGROUND_IDX);
            std::vector<size_t> out;
            for (size_t i = 12613; i < 12621; ++i)
            {
                REQUIRE(gene_set[1].has_overlap(i, out));
                REQUIRE_THAT(out, Catch::UnorderedEquals<size_t>(idx));
            }
            REQUIRE_FALSE(gene_set[1].has_overlap(12612, out));
            REQUIRE_FALSE(gene_set[1].has_overlap(12622, out));
            idx = msigdb_list["ENSG00000175084"];
            if (!has_background) idx.push_back(BACKGROUND_IDX);
            for (size_t i = 220288499; i < 220288542; ++i)
            {
                REQUIRE(gene_set[2].has_overlap(i, out));
                REQUIRE_THAT(out, Catch::UnorderedEquals<size_t>(idx));
            }
            REQUIRE_FALSE(gene_set[2].has_overlap(220288498, out));
            REQUIRE_FALSE(gene_set[2].has_overlap(220288543, out));
            if (!has_background)
            {

                for (size_t i = 179497253; i < 179497299; ++i)
                {
                    REQUIRE(gene_set[2].has_overlap(i, out));
                    REQUIRE_THAT(out, Catch::UnorderedEquals<size_t>({1}));
                }
                REQUIRE_FALSE(gene_set[2].has_overlap(179497252, out));
                REQUIRE_FALSE(gene_set[2].has_overlap(179497300, out));
            }
        }
    }
}
TEST_CASE("Full load bed file")
{
    mock_region region;
    Reporter reporter("log", 60, true);
    region.set_reporter(&reporter);
    region.set_wind(1, 2);
    std::ofstream bed_file("full_bed.test");
    bed_file << "chr3 123 148" << std::endl;
    bed_file << "chr3 144 152" << std::endl;
    bed_file.close();
    auto set_idx = GENERATE(take(1, random(1ul, 1026ul)));
    auto ori_idx = set_idx;
    using record = std::tuple<std::string, std::string>;
    auto input = GENERATE(table<std::string, std::string>(
        {record {"full_bed.test", "full_bed.test"},
         record {"full_bed.test:SetA", "SetA"}}));
    REQUIRE(region.test_load_bed_regions(std::get<0>(input), set_idx, 22));
    region.index();
    // we don't advance the set idx here
    REQUIRE(ori_idx == set_idx);
    auto name = region.get_names();
    REQUIRE_THAT(name, Catch::Equals<std::string>({std::get<1>(input)}));
    auto set = region.get_gene_sets();
    std::vector<size_t> out;
    for (size_t i = 123; i < 154; ++i)
    {
        REQUIRE(set[3].has_overlap(i, out));
        REQUIRE_THAT(out, Catch::Equals<size_t>({ori_idx}));
    }
    REQUIRE_FALSE(set[3].has_overlap(122, out));
    REQUIRE_FALSE(set[3].has_overlap(155, out));
}
TEST_CASE("load msigdb")
{
    mock_region region;
    Reporter reporter("log", 60, true);
    region.set_reporter(&reporter);
    std::unordered_map<std::string, std::vector<size_t>> msigdb_list;
    auto set_idx = GENERATE(take(1, random(1ul, 1025ul)));
    size_t ori_idx = set_idx;
    SECTION("Invalid msig")
    {
        auto input = std::make_unique<std::istringstream>("Set1");
        REQUIRE_THROWS(
            region.test_load_msigdb(msigdb_list, std::move(input), set_idx));
        REQUIRE(set_idx == ori_idx);
    }
    SECTION("valid msigdb")
    {
        auto input = std::make_unique<std::istringstream>("Set1 Gene1 Gene2\n"
                                                          "Set2 Hello world\n"
                                                          "Set1 Duplicate set\n"
                                                          "Set3 Ok Gene1 on\n");
        REQUIRE_NOTHROW(
            region.test_load_msigdb(msigdb_list, std::move(input), set_idx));
        // +3 because there is one duplicate

        REQUIRE(set_idx == ori_idx + 3);
        using record = std::tuple<std::string, std::vector<size_t>>;
        std::vector<record> expected = {{"Gene1", {ori_idx, ori_idx + 2}},
                                        {"Gene2", {ori_idx}},
                                        {"Hello", {ori_idx + 1}},
                                        {"world", {ori_idx + 1}},
                                        {"Ok", {ori_idx + 2}},
                                        {"on", {ori_idx + 2}}};

        for (auto&& exp : expected)
        {
            auto res = msigdb_list.find(std::get<0>(exp));
            REQUIRE(res != msigdb_list.end());
            REQUIRE_THAT(res->second,
                         Catch::Equals<size_t>({std::get<1>(exp)}));
        }
    }
}


TEST_CASE("Load snp sets")
{
    std::vector<SNP> snp_list;
    std::unordered_map<std::string, size_t> snp_list_idx;
    // generate 1000 fake SNPs
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> chr(1, 22);
    std::uniform_int_distribution<> bp(1, std::numeric_limits<int>::max());
    for (size_t i = 0; i < 1000; ++i)
    {
        snp_list_idx[std::to_string(i)] = i;
        snp_list.emplace_back(
            SNP(std::to_string(i), chr(gen), bp(gen), "A", "C", 0, 1));
    }
    mock_region region;
    Reporter reporter("log", 60, true);
    region.set_reporter(&reporter);
    auto set_idx = GENERATE(take(1, random(1ul, 10026ul)));
    SECTION("SNP set")
    {
        SECTION("Transverse file")
        {
            size_t ori_idx = set_idx;
            auto input =
                std::make_unique<std::istringstream>("Test 1 2 3\n"
                                                     "Test 2 3 4\n"
                                                     "Base 7 8 1\n"
                                                     "Set2 5 6 7 8 9\n");
            const bool is_snp_set = true;
            region.test_transverse_snp_file(snp_list_idx, snp_list, is_snp_set,
                                            std::move(input), set_idx);

            // + 3 because the second Test are ignored (usually we also
            // ignore base but as we haven;t gone through the proper
            // generate_region function, that wasn't set up yet
            REQUIRE(set_idx == ori_idx + 3);
            auto names = region.get_names();
            REQUIRE_THAT(names,
                         Catch::Equals<std::string>({"Test", "Base", "Set2"}));
            size_t max_chr = 0;
            for (auto&& i : {1, 2, 3, 5, 6, 7, 8, 9})
            {
                if (snp_list[i].chr() > max_chr) max_chr = snp_list[i].chr();
            }
            auto gene_sets = region.get_gene_sets();
            for (auto&& tree : gene_sets) { tree.index(); }
            REQUIRE(gene_sets.size() == max_chr + 1);
            std::vector<size_t> out;
            for (auto&& i : {1, 2, 3, 5, 6, 7, 8, 9})
            {
                auto&& cur_snp = snp_list[i];
                REQUIRE(
                    gene_sets[cur_snp.chr()].has_overlap(cur_snp.loc(), out));
                // don't test the not cases as we don't know if by chance
                // there are other SNPs that were randomly selected are
                // located in those space
                switch (i)
                {

                case 2:
                case 3:
                    CHECK_THAT(out, Catch::UnorderedEquals<size_t>({ori_idx}));
                    break;
                case 1:
                    CHECK_THAT(out, Catch::UnorderedEquals<size_t>(
                                        {ori_idx, ori_idx + 1}));
                    break;
                case 7:
                case 8:
                    CHECK_THAT(out, Catch::UnorderedEquals<size_t>(
                                        {ori_idx + 1, ori_idx + 2}));
                    break;
                case 5:
                case 6:
                case 9:
                    CHECK_THAT(out,
                               Catch::UnorderedEquals<size_t>({ori_idx + 2}));
                    break;
                }
            }
        }
        SECTION("Read through file")
        {
            std::ofstream snp_set("snp_set.test");
            snp_set << "Test 1 2 3\n"
                       "Test 2 3 4\n"
                       "Base 7 8 1\n"
                       "Set2 5 6 7 8 9\n";
            snp_set.close();
            SECTION("valid SNP Set file")
            {
                size_t ori_idx = set_idx;
                REQUIRE_NOTHROW(region.test_load_snp_sets(
                    snp_list_idx, snp_list, "snp_set.test", set_idx));
                REQUIRE(set_idx == ori_idx + 3);
                auto names = region.get_names();
                REQUIRE_THAT(names, Catch::Equals<std::string>(
                                        {"Test", "Base", "Set2"}));
                size_t max_chr = 0;
                for (auto&& i : {1, 2, 3, 5, 6, 7, 8, 9})
                {
                    if (snp_list[i].chr() > max_chr)
                        max_chr = snp_list[i].chr();
                }
                auto gene_sets = region.get_gene_sets();
                for (auto&& tree : gene_sets) { tree.index(); }
                REQUIRE(gene_sets.size() == max_chr + 1);
                std::vector<size_t> out;
                for (auto&& i : {1, 2, 3, 5, 6, 7, 8, 9})
                {
                    auto&& cur_snp = snp_list[i];
                    REQUIRE(gene_sets[cur_snp.chr()].has_overlap(cur_snp.loc(),
                                                                 out));
                    // don't test the not cases as we don't know if by
                    // chance there are other SNPs that were randomly
                    // selected are located in those space
                    switch (i)
                    {

                    case 2:
                    case 3:
                        CHECK_THAT(out,
                                   Catch::UnorderedEquals<size_t>({ori_idx}));
                        break;
                    case 1:
                        CHECK_THAT(out, Catch::UnorderedEquals<size_t>(
                                            {ori_idx, ori_idx + 1}));
                        break;
                    case 7:
                    case 8:
                        CHECK_THAT(out, Catch::UnorderedEquals<size_t>(
                                            {ori_idx + 1, ori_idx + 2}));
                        break;
                    case 5:
                    case 6:
                    case 9:
                        CHECK_THAT(
                            out, Catch::UnorderedEquals<size_t>({ori_idx + 2}));
                        break;
                    }
                }
            }
            SECTION("invalid SNP set input format")
            {
                // we don't allow user defined set name when snp set is
                // provided
                REQUIRE_THROWS(region.test_load_snp_sets(
                    snp_list_idx, snp_list, "snp_set.test:Test", set_idx));
            }
        }
    }
    SECTION("SNP list")
    {
        SECTION("load from file")
        {
            std::ofstream snp_set("snp_list.test");
            snp_set << "1\n"
                       "3\n"
                       "4\n"
                       "6\n";
            snp_set.close();
            size_t ori_idx = set_idx;
            using record = std::tuple<std::string, std::string>;
            auto expected = GENERATE(table<std::string, std::string>(
                {record {"snp_list.test", "snp_list.test"},
                 record {"snp_list.test:Test", "Test"}}));
            REQUIRE_NOTHROW(region.test_load_snp_sets(
                snp_list_idx, snp_list, std::get<0>(expected), set_idx));
            REQUIRE(set_idx == ori_idx + 1);
            auto names = region.get_names();
            REQUIRE_THAT(names,
                         Catch::Equals<std::string>({std::get<1>(expected)}));
            auto gene_sets = region.get_gene_sets();
            size_t max_chr = 0;
            std::unordered_set<size_t> expected_loc;
            for (auto&& i : {1, 3, 4, 6})
            {
                if (snp_list[i].chr() > max_chr) max_chr = snp_list[i].chr();
                expected_loc.insert(snp_list[i].loc());
            }
            REQUIRE(gene_sets.size() == max_chr + 1);
            std::vector<size_t> out;
            for (auto&& i : {1, 3, 4, 6})
            {
                auto&& cur_snp = snp_list[i];
                REQUIRE(
                    gene_sets[cur_snp.chr()].has_overlap(cur_snp.loc(), out));
                REQUIRE_THAT(out, Catch::Equals<size_t>({ori_idx}));
                if (expected_loc.find(cur_snp.loc() - 1) == expected_loc.end())
                {
                    REQUIRE_FALSE(gene_sets[cur_snp.chr()].has_overlap(
                        cur_snp.loc() - 1, out));
                }
                if (expected_loc.find(cur_snp.loc() + 1) == expected_loc.end())
                {
                    REQUIRE_FALSE(gene_sets[cur_snp.chr()].has_overlap(
                        cur_snp.loc() + 1, out));
                }
            }
        }
        SECTION("Transverse file")
        {
            size_t ori_idx = set_idx;
            auto input = std::make_unique<std::istringstream>("1\n"
                                                              "3\n"
                                                              "4\n"
                                                              "6\n");
            const bool is_snp_set = false;
            region.test_transverse_snp_file(snp_list_idx, snp_list, is_snp_set,
                                            std::move(input), set_idx);
            // won't change as we need to do that manually
            REQUIRE(set_idx == ori_idx);
            auto names = region.get_names();
            REQUIRE(names.empty());
            auto gene_sets = region.get_gene_sets();
            size_t max_chr = 0;
            for (auto&& i : {1, 3, 4, 6})
            {
                if (snp_list[i].chr() > max_chr) max_chr = snp_list[i].chr();
            }
            REQUIRE(gene_sets.size() == max_chr + 1);
            std::vector<size_t> out;
            for (auto&& i : {1, 3, 4, 6})
            {
                auto&& cur_snp = snp_list[i];
                REQUIRE(
                    gene_sets[cur_snp.chr()].has_overlap(cur_snp.loc(), out));
                REQUIRE_THAT(out, Catch::Equals<size_t>({ori_idx}));
            }
        }
    }
}

TEST_CASE("Read bed file")
{
    std::vector<IITree<size_t, size_t>> gene_set;
    bool warning = false;
    size_t wind_5 = 0, wind_3 = 0, max_chr = 22;
    auto set_idx = GENERATE(take(2, random(0ul, 1025ul)));
    auto zero_based = GENERATE(true, false);
    mock_region region;
    SECTION("invalid bed file")
    {
        auto input = std::make_unique<std::istringstream>("chr1 10");
        REQUIRE_THROWS(region.test_read_bed(std::move(input), gene_set, warning,
                                            wind_5, wind_3, max_chr, set_idx,
                                            zero_based));
    }
    SECTION("invalid start end")
    {
        auto input = std::make_unique<std::istringstream>("chr1 10 -10");
        REQUIRE_THROWS(region.test_read_bed(std::move(input), gene_set, warning,
                                            wind_5, wind_3, max_chr, set_idx,
                                            zero_based));
    }
    SECTION("chromosome too big")
    {
        size_t chr = 100;
        size_t start = 1;
        size_t end = 1000;
        auto input = std::make_unique<std::istringstream>(
            "chr" + std::to_string(chr) + " " + std::to_string(start) + " "
            + std::to_string(end));
        region.test_read_bed(std::move(input), gene_set, warning, wind_5,
                             wind_3, max_chr, set_idx, zero_based);
        // empty because we never read in any chromosome
        REQUIRE(gene_set.empty());
    }
    SECTION("overlapping regions")
    {
        auto input = std::make_unique<std::istringstream>("chr1 123 346\n"
                                                          "chr1 246 348\n");
        region.test_read_bed(std::move(input), gene_set, warning, wind_5,
                             wind_3, max_chr, set_idx, zero_based);
        for (auto&& tree : gene_set) { tree.index(); }
        size_t exp_start = 123 + zero_based - wind_5;
        size_t exp_end = 348 + 1 + wind_3;
        size_t chr = 1;
        std::vector<size_t> out;
        for (size_t i = exp_start; i < exp_end; ++i)
        {
            out.clear();
            REQUIRE(gene_set[chr].has_overlap(i, out));
            REQUIRE_THAT(out, Catch::Equals<size_t>({set_idx}));
        }
        out.clear();
        REQUIRE_FALSE(gene_set[chr].has_overlap(exp_start - 1, out));
        out.clear();
        REQUIRE_FALSE(gene_set[chr].has_overlap(exp_end, out));
    }
    SECTION("valid input")
    {

        size_t chr = GENERATE(take(1, random(1ul, 22ul)));
        size_t start = GENERATE(
            take(1, random(1ul, std::numeric_limits<int>::max() / 2ul)));
        size_t end = start + GENERATE(take(1, random(1ul, 50ul)));

        // if zero based we want any SNP with coord of 1024 - 1027 to be
        // included. otherwise, we want any SNP with coord of 1023-1027 to
        // be included
        auto equal_wind = GENERATE(true, false);
        // don't want the window too big as it will take forever to run the
        // unit test
        auto wind_5 = GENERATE(take(1, random(1ul, 20ul)));
        size_t wind_3;
        if (equal_wind) { wind_3 = wind_5; }
        else
        {
            wind_3 = wind_5 + GENERATE(take(1, random(1ul, 10ul)));
        }
        SECTION("no strand info")
        {
            auto input = std::make_unique<std::istringstream>(
                "chr" + std::to_string(chr) + " " + std::to_string(start) + " "
                + std::to_string(end));
            region.test_read_bed(std::move(input), gene_set, warning, wind_5,
                                 wind_3, max_chr, set_idx, zero_based);
            // now index the vector
            for (auto&& tree : gene_set) { tree.index(); }
            if (!equal_wind) { REQUIRE(warning); }
            REQUIRE(gene_set.size() == chr + 1);
            size_t exp_start = start + zero_based - wind_5;
            size_t exp_end = end + 1 + wind_3;
            std::vector<size_t> out;
            for (size_t i = exp_start; i < exp_end; ++i)
            {
                out.clear();
                REQUIRE(gene_set[chr].has_overlap(i, out));
                REQUIRE_THAT(out, Catch::Equals<size_t>({set_idx}));
            }
            out.clear();
            REQUIRE_FALSE(gene_set[chr].has_overlap(exp_start - 1, out));
            out.clear();
            REQUIRE_FALSE(gene_set[chr].has_overlap(exp_end, out));
        }
        SECTION("with strand info")
        {
            std::string strand = GENERATE("-", "+", ".");

            auto input = std::make_unique<std::istringstream>(
                "chr" + std::to_string(chr) + " " + std::to_string(start) + " "
                + std::to_string(end) + " test 0 " + strand);
            region.test_read_bed(std::move(input), gene_set, warning, wind_5,
                                 wind_3, max_chr, set_idx, zero_based);
            // warning should only be issued when file doesn't even containt
            // the strand info
            REQUIRE(gene_set.size() == chr + 1);
            size_t exp_start = start + zero_based - wind_5;
            if (start + zero_based <= wind_5) { exp_start = 1; }

            size_t exp_end = end + 1 + wind_3;
            if (strand == "-")
            {
                exp_start = start + zero_based - wind_3;
                if (start + zero_based <= wind_3) { exp_start = 1; }
                exp_end = end + 1 + wind_5;
            }
            std::vector<size_t> out;
            for (size_t i = exp_start; i < exp_end; ++i)
            {
                REQUIRE(gene_set[chr].has_overlap(i, out));
                REQUIRE_THAT(out, Catch::Equals<size_t>({set_idx}));
            }
            out.clear();
            REQUIRE_FALSE(gene_set[chr].has_overlap(exp_start - 1, out));
            out.clear();
            REQUIRE_FALSE(gene_set[chr].has_overlap(exp_end, out));
        }
    }
}
