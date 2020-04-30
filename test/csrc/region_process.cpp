#include "IITree.h"
#include "catch.hpp"
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

TEST_CASE("Load background") {}
TEST_CASE("Load GTF") {}
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

            // + 3 because the second Test are ignored (usually we also ignore
            // base but as we haven;t gone through the proper generate_region
            // function, that wasn't set up yet
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
                // don't test the not cases as we don't know if by chance there
                // are other SNPs that were randomly selected are located in
                // those space
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
                    // don't test the not cases as we don't know if by chance
                    // there are other SNPs that were randomly selected are
                    // located in those space
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
                // we don't allow user defined set name when snp set is provided
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
