#include "IITree.h"
#include "catch.hpp"
#include "mock_region.hpp"
#include "region.hpp"


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
            if (!gene_set[chr].has_overlap(i, out))
            { std::cerr << "Check: " << i << std::endl; }
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
        // included. otherwise, we want any SNP with coord of 1023-1027 to be
        // included
        auto equal_wind = GENERATE(true, false);
        // don't want the window too big as it will take forever to run the unit
        // test
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
            // warning should only be issued when file doesn't even containt the
            // strand info
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
    // chr >= max chr
}
