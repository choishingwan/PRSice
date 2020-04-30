#include "catch.hpp"
#include "genotype.hpp"
#include "region.hpp"

TEST_CASE("Exclusion test")
{
    std::vector<IITree<size_t, size_t>> gene_set;
    SECTION("Bed file")
    {
        std::string input = "exclude.bed";
        std::ofstream bed(input);
        bed << "chr1 123 145 test 0 +";
        bed.close();
        Region::generate_exclusion(gene_set, input);
        REQUIRE(gene_set.size() == 2);
        std::vector<size_t> out;
        for (size_t i = 124; i <= 145; ++i)
        { REQUIRE(gene_set[1].has_overlap(i, out)); }
        REQUIRE_FALSE(gene_set[1].has_overlap(123, out));
        REQUIRE_FALSE(gene_set[1].has_overlap(146, out));
    }
    SECTION("combination")
    {
        std::ofstream bed("multi_exclude.bed");
        bed << "chr1 123 145 test 0 +";
        bed.close();
        Region::generate_exclusion(
            gene_set, "chr22:124387-124390,multi_exclude.bed,chr1:456");
        REQUIRE(gene_set.size() == 23);
        std::vector<size_t> out;
        // first range
        for (size_t i = 124387; i <= 124390; ++i)
        { REQUIRE(gene_set[22].has_overlap(i, out)); }
        REQUIRE_FALSE(gene_set[22].has_overlap(124386, out));
        REQUIRE_FALSE(gene_set[22].has_overlap(124391, out));
        // bed file (that's why we need to +1)
        for (size_t i = 124; i <= 145; ++i)
        { REQUIRE(gene_set[1].has_overlap(i, out)); }
        REQUIRE_FALSE(gene_set[1].has_overlap(123, out));
        REQUIRE_FALSE(gene_set[1].has_overlap(146, out));

        // final single base
        REQUIRE_FALSE(gene_set[1].has_overlap(455, out));
        REQUIRE(gene_set[1].has_overlap(456, out));
        REQUIRE_FALSE(gene_set[1].has_overlap(457, out));
    }
    SECTION("Direct input")
    {
        SECTION("invalid input")
        {
            auto invalid = GENERATE("chr1:123-1234-123", "chr1:1:2:3",
                                    "chr1:1234--12", "chr12-123:33");
            REQUIRE_THROWS(Region::generate_exclusion(gene_set, invalid));
        }
        SECTION("single base")
        {
            using record = std::tuple<std::string, size_t, size_t>;
            auto input = GENERATE(table<std::string, size_t, size_t>(
                {record {"chr1:123", 1, 123}, record {"chr22:446", 22, 446}}));
            Region::generate_exclusion(gene_set, std::get<0>(input));
            std::vector<size_t> out;
            REQUIRE(gene_set.size() == std::get<1>(input) + 1);
            REQUIRE(gene_set[std::get<1>(input)].has_overlap(std::get<2>(input),
                                                             out));
            REQUIRE_FALSE(gene_set[std::get<1>(input)].has_overlap(
                std::get<2>(input) - 1, out));
            REQUIRE_FALSE(gene_set[std::get<1>(input)].has_overlap(
                std::get<2>(input) + 1, out));
        }
        SECTION("range")
        {
            using record = std::tuple<std::string, size_t, size_t, size_t>;
            auto input = GENERATE(table<std::string, size_t, size_t, size_t>(
                {record {"chr1:123-143", 1, 123, 143},
                 record {"chr22:446-450", 22, 446, 450}}));
            Region::generate_exclusion(gene_set, std::get<0>(input));
            std::vector<size_t> out;
            REQUIRE(gene_set.size() == std::get<1>(input) + 1);
            for (size_t i = std::get<2>(input); i <= std::get<3>(input); ++i)
            { REQUIRE(gene_set[std::get<1>(input)].has_overlap(i, out)); }
            REQUIRE_FALSE(gene_set[std::get<1>(input)].has_overlap(
                std::get<2>(input) - 1, out));
            REQUIRE_FALSE(gene_set[std::get<1>(input)].has_overlap(
                std::get<3>(input) + 1, out));
        }
    }
}
