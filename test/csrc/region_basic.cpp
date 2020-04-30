#include "catch.hpp"
#include "genotype.hpp"
#include "mock_region.hpp"
#include "region.hpp"

TEST_CASE("Feature found")
{
    std::vector<std::string> feature;
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()};
    std::uniform_int_distribution<size_t> rand_name {1, 20};
    auto n_features = GENERATE(take(2, random(1ul, 10ul)));
    for (size_t i = 0; i < n_features; ++i)
    { feature.push_back(std::to_string(rand_name(mersenne_engine))); }
    mock_region region;
    SECTION("Expect found")
    {
        std::uniform_int_distribution<size_t> rand_idx {0, n_features - 1};
        for (size_t i = 0; i < 10; ++i)
        {
            auto search = feature[rand_idx(mersenne_engine)];
            REQUIRE(region.test_in_region(search, feature));
        }
    }
    SECTION("Expect not found")
    {
        std::uniform_int_distribution<size_t> out_range {
            21, std::numeric_limits<size_t>::max()};
        for (size_t i = 0; i < 10; ++i)
        {
            REQUIRE_FALSE(region.test_in_region(
                std::to_string(out_range(mersenne_engine)), feature));
        }
    }
}

TEST_CASE("Strand check")
{
    SECTION("Valid strand")
    {
        auto input = GENERATE("+", "-", ".");
        REQUIRE(mock_region::test_valid_strand(input));
    }
    SECTION("Invalid strand")
    {
        REQUIRE_FALSE(mock_region::test_valid_strand("@"));
    }
}

TEST_CASE("get set name")
{
    mock_region region;
    SECTION("single input")
    {
        auto [file_name, set_name, user] = region.test_get_set_name("Input");
        REQUIRE(file_name == set_name);
        REQUIRE(file_name == "Input");
        REQUIRE_FALSE(user);
    }
    SECTION("with user input")
    {
        auto [file_name, set_name, user] =
            region.test_get_set_name("Input:Name");
        REQUIRE(user);
        REQUIRE(file_name == "Input");
        REQUIRE(set_name == "Name");
    }
    SECTION("invalid input")
    {
        REQUIRE_THROWS(region.test_get_set_name("Input:Name:fun"));
    }
}

TEST_CASE("region extension")
{
    SECTION("Positive strand")
    {
        // . is considered as +
        std::string strand = GENERATE("+", ".");
        size_t wind_5 = 100, wind_3 = 20, end = 1000;
        SECTION("wind bigger than start")
        {
            size_t start = 10;
            REQUIRE_NOTHROW(mock_region::test_extend_region(
                strand, wind_5, wind_3, start, end));
            REQUIRE(start == 1);
            REQUIRE(end == 1000 + wind_3);
        }
        SECTION("wind smaller than start")
        {
            size_t start = 500;
            REQUIRE_NOTHROW(mock_region::test_extend_region(
                strand, wind_5, wind_3, start, end));
            REQUIRE(start == 500 - wind_5);
            REQUIRE(end == 1000 + wind_3);
        }
    }
    SECTION("Negative strand")
    {
        std::string strand = "-";
        size_t wind_5 = 10, wind_3 = 20, end = 1000;
        SECTION("wind bigger than start")
        {
            size_t start = 10;
            REQUIRE_NOTHROW(mock_region::test_extend_region(
                strand, wind_5, wind_3, start, end));
            REQUIRE(start == 1);
            REQUIRE(end == 1000 + wind_5);
        }
        SECTION("wind smaller than start")
        {
            size_t start = 100;
            REQUIRE_NOTHROW(mock_region::test_extend_region(
                strand, wind_5, wind_3, start, end));
            REQUIRE(start == 100 - wind_3);
            REQUIRE(end == 1000 + wind_5);
        }
    }
    SECTION("invalid inputs")
    {
        std::string strand = "@";
        size_t wind_5 = 10, wind_3 = 10, start = 100, end = 10000;
        REQUIRE_THROWS(mock_region::test_extend_region(strand, wind_5, wind_3,
                                                       start, end));
    }
}

TEST_CASE("Duplicate set check")
{
    mock_region region;
    Reporter reporter("log", 60, true);
    region.set_reporter(&reporter);
    // first region will always be empty
    std::string region_name = "Testing";
    SECTION("string")
    {
        REQUIRE_FALSE(region.test_duplicated_set(region_name));
        REQUIRE(region.test_duplicated_set(region_name));
    }
    SECTION("string view")
    {
        std::string_view r = region_name;
        REQUIRE_FALSE(region.test_duplicated_set(r));
        REQUIRE(region.test_duplicated_set(r));
    }
}

TEST_CASE("Start end parsing")
{
    SECTION("invalid string input")
    {
        auto invalid = GENERATE("-123", "ABC", "0.05");
        std::string_view invalid_in = invalid;
        auto valid = "12345";
        std::string_view valid_in = valid;
        auto zero = GENERATE(true, false);
        SECTION("invalid start")
        {
            REQUIRE_THROWS_WITH(
                mock_region::test_start_end(invalid_in, valid_in, zero),
                Catch::Contains("start"));
        }
        SECTION("invalid end")
        {
            REQUIRE_THROWS_WITH(
                mock_region::test_start_end(valid_in, invalid_in, zero),
                Catch::Contains("end"));
        }
    }
    SECTION("invalid size")
    {
        auto zero = GENERATE(true, false);
        REQUIRE_THROWS_WITH(
            mock_region::test_start_end(std::string("12345"),
                                        std::string("1234"), zero),
            Catch::Contains("Error: Start coordinate should be smaller"));
    }
    SECTION("Zero based")
    {
        std::string start = "123", end = "12345";
        auto [s, e] = mock_region::test_start_end(start, end, true);
        REQUIRE(s == 124);
        REQUIRE(e == 12345);
    }
    SECTION("Not zero based")
    {
        std::string start = "123", end = "12345";
        auto [s, e] = mock_region::test_start_end(start, end, false);
        REQUIRE(s == 123);
        REQUIRE(e == 12345);
    }
}


TEST_CASE("bed header check")
{
    mock_region region;
    Reporter reporter("log", 60, true);
    region.set_reporter(&reporter);
    SECTION("is headers")
    {
        auto header = GENERATE("track number is something", "browser firefox",
                               "#Chr size something");
        auto token = misc::tokenize(header);
        auto column_size = GENERATE(take(1, random(1ul, 1000ul)));
        auto ori_csize = column_size;
        REQUIRE(region.test_is_bed_header(token, column_size));
        // shouldn't update the column size
        REQUIRE(ori_csize == column_size);
    }
    SECTION("invalid bed file")
    {
        auto line = "A B";
        auto token = misc::tokenize(line);
        size_t column_size = 0;
        REQUIRE_THROWS(region.test_is_bed_header(token, column_size));
    }
    SECTION("Not header")
    {
        SECTION("column number check")
        {
            auto line = "A B C D";
            auto token = misc::tokenize(line);
            SECTION("valid column size")
            {
                size_t column_size = 4;
                REQUIRE_FALSE(region.test_is_bed_header(token, column_size));
                REQUIRE(column_size == 4);
            }
            SECTION("invalid column size")
            {
                auto column_size = GENERATE(3ul, 5ul);
                REQUIRE_THROWS(region.test_is_bed_header(token, column_size));
            }
            SECTION("without column")
            {
                auto column_size = 0ul;
                REQUIRE_FALSE(region.test_is_bed_header(token, column_size));
                REQUIRE(column_size == token.size());
            }
        }
        SECTION("Invalid strand information")
        {
            auto line = "1 102 123 Test 0 @";
            auto token = misc::tokenize(line);
            size_t column_size = 0;
            REQUIRE_THROWS(region.test_is_bed_header(token, column_size));
        }
    }
}

TEST_CASE("Parse gene id")
{
    mock_region region;
    SECTION("malformed")
    {
        std::string input = "gene_id:test";
        REQUIRE_THROWS(region.test_parse_gene_id(input, "gene_id"));
    }
    SECTION("Find gene id")
    {
        SECTION("found")
        {
            std::string input =
                GENERATE("gene_id test", " gene_id test", "gene_id \"test\"");
            auto [s, b] = region.test_parse_gene_id(input, "gene_id");
            REQUIRE(b);
            REQUIRE(s == "test");
        }
        SECTION("not found")
        {
            std::string input = "gene_name:test";
            auto [s, b] = region.test_parse_gene_id(input, "gene_id");
            REQUIRE_FALSE(b);
            REQUIRE(s.empty());
        }
    }
    SECTION("find gene name")
    {
        SECTION("found")
        {
            std::string input = GENERATE("gene_name test", " gene_name test",
                                         "gene_name \"test\"");
            auto [s, b] = region.test_parse_gene_id(input, "gene_name");
            REQUIRE(b);
            REQUIRE(s == "test");
        }
        SECTION("not found")
        {
            std::string input = "gene_id test";
            auto [s, b] = region.test_parse_gene_id(input, "gene_name");
            REQUIRE_FALSE(b);
            REQUIRE(s.empty());
        }
    }
}

TEST_CASE("find gene info")
{
    mock_region region;
    std::string gene_id, gene_name;
    bool found_id = false, found_name = false;
    SECTION("not found")
    {
        std::string input = "gene_source \"havana\"";
        REQUIRE_FALSE(region.test_find_gene_info(input, gene_id, gene_name,
                                                 found_id, found_name));
        REQUIRE_FALSE(found_id);
        REQUIRE_FALSE(found_name);
    }
    SECTION("Find gene id")
    {

        std::string input = "gene_id \"ENSG00000223972\"";
        // false because gene name not found yet
        REQUIRE_FALSE(region.test_find_gene_info(input, gene_id, gene_name,
                                                 found_id, found_name));
        REQUIRE(found_id);
        REQUIRE(gene_id == "ENSG00000223972");
        REQUIRE_FALSE(found_name);
        SECTION("Find gene name")
        {
            input = "gene_name \"DDX11L1\"";
            REQUIRE(region.test_find_gene_info(input, gene_id, gene_name,
                                               found_id, found_name));
            REQUIRE(found_id);
            REQUIRE(gene_id == "ENSG00000223972");
            REQUIRE(found_name);
            REQUIRE(gene_name == "DDX11L1");
        }
    }
    SECTION("find gene name")
    {
        std::string input = "gene_name \"DDX11L1\"";
        REQUIRE_FALSE(region.test_find_gene_info(input, gene_id, gene_name,
                                                 found_id, found_name));
        REQUIRE_FALSE(found_id);
        REQUIRE(found_name);
        REQUIRE(gene_name == "DDX11L1");
        SECTION("Find gene id")
        {
            input = "gene_id \"ENSG00000223972\"";
            REQUIRE(region.test_find_gene_info(input, gene_id, gene_name,
                                               found_id, found_name));
            REQUIRE(found_id);
            REQUIRE(gene_id == "ENSG00000223972");
            REQUIRE(found_name);
            REQUIRE(gene_name == "DDX11L1");
        }
    }
}

TEST_CASE("parse attribute")
{
    mock_region region;
    std::string gene_id, gene_name;
    SECTION("Found")
    {
        std::string input =
            "gene_id \"ENSG00000223972\"; gene_name \"DDX11L1\"; gene_source "
            "\"havana\"; gene_biotype \"transcribed_unprocessed_pseudogene\"";
        region.test_parse_attribute(input, gene_id, gene_name);
        REQUIRE(gene_id == "ENSG00000223972");
        REQUIRE(gene_name == "DDX11L1");
    }
    SECTION("not found")
    {
        std::string input = "gene_source \"havana\"; gene_biotype "
                            "\"transcribed_unprocessed_pseudogene\"";
        region.test_parse_attribute(input, gene_id, gene_name);
        REQUIRE(gene_id.empty());
        REQUIRE(gene_name.empty());
    }
}
