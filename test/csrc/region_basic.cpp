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
        std::uniform_int_distribution<size_t> rand_idx {0, n_features};
        for (size_t i = 0; i < 10; ++i)
        {
            REQUIRE(region.test_in_region(feature[rand_idx(mersenne_engine)],
                                          feature));
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
}
