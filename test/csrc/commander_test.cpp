#include "catch.hpp"
#include "commander.hpp"
#include "mock_commander.hpp"

void transform_test(const std::string& cov,
                    const std::vector<std::string>& expected,
                    bool expect_success)
{
    std::vector<std::string> result;
    if (!expect_success)
    { REQUIRE_FALSE(mockCommander::transform_wrapper(cov, result)); }
    else
    {
        REQUIRE(mockCommander::transform_wrapper(cov, result));
        REQUIRE_THAT(result, Catch::UnorderedEquals<std::string>(expected));
    }
}

TEST_CASE("Covariate Transformation")
{
    SECTION("find_first_end")
    {
        std::string cov = "PC[1-5]";
        size_t res;
        SECTION("wrong start index")
        {
            REQUIRE_FALSE(mockCommander::find_first_end_wrapper(cov, 0, res));
        }
        SECTION("correct start")
        {
            REQUIRE(mockCommander::find_first_end_wrapper(cov, 2, res));
        }
        SECTION("nested list not allowed")
        {
            cov = "PC[1-5[1-5]]";
            SECTION("wrong start index")
            {
                REQUIRE_FALSE(
                    mockCommander::find_first_end_wrapper(cov, 0, res));
            }
            SECTION("correct start")
            {
                REQUIRE_FALSE(
                    mockCommander::find_first_end_wrapper(cov, 2, res));
            }
        }
    }
    SECTION("parse_range")
    {
        std::vector<size_t> res;
        SECTION("Expect [] removed")
        {
            REQUIRE_FALSE(mockCommander::parse_range_wrapper("[1-5]", res));
        }
        SECTION("valid ranges")
        {
            int start = GENERATE(take(5, random(1, 1000)));
            int end = GENERATE(take(5, random(1, 1000)));
            REQUIRE(mockCommander::parse_range_wrapper(
                std::to_string(start) + "-" + std::to_string(end), res));
            if (start > end) { std::swap(start, end); }
            std::vector<size_t> expected(end - start + 1, start);
            std::iota(expected.begin(), expected.end(), start);
            REQUIRE_THAT(res, Catch::Equals<size_t>(expected));
        }
        SECTION("start == end")
        {
            REQUIRE(mockCommander::parse_range_wrapper("1-1", res));
            REQUIRE_THAT(res, Catch::Equals<size_t>({1}));
        }
        SECTION("single value")
        {
            auto i = GENERATE(take(5, random(-100, 100)));
            if (i < 0)
            {
                REQUIRE_FALSE(
                    mockCommander::parse_range_wrapper(std::to_string(i), res));
            }
            else
            {
                REQUIRE(
                    mockCommander::parse_range_wrapper(std::to_string(i), res));
                REQUIRE_THAT(res,
                             Catch::Equals<size_t>({static_cast<size_t>(i)}));
            }
        }
        SECTION("not allow negative range")
        {
            REQUIRE_FALSE(mockCommander::parse_range_wrapper("1--5", res));
        }
        SECTION("not allow full negative")
        {
            REQUIRE_FALSE(mockCommander::parse_range_wrapper("-1--5", res));
        }
        SECTION("not allow comma")
        {
            REQUIRE_FALSE(mockCommander::parse_range_wrapper("1,5", res));
        }
    }
    SECTION("get_range")
    {
        std::vector<size_t> res;
        SECTION("invalid start")
        {
            REQUIRE_FALSE(
                mockCommander::get_range_wrapper("PC[1-5]", 0, 6, res));
        }
        SECTION("invalid end")
        {
            REQUIRE_FALSE(
                mockCommander::get_range_wrapper("PC[1-5]", 2, 4, res));
        }
        SECTION("out of bound end")
        {
            REQUIRE_FALSE(
                mockCommander::get_range_wrapper("PC[1-5]", 2, 7, res));
        }
        SECTION("Correct parse")
        {
            REQUIRE(mockCommander::get_range_wrapper("PC[1-5]", 2, 6, res));
            REQUIRE_THAT(res, Catch::Equals<size_t>({1, 2, 3, 4, 5}));
        }
        SECTION("mixed parse")
        {
            REQUIRE(
                mockCommander::get_range_wrapper("PC[1-5.8.7-10]", 2, 13, res));
            REQUIRE_THAT(res, Catch::UnorderedEquals<size_t>(
                                  {1, 2, 3, 4, 5, 7, 8, 9, 10}));
        }
        SECTION("double parse")
        {
            // should be disallowed, but can't detect, need to check log
            REQUIRE(mockCommander::get_range_wrapper("PC[1-5.6.0.005]", 2, 14,
                                                     res));
            REQUIRE_THAT(res,
                         Catch::UnorderedEquals<size_t>({0, 1, 2, 3, 4, 5, 6}));
        }
    }
    SECTION("update_covariate_ranges")
    {
        std::vector<std::string> result;
        std::vector<size_t> range;
        SECTION("empty vectors")
        {
            REQUIRE_FALSE(
                mockCommander::update_covariate_ranges_wrapper(result, range));
        }
        range = {1, 3, 5, 7, 9};
        SECTION("valid range")
        {
            REQUIRE(
                mockCommander::update_covariate_ranges_wrapper(result, range));
            std::vector<std::string> expected = {"1", "3", "5", "7", "9"};
            REQUIRE_THAT(result, Catch::UnorderedEquals<std::string>(expected));
        }
        SECTION("append range")
        {
            result = {"PC1AB", "PC2CD"};
            REQUIRE(
                mockCommander::update_covariate_ranges_wrapper(result, range));
            std::vector<std::string> expected = {
                "PC1AB1", "PC1AB3", "PC1AB5", "PC1AB7", "PC1AB9",
                "PC2CD1", "PC2CD3", "PC2CD5", "PC2CD7", "PC2CD9"};
            REQUIRE_THAT(result, Catch::UnorderedEquals<std::string>(expected));
        }
    }
    SECTION("transformation")
    {
        SECTION("without parsing")
        {
            transform_test("@PC", std::vector<std::string> {"PC"}, true);
        }
        SECTION("double @")
        {
            transform_test("@@PC", std::vector<std::string> {"@PC"}, true);
        }
        SECTION("without @")
        {
            // this is ok, as we'd removed the @ when we start transformation
            transform_test("PC1", std::vector<std::string> {"PC1"}, true);
        }
        SECTION("Range")
        {
            transform_test("PC1-5", std::vector<std::string> {"PC1-5"}, true);
        }
        SECTION("Internal @")
        {
            transform_test("PC1@5", std::vector<std::string> {"PC1@5"}, true);
        }
        SECTION("Without []")
        {
            // we will always remove first @
            transform_test("@PC1-5", std::vector<std::string> {"PC1-5"}, true);
        }
        SECTION("Normal transformation")
        {
            transform_test(
                "@PC[1-5]",
                std::vector<std::string> {"PC1", "PC2", "PC3", "PC4", "PC5"},
                true);
        }
        SECTION("Mixed tranform")
        {
            transform_test("@PC[1-2.5]",
                           std::vector<std::string> {"PC1", "PC2", "PC5"},
                           true);
        }
        SECTION("Complex transform")
        {
            transform_test("@PC[1-2.4.3-6]",
                           std::vector<std::string> {"PC1", "PC2", "PC3", "PC4",
                                                     "PC5", "PC6"},
                           true);
        }
        SECTION("Parse internal")
        {
            transform_test("@PC[1-2]A",
                           std::vector<std::string> {"PC1A", "PC2A"}, true);
        }
        SECTION("Multiple range")
        {
            transform_test(
                "@PC[1-2]A[1-2]",
                std::vector<std::string> {"PC1A1", "PC1A2", "PC2A1", "PC2A2"},
                true);
        }
    }
}

TEST_CASE("Unit parsing")
{
    mockCommander commander;
    size_t value = 0;
    SECTION("memory without number")
    {
        auto i = GENERATE("m", "b", "mb", "hi", "TB");
        REQUIRE_FALSE(commander.check_parse_unit_value(i, "--mem", 0, value));
    }
    SECTION("out of bound")
    {
        SECTION("high power")
        {
            REQUIRE_FALSE(commander.check_parse_unit_value("1", "", 7, value));
        }
        SECTION("negative input")
        {
            REQUIRE_FALSE(
                commander.check_parse_unit_value("1000000000tb", "", 1, value));
        }
    }
    SECTION("negative input")
    {
        SECTION("without unit")
        {
            REQUIRE_FALSE(commander.check_parse_unit_value("-1", "", 1, value));
        }
        SECTION("with unit")
        {
            REQUIRE_FALSE(
                commander.check_parse_unit_value("-1tb", "", 1, value));
        }
    }
    SECTION("unit provided")
    {
        SECTION("ignore default")
        {
            auto i = GENERATE(range(0, 7));
            REQUIRE(commander.check_parse_unit_value(
                "1b", "", static_cast<size_t>(i), value));
            REQUIRE(value == 1);
        }
        SECTION("different units")
        {
            using record = std::tuple<std::string, size_t>;
            auto expected = GENERATE(table<std::string, size_t>(
                {record {"1k", 1000}, record {"1kb", 1000},
                 record {"1m", 1000000}, record {"1mb", 1000000},
                 record {"1g", 1000000000}, record {"1gb", 1000000000},
                 record {"1t", 1000000000000}, record {"1tb", 1000000000000}}));
            auto input = std::get<0>(expected);
            auto res = std::get<1>(expected);
            REQUIRE(commander.check_parse_unit_value(input, "", 1, value));
            REQUIRE(value == res);
        }
        SECTION("non-integer")
        {
            using record = std::tuple<std::string, size_t>;
            auto expected = GENERATE(table<std::string, size_t>(
                {record {"1.5k", 1500}, record {"1.004kb", 1004}}));
            auto input = std::get<0>(expected);
            auto res = std::get<1>(expected);
            REQUIRE(commander.check_parse_unit_value(input, "", 1, value));
            REQUIRE(value == res);
            SECTION("result can't be double")
            {
                REQUIRE_FALSE(
                    commander.check_parse_unit_value("1.5b", "", 1, value));
            }
        }
    }
    SECTION("default value")
    {
        auto power = GENERATE(range(0, 6));
        auto exp = std::pow(1000, power);
        REQUIRE(commander.check_parse_unit_value(
            "1", "", static_cast<size_t>(power), value));
        REQUIRE(value == static_cast<size_t>(exp));
    }
}
