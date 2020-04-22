#include "catch.hpp"
#include "misc.hpp"


TEST_CASE("Convertor")
{
    SECTION("size_t")
    {
        SECTION("valid")
        {
            auto i = GENERATE(take(100, random(-10000, 10000)));
            if (i < 0)
            {
                REQUIRE_THROWS(
                    misc::Convertor::convert<size_t>(std::to_string(i)));
            }
            else
            {
                REQUIRE(misc::Convertor::convert<size_t>(std::to_string(i))
                        == static_cast<size_t>(i));
            }
        }
    }
    SECTION("int")
    {
        auto i = GENERATE(take(100, random(-10000, 10000)));
        REQUIRE(misc::Convertor::convert<int>(std::to_string(i)) == i);
    }
    // with double there will be float point exception so we can't use generator
    // directly
    SECTION("double")
    {
        using record = std::tuple<std::string, double>;
        auto extent = GENERATE(table<std::string, double>(
            {record {"0.5", 0.5}, record {"1e-10", 1e-10},
             record {"-1e20", -1e20}}));
        REQUIRE(misc::Convertor::convert<double>(std::get<0>(extent))
                == Approx(std::get<1>(extent)));
    }
    SECTION("double overflow")
    {
        REQUIRE_THROWS(misc::Convertor::convert<double>("1e-400"));
        REQUIRE_THROWS(misc::Convertor::convert<double>("1e400"));
    }
}
TEST_CASE("stringview trimming")
{
    std::string ref = " testing \n";
    std::string_view s = ref;
    std::string_view t = ref;
    s = misc::trimmed(s);
    REQUIRE(s == "testing");
    REQUIRE(ref == " testing \n");
    misc::trim(t);
    REQUIRE(t == "testing");
}


TEST_CASE("hasEnding")
{
    std::string full_text;
    std::string ending;
    SECTION("empty text")
    {
        REQUIRE_THROWS(misc::hasEnding(full_text, ending));
    }
    full_text = "Hello";
    SECTION("empty ending")
    {
        REQUIRE_THROWS(misc::hasEnding(full_text, ending));
    }
    SECTION("found ending")
    {
        auto i = GENERATE("ello", "lo", "o", "Hello");
        REQUIRE(misc::hasEnding(full_text, i));
    }
    SECTION("ending not found")
    {
        auto i = GENERATE("Lo", "O", "aHello");
        REQUIRE_FALSE(misc::hasEnding(full_text, i));
    }
}


TEST_CASE("split")
{
    SECTION("stringview")
    {
        std::string input = "test,the,splitter,and,make,sure-it,works well";
        std::vector<std::string_view> token = misc::tokenize(input, ",");
        REQUIRE_THAT(token, Catch::Equals<std::string_view>(
                                {"test", "the", "splitter", "and", "make",
                                 "sure-it", "works well"}));
        input = "sure-it\tworks well ok";
        token = misc::tokenize(input);
        REQUIRE_THAT(token, Catch::Equals<std::string_view>(
                                {"sure-it", "works", "well", "ok"}));
        /*
        input = "what\tif\twe\tgot\tempty\t\tinput";
        token = misc::tokenize(input);
        REQUIRE_THAT(token,
                     Catch::Equals<std::string_view>(
                         {"what", "if", "we", "got", "empty", "", "input"}));
        input = "back\tempty\t";
        token = misc::tokenize(input);
        REQUIRE_THAT(token,
                     Catch::Equals<std::string_view>({"back", "empty", ""}));
        input = "\tfront\tempty";
        token = misc::tokenize(input);
        REQUIRE_THAT(token,
                     Catch::Equals<std::string_view>({"", "front", "empty"}));*/
    }
    SECTION("string")
    {
        std::string input = "test,the,splitter,and,make,sure-it,works well";
        std::vector<std::string> token = misc::split(input, ",");
        std::vector<std::string> alt;
        misc::split(alt, input, ",");
        REQUIRE_THAT(token, Catch::Equals<std::string>(
                                {"test", "the", "splitter", "and", "make",
                                 "sure-it", "works well"}));
        REQUIRE_THAT(
            alt, Catch::Equals<std::string>({"test", "the", "splitter", "and",
                                             "make", "sure-it", "works well"}));
        input = "sure-it\tworks well ok";
        token = misc::split(input);
        misc::split(alt, input);
        REQUIRE_THAT(token, Catch::Equals<std::string>(
                                {"sure-it", "works", "well", "ok"}));
        REQUIRE_THAT(alt, Catch::Equals<std::string>(
                              {"sure-it", "works", "well", "ok"}));
        /*
        input = "what\tif\twe\tgot\tempty\t\tinput";
        token = misc::split(input);
        misc::split(alt, input);
        REQUIRE_THAT(token,
                     Catch::Equals<std::string>(
                         {"what", "if", "we", "got", "empty", "", "input"}));
        REQUIRE_THAT(alt, Catch::Equals<std::string>({"what", "if", "we", "got",
                                                      "empty", "", "input"}));
        input = "back\tempty\t";
        token = misc::split(input);
        misc::split(alt, input);
        REQUIRE_THAT(token, Catch::Equals<std::string>({"back", "empty", ""}));
        REQUIRE_THAT(alt, Catch::Equals<std::string>({"back", "empty", ""}));
        input = "\tfront\tempty";
        token = misc::split(input);
        misc::split(alt, input);
        REQUIRE_THAT(token, Catch::Equals<std::string>({"", "front", "empty"}));
        REQUIRE_THAT(alt, Catch::Equals<std::string>({"", "front", "empty"}));*/
    }
}
