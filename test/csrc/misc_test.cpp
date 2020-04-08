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
