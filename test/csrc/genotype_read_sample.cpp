#include "catch.hpp"
#include "genotype.hpp"
#include "mock_genotype.hpp"

TEST_CASE("initialize sample vectors")
{
    mockGenotype geno;
    SECTION("without any samples")
    {
        REQUIRE_THROWS(geno.init_sample_vectors());
    }
    SECTION("with samples")
    {
        auto n_sample = GENERATE(as<uintptr_t> {}, 1023, 1025);
        geno.set_sample(n_sample);
        REQUIRE_NOTHROW(geno.init_sample_vectors());
    }
}
