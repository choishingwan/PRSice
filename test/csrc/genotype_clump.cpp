#include "catch.hpp"
#include "genotype.hpp"
#include "mock_genotype.hpp"

TEST_CASE("Sort by p")
{
    mockGenotype geno;
    Reporter reporter("log", 60, true);
    geno.set_reporter(&reporter);
    std::vector<SNP> input = {SNP("rs4567", 1, 1000, "A", "C", 0, 0.06, 0, 0),
                              SNP("rs3456", 2, 1000, "A", "C", 0, 0.05, 0, 0),
                              SNP("rs3455", 2, 1000, "A", "C", 0, 0.06, 0, 0),
                              SNP("rs3452", 2, 1003, "A", "C", 0, 0.06, 0, 0),
                              SNP("rs3457", 2, 1003, "A", "C", 0, 0.06, 0, 0)

    };
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(input.begin(), input.end(), g);
    for (auto&& snp : input) { geno.load_snp(snp); }
    REQUIRE_NOTHROW(geno.sort_by_p());
    auto idx = geno.sorted_p_index();
    std::vector<std::string> expected_order = {"rs4567", "rs3456", "rs3455",
                                               "rs3452", "rs3457"};
    auto res = geno.existed_snps();
    REQUIRE(res.size() == expected_order.size());
    for (size_t i = 0; i < res.size(); ++i)
    { REQUIRE(res[idx[i]].rs() == expected_order[i]); }
}
TEST_CASE("Build Clump window")
{
    mockGenotype geno;
    Reporter reporter("log", 60, true);
    geno.set_reporter(&reporter);
    std::vector<SNP> input = {SNP("rs1", 1, 10, "A", "C", 0, 10),
                              SNP("rs5", 1, 10, "A", "C", 0, 13),
                              SNP("rs3", 1, 10, "A", "C", 1, 10),
                              SNP("rs4", 1, 20, "A", "C", 1, 10),
                              SNP("rs10", 1, 25, "A", "C", 1, 10),
                              SNP("rs6", 1, 40, "A", "C", 1, 10),
                              SNP("rs7", 1, 70, "A", "C", 1, 10),
                              SNP("rs8", 3, 70, "A", "C", 1, 10)};
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(input.begin(), input.end(), g);
    for (auto&& snp : input) { geno.load_snp(snp); }
    REQUIRE_NOTHROW(geno.build_clump_windows(10));
    auto res = geno.existed_snps();
    std::vector<std::string> expected_order = {"rs1",  "rs5", "rs3", "rs4",
                                               "rs10", "rs6", "rs7", "rs8"};
    // window size only need to cover # SNP before
    REQUIRE(geno.max_window() == 4);
    REQUIRE(res.size() == expected_order.size());
    for (size_t i = 0; i < res.size(); ++i)
    {
        REQUIRE(res[i].rs() == expected_order[i]);
        switch (i)
        {
        case 0:
        case 1:
        case 2:
            REQUIRE(res[i].up_bound() == 4);
            REQUIRE(res[i].low_bound() == 0);
            break;
        case 3:
            REQUIRE(res[i].up_bound() == 5);
            REQUIRE(res[i].low_bound() == 0);
            break;
        case 4:
            REQUIRE(res[i].up_bound() == 5);
            REQUIRE(res[i].low_bound() == 3);
            break;
        case 5:
            REQUIRE(res[i].up_bound() == 6);
            REQUIRE(res[i].low_bound() == 5);
            break;
        case 6:
            REQUIRE(res[i].up_bound() == 7);
            REQUIRE(res[i].low_bound() == 6);
            break;
        case 7:
            REQUIRE(res[i].up_bound() == 8);
            REQUIRE(res[i].low_bound() == 7);
            break;
        }
    }
}

// Things that we need to test
// Read_Genotype function for both binaryplink and binarygen
// update_index_tot function (might need to use plink and predefined data)
// get_r2 function (again, might need to use predefined data and reference to
// plink)
// clump from SNP
