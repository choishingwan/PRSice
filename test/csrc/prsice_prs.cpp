#include "catch.hpp"
#include "mock_prsice.hpp"
#include "prsice.hpp"
#include "storage.hpp"

TEST_CASE("Initialize progress bar")
{
    // assume 10 sets
    std::random_device rd;
    std::mt19937 mersenne_engine(rd());
    std::uniform_real_distribution<double> statistic(0.0, 1.0);
    std::uniform_int_distribution<size_t> num_thres(0, 10);
    const size_t num_set = 10;
    std::vector<std::set<double>> sets(num_set);
    // generate random number of thresholds
    size_t total_thres = 0;
    for (size_t i_set = 0; i_set < num_set; ++i_set)
    {
        auto thres = num_thres(mersenne_engine);
        auto start_stat = statistic(mersenne_engine);
        for (size_t i = 0; i < thres; ++i)
        { sets[i_set].insert(start_stat + i); }
        if (i_set != 1) total_thres += thres;
    }
    CalculatePRS prs_info;
    PThresholding thres;
    Phenotype pheno;
    Permutations perm;
    using record = std::tuple<bool, bool>;
    auto perm_type = GENERATE(table<bool, bool>(
        {record {false, false}, record {true, false}, record {false, true}}));
    perm.run_perm = std::get<0>(perm_type);
    perm.run_set_perm = std::get<1>(perm_type);
    perm.num_permutation = GENERATE(take(1, random(1ul, 10000ul)));
    std::string output = "PRSice";
    Reporter reporter("log", 60, true);
    mock_prsice prsice(prs_info, thres, perm, output, true, &reporter);
    size_t exp_total = 0;
    if (perm.run_perm)
        exp_total = total_thres * (perm.num_permutation + 1);
    else
        exp_total = total_thres;
    size_t exp_competitive = 0;
    if (perm.run_set_perm)
    { exp_competitive = perm.num_permutation * (num_set - 2); }
    auto analysis_done = GENERATE(take(1, random(1ul, 1000ul)));
    auto competitive_done = GENERATE(take(1, random(1ul, 1000ul)));
    prsice.set_progress(analysis_done, competitive_done);
    auto [ad, cd] = prsice.get_current_progress();
    REQUIRE(ad == analysis_done);
    REQUIRE(cd == competitive_done);
    prsice.init_progress_count(sets);
    auto [total, com_total] = prsice.get_progress();
    REQUIRE(total == exp_total);
    REQUIRE(com_total == exp_competitive);
    std::tie(ad, cd) = prsice.get_current_progress();
    REQUIRE(ad == 0);
    REQUIRE(cd == 0);
}
