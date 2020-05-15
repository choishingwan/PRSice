#include "catch.hpp"
#include "mock_genotype.hpp"
#include "mock_prsice.hpp"
#include "prsice.hpp"
#include <algorithm>
#include <functional>
TEST_CASE("Simple Phenotype check")
{
    Reporter reporter("log", 60, true);
    Phenotype pheno;
    SECTION("Get Phenotype index")
    {
        std::string input = "FID IID SCZ MDD Bipolar SCZ";
        std::vector<std::string_view> col_sv = misc::tokenize(input);
        SECTION("With duplicated target phenotype column")
        {
            REQUIRE_THROWS(
                mock_prsice::test_get_pheno_idx(col_sv, pheno, "SCZ"));
        }
        SECTION("Phenotyep not found in header")
        {
            auto [idx, found] =
                mock_prsice::test_get_pheno_idx(col_sv, pheno, "ANX");
            REQUIRE(idx == 0);
            REQUIRE_FALSE(found);
        }
        SECTION("Phenotype found")
        {
            auto [idx, found] =
                mock_prsice::test_get_pheno_idx(col_sv, pheno, "MDD");
            REQUIRE(found);
            REQUIRE(idx == 3);
        }
    }
    SECTION("Parse phenotype file header")
    {
        PThresholding pthres;
        Permutations perm;
        CalculatePRS prs_info;
        Phenotype pheno;
        pheno.ignore_fid = GENERATE(true, false);
        prs_info.no_regress = false;
        SECTION("Invalid phenotype file")
        {
            auto invalid = GENERATE("\nFID IID Ok\n", "Invalid");
            auto input = std::make_unique<std::istringstream>(invalid);
            REQUIRE_THROWS(mock_prsice::test_parse_pheno_header(
                std::move(input), pheno, reporter));
        }
        SECTION("User did not provide a phenotype header")
        {
            using record = std::tuple<std::string, std::string>;
            auto headers = GENERATE(table<std::string, std::string>(
                {record {"FID IID SCZ MDD Bipolar SCZ", "SCZ"},
                 record {"FID IID 1 2 3", "Phenotype"}}));
            auto input =
                std::make_unique<std::istringstream>(std::get<0>(headers));
            Phenotype original = pheno;
            REQUIRE_NOTHROW(mock_prsice::test_parse_pheno_header(
                std::move(input), pheno, reporter));
            if (original.ignore_fid)
            {
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"IID"}));
            }
            else
            {
                REQUIRE_THAT(pheno.pheno_col, Catch::Equals<std::string>(
                                                  {std::get<1>(headers)}));
            }
            size_t exp_idx = !pheno.ignore_fid + 1;
            REQUIRE_THAT(pheno.pheno_col_idx, Catch::Equals<size_t>({exp_idx}));
        }
        SECTION("User provided multiple header")
        {
            auto input = std::make_unique<std::istringstream>(
                "FID IID SCZ MDD Bipolar SCZ");
            SECTION("Non of the column were found")
            {
                pheno.pheno_col = {"ALZ", "ADD", "ANX"};
                REQUIRE_THROWS(mock_prsice::test_parse_pheno_header(
                    std::move(input), pheno, reporter));
            }
            SECTION("Some column were found")
            {
                pheno.pheno_col = {"ALZ", "MDD"};
                REQUIRE_NOTHROW(mock_prsice::test_parse_pheno_header(
                    std::move(input), pheno, reporter));
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"ALZ", "MDD"}));
                REQUIRE_THAT(pheno.skip_pheno,
                             Catch::Equals<bool>({true, false}));
                REQUIRE_THAT(pheno.pheno_col_idx,
                             Catch::Equals<size_t>({0, 3}));
            }
        }
    }
    SECTION("Full pheno check")
    {
        // ignore fid don't have any effect on  when phenotype file isn't
        // provided
        pheno.ignore_fid = GENERATE(true, false);
        pheno.binary = {true};
        SECTION("Did not provide phenotype file")
        {
            auto no_regress = GENERATE(true, false);
            REQUIRE_NOTHROW(PRSice::pheno_check(no_regress, pheno, reporter));
            if (no_regress)
            {
                // with no regress, user are not required to provide binary
                // flag, thus we will just add one
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"PlaceHolder"}));
                REQUIRE_THAT(pheno.pheno_col_idx,
                             Catch::Equals<size_t>({~size_t(0)}));
                REQUIRE(pheno.binary.size() == pheno.pheno_col_idx.size());
                REQUIRE_THAT(pheno.skip_pheno, Catch::Equals<bool>({false}));
            }
            else
            {
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"Phenotype"}));
                REQUIRE_THAT(pheno.pheno_col_idx,
                             Catch::Equals<size_t>({~size_t(0)}));
                REQUIRE_THAT(pheno.skip_pheno, Catch::Equals<bool>({false}));
            }
        }
        SECTION("Provided a phenotype file")
        {
            std::ofstream pheno_file("pheno_check.prs");
            pheno_file << "FID IID SCZ MDD" << std::endl;
            pheno_file << "ID1 ID1 0 1" << std::endl;
            pheno_file.close();
            pheno.pheno_file = "pheno_check.prs";
            SECTION("Did not provide a header")
            {
                PRSice::pheno_check(false, pheno, reporter);
                REQUIRE_THAT(pheno.pheno_col_idx,
                             Catch::Equals<size_t>({1ul + !pheno.ignore_fid}));
                // because header isn't digit, we will use the column header as
                // the phenotype name
                if (pheno.ignore_fid)
                {
                    REQUIRE_THAT(pheno.pheno_col,
                                 Catch::Equals<std::string>({"IID"}));
                }
                else
                {
                    REQUIRE_THAT(pheno.pheno_col,
                                 Catch::Equals<std::string>({"SCZ"}));
                }
                REQUIRE_THAT(pheno.skip_pheno, Catch::Equals<bool>({false}));
            }
            SECTION("provided a header")
            {
                pheno.pheno_col = {"MDD"};
                PRSice::pheno_check(false, pheno, reporter);
                REQUIRE_THAT(pheno.pheno_col_idx, Catch::Equals<size_t>({3}));
                REQUIRE_THAT(pheno.pheno_col,
                             Catch::Equals<std::string>({"MDD"}));
                REQUIRE_THAT(pheno.skip_pheno, Catch::Equals<bool>({false}));
            }
        }
    }
}

TEST_CASE("load_pheno_map")
{
    std::string delim = " ";
    auto ignore_fid = GENERATE(true, false);
    Reporter reporter("log", 60, true);
    mock_prsice prsice;
    size_t idx = 2;
    SECTION("Without duplicated samples")
    {
        auto pheno_file = std::make_unique<std::istringstream>("FID IID Pheno\n"
                                                               "ID1 ID1 1\n"
                                                               "ID2 ID2 2\n"
                                                               "ID3 ID3 3\n"
                                                               "ID4 ID4 4\n"
                                                               "ID5 ID5 5\n"
                                                               "ID6 ID6 6\n");

        SECTION("Wrong idx size")
        {
            REQUIRE_THROWS(prsice.test_load_pheno_map(delim, 4, ignore_fid,
                                                      std::move(pheno_file)));
        }
        SECTION("Correct idx")
        {
            auto res = prsice.test_load_pheno_map(delim, idx, ignore_fid,
                                                  std::move(pheno_file));
            std::unordered_map<std::string, std::string> expected;
            for (size_t i = 0; i < 6; ++i)
            {
                std::string id = "ID" + std::to_string(i + 1);
                if (!ignore_fid) id = id + delim + id;
                expected[id] = std::to_string(i + 1);
            }
            // we also read in the header as we don't bother differentiating the
            // header line from the phenotype line
            auto header = "FID" + (ignore_fid ? "" : delim + "IID");
            expected[header] = "Pheno";
            REQUIRE(expected.size() == res.size());
            for (auto&& e : expected)
            {
                auto check = res.find(e.first);
                REQUIRE(check != res.end());
                REQUIRE(check->second == e.second);
            }
        }
    }
    SECTION("With duplicated FID")
    {

        // should be valid if we don't ignore FID
        auto pheno_file = std::make_unique<std::istringstream>("FID IID Pheno\n"
                                                               "ID1 ID1 1\n"
                                                               "ID2 ID2 2\n"
                                                               "ID3 ID3 3\n"
                                                               "ID4 ID4 4\n"
                                                               "ID5 ID5 5\n"
                                                               "ID3 ID6 6\n");
        if (ignore_fid)
        {
            REQUIRE_THROWS(prsice.test_load_pheno_map(delim, idx, ignore_fid,
                                                      std::move(pheno_file)));
        }
        else
        {
            // don't bother to do additional check
            REQUIRE_NOTHROW(prsice.test_load_pheno_map(delim, idx, ignore_fid,
                                                       std::move(pheno_file)));
        }
    }
}
TEST_CASE("test phenotype file processes")
{

    Reporter reporter("log", 60, true);
    auto ignore_fid = GENERATE(true, false);
    mockGenotype geno;
    SECTION("binary trait")
    {
        // two case, two control, one invalid pheno, one -9, one nan, one NA
        // all with regression = true
        std::vector<Sample_ID> samples {
            Sample_ID("Case1", "Case1", "1", true),
            Sample_ID("Case2", "Case2", "1", true),
            Sample_ID("Control1", "Control1", "0", true),
            Sample_ID("Case3", "Case3", "1", true),
            Sample_ID("Control2", "Control2", "0", true),
            Sample_ID("Missing", "Missing", "-9", true),
            Sample_ID("na1", "na1", "nan", true),
            Sample_ID("invalid", "invalid", "string", true),
            Sample_ID("na2", "na2", "NA", true)};

        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(samples.begin(), samples.end(), g);
        for (auto&& s : samples) { geno.add_sample(s); }
        geno.set_sample(samples.size());
        mock_prsice prsice(true, &reporter);
        SECTION("From phenotype file")
        {
            // also need one not found
            std::ofstream pheno_file("pheno_test");
            for (auto&& s : samples)
            {
                if (s.FID != "Control2")
                {
                    pheno_file << s.FID << "\t" << s.IID << "\t" << s.pheno
                               << std::endl;
                }
            }
            pheno_file << "NotFound\tNotFound\t1" << std::endl;
            pheno_file.close();
            auto [pheno_store, num_not_found, invalid, max_code] =
                prsice.test_process_phenotype_file("pheno_test", " ", 2,
                                                   ignore_fid, geno);
            REQUIRE(num_not_found == 1);
            REQUIRE(invalid == 1);
            REQUIRE(max_code == 1);
            auto pheno_map = prsice.sample_with_phenotypes();
            size_t valid_idx = 0;
            std::vector<double> expected;
            for (size_t i = 0; i < samples.size(); ++i)
            {
                if (samples[i].FID != "Missing" && samples[i].FID != "na1"
                    && samples[i].FID != "na2" && samples[i].FID != "invalid"
                    && samples[i].FID != "Control2")
                {
                    auto id = ignore_fid
                                  ? samples[i].IID
                                  : samples[i].FID + " " + samples[i].IID;
                    auto loc = pheno_map.find(id);
                    REQUIRE(loc != pheno_map.end());
                    REQUIRE(loc->second == valid_idx);
                    expected.push_back(misc::convert<double>(samples[i].pheno));
                    ++valid_idx;
                }
            }
            REQUIRE_THAT(pheno_store, Catch::Equals<double>(expected));
        }
        SECTION("Directly from fam")
        {
            auto [pheno_store, invalid, max_code] =
                prsice.test_process_phenotype_info(" ", ignore_fid, geno);
            REQUIRE(invalid == 1);
            REQUIRE(max_code == 1);
            auto pheno_map = prsice.sample_with_phenotypes();
            size_t valid_idx = 0;
            std::vector<double> expected;
            for (size_t i = 0; i < samples.size(); ++i)
            {
                if (samples[i].FID != "Missing" && samples[i].FID != "na1"
                    && samples[i].FID != "na2" && samples[i].FID != "invalid")
                {
                    auto id = ignore_fid
                                  ? samples[i].IID
                                  : samples[i].FID + " " + samples[i].IID;
                    auto loc = pheno_map.find(id);
                    REQUIRE(loc != pheno_map.end());
                    REQUIRE(loc->second == valid_idx);
                    expected.push_back(misc::convert<double>(samples[i].pheno));
                    ++valid_idx;
                }
            }
            REQUIRE_THAT(pheno_store, Catch::Equals<double>(expected));
        }
    }
    SECTION("quantitative trait")
    {

        // 4 valid , one invalid pheno, one -9, one nan, one NA
        // all with regression = true
        std::vector<Sample_ID> samples {
            Sample_ID("ID1", "ID1", "1", true),
            Sample_ID("ID2", "ID2", "2", true),
            Sample_ID("ID3", "ID3", "3", true),
            Sample_ID("ID4", "ID4", "4", true),
            Sample_ID("Missing", "Missing", "-9", true),
            Sample_ID("na1", "na1", "nan", true),
            Sample_ID("invalid", "invalid", "string", true),
            Sample_ID("na2", "na2", "NA", true)};
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(samples.begin(), samples.end(), g);
        for (auto&& s : samples) { geno.add_sample(s); }
        geno.set_sample(samples.size());

        mock_prsice prsice(false, &reporter);
        SECTION("Directly from fam")
        {
            auto [pheno_store, invalid, max_code] =
                prsice.test_process_phenotype_info(" ", ignore_fid, geno);
            REQUIRE(invalid == 1);
            // doesn't matter with max code
            auto pheno_map = prsice.sample_with_phenotypes();
            size_t valid_idx = 0;
            std::vector<double> expected;
            for (size_t i = 0; i < samples.size(); ++i)
            {
                if (samples[i].FID != "na1" && samples[i].FID != "na2"
                    && samples[i].FID != "invalid")
                {
                    auto id = ignore_fid
                                  ? samples[i].IID
                                  : samples[i].FID + " " + samples[i].IID;
                    auto loc = pheno_map.find(id);
                    REQUIRE(loc != pheno_map.end());
                    REQUIRE(loc->second == valid_idx);
                    expected.push_back(misc::convert<double>(samples[i].pheno));
                    ++valid_idx;
                }
            }
            REQUIRE_THAT(pheno_store, Catch::Equals<double>(expected));
        }

        SECTION("From phenotype file")
        {
            // also need one not found
            std::ofstream pheno_file("pheno_test");
            for (auto&& s : samples)
            {
                if (s.FID != "ID4")
                {
                    pheno_file << s.FID << "\t" << s.IID << "\t" << s.pheno
                               << std::endl;
                }
            }
            pheno_file << "NotFound\tNotFound\t1" << std::endl;
            pheno_file.close();
            auto [pheno_store, num_not_found, invalid, max_code] =
                prsice.test_process_phenotype_file("pheno_test", " ", 2,
                                                   ignore_fid, geno);
            REQUIRE(num_not_found == 1);
            REQUIRE(invalid == 1);
            auto pheno_map = prsice.sample_with_phenotypes();
            size_t valid_idx = 0;
            std::vector<double> expected;
            for (size_t i = 0; i < samples.size(); ++i)
            {
                if (samples[i].FID != "na1" && samples[i].FID != "na2"
                    && samples[i].FID != "invalid" && samples[i].FID != "ID4")
                {
                    auto id = ignore_fid
                                  ? samples[i].IID
                                  : samples[i].FID + " " + samples[i].IID;
                    auto loc = pheno_map.find(id);
                    REQUIRE(loc != pheno_map.end());
                    REQUIRE(loc->second == valid_idx);
                    expected.push_back(misc::convert<double>(samples[i].pheno));
                    ++valid_idx;
                }
            }
            REQUIRE_THAT(pheno_store, Catch::Equals<double>(expected));
        }
    }
}
TEST_CASE("Phenotype checking")
{
    Reporter reporter("log", 60, true);

    std::vector<double> pheno_store;
    SECTION("Quantitative trait")
    {
        mock_prsice prsice(false, &reporter);
        SECTION("One phenotype only")
        {
            pheno_store.resize(1000, -9);
            REQUIRE_FALSE(prsice.test_quantitative_pheno_is_valid(pheno_store));
        }
        SECTION("Valid")
        {
            std::random_device rnd_device;
            std::mt19937 mersenne_engine {rnd_device()};
            std::uniform_real_distribution<double> dist {-1.0, 1.0};

            auto gen = [&dist, &mersenne_engine]() {
                return dist(mersenne_engine);
            };
            pheno_store.resize(1000);
            std::generate(std::begin(pheno_store), std::end(pheno_store), gen);
            REQUIRE(prsice.test_quantitative_pheno_is_valid(pheno_store));
        }
    }
    SECTION("Binary trait")
    {
        mock_prsice prsice(false, &reporter);
        SECTION("Invalid pheno ")
        {
            // mixing 0 1 2
            std::random_device rnd_device;
            std::mt19937 mersenne_engine {rnd_device()};
            std::uniform_int_distribution<size_t> dist {0, 2};

            auto gen = [&dist, &mersenne_engine]() {
                return dist(mersenne_engine);
            };
            pheno_store.resize(1000);
            std::generate(std::begin(pheno_store), std::end(pheno_store), gen);
            auto [valid, ncase, ncontrol] =
                prsice.test_binary_pheno_is_valid(2, pheno_store);
            REQUIRE_FALSE(valid);
            // don't need to check case and control
        }
        SECTION("Valid")
        {
            // just 1 and 2 or 0 /1
            auto base = GENERATE(0ul, 1ul);
            std::random_device rnd_device;
            std::mt19937 mersenne_engine {rnd_device()};
            std::uniform_int_distribution<size_t> dist {0, 1};

            auto gen = [&dist, &mersenne_engine]() {
                return dist(mersenne_engine);
            };
            pheno_store.resize(1000);
            std::generate(std::begin(pheno_store), std::end(pheno_store), gen);
            std::vector<double> expected = pheno_store;
            auto expect_case =
                std::count(expected.begin(), expected.end(), 1.0);
            auto expect_control =
                std::count(expected.begin(), expected.end(), 0.0);
            if (base == 1)
            {
                for (auto&& p : pheno_store) ++p;
            }
            auto [valid, ncase, ncontrol] =
                prsice.test_binary_pheno_is_valid(1 + base, pheno_store);
            REQUIRE(valid);
            REQUIRE(ncase == Approx(expect_case));
            REQUIRE(ncontrol == Approx(expect_control));
            REQUIRE_THAT(pheno_store, Catch::Equals<double>(expected));
            // don't need to check case and control
        }
    }
}
TEST_CASE("Parse phenotype")
{
    std::vector<double> pheno_store;
    Reporter reporter("log", 60, true);
    int max_pheno = 0;
    SECTION("bianry trait")
    {
        // should have filtered out na nan and -9 before
        mock_prsice prsice(true, &reporter);
        SECTION("Invalid inputs")
        {
            std::string pheno = GENERATE("3", "String", "Case", "Control");
            REQUIRE_THROWS(
                prsice.test_parse_pheno(pheno, pheno_store, max_pheno));
        }
        SECTION("Valid inputs")
        {
            std::string pheno = GENERATE("0", "1", "2");
            REQUIRE_NOTHROW(
                prsice.test_parse_pheno(pheno, pheno_store, max_pheno));
            REQUIRE(max_pheno == pheno_store.back());
        }
    }
    SECTION("Quantitative trait")
    {
        mock_prsice prsice(false, &reporter);
        SECTION("Valid inputs")
        {
            std::string pheno = GENERATE("-1.96", "0", "1", "2");
            REQUIRE_NOTHROW(
                prsice.test_parse_pheno(pheno, pheno_store, max_pheno));
            REQUIRE(misc::convert<double>(pheno) == pheno_store.back());
        }
        SECTION("Invalid input")
        {
            std::string pheno = GENERATE("string", "nan");
            REQUIRE_THROWS(
                prsice.test_parse_pheno(pheno, pheno_store, max_pheno));
        }
    }
}
TEST_CASE("gen_pheno_vec")
{

    auto pheno_name = "Phenotype";
    const size_t sample_size = 100;
    std::vector<Sample_ID> samples;
    samples.reserve(sample_size);
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()};
    std::uniform_real_distribution<double> dist {-1.0, 1.0};

    auto gen = [&dist, &mersenne_engine]() { return dist(mersenne_engine); };
    for (size_t i = 0; i < sample_size; ++i)
    {
        samples.emplace_back(Sample_ID(std::to_string(i), std::to_string(i),
                                       std::to_string(gen()), true));
    }
    mockGenotype geno;
    std::shuffle(samples.begin(), samples.end(), mersenne_engine);
    for (auto&& s : samples) { geno.add_sample(s); }
    geno.set_sample(samples.size());
    Reporter reporter("log", 60, true);
    mock_prsice prsice(false, &reporter);
    SECTION("From file")
    {
        std::ofstream pheno_file("pheno_full_load");
        for (auto&& s : samples)
        {
            pheno_file << s.FID << "\t" << s.IID << "\t" << s.pheno
                       << std::endl;
        }
        pheno_file.close();
        REQUIRE_NOTHROW(prsice.test_gen_pheno_vec("pheno_full_load", pheno_name,
                                                  " ", 2, false, geno));
        Eigen::VectorXd res = prsice.phenotype_matrix();
        for (size_t i = 0; i < samples.size(); ++i)
        {
            REQUIRE(res(i, 0)
                    == Approx(misc::convert<double>(samples[i].pheno)));
        }
    }
    SECTION("From fam")
    {
        REQUIRE_NOTHROW(
            prsice.test_gen_pheno_vec("", pheno_name, " ", 2, false, geno));
        Eigen::VectorXd res = prsice.phenotype_matrix();
        for (size_t i = 0; i < samples.size(); ++i)
        {
            REQUIRE(res(i, 0)
                    == Approx(misc::convert<double>(samples[i].pheno)));
        }
    }
}
TEST_CASE("Print Phenotype log")
{
    auto ignore_fid = GENERATE(true, false);
    auto sample_ct = 1000ul;
    auto pheno_name = "Phenotype";

    Reporter reporter("log", 60, true);

    std::vector<double> pheno_store(sample_ct);
    SECTION("Empty pheno_store")
    {
        pheno_store.clear();
        mock_prsice prsice(true, &reporter);
        REQUIRE_THROWS_WITH(prsice.test_print_pheno_log(pheno_name, sample_ct,
                                                        0, 0, 1, ignore_fid,
                                                        pheno_store),
                            Catch::Contains("No phenotype presented"));
    }
    SECTION("All not found")
    {
        std::string expected_message =
            ignore_fid
                ? "Maybe the first column of your phenotype file "
                  "is the FID?"
                : "Maybe your phenotype file doesn not contain the FID?\n"
                  "Might want to consider using --ignore-fid\n";
        mock_prsice prsice(true, &reporter);
        REQUIRE_THROWS_WITH(
            prsice.test_print_pheno_log(pheno_name, sample_ct, sample_ct, 0, 1,
                                        ignore_fid, pheno_store),
            Catch::Contains("Error: No sample left")
                && Catch::Contains(expected_message));
    }
    SECTION("All invalid")
    {
        mock_prsice prsice(true, &reporter);
        REQUIRE_THROWS_WITH(
            prsice.test_print_pheno_log(pheno_name, sample_ct, 0, sample_ct, 1,
                                        ignore_fid, pheno_store),
            Catch::Contains("Error: No sample left")
                && Catch::Contains("All sample has invalid phenotypes"));
    }
    SECTION("Problematic Binary trait ")
    {
        mock_prsice prsice(true, &reporter);
        SECTION("Singluar situation")
        {
            using record = std::tuple<size_t, std::string>;
            auto input = GENERATE(table<size_t, std::string>(
                {record {0, "cases"}, record {1, "control samples"}}));
            std::fill(pheno_store.begin(), pheno_store.end(),
                      std::get<0>(input));
            REQUIRE_THROWS_WITH(
                prsice.test_print_pheno_log(pheno_name, sample_ct, 0, 0, 1,
                                            ignore_fid, pheno_store),
                Catch::Contains("There are no " + std::get<1>(input)));
        }
        SECTION("Mixed encoding")
        {
            std::random_device rnd_device;
            std::mt19937 mersenne_engine {rnd_device()};
            std::uniform_int_distribution<size_t> dist {0, 2};
            auto gen = [&dist, &mersenne_engine]() {
                return dist(mersenne_engine);
            };
            std::generate(std::begin(pheno_store), std::end(pheno_store), gen);
            int max_pheno_code = GENERATE(1, 2);
            REQUIRE_THROWS_WITH(prsice.test_print_pheno_log(
                                    pheno_name, sample_ct, 0, 0, max_pheno_code,
                                    ignore_fid, pheno_store),
                                Catch::Contains("Mixed encoding"));
        }
    }
    SECTION("Single Quantitative trait")
    {
        mock_prsice prsice(false, &reporter);
        SECTION("All -9")
        {
            std::fill(pheno_store.begin(), pheno_store.end(), -9);
            REQUIRE_THROWS_WITH(
                prsice.test_print_pheno_log(pheno_name, sample_ct, 0, 0, 2,
                                            ignore_fid, pheno_store),
                Catch::Contains("they are all -9")
                    && Catch::Contains("Not enough valid phenotype"));
        }
        SECTION("Not -9")
        {
            std::fill(pheno_store.begin(), pheno_store.end(), 1.96);
            REQUIRE_THROWS_WITH(
                prsice.test_print_pheno_log(pheno_name, sample_ct, 0, 0, 2,
                                            ignore_fid, pheno_store),
                !Catch::Contains("they are all -9")
                    && Catch::Contains("Not enough valid phenotype"));
        }
    }
}
