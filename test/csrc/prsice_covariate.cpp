#include "catch.hpp"
#include "mock_genotype.hpp"
#include "mock_prsice.hpp"

TEST_CASE("Check if covariate is valid")
{
    std::vector<std::string> covariates = {"ID1", "ID1", "1",     "1", "str",
                                           "1",   "nan", "TypeA", "1", "NA",
                                           "1",   "1",   "EUR"};
    std::vector<size_t> missing_count;
    SECTION("All valid")
    {
        std::vector<size_t> cov_idx = {2, 3, 5, 7, 11, 12};
        std::set<size_t> factor_idx = {7, 12};
        missing_count.resize(cov_idx.size(), 0);
        REQUIRE(mock_prsice::test_is_valid_covariate(
            factor_idx, cov_idx, covariates, missing_count));
        REQUIRE_THAT(missing_count, Catch::Equals<size_t>(std::vector<size_t>(
                                        cov_idx.size(), 0)));
    }
    SECTION("With missing or invalid")
    {
        auto invalid_idx = GENERATE(4ul, 6ul, 9ul);
        std::vector<size_t> cov_idx = {2, 3, 5, 7, 11, 12, invalid_idx};
        std::sort(cov_idx.begin(), cov_idx.end());
        missing_count.resize(cov_idx.size(), 0);
        std::set<size_t> factor_idx = {7, 12};
        REQUIRE_FALSE(mock_prsice::test_is_valid_covariate(
            factor_idx, cov_idx, covariates, missing_count));
        std::vector<size_t> expected(cov_idx.size(), 0);
        for (size_t i = 0; i < expected.size(); ++i)
        {
            if (cov_idx[i] == invalid_idx) { ++expected[i]; }
        }
        REQUIRE_THAT(missing_count, Catch::Equals<size_t>(expected));
    }
    SECTION("With multiple problem")
    {
        std::vector<size_t> cov_idx = {2, 3, 4, 5, 6, 7, 11, 12};
        std::set<size_t> factor_idx = {7, 12};
        missing_count.resize(cov_idx.size(), 0);
        REQUIRE_FALSE(mock_prsice::test_is_valid_covariate(
            factor_idx, cov_idx, covariates, missing_count));
        std::vector<size_t> expected(cov_idx.size(), 0);
        expected[2] = 1;
        expected[4] = 1;
        REQUIRE_THAT(missing_count, Catch::Equals<size_t>(expected));
    }
}

TEST_CASE("covarience check and factor level count")
{
    // check detect duplicated ID (shouldn't check if no valid pheno)
    // check resulting levels
    // proper termination if no valid samples
    // need to initialize m_sample_with_phenotypes
    const size_t num_sample = 1000;
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()};
    // select ~70% of samples
    std::uniform_int_distribution<size_t> dist {1, 10};
    std::uniform_int_distribution<size_t> sex_dist {0, 1};
    std::uniform_int_distribution<size_t> batch_dist {0, 100};
    auto in_reg = [&dist, &mersenne_engine]() {
        return dist(mersenne_engine) > 3;
    };
    auto batch = [&batch_dist, &mersenne_engine] {
        return "b" + std::to_string(batch_dist(mersenne_engine));
    };
    auto sex = [&sex_dist, &mersenne_engine] {
        switch (sex_dist(mersenne_engine))
        {
        case 0: return "F";
        case 1: return "M";
        default: return "NA";
        }
    };
    std::normal_distribution pheno_dist;
    auto pheno = [&pheno_dist, &mersenne_engine]() {
        return pheno_dist(mersenne_engine);
    };
    mockGenotype geno;
    Reporter reporter("log", 60, true);
    mock_prsice prsice;
    prsice.set_reporter(&reporter);
    geno.set_reporter(&reporter);
    std::vector<double> sample_pheno;
    std::vector<bool> sample_in_regression(num_sample, false);
    auto&& sample_with_pheno = prsice.sample_with_phenotypes();
    const std::string delim = " ";
    auto ignore_fid = GENERATE(true, false);
    for (size_t i = 0; i < num_sample; ++i)
    {
        // auto cur_pheno = pheno();
        auto cur_pheno = i;
        auto in_regression = in_reg();
        geno.add_sample(Sample_ID(std::to_string(i), std::to_string(i),
                                  std::to_string(cur_pheno), in_regression));
        geno.update_valid_sample(i, in_regression);
        if (in_regression)
        {
            sample_in_regression[i] = true;
            std::string id = std::to_string(i);
            if (!ignore_fid) id.append(delim + std::to_string(i));
            sample_with_pheno[id] = i;
            sample_pheno.push_back(
                misc::convert<double>(std::to_string(cur_pheno)));
        }
    }
    prsice.phenotype_matrix() = Eigen::Map<Eigen::VectorXd>(
        sample_pheno.data(), static_cast<Eigen::Index>(sample_pheno.size()));
    std::vector<std::string> cov_names = {"PC1", "PC2", "Sex", "Age", "Batch"};
    std::vector<size_t> cov_idx = {2, 4, 5, 6, 7};
    std::set<size_t> factor_idx = {5, 7};
    std::string cov_header = "FID IID PC1 Something PC2 Sex Age Batch Centre\n";
    SECTION("Valid input")
    {
        const size_t sample_in_file = 2000;
        std::string cov_file = cov_header;
        std::unordered_map<std::string, size_t> expected_batch, expected_sex;
        size_t batch_count = 0, sex_count = 0;
        std::vector<double> expected_pheno;
        auto&& sample_vec = geno.get_sample_vec();
        for (size_t i = 0; i < sample_in_file; ++i)
        {
            if (i >= num_sample || !sample_vec[i].in_regression)
            {
                // doesn't matter what we sim
                cov_file.append(std::to_string(i) + " " + std::to_string(i)
                                + " " + std::to_string(pheno()) + " "
                                + std::to_string(pheno()) + " "
                                + std::to_string(pheno()) + " " + sex() + " "
                                + std::to_string(pheno()) + " " + batch() + " "
                                + batch() + "\n");
            }
            else
            {
                auto valid_cov = in_reg();
                auto pc1 = valid_cov ? std::to_string(pheno()) : "NA";
                auto cur_batch = batch();
                auto cur_sex = sex();
                if (valid_cov)
                {
                    expected_pheno.push_back(
                        misc::convert<double>(sample_vec[i].pheno));
                    if (expected_sex.find(cur_sex) == expected_sex.end())
                    { expected_sex[cur_sex] = sex_count++; }
                    if (expected_batch.find(cur_batch) == expected_batch.end())
                    { expected_batch[cur_batch] = batch_count++; }
                    if (geno.sample_valid_for_regress(i))
                    { geno.update_valid_sample(i, true); }
                }
                cov_file.append(std::to_string(i) + " " + std::to_string(i)
                                + " " + pc1 + " " + std::to_string(pheno())
                                + " " + std::to_string(pheno()) + " " + cur_sex
                                + " " + std::to_string(pheno()) + " "
                                + cur_batch + " " + batch() + "\n");
            }
        }
        std::unique_ptr<std::istream> input_file =
            std::make_unique<std::istringstream>(cov_file);

        auto res = prsice.test_cov_check_and_factor_level_count(
            factor_idx, cov_names, cov_idx, delim, ignore_fid, input_file,
            geno);
        auto sex_factor = res.front();
        auto batch_factor = res.back();
        REQUIRE(sex_factor.size() == expected_sex.size());
        REQUIRE(batch_factor.size() == expected_batch.size());
        for (auto&& b : expected_batch)
        {
            REQUIRE(batch_factor.find(b.first) != batch_factor.end());
            REQUIRE(batch_factor[b.first] == b.second);
        }
        for (auto&& s : expected_sex)
        {
            REQUIRE(sex_factor.find(s.first) != sex_factor.end());
            REQUIRE(sex_factor[s.first] == s.second);
        }
        Eigen::VectorXd pheno_matrix = prsice.phenotype_matrix();
        REQUIRE(static_cast<size_t>(pheno_matrix.rows())
                == expected_pheno.size());
        for (size_t i = 0; i < expected_pheno.size(); ++i)
        { REQUIRE(pheno_matrix(i, 0) == Approx(expected_pheno[i])); }
    }
    SECTION("Invalid cov file format")
    {
        std::unique_ptr<std::istream> cov_file =
            std::make_unique<std::istringstream>(cov_header
                                                 + "1 1 0.1 NA M 10 B1 c1\n"
                                                   "2 2 0.1 NA M 10\n");

        REQUIRE_THROWS_WITH(prsice.test_cov_check_and_factor_level_count(
                                factor_idx, cov_names, cov_idx, delim,
                                ignore_fid, cov_file, geno),
                            Catch::Contains("Error: Malformed covariate file, "
                                            "should have at least 8 columns"));
    }
    SECTION("Duplicated covariate ID")
    {
        std::string input = cov_header;
        for (size_t i = 0; i < num_sample; ++i)
        {
            input.append(
                std::to_string(i) + " " + std::to_string(i) + " "
                + std::to_string(pheno()) + " " + std::to_string(pheno()) + " "
                + std::to_string(pheno()) + " " + std::to_string(pheno()) + " "
                + std::to_string(pheno()) + " " + std::to_string(pheno())
                + "\n");
            input.append(
                std::to_string(i) + " " + std::to_string(i) + " "
                + std::to_string(pheno()) + " " + std::to_string(pheno()) + " "
                + std::to_string(pheno()) + " " + std::to_string(pheno()) + " "
                + std::to_string(pheno()) + " " + std::to_string(pheno())
                + "\n");
        }
        std::unique_ptr<std::istream> cov_file =
            std::make_unique<std::istringstream>(input);
        REQUIRE_THROWS_WITH(
            prsice.test_cov_check_and_factor_level_count(
                factor_idx, cov_names, cov_idx, delim, ignore_fid, cov_file,
                geno),
            Catch::Contains("duplicated IDs in covariate file"));
    }
    SECTION("No valid samples")
    {
        std::unique_ptr<std::istream> cov_file =
            std::make_unique<std::istringstream>(
                cov_header
                + "1 1 NA 0.1 0.1 M 10 b1 c1\n"
                  "2 2 0.1 0.1 NA M 10 b1 c1\n"
                  "3 3 0.1 0.1 0.1 NA 10 b1 c1\n"
                  "4 4 0.1 0.1 0.1 M NA b1 c1\n"
                  "4 4 0.1 0.1 0.1 M 10 NA c1\n");
        REQUIRE_THROWS_WITH(prsice.test_cov_check_and_factor_level_count(
                                factor_idx, cov_names, cov_idx, delim,
                                ignore_fid, cov_file, geno),
                            Catch::Contains("Error: All samples removed due to "
                                            "missingness in covariate file!"));
    }
}
TEST_CASE("Update valid samples from m_phenotype")
{
    // we need a genotype object with sample initialized (SampleID)
    // we also need to initialize the m_phenotype matrix
    // first generate samples
    const size_t num_sample = 1000;
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()};
    // select ~70% of samples
    std::uniform_int_distribution<size_t> dist {1, 10};
    auto valid = [&dist, &mersenne_engine]() {
        return dist(mersenne_engine) > 7;
    };
    std::normal_distribution pheno_dist;
    auto pheno = [&pheno_dist, &mersenne_engine]() {
        return pheno_dist(mersenne_engine);
    };
    mockGenotype geno;
    Reporter reporter("log", 60, true);
    geno.set_reporter(&reporter);
    std::vector<double> sample_pheno;
    std::vector<bool> valid_after_covariate(num_sample, false);
    size_t num_cov_valid = 0;
    std::vector<double> expected_pheno;
    const std::string delim = " ";
    auto ignore_fid = GENERATE(true, false);
    std::unordered_map<std::string, size_t> expected_id_map;
    for (size_t i = 0; i < num_sample; ++i)
    {
        auto cur_pheno = pheno();
        auto valid_sample = valid();
        geno.add_sample(Sample_ID(std::to_string(i), std::to_string(i),
                                  std::to_string(cur_pheno), valid_sample));
        geno.update_valid_sample(i, false);
        if (valid_sample)
        {
            sample_pheno.push_back(cur_pheno);
            geno.update_valid_sample(i, true);
            valid_after_covariate[i] = valid();
            if (valid_after_covariate[i])
            {
                auto id = std::to_string(i);
                if (!ignore_fid) id.append(delim + std::to_string(i));
                expected_id_map[id] = num_cov_valid;
                expected_pheno.push_back(cur_pheno);
                ++num_cov_valid;
            }
        }
    }

    mock_prsice prsice;
    prsice.phenotype_matrix() = Eigen::Map<Eigen::VectorXd>(
        sample_pheno.data(), static_cast<Eigen::Index>(sample_pheno.size()));
    prsice.test_update_phenotype_matrix(valid_after_covariate, delim,
                                        num_cov_valid, ignore_fid, geno);
    Eigen::VectorXd res = prsice.phenotype_matrix();
    REQUIRE(static_cast<size_t>(res.rows()) == expected_pheno.size());
    for (size_t i = 0; i < expected_pheno.size(); ++i)
    { REQUIRE(res(i, 0) == Approx(expected_pheno[i])); }
    auto sample_map = prsice.sample_with_phenotypes();
    REQUIRE(sample_map.size() == expected_id_map.size());
    for (auto&& id : expected_id_map)
    {
        auto loc = sample_map.find(id.first);
        REQUIRE(loc != sample_map.end());
        REQUIRE(loc->second == id.second);
    }
}

TEST_CASE("Get covariate start position")
{
    mock_prsice prsice;
    Reporter reporter("log", 60, true);
    prsice.set_reporter(&reporter);
    std::vector<std::unordered_map<std::string, size_t>> factor_levels = {
        {{"A", 0}, {"B", 1}, {"C", 2}, {"D", 3}, {"E", 4}},
        {{"M", 0}, {"F", 1}},
        {{"batch1", 0}, {"batch3", 1}, {"batch2", 2}}};
    std::set<size_t> is_factor = {3, 4, 9};
    std::vector<size_t> cov_idx = {1, 3, 4, 5, 6, 9};
    auto [cov_start, num_col] =
        prsice.test_get_cov_start(factor_levels, is_factor, cov_idx);
    REQUIRE(num_col == 12);
    std::vector<size_t> expected_start = {2, 3, 7, 8, 9, 10};
    REQUIRE_THAT(cov_start, Catch::Equals<size_t>(expected_start));
}

TEST_CASE("Generate covariate matrix")
{
    // need to get the number of column, number of covariate start the factor
    // levels and a cov file
    size_t num_sample = 1000;
    // one numeric covariate and one factor
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()};
    std::uniform_int_distribution<size_t> dist {1, 4};
    std::uniform_int_distribution<size_t> valid_cov {1, 10};
    std::normal_distribution pheno_dist;
    auto pheno = [&pheno_dist, &mersenne_engine] {
        return pheno_dist(mersenne_engine);
    };
    auto is_valid = [&valid_cov, &mersenne_engine] {
        return valid_cov(mersenne_engine) > 3;
    };
    auto factor = [&dist, &mersenne_engine] {
        return "b" + std::to_string(dist(mersenne_engine));
    };
    const std::string cov_header = "FID IID PC1 Batch\n";
    std::vector<std::string> cov_input;
    mock_prsice prsice;
    Reporter reporter("log", 60, true);
    prsice.set_reporter(&reporter);
    auto&& sample_with_pheno = prsice.sample_with_phenotypes();
    auto ignore_fid = GENERATE(true, false);
    const std::string delim = " ";
    std::unordered_set<std::string> valid_factor;
    size_t valid_idx = 0;
    mockGenotype geno;
    geno.set_reporter(&reporter);
    size_t n_fam_sample = 0;
    std::vector<double> pheno_value;
    for (size_t i = 0; i < num_sample; ++i)
    {
        auto cur_pheno = pheno();
        auto cur_factor = factor();

        auto id = std::to_string(i);
        cov_input.push_back(id + " " + id + " " + std::to_string(cur_pheno)
                            + " " + cur_factor + "\n");

        if (is_valid())
        {
            if (!ignore_fid) { id.append(delim + std::to_string(i)); }
            sample_with_pheno[id] = valid_idx++;
            valid_factor.insert(cur_factor);
            geno.add_sample(Sample_ID(std::to_string(i), std::to_string(i),
                                      std::to_string(cur_pheno), true));
            pheno_value.push_back(
                misc::convert<double>(std::to_string(cur_pheno)));
            ++n_fam_sample;
        }
    }
    prsice.phenotype_matrix() = Eigen::Map<Eigen::VectorXd>(
        pheno_value.data(), static_cast<Eigen::Index>(pheno_value.size()));
    geno.set_sample_vector(n_fam_sample);
    geno.set_sample(n_fam_sample);
    std::shuffle(cov_input.begin(), cov_input.end(), mersenne_engine);
    size_t num_factor = valid_factor.size();
    Eigen::MatrixXd expected_independent =
        Eigen::MatrixXd::Zero(sample_with_pheno.size(), 2 + num_factor);
    expected_independent.col(0).setOnes();
    expected_independent.col(1).setOnes();
    std::string cov_input_str = cov_header;
    std::vector<std::string> token;
    std::vector<std::unordered_map<std::string, size_t>> factor_level(1);
    auto cur_level = factor_level.front();
    size_t num_level = 0;
    std::set<size_t> is_factor = {3};
    std::vector<size_t> cov_idx = {2, 3}, factor_idx = {3};
    for (size_t i = 0; i < cov_input.size(); ++i)
    {
        cov_input_str.append(cov_input[i]);
        misc::trim(cov_input[i]);
        token = misc::split(cov_input[i]);
        auto id = token[0];
        if (!ignore_fid) id.append(delim + token[1]);
        auto found = sample_with_pheno.find(id);
        if (found == sample_with_pheno.end()) continue;
        const auto idx = found->second;
        expected_independent(idx, 2) = misc::convert<double>(token[2]);

        if (factor_level[0].find(token[3]) != factor_level[0].end())
        {
            if (factor_level[0][token[3]] != 0)
                expected_independent(idx, 2 + factor_level[0][token[3]]) = 1;
        }
        else
        {
            factor_level[0][token[3]] = num_level++;
            if (factor_level[0][token[3]] != 0)
                expected_independent(idx, 2 + factor_level[0][token[3]]) = 1;
        }
    }
    std::unique_ptr<std::istream> cov_file =
        std::make_unique<std::istringstream>(cov_input_str);
    SECTION("Propagate independent matrix")
    {
        std::vector<size_t> cov_start = {2, 3};
        prsice.init_independent(valid_idx, num_factor + 2);
        prsice.test_propagate_independent_matrix(
            factor_level, is_factor, cov_idx, cov_start, delim, ignore_fid,
            std::move(cov_file));
        Eigen::MatrixXd result = prsice.get_independent();
        REQUIRE(result == expected_independent);
    }
    SECTION("Read from covariate file")
    {
        SECTION("no file")
        {
            REQUIRE_NOTHROW(prsice.test_gen_cov_matrix(
                std::vector<std::string> {}, std::vector<size_t> {},
                std::vector<size_t> {}, "", delim, ignore_fid, geno));
            REQUIRE(prsice.get_independent()
                    == Eigen::MatrixXd::Ones(n_fam_sample, 2));
        }
        SECTION("From file")
        {
            std::ofstream cov_out("cov_file_test");
            cov_out << cov_input_str << std::endl;
            cov_out.close();
            REQUIRE_NOTHROW(prsice.test_gen_cov_matrix(
                std::vector<std::string> {"PC1", "Batch"}, cov_idx, factor_idx,
                "cov_file_test", delim, ignore_fid, geno));
            REQUIRE(prsice.get_independent() == expected_independent);
        }
    }
}


TEST_CASE("Full init_matrix test")
{
    SECTION("No regress") {}
}
