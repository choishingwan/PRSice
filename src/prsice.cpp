// This file is part of PRSice-2, copyright (C) 2016-2019
// Shing Wan Choi, Paul F. O’Reilly
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "prsice.hpp"

std::mutex PRSice::lock_guard;
std::tuple<size_t, bool>
PRSice::get_pheno_idx(const std::vector<std::string_view>& column,
                      const Phenotype& pheno_info, const std::string& pheno)
{
    if (pheno.empty()) return {1 + !pheno_info.ignore_fid, false};
    // there should not be nay duplicated column in the pheno_col
    bool duplicated_file_column = false;
    size_t col_idx = 0;
    // don't start from 0 as we expect those to be the FID / IID
    for (size_t i_col = 1 + !pheno_info.ignore_fid; i_col < column.size();
         ++i_col)
    {
        if (column[i_col] == pheno)
        {
            if (duplicated_file_column)
            {
                throw std::runtime_error(
                    "Error: Multiple column in the phenotype file mataches "
                    "with the require pehnotyp name: "
                    + pheno);
            }
            duplicated_file_column = true;
            col_idx = i_col;
        }
    }
    return {col_idx, duplicated_file_column};
}
void PRSice::parse_pheno_header(std::unique_ptr<std::istream> pheno_file,
                                Phenotype& pheno_info, Reporter& reporter)
{
    std::string line;
    std::getline(*pheno_file, line);
    if (line.empty())
    {
        throw std::runtime_error(
            "Error: First line of phenotype file cannot be empty");
    }
    pheno_file.reset();
    misc::trim(line);
    std::vector<std::string_view> col = misc::tokenize(line);
    if (col.size() < 2ul + !pheno_info.ignore_fid)
    {
        throw std::runtime_error(
            "Error: Not enough column in Phenotype file. "
            "If the phenotype does not contain the FID, use --ignore-fid");
    }
    std::string message = "";
    std::string sample_id = std::string(col[0]);
    if (!pheno_info.ignore_fid) { sample_id.append("+" + std::string(col[1])); }
    message.append("Phenotype file: " + pheno_info.pheno_file + "\n");
    message.append("Column Name of Sample ID: " + sample_id + "\n");
    message.append("Note: If the phenotype file does not contain a header, "
                   "the column name will be displayed as the Sample ID "
                   "which is expected.\n");
    // no user input
    if (pheno_info.pheno_col.empty())
    {
        pheno_info.pheno_col_idx.push_back(1 + !pheno_info.ignore_fid);
        pheno_info.pheno_col = {std::string(col[1 + !pheno_info.ignore_fid])};
        if (isdigit(col[1 + !pheno_info.ignore_fid].at(0)))
        { pheno_info.pheno_col = {"Phenotype"}; }
        pheno_info.skip_pheno = {false};
    }
    else
    {
        pheno_info.skip_pheno.resize(pheno_info.pheno_col.size());
        bool has_valid_pheno = false;
        for (size_t i_pheno = 0; i_pheno < pheno_info.pheno_col.size();
             ++i_pheno)
        {
            auto [idx, found] =
                get_pheno_idx(col, pheno_info, pheno_info.pheno_col[i_pheno]);
            if (found) { has_valid_pheno = true; }
            else
            {
                message.append("Warning: Phenotype: "
                               + pheno_info.pheno_col[i_pheno]
                               + " cannot be found in the phenotype file\n");
                pheno_info.skip_pheno[i_pheno] = true;
            }
            pheno_info.pheno_col_idx.push_back(idx);
        }
        if (!has_valid_pheno)
        {
            message.append("Error: None of the phenotype(s) can be found "
                           "in the phenotype file!\n");
            throw std::runtime_error(message);
        }
    }
    reporter.report(message);
}

void PRSice::pheno_check(const bool no_regress, Phenotype& pheno,
                         Reporter& reporter)
{
    // don't bother to check anything, just add a place holder in binary
    if (no_regress)
    {
        pheno.binary = {true};
        pheno.pheno_col = {"PlaceHolder"};
        pheno.skip_pheno = {false};
        pheno.pheno_col_idx = {~size_t(0)};
        return;
    }
    if (pheno.binary.empty())
    { throw std::runtime_error("Error: No phenotype provided"); }

    // want to update the binary and prevalence vector by removing
    // phenotypes not found / duplicated
    pheno.skip_pheno.resize(pheno.binary.size(), false);
    if (!pheno.pheno_file.empty())
    {
        auto pheno_file = misc::load_stream(pheno.pheno_file);
        parse_pheno_header(std::move(pheno_file), pheno, reporter);
    }
    else
    {
        pheno.pheno_col = {"Phenotype"};
        pheno.pheno_col_idx = {~size_t(0)};
    }
    assert(m_pheno_info.pheno_col_idx.size() == m_pheno_info.pheno_col.size());
    auto message = "There are a total of "
                   + std::to_string(pheno.pheno_col.size())
                   + " phenotype to process\n";
    reporter.report(message);
}

void PRSice::init_matrix(const Phenotype& pheno_info, const std::string& delim,
                         const size_t pheno_idx, Genotype& target)
{
    const auto file_idx = pheno_info.pheno_col_idx[pheno_idx];
    const auto file_name = pheno_info.pheno_file;
    const auto pheno_name = pheno_info.pheno_col[pheno_idx];
    const auto cov_names = pheno_info.cov_colname;
    const auto cov_idx = pheno_info.col_index_of_cov;
    const auto factor_idx = pheno_info.col_index_of_factor_cov;
    const auto cov_file_name = pheno_info.cov_file;
    const auto ignore_fid = pheno_info.ignore_fid;
    const auto no_regress = m_prs_info.no_regress;
    if (m_binary_trait && m_prs_info.scoring_method == SCORING::CONTROL_STD)
    { target.reset_std_flag(); }
    // we need genotype for no-regress if we are trying to do control std
    gen_pheno_vec(file_name, pheno_name, delim, file_idx, ignore_fid, target);
    if (!no_regress)
    {
        // won't use covariate when no regression is performed
        gen_cov_matrix(cov_names, cov_idx, factor_idx, cov_file_name, delim,
                       ignore_fid, target);
    }
    if (m_binary_trait && m_prs_info.scoring_method == SCORING::CONTROL_STD)
        set_std_exclusion_flag(delim, ignore_fid, target);
    m_matrix_index = get_matrix_idx(delim, ignore_fid, target);
    if (no_regress) return;
    double null_r2_adjust = 0.0;
    bool has_covariate = m_independent_variables.cols() > 2;
    if (has_covariate)
    {
        auto n_thread = m_prs_info.thread;
        // only do it if we have the correct number of sample
        assert(m_independent_variables.rows() == m_phenotype.rows());
        if (m_binary_trait)
        {
            // ignore the first column
            // this is ok as both the first column (intercept) and the
            // second column (PRS) is currently 1
            Regression::glm(m_phenotype,
                            m_independent_variables.topRightCorner(
                                m_independent_variables.rows(),
                                m_independent_variables.cols() - 1),
                            m_null_p, m_null_r2, m_null_coeff, m_null_se,
                            n_thread);
        }
        else
        {
            // ignore the first column
            // and perform linear regression
            Regression::fastLm(m_phenotype,
                               m_independent_variables.topRightCorner(
                                   m_independent_variables.rows(),
                                   m_independent_variables.cols() - 1),
                               m_null_p, m_null_r2, null_r2_adjust,
                               m_null_coeff, m_null_se, n_thread, true);
        }
    }
    m_best_sample_score.resize(target.num_sample());
}

void PRSice::parse_pheno(const std::string& pheno,
                         std::vector<double>& pheno_store, int& max_pheno_code)
{
    if (m_binary_trait)
    {
        int temp = misc::convert<int>(pheno);
        if (temp >= 0 && temp <= 2)
        {
            pheno_store.push_back(temp);
            if (max_pheno_code < temp) max_pheno_code = temp;
        }
        else
        {
            throw std::runtime_error("Invalid binary phenotype format!");
        }
    }
    else
    {
        pheno_store.push_back(misc::convert<double>(pheno));
    }
}

std::unordered_map<std::string, std::string>
PRSice::load_pheno_map(const std::string& delim, const size_t idx,
                       const bool ignore_fid,
                       std::unique_ptr<std::istream> pheno_file)
{
    std::unordered_map<std::string, std::string> phenotype_info;
    std::vector<std::string_view> token;
    std::string line, id;
    while (std::getline(*pheno_file, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::tokenize(line);
        // Check if we have the minimal required column number
        if (token.size() < idx + 1)
        {
            throw std::runtime_error(
                "Malformed pheno file, should contain at least "
                + misc::to_string(idx + 1)
                + " columns. "
                  "Have you use the --ignore-fid option?");
        }
        id = (ignore_fid)
                 ? std::string(token[0])
                 : std::string(token[0]) + delim + std::string(token[1]);
        if (phenotype_info.find(id) != phenotype_info.end())
        {
            throw std::runtime_error("Error: Duplicated sample ID in "
                                     "phenotype file: "
                                     + id
                                     + ". Please "
                                       "check if your input is correct!");
        }
        phenotype_info[id] = std::string(token[idx]);
    }
    pheno_file.reset();
    return phenotype_info;
}

std::vector<size_t> PRSice::get_matrix_idx(const std::string& delim,
                                           const bool ignore_fid,
                                           Genotype& target)
{
    std::vector<size_t> matrix_idx;
    for (size_t i_sample = 0; i_sample < target.num_sample(); ++i_sample)
    {
        auto id = ignore_fid ? target.iid(i_sample)
                             : target.sample_id(i_sample, delim);
        auto&& pheno_idx = m_sample_with_phenotypes.find(id);
        if (pheno_idx == m_sample_with_phenotypes.end()) { continue; }
        matrix_idx.push_back(i_sample);
    }
    return matrix_idx;
}
void PRSice::set_std_exclusion_flag(const std::string& delim,
                                    const bool ignore_fid, Genotype& target)
{
    for (size_t i_sample = 0; i_sample < target.num_sample(); ++i_sample)
    {
        auto id = ignore_fid ? target.iid(i_sample)
                             : target.sample_id(i_sample, delim);
        auto&& pheno_idx = m_sample_with_phenotypes.find(id);
        if (pheno_idx == m_sample_with_phenotypes.end()) { continue; }
        if (!misc::logically_equal(m_phenotype(pheno_idx->second), 0))
        { target.exclude_from_std(i_sample); }
    }
}
std::tuple<bool, size_t, size_t>
PRSice::binary_pheno_is_valid(const int max_pheno_code,
                              std::vector<double>& pheno_store)
{
    size_t min_pheno = max_pheno_code == 1 ? 0 : 1;
    size_t num_case = 0;
    size_t num_control = 0;
    bool valid = true;
    for (auto&& pheno : pheno_store)
    {
        if (pheno < min_pheno || pheno > max_pheno_code)
        {
            valid = false;
            break;
        }
        pheno -= min_pheno;
        misc::logically_equal(pheno, 1) ? ++num_case : ++num_control;
    }
    return {valid, num_case, num_control};
}
bool PRSice::quantitative_pheno_is_valid(const std::vector<double>& pheno_store)
{
    // assume pheno_store is valid
    double prev = pheno_store.front();
    for (auto&& p : pheno_store)
    {
        if (!misc::logically_equal(p, prev)) { return true; }
    }
    return false;
}
std::tuple<std::vector<double>, size_t, int>
PRSice::process_phenotype_info(const std::string& delim, const bool ignore_fid,
                               Genotype& target)
{
    const size_t sample_ct = target.num_sample();
    std::vector<double> pheno_store;
    pheno_store.reserve(sample_ct);
    size_t invalid_pheno = 0;
    int max_pheno_code = 0;
    size_t pheno_matrix_idx = 0;
    for (size_t i_sample = 0; i_sample < sample_ct; ++i_sample)
    {
        target.update_valid_sample(i_sample, false);
        if (target.pheno_is_na(i_sample)
            || !target.sample_selected_for_prs(i_sample)
            || (m_binary_trait && target.pheno(i_sample) == "-9"))
        { continue; }
        try
        {
            parse_pheno(target.pheno(i_sample), pheno_store, max_pheno_code);
            auto id = ignore_fid ? target.iid(i_sample)
                                 : target.sample_id(i_sample, delim);
            m_sample_with_phenotypes[id] = pheno_matrix_idx++;
            target.update_valid_sample(i_sample, true);
        }
        catch (const std::runtime_error&)
        {
            ++invalid_pheno;
        }
    }
    return {pheno_store, invalid_pheno, max_pheno_code};
}


std::tuple<std::vector<double>, size_t, size_t, int>
PRSice::process_phenotype_file(const std::string& file_name,
                               const std::string& delim,
                               const std::size_t pheno_idx,
                               const bool ignore_fid, Genotype& target)
{
    auto pheno_stream = misc::load_stream(file_name);
    auto phenotype_info =
        load_pheno_map(delim, pheno_idx, ignore_fid, std::move(pheno_stream));
    const size_t sample_ct = target.num_sample();
    std::string id;
    int max_pheno_code = 0;
    size_t invalid_pheno = 0;
    size_t num_not_found = 0;
    std::vector<double> pheno_store;
    pheno_store.reserve(sample_ct);
    size_t pheno_matrix_idx = 0;
    for (size_t i_sample = 0; i_sample < sample_ct; ++i_sample)
    {
        target.update_valid_sample(i_sample, false);
        id = (ignore_fid) ? target.iid(i_sample)
                          : target.sample_id(i_sample, delim);
        auto pheno_in_file = phenotype_info.find(id);
        if (pheno_in_file != phenotype_info.end())
        {
            auto pheno_tmp = pheno_in_file->second;
            misc::to_lower(pheno_tmp);
            if (pheno_tmp != "na" && phenotype_info[id] != "nan"
                && target.sample_selected_for_prs(i_sample)
                && !(m_binary_trait && target.pheno(i_sample) == "-9"))
            {
                try
                {
                    parse_pheno(phenotype_info[id], pheno_store,
                                max_pheno_code);
                    m_sample_with_phenotypes[id] = pheno_matrix_idx++;
                    target.update_valid_sample(i_sample, true);
                }
                catch (...)
                {
                    ++invalid_pheno;
                }
            }
        }
        else
        {
            ++num_not_found;
        }
    }
    return {pheno_store, num_not_found, invalid_pheno, max_pheno_code};
}
void PRSice::print_pheno_log(const std::string& name, const size_t sample_ct,
                             const size_t num_not_found,
                             const size_t invalid_pheno,
                             const int max_pheno_code, const bool ignore_fid,
                             std::vector<double>& pheno_store)
{
    std::string message = name + " is a ";
    if (m_binary_trait) { message.append("binary phenotype\n"); }
    else
    {
        message.append("continuous phenotype\n");
    }
    if (num_not_found != 0)
    {
        message.append(std::to_string(num_not_found)
                       + " sample(s) without phenotype\n");
    }
    if (invalid_pheno != 0)
    {
        message.append(std::to_string(invalid_pheno)
                       + " sample(s) with invalid phenotype\n");
    }
    if (num_not_found == sample_ct)
    {
        message.append("None of the target samples were found in the "
                       "phenotype file. ");
        if (ignore_fid)
        {
            message.append("Maybe the first column of your phenotype file "
                           "is the FID?");
        }
        else
        {
            message.append(
                "Maybe your phenotype file doesn not contain the FID?\n"
                "Might want to consider using --ignore-fid\n");
        }
        message.append("Or it is possible that only non-founder sample have "
                       "phenotype information "
                       " and you did not use "
                       "--nonfounders?\n");
        throw std::runtime_error(message + "Error: No sample left");
    }
    if (invalid_pheno == sample_ct)
    {
        throw std::runtime_error("Error: All sample has invalid phenotypes!\n"
                                 "Error: No sample left");
    }
    if (pheno_store.empty())
    { throw std::runtime_error("No phenotype presented"); }
    if (m_binary_trait)
    {
        auto [valid, num_case, num_control] =
            binary_pheno_is_valid(max_pheno_code, pheno_store);
        if (!valid)
        {
            throw std::runtime_error(
                "Mixed encoding! Both 0/1 and 1/2 encoding found!");
        }
        message.append(std::to_string(num_control) + " control(s)\n");
        message.append(std::to_string(num_case) + " case(s)\n");
        if (num_control == 0)
        {
            m_reporter->report(message);
            throw std::runtime_error("There are no control samples");
        }
        if (num_case == 0)
        {
            m_reporter->report(message);
            throw std::runtime_error("There are no cases");
        }
    }
    else
    {
        bool valid = quantitative_pheno_is_valid(pheno_store);
        if (!valid)
        {
            m_reporter->report(message);
            std::string error_message = "Only one phenotype value detected";
            if (misc::logically_equal(pheno_store.front(), -9))
            { error_message.append(" and they are all -9"); }
            throw std::runtime_error(error_message
                                     + ". Not enough valid phenotype");
        }
        else
        {
            message.append(std::to_string(m_phenotype.rows())
                           + " sample(s) with valid phenotype\n");
        }
    }
    m_reporter->report(message);
}

void PRSice::gen_pheno_vec(const std::string& pheno_file,
                           const std::string& pheno_name,
                           const std::string& delim,
                           const size_t pheno_file_idx, const bool ignore_fid,
                           Genotype& target)
{
    // we will first store the phenotype into the double vector and then
    // later use this to construct the matrix
    const size_t sample_ct = target.num_sample();
    size_t num_not_found = 0;
    size_t invalid_pheno = 0;
    int max_pheno_code = 0;
    std::vector<double> pheno_store;
    // m_sample_with_phenotype will have ID as key and idx on the vector
    // SampleID as value because we know the phenobtype are stored w.r.t the
    // order of vector SampleID
    if (!pheno_file.empty()) // use phenotype file
    {
        std::tie(pheno_store, num_not_found, invalid_pheno, max_pheno_code) =
            process_phenotype_file(pheno_file, delim, pheno_file_idx,
                                   ignore_fid, target);
    }
    else
    {
        // No phenotype file is provided
        // Use information from the fam file directly
        std::tie(pheno_store, invalid_pheno, max_pheno_code) =
            process_phenotype_info(delim, ignore_fid, target);
    }
    print_pheno_log(pheno_name, sample_ct, num_not_found, invalid_pheno,
                    max_pheno_code, ignore_fid, pheno_store);
    m_phenotype = Eigen::Map<Eigen::VectorXd>(
        pheno_store.data(), static_cast<Eigen::Index>(pheno_store.size()));
}
bool PRSice::is_valid_covariate(const std::set<size_t>& factor_idx,
                                const std::vector<size_t>& cov_idx,
                                std::vector<std::string>& cov_line,
                                std::vector<size_t>& missing_count)
{
    bool valid = true;
    for (size_t idx = 0; idx < cov_idx.size(); ++idx)
    {
        auto cur_cov_idx = cov_idx[idx];
        auto cur_cov = cov_line[cur_cov_idx];
        misc::to_upper(cur_cov);
        if (cur_cov == "NAN" || cur_cov == "NA")
        {
            ++missing_count[idx];
            valid = false;
        }
        else if (factor_idx.find(cur_cov_idx) == factor_idx.end())
        {
            // not factor
            try
            {
                misc::convert<double>(cur_cov);
            }
            catch (...)
            {
                ++missing_count[idx];
                valid = false;
            }
        }
    }
    return valid;
}
std::string PRSice::output_missing(const std::set<size_t>& factor_idx,
                                   const std::vector<std::string>& cov_names,
                                   const std::vector<size_t>& cov_idx,
                                   const std::vector<size_t>& factor_levels,
                                   const std::vector<size_t>& missing_count)
{
    std::string message =
        "Include Covariates:\nName\tMissing\tNumber of levels\n";
    size_t cur_name_idx = 0;
    size_t i_factor = 0;
    for (auto&& cov : cov_idx)
    {
        if (factor_idx.find(cov) == factor_idx.end())
        {
            // non-factor
            message.append(cov_names[cur_name_idx] + "\t"
                           + std::to_string(missing_count[cur_name_idx])
                           + "\t-\n");
        }
        else
        {
            message.append(cov_names[cur_name_idx] + "\t"
                           + std::to_string(missing_count[cur_name_idx]) + "\t"
                           + std::to_string(factor_levels[i_factor++]) + "\n");
        }
        ++cur_name_idx;
    }
    return message;
}
void PRSice::update_phenotype_matrix(const std::vector<bool>& valid_samples,
                                     const std::string& delim,
                                     const size_t num_valid,
                                     const bool ignore_fid, Genotype& target)
{
    Eigen::VectorXd new_pheno = Eigen::VectorXd::Zero(num_valid);
    const size_t num_sample = target.num_sample();
    size_t new_matrix_idx = 0, pheno_matrix_idx = 0;
    for (size_t i = 0; i < num_sample; ++i)
    {
        if (target.sample_valid_for_regress(i))
        {
            auto id = ignore_fid ? target.iid(i) : target.sample_id(i, delim);
            if (valid_samples[i])
            {
                m_sample_with_phenotypes[id] = new_matrix_idx;
                new_pheno(new_matrix_idx, 0) = m_phenotype(pheno_matrix_idx, 0);
                ++new_matrix_idx;
            }
            else
            {
                m_sample_with_phenotypes.erase(id);
                target.update_valid_sample(i, false);
            }
            ++pheno_matrix_idx;
        }
    }
    m_phenotype = new_pheno;
}

std::vector<std::unordered_map<std::string, size_t>>
PRSice::cov_check_and_factor_level_count(
    const std::set<size_t>& factor_idx,
    const std::vector<std::string>& cov_names,
    const std::vector<size_t>& cov_idx, const std::string& delim,
    const bool ignore_fid, std::unique_ptr<std::istream>& cov_file,
    Genotype& target)
{
    const size_t max_idx = cov_idx.back() + 1;
    std::vector<size_t> missing_count(cov_idx.size(), 0);
    std::vector<size_t> current_factor_level(factor_idx.size(), 0);
    std::unordered_set<std::string> duplicated_id;
    std::string line, id;
    std::vector<std::string> token;
    size_t num_duplicated_id = 0;
    // indicate if the sample is valid after covariate read
    // we need this as the covariate file might not follow order of the
    // fam file. Size should be total number of samples in vector not in the
    // matrix as we don't know the order
    std::vector<bool> valid_samples(target.num_sample(), false);
    bool valid;
    size_t num_valid = 0;
    std::vector<std::unordered_map<std::string, size_t>> factor_levels(
        factor_idx.size());
    while (std::getline(*cov_file, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() < max_idx)
        {
            throw std::runtime_error(
                "Error: Malformed covariate file, should have at least "
                + std::to_string(max_idx) + " columns");
        }
        id = (ignore_fid) ? token[0] : token[0] + delim + token[1];
        auto sample_idx = m_sample_with_phenotypes.find(id);
        if (sample_idx != m_sample_with_phenotypes.end())
        {
            valid =
                is_valid_covariate(factor_idx, cov_idx, token, missing_count);
            if (!valid) continue;
            if (duplicated_id.find(id) != duplicated_id.end())
            {
                ++num_duplicated_id;
                continue;
            }
            duplicated_id.insert(id);
            // sample is valid
            valid_samples[sample_idx->second] = true;
            ++num_valid;
            size_t i_factor = 0;
            for (auto&& f : factor_idx)
            {
                auto&& cur_level = factor_levels[i_factor];
                auto&& cur_cov = token[f];
                if (cur_level.find(cur_cov) == cur_level.end())
                { cur_level[cur_cov] = current_factor_level[i_factor]++; }
                ++i_factor;
            }
        }
    }
    if (num_duplicated_id != 0)
    {
        throw std::runtime_error("Error: " + std::to_string(num_duplicated_id)
                                 + " duplicated IDs in covariate file!\n");
    }
    else if (num_valid == 0)
    {
        throw std::runtime_error("Error: All samples removed due to "
                                 "missingness in covariate file!");
    }
    auto message = output_missing(factor_idx, cov_names, cov_idx,
                                  current_factor_level, missing_count);
    double ratio = static_cast<double>(num_valid)
                   / static_cast<double>(m_sample_with_phenotypes.size());
    if (ratio < 0.95)
    {
        message.append(
            "Warning: More than " + std::to_string((1.0 - ratio) * 100)
            + "% of your samples were removed! "
              "You should check if your covariate file is correct\n");
    }
    m_reporter->report(message);
    update_phenotype_matrix(valid_samples, delim, num_valid, ignore_fid,
                            target);
    // update phenotype matrix here
    cov_file->clear();
    cov_file->seekg(0, cov_file->beg);
    return factor_levels;
}
std::tuple<std::vector<size_t>, size_t> PRSice::get_cov_start(
    const std::vector<std::unordered_map<std::string, size_t>>& factor_levels,
    const std::set<size_t>& is_factor, const std::vector<size_t>& cov_idx)
{
    // start at 2 because 0 = intercept, 1 = PRS
    std::vector<size_t> cov_start(cov_idx.size(), 2);
    size_t num_matrix_col = 2;
    size_t i_factor = 0;
    for (size_t i = 0; i < cov_idx.size(); ++i)
    {
        auto cur_cov_idx = cov_idx[i];
        cov_start[i] = num_matrix_col;
        if (is_factor.find(cur_cov_idx) != is_factor.end())
        {
            // this is a factor
            num_matrix_col += factor_levels[i_factor].size() - 1;
            ++i_factor;
        }
        else
        {
            ++num_matrix_col;
        }
    }
    return {cov_start, num_matrix_col};
}

void PRSice::propagate_independent_matrix(
    const std::vector<std::unordered_map<std::string, size_t>>& factor_levels,
    const std::set<size_t>& is_factor, const std::vector<size_t>& cov_idx,
    const std::vector<size_t>& cov_start, const std::string& delim,
    const bool ignore_fid, std::unique_ptr<std::istream> cov_file)
{
    std::string line, id;
    // m_sample_with_phenotypes will tell us which row should we add the
    // covariate to
    std::vector<std::string> token;
    size_t i_factor = 0;
    while (std::getline(*cov_file, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        // don't need to check size, as we have done it when we check for valid
        // samples
        token = misc::split(line);
        id = ignore_fid ? token[0] : token[0] + delim + token[1];
        auto cur_sample = m_sample_with_phenotypes.find(id);
        if (cur_sample != m_sample_with_phenotypes.end())
        {
            auto row_idx = cur_sample->second;
            i_factor = 0;
            for (size_t i_col = 0; i_col < cov_idx.size(); ++i_col)
            {
                auto cur_cov_idx = cov_idx[i_col];
                if (is_factor.find(cur_cov_idx) != is_factor.end())
                {
                    // -1 so that the first non-reference level will start at
                    // the first column
                    auto cur_factor_level =
                        factor_levels[i_factor].at(token[cur_cov_idx]);
                    if (cur_factor_level != 0)
                    {
                        --cur_factor_level;
                        m_independent_variables(
                            row_idx, cov_start[i_col] + cur_factor_level) = 1;
                    }
                    ++i_factor;
                }
                else
                {
                    m_independent_variables(row_idx, cov_start[i_col]) =
                        misc::convert<double>(token[cur_cov_idx]);
                }
            }
        }
    }
    cov_file.reset();
}
void PRSice::gen_cov_matrix(const std::vector<std::string>& cov_names,
                            const std::vector<size_t>& cov_idx,
                            const std::vector<size_t>& factor_idx,
                            const std::string& cov_file_name,
                            const std::string& delim, const bool ignore_fid,
                            Genotype& target)
{
    Eigen::Index num_sample =
        static_cast<Eigen::Index>(m_sample_with_phenotypes.size());
    if (cov_file_name.empty())
    {
        m_independent_variables = Eigen::MatrixXd::Ones(num_sample, 2);
        return;
    }
    // convert factor idx to unordered_set for more elegant way of determining
    // if idx is factor
    std::set<size_t> is_factor;
    for (auto&& f : factor_idx) { is_factor.insert(f); }
    auto cov_file = misc::load_stream(cov_file_name);
    m_reporter->report("Processing the covariate file: " + cov_file_name
                       + "\n==============================\n");
    auto factor_levels = cov_check_and_factor_level_count(
        is_factor, cov_names, cov_idx, delim, ignore_fid, cov_file, target);
    auto [cov_start, num_matrix_col] =
        get_cov_start(factor_levels, is_factor, cov_idx);
    // update size of independent variable
    m_independent_variables =
        Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(num_sample),
                              static_cast<Eigen::Index>(num_matrix_col));
    m_independent_variables.col(0).setOnes();
    m_independent_variables.col(1).setOnes();
    // now propagate the matrix
    propagate_independent_matrix(factor_levels, is_factor, cov_idx, cov_start,
                                 delim, ignore_fid, std::move(cov_file));
    m_reporter->report("After reading the covariate file, "
                       + std::to_string(m_sample_with_phenotypes.size())
                       + " sample(s) included in the analysis\n");
}

void PRSice::reset_result_containers(const Genotype& target,
                                     const size_t region_idx)
{
    m_best_index = -1;
    m_num_snp_included = 0;
    // m_perm_result stores the result (Z scores) from each permutation and
    // is then used for calculation of empirical p value
    m_perm_result.resize(m_perm_info.num_permutation, 0);
    m_prs_results.resize(target.num_threshold(region_idx), prsice_result());
    std::fill(m_best_sample_score.begin(), m_best_sample_score.end(), 0);
}

void PRSice::print_all_score(const size_t num_sample,
                             std::unique_ptr<std::ostream>& all_score_file,
                             Genotype& target)
{
    for (size_t i_sample = 0; i_sample < num_sample; ++i_sample)
    {
        // we will calculate the the number of white space we need
        // to skip to reach the current sample + threshold's output
        // position
        const long long loc =
            m_all_file.header_length
            + static_cast<long long>(i_sample)
                  * (m_all_file.line_width + NEXT_LENGTH)
            + NEXT_LENGTH + m_all_file.skip_column_length
            + m_all_file.processed_threshold
            + m_all_file.processed_threshold * m_numeric_width;
        all_score_file->seekp(loc);
        // then we will output the score
        (*all_score_file) << std::setprecision(static_cast<int>(m_precision))
                          << target.calculate_score(i_sample);
    }
    ++m_all_file.processed_threshold;
}
void PRSice::run_prsice(const std::vector<size_t>& set_snp_idx,
                        const std::vector<std::string>& region_names,
                        const std::string& pheno_name, const double prevalence,
                        const size_t pheno_idx, const size_t region_idx,
                        const bool all_scores, const bool has_prevalence,
                        std::unique_ptr<std::ostream>& prsice_out,
                        std::unique_ptr<std::ostream>& best_score_file,
                        std::unique_ptr<std::ostream>& all_score_file,
                        Genotype& target)
{
    const bool print_all_scores = all_scores && pheno_idx == 0;
    const bool no_regress = m_prs_info.no_regress;
    const auto num_thread = m_prs_info.thread;
    const size_t num_sample = target.num_sample();
    if (set_snp_idx.empty()) return;
    Eigen::initParallel();
    Eigen::setNbThreads(m_prs_info.thread);
    reset_result_containers(target, region_idx);
    size_t prs_result_idx = 0;
    double cur_threshold = 0.0;
    print_progress();
    bool first_run = true;
    std::vector<size_t>::const_iterator start = set_snp_idx.begin();
    double top, bot;
    if (prevalence <= 1.0)
    { std::tie(top, bot) = lee_adjustment_factor(prevalence); }
    while (target.get_score(start, set_snp_idx.cend(), cur_threshold,
                            m_num_snp_included, first_run))
    {
        ++m_analysis_done;
        print_progress();
        if (print_all_scores && pheno_idx == 0)
        { print_all_score(num_sample, all_score_file, target); }
        if (!no_regress)
        {
            regress_score(target, cur_threshold, num_thread, prs_result_idx);
            print_prsice_output(m_prs_results[prs_result_idx], pheno_name,
                                region_names[region_idx], cur_threshold, top,
                                bot, has_prevalence, prsice_out);
            if (m_perm_info.run_perm) { permutation(num_thread); }
        }
        else
        {
            (*prsice_out) << pheno_name << "\t" << region_names[region_idx]
                          << "\t" << cur_threshold << "\t" << m_num_snp_included
                          << "\n";
        }
        ++prs_result_idx;
        first_run = false;
    }

    if (m_quick_best && !no_regress)
    {
        // if we can, store all best score in a matrix and output once to speed
        // things up
        m_fast_best_output.col(static_cast<Eigen::Index>(region_idx)) =
            Eigen::Map<Eigen::VectorXd>(
                m_best_sample_score.data(),
                static_cast<Eigen::Index>(m_best_sample_score.size()));
    }
    else if (!no_regress)
    {
        slow_print_best(best_score_file, target);
    }
    // we need to process the permutation result if permutation is required
    if (m_perm_info.run_perm) process_permutations();
    if (!no_regress && !(m_best_index < 0))
    {
        // postpone summary output until we have finished competitive
        // permutation
        auto&& best_info = m_prs_results[static_cast<size_t>(m_best_index)];
        m_prs_summary.push_back(prsice_summary(
            best_info, region_names[region_idx], (region_idx == 0)));
        if (best_info.p > 0.1)
            ++m_significant_store[0];
        else if (best_info.p > 1e-5)
            ++m_significant_store[1];
        else
            ++m_significant_store[2];
    }
}


void PRSice::slow_print_best(std::unique_ptr<std::ostream>& best_file,
                             Genotype& target)
{
    if (m_quick_best) return;
    if (m_best_index < 0)
    {
        // no best threshold
        m_reporter->report("Error: No best score obtained\nCannot output the "
                           "best PRS score\n");
        return;
    }
    auto&& best_info = m_prs_results[static_cast<size_t>(m_best_index)];
    size_t best_snp_size = best_info.num_snp;
    if (best_snp_size == 0)
    {
        m_reporter->report("Error: Best R2 obtained when no SNPs were "
                           "included\n"
                           "Cannot output the best PRS score\n");
    }
    else
    {
        for (size_t i_sample = 0; i_sample < target.num_sample(); ++i_sample)
        {
            long long loc = m_best_file.header_length
                            + static_cast<long long>(i_sample)
                                  * (m_best_file.line_width + NEXT_LENGTH)
                            + NEXT_LENGTH + m_best_file.skip_column_length
                            + m_best_file.processed_threshold
                            + m_best_file.processed_threshold * m_numeric_width;
            best_file->seekp(loc);
            (*best_file) << std::setprecision(static_cast<int>(m_precision))
                         << m_best_sample_score[i_sample];
        }
    }
    ++m_best_file.processed_threshold;
}

void PRSice::print_best(
    const std::vector<std::vector<std::size_t>>& region_membership,
    std::unique_ptr<std::ostream> best_file, Genotype& target)
{
    if (!m_quick_best)
    {
        // finished printing
        best_file.reset();
        return;
    }
    if (m_best_index < 0)
    {
        // no best threshold
        m_reporter->report("Error: No best score obtained\nCannot output the "
                           "best PRS score\n");
        return;
    }

    auto&& best_info = m_prs_results[static_cast<size_t>(m_best_index)];
    size_t best_snp_size = best_info.num_snp;
    if (best_snp_size == 0)
    {
        m_reporter->report("Error: Best R2 obtained when no SNPs were "
                           "included\nCannot output the best PRS score\n");
    }
    else
    {
        for (size_t i_sample = 0; i_sample < target.num_sample(); ++i_sample)
        {
            (*best_file) << target.sample_id(i_sample, " ") << " "
                         << ((target.sample_valid_for_regress(i_sample)) ? "Yes"
                                                                         : "No")
                         << std::setprecision(static_cast<int>(m_precision));
            for (Eigen::Index i = 0; i < m_fast_best_output.cols(); ++i)
            {
                if (i == 1 || region_membership[i].empty()) continue;
                (*best_file) << " " << m_fast_best_output(i_sample, i);
            }
            (*best_file) << "\n";
        }
        // can just close it as we assume we only need to do it once.
        best_file.reset();
    }
}

void PRSice::regress_score(Genotype& target, const double threshold,
                           const int thread, const size_t prs_result_idx)
{
    double r2 = 0.0, r2_adjust = 0.0, p_value = 0.0, coefficient = 0.0,
           se = 0.0;
    const auto num_regress_samples = m_matrix_index.size();
    // should never have m_num_snp_included == 0
    assert(!m_prs_results.empty());
    if (m_num_snp_included == m_prs_results[prs_result_idx].num_snp
        && !m_prs_info.non_cumulate)
    { return; }

    for (size_t sample_id = 0; sample_id < num_regress_samples; ++sample_id)
    {
        m_independent_variables(sample_id, 1) =
            target.calculate_score(m_matrix_index[sample_id]);
    }

    if (m_binary_trait)
    {
        try
        {
            Regression::glm(m_phenotype, m_independent_variables, p_value, r2,
                            coefficient, se, thread);
        }
        catch (const std::runtime_error& error)
        {
            fprintf(stderr, "Error: GLM model did not converge!\n");
            fprintf(stderr,
                    "       This is usually caused by small sample\n"
                    "       size or caused by problem in the input file\n");
            fprintf(stderr, "Error: %s\n", error.what());
        }
    }
    else
    {
        // we can run the linear regression
        Regression::fastLm(m_phenotype, m_independent_variables, p_value, r2,
                           r2_adjust, coefficient, se, thread, true);
    }
    // If this is the best r2, then we will add it
    int best_index = m_best_index;
    if (prs_result_idx == 0 || best_index < 0
        || m_prs_results[static_cast<size_t>(best_index)].r2 < r2)
    {
        m_best_index = static_cast<int>(prs_result_idx);
        const size_t num_include_samples = target.num_sample();
        // load all sample, including those that are not used for regression
        for (size_t s = 0; s < num_include_samples; ++s)
        { m_best_sample_score[s] = target.calculate_score(s); }
    }
    // we can now store the prsice_result

    m_prs_results[prs_result_idx] =
        prsice_result(threshold, r2, r2_adjust, coefficient, p_value, -1, se,
                      -1, m_num_snp_included);
}


void PRSice::process_permutations()
{
    // can't generate an empirical p-value if there is no observed p-value
    if (m_best_index == -1) return;
    size_t best_index = static_cast<size_t>(m_best_index);
    const double best_t = std::fabs(m_prs_results[best_index].coefficient
                                    / m_prs_results[best_index].se);
    const auto num_better =
        std::count_if(m_perm_result.begin(), m_perm_result.end(),
                      [&best_t](double t) { return t > best_t; });
    m_prs_results[best_index].emp_p =
        (num_better + 1.0) / (m_perm_info.num_permutation + 1.0);
}

void PRSice::permutation(const int n_thread)
{
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm_matrix(
        m_phenotype.rows());
    Eigen::setNbThreads(n_thread);
    Eigen::Index rank = 0;
    // logit_perm can only be true if it is binary trait and user used the
    // --logit-perm flag
    // can always do the following if
    // 1. QT trait (!is_binary)
    // 2. Not require logit perm
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat;
    Eigen::MatrixXd R;
    bool run_glm = true;
    if (!m_binary_trait || !m_perm_info.logit_perm)
    {
        //  first decompose the independent variable to speed up the other
        //  processes
        PQR.compute(m_independent_variables);
        Pmat = PQR.colsPermutation();
        rank = PQR.rank();
        if (rank != m_independent_variables.cols())
        {
            PQR.matrixQR()
                .topLeftCorner(rank, rank)
                .triangularView<Eigen::Upper>()
                .solve(Eigen::MatrixXd::Identity(rank, rank));
        }
        run_glm = false;
    }
    if (n_thread == 1)
    {
        // we will run the single thread function to reduce overhead
        run_null_perm_no_thread(PQR, Pmat, R, run_glm);
    }
    else
    {
        Thread_Queue<std::pair<Eigen::VectorXd, size_t>> set_perm_queue;
        std::thread producer(&PRSice::gen_null_pheno, this,
                             std::ref(set_perm_queue), n_thread - 1);
        std::vector<std::thread> consume_store;
        // we have used one thread as the producer, therefore we need to
        // reduce the number of available thread by 1
        for (int i = 0; i < n_thread - 1; ++i)
        {
            consume_store.push_back(std::thread(
                &PRSice::consume_null_pheno, this, std::ref(set_perm_queue),
                std::ref(PQR), std::ref(Pmat), std::ref(R), run_glm));
        }
        // wait for all the threads to complete their job
        producer.join();
        for (auto&& consume : consume_store) consume.join();
    }
}

void PRSice::run_null_perm_no_thread(
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& PQR,
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType& Pmat,
    const Eigen::MatrixXd& R, const bool run_glm)
{
    // reset the seed for each new threshold such that we will always
    // generate the same phenotpe for each threhsold without us needing to
    // regenerate the PRS
    std::mt19937 rand_gen {m_perm_info.seed};
    // we want to count the number of samples included in the analysis
    const Eigen::Index num_regress_sample = m_phenotype.rows();
    const Eigen::Index p = m_independent_variables.cols();
    const Eigen::Index rank = PQR.rank();
    // we will copy the phenotype into a new vector (Maybe not necessary,
    // but better safe than sorry)
    Eigen::VectorXd perm_pheno = m_phenotype;
    // count the number of loop we've finished so far
    size_t processed = 0;
    // pre-initialize these parameters
    double coefficient, standard_error, r2, obs_p;
    double obs_t = -1;

    Eigen::VectorXd beta, se, effects, fitted, resid;
    Eigen::Index df;
    while (processed < m_perm_info.num_permutation)
    {
        // for quantitative trait, we can directly compute the results
        // without re-computing the decomposition
        perm_pheno = m_phenotype;
        std::shuffle(perm_pheno.data(), perm_pheno.data() + num_regress_sample,
                     rand_gen);
        m_analysis_done++;
        print_progress();
        if (run_glm)
        {
            Regression::glm(perm_pheno, m_independent_variables, obs_p, r2,
                            coefficient, standard_error, 1);
        }
        else
        {
            // directly solve the current phenotype to obtain the required
            // beta
            if (p == rank)
            {
                beta = PQR.solve(perm_pheno);
                fitted = m_independent_variables * beta;
                se = Pmat
                     * PQR.matrixQR()
                           .topRows(p)
                           .triangularView<Eigen::Upper>()
                           .solve(lm::I_p(p))
                           .rowwise()
                           .norm();
            }
            else
            {
                effects = PQR.householderQ().adjoint() * perm_pheno;
                beta = Eigen::VectorXd::Constant(
                    p, std::numeric_limits<double>::quiet_NaN());
                beta.head(rank) = R * effects.head(rank);
                beta = Pmat * beta;
                // create fitted values from effects
                // (can't use X*m_coef if X is rank-deficient)
                effects.tail(num_regress_sample - rank).setZero();
                fitted = PQR.householderQ() * effects;
                se = Eigen::VectorXd::Constant(
                    p, std::numeric_limits<double>::quiet_NaN());
                se.head(rank) = R.rowwise().norm();
                se = Pmat * se;
            }
            // we take the absolute of the T-value as we only concern about
            // the magnitude
            resid = perm_pheno - fitted;
            df = (rank >= 0) ? num_regress_sample - p
                             : num_regress_sample - rank;
            double s = resid.norm() / std::sqrt(double(df));
            se = s * se;
            coefficient = beta(1);
            standard_error = se(1);
        }
        obs_t = std::fabs(coefficient / standard_error);
        m_perm_result[processed] = std::max(obs_t, m_perm_result[processed]);
        // we have finished the current analysis.
        ++processed;
    }
}


void PRSice::gen_null_pheno(Thread_Queue<std::pair<Eigen::VectorXd, size_t>>& q,
                            size_t num_consumer)
{
    size_t processed = 0;
    // we need to reset the seed for each threshold so that the phenotype
    // generated should always be the same for each threshold. This help us
    // avoid needing to repeat reading in the PRS for each permutation
    std::mt19937 rand_gen {m_perm_info.seed};
    Eigen::setNbThreads(1);
    const Eigen::Index num_regress_sample = m_phenotype.rows();
    // this matrix is use for storing the shuffle information
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm_matrix(
        m_phenotype.rows());

    while (processed < m_perm_info.num_permutation)
    {
        // I am not completely sure whether we pass by reference or pass by
        // value to the queue (or if we move this out and cause undefined
        // behaviour). To avoid the aforementioned problems, we reinitialize
        // a new vector for each permutation
        Eigen::VectorXd null_pheno = m_phenotype;
        // then we shuffle it
        std::shuffle(null_pheno.data(), null_pheno.data() + num_regress_sample,
                     rand_gen);
        // and the we will push it to the queue where the consumers will
        // pick up and work on it
        q.emplace(std::make_pair(null_pheno, processed), num_consumer);
        ++m_analysis_done;
        print_progress();
        ++processed;
    }
    // send termination signal to the consumers
    q.completed();
}

void PRSice::consume_null_pheno(
    Thread_Queue<std::pair<Eigen::VectorXd, size_t>>& q,
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& PQR,
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType& Pmat,
    const Eigen::MatrixXd& R, bool run_glm)
{
    const Eigen::Index n = m_phenotype.rows();
    const Eigen::Index p = m_independent_variables.cols();
    const Eigen::Index rank = PQR.rank();
    // to avoid false sharing, all consumer will first store their
    // permutation result in their own vector and only update the master
    // vector at the end of permutation temp_store stores all the T-value
    std::vector<double> temp_store;
    // and temp_index indicate which permutation the T-value corresponse to
    // it is important we don't mix up the permutation order as we're
    // supposed to mimic re-running PRSice N times with different
    // permutation
    std::vector<size_t> temp_index;
    Eigen::VectorXd beta, se, effects, fitted, resid;
    Eigen::Index df;
    std::pair<Eigen::VectorXd, size_t> input;
    double coefficient, standard_error, r2, obs_p;
    double obs_t = -1;
    while (!q.pop(input))
    {
        // as long as we have not received a termination signal, we will
        // continue our processing and should read from the queue
        if (run_glm)
        {
            // the first entry from the queue should be the permuted
            // phenotype and the second entry is the index. We will pass the
            // phenotype for GLM analysis if required
            Regression::glm(std::get<0>(input), m_independent_variables, obs_p,
                            r2, coefficient, standard_error, 1);
        }
        else
        {
            if (p == rank)
            {
                beta = PQR.solve(std::get<0>(input));
                fitted = m_independent_variables * beta;
                se = Pmat
                     * PQR.matrixQR()
                           .topRows(p)
                           .triangularView<Eigen::Upper>()
                           .solve(lm::I_p(p))
                           .rowwise()
                           .norm();
            }
            else
            {
                beta = Eigen::VectorXd::Constant(
                    p, std::numeric_limits<double>::quiet_NaN());
                effects = PQR.householderQ().adjoint() * std::get<0>(input);
                beta.head(rank) = R * effects.head(rank);
                beta = Pmat * beta;
                // create fitted values from effects
                // (can't use X*m_coef if X is rank-deficient)
                effects.tail(n - rank).setZero();
                fitted = PQR.householderQ() * effects;
                se = Eigen::VectorXd::Constant(
                    p, std::numeric_limits<double>::quiet_NaN());
                se.head(rank) = R.rowwise().norm();
                se = Pmat * se;
            }
            coefficient = beta(1);
            resid = std::get<0>(input) - fitted;
            df = (rank >= 0) ? n - p : n - rank;
            double s = resid.norm() / std::sqrt(double(df));
            se = s * se;
            standard_error = se(1);
        }
        obs_t = std::fabs(coefficient / standard_error);
        temp_store.push_back(obs_t);
        temp_index.push_back(std::get<1>(input));
    }
    // once we received the termination signal, we can start propagating the
    // master vector with out content
    std::lock_guard<std::mutex> lock(lock_guard);
    for (size_t i = 0; i < temp_store.size(); ++i)
    {
        double obs_t = temp_store[i];
        auto&& index = temp_index[i];
        // if the t-value in the master vector is lower than our observed t,
        // update it
        if (m_perm_result[index] < obs_t) { m_perm_result[index] = obs_t; }
    }
}

void PRSice::prep_best_output(
    const Genotype& target,
    const std::vector<std::vector<size_t>>& region_membership,
    const std::vector<std::string>& region_name, const size_t max_fid,
    const size_t max_iid, std::unique_ptr<std::ostream>& best_file)
{
    const size_t num_region = region_name.size();
    const size_t num_samples = target.num_sample();
    const long long begin_byte = best_file->tellp();
    (*best_file) << "FID IID In_Regression";
    if (!(num_region > 2)) { (*best_file) << " PRS\n"; }
    for (size_t i = 0; i < region_name.size(); ++i)
    {
        if (i == 1 || region_membership[i].empty()) continue;
        (*best_file) << " " << region_name[i];
    }
    (*best_file) << "\n";
    try
    {
        m_fast_best_output =
            Eigen::MatrixXd::Zero(num_samples, region_name.size());
        m_quick_best = true;
    }
    catch (...)
    {
        m_reporter->report(
            "Warning: Not enough memory to store all best scores "
            "into the memory, will use a slower method to output "
            "the best score file");
        // not enough memory for fast best use the seekg method

        m_quick_best = num_region <= 2;
        const long long end_byte = best_file->tellp();
        assert(end_byte >= begin_byte);
        m_best_file.header_length = end_byte - begin_byte;
        m_best_file.processed_threshold = 0;
        // each numeric output took 12 spaces, then for each output,
        // there is one space next to each
        m_best_file.line_width =
            max_fid /* FID */ + 1LL                    /* space */
            + max_iid                                  /* IID */
            + 1LL /* space */ + 3LL /* Yes/No */ + 1LL /* space */
            + (num_region - 1) /* each region -1 to remove background*/
                  * (m_numeric_width + 1LL /* space */)
            + 1LL /* new line */;
        m_best_file.skip_column_length =
            max_fid + 1LL + max_iid + 1LL + 3LL + 1LL;
        if (!m_quick_best)
        {
            // need to pre-write some info
            std::string best_line;
            for (size_t i_sample = 0; i_sample < num_samples; ++i_sample)
            {
                best_line =
                    target.sample_id(i_sample, " ") + " "
                    + ((target.sample_valid_for_regress(i_sample)) ? "Yes"
                                                                   : "No");
                // TODO: Bug if line width is bigger than what setw can handle
                (*best_file)
                    << std::setfill(' ') << std::setw(m_best_file.line_width)
                    << std::left << best_line << "\n";
            }
            // seek back to front right after header
            best_file->seekp(end_byte, best_file->beg);
        }
    }
    ++m_best_file.line_width;
}

void PRSice::prep_all_score_output(
    const Genotype& target,
    const std::vector<std::vector<size_t>>& region_membership,
    const std::vector<std::string>& region_name, const size_t max_fid,
    const size_t max_iid, std::unique_ptr<std::ostream>& all_score_file)
{
    auto set_thresholds = target.get_set_thresholds();
    unsigned long long total_set_thresholds = 0;
    for (size_t thres = 0; thres < set_thresholds.size(); ++thres)
    {
        // skip bakcground and empty region
        if (thres == 1 || region_membership[thres].empty()) continue;
        if (total_set_thresholds > std::numeric_limits<long long>::max()
                                       - set_thresholds[thres].size())
        {
            throw std::runtime_error(
                "Error: Too many combinations of number of regions and "
                "number "
                "of thresholds, will cause integer overflow.");
        }
        total_set_thresholds += set_thresholds[thres].size();
    }
    const long long begin_byte = all_score_file->tellp();
    (*all_score_file) << "FID IID";
    if (!(region_name.size() > 2))
    {
        // add character in front so that when R parse it, it doesn't add the
        // annoying X and properly treat it as a header
        for (auto& thres : set_thresholds.front())
        { (*all_score_file) << " Pt_" << thres; }
    }
    else
    {
        // don't output all score for background
        // not all set has snps in all threshold
        for (size_t i = 0; i < region_name.size(); ++i)
        {
            if (i == 1 || region_membership[i].empty()) continue;
            for (auto& thres : set_thresholds[i])
            { (*all_score_file) << " " << region_name[i] << "_" << thres; }
        }
    }
    const long long end_byte = all_score_file->tellp();
    // if the line is too long, we might encounter overflow
    assert(end_byte >= begin_byte);
    m_all_file.header_length = end_byte - begin_byte;
    m_all_file.processed_threshold = 0;
    m_all_file.line_width =
        max_fid + 1LL + max_iid + 1LL
        + static_cast<long long>(total_set_thresholds) * (m_numeric_width + 1LL)
        + 1LL;
    m_all_file.skip_column_length = max_fid + max_iid + 2;
    size_t num_sample = target.num_sample();
    // now print out all the empty lines
    std::string name;
    for (size_t i_sample = 0; i_sample < num_sample; ++i_sample)
    {
        name = target.sample_id(i_sample, " ");
        // TODO: Bug if line_width is bigger than what setw can handle
        (*all_score_file) << std::setfill(' ')
                          << std::setw(m_all_file.line_width) << std::left
                          << name << "\n";
    }
    ++m_all_file.line_width;
    all_score_file->seekp(0, all_score_file->beg);
}

std::tuple<double, double>
PRSice::lee_adjustment_factor(const double prevalence)
{
    double num_case = m_phenotype.sum();
    double case_ratio = num_case / static_cast<double>(m_phenotype.rows());
    // the following is from Lee et al A better coefficient paper
    double x = misc::qnorm(1 - prevalence);
    double z = misc::dnorm(x);
    double i2 = z / prevalence;
    double cc = prevalence * (1 - prevalence) * prevalence * (1 - prevalence)
                / (z * z * case_ratio * (1 - case_ratio));
    double theta = i2 * ((case_ratio - prevalence) / (1 - prevalence))
                   * (i2 * ((case_ratio - prevalence) / (1 - prevalence)) - x);
    double e = 1
               - pow(case_ratio, (2 * case_ratio))
                     * pow((1 - case_ratio), (2 * (1 - case_ratio)));
    return {cc * e, cc * e * theta};
}
void PRSice::adjustment_factor(const double prevalence, double& top,
                               double& bottom)
{
    double num_case = m_phenotype.sum();
    double case_ratio = num_case / static_cast<double>(m_phenotype.rows());
    // the following is from Lee et al A better coefficient paper
    double x = misc::qnorm(1 - prevalence);
    double z = misc::dnorm(x);
    double i2 = z / prevalence;
    double cc = prevalence * (1 - prevalence) * prevalence * (1 - prevalence)
                / (z * z * case_ratio * (1 - case_ratio));
    double theta = i2 * ((case_ratio - prevalence) / (1 - prevalence))
                   * (i2 * ((case_ratio - prevalence) / (1 - prevalence)) - x);
    double e = 1
               - pow(case_ratio, (2 * case_ratio))
                     * pow((1 - case_ratio), (2 * (1 - case_ratio)));
    top = cc * e;
    bottom = cc * e * theta;
}

void PRSice::print_summary(const std::string& pheno_name,
                           const double prevalence, const bool has_prevalence,
                           std::vector<size_t>& significant_count,
                           std::unique_ptr<std::ostream>& summary_file)
{
    assert(significant_count.size() == m_significant_store.size());
    for (size_t i = 0; i < significant_count.size(); ++i)
    { significant_count[i] += m_significant_store[i]; }
    double top, bot;
    if (prevalence < 1.0)
    { std::tie(top, bot) = lee_adjustment_factor(prevalence); }
    for (auto&& sum : m_prs_summary)
    {
        (*summary_file) << pheno_name << "\t" << sum.set << "\t"
                        << sum.result.threshold << "\t"
                        << sum.result.r2 - m_null_r2;
        if (has_prevalence && m_binary_trait)
        {
            (*summary_file)
                << "\t"
                << get_adjusted_r2(sum.result.r2, top, bot)
                       - get_adjusted_r2(m_null_r2, top, bot)
                << "\t" << get_adjusted_r2(sum.result.r2, top, bot) << "\t"
                << get_adjusted_r2(m_null_r2, top, bot) << "\t" << prevalence;
        }
        else if (has_prevalence)
        {
            (*summary_file)
                << "\tNA\t" << sum.result.r2 << "\t" << m_null_r2 << "\t-";
        }
        else
        {
            (*summary_file)
                << "\t" << sum.result.r2 << "\t" << m_null_r2 << "\t-";
        }
        // now generate the rest of the output
        (*summary_file) << "\t" << sum.result.coefficient << "\t"
                        << sum.result.se << "\t" << sum.result.p << "\t"
                        << sum.result.num_snp;
        if (m_perm_info.run_set_perm && (sum.result.competitive_p >= 0.0))
        { (*summary_file) << "\t" << sum.result.competitive_p; }
        else if (m_perm_info.run_set_perm)
        {
            (*summary_file) << "\tNA";
        }
        if (m_perm_info.run_perm) (*summary_file) << "\t" << sum.result.emp_p;
        (*summary_file) << "\n";
    }
}


PRSice::~PRSice() {}

void PRSice::get_se_matrix(
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& PQR,
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType& Pmat,
    const Eigen::MatrixXd& Rinv, const Eigen::Index p, const Eigen::Index rank,
    Eigen::VectorXd& se)
{
    if (p == rank)
    {
        se = Pmat
             * PQR.matrixQR()
                   .topRows(p)
                   .triangularView<Eigen::Upper>()
                   .solve(lm::I_p(p))
                   .rowwise()
                   .norm();
    }
    else
    {
        se = Eigen::VectorXd::Constant(
            p, std::numeric_limits<double>::quiet_NaN());
        se.head(rank) = Rinv.rowwise().norm();
        se = Pmat * se;
    }
}

double PRSice::get_coeff_resid_norm(const Regress& decomposed,
                                    const Eigen::VectorXd& prs,
                                    Eigen::VectorXd& beta,
                                    Eigen::VectorXd effects)
{
    const Eigen::Index p = m_independent_variables.cols();
    if (decomposed.rank == p)
    {
        // full rank case
        beta = decomposed.PQR.solve(prs);
        return (prs - decomposed.YCov * beta).norm();
    }
    effects = decomposed.PQR.householderQ().adjoint() * prs;
    beta =
        Eigen::VectorXd::Constant(p, std::numeric_limits<double>::quiet_NaN());
    beta.head(decomposed.rank) =
        decomposed.Rinv * effects.head(decomposed.rank);
    beta = decomposed.Pmat * beta;
    const Eigen::Index num_regress_sample =
        static_cast<Eigen::Index>(m_matrix_index.size());
    effects.tail(num_regress_sample - decomposed.rank).setZero();
    return (prs - decomposed.PQR.householderQ() * effects).norm();
}

double PRSice::get_t_value(const Regress& decomposed,
                           const Eigen::VectorXd& prs, Eigen::VectorXd& beta,
                           Eigen::VectorXd& effects, double& coefficient,
                           double& standard_error)
{
    const Eigen::Index rank = decomposed.rank;
    const Eigen::Index num_regress_sample = decomposed.YCov.rows();
    const Eigen::Index p = m_independent_variables.cols();
    const Eigen::Index df = (rank >= 0) ? num_regress_sample - p
                                        : num_regress_sample - decomposed.rank;
    const double s = get_coeff_resid_norm(decomposed, prs, beta, effects)
                     / std::sqrt(double(df));
    coefficient = beta(1);
    standard_error = (s * decomposed.se)(1);
    return coefficient / standard_error;
}
