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

void PRSice::new_phenotype(Genotype& target)
{
    // remove all content from phenotype vector and independent matrix so
    // that we don't get into trouble also clean out the dictionary which
    // indicate which sample contain valid phenotype
    if (m_prsice_out.is_open()) m_prsice_out.close();
    if (m_all_out.is_open()) m_all_out.close();
    if (m_best_out.is_open()) m_best_out.close();
    m_prsice_out.clear();
    m_all_out.clear();
    m_best_out.clear();
    m_null_r2 = 0.0;
    m_phenotype = Eigen::VectorXd::Zero(0);
    m_independent_variables.resize(0, 0);
    m_sample_with_phenotypes.clear();
    target.reset_in_regression_flag();
    target.reset_std_flag();
    // As m_sample_with_phenotypes is empty, we won't go into the for loop
    // with update to the regression flag or the exclude_std flag, and we
    // will not need the delimitor (as none will be found anyway), so we
    // only need to give the target file to calculate the FID and IID length
    if (m_prs_info.no_regress) { update_sample_included("", false, target); }
}

void PRSice::init_matrix(const std::string& file_name,
                         const std::string& pheno_name,
                         const std::string& delim, const size_t file_idx,
                         const bool ignore_fid, Genotype& target)
{
    gen_pheno_vec(file_name, pheno_name, delim, file_idx, ignore_fid, target);
}
void PRSice::init_matrix(const size_t pheno_index, const std::string& delim,
                         Genotype& target)
{
    // this reset the in_regression flag of all samples
    // don't need to do anything if we don't need to do regression

    // read in phenotype vector
    gen_pheno_vec(m_pheno_info.pheno_file, m_pheno_info.pheno_col[pheno_index],
                  delim, m_pheno_info.pheno_col_idx[pheno_index],
                  m_pheno_info.ignore_fid, target);
    // now that we've got the phenotype, we can start processing the more
    // complicated covariate
    gen_cov_matrix(delim);
    // Update has pheno flag, as some sample might have missing covariates
    update_sample_included(delim, m_pheno_info.binary[pheno_index], target);

    // now we want to calculate the null R2 (if covariates are included)
    double null_r2_adjust = 0.0;
    // get the number of thread available
    const int n_thread = m_prs_info.thread;
    if (m_independent_variables.cols() > 2 && !m_prs_info.no_regress)
    {
        // only do it if we have the correct number of sample
        assert(m_independent_variables.rows() == m_phenotype.rows());
        if (m_pheno_info.binary[pheno_index])
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
}

void PRSice::update_sample_included(const std::string& delim, const bool binary,
                                    Genotype& target)
{
    // we also want to avoid always having to search from
    // m_sample_with_phenotypes. Therefore, we push in the sample index to
    // m_matrix_index
    // as our phenotype vector and independent variable matrix all follow
    // the order of samples appear in the target genotype object, we can
    // safely store the matrix sequentially and the m_matrix should still
    // correspond to the phenotype vector and covariance matrix's order
    // correctly
    m_matrix_index.clear();
    const bool ctrl_std =
        binary && m_prs_info.scoring_method == SCORING::CONTROL_STD;
    for (size_t i_sample = 0; i_sample < target.num_sample(); ++i_sample)
    {
        // got through each sample
        // will be problematic if the fid and iid are too long
        if (target.fid(i_sample).length()
                > std::numeric_limits<long long>::max()
            || target.iid(i_sample).length()
                   > std::numeric_limits<long long>::max())
        {
            throw std::runtime_error(
                "Error: FID / IID are pathologically long");
        }
        // update the in regression flag according to covariate
        if (!m_sample_with_phenotypes.empty())
        {
            auto&& pheno_idx = m_sample_with_phenotypes.find(
                target.sample_id(i_sample, delim));
            if (pheno_idx != m_sample_with_phenotypes.end())
            {
                m_matrix_index.push_back(i_sample);
                // the in regression flag is only use for output
                target.set_in_regression(i_sample);
                // so only standardize samples if they are with valid
                // phenotype
                if (ctrl_std
                    && !misc::logically_equal(
                        m_phenotype[static_cast<Eigen::Index>(
                            pheno_idx->second)],
                        0))
                {
                    // this will not be used for standardization
                    target.exclude_from_std(i_sample);
                }
            }
        }
    }
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
            m_sample_with_phenotypes[id] = i_sample;
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
                    m_sample_with_phenotypes[id] = i_sample;
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

bool PRSice::validate_covariate(const std::string& covariate,
                                const size_t num_factors, const size_t idx,
                                size_t& factor_level_idx,
                                std::vector<size_t>& missing_count

)
{
    if (covariate == "NA" || covariate == "NAN")
    {
        ++missing_count[idx];
        return false;
    }
    // we first check if the factor_level_index is larger than the
    // number of factor. If that is the case, this must not be a
    // factor covaraite.
    // If not, then we check if the current index corresponds to a
    // factor index.
    else if (factor_level_idx >= num_factors
             || idx != m_pheno_info.col_index_of_factor_cov[factor_level_idx])
    {
        // not a factor. Check if this is a valid covariate
        try
        {
            misc::convert<double>(covariate);
        }
        catch (const std::runtime_error&)
        {
            ++missing_count[idx];
            return false;
        }
    }
    // we will iterate the factor_level only if this a factor
    if (factor_level_idx < num_factors)
    {
        factor_level_idx +=
            (idx == m_pheno_info.col_index_of_factor_cov[factor_level_idx]);
    }
    return true;
}

void PRSice::update_sample_matrix(
    const std::vector<size_t>& missing_count,
    std::vector<std::pair<std::string, size_t>>& valid_sample_index)
{
    // helpful to give the overview
    const size_t num_sample = m_sample_with_phenotypes.size();
    const int removed = static_cast<int>(num_sample)
                        - static_cast<int>(valid_sample_index.size());
    std::string message =
        std::to_string(removed) + " sample(s) with invalid covariate:\n\n";
    double portion =
        static_cast<double>(removed) / static_cast<double>(num_sample);

    if (valid_sample_index.size() == 0)
    {
        // if all samples are removed
        size_t cur_cov_index = 0;
        for (auto&& cov : m_pheno_info.col_index_of_cov)
        {
            if (missing_count[cov] == num_sample)
            {
                // inform user which covariate is the culprits
                message.append("Error: "
                               + m_pheno_info.cov_colname[cur_cov_index]
                               + " is invalid, please check it is of the "
                                 "correct format\n");
            }
            ++cur_cov_index;
        }
        m_reporter->report(message);
        throw std::runtime_error("Error: All samples removed due to "
                                 "missingness in covariate file!");
    }
    // provide a warning if too many samples were removed due to covariate
    if (portion > 0.05)
    {
        message.append(
            "Warning: More than " + std::to_string(portion * 100)
            + "% of your samples were removed! "
              "You should check if your covariate file is correct\n");
    }
    m_reporter->report(message);
    // sort the sample index
    // Sorting is required because our ordering follows the covariate
    // file, which does not need to have the same ordering as the
    // target file.
    // Also, this does means that the base factor will be the
    // first factor observed within the covariate file, not
    // the one observed in the first sample in the genotype
    // TODO: If I have time, maybe allow users to select the
    //       base factor? (Would be a pain though)
    std::sort(begin(valid_sample_index), end(valid_sample_index),
              [](std::pair<std::string, size_t> const& t1,
                 std::pair<std::string, size_t> const& t2) {
                  if (std::get<1>(t1) == std::get<1>(t2))
                      return std::get<0>(t1).compare(std::get<0>(t2)) < 0;
                  else
                      return std::get<1>(t1) < std::get<1>(t2);
              });


    // update the m_phenotype and m_independent
    m_sample_with_phenotypes.clear();
    // vector contains the name of samples that we keep
    // and also their original index on m_phenotype
    for (size_t cur_index = 0; cur_index < valid_sample_index.size();
         ++cur_index)
    {
        const std::string name = std::get<0>(valid_sample_index[cur_index]);
        // update sample's index on matrix
        m_sample_with_phenotypes[name] = cur_index;
        size_t original_index = std::get<1>(valid_sample_index[cur_index]);
        if (original_index != cur_index)
        {
            // update the content of the phenotype matrix
            m_phenotype(static_cast<Eigen::Index>(cur_index), 0) =
                m_phenotype(static_cast<Eigen::Index>(original_index), 0);
        }
    }

    m_phenotype.conservativeResize(
        static_cast<Eigen::Index>(valid_sample_index.size()), 1);
}


void PRSice::process_cov_file(
    std::vector<size_t>& cov_start_index,
    std::vector<std::unordered_map<std::string, size_t>>& factor_levels,
    Eigen::Index& num_column, const std::string& delim)
{
    // first, go through the covariate and generate the factor level vector
    std::ifstream cov;
    // we will generate a vector containing the information of all samples
    // with valid covariate. The pair contain Sample Name and the index
    // before removal
    std::vector<std::pair<std::string, size_t>> valid_sample_index;
    // is the token to store the tokenized string
    std::vector<std::string> token;
    // contain the current level of factor
    // at the end, this = number of levels in each factor covariate -1
    std::vector<size_t> current_factor_level(m_pheno_info.factor_cov.size(), 0);
    // is the number of missingness in each covariate
    std::vector<size_t> missing_count(m_pheno_info.col_index_of_cov.back() + 1,
                                      0);
    std::string line, id;
    std::unordered_set<std::string> dup_id_check;
    // is the maximum column index required
    const size_t max_index = m_pheno_info.col_index_of_cov.back() + 1;
    // This is the index for iterating the current_vector_level (reset after
    // each line)
    size_t index = 0, factor_level_index = 0;
    int num_valid = 0, dup_id_count = 0;
    bool valid = true;
    // we initialize the storage facility for the factor levels
    const size_t num_factors = m_pheno_info.factor_cov.size();
    factor_levels.resize(num_factors);
    // open the covariate file
    cov.open(m_pheno_info.cov_file.c_str());
    if (!cov.is_open())
    {
        throw std::runtime_error("Error: Cannot open covariate file: "
                                 + m_pheno_info.cov_file);
    }
    // the number of factor is used to guard against array out of bound


    while (std::getline(cov, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        // we don't need to remove header as we will use the FID/IID to map
        // the samples and unless there's a sample called FID or IID, we
        // should be ok
        token = misc::split(line);
        if (token.size() < max_index)
        {
            throw std::runtime_error(
                "Error: Malformed covariate file, should have at least "
                + std::to_string(max_index) + " columns");
        }
        // check if this sample has a valid phenotype
        id = (m_pheno_info.ignore_fid) ? token[0] : token[0] + delim + token[1];
        if (m_sample_with_phenotypes.find(id) != m_sample_with_phenotypes.end())
        {
            valid = true;
            factor_level_index = 0;
            for (auto&& header : m_pheno_info.col_index_of_cov)
            {
                misc::to_upper(token[header]);
                valid &= validate_covariate(token[header], num_factors, header,
                                            factor_level_index, missing_count);
            }
            if (valid)
            {
                if (dup_id_check.find(id) != dup_id_check.end())
                {
                    // check if there are duplicated ID in the covariance
                    // file
                    ++dup_id_count;
                    continue;
                }
                dup_id_check.insert(id);
                // this is a valid sample, so we want to keep its
                // information in the valid_sample_index first, obtain its
                // current index on the phenotype vector
                index = m_sample_with_phenotypes[id];
                // store the index information
                valid_sample_index.push_back(
                    std::pair<std::string, size_t>(id, index));
                // we reset the factor level index to 0
                factor_level_index = 0;
                ++num_valid;
                for (auto&& factor : m_pheno_info.col_index_of_factor_cov)
                {
                    // now we go through each factor covariate and check if
                    // we have a new level
                    auto&& cur_level = factor_levels[factor_level_index];
                    if (cur_level.find(token[factor]) == cur_level.end())
                    {
                        // if this input is a new level, we will add it to
                        // our factor map
                        cur_level[token[factor]] =
                            current_factor_level[factor_level_index]++;
                    }
                    ++factor_level_index;
                }
            }
        }
    }
    cov.close();

    if (dup_id_count != 0)
    {
        throw std::runtime_error("Error: " + std::to_string(dup_id_count)
                                 + " duplicated IDs in covariate file!\n");
    }
    // Here, we should know the identity of the valid sample and also
    // the factor levels
    // now calculate the number of column required
    // 1 for intercept, 1 for PRS
    // first, we generate the output regarding the covariate inclusion
    std::vector<std::string> all_missing_cov;
    std::string message =
        "Include Covariates:\nName\tMissing\tNumber of levels\n";
    uint32_t total_column = 2;
    const size_t num_sample = m_sample_with_phenotypes.size();
    // reset the factor index
    factor_level_index = 0;
    size_t cur_cov_index = 0;
    size_t num_level = 0;
    // iterate through each covariate
    for (auto&& cov : m_pheno_info.col_index_of_cov)
    {
        cov_start_index.push_back(total_column);
        if (factor_level_index >= m_pheno_info.col_index_of_factor_cov.size()
            || cov != m_pheno_info.col_index_of_factor_cov[factor_level_index])
        {
            // if this is not a factor we will just output the number of
            // missingness and use - for number of factor
            ++total_column;
            message.append(m_pheno_info.cov_colname[cur_cov_index] + "\t"
                           + std::to_string(missing_count[cov]) + "\t-\n");
        }
        else
        {
            // this is a factor
            num_level = factor_levels[factor_level_index++].size();
            // need to add number of level - 1 (as reference level doesn't
            // require additional column) to the total column required
            total_column += num_level - 1;
            // output the missing information
            message.append(m_pheno_info.cov_colname[cur_cov_index] + "\t"
                           + std::to_string(missing_count[cov]) + "\t"
                           + std::to_string(num_level) + "\n");
        }
        ++cur_cov_index;
    }
    // we now output the covariate information
    m_reporter->report(message);
    // now update the m_phenotype vector, removing any sample with missing
    // covariates
    if (valid_sample_index.size() != num_sample && num_sample != 0)
    { update_sample_matrix(missing_count, valid_sample_index); }
    num_column = total_column;
}
/*
void PRSice::gen_cov_matrix(const std::string& cov_file,
                            const std::string& delim)
{
    Eigen::Index num_sample =
        static_cast<Eigen::Index>(m_sample_with_phenotypes.size());
    if (cov_file.empty())
    {
        m_independent_variables = Eigen::MatrixXd::Ones(num_sample, 2);
        return;
    }
}*/

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
                                     const size_t num_valid, Genotype& target)
{
    Eigen::VectorXd new_pheno = Eigen::VectorXd::Zero(num_valid);
    const size_t num_sample = target.num_sample();
    size_t new_matrix_idx = 0, pheno_matrix_idx = 0;
    for (size_t i = 0; i < num_sample; ++i)
    {
        if (target.sample_valid_for_regress(i))
        {
            if (valid_samples[i])
            {
                new_pheno(new_matrix_idx, 0) = m_phenotype(pheno_matrix_idx, 0);
                ++new_matrix_idx;
            }
            else
            {
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
    update_phenotype_matrix(valid_samples, num_valid, target);
    // update phenotype matrix here
    cov_file->clear();
    cov_file->seekg(0, cov_file->beg);
    return factor_levels;
}

void PRSice::gen_cov_matrix(const std::string& delim)
{
    // The size of the map should be informative of the number of sample
    // currently included in the data
    Eigen::Index num_sample =
        static_cast<Eigen::Index>(m_sample_with_phenotypes.size());
    if (m_pheno_info.cov_file.empty())
    {
        m_independent_variables = Eigen::MatrixXd::Ones(num_sample, 2);
        return;
    }
    // convert factor idx to unordered_set for more elegant way of determining
    // if idx is factor
    std::set<size_t> is_factor;
    for (auto&& f : m_pheno_info.col_index_of_factor_cov)
    { is_factor.insert(f); }

    // obtain the index of each covariate the key is the
    // variable name and the value is the index on the matrix

    // need to account for the situation where the same variable name can
    // occur in different covariates
    // As the index are sorted, we can use vector

    // the index of the factor_list is the index of the covariate
    // the key of the nested unorder map is the factor and the value is the
    // factor level (similar to column index)
    std::vector<std::unordered_map<std::string, size_t>> factor_list;

    // an indexor to indicate whcih column should each covariate start from
    // (as there're factor covariates, the some covariates might take up
    // more than one column
    std::vector<size_t> cov_start_index;
    // by default the required number of column for the matrix is
    // intercept+PRS+number of covariate (when there're no factor input)
    Eigen::Index num_column =
        2 + static_cast<Eigen::Index>(m_pheno_info.cov_colname.size());
    // we will perform the first pass to the covariate file which will
    // remove all samples with missing covariate and will also generate the
    // factor level
    std::string message =
        "Processing the covariate file: " + m_pheno_info.cov_file + "\n";
    message.append("==============================\n");
    m_reporter->report(message);
    process_cov_file(cov_start_index, factor_list, num_column, delim);
    // update the number of sample to account for missing covariates
    num_sample = static_cast<Eigen::Index>(m_sample_with_phenotypes.size());
    // initalize the matrix to the desired size
    m_independent_variables =
        Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(num_sample),
                              static_cast<Eigen::Index>(num_column));
    m_independent_variables.col(0).setOnes();
    m_independent_variables.col(1).setOnes();
    // now we only need to fill in the independent matrix without worry
    // about other stuff
    std::ifstream cov;
    cov.open(m_pheno_info.cov_file.c_str());
    if (!cov.is_open())
    {
        throw std::runtime_error("Error: Cannot open covariate file: "
                                 + m_pheno_info.cov_file);
    }
    std::vector<std::string> token;
    std::string line, id;
    size_t max_index = m_pheno_info.col_index_of_cov.back() + 1;
    Eigen::Index cur_index, index, f_level;
    uint32_t cur_factor_index = 0;
    size_t num_factor = m_pheno_info.col_index_of_factor_cov.size(),
           num_cov = m_pheno_info.col_index_of_cov.size();
    while (std::getline(cov, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() < max_index)
        {
            throw std::runtime_error(
                "Error: Malformed covariate file, should contain at least "
                + std::to_string(max_index) + " column!");
        }
        id = (m_pheno_info.ignore_fid) ? token[0] : token[0] + delim + token[1];
        if (m_sample_with_phenotypes.find(id) != m_sample_with_phenotypes.end())
        {
            // Only valid samples will be found in the
            // m_sample_with_phenotypes map structure reset the index to 0
            cur_factor_index = 0;
            // get the row number
            index = static_cast<Eigen::Index>(m_sample_with_phenotypes[id]);
            std::string covariate;
            for (size_t i_cov = 0; i_cov < num_cov; ++i_cov)
            {
                covariate = token[m_pheno_info.col_index_of_cov[i_cov]];
                if (cur_factor_index >= num_factor
                    || m_pheno_info.col_index_of_cov[i_cov]
                           != m_pheno_info
                                  .col_index_of_factor_cov[cur_factor_index])
                {
                    // noraml covariate
                    // we don't need to deal with invalid conversion
                    // situation as those should be taken cared of by the
                    // process_cov_file function
                    m_independent_variables(
                        index,
                        static_cast<Eigen::Index>(cov_start_index[i_cov])) =
                        misc::convert<double>(covariate);
                }
                else
                {
                    // this is a factor
                    // and the level of the current factor is f_level
                    f_level = static_cast<Eigen::Index>(
                        factor_list[cur_factor_index][covariate]);
                    if (f_level != 0)
                    {
                        // if this is not the reference level, we will add 1
                        // to the matrix we need to -1 as the reference
                        // level = 0 and the second level = 1 but for the
                        // second level, it should propagate the first
                        // column
                        cur_index =
                            static_cast<Eigen::Index>(cov_start_index[i_cov])
                            + f_level - 1;
                        m_independent_variables(
                            index, static_cast<Eigen::Index>(cur_index)) = 1;
                    }
                    ++cur_factor_index;
                }
            }
        }
    }

    message = "After reading the covariate file, "
              + std::to_string(m_sample_with_phenotypes.size())
              + " sample(s) included in the analysis\n";
    m_reporter->report(message);
}

void PRSice::reset_result_containers(const Genotype& target,
                                     const size_t region_idx)
{
    m_best_index = -1;
    // m_num_snp_included will store the current number of SNP included.
    // This is the global number and the true number per sample might differ
    // due to missingness. This is only use for display
    m_num_snp_included = 0;
    // m_perm_result stores the result (T-value) from each permutation and
    // is then used for calculation of empirical p value
    m_perm_result.resize(m_perm_info.num_permutation, 0);
    m_prs_results.resize(target.num_threshold(region_idx));
    // set to -1 to indicate not done
    for (auto&& p : m_prs_results)
    {
        p.threshold = -1;
        p.r2 = 0.0;
        p.num_snp = 0;
    }
    // if cur_start_idx == cur_end_idx, this is an empty region
    // initialize score vector
    // this stores the best score for each sample. Ideally, we will only do
    // it once per phenotype as subsequent region should all have the same
    // number of samples
    if (region_idx == 0) m_best_sample_score.resize(target.num_sample());
    std::fill(m_best_sample_score.begin(), m_best_sample_score.end(), 0);
}
bool PRSice::run_prsice(
    const size_t pheno_index, const size_t region_index,
    const std::vector<std::vector<size_t>>& region_membership,
    const bool all_scores, Genotype& target)
{

    // only print out all scores if this is the first phenotype
    const bool print_all_scores = all_scores && pheno_index == 0;
    const size_t num_samples_included = target.num_sample();
    auto set_snp_idx = region_membership[region_index];
    if (set_snp_idx.empty()) { return false; }

    Eigen::initParallel();
    Eigen::setNbThreads(m_prs_info.thread);
    reset_result_containers(target, region_index);
    // now prepare all score
    // current threshold iteration
    // must iterate after each threshold even if no-regress is called
    size_t prs_result_idx = 0;
    double cur_threshold = 0.0;
    // we want to know if user want to obtain the standardized PRS
    // print the progress bar
    print_progress();
    // indicate if this is the first run. If this is the first run,
    // get_score will perform assignment instead of addition
    bool first_run = true;
    std::vector<size_t>::const_iterator start = set_snp_idx.begin();
    while (target.get_score(start, set_snp_idx.cend(), cur_threshold,
                            m_num_snp_included, first_run))
    {
        ++m_analysis_done;
        print_progress();
        if (print_all_scores && pheno_index == 0)
        {
            for (size_t sample = 0; sample < num_samples_included; ++sample)
            {
                // we will calculate the the number of white space we need
                // to skip to reach the current sample + threshold's output
                // position
                const long long loc =
                    m_all_file.header_length
                    + static_cast<long long>(sample)
                          * (m_all_file.line_width + NEXT_LENGTH)
                    + NEXT_LENGTH + m_all_file.skip_column_length
                    + m_all_file.processed_threshold
                    + m_all_file.processed_threshold * m_numeric_width;
                m_all_out.seekp(loc);
                // then we will output the score
                m_all_out << std::setprecision(static_cast<int>(m_precision))
                          << target.calculate_score(sample);
            }
        }
        // we need to then tell the file that we have finish processing one
        // threshold. Next time we output another PRS, it should be output
        // in the column of the next threshold
        ++m_all_file.processed_threshold;
        if (!m_prs_info.no_regress)
        {
            regress_score(target, cur_threshold, m_prs_info.thread, pheno_index,
                          prs_result_idx);
            if (m_perm_info.run_perm)
            {
                permutation(m_prs_info.thread,
                            m_pheno_info.binary[pheno_index]);
            }
        }
        else
        {
            prsice_result cur_result;
            cur_result.threshold = cur_threshold;
            cur_result.num_snp = m_num_snp_included;
            cur_result.p = -1;
            m_prs_results[prs_result_idx] = cur_result;
        }
        ++prs_result_idx;
        first_run = false;
    }
    if (m_quick_best && !m_prs_info.no_regress)
    {
        m_fast_best_output.col(static_cast<Eigen::Index>(region_index)) =
            Eigen::Map<Eigen::VectorXd>(
                m_best_sample_score.data(),
                static_cast<Eigen::Index>(m_best_sample_score.size()));
    }

    // we need to process the permutation result if permutation is required
    if (m_perm_info.run_perm) process_permutations();
    if (!m_prs_info.no_regress & !m_quick_best)
    { slow_print_best(target, pheno_index); }
    return true;
}

void PRSice::slow_print_best(Genotype& target, const size_t pheno_index)
{
    if (m_quick_best) return;
    std::string pheno_name = "";
    if (m_pheno_info.pheno_col.size() > 1)
        pheno_name = m_pheno_info.pheno_col[pheno_index];
    std::string output_prefix = m_prefix;
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    output_prefix.append(".best");
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
        for (size_t sample = 0; sample < target.num_sample(); ++sample)
        {
            long long loc = m_best_file.header_length
                            + static_cast<long long>(sample)
                                  * (m_best_file.line_width + NEXT_LENGTH)
                            + NEXT_LENGTH + m_best_file.skip_column_length
                            + m_best_file.processed_threshold
                            + m_best_file.processed_threshold * m_numeric_width;
            m_best_out.seekp(loc);
            m_best_out << std::setprecision(static_cast<int>(m_precision))
                       << m_best_sample_score[sample];
        }
    }
    // once we finish outputing the result, we need to increment the
    // processed_threshold index such that when we process the next region,
    // we will be writing to the next column instead of overwriting the
    // current column
    ++m_best_file.processed_threshold;
};
void PRSice::print_best(Genotype& target,
                        const std::vector<std::string>& region_name,
                        const size_t pheno_index)
{
    // read in the name of the phenotype. If there's only one phenotype
    // name, we'll do assign an empty string to phenotyp name
    if (!m_quick_best) return;
    std::string pheno_name = "";
    if (m_pheno_info.pheno_col.size() > 1)
        pheno_name = m_pheno_info.pheno_col[pheno_index];
    std::string output_prefix = m_prefix;
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    output_prefix.append(".best");
    // we generate one best score file per phenotype. The reason for this is
    // to not make the best score file too big, and because each phenotype
    // might have different set of sample included in the regression due to
    // missingness
    // we have to overwrite the white spaces with the desired values
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

        m_best_out.close();
        m_best_out.clear();
        m_best_out.open(output_prefix.c_str());
        if (!m_best_out.is_open())
        {
            throw std::runtime_error("Error: Cannot open best file for output: "
                                     + output_prefix);
        }
        m_best_out << "FID IID In_Regression";
        if (region_name.size() <= 2) { m_best_out << " PRS"; }
        else
        {
            for (size_t i = 0; i < region_name.size(); ++i)
            {
                if (i == 1 || region_name[i].empty()) continue;
                m_best_out << " " << region_name[i];
            }
        }
        m_best_out << std::endl;
        for (size_t sample = 0; sample < target.num_sample(); ++sample)
        {
            m_best_out << target.fid(sample) << " " << target.iid(sample) << " "
                       << ((target.sample_in_regression(sample)) ? "Yes" : "No")
                       << std::setprecision(static_cast<int>(m_precision));
            for (Eigen::Index i = 0; i < m_fast_best_output.cols(); ++i)
            {
                if (i == 1 || region_name[i].empty()) continue;
                m_best_out << " " << m_fast_best_output(sample, i);
            }
            m_best_out << "\n";
        }
        // can just close it as we assume we only need to do it once.
        m_best_out.close();
    }
    m_fast_best_output.resize(0, 0);
}

void PRSice::regress_score(Genotype& target, const double threshold,
                           const int thread, const size_t pheno_index,
                           const size_t prs_result_idx)
{
    double r2 = 0.0, r2_adjust = 0.0, p_value = 0.0, coefficient = 0.0,
           se = 0.0;
    const Eigen::Index num_regress_samples =
        static_cast<Eigen::Index>(m_matrix_index.size());
    // should never have m_num_snp_included == 0
    assert(!m_prs_results.empty());
    if (m_num_snp_included == m_prs_results[prs_result_idx].num_snp
        && !m_prs_info.non_cumulate)
    { return; }

    for (Eigen::Index sample_id = 0; sample_id < num_regress_samples;
         ++sample_id)
    {
        // we can directly read in the matrix index from m_matrix_index
        // vector and assign the PRS directly to the indep variable matrix
        m_independent_variables(sample_id, 1) = target.calculate_score(
            m_matrix_index[static_cast<size_t>(sample_id)]);
    }

    if (m_pheno_info.binary[pheno_index])
    {
        // if this is a binary phenotype, we will perform the GLM model
        try
        {
            Regression::glm(m_phenotype, m_independent_variables, p_value, r2,
                            coefficient, se, thread);
        }
        catch (const std::runtime_error& error)
        {
            // This should only happen when the glm doesn't converge.
            // And it actually happen quite often
            fprintf(stderr, "Error: GLM model did not converge!\n");
            fprintf(stderr,
                    "       This is usually caused by small sample\n"
                    "       size or caused by problem in the input file\n"
                    "       If you are certain it is not due to small\n"
                    "       sample size and problematic input, please\n"
                    "       send me the DEBUG files\n");
            std::ofstream debug;
            debug.open("DEBUG");
            debug << m_independent_variables << "\n";
            debug.close();
            debug.open("DEBUG.y");
            debug << m_phenotype << "\n";
            debug.close();
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
        size_t num_include_samples = target.num_sample();
        for (size_t s = 0; s < num_include_samples; ++s)
        {
            // we will have to store the best scores. we cannot directly
            // copy from the m_independent_variable as some samples which
            // might have excluded from the regression model but we still
            // want their PRS.
            m_best_sample_score[s] = target.calculate_score(s);
        }
    }
    // we can now store the prsice_result
    prsice_result cur_result;
    cur_result.threshold = threshold;
    cur_result.r2 = r2;
    cur_result.r2_adj = r2_adjust;
    cur_result.coefficient = coefficient;
    cur_result.p = p_value;
    cur_result.emp_p = -1.0;
    cur_result.num_snp = m_num_snp_included;
    cur_result.se = se;
    cur_result.competitive_p = -1.0;
    m_prs_results[prs_result_idx] = cur_result;
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

void PRSice::permutation(const int n_thread, const bool is_binary)
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
    if (!is_binary || !m_perm_info.logit_perm)
    {
        // if our trait isn't binary or if we don't need to perform logistic
        // regression in our permutation, we will first decompose the
        // independent variable once, therefore speed up the other processes
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
        // we will run teh thread queue where one thread is responsible for
        // generating the shuffled phenotype whereas other threads are
        // responsible for calculating the t statistics
        Thread_Queue<std::pair<Eigen::VectorXd, size_t>> set_perm_queue;
        // For multi-threading we use the producer consumer pattern where
        // the producer will keep random shuffle the phenotypes and the
        // consumers will calculate the t-values All we need to provide to
        // the producer is the number of consumers and the permutation queue
        // use for
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
        m_analysis_done++;
        print_progress();
        processed++;
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

void PRSice::prep_output(const Genotype& target,
                         const std::vector<std::string>& region_name,
                         const size_t pheno_index, const bool all_score)
{
    // As R has a default precision of 7, we will go a bit
    // higher to ensure we use up all precision
    std::string pheno_name = "";
    // > 1 because we don't need to specify the name when there's only one
    // phenotype
    if (m_pheno_info.pheno_col.size() > 1)
        pheno_name = m_pheno_info.pheno_col[pheno_index];
    std::string output_prefix = m_prefix;
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    const std::string out_prsice = output_prefix + ".prsice";
    // all score will be the same for all phenotype anyway, so we will only
    // generate it once
    const std::string out_all = m_prefix + ".all.score";
    const std::string out_best = output_prefix + ".best";
    const size_t num_samples_included = target.num_sample();
    const long long empty_sets =
        std::count_if(region_name.begin(), region_name.end(), empty_name);
    const long long num_region =
        static_cast<long long>(region_name.size()) - empty_sets;
    if (num_region > std::numeric_limits<long long>::max())
    {
        throw std::runtime_error("Error: Too many regions, will cause integer "
                                 "overflow when generating the best file");
    }
    // .prsice output
    // we only need to generate the header for it
    m_prsice_out.open(out_prsice.c_str());
    if (!m_prsice_out.is_open())
    {
        throw std::runtime_error("Error: Cannot open file: " + out_prsice
                                 + " to write");
    }
    // we won't store the empirical p and competitive p output in the
    // prsice file now as that seems like a waste (only one threshold
    // will contain that information, storing that in the summary file
    // should be enough)
    m_prsice_out << "Set\tThreshold\tR2\t";
    // but generate the adjusted R2 if prevalence is provided
    if (!m_pheno_info.prevalence.empty()) m_prsice_out << "R2.adj\t";
    m_prsice_out << "P\tCoefficient\tStandard.Error\tNum_SNP\n";
    if (!m_prs_info.no_regress)
    {
        // .best output
        try
        {
            m_fast_best_output =
                Eigen::MatrixXd::Zero(num_samples_included, region_name.size());
            m_quick_best = true;
        }
        catch (...)
        {
            m_reporter->report(
                "Warning: Not enough memory to store all best scores "
                "into the memory, will use a slower method to output "
                "the best score file");
            // not enough memory for fast best use the seekg method
            m_best_out.open(out_best.c_str());
            if (!m_best_out.is_open())
            {
                throw std::runtime_error("Error: Cannot open file: " + out_best
                                         + " to write");
            }
            std::string header_line = "FID IID In_Regression";
            // The default name of the output should be PRS, but if we are
            // running PRSet, it should be call Base
            if (!(num_region > 2))
                header_line.append(" PRS");
            else
            {
                for (long long i = 0;
                     i < static_cast<long long>(region_name.size()); ++i)
                {
                    if (i == 1 || region_name[static_cast<size_t>(i)].empty())
                        continue;
                    header_line.append(" "
                                       + region_name[static_cast<size_t>(i)]);
                }
                m_quick_best = num_region <= 2;
            }
            // the safetest way to calculate the length we need to skip is
            // to directly count the number of byte involved
            const long long begin_byte = m_best_out.tellp();
            m_best_out << header_line << "\n";
            const long long end_byte = m_best_out.tellp();
            // we now know the exact number of byte the header contain and
            // can correctly skip it acordingly
            assert(end_byte >= begin_byte);
            m_best_file.header_length = end_byte - begin_byte;
            // we will set the processed_threshold information to 0
            m_best_file.processed_threshold = 0;

            // each numeric output took 12 spaces, then for each output,
            // there is one space next to each
            m_best_file.line_width =
                m_max_fid_length /* FID */ + 1LL           /* space */
                + m_max_iid_length                         /* IID */
                + 1LL /* space */ + 3LL /* Yes/No */ + 1LL /* space */
                + (num_region - 1) /* each region -1 to remove background*/
                      * (m_numeric_width + 1LL /* space */)
                + 1LL /* new line */;
            m_best_file.skip_column_length =
                m_max_fid_length + 1LL + m_max_iid_length + 1LL + 3LL + 1LL;
        }
    }

    // also handle all score here
    // but we will only try and generate the all score file when we are
    // dealing with the first phenotype (pheno_index == 0)
    const bool all_scores = all_score && !pheno_index;
    if (all_scores)
    {
        m_all_out.open(out_all.c_str());
        if (!m_all_out.is_open())
        {
            throw std::runtime_error("Cannot open file " + out_all
                                     + " for write");
        }
        // we need to know the number of available thresholds so that we can
        // know how many white spaces we need to pad in
        std::vector<std::set<double>> set_thresholds =
            target.get_set_thresholds();
        unsigned long long total_set_thresholds = 0;
        for (size_t thres = 0; thres < set_thresholds.size(); ++thres)
        {
            if (thres == 1) continue;
            if (total_set_thresholds
                    == std::numeric_limits<unsigned long long>::max()
                || total_set_thresholds > std::numeric_limits<long long>::max())
            {
                throw std::runtime_error(
                    "Error: Too many combinations of number of regions and "
                    "number "
                    "of thresholds, will cause integer overflow.");
            }
            total_set_thresholds += set_thresholds[thres].size();
        }
        // we want the threshold to be in sorted order as we will process
        // the SNPs from the smaller threshold to the highest (therefore,
        // the processed_threshold index should be correct)

        const long long begin_byte = m_all_out.tellp();
        m_all_out << "FID IID";
        // size_t header_length = 3+1+3;
        if (!(region_name.size() > 2))
        {
            for (auto& thres : set_thresholds.front())
            { m_all_out << " " << thres; }
        }
        else
        {
            // don't output all score for background
            // not all set has snps in all threshold
            for (size_t i = 0; i < region_name.size(); ++i)
            {
                if (i == 1) continue;
                for (auto& thres : set_thresholds[i])
                { m_all_out << " " << region_name[i] << "_" << thres; }
            }
        }
        m_all_out << "\n";
        const long long end_byte = m_all_out.tellp();
        // if the line is too long, we might encounter overflow
        assert(end_byte >= begin_byte);
        m_all_file.header_length = end_byte - begin_byte;
        m_all_file.processed_threshold = 0;
        m_all_file.line_width = m_max_fid_length + 1LL + m_max_iid_length + 1LL
                                + static_cast<long long>(total_set_thresholds)
                                      * (m_numeric_width + 1LL)
                                + 1LL;
        m_all_file.skip_column_length = m_max_fid_length + m_max_iid_length + 2;
    }

    // output sample IDs
    std::string best_line;
    std::string name;
    if (all_scores || (!m_prs_info.no_regress && !m_quick_best))
    {
        for (size_t i_sample = 0; i_sample < num_samples_included; ++i_sample)
        {
            name = target.sample_id(i_sample, " ");
            // when we print the best file, we want to also print whether
            // the sample is used in regression or not (so that user can
            // easily reproduce their results)
            if (!m_prs_info.no_regress && !m_quick_best)
            {
                best_line =
                    name + " "
                    + ((target.sample_in_regression(i_sample)) ? "Yes" : "No");
                // we print a line containing m_best_file.line_width white
                // space characters, which we can then overwrite later on,
                // therefore achieving a vertical output
                m_best_out << std::setfill(
                    ' ') << std::setw(static_cast<int>(m_best_file.line_width))
                           << std::left << best_line << "\n";
            }
            if (all_scores)
            {
                m_all_out << std::setfill(' ')
                          << std::setw(static_cast<int>(m_all_file.line_width))
                          << std::left << name << "\n";
            }
        }
    }

    // another one spacing for new line (just to be safe)
    ++m_all_file.line_width;
    ++m_best_file.line_width;
    // don't need to close the files as they will automatically be closed
    // when we move out of the function
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

void PRSice::print_na(const std::string& region_name, const double threshold,
                      const size_t num_snp, const bool has_prevalence)
{
    assert(m_prsice_out.is_open());
    m_prsice_out << region_name << "\t" << threshold << "\tNA\tNA\tNA\tNA\t";
    if (has_prevalence) m_prsice_out << "NA\t";
    m_prsice_out << num_snp << "\n";
}

void PRSice::store_best(const std::string& pheno_name,
                        const std::string& region_name, const double top,
                        const double bottom, const double prevalence,
                        const bool is_base)
{
    if (m_best_index < 0) { return; }
    auto&& best_info =
        m_prs_results[static_cast<std::vector<prsice_result>::size_type>(
            m_best_index)];
    // we will extract the information of the best threshold, store it and
    // use it to generate the summary file in theory though, I should be
    // able to start generating the summary file
    prsice_summary prs_sum;
    prs_sum.pheno = pheno_name;
    prs_sum.set = region_name;
    prs_sum.result = best_info;
    prs_sum.r2_null = m_null_r2;
    prs_sum.top = top;
    prs_sum.bottom = bottom;
    prs_sum.prevalence = prevalence;
    // we don't run competitive testing on the base region
    // therefore we skip region_index == 0 (base is always
    // the first region)
    prs_sum.has_competitive = is_base;
    m_prs_summary.push_back(prs_sum);
    if (best_info.p > 0.1)
        ++m_significant_store[0];
    else if (best_info.p > 1e-5)
        ++m_significant_store[1];
    else
        ++m_significant_store[2];
}
void PRSice::output(const std::vector<std::string>& region_names,
                    const size_t pheno_index, const size_t region_index)
{
    // if prevalence is provided, we'd like to generate calculate the
    // adjusted R2

    const bool has_prevalence = !m_pheno_info.prevalence.empty();
    const bool is_binary = m_pheno_info.binary[pheno_index];
    double top = 1.0, bottom = 1.0, prevalence = -1;
    if (has_prevalence && is_binary)
    {
        size_t num_prev_binary = 0;
        for (size_t i = 0; i < pheno_index; ++i)
        {
            if (m_pheno_info.binary[pheno_index]) ++num_prev_binary;
        }
        adjustment_factor(m_pheno_info.prevalence[num_prev_binary], top,
                          bottom);
    }
    const std::string pheno_name = (m_pheno_info.pheno_col.size() > 1)
                                       ? m_pheno_info.pheno_col[pheno_index]
                                       : "";
    std::string output_prefix = m_prefix;
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    if (m_best_index < 0 && !m_prs_info.no_regress)
    {
        m_reporter->report("Error: No valid PRS for "
                           + region_names[region_index] + "!");
        return;
    }
    // now we know we can generate the prsice file
    // go through every result and output
    for (size_t i = 0; i < m_prs_results.size(); ++i)
    {
        if (m_prs_results[i].threshold < 0 || m_prs_results[i].p < 0)
        {
            print_na(region_names[region_index], m_prs_results[i].threshold,
                     m_prs_results[i].num_snp, has_prevalence);
            continue;
        }
        double full = m_prs_results[i].r2;
        double null = m_null_r2;
        double full_adj = full;
        double null_adj = null;
        if (has_prevalence)
        {
            full_adj = top * full / (1 + bottom * full);
            null_adj = top * null / (1 + bottom * null);
        }
        double r2 = full - null;
        m_prsice_out << region_names[region_index] << "\t"
                     << m_prs_results[i].threshold << "\t" << r2 << "\t";
        if (has_prevalence)
        {
            if (is_binary)
                m_prsice_out << full_adj - null_adj << "\t";
            else
                m_prsice_out << "NA\t";
        }
        m_prsice_out << m_prs_results[i].p << "\t"
                     << m_prs_results[i].coefficient << "\t"
                     << m_prs_results[i].se << "\t" << m_prs_results[i].num_snp
                     << "\n";
        // the empirical p-value will now be excluded from the .prsice
        // output (the "-" isn't that helpful anyway)
    }
    store_best(pheno_name, region_names[region_index], top, bottom, prevalence,
               region_index == 0);
}

void PRSice::summarize()
{
    // we need to know if we are going to write "and" in the output, thus
    // need a flag to indicate if there are any previous outputs
    bool has_previous_output = false;
    // we will output a short summary file
    std::string message = "There are ";
    if (m_significant_store[0] != 0)
    {
        message.append(
            misc::to_string(m_significant_store[0])
            + " region(s)/phenotype(s) with p-value > 0.1 (\033[1;31mnot "
              "significant\033[0m);");
        has_previous_output = true;
    }
    if (m_significant_store[1] != 0)
    {
        if (m_significant_store[2] == 0 && has_previous_output)
        { message.append(" and "); }
        message.append(
            misc::to_string(m_significant_store[1])
            + " region(s) with p-value between "
              "0.1 and 1e-5 (\033[1;31mmay not be significant\033[0m);");
        has_previous_output = true;
    }
    if (m_significant_store[2] != 0)
    {
        if (has_previous_output) message.append(" and ");
        message.append(std::to_string(m_significant_store[2])
                       + " region(s) with p-value less than 1e-5.");
    }
    if (!has_previous_output)
    {
        message.append(
            " Please note that these results are inflated due to the "
            "overfitting inherent in finding the best-fit "
            "PRS (but it's still best to find the best-fit PRS!).\n"
            "You can use the --perm option (see manual) to calculate "
            "an empirical P-value.");
    }
    m_reporter->report(message);
    // now we generate the output file
    std::string out_name = m_prefix + ".summary";
    std::ofstream out;
    out.open(out_name.c_str());
    if (!out.is_open())
    {
        std::string error_message =
            "Error: Cannot open file: " + out_name + " to write";
        throw std::runtime_error(error_message);
    }
    const bool has_prevalence = !m_pheno_info.prevalence.empty();
    out << "Phenotype\tSet\tThreshold\tPRS.R2";
    if (has_prevalence) { out << "\tPRS.R2.adj"; }
    out << "\tFull.R2\tNull."
           "R2\tPrevalence\tCoefficient\tStandard.Error\tP\tNum_SNP";
    if (m_perm_info.run_set_perm) out << "\tCompetitive.P";
    if (m_perm_info.run_perm) out << "\tEmpirical-P";
    out << "\n";
    for (auto&& sum : m_prs_summary)
    {
        out << ((sum.pheno.empty()) ? "-" : sum.pheno) << "\t" << sum.set
            << "\t" << sum.result.threshold << "\t"
            << sum.result.r2 - sum.r2_null;
        // by default, phenotypethat doesn't have the prevalence
        // information will have a prevalence of -1
        if (sum.prevalence > 0)
        {
            // calculate the adjusted R2 for binary traits
            double full = sum.result.r2;
            double null = sum.r2_null;
            full = sum.top * full / (1 + sum.bottom * full);
            null = sum.top * null / (1 + sum.bottom * null);
            out << "\t" << full - null << "\t" << full << "\t" << null << "\t"
                << sum.prevalence;
        }
        else if (has_prevalence)
        {
            // and replace the R2 adjust by NA if the sample doesn't have
            // prevalence (i.e. quantitative trait)
            out << "\tNA\t" << sum.result.r2 << "\t" << sum.r2_null << "\t"
                << sum.prevalence;
        }
        else
        {
            out << "\t" << sum.result.r2 << "\t" << sum.r2_null << "\t-";
        }
        // now generate the rest of the output
        out << "\t" << sum.result.coefficient << "\t" << sum.result.se << "\t"
            << sum.result.p << "\t" << sum.result.num_snp;
        if (m_perm_info.run_set_perm && (sum.result.competitive_p >= 0.0))
        { out << "\t" << sum.result.competitive_p; }
        else if (m_perm_info.run_set_perm)
        {
            out << "\tNA";
        }
        if (m_perm_info.run_perm) out << "\t" << sum.result.emp_p;
        out << "\n";
    }
    out.close();
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
