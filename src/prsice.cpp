// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. Oâ€™Reilly
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
void PRSice::pheno_check(const std::string& pheno_file,
                         const std::vector<std::string>& pheno_header,
                         const std::vector<bool>& is_binary, Reporter& reporter)
{
    std::string message = "";
    if (pheno_file.empty()) {
        // user did not provide a phenotype file, will use the information on
        // the fam file or sample file as input
        pheno_info.use_pheno = false;
        // we will still need to know if the sample is binary or not. As no
        // phenotype file is provided, there can only be one phenotype, thus the
        // first entry of the binary_target should correspond to the binary
        // status of our phenotype
        pheno_info.binary.push_back(is_binary[0]);
    }
    else
    {
        // user provided the phenotype file
        std::ifstream pheno;
        pheno.open(pheno_file.c_str());
        if (!pheno.is_open()) {
            std::string error_message =
                "Cannot open phenotype file: " + pheno_file;
            throw std::runtime_error(error_message);
        }
        std::string line;
        // read in the header line to check if the phenotype is here
        std::getline(pheno, line);
        if (line.empty()) {
            throw std::runtime_error(
                "Cannot have empty header line for phenotype file!");
        }
        pheno.close();
        misc::trim(line);
        std::vector<std::string> col = misc::split(line);
        // we need at least 2 columns (IID + Phenotype)
        // or 3 if m_ignore_fid != T
        if (col.size() < static_cast<std::vector<std::string>::size_type>(
                             2 + !m_ignore_fid))
        {
            throw std::runtime_error(
                "Error: Not enough column in Phenotype file. "
                "Have you use the --ignore-fid option");
        }
        // the first column must be considered as the ID
        std::string sample_id = col[0];
        // if we should not ignore the fid, we should then use both the first
        // and second column as our sample ID
        if (!m_ignore_fid && col.size() > 1) sample_id.append("+" + col[1]);
        message.append("Check Phenotype file: " + pheno_file + "\n");
        message.append("Column Name of Sample ID: " + sample_id + "\n");
        message.append("Note: If the phenotype file does not contain a header, "
                       "the column name will be displayed as the Sample ID "
                       "which is ok.\n");
        bool found = false;
        // now we want to check if the file contain the phenotype header
        std::unordered_map<std::string, bool> dup_col;
        if (pheno_header.size() == 0) {
            // user did not provide a phenotype name. We will therefore simply
            // use the first entry
            pheno_info.use_pheno = true;
            // The index of the phenotype will be 1 (IID Phenotype) or 2 (FID
            // IID Pheno)
            pheno_info.col.push_back(1 + !m_ignore_fid);
            // And we will use the default name (as file without header is
            // still consider as a valid input)
            pheno_info.name.push_back("Phenotype");
            // phenotype order correspond to the phenotype name in the phenotype
            // name vector. Here we place 0 as a place holder
            pheno_info.order.push_back(0);
            pheno_info.binary.push_back(is_binary[0]);
            message.append(
                "Phenotype Name: "
                + col[static_cast<std::vector<std::string>::size_type>(
                      pheno_info.col.back())]
                + "\n");
        }
        else
        {
            // if user provide the phenotype names, we will go through the
            // phenotype header and try to identify the corresponding phenotype
            // index
            for (std::vector<std::string>::size_type i_pheno = 0;
                 i_pheno < pheno_header.size(); ++i_pheno)
            {
                if (dup_col.find(pheno_header[i_pheno]) == dup_col.end()) {
                    // we will ignore any duplicate phenotype input.
                    // it should still be ok as the binary_target should have
                    // the same length as the phenotype column name and we will
                    // store the binary target information
                    found = false;
                    dup_col[pheno_header[i_pheno]] = true;
                    // start from 1+!m_ignore_fid to skip the iid and fid part
                    for (std::vector<std::string>::size_type i_column =
                             1 + !m_ignore_fid;
                         i_column < col.size(); ++i_column)
                    {
                        // now go through each column of the input file to
                        // identify the index
                        // NOTE: If there are multiple column with the same
                        // column name that the user required, we will terminate
                        // as we don't know which one to use
                        if (col[i_column] == (pheno_header[i_pheno])) {
                            if (found) {
                                std::string error_message =
                                    "Error: Multiple Column of your phenotype "
                                    "file matches with the required phenotype "
                                    "name: "
                                    + pheno_header[i_pheno];
                                throw std::runtime_error(error_message);
                            }
                            found = true;
                            pheno_info.use_pheno = true;
                            // store the column index
                            pheno_info.col.push_back(
                                static_cast<int>(i_column));
                            // store the phenotype name
                            pheno_info.name.push_back(pheno_header[i_pheno]);
                            // store the order of the phenotype name in the
                            // pheno-col
                            pheno_info.order.push_back(
                                static_cast<int>(i_pheno));
                            // store the binary information of teh phenotype
                            pheno_info.binary.push_back(is_binary[i_pheno]);
                            break;
                        }
                    }
                    if (!found) {
                        message.append(
                            "Phenotype: " + pheno_header[i_pheno]
                            + " cannot be found in phenotype file\n");
                    }
                }
            }
        }
    }
    size_t num_pheno = (pheno_info.use_pheno) ? pheno_info.col.size() : 1;
    message.append("There are a total of " + std::to_string(num_pheno)
                   + " phenotype to process\n");
    reporter.report(message);
}

void PRSice::init_matrix(const Commander& c_commander, const size_t pheno_index,
                         const std::string& delim, Genotype& target,
                         Reporter& reporter)
{
    // reset the null R2 to 0 (different phenotype might lead to different
    // missingness, thus a different null model R2)
    m_null_r2 = 0.0;
    // remove all content from phenotype vector and independent matrix so that
    // we don't get into trouble
    m_phenotype = Eigen::VectorXd::Zero(0);
    m_independent_variables.resize(0, 0);
    // also clean out the dictionary which indicate which sample contain valid
    // phenotype
    m_sample_with_phenotypes.clear();
    // check if we want to perform regression
    const bool no_regress = c_commander.no_regress();
    // get the phenotype file name
    const std::string pheno_file = c_commander.pheno_file();
    // also get the output file prefix
    const std::string output_name = c_commander.out();

    // this reset the in_regression flag of all samples
    target.reset_in_regression_flag();
    // don't need to do anything if we don't need to do regression
    if (no_regress) {
        update_sample_included(delim, target);
        return;
    }
    // read in phenotype vector
    gen_pheno_vec(target, pheno_file, pheno_index, delim, reporter);
    // now that we've got the phenotype, we can start processing the more
    // complicated covariate
    gen_cov_matrix(c_commander.get_cov_file(), c_commander.get_cov_name(),
                   c_commander.get_cov_index(),
                   c_commander.get_factor_cov_index(), delim, reporter);
    // NOTE: After gen_cov_matrix, the has_pheno flag in m_sample_names is no
    // longer correct as we have not updated that to account for invalid
    // covariates.

    // now inform PRSice which samples should be included in the regression
    // model (mainly for best score output)
    update_sample_included(delim, target);

    // now we want to calculate the null R2 (if covariates are included)
    double null_r2_adjust = 0.0;
    // get the number of thread available
    size_t n_thread = c_commander.thread();
    if (m_independent_variables.cols() > 2 && !no_regress) {
        // only do it if we have the correct number of sample
        assert(m_independent_variables.rows() == m_phenotype.rows());
        if (c_commander.is_binary(static_cast<size_t>(pheno_index))) {
            // ignore the first column
            // this is ok as both the first column (intercept) and the
            // second column (PRS) is currently 1
            Regression::glm(m_phenotype,
                            m_independent_variables.topRightCorner(
                                m_independent_variables.rows(),
                                m_independent_variables.cols() - 1),
                            m_null_p, m_null_r2, m_null_coeff, m_null_se, 25,
                            n_thread, true);
        }
        else
        {
            // ignore the first column
            // and perform linear regression
            Regression::linear_regression(
                m_phenotype,
                m_independent_variables.topRightCorner(
                    m_independent_variables.rows(),
                    m_independent_variables.cols() - 1),
                m_null_p, m_null_r2, null_r2_adjust, m_null_coeff, m_null_se,
                n_thread, true);
        }
    }
}

void PRSice::update_sample_included(const std::string& delim, Genotype& target)
{
    // this is a bit tricky. The reason we need to calculate the max fid and iid
    // length is so that we can generate the best file and all score file
    // vertically by replacing pre-placed space characters with the desired
    // number
    m_max_fid_length = 3;
    m_max_iid_length = 3;
    // we also want to avoid always having to search from
    // m_sample_with_phenotypes. Therefore, we push in the sample index to
    // m_matrix_index
    // as our phenotype vector and independent variable matrix all follow the
    // order of samples appear in the target genotype object, we can safely
    // store the matrix sequentially and the m_matrix should still correspond to
    // the phenotype vector and covariance matrix's order correctly
    m_matrix_index.clear();
    int32_t fid_length, iid_length;
    for (size_t i_sample = 0; i_sample < target.num_sample(); ++i_sample) {
        // got through each sample
        fid_length = static_cast<int32_t>(target.fid(i_sample).length());
        iid_length = static_cast<int32_t>(target.iid(i_sample).length());
        if (m_max_fid_length < fid_length) m_max_fid_length = fid_length;
        if (m_max_iid_length < iid_length) m_max_iid_length = iid_length;
        // update the in regression flag according to covariate
        if (m_sample_with_phenotypes.find(target.sample_id(i_sample, delim))
            != m_sample_with_phenotypes.end())
        {
            m_matrix_index.push_back(i_sample);
            // the in regression flag is only use for output
            target.set_in_regression(i_sample);
        }
    }
}


void PRSice::gen_pheno_vec(Genotype& target, const std::string& pheno_file_name,
                           const size_t pheno_index, const std::string& delim,
                           Reporter& reporter)
{

    // reserve the maximum size (All samples)
    // check if the phenotype is binary or not
    const bool binary = pheno_info.binary[pheno_index];
    const size_t sample_ct = target.num_sample();
    std::string line;
    int max_pheno_code = 0;
    int num_case = 0;
    int num_control = 0;
    size_t invalid_pheno = 0;
    size_t num_not_found = 0;
    size_t sample_index_ct = 0;
    // we will first store the phenotype into the double vector and then later
    // use this to construct the matrix
    std::vector<double> pheno_store;
    pheno_store.reserve(sample_ct);
    std::string pheno_name = "Phenotype";
    std::string id;

    // check if input is sensible
    double first_pheno = 0.0;
    bool more_than_one_pheno = false;

    if (pheno_info.use_pheno) // use phenotype file
    {
        // read in the phenotype index
        const size_t pheno_col_index =
            static_cast<size_t>(pheno_info.col[pheno_index]);
        // and get the phenotype name
        pheno_name = pheno_info.name[pheno_index];
        // now read in the phenotype file
        std::ifstream pheno_file;
        // check if the file is open
        pheno_file.open(pheno_file_name.c_str());
        if (!pheno_file.is_open()) {
            std::string error_message =
                "Cannot open phenotype file: " + pheno_file_name;
            throw std::runtime_error(error_message);
        }

        // we first store everything into a map. This allow the phenotype and
        // genotype file to have completely different ordering and allow
        // different samples to be included in each file
        std::unordered_map<std::string, std::string> phenotype_info;
        std::vector<std::string> token;
        // do not remove header line as that won't match anyway
        while (std::getline(pheno_file, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            // Check if we have the minimal required column number
            if (token.size() < pheno_col_index + 1) {
                std::string error_message =
                    "Malformed pheno file, should contain at least "
                    + misc::to_string(pheno_col_index + 1)
                    + " columns. "
                      "Have you use the --ignore-fid option?";
                throw std::runtime_error(error_message);
            }
            // read in the sample ID
            // TODO: potential problem with BGEN. Might want to allow for
            //       delim here
            id = (m_ignore_fid) ? token[0] : token[0] + delim + token[1];
            // and store the information into the map
            if (phenotype_info.find(id) != phenotype_info.end()) {
                std::string error_message = "Error: Duplicated sample ID in "
                                            "phenotype file: "
                                            + id
                                            + ". Please "
                                              "check if your input is correct!";
                throw std::runtime_error(error_message);
            }
            phenotype_info[id] = token[pheno_col_index];
        }
        pheno_file.close();
        for (size_t i_sample = 0; i_sample < sample_ct; ++i_sample) {
            // now we go through all the samples
            // get the sample ID from the genotype object
            id = target.sample_id(i_sample, delim);
            if (phenotype_info.find(id) != phenotype_info.end()
                && phenotype_info[id] != "NA" && target.is_founder(i_sample))
            {
                // if this sample is found in the phenotype file, and the
                // phenotype isn't NA, and the sample is founder (or
                // keep-founder is used)
                try
                {
                    if (binary) {
                        // if trait is binary
                        // we first convert it to a temporary
                        int temp = misc::convert<int>(phenotype_info[id]);
                        // so taht we can check if the input is valid
                        if (temp >= 0 && temp <= 2) {
                            pheno_store.push_back(temp);
                            // we will also check what is the maximum phenotype
                            // code (1 or 2)
                            // this should happen relatively infrequently so
                            // branch prediction should be rather accurate?
                            if (max_pheno_code < temp) max_pheno_code = temp;
                            // for now, we assume the coding is 0/1
                            num_case += (temp == 1);
                            num_control += (temp == 0);
                        }
                        else
                        {
                            // the phenotype of this sample is invalid
                            throw std::runtime_error(
                                "Invalid binary phenotype format!");
                        }
                    }
                    else
                    {
                        // we will directly push_back the phenotype, if it is
                        // not convertable, it will go into the catch
                        pheno_store.push_back(
                            misc::convert<double>(phenotype_info[id]));
                        if (pheno_store.size() == 1) {
                            // this is the first entrance
                            first_pheno = pheno_store[0];
                        }
                        else if (!more_than_one_pheno
                                 && !misc::logically_equal(first_pheno,
                                                           pheno_store.back()))
                        {
                            // if we found something different from previous
                            // input, then we will set more than one pheno as
                            // true
                            more_than_one_pheno = true;
                        }
                    }
                    // we indicate we have this phenotype
                    // sample_index_ct should currently be representative of the
                    // sample index because we use the suffix++
                    m_sample_with_phenotypes[id] = sample_index_ct++;
                }
                catch (...)
                {
                    invalid_pheno++;
                }
            }
            else
            {
                // we cannot find this sample in the phenotype file / or that we
                // don't want to
                // TODO: Differentiate not include & not include for regression
                num_not_found++;
            }
        }
    }
    else
    {
        // No phenotype file is provided
        // Use information from the fam file directly
        for (size_t i_sample = 0; i_sample < sample_ct; ++i_sample) {

            if (target.pheno_is_na(i_sample) || !target.is_founder(i_sample)) {
                // it is ok to skip NA as default = sample.has_pheno = false
                continue;
            }
            try
            {
                if (binary) {
                    // try to convert the input to int (we stored it as string)
                    int temp = misc::convert<int>(target.pheno(i_sample));
                    if (temp >= 0 && temp <= 2) {
                        // again, check if the input is within reasonable range
                        pheno_store.push_back(temp);
                        if (max_pheno_code < temp) max_pheno_code = temp;
                        // assume 0/1 encoding first
                        num_case += (temp == 1);
                        num_control += (temp == 0);
                    }
                    else
                    {
                        // this isn't a valid binary phenotype
                        throw std::runtime_error(
                            "Invalid binary phenotype format!");
                    }
                }
                else
                {
                    pheno_store.push_back(
                        misc::convert<double>(target.pheno(i_sample)));
                    if (pheno_store.size() == 1) {
                        // this is the first entrance
                        first_pheno = pheno_store[0];
                    }
                    else if (!more_than_one_pheno
                             && !misc::logically_equal(first_pheno,
                                                       pheno_store.back()))
                    {
                        // if we found something different from previous
                        // input, then we will set more than one pheno as
                        // true
                        more_than_one_pheno = true;
                    }
                }
                // we indicate we have this phenotype
                // sample_index_ct should currently be representative of the
                // sample index because we use the suffix++
                m_sample_with_phenotypes[target.sample_id(i_sample, delim)] =
                    sample_index_ct++;
            }
            catch (const std::runtime_error&)
            {
                invalid_pheno++;
            }
        }
    }

    std::string message = "";
    message = pheno_name + " is a ";
    if (binary) {
        message.append("binary phenotype\n");
    }
    else
    {
        message.append("continuous phenotype\n");
    }
    if (num_not_found != 0) {
        message.append(std::to_string(num_not_found)
                       + " sample(s) without phenotype\n");
    }
    if (invalid_pheno != 0) {
        message.append(std::to_string(invalid_pheno)
                       + " sample(s) with invalid phenotype\n");
    }

    if (num_not_found == sample_ct) {
        // it is also possible that the only sample that were found in the
        // phenotype file are the non-founder
        message.append(
            "None of the target samples were found in the phenotype file. ");
        if (m_ignore_fid) {
            message.append(
                "Maybe the first column of your phenotype file is the FID?");
        }
        else
        {
            message.append(
                "Maybe your phenotype file doesn not contain the FID?\n");
            message.append("Might want to consider using --ignore-fid\n");
        }
        message.append("Or it is possible that only non-founder sample contain "
                       "the phenotype information and you did not use "
                       "--nonfounders?\n");
        reporter.report(message);
        throw std::runtime_error("Error: No sample left");
    }
    if (invalid_pheno == sample_ct) {
        message.append("Error: All sample has invalid phenotypes!");
        reporter.report(message);
        throw std::runtime_error("Error: No sample left");
    }
    if (!binary && !more_than_one_pheno) {
        message.append("Only one phenotype value detected");
        if (misc::logically_equal(first_pheno, -9)) {
            message.append(" and they are all -9");
        }
        reporter.report(message);
        throw std::runtime_error("Not enough valid phenotype");
    }
    // finished basic logs
    // we now check if the binary encoding is correct
    bool error = false;
    if (max_pheno_code > 1 && binary) {
        // this is likely code in 1/2
        num_case = 0;
        num_control = 0;
        for (auto&& pheno : pheno_store) {
            pheno--;
            if (pheno < 0) {
                error = true;
            }
            else
                (misc::logically_equal(pheno, 1)) ? num_case++ : num_control++;
        }
    }
    if (error) {
        reporter.report(message);
        throw std::runtime_error(
            "Mixed encoding! Both 0/1 and 1/2 encoding found!");
    }
    if (pheno_store.size() == 0) {
        reporter.report(message);
        throw std::runtime_error("No phenotype presented");
    }
    // now store the vector into the m_phenotype vector
    m_phenotype = Eigen::Map<Eigen::VectorXd>(
        pheno_store.data(), static_cast<Eigen::Index>(pheno_store.size()));


    if (binary) {
        message.append(std::to_string(num_control) + " control(s)\n");
        message.append(std::to_string(num_case) + " case(s)\n");
        if (num_control == 0)
            throw std::runtime_error("There are no control samples");
        if (num_case == 0) throw std::runtime_error("There are no cases");
    }
    else
    {
        message.append(std::to_string(m_phenotype.rows())
                       + " sample(s) with valid phenotype\n");
    }
    reporter.report(message);
}


void PRSice::process_cov_file(
    const std::string& cov_file, const std::vector<uint32_t>& factor_cov_index,
    std::vector<size_t>& cov_start_index,
    const std::vector<uint32_t>& cov_index,
    const std::vector<std::string>& cov_name,
    std::vector<std::unordered_map<std::string, size_t>>& factor_levels,
    size_t& num_column, const std::string& delim, Reporter& reporter)
{
    // first, go through the covariate and generate the factor level vector
    std::ifstream cov;
    // we will generate a vector containing the information of all samples with
    // valid covariate. The pair contain Sample Name and the index before
    // removal
    std::vector<std::pair<std::string, size_t>> valid_sample_index;
    // is the token to store the tokenized string
    std::vector<std::string> token;
    // contain the current level of factor
    // at the end, this = number of levels in each factor covariate -1
    std::vector<size_t> current_factor_level(factor_cov_index.size(), 0);
    // is the number of missingness in each covariate
    std::vector<size_t> missing_count(cov_index.back() + 1, 0);
    std::string line, id;
    std::unordered_set<std::string> dup_id_check;
    // is the maximum column index required
    const size_t max_index = cov_index.back() + 1;
    // This is the index for iterating the current_vector_level (reset after
    // each line)
    size_t factor_level_index = 0;
    int num_valid = 0, index = 0, dup_id_count = 0;
    bool valid = true;
    // we initialize the storage facility for the factor levels
    factor_levels.resize(factor_cov_index.size());
    // open the covariate file
    cov.open(cov_file.c_str());
    if (!cov.is_open()) {
        throw std::runtime_error("Error: Cannot open covariate file: "
                                 + cov_file);
    }
    // the number of factor is used to guard against array out of bound
    const size_t num_factors = factor_cov_index.size();

    while (std::getline(cov, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        // we don't need to remove header as we will use the FID/IID to map
        // the samples and unless there's a sample called FID or IID, we should
        // be ok
        token = misc::split(line);
        if (token.size() < max_index) {
            throw std::runtime_error(
                "Error: Malformed covariate file, should have at least "
                + std::to_string(max_index) + " columns");
        }
        // check if this sample has a valid phenotype
        id = (m_ignore_fid) ? token[0] : token[0] + delim + token[1];
        if (m_sample_with_phenotypes.find(id) != m_sample_with_phenotypes.end())
        {
            valid = true;
            factor_level_index = 0;
            for (auto&& header : cov_index) {
                std::transform(token[header].begin(), token[header].end(),
                               token[header].begin(), ::toupper);
                if (token[header] == "NA") {
                    // this sample has a missing covariate
                    valid = false;
                    ++missing_count[header];
                }
                // we first check if the factor_level_index is larger than the
                // number of factor. If that is the case, this must not be a
                // factor covaraite.
                // If not, then we check if the current index corresponds to a
                // factor index.
                else if (factor_level_index >= num_factors
                         || header != factor_cov_index[factor_level_index])
                {
                    // not a factor. Check if this is a valid covariate
                    try
                    {
                        misc::convert<double>(token[header]);
                    }
                    catch (const std::runtime_error&)
                    {
                        valid = false;
                        ++missing_count[header];
                    }
                }
                // we will iterate the factor_level only if this a factor
                if (factor_level_index < num_factors) {
                    factor_level_index +=
                        (header == factor_cov_index[factor_level_index]);
                }
            }
            if (valid) {
                if (dup_id_check.find(id) != dup_id_check.end()) {
                    // check if there are duplicated ID in the covariance file
                    dup_id_count++;
                    continue;
                }
                dup_id_check.insert(id);
                // this is a valid sample, so we want to keep its information in
                // the valid_sample_index
                // first, obtain its current index on the phenotype vector
                index = static_cast<int>(m_sample_with_phenotypes[id]);
                // store the index information
                valid_sample_index.push_back(
                    std::pair<std::string, size_t>(id, index));
                // we reset the factor level index to 0
                factor_level_index = 0;
                ++num_valid;
                for (auto&& factor : factor_cov_index) {
                    // now we go through each factor covariate and check if we
                    // have a new level
                    auto&& cur_level = factor_levels[factor_level_index];
                    if (cur_level.find(token[factor]) == cur_level.end()) {
                        // if this input is a new level, we will add it to our
                        // factor map
                        cur_level[token[factor]] =
                            current_factor_level[factor_level_index]++;
                    }
                    ++factor_level_index;
                }
            }
        }
    }
    cov.close();

    if (dup_id_count != 0) {
        std::string message = "Error: " + std::to_string(dup_id_count)
                              + " duplicated IDs in covariate file!\n";
        throw std::runtime_error(message);
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
    size_t num_sample = m_sample_with_phenotypes.size();
    // reset the factor index
    factor_level_index = 0;
    size_t cur_cov_index = 0;
    size_t num_level = 0;
    // iterate through each covariate
    for (auto&& cov : cov_index) {
        cov_start_index.push_back(total_column);
        if (factor_level_index >= factor_cov_index.size()
            || cov != factor_cov_index[factor_level_index])
        {
            // if this is not a factor we will just output the number of
            // missingness and use - for number of factor
            ++total_column;
            message.append(cov_name[cur_cov_index] + "\t"
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
            message.append(cov_name[cur_cov_index] + "\t"
                           + std::to_string(missing_count[cov]) + "\t"
                           + std::to_string(num_level) + "\n");
        }
        ++cur_cov_index;
    }
    // we now output the covariate information
    reporter.report(message);
    // now update the m_phenotype vector, removing any sample with missing
    // covariates
    if (valid_sample_index.size() != num_sample && num_sample != 0) {
        // helpful to give the overview
        int removed = static_cast<int>(num_sample)
                      - static_cast<int>(valid_sample_index.size());
        message =
            std::to_string(removed) + " sample(s) with invalid covariate:\n\n";
        double portion =
            static_cast<double>(removed) / static_cast<double>(num_sample);
        if (valid_sample_index.size() == 0) {
            // if all samples are removed
            cur_cov_index = 0;
            for (auto&& cov : cov_index) {
                if (missing_count[cov] == num_sample) {
                    // inform user which covariate is the culprits
                    message.append("Error: " + cov_name[cur_cov_index]
                                   + " is invalid, please check it is of the "
                                     "correct format\n");
                }
                ++cur_cov_index;
            }
            reporter.report(message);
            throw std::runtime_error("Error: All samples removed due to "
                                     "missingness in covariate file!");
        }
        // provide a warning if too many samples were removed due to covariate
        if (portion > 0.05) {
            message.append(
                "Warning: More than " + std::to_string(portion * 100)
                + "% of your samples were removed! "
                  "You should check if your covariate file is correct\n");
        }
        reporter.report(message);
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
            std::string name = std::get<0>(valid_sample_index[cur_index]);
            m_sample_with_phenotypes[name] = cur_index;
            size_t original_index = std::get<1>(valid_sample_index[cur_index]);
            if (original_index != cur_index) {
                m_phenotype(static_cast<Eigen::Index>(cur_index), 0) =
                    m_phenotype(static_cast<Eigen::Index>(original_index), 0);
            }
        }

        m_phenotype.conservativeResize(
            static_cast<Eigen::Index>(valid_sample_index.size()), 1);
    }
    num_column = total_column;
}

void PRSice::gen_cov_matrix(const std::string& c_cov_file,
                            const std::vector<std::string> cov_header_name,
                            const std::vector<uint32_t> cov_header_index,
                            const std::vector<uint32_t> factor_cov_index,
                            const std::string& delim, Reporter& reporter)
{
    // The size of the map should be informative of the number of sample
    // currently included in the data
    size_t num_sample = m_sample_with_phenotypes.size();
    if (c_cov_file.empty()) {
        // if no covariates, just return a matrix of 1 with two column, one for
        // the intercept and the other for storing PRS. Both is 1 so that when
        // we need to calculate null, we can simply remove the first without
        // trouble
        m_independent_variables =
            Eigen::MatrixXd::Ones(static_cast<Eigen::Index>(num_sample), 2);
        return;
    }
    // obtain the index of each covariate
    // the key is the variable name and the value is the index on the matrix

    // need to account for the situation where the same variable name can
    // occur in different covariates
    // As the index are sorted, we can use vector

    // the index of the factor_list is the index of the covariate
    // the key of the nested unorder map is the factor and the value is the
    // factor level (similar to column index)
    std::vector<std::unordered_map<std::string, size_t>> factor_list;

    // an indexor to indicate whcih column should each covariate start from (as
    // there're factor covariates, the some covariates might take up more than
    // one column
    std::vector<size_t> cov_start_index;
    // by default the required number of column for the matrix is
    // intercept+PRS+number of covariate (when there're no factor input)
    size_t num_column = 2 + cov_header_index.size();
    // we will perform the first pass to the covariate file which will remove
    // all samples with missing covariate and will also generate the factor
    // level
    process_cov_file(c_cov_file, factor_cov_index, cov_start_index,
                     cov_header_index, cov_header_name, factor_list, num_column,
                     delim, reporter);
    std::string message = "Processing the covariate file: " + c_cov_file + "\n";
    message.append("==============================\n");
    reporter.report(message);
    // update the number of sample to account for missing covariates
    num_sample = m_sample_with_phenotypes.size();
    // initalize the matrix to the desired size
    m_independent_variables =
        Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(num_sample),
                              static_cast<Eigen::Index>(num_column));
    m_independent_variables.col(0).setOnes();
    m_independent_variables.col(1).setOnes();
    // now we only need to fill in the independent matrix without worry
    // about other stuff
    std::ifstream cov;
    cov.open(c_cov_file.c_str());
    if (!cov.is_open()) {
        std::string error_message =
            "Error: Cannot open covariate file: " + c_cov_file;
        throw std::runtime_error(error_message);
    }
    std::vector<std::string> token;
    std::string line, id;
    size_t max_index = cov_header_index.back() + 1, index, cur_index, f_level;
    uint32_t cur_factor_index = 0;
    size_t num_factor = factor_cov_index.size(),
           num_cov = cov_header_index.size();
    while (std::getline(cov, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() < max_index) {
            std::string error_message =
                "Error: Malformed covariate file, should contain at least "
                + std::to_string(max_index) + " column!";
            throw std::runtime_error(error_message);
        }
        id = (m_ignore_fid) ? token[0] : token[0] + delim + token[1];
        if (m_sample_with_phenotypes.find(id) != m_sample_with_phenotypes.end())
        {
            // Only valid samples will be found in the m_sample_with_phenotypes
            // map structure
            // reset the index to 0
            cur_factor_index = 0;
            // get the row number
            index = m_sample_with_phenotypes[id];
            for (size_t i_cov = 0; i_cov < num_cov; ++i_cov) {
                if (cur_factor_index >= num_factor
                    || cov_header_index[i_cov]
                           != factor_cov_index[cur_factor_index])
                {
                    // noraml covariate
                    // we don't need to deal with invalid conversion situation
                    // as those should be taken cared of by the process_cov_file
                    // function
                    m_independent_variables(
                        static_cast<Eigen::Index>(index),
                        static_cast<Eigen::Index>(cov_start_index[i_cov])) =
                        misc::convert<double>(token[cov_header_index[i_cov]]);
                }
                else
                {
                    // this is a factor
                    // and the level of the current factor is f_level
                    f_level = factor_list[cur_factor_index]
                                         [token[cov_header_index[i_cov]]];
                    if (f_level != 0) {
                        // if this is not the reference level, we will add 1 to
                        // the matrix
                        // we need to -1 as the reference level = 0 and the
                        // second level = 1 but for the second level, it should
                        // propagate the first column
                        cur_index = cov_start_index[i_cov] + f_level - 1;
                        m_independent_variables(
                            static_cast<Eigen::Index>(index),
                            static_cast<Eigen::Index>(cur_index)) = 1;
                    }
                    ++cur_factor_index;
                }
            }
        }
    }

    message = "After reading the covariate file, "
              + std::to_string(m_sample_with_phenotypes.size())
              + " sample(s) included in the analysis\n";
    reporter.report(message);
}

bool PRSice::run_prsice(const Commander& c_commander, const size_t pheno_index,
                        const size_t region_index,
                        const std::vector<size_t>& region_membership,
                        const std::vector<size_t>& region_start_idx,
                        Genotype& target)
{
    const bool no_regress = c_commander.no_regress();
    // only print out all scores if this is the first phenotype
    const bool print_all_scores = c_commander.all_scores() && pheno_index == 0;
    const size_t num_thread = c_commander.thread();
    const bool non_cumulate = c_commander.non_cumulate();
    const bool use_ref_maf = c_commander.use_ref_maf();
    const size_t num_samples_included = target.num_sample();

    std::vector<size_t>::const_iterator cur_start_idx =
        region_membership.cbegin();
    std::advance(cur_start_idx, region_start_idx[region_index]);
    std::vector<size_t>::const_iterator cur_end_idx =
        region_membership.cbegin();
    if (region_index + 1 >= region_start_idx.size()) {
        cur_end_idx = region_membership.cend();
    }
    else
    {
        std::advance(cur_end_idx, region_start_idx[region_index + 1]);
    }

    Eigen::initParallel();
    Eigen::setNbThreads(static_cast<int>(num_thread));
    m_best_index = -1;
    // m_num_snp_included will store the current number of SNP included. This is
    // the global number and the true number per sample might differ due to
    // missingness. This is only use for display
    m_num_snp_included = 0;
    // m_perm_result stores the result (T-value) from each permutation and is
    // then used for calculation of empirical p value
    m_perm_result.resize(m_num_perm, 0);
    m_best_sample_score.clear();
    m_prs_results.resize(target.num_threshold());
    // set to -1 to indicate not done
    for (auto&& p : m_prs_results) {
        p.threshold = -1;
        p.r2 = 0.0;
        p.num_snp = 0;
    }
    // if cur_start_idx == cur_end_idx, this is an empty region
    if (cur_start_idx == cur_end_idx) {
        return false;
    }
    // initialize score vector
    // this stores the best score for each sample. Ideally, we will only do it
    // once per phenotype as subsequent region should all have the same number
    // of samples
    if (region_index == 0) m_best_sample_score.resize(target.num_sample());

    // now prepare all score
    std::fstream all_out;
    if (print_all_scores && pheno_index == 0) {
        const std::string all_out_name = c_commander.out() + ".all.score";
        all_out.open(all_out_name.c_str(),
                     std::fstream::out | std::fstream::in | std::fstream::ate);
        if (!all_out.is_open()) {
            std::string error_message =
                "Cannot open file " + all_out_name + " for write";
            throw std::runtime_error(error_message);
        }
    }


    // current threshold iteration
    // must iterate after each threshold even if no-regress is called
    size_t prs_result_idx = 0;
    double cur_threshold = 0.0;
    // we want to know if user want to obtain the standardized PRS
    const bool require_standardize = (m_score == SCORING::STANDARDIZE);
    // print the progress bar
    print_progress();
    // indicate if this is the first run. If this is the first run, get_score
    // will perform assignment instead of addition
    bool first_run = true;
    // we will call the read score function from the target genotype, which will
    // then proceed to read in and calculate the PRS for the given category
    // (defined by the cur_index, which points to the first SNP of the p-value
    // threshold)


    while (target.get_score(cur_start_idx, cur_end_idx, cur_threshold,
                            m_num_snp_included, non_cumulate,
                            require_standardize, first_run, use_ref_maf))
    {
        m_analysis_done++;
        print_progress();
        if (print_all_scores) {
            for (size_t sample = 0; sample < num_samples_included; ++sample) {
                // we will calculate the the number of white space we need to
                // skip to reach the current sample + threshold's output
                // position
                int loc = m_all_file.header_length
                          + static_cast<int>(sample)
                                * (m_all_file.line_width + NEXT_LENGTH)
                          + NEXT_LENGTH + m_all_file.skip_column_length
                          + m_all_file.processed_threshold
                          + m_all_file.processed_threshold * m_numeric_width;
                all_out.seekp(loc);
                // then we will output the score
                all_out << std::setprecision(m_precision)
                        << target.calculate_score(m_score, sample);
            }
        }
        // we need to then tell the file that we have finish processing one
        // threshold. Next time we output another PRS, it should be output in
        // the column of the next threshold
        m_all_file.processed_threshold++;
        if (!no_regress) {
            // We only perform the regression analysis if we would like to
            // perform the regresswion
            regress_score(target, cur_threshold, num_thread, pheno_index,
                          prs_result_idx);

            if (m_perform_perm) {
                // perform permutation if required
                permutation(num_thread, m_target_binary[pheno_index]);
            }
        }
        prs_result_idx++;
        first_run = false;
    }

    // we need to process the permutation result if permutation is required
    if (m_perform_perm) process_permutations();
    if (!no_regress) {
        // if regression was performed, we will also generate the best score
        // output
        print_best(target, pheno_index, c_commander);
    }
    return true;
}

void PRSice::print_best(Genotype& target, const size_t pheno_index,
                        const Commander& commander)
{

    // read in the name of the phenotype. If there's only one phenotype name,
    // we'll do assign an empty string to phenotyp name
    std::string pheno_name = "";
    if (pheno_info.name.size() > 1) pheno_name = pheno_info.name[pheno_index];
    std::string output_prefix = commander.out();
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    // we generate one best score file per phenotype. The reason for this is to
    // not make the best score file too big, and because each phenotype might
    // have different set of sample included in the regression due to
    // missingness
    const std::string out_best = output_prefix + ".best";
    // we have to overwrite the white spaces with the desired values
    std::fstream best_out(out_best.c_str(), std::fstream::out | std::fstream::in
                                                | std::fstream::ate);
    auto&& best_info = m_prs_results[static_cast<size_t>(m_best_index)];
    size_t best_snp_size = best_info.num_snp;
    if (best_snp_size == 0) {
        fprintf(stderr, "Error: Best R2 obtained when no SNPs were included\n");
        fprintf(stderr, "       Cannot output the best PRS score\n");
    }
    else
    {
        for (int sample = 0; sample < static_cast<int>(target.num_sample());
             ++sample)
        {
            // samples that are extracted are ignored
            // sample excluded will not be output here
            /* std::string has_pheno =
                target.sample_in_regression(static_cast<size_t>(sample)) ? "Yes"
                                                                         : "No";
                                                                         */
            int loc = m_best_file.header_length
                      + sample * (m_best_file.line_width + NEXT_LENGTH)
                      + NEXT_LENGTH + m_best_file.skip_column_length
                      + m_best_file.processed_threshold
                      + m_best_file.processed_threshold * m_numeric_width;

            best_out.seekp(loc);
            best_out << std::setprecision(m_precision)
                     << m_best_sample_score[static_cast<size_t>(sample)];
        }
    }
    best_out.close();
    // once we finish outputing the result, we need to increment the
    // processed_threshold index such that when we process the next region, we
    // will be writing to the next column instead of overwriting the current
    // column
    m_best_file.processed_threshold++;
}

void PRSice::regress_score(Genotype& target, const double threshold,
                           size_t thread, const size_t pheno_index,
                           const size_t prs_result_idx)
{
    double r2 = 0.0, r2_adjust = 0.0, p_value = 0.0, coefficient = 0.0,
           se = 0.0;
    const Eigen::Index num_regress_samples =
        static_cast<Eigen::Index>(m_matrix_index.size());
    if (m_num_snp_included == 0
        || (m_num_snp_included == m_prs_results[prs_result_idx].num_snp))
    {
        // if we haven't read in any SNP, or that we have the same number of SNP
        // as the previous threshold, we will skip (normally this should not
        // happen, as we have removed any threshold that doesn't contain any
        // SNPs)
        return;
    }

    for (Eigen::Index sample_id = 0; sample_id < num_regress_samples;
         ++sample_id)
    {
        // we can directly read in the matrix index from m_matrix_index vector
        // and assign the PRS directly to the indep variable matrix
        m_independent_variables(sample_id, 1) = target.calculate_score(
            m_score, m_matrix_index[static_cast<size_t>(sample_id)]);
    }

    if (m_target_binary[pheno_index]) {
        // if this is a binary phenotype, we will perform the GLM model
        try
        {
            Regression::glm(m_phenotype, m_independent_variables, p_value, r2,
                            coefficient, se, 25, thread, true);
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
        Regression::linear_regression(m_phenotype, m_independent_variables,
                                      p_value, r2, r2_adjust, coefficient, se,
                                      thread, true);
    }

    // If this is the best r2, then we will add it
    int best_index = m_best_index;
    if (prs_result_idx == 0 || best_index < 0
        || m_prs_results[static_cast<size_t>(best_index)].r2 < r2)
    {
        m_best_index = static_cast<int>(prs_result_idx);
        size_t num_include_samples = target.num_sample();
        for (size_t s = 0; s < num_include_samples; ++s) {
            // we will have to store the best scores. we cannot directly copy
            // from the m_independent_variable as some samples which might have
            // excluded from the regression model but we still want their PRS.
            m_best_sample_score[s] = target.calculate_score(m_score, s);
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
    const double best_t = std::abs(m_prs_results[best_index].coefficient
                                   / m_prs_results[best_index].se);
    const auto num_better =
        std::count_if(m_perm_result.begin(), m_perm_result.end(),
                      [&best_t](double t) { return t > best_t; });
    m_prs_results[best_index].emp_p = (num_better + 1.0) / (m_num_perm + 1.0);
}

void PRSice::permutation(const size_t n_thread, bool is_binary)
{

    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm_matrix(
        m_phenotype.rows());
    Eigen::setNbThreads(static_cast<int>(n_thread));
    Eigen::Index rank = 0;
    // logit_perm can only be true if it is binary trait and user used the
    // --logit-perm flag
    // can always do the following if
    // 1. QT trait (!is_binary)
    // 2. Not require logit perm
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomposed;
    Eigen::VectorXd pre_se_calulated;
    bool run_glm = true;
    if (!is_binary || !m_logit_perm) {
        // if our trait isn't binary or if we don't need to perform logistic
        // regression in our permutation, we will first decompose the
        // independent variable once, therefore speed up the other processes
        decomposed.compute(m_independent_variables);
        rank = decomposed.rank();
        // we also pre-calculate the matrix required for SE calculation as this
        // can be reused
        Eigen::MatrixXd R = decomposed.matrixR()
                                .topLeftCorner(rank, rank)
                                .triangularView<Eigen::Upper>();
        pre_se_calulated = (R.transpose() * R).inverse().diagonal();
        run_glm = false;
    }
    if (n_thread == 1) {
        // we will run the single thread function to reduce overhead
        run_null_perm_no_thread(decomposed, rank, pre_se_calulated, run_glm);
    }
    else
    {
        // we will run teh thread queue where one thread is responsible for
        // generating the shuffled phenotype whereas other threads are
        // responsible for calculating the t statistics
        Thread_Queue<std::pair<Eigen::VectorXd, size_t>> set_perm_queue;
        // For multi-threading we use the producer consumer pattern where the
        // producer will keep random shuffle the phenotypes and the consumers
        // will calculate the t-values
        // All we need to provide to the producer is the number of consumers and
        // the permutation queue use for
        std::thread producer(&PRSice::gen_null_pheno, this,
                             std::ref(set_perm_queue), n_thread - 1);
        std::vector<std::thread> consume_store;
        // we have used one thread as the producer, therefore we need to reduce
        // the number of available thread by 1
        for (size_t i = 0; i < n_thread - 1; ++i) {
            consume_store.push_back(
                std::thread(&PRSice::consume_null_pheno, this,
                            std::ref(set_perm_queue), std::ref(decomposed),
                            rank, std::cref(pre_se_calulated), run_glm));
        }
        // wait for all the threads to complete their job
        producer.join();
        for (auto&& consume : consume_store) consume.join();
    }
}

void PRSice::run_null_perm_no_thread(
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& decomposed,
    const Eigen::Index rank, const Eigen::VectorXd& pre_se, const bool run_glm)
{

    // reset the seed for each new threshold such that we will always generate
    // the same phenotpe for each threhsold without us needing to regenerate the
    // PRS
    std::mt19937 rand_gen{m_seed};
    // we want to count the number of samples included in the analysis
    const Eigen::Index num_regress_sample = m_phenotype.rows();
    // we will always include the intercept. If there's a demand, we can also
    // try to remove this requirement
    const bool intercept = true;
    // we will copy the phenotype into a new vector (Maybe not necessary, but
    // better safe than sorry)
    Eigen::VectorXd perm_pheno = m_phenotype;
    // count the number of loop we've finished so far
    size_t processed = 0;
    // pre-initialize these parameters
    double coefficient, se, r2, obs_p;
    double obs_t = -1;

    if (run_glm) {
        // we we want to use the logistic regression
        while (processed < m_num_perm) {
            // reassign the phenotype matrix. This is to ensure single threading
            // will produce the same result as multithreading given the same
            // seed
            perm_pheno = m_phenotype;
            // random shuffle the phenotype matrix
            std::shuffle(perm_pheno.data(),
                         perm_pheno.data() + num_regress_sample, rand_gen);
            // update the progress bar
            m_analysis_done++;
            print_progress();
            // run the logistic regression on the permuted phenotype
            Regression::glm(perm_pheno, m_independent_variables, obs_p, r2,
                            coefficient, se, 25, 1, true);
            obs_t = std::abs(coefficient / se);
            // note that for us to calculate the p-value from logistic
            // regression, we take the square of obs_t, but abs should give us
            // similar result and has the added benefit of not needing to worry
            // about binary or QT when we process the permutation results. We
            // therefore square the obs_t for the comparison (the higher the t,
            // the most significant). We use T for it is slower to reach bound
            // than p-value
            //  obs_t *= obs_t;
            // we then store the best t-value in our permutation results
            m_perm_result[processed] =
                std::max(obs_t, m_perm_result[processed]);
            // we have finished the current analysis.
            processed++;
        }
    }
    else
    {
        Eigen::VectorXd beta;
        Eigen::VectorXd se;
        Eigen::Index rdf;
        double rss, resvar;
        int se_index;
        while (processed < m_num_perm) {
            // for quantitative trait, we can directly compute the results
            // without re-computing the decomposition
            perm_pheno = m_phenotype;
            std::shuffle(perm_pheno.data(),
                         perm_pheno.data() + num_regress_sample, rand_gen);
            m_analysis_done++;
            print_progress();
            // directly solve the current phenotype to obtain the required beta
            beta = decomposed.solve(perm_pheno);
            rdf = num_regress_sample - rank;
            rss = (m_independent_variables * beta - perm_pheno).squaredNorm();
            se_index = intercept;
            // the decomposition might have moved our column order, thus we need
            // to know which column contain the results for our PRS
            for (int ind = 0; ind < beta.rows(); ++ind) {
                if (decomposed.colsPermutation().indices()(ind) == intercept) {
                    se_index = ind;
                    break;
                }
            }
            resvar = rss / static_cast<double>(rdf);
            // calculate the SE and also the T-value
            se = (pre_se * resvar).array().sqrt();
            // we take the absolute of the T-value as we only concern about the
            // magnitude
            obs_t = std::abs(beta(intercept) / se(se_index));
            // for T-value, we need the maximum T, not the smallest
            m_perm_result[processed] =
                std::max(obs_t, m_perm_result[processed]);
            processed++;
        }
    }
}

void PRSice::gen_null_pheno(Thread_Queue<std::pair<Eigen::VectorXd, size_t>>& q,
                            size_t num_consumer)
{
    size_t processed = 0;
    // we need to reset the seed for each threshold so that the phenotype
    // generated should always be the same for each threshold. This help us
    // avoid needing to repeat reading in the PRS for each permutation
    std::mt19937 rand_gen{m_seed};
    Eigen::setNbThreads(1);
    const Eigen::Index num_regress_sample = m_phenotype.rows();
    // this matrix is use for storing the shuffle information
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm_matrix(
        m_phenotype.rows());

    while (processed < m_num_perm) {
        // I am not completely sure whether we pass by reference or pass by
        // value to the queue (or if we move this out and cause undefined
        // behaviour). To avoid the aforementioned problems, we reinitialize a
        // new vector for each permutation
        Eigen::VectorXd null_pheno = m_phenotype;
        // then we shuffle it
        std::shuffle(null_pheno.data(), null_pheno.data() + num_regress_sample,
                     rand_gen);
        // and the we will push it to the queue where the consumers will pick up
        // and work on it
        q.emplace(std::make_pair(null_pheno, processed), num_consumer);
        m_analysis_done++;
        print_progress();
        processed++;
    }
    // send termination signal to the consumers
    for (size_t i = 0; i < num_consumer; ++i) {
        q.push(std::pair<Eigen::VectorXd, size_t>(Eigen::VectorXd(1), 0),
               num_consumer);
    }
}

void PRSice::consume_null_pheno(
    Thread_Queue<std::pair<Eigen::VectorXd, size_t>>& q,
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& decomposed, int rank,
    const Eigen::VectorXd& pre_se, bool run_glm)
{
    const Eigen::Index n = m_phenotype.rows();
    const bool intercept = true;
    // to avoid false sharing, all consumer will first store their permutation
    // result in their own vector and only update the master vector at the end
    // of permutation
    // temp_store stores all the T-value
    std::vector<double> temp_store;
    // and temp_index indicate which permutation the T-value corresponse to
    // it is important we don't mix up the permutation order as we're supposed
    // to mimic re-running PRSice N times with different permutation
    std::vector<size_t> temp_index;
    Eigen::VectorXd beta, se;
    Eigen::Index rdf = n - rank;
    std::pair<Eigen::VectorXd, size_t> input;
    double coefficient, se_res, r2, obs_p, rss, resvar;
    double obs_t = -1;
    Eigen::Index se_index;
    while (true) {
        // as long as we have not received a termination signal, we will
        // continue our processing and should read from the queue
        q.pop(input);
        // the termination signal is represented by an vector with only 1 row
        if (std::get<0>(input).rows() == 1) break;
        if (run_glm) {
            // the first entry from the queue should be the permuted phenotype
            // and the second entry is the index. We will pass the phenotype for
            // GLM analysis if required
            Regression::glm(std::get<0>(input), m_independent_variables, obs_p,
                            r2, coefficient, se_res, 25, 1, true);
            obs_t = std::abs(coefficient / se_res);
            // although the obs_t^2 is used for calculation of p-value, abs
            // should give us similar result and allow us to use the same
            // comparison method to process the permutation results obs_t *=
            // obs_t;
        }
        else
        {
            beta = decomposed.solve(std::get<0>(input));
            // Eigen::MatrixXd fitted = m_independent_variables * beta;
            rss = (m_independent_variables * beta - std::get<0>(input))
                      .squaredNorm();
            se_index = intercept;
            // need to account for situation where the matrix is being shuffled
            // (PRS no longer at the second column)
            for (Eigen::Index ind = 0; ind < beta.rows(); ++ind) {
                if (decomposed.colsPermutation().indices()(ind) == intercept) {
                    se_index = ind;
                    break;
                }
            }
            resvar = rss / static_cast<double>(rdf);
            se = (pre_se * resvar).array().sqrt();
            obs_t = std::abs(beta(intercept) / se(se_index));
        }
        temp_store.push_back(obs_t);
        temp_index.push_back(std::get<1>(input));
    }
    // once we received the termination signal, we can start propagating the
    // master vector with out content
    std::lock_guard<std::mutex> lock(lock_guard);
    for (size_t i = 0; i < temp_store.size(); ++i) {
        double obs_t = temp_store[i];
        auto&& index = temp_index[i];
        // if the t-value in the master vector is lower than our observed t,
        // update it
        if (m_perm_result[index] < obs_t) {
            m_perm_result[index] = obs_t;
        }
    }
}

void PRSice::prep_output(const std::string& out, const bool all_score,
                         const bool has_prev, const Genotype& target,
                         const std::vector<std::string>& region_name,
                         const size_t pheno_index)
{
    // As R has a default precision of 7, we will go a bit
    // higher to ensure we use up all precision
    std::string pheno_name = "";
    if (pheno_info.name.size() > 1)
        pheno_name = pheno_info.name[static_cast<size_t>(pheno_index)];
    std::string output_prefix = out;
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    const std::string output_name = output_prefix;
    const std::string out_prsice = output_name + ".prsice";
    // all score will be the same for all phenotype anyway, so we will only
    // generate it once
    const std::string out_all = out + ".all.score";
    const std::string out_best = output_name + ".best";
    std::ofstream prsice_out, best_out, all_out;

    // .prsice output
    // we only need to generate the header for it
    prsice_out.open(out_prsice.c_str());
    if (!prsice_out.is_open()) {
        std::string error_message =
            "Error: Cannot open file: " + out_prsice + " to write";
        throw std::runtime_error(error_message);
    }
    // we won't store the empirical p and competitive p output in the prsice
    // file now as that seems like a waste (only one threshold will contain that
    // information, storing that in the summary file should be enough)
    prsice_out << "Set\tThreshold\tR2\t";
    // but generate the adjusted R2 if prevalence is provided
    if (has_prev) prsice_out << "R2.adj\t";
    prsice_out << "P\tCoefficient\tStandard.Error\tNum_SNP\n";
    prsice_out.close();

    // .best output
    best_out.open(out_best.c_str());
    if (!best_out.is_open()) {
        std::string error_message =
            "Error: Cannot open file: " + out_best + " to write";
        throw std::runtime_error(error_message);
    }
    std::string header_line = "FID IID In_Regression";
    // The default name of the output should be PRS, but if we are running
    // PRSet, it should be call Base
    if (!m_perform_prset)
        header_line.append(" PRS");
    else
    {
        for (size_t i = 0; i < region_name.size(); ++i) {
            // the second item is always the background, which we will not
            // output
            if (i == 1) continue;
            header_line.append(" " + region_name[i]);
        }
    }
    // the safetest way to calculate the length we need to speed is to directly
    // count the number of byte involved
    auto begin_byte = best_out.tellp();
    best_out << header_line << "\n";
    auto end_byte = best_out.tellp();
    // we now know the exact number of byte the header contain and can correctly
    // skip it acordingly
    m_best_file.header_length = static_cast<int>(end_byte - begin_byte);
    // we will set the processed_threshold information to 0
    m_best_file.processed_threshold = 0;

    // each numeric output took 12 spaces, then for each output, there is one
    // space next to each
    m_best_file.line_width =
        m_max_fid_length /* FID */ + 1 /* space */ + m_max_iid_length /* IID */
        + 1 /* space */ + 3 /* Yes/No */ + 1   /* space */
        + static_cast<int>(region_name.size()) /* each region */
              * (m_numeric_width + 1 /* space */)
        + 1 /* new line */;

    m_best_file.skip_column_length =
        m_max_fid_length + 1 + m_max_iid_length + 1 + 3 + 1;


    // also handle all score here
    // but we will only try and generate the all score file when we are dealing
    // with the first phenotype (pheno_index == 0)
    const bool all_scores = all_score && !pheno_index;
    const bool print_background = (region_name.size() != 2);
    if (all_scores) {
        all_out.open(out_all.c_str());
        if (!all_out.is_open()) {
            std::string error_message =
                "Cannot open file " + out_all + " for write";
            throw std::runtime_error(error_message);
        }
        // we need to know the number of available thresholds so that we can
        // know how many white spaces we need to pad in
        std::vector<double> avail_thresholds = target.get_thresholds();
        // we want the threshold to be in sorted order as we will process the
        // SNPs from the smaller threshold to the highest (therefore, the
        // processed_threshold index should be correct)
        std::sort(avail_thresholds.begin(), avail_thresholds.end());
        int num_thresholds = static_cast<int>(avail_thresholds.size());
        begin_byte = all_out.tellp();
        all_out << "FID IID";
        // size_t header_length = 3+1+3;
        if (!m_perform_prset) {
            for (auto& thres : avail_thresholds) {
                all_out << " " << thres;
                // if we are not performing PRSet, it is easy, just one
                // column per threshold
            }
        }
        else
        {
            // but if we are performing PRSet, we will need to have region
            // number * threshold number thresholds
            for (size_t i = 0; i < region_name.size() - 1; ++i) {
                for (auto& thres : avail_thresholds) {
                    all_out << " " << region_name[i] << "_" << thres;
                }
            }
        }
        all_out << "\n";
        end_byte = all_out.tellp();
        // if the line is too long, we might encounter overflow
        m_all_file.header_length = static_cast<int>(end_byte - begin_byte);
        m_all_file.processed_threshold = 0;
        m_all_file.line_width =
            m_max_fid_length + 1 + m_max_iid_length + 1
            + num_thresholds
                  * static_cast<int>(region_name.size() - !print_background)
                  * (m_numeric_width + 1)
            + 1;
        m_all_file.skip_column_length = m_max_fid_length + m_max_iid_length + 2;
        // all_out << header_line << "\n";
    }

    // output sample IDs
    size_t num_samples_included = target.num_sample();
    std::string best_line;
    std::string name;
    for (size_t i_sample = 0; i_sample < num_samples_included; ++i_sample) {
        name = target.fid(i_sample) + " " + target.iid(i_sample);
        // when we print the best file, we want to also print whether the sample
        // is used in regression or not (so that user can easily reproduce their
        // results)
        best_line = name + " "
                    + ((target.sample_in_regression(i_sample)) ? "Yes" : "No");
        // we print a line containing m_best_file.line_width white space
        // characters, which we can then overwrite later on, therefore achieving
        // a vertical output
        best_out << std::setfill(' ') << std::setw(m_best_file.line_width)
                 << std::left << best_line << "\n";
        if (all_scores) {
            all_out << std::setfill(' ') << std::setw(m_all_file.line_width)
                    << std::left << name << "\n";
        }
    }

    // another one spacing for new line (just to be safe)
    m_all_file.line_width++;
    m_best_file.line_width++;
    // don't need to close the files as they will automatically be closed when
    // we move out of the function
}

void PRSice::output(const Commander& c_commander,
                    const std::vector<std::string>& region_names,
                    const size_t pheno_index, const size_t region_index)
{
    // if prevalence is provided, we'd like to generate calculate the adjusted
    // R2
    std::vector<double> prev = c_commander.prevalence();
    bool has_prevalence = c_commander.has_prevalence();
    const bool is_binary =
        c_commander.is_binary(static_cast<size_t>(pheno_index));
    double top = 1.0, bottom = 1.0, prevalence = -1;
    if (has_prevalence && is_binary) {
        size_t num_binary = 0;
        for (size_t i = 0; i < static_cast<size_t>(pheno_index); ++i) {
            if (c_commander.is_binary(i))
                num_binary++; // this is the number of previous binary traits
        }
        int num_case = static_cast<int>(m_phenotype.sum());
        double case_ratio = static_cast<double>(num_case)
                            / static_cast<double>(m_phenotype.rows());
        prevalence = prev[num_binary];
        // the following is from Lee et al A better coefficient paper
        double x = misc::qnorm(1 - prevalence);
        double z = misc::dnorm(x);
        double i2 = z / prevalence;
        double cc = prevalence * (1 - prevalence) * prevalence
                    * (1 - prevalence)
                    / (z * z * case_ratio * (1 - case_ratio));
        double theta =
            i2 * ((case_ratio - prevalence) / (1 - prevalence))
            * (i2 * ((case_ratio - prevalence) / (1 - prevalence)) - x);
        double e = 1
                   - pow(case_ratio, (2 * case_ratio))
                         * pow((1 - case_ratio), (2 * (1 - case_ratio)));
        top = cc * e;
        bottom = cc * e * theta;
    }

    const std::string pheno_name =
        (pheno_info.name.size() > 1)
            ? pheno_info.name[static_cast<std::vector<std::string>::size_type>(
                  pheno_index)]
            : "";
    std::string output_prefix = c_commander.out();
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);

    // check if this is a valid phenotyep
    if (m_best_index == -1) {
        // when m_best_index == -1, we don't have any valid PRS output
        fprintf(stderr, "Error: No valid PRS ");
        if (m_perform_prset)
            fprintf(stderr, "for %s", region_names[region_index].c_str());
        fprintf(stderr, "!\n");
        return;
    }
    // now we know can generate the prsice file
    std::string out_prsice = output_prefix + ".prsice";
    std::ofstream prsice_out;
    prsice_out.open(out_prsice.c_str(), std::fstream::app);
    if (!prsice_out.is_open()) {
        std::string error_message =
            "Error: Cannot open file: " + out_prsice + " to write";
        throw std::runtime_error(error_message);
    }
    // go through every result and output
    for (size_t i = 0; i < m_prs_results.size(); ++i) {
        if (m_prs_results[i].threshold < 0 || m_prs_results[i].p < 0) continue;
        double full = m_prs_results[i].r2;
        double null = m_null_r2;
        double full_adj = full;
        double null_adj = null;
        if (has_prevalence) {
            full_adj = top * full / (1 + bottom * full);
            null_adj = top * null / (1 + bottom * null);
        }

        double r2 = full - null;
        prsice_out << region_names[region_index] << "\t"
                   << m_prs_results[i].threshold << "\t" << r2 << "\t";
        if (has_prevalence) {
            if (is_binary)
                prsice_out << full_adj - null_adj << "\t";
            else
                prsice_out << "NA\t";
        }
        prsice_out << m_prs_results[i].p << "\t" << m_prs_results[i].coefficient
                   << "\t" << m_prs_results[i].se << "\t"
                   << m_prs_results[i].num_snp << "\n";
        // the empirical p-value will now be excluded from the .prsice output
        // (the "-" isn't that helpful anyway)
    }
    prsice_out.close();
    auto&& best_info =
        m_prs_results[static_cast<std::vector<prsice_result>::size_type>(
            m_best_index)];


    // we will extract the information of the best threshold, store it and use
    // it to generate the summary file
    // in theory though, I should be able to start generating the summary file
    prsice_summary prs_sum;
    prs_sum.pheno = pheno_name;
    prs_sum.set = region_names[region_index];
    prs_sum.result = best_info;
    prs_sum.r2_null = m_null_r2;
    prs_sum.top = top;
    prs_sum.bottom = bottom;
    prs_sum.prevalence = prevalence;
    // we don't run competitive testing on the base region
    // therefore we skip region_index == 0 (base is always
    // the first region)
    prs_sum.has_competitive = (region_index == 0);
    m_prs_summary.push_back(prs_sum);
    if (best_info.p > 0.1)
        m_significant_store[0]++;
    else if (best_info.p > 1e-5)
        m_significant_store[1]++;
    else
        m_significant_store[2]++;
}

void PRSice::summarize(const Commander& commander, Reporter& reporter)
{
    // we need to know if we are going to write "and" in the output, thus need a
    // flag to indicate if there are any previous outputs

    bool has_previous_output = false;
    // we will output a short summary file
    std::string message = "There are ";
    if (m_significant_store[0] != 0) {
        message.append(
            misc::to_string(m_significant_store[0])
            + " region(s)/phenotype(s) with p-value > 0.1 (\033[1;31mnot "
              "significant\033[0m);");
        has_previous_output = true;
    }
    if (m_significant_store[1] != 0) {
        if (m_significant_store[2] == 0 && has_previous_output) {
            message.append(" and ");
        }
        message.append(
            misc::to_string(m_significant_store[1])
            + " region(s) with p-value between "
              "0.1 and 1e-5 (\033[1;31mmay not be significant\033[0m);");
        has_previous_output = true;
    }
    if (m_significant_store[2] != 0) {
        if (has_previous_output) message.append(" and ");
        message.append(std::to_string(m_significant_store[2])
                       + " region(s) with p-value less than 1e-5.");
    }
    if (!has_previous_output) {
        message.append(
            " Please note that these results are inflated due to the "
            "overfitting inherent in finding the best-fit "
            "PRS (but it's still best to find the best-fit PRS!).\n"
            "You can use the --perm option (see manual) to calculate "
            "an empirical P-value.");
    }
    reporter.report(message);
    // now we generate the output file
    std::string out_name = commander.out() + ".summary";
    std::ofstream out;
    out.open(out_name.c_str());
    if (!out.is_open()) {
        std::string error_message =
            "Error: Cannot open file: " + out_name + " to write";
        throw std::runtime_error(error_message);
    }
    const bool has_prevalence = commander.has_prevalence();
    out << "Phenotype\tSet\tThreshold\tPRS.R2";
    if (has_prevalence) {
        // if we have the prevalence adjustment, we would also like to output
        // the adjusted R2 together with the unadjusted (just in case)
        out << "\tPRS.R2.adj";
    }
    out << "\tFull.R2\tNull."
           "R2\tPrevalence\tCoefficient\tStandard.Error\tP\tNum_SNP";
    if (m_perform_competitive) out << "\tCompetitive.P";
    if (m_perform_perm) out << "\tEmpirical-P";
    out << "\n";
    for (auto&& sum : m_prs_summary) {
        out << ((sum.pheno.empty()) ? "-" : sum.pheno) << "\t" << sum.set
            << "\t" << sum.result.threshold;
        // by default, phenotype that doesn't have the prevalence information
        // will have a prevalence of -1
        if (sum.prevalence > 0) {
            // calculate the adjusted R2 for binary traits
            double full = sum.result.r2;
            double null = sum.r2_null;
            full = sum.top * full / (1 + sum.bottom * full);
            null = sum.top * null / (1 + sum.bottom * null);
            out << "\t" << sum.result.r2 - sum.r2_null << "\t" << full - null
                << "\t" << full << "\t" << null << "\t" << sum.prevalence;
        }
        else if (has_prevalence)
        {
            // and replace the R2 adjust by NA if the sample doesn't have
            // prevalence (i.e. quantitative trait)
            out << "\t" << sum.result.r2 - sum.r2_null << "\tNA\t"
                << sum.result.r2 << "\t" << sum.r2_null << "\t"
                << sum.prevalence;
        }
        else
        {
            // if the prevalence is never provided, we don't need to calculate
            // the adjusted R2
            out << "\t" << sum.result.r2 - sum.r2_null << "\t" << sum.result.r2
                << "\t" << sum.r2_null << "\t-";
        }
        // now generate the rest of the output
        out << "\t" << sum.result.coefficient << "\t" << sum.result.se << "\t"
            << sum.result.p << "\t" << sum.result.num_snp;
        // As we never run competitive analysis on the base data set, we need to
        // account for that (default will have a p-value less than 0)
        if (m_perform_competitive && (sum.result.competitive_p >= 0.0)) {
            out << "\t" << sum.result.competitive_p;
        }
        else if (m_perform_competitive)
        {
            // this is the base. While it look nicer to have - to represent not
            // available, it become nightmarish for R and other downstream
            // analysis where we need to as.numeric(as.character()). Therefore
            // we now use NA instead
            out << "\tNA";
        }
        if (m_perform_perm) out << "\t" << sum.result.emp_p;
        out << "\n";
    }
    out.close();
}

PRSice::~PRSice()
{
    // dtor
}

void PRSice::null_set_no_thread(
    Genotype& target, const std::vector<size_t>::const_iterator& bk_start_idx,
    const std::vector<size_t>::const_iterator& bk_end_idx,
    const std::map<size_t, std::vector<size_t>>& set_index,
    std::vector<double>& obs_t_value, std::vector<size_t>& set_perm_res,
    const size_t num_perm, const bool is_binary, const bool require_standardize,
    const bool use_ref_maf)
{
    // we need to know the size of the largest gene set (excluding the base and
    // background)
    // it's a map, last element should be the largest
    const size_t max_size = set_index.rbegin()->first;
    // need to count the number of permutation done
    size_t processed = 0;
    const Eigen::Index num_sample =
        static_cast<Eigen::Index>(m_matrix_index.size());
    double coefficient, se, r2, r2_adjust, obs_p, t_value;
    std::mt19937 g(m_seed);
    // we need to know how many background SNPs are there
    const size_t num_background =
        static_cast<size_t>(std::distance(bk_start_idx, bk_end_idx));
    std::vector<size_t> background(bk_start_idx, bk_end_idx);
    // a boolean to tell the genotype class whether the PRS should be reset
    bool first_run = true;
    while (processed < num_perm) {
        size_t begin = 0;
        // we will shuffle n where n is the set with the largest size
        // this is the Fisher-Yates shuffle algorithm for random selection
        // without replacement
        size_t num_snp = max_size;
        while (num_snp--) {
            std::uniform_int_distribution<int> dist(
                static_cast<int>(begin), static_cast<int>(num_background) - 1);
            int advance_index = dist(g);
            std::swap(background[static_cast<size_t>(begin)],
                      background[static_cast<size_t>(advance_index)]);
            ++begin;
        }
        //  we have now selected N SNPs from the background. We can then
        //  construct the PRS based on these index
        first_run = true;
        size_t prev_size = 0;
        for (auto&& set_size : set_index) {
            // now we iterate through each set size
            // in theory this will reduce our I/O. If the set sizes
            // are 10, 100 and 1000, then the number of SNPs we read
            // per set will be 10, 90, 900. Which will help to reduce
            // a lot of reading if the set sizes are very similar

            // read in genotype here
            target.get_null_score(set_size.first, prev_size, background,
                                  first_run, require_standardize, use_ref_maf);

            // now that we've constructed the PRS, for any subsequent set sizes,
            // we should not reset the PRS calculation until the next
            // permutation, thus we need to set first_run to true
            first_run = false;
            prev_size = set_size.first;
            // we have now read in the PRS and should now assign it to the
            // m_independent_variables matrix
            for (Eigen::Index sample_id = 0; sample_id < num_sample;
                 ++sample_id)
            {
                m_independent_variables(sample_id, 1) = target.calculate_score(
                    m_score, m_matrix_index[static_cast<size_t>(sample_id)]);
            }
            m_analysis_done++;
            print_progress();
            //  we can now perform the glm or linear regression analysis
            if (is_binary) {
                Regression::glm(m_phenotype, m_independent_variables, obs_p, r2,
                                coefficient, se, 25, 1, true);
                t_value = std::abs(coefficient / se);
                // in GLM, p_value is calculated by t_value^2, so to mimic that,
                // we also do t_value^2, though this should give the same result
                // as abs(t_value) so we will just use abs(t-value) for both QT
                // and binary traits t_value *= t_value;
            }
            else
            {
                Regression::linear_regression(
                    m_phenotype, m_independent_variables, obs_p, r2, r2_adjust,
                    coefficient, se, 1, true);
                t_value = std::abs(coefficient / se);
            }
            // set_size second contain the indexs to each set with this size
            for (auto&& set_index : set_size.second) {
                set_perm_res[set_index] += (obs_t_value[set_index] < t_value);
            }
        }
        processed++;
    }
}

void PRSice::produce_null_prs(
    Thread_Queue<std::pair<std::vector<double>, size_t>>& q, Genotype& target,
    const std::vector<size_t>::const_iterator& bk_start_idx,
    const std::vector<size_t>::const_iterator& bk_end_idx, size_t num_consumer,
    std::map<size_t, std::vector<size_t>>& set_index, const size_t num_perm,
    const bool require_standardize, const bool use_ref_maf)
{
    // we need to know the size of the biggest set
    const size_t max_size = set_index.rbegin()->first;
    const size_t num_sample = m_matrix_index.size();
    const size_t num_regress_sample =
        static_cast<size_t>(m_independent_variables.rows());
    const size_t num_background =
        static_cast<size_t>(std::distance(bk_start_idx, bk_end_idx));
    size_t processed = 0;
    size_t prev_size = 0;
    size_t r;
    // we seed the random number generator
    std::mt19937 g(m_seed);
    // get the SNP index for background
    std::vector<size_t> background(bk_start_idx, bk_end_idx);
    bool first_run = true;
    std::vector<size_t>::size_type advance_index, begin;
    while (processed < num_perm) {
        // here we perform random sampling without replacement using the
        // Fisher-Yates shuffle algorithm
        begin = 0;
        size_t num_snp = max_size;
        while (num_snp--) {
            std::uniform_int_distribution<int> dist(
                static_cast<int>(begin), static_cast<int>(num_background) - 1);
            r = background[begin];
            advance_index = static_cast<size_t>(dist(g));
            background[begin] = background[advance_index];
            background[advance_index] = r;
            ++begin;
        }
        first_run = true;
        prev_size = 0;
        for (auto&& set_size : set_index) {
            // for each gene sets size, we calculate the PRS
            target.get_null_score(set_size.first, prev_size, background,
                                  first_run, require_standardize, use_ref_maf);
            first_run = false;
            // we need to know how many SNPs we have already read, such that we
            // can skip reading this number of SNPs for the next set
            prev_size = set_size.first;
            // we store the PRS in a new vector to avoid crazy error with move
            // semetics and stuff which I have not fully understand
            std::vector<double> prs(num_regress_sample, 0);
            for (size_t sample_id = 0; sample_id < num_sample; ++sample_id) {
                // propagate the prs vector
                prs[sample_id] =
                    target.calculate_score(m_score, m_matrix_index[sample_id]);
            }
            // then we push the result prs to the queue, which can then picked
            // up by the consumers
            q.emplace(std::make_pair(prs, set_size.first), num_consumer);
            m_analysis_done++;
            print_progress();
        }
        processed++;
    }
    // send termination signal to the consumers
    for (size_t i = 0; i < num_consumer; ++i) {
        // termination signal is represented by an empty vector
        q.emplace(std::make_pair(std::vector<double>(), 0), num_consumer);
    }
}


void PRSice::consume_prs(
    Thread_Queue<std::pair<std::vector<double>, size_t>>& q,
    std::map<size_t, std::vector<size_t>>& set_index,
    std::vector<double>& obs_t_value, std::vector<size_t>& set_perm_res,
    const bool is_binary)
{
    // we first make a local copy of the independent matrix to ensure thread
    // safety
    Eigen::MatrixXd independent = m_independent_variables;
    const Eigen::Index num_regress_sample =
        static_cast<Eigen::Index>(m_matrix_index.size());
    // to avoid false sharing and frequent lock, we wil first store all
    // permutation results within a temporary vector
    std::vector<uint32_t> temp_perm_res(set_perm_res.size(), 0);
    double coefficient, se, r2, r2_adjust;
    double obs_p = 2.0; // for safety reason, make sure it is out bound
    // results from queue will be stored in the prs_info
    std::pair<std::vector<double>, size_t> prs_info;
    // now listen for producer
    while (true) {
        q.pop(prs_info);
        if (std::get<0>(prs_info).empty()) {
            // all job finished as the termination signal = empty vector
            break;
        }
        // update the independent variable matrix with the new PRS
        for (Eigen::Index i_sample = 0; i_sample < num_regress_sample;
             ++i_sample)
        {
            independent(i_sample, 1) =
                std::get<0>(prs_info)[static_cast<size_t>(i_sample)];
        }
        // then perform regression analysis to obtain the t-value
        if (is_binary) {
            Regression::glm(m_phenotype, independent, obs_p, r2, coefficient,
                            se, 25, 1, true);
        }
        else
        {
            Regression::linear_regression(m_phenotype, independent, obs_p, r2,
                                          r2_adjust, coefficient, se, 1, true);
        }
        double t_value = std::abs(coefficient / se);
        auto&& index = set_index[std::get<1>(prs_info)];
        // we register the number of time a more significant / bigger t-value is
        // obtained when compared to the observed t-value
        for (auto&& ref : index) {
            temp_perm_res[ref] += (obs_t_value[ref] < t_value);
        }
    }

    {
        // keep mutex lock within this scope
        std::unique_lock<std::mutex> locker(m_thread_mutex);
        size_t num_sets = temp_perm_res.size();
        // once everything is done, we go through the master copy of the
        // set_perm_res and add up the results
        for (size_t i = 0; i < num_sets; ++i) {
            set_perm_res[i] += temp_perm_res[i];
        }
    }
}

void PRSice::run_competitive(
    Genotype& target, const std::vector<size_t>::const_iterator& bk_start_idx,
    const std::vector<size_t>::const_iterator& bk_end_idx,
    const Commander& commander, const size_t pheno_index, Reporter& reporter)
{
    m_perform_competitive = true;
    fprintf(stderr, "\nStart competitive permutation\n");
    size_t num_perm;
    if (!commander.set_perm(num_perm)) {
        // false when we don't want to perform the competitive analysis
        return;
    }
    const bool require_standardize = (m_score == SCORING::STANDARDIZE);
    const bool is_binary = m_target_binary[pheno_index];
    const bool use_ref_maf = commander.use_ref_maf();
    const size_t num_bk_snps =
        static_cast<size_t>(std::distance(bk_start_idx, bk_end_idx));
    // the number of items to skip from the front of prs_summary
    size_t pheno_start_idx = 0;
    // obs_t_value stores the observed t-value
    std::vector<double> obs_t_value;
    // set_perm_res stores number of perm where a more sig result is obtained
    std::vector<size_t> set_perm_res;
    // set_index stores the index of sets with "key" size
    std::map<size_t, std::vector<size_t>> set_index;
    const size_t num_prs_res = m_prs_summary.size();
    bool started = false;
    // start at 1 to avoid the base set
    size_t cur_set_index = 0;
    size_t max_set_size = 0;
    for (size_t i = 0; i < num_prs_res; ++i) {
        // if we have already calculated the competitive p-value for the set, we
        // will just skip them. This help us to handle multiple-phenotype
        // without too much additional coding
        if (m_prs_summary[i].has_competitive || m_prs_summary[i].set == "Base")
            continue;
        if (!started) {
            // remembering the index of the first set that need to perform the
            // competitive p-value calculation. This allow us to later reassign
            // results to the sets
            pheno_start_idx = i;
            started = true;
        }
        auto&& res = m_prs_summary[i].result;
        // store the location of this set w.r.t number of SNPs in set
        set_index[res.num_snp].push_back(cur_set_index++);
        if (res.num_snp > max_set_size) max_set_size = res.num_snp;
        // ori_t_value will contain the obesrved t-value
        obs_t_value.push_back(std::abs(res.coefficient / res.se));
        set_perm_res.push_back(0);
    }
    if (max_set_size > num_bk_snps) {
        std::string error_messgae =
            "Error: Insufficient background SNPs for "
            "competitive analysis. Please ensure you have "
            "use the correct background. Will now generate skip "
            "the competitive analysis\n";
        for (size_t i = pheno_start_idx; i < num_prs_res; ++i) {
            // set them to true so that we will skip them for the next
            // phenotype (though in reality, they will all encounter the
            // same error. Need a better structure here)
            m_prs_summary[i].has_competitive = true;
        }
        reporter.report(error_messgae);
        return;
    }
    // now we can run the competitive testing
    // know how many thread we are allowed to use
    size_t num_thread = commander.thread();
    // find out the number of samples involved so that we can estimate the
    // required memory
    // The reason we need to go through the next set of calculation is that we
    // want to be certain that we've enough memory for the competitive analysis,
    // which will generate one new indepenendet variable matrix for each thread.
    // To avoid crashing due to insufficient memory, we will need to limit the
    // number of thread used
    const size_t mb = 1048576;
    const size_t num_regress_sample =
        static_cast<size_t>(m_independent_variables.rows());
    // the required memory is roughly number of Sample * number of covaraite.
    // But then for GLM, we might do iterative reweighting which also generate a
    // bunch of intermediate that requires memory. As a result of that, we
    // provide an over-estimation of required memory here to ensure we have
    // sufficient memory for our analysis
    const size_t basic_memory_required_per_thread =
        num_regress_sample * sizeof(double)
        * (static_cast<size_t>(m_independent_variables.cols()) * 6 + 15);
    // we need to calculate the total amount of memory available right now
    const size_t total_memory = misc::total_ram_available();
    // and the calculate the amount of memory we are allowed to use
    const size_t valid_memory = commander.max_memory(total_memory);
    // then calculate the amount of memory we've already used
    const size_t used_memory = misc::current_ram_usage();
    if (valid_memory <= used_memory) {
        // if we have used up all memory, we will exit
        std::string error_message =
            "Error: Not enough memory left for permutation. "
            "User allowed "
            + std::to_string(valid_memory / mb)
            + " Mb of memory but already used " + std::to_string(used_memory)
            + " Mb";
        fprintf(stderr, "\n");
        throw std::runtime_error(error_message);
    }
    // artificially reduce available memory to avoid memory overflow
    // ideally, if whole PRSice is within the same memory pool, we will
    // not need to do this reduction
    const size_t available_memory =
        static_cast<size_t>((valid_memory - used_memory) * 0.5);

    if (available_memory < basic_memory_required_per_thread) {
        std::string error_message =
            "Error: Not enough memory left for permutation. "
            "System allowed "
            + std::to_string(available_memory / mb)
            + " Mb of memory but required at least "
            + std::to_string(basic_memory_required_per_thread) + " Mb";
        fprintf(stderr, "\n");
        throw std::runtime_error(error_message);
    }
    // reduce number of threads to account for memory available
    if (available_memory / basic_memory_required_per_thread < num_thread) {
        num_thread = available_memory / basic_memory_required_per_thread;
    }
    // now, we should be safe to run the competitive p-value analysis wiht
    // num_thread without worry about insufficient memory

    if (num_thread > 1) {
        //  similar to permutation for empirical p-value calculation, we employ
        //  the producer consumer pattern where one thread is responsible for
        //  reading in the PRS and construct the required independent variable
        //  and other threads are responsible for the calculation
        Thread_Queue<std::pair<std::vector<double>, size_t>> set_perm_queue;
        std::thread producer(&PRSice::produce_null_prs, this,
                             std::ref(set_perm_queue), std::ref(target),
                             std::cref(bk_start_idx), std::cref(bk_end_idx),
                             num_thread - 1, std::ref(set_index), num_perm,
                             require_standardize, use_ref_maf);
        std::vector<std::thread> consumer_store;
        for (size_t i_thread = 0; i_thread < num_thread - 1; ++i_thread) {
            consumer_store.push_back(std::thread(
                &PRSice::consume_prs, this, std::ref(set_perm_queue),
                std::ref(set_index), std::ref(obs_t_value),
                std::ref(set_perm_res), is_binary));
        }

        producer.join();
        for (auto&& thread : consumer_store) thread.join();
    }
    else
    {
        // alternatively, if we only got one thread, we will use the no thread
        // function to reduce threading overhead
        null_set_no_thread(target, bk_start_idx, bk_end_idx, set_index,
                           obs_t_value, set_perm_res, num_perm, is_binary,
                           require_standardize, use_ref_maf);
    }
    // start_index is the index of m_prs_summary[i], not the actual index
    // on set_perm_res.
    // this will iterate all sets from beginning of current phenotype
    // to the last set within the current phenotype
    // because of the sequence of how set_perm_res is contructed,
    // the results for each set should be sequentially presented in
    // set_perm_res. Index for set_perm_res results are therefore
    // i - start_index
    for (size_t i = pheno_start_idx; i < num_prs_res; ++i) {
        auto&& res = m_prs_summary[i].result;
        // we need to minus out the start index from i such that our index start
        // at 0, which is the assumption of set_perm_res
        res.competitive_p =
            (static_cast<double>(set_perm_res[(i - pheno_start_idx)]) + 1.0)
            / (static_cast<double>(num_perm) + 1.0);
        m_prs_summary[i].has_competitive = true;
    }
}
