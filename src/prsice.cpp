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


// ARG!! THIS IS UGLY!!
// But I don't want to pass a void package to the function containing all this
// that'd be a nigtmare to backware convert them back...
bool PRSice::g_logit_perm = false;
Eigen::MatrixXd PRSice::g_independent_variables;
std::vector<double> PRSice::g_perm_result;
std::unordered_map<uintptr_t, PRSice::perm_info> PRSice::g_perm_range;
Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PRSice::g_perm_pre_decomposed;
std::vector<double> PRSice::g_permuted_pheno;
Eigen::VectorXd PRSice::g_pre_se_calulated;

void PRSice::pheno_check(const Commander& c_commander, Reporter& reporter)
{
    std::vector<std::string> pheno_header = c_commander.pheno_col();
    std::string pheno_file = c_commander.pheno_file();
    std::string message = "";
    if (pheno_file.empty()) {
        pheno_info.use_pheno = false;
        pheno_info.binary.push_back(c_commander.is_binary(0));
    }
    else
    {
        std::ifstream pheno;
        pheno.open(pheno_file.c_str());
        if (!pheno.is_open()) {
            std::string error_message =
                "Cannot open phenotype file: " + pheno_file;
            throw std::runtime_error(error_message);
        }
        std::string line;
        std::getline(pheno, line); // assume header line
        if (line.empty()) {
            throw std::runtime_error(
                "Cannot have empty header line for phenotype file!");
        }
        pheno.close();
        misc::trim(line);
        std::vector<std::string> col = misc::split(line);
        if (col.size() < (size_t)(1 + !m_ignore_fid)) {
            throw std::runtime_error(
                "Error: Not enough column in Phenotype file."
                "Have you use the --ignore-fid option");
        }
        std::string sample_id = col[0];
        if (!m_ignore_fid && col.size() > 1) sample_id.append("+" + col[1]);
        message.append("Check Phenotype file: " + pheno_file + "\n");
        message.append("Column Name of Sample ID: " + sample_id + "\n");
        message.append("Note: If the phenotype file does not contain a header, "
                       "the column name will be displayed as the Sample ID "
                       "which is ok.\n");
        bool found = false;
        std::unordered_map<std::string, bool> dup_col;
        if (pheno_header.size() == 0) {
            pheno_info.use_pheno = true;
            pheno_info.col.push_back(1 + !m_ignore_fid);

            pheno_info.name.push_back("");
            pheno_info.order.push_back(0);
            pheno_info.binary.push_back(c_commander.is_binary(0));
            message.append("Phenotype Name: " + col[pheno_info.col.back()]
                           + "\n");
        }
        else
        {
            for (size_t i_pheno = 0; i_pheno < pheno_header.size(); ++i_pheno) {
                if (dup_col.find(pheno_header[i_pheno]) == dup_col.end()) {
                    found = false;
                    dup_col[pheno_header[i_pheno]] = true;
                    // start from 1+!m_ignore_fid to skip the iid and fid part
                    for (size_t i_column = 1 + !m_ignore_fid;
                         i_column < col.size(); ++i_column)
                    {
                        if (col[i_column].compare(pheno_header[i_pheno]) == 0) {
                            found = true;
                            pheno_info.use_pheno = true;
                            pheno_info.col.push_back(i_column);
                            pheno_info.name.push_back(pheno_header[i_pheno]);
                            pheno_info.order.push_back(i_pheno);
                            pheno_info.binary.push_back(
                                c_commander.is_binary(i_pheno));
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
    // TODO: Might want to error out when duplicated column is detected within
    // the phenotype file
    size_t num_pheno = (pheno_info.use_pheno) ? pheno_info.col.size() : 1;
    message.append("There are a total of " + std::to_string(num_pheno)
                   + " phenotype to process\n");
    reporter.report(message);
}

void PRSice::init_matrix(const Commander& c_commander, const size_t pheno_index,
                         Genotype& target, Reporter& reporter,
                         const bool prslice)
{
    m_null_r2 = 0.0;
    m_phenotype = Eigen::VectorXd::Zero(0);
    g_independent_variables.resize(0, 0);
    m_sample_with_phenotypes.clear();


    const bool no_regress = c_commander.no_regress();
    const std::string pheno_file = c_commander.pheno_file();
    const std::string output_name = c_commander.out();

    // this reset the in_regression flag of all samples
    target.reset_sample_pheno();
    // this includes all samples

    gen_pheno_vec(target, pheno_file, pheno_index, !no_regress, reporter);
    if (!no_regress) {
        std::vector<std::string> cov_header = c_commander.get_cov_header();
        gen_cov_matrix(c_commander.get_cov_file(), cov_header, reporter);
    }
    // NOTE: After gen_cov_matrix, the has_pheno flag in m_sample_names is no
    // longer correct

    // now inform PRSice which samples should be included
    update_sample_included(target);

    // get the null r2
    double null_r2_adjust = 0.0;
    int n_thread = c_commander.thread();
    if (g_independent_variables.cols() > 2 && !no_regress) {
        assert(g_independent_variables.rows() == m_phenotype.rows());
        if (c_commander.is_binary(pheno_index)) {
            // ignore the first column
            // this is ok as both the first column (intercept) and the
            // second column (PRS) is currently 1
            Regression::glm(m_phenotype,
                            g_independent_variables.topRightCorner(
                                g_independent_variables.rows(),
                                g_independent_variables.cols() - 1),
                            m_null_p, m_null_r2, m_null_coeff, m_null_se, 25,
                            n_thread, true);
        }
        else
        {
            // ignore the first column
            Regression::linear_regression(
                m_phenotype,
                g_independent_variables.topRightCorner(
                    g_independent_variables.rows(),
                    g_independent_variables.cols() - 1),
                m_null_p, m_null_r2, null_r2_adjust, m_null_coeff, m_null_se,
                n_thread, true);
        }
    }
}

void PRSice::update_sample_included(Genotype& target)
{
    m_max_fid_length = 3;
    m_max_iid_length = 3;
    // anyone that's included in the study are considered
    // therefore, it should work even for multiple different
    // phenotypes
    for (size_t i_sample = 0; i_sample < target.num_sample(); ++i_sample) {
        bool included = target.is_include(i_sample);
        if (!included) continue;
        m_max_fid_length = (m_max_fid_length > target.fid(i_sample).length())
                               ? m_max_fid_length
                               : target.fid(i_sample).length();
        m_max_iid_length = (m_max_iid_length > target.iid(i_sample).length())
                               ? m_max_iid_length
                               : target.iid(i_sample).length();

        if (m_sample_with_phenotypes.find(target.sample_id(i_sample))
            == m_sample_with_phenotypes.end())
        {
            target.set_in_regression(i_sample, false);
        }
    }
}

void PRSice::gen_pheno_vec(Genotype& target, const std::string& pheno_file_name,
                           const int pheno_index, bool regress,
                           Reporter& reporter)
{
    std::vector<double> pheno_store;
    // reserve the maximum size (All samples)
    pheno_store.reserve(target.num_sample());
    const bool binary = pheno_info.binary[pheno_index];
    int max_num = 0;
    int num_case = 0;
    int num_control = 0;
    size_t invalid_pheno = 0;
    size_t num_not_found = 0;
    std::string line;
    size_t sample_index_ct = 0;
    size_t num_included = 0;
    std::unordered_set<double> input_sanity_check; // check if input is sensible
    if (pheno_info.use_pheno)                      // use phenotype file
    {
        int pheno_col_index =
            pheno_info.col[pheno_index]; // obtain the phenotype index
        std::ifstream pheno_file;
        pheno_file.open(pheno_file_name.c_str());
        if (!pheno_file.is_open()) {
            std::string error_message =
                "Cannot open phenotype file: " + pheno_file_name;
            throw std::runtime_error(error_message);
        }

        // Read in phenotype from phenotype file
        std::unordered_map<std::string, std::string> phenotype_info;
        // do not remove header line as that won't match anyway
        while (std::getline(pheno_file, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            std::vector<std::string> token = misc::split(line);
            if (token.size()
                <= (size_t)(pheno_index + 1
                            + !m_ignore_fid)) // need to check the range
            {
                std::string error_message =
                    "Malformed pheno file, should contain at least "
                    + std::to_string(pheno_index + 2 + !m_ignore_fid)
                    + " columns. "
                      "Have you use the --ignore-fid option?";
                throw std::runtime_error(error_message);
            }
            std::string id =
                (m_ignore_fid) ? token[0] : token[0] + "_" + token[1];
            phenotype_info[id] = token[pheno_col_index];
        }
        pheno_file.close();
        for (size_t i_sample = 0; i_sample < target.num_sample(); ++i_sample) {
            std::string id = target.sample_id(i_sample);
            bool included = target.is_include(i_sample);
            if (included) num_included++;
            if (phenotype_info.find(id) != phenotype_info.end() && included
                && phenotype_info[id].compare("NA") != 0)
            {
                try
                {
                    if (binary) {
                        int temp = misc::convert<int>(phenotype_info[id]);
                        if (temp >= 0 && temp <= 2) {
                            pheno_store.push_back(temp);
                            max_num = (temp > max_num) ? temp : max_num;
                            num_case += (temp == 1);
                            num_control += (temp == 0);
                        }
                        else
                        {
                            // so that it will add invalid
                            throw std::runtime_error(
                                "Invalid binary phenotype format!");
                        }
                    }
                    else
                    {
                        pheno_store.push_back(
                            misc::convert<double>(phenotype_info[id]));
                        if (input_sanity_check.size() < 2) {
                            input_sanity_check.insert(pheno_store.back());
                        }
                    }
                    m_sample_with_phenotypes[id] = sample_index_ct++;
                    target.set_in_regression(i_sample, true);
                }
                catch (const std::runtime_error& error)
                {
                    invalid_pheno++;
                }
            }
            else
            {
                num_not_found++;
            }
        }
    }
    else
    {
        // No phenotype file is provided
        // Use information from the fam file directly
        for (size_t i_sample = 0; i_sample < target.num_sample(); ++i_sample) {
            bool included = target.is_include(i_sample);
            if (included) num_included++;
            if (target.pheno_is_na(i_sample) || !included) {
                // it is ok to skip NA as default = sample.has_pheno = false
                continue;
            }
            try
            {
                if (binary) {
                    int temp = misc::convert<int>(target.pheno(i_sample));
                    if (temp >= 0 && temp <= 2) {
                        pheno_store.push_back(temp);
                        max_num = (temp > max_num) ? temp : max_num;
                        num_case += (temp == 1);
                        num_control += (temp == 0);
                    }
                    else
                    {
                        throw std::runtime_error(
                            "Invalid binary phenotype format!");
                    }
                }
                else
                {
                    pheno_store.push_back(
                        misc::convert<double>(target.pheno(i_sample)));
                    if (input_sanity_check.size() < 2) {
                        input_sanity_check.insert(pheno_store.back());
                    }
                }
                m_sample_with_phenotypes[target.sample_id(i_sample)] =
                    sample_index_ct++;
                target.set_in_regression(i_sample, true);
            }
            catch (const std::runtime_error& error)
            {
                invalid_pheno++;
            }
        }
    }

    std::string message = "";
    if (num_not_found != 0) {
        message.append(std::to_string(num_not_found)
                       + " sample(s) without phenotype\n");
    }
    if (invalid_pheno != 0) {
        message.append(std::to_string(invalid_pheno)
                       + " sample(s) with invalid phenotype\n");
    }
    if (num_not_found == num_included && regress) {
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
        reporter.report(message);
        throw std::runtime_error("Error: No sample left");
    }
    if (invalid_pheno == num_included && regress) {
        message.append("Error: All sample has invalid phenotypes!");
        reporter.report(message);
        throw std::runtime_error("Error: No sample left");
    }
    if (input_sanity_check.size() < 2 && !binary && regress) {
        message.append("Only one phenotype value detected");
        auto itr = input_sanity_check.begin();
        if ((*itr) == -9) {
            message.append(" and they are all -9");
        }
        reporter.report(message);
        throw std::runtime_error("Not enough valid phenotype");
    }
    bool error = false;
    if (max_num > 1 && binary) {
        num_case = 0;
        num_control = 0;
        size_t check = 0;
        for (auto&& pheno : pheno_store) {
            pheno--;
            if (pheno < 0) {
                error = true;
            }
            else
                (pheno == 1) ? num_case++ : num_control++;
            check++;
        }
    }
    if (error && regress) {
        reporter.report(message);
        throw std::runtime_error(
            "Mixed encoding! Both 0/1 and 1/2 encoding found!");
    }
    if (pheno_store.size() == 0 && regress) {
        reporter.report(message);
        throw std::runtime_error("No phenotype presented");
    }
    // now store the vector into the m_phenotype vector
    m_phenotype =
        Eigen::Map<Eigen::VectorXd>(pheno_store.data(), pheno_store.size());


    if (binary) {
        message.append(std::to_string(num_control) + " control(s)\n");
        message.append(std::to_string(num_case) + " case(s)\n");
        if (regress) {
            if (num_control == 0)
                throw std::runtime_error("There are no control samples");
            if (num_case == 0) throw std::runtime_error("There are no cases");
        }
    }
    else
    {
        message.append(std::to_string(m_phenotype.rows())
                       + " sample(s) with valid phenotype\n");
    }
    reporter.report(message);
}


std::vector<size_t> PRSice::get_cov_index(const std::string& c_cov_file,
                                          std::vector<std::string>& cov_header,
                                          Reporter& reporter)
{
    std::vector<size_t> cov_index;
    std::ifstream cov;
    cov.open(c_cov_file.c_str());
    if (!cov.is_open()) {
        std::string error_message =
            "Error: Cannot open covariate file: " + c_cov_file;
        throw std::runtime_error(error_message);
    }
    std::string line;
    std::getline(cov, line);
    // obtain the header information of the covariate file
    if (line.empty())
        throw std::runtime_error("First line of covariate file is empty!");
    std::vector<std::string> token = misc::split(line);
    if (cov_header.size() == 0) {
        cov_header = token;
        // if no header is provided, we will use all the covariates included
        for (size_t i = 1 + !m_ignore_fid; i < token.size(); ++i)
            cov_index.push_back(i); // FID, therefore+1
    }
    else
    {
        // The command class should have phrased the covariate regrex.
        // we just need to perform the matching
        std::unordered_set<std::string> included;
        for (auto cov : cov_header) {
            if (cov.empty()) continue;
            // to avoid duplicated covariance
            if (included.find(cov) == included.end()) {
                included.insert(cov);
            }
        }
        //+1 when fid is include
        for (size_t i_header = 1 + !m_ignore_fid; i_header < token.size();
             ++i_header)
        {
            if (included.find(token[i_header]) != included.end()) {
                cov_index.push_back(i_header);
            }
        }
    }

    // while the cov_index is sorted, this index corresponds to the header line
    // so that is ok?
    std::sort(cov_index.begin(), cov_index.end());
    if (cov_index.size() == 0) {
        throw std::runtime_error("Error: No valid covariates!");
    }
    else
    {
        std::string message = "";
        if (cov_index.size() == 1) {
            message.append("1 valid covariate included\n");
        }
        else
            message.append(std::to_string(cov_index.size())
                           + " valid covariates included\n");

        reporter.report(message);
    }
    cov_header = token;
    return cov_index;
}

// Funcion to get the factors from the covariate file
// This function won't go live until we have good way to handle the factors
// e.g. determining the base factor & detecting if covariate is actually a
// factor
void PRSice::check_factor_cov(
    const std::string& c_cov_file, const std::vector<std::string>& c_cov_header,
    const std::vector<size_t>& cov_index,
    std::vector<std::unordered_map<std::string, int>>& factor_levels)
{
    std::ifstream cov;
    cov.open(c_cov_file.c_str());
    if (!cov.is_open()) {
        std::string error_message =
            "Error: Cannot open covariate file: " + c_cov_file;
        throw std::runtime_error(error_message);
    }
    std::string line;
    std::getline(cov, line); // remove header
    std::vector<std::unordered_map<std::string, int>> current_factors(
        cov_index.size());
    std::vector<size_t> convertable(cov_index.size(), 0);
    size_t max_index = cov_index.back() + 1;
    while (std::getline(cov, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (token.size() < max_index) {
            std::string error_message =
                "Error: Malformed covariate file, should contain at least "
                + std::to_string(max_index) + " column!";
            throw std::runtime_error(error_message);
        }
        std::string id = (m_ignore_fid) ? token[0] : token[0] + "_" + token[1];
        if (m_sample_with_phenotypes.find(id) != m_sample_with_phenotypes.end())
        { // sample is found in the phenotype vector
            for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov) {
                size_t covar_index = cov_index[i_cov];
                if (current_factors[i_cov].find(token[covar_index])
                    != current_factors[i_cov].end())
                {
                    current_factors[i_cov][token[covar_index]]++;
                }
                else
                {
                    current_factors[i_cov][token[covar_index]] = 1;
                }
                try
                {
                    // this is for catching unconvertable covariate
                    misc::convert<double>(token[covar_index]);
                    convertable[i_cov]++;
                }
                catch (const std::runtime_error& error)
                {
                    std::string str = token[covar_index];
                    std::transform(str.begin(), str.end(), str.begin(),
                                   ::toupper);
                    // we also consider missing as convertable
                    if (str.compare("NA") == 0 || str.compare("NULL") == 0)
                        convertable[i_cov]++;
                }
            }
        }
    }
    cov.close();
    factor_levels.resize(cov_index.size()); // make sure size is ok
    size_t num_sample = m_sample_with_phenotypes.size();
    std::ofstream log_file_stream;
    log_file_stream.open(m_log_file.c_str(), std::ofstream::app);
    if (!log_file_stream.is_open()) {
        std::string error_message =
            "Error: Cannot open log file: " + m_log_file;
        throw std::runtime_error(error_message);
    }

    for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov) {
        // if all convertable,then it is not a factor
        if (convertable[i_cov] == num_sample) continue;
        factor_levels[i_cov] = current_factors[i_cov];
        log_file_stream << c_cov_header[cov_index[i_cov]]
                        << " is a factor with " << factor_levels[i_cov].size()
                        << " levels" << std::endl;
    }
    log_file_stream << std::endl;
    log_file_stream.close();
}

void PRSice::gen_cov_matrix(const std::string& c_cov_file,
                            std::vector<std::string>& cov_header,
                            Reporter& reporter)
{
    // The size of the map should be informative of the number of sample
    size_t num_sample = m_sample_with_phenotypes.size();
    if (c_cov_file.empty()) {
        // if no covariates, just return a matrix of 1
        g_independent_variables = Eigen::MatrixXd::Ones(num_sample, 2);
        return;
    }
    // obtain the index of each covariate

    std::vector<size_t> cov_index =
        get_cov_index(c_cov_file, cov_header, reporter);

    std::string message = "Processing the covariate file: " + c_cov_file + "\n";
    message.append("==============================\n");
    reporter.report(message);
    std::vector<std::pair<std::string, size_t>> valid_sample_index;
    // Initialize the independent variables matrix with 1s
    // might be worth while to check the covariates
    // not live yet. Wait till we have time to debug it
    // std::vector<std::unordered_map<std::string, int>> factor_levels;
    // check_factor_cov( c_cov_file,c_cov_header, cov_index, factor_levels);

    g_independent_variables =
        Eigen::MatrixXd::Ones(num_sample, cov_index.size() + 2);
    bool valid = true;
    size_t num_valid = 0;
    std::ifstream cov;
    cov.open(c_cov_file.c_str());
    if (!cov.is_open()) {
        std::string error_message =
            "Error: Cannot open covariate file: " + c_cov_file;
        throw std::runtime_error(error_message);
    }
    std::string line;
    std::getline(cov, line); // remove header
    size_t max_index = cov_index.back() + 1;
    // Check the number of missingness for each covariates
    std::vector<int> missing_count(cov_index.size(), 0);

    while (std::getline(cov, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        valid = true;
        std::vector<std::string> token = misc::split(line);
        if (token.size() < max_index) {
            std::string error_message =
                "Error: Malformed covariate file, should contain at least "
                + std::to_string(max_index) + " column!";
            throw std::runtime_error(error_message);
        }
        std::string id = (m_ignore_fid) ? token[0] : token[0] + "_" + token[1];
        if (m_sample_with_phenotypes.find(id) != m_sample_with_phenotypes.end())
        {
            // sample is found in the phenotype vector
            int index = m_sample_with_phenotypes[id]; // index on vector
            for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov) {
                if (token[cov_index[i_cov]].compare("NA") == 0
                    || token[cov_index[i_cov]].compare("Na") == 0
                    || token[cov_index[i_cov]].compare("na") == 0
                    || token[cov_index[i_cov]].compare("nA") == 0)
                {
                    valid = false;
                    g_independent_variables(index, i_cov + 2) = 0;
                    missing_count[i_cov]++;
                }
                else
                {
                    try
                    {
                        double temp =
                            misc::convert<double>(token[cov_index[i_cov]]);
                        g_independent_variables(index, i_cov + 2) = temp;
                        // + 2 because first line = intercept, second
                        // line = PRS
                    }
                    catch (const std::runtime_error& error)
                    {
                        valid = false;
                        // place holder as 0, will remove it later
                        g_independent_variables(index, i_cov + 2) = 0;
                        missing_count[i_cov]++;
                    }
                }
            }
            if (valid) {
                valid_sample_index.push_back(
                    std::pair<std::string, size_t>(id, index));
                num_valid++;
            }
        }
    }

    // update the phenotype and independent variable matrix to
    // remove missing samples
    // we don't bother imputing the value for the user as of now

    if (valid_sample_index.size() != num_sample && num_sample != 0) {
        // helpful to give the overview
        int removed = num_sample - valid_sample_index.size();
        message =
            std::to_string(removed) + " sample(s) with invalid covariate:\n\n";
        message.append("Covariate\tNumber of Missing Samples\n");
        for (size_t miss = 0; miss < missing_count.size(); ++miss) {
            message.append(cov_header[cov_index[miss]] + "\t"
                           + std::to_string(missing_count[miss]) + "\n");
        }
        double portion = (double) removed / (double) num_sample;
        if (valid_sample_index.size() == 0) {
            // if all samples are removed
            for (size_t miss = 0; miss < missing_count.size(); ++miss) {
                if (missing_count[miss] == num_sample) {
                    // we sorted the column index so we can't tell what the
                    // column name is useless we also store the head of the file
                    // (too troublesome)
                    message.append("Error: Column " + std::to_string(miss)
                                   + " is invalid, please check it is of the "
                                     "correct format\n");
                }
            }
            reporter.report(message);
            throw std::runtime_error("Error: All samples removed due to "
                                     "missingness in covariate file!");
        }
        if (portion > 0.05) {
            message.append(
                "Warning: More than " + std::to_string(portion * 100)
                + "% of your samples were removed! "
                  "You should check if your covariate file is correct\n");
        }
        reporter.report(message);
        // sort the sample index
        std::sort(begin(valid_sample_index), end(valid_sample_index),
                  [](std::pair<std::string, size_t> const& t1,
                     std::pair<std::string, size_t> const& t2) {
                      if (std::get<1>(t1) == std::get<1>(t2))
                          return std::get<0>(t1).compare(std::get<0>(t2)) < 0;
                      else
                          return std::get<1>(t1) < std::get<1>(t2);
                  });

        // update the m_phenotype and g_independent
        m_sample_with_phenotypes.clear();
        for (size_t cur_index = 0; cur_index < valid_sample_index.size();
             ++cur_index)
        {
            std::string name = std::get<0>(valid_sample_index[cur_index]);
            m_sample_with_phenotypes[name] = cur_index;
            size_t original_index = std::get<1>(valid_sample_index[cur_index]);
            if (original_index != cur_index) {
                m_phenotype(cur_index, 0) = m_phenotype(original_index, 0);
                for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov) {
                    g_independent_variables(cur_index, i_cov + 2) =
                        g_independent_variables(original_index, i_cov + 2);
                }
            }
        }
        g_independent_variables.conservativeResize(
            valid_sample_index.size(), g_independent_variables.cols());
        m_phenotype.conservativeResize(valid_sample_index.size(), 1);
    }
    message = "After reading the covariate file, "
              + std::to_string(valid_sample_index.size())
              + " sample(s) included in the analysis\n";
    reporter.report(message);
}

void PRSice::run_prsice(const Commander& c_commander,
                        const std::string& region_name,
                        const size_t pheno_index, const size_t region_index,
                        Genotype& target)
{

    target.reset_sample_prs();
    // prslice can easily be implemented using PRSet functionality
    // so maybe remove prslice from this function
    const bool no_regress = c_commander.no_regress();
    const bool print_all_scores = c_commander.all_scores();
    const int num_thread = c_commander.thread();
    const bool multi = pheno_info.name.size() > 1;
    const size_t num_samples_included = target.num_sample();
    Eigen::initParallel();
    Eigen::setNbThreads(num_thread);
    m_best_index = -1;
    m_num_snp_included = 0;
    g_perm_result.resize(m_num_perm, 2);
    m_prs_results.clear();
    m_best_sample_score.clear();
    m_prs_results.resize(target.num_threshold());
    // set to -1 to indicate not done
    for (auto&& p : m_prs_results) p.threshold = -1;
    // initialize score vector
    m_best_sample_score.resize(target.num_sample());

    // now prepare all score
    // in theory, we only need to calulate it once for every phenotype + sets
    // but it is easier to do it this way
    std::fstream all_out;
    if (print_all_scores) {
        std::string all_out_name = c_commander.out();
        if (multi) {
            all_out_name.append("." + pheno_info.name[pheno_index]);
        }
        all_out_name.append(".all.score");
        all_out.open(all_out_name.c_str(), std::fstream::out | std::fstream::in
                                               | std::fstream::ate);
        if (!all_out.is_open()) {
            std::string error_message =
                "Cannot open file " + all_out_name + " for write";
            throw std::runtime_error(error_message);
        }
    }


    // current threshold iteration
    size_t iter_threshold = 0;
    // +1 such that only 100% when finished
    size_t max_category = target.max_category() + 1;
    int cur_category = 0, cur_index = -1;
    double cur_threshold = 0.0, prev_progress = 0.0;
    bool require_standardize = (m_score == SCORING::STANDARDIZE);
    while (target.get_score(cur_index, cur_category, cur_threshold,
                            m_num_snp_included, region_index,
                            require_standardize))
    {
        double progress =
            (double) cur_category / (double) (max_category) *100.0;
        if (progress - prev_progress > 0.01 && !m_prset) {
            fprintf(stderr, "\rProcessing %03.2f%%", progress);
            prev_progress = progress;
        }

        if (print_all_scores) {
            for (size_t sample = 0; sample < num_samples_included; ++sample) {
                double score =
                    target.calculate_score(m_score, sample);
                size_t loc = m_all_file.header_length+sample*m_all_file.line_width+
                		m_all_file.skip_column_length+m_all_file.processed_threshold+m_all_file.processed_threshold*m_numeric_width;
                all_out.seekp(loc);
                all_out << std::setprecision(m_precision) << score;
            }
        }
        if (no_regress) {
            iter_threshold++;
            continue;
        }
        regress_score(target, cur_threshold, num_thread, pheno_index,
                      iter_threshold);

        if (c_commander.permutation() != 0) {
            permutation(target, num_thread,
                        c_commander.logit_perm()
                            && m_target_binary[pheno_index]);
        }
        iter_threshold++;
        m_all_file.processed_threshold++;
    }
    if (all_out.is_open()) all_out.close();
    if (!m_prset) fprintf(stderr, "\rProcessing %03.2f%%\n", 100.0);
    process_permutations();
    if(!no_regress){
    	print_best(target, pheno_index, c_commander);
    	if(m_prset){
    		// if it is prset based permutation, we wouldn't have
    		// reserved the memory and due to the special usaage,
    		// we'd like to have a separate memory pool for each thread.
    		// That can be done by having each thread allocating their
    		// own memory, though that might be more complicated than I'd
    		// like as we might also need to take into account the number
    		// of situation that can be done in one go
    	}
    }


}

void PRSice::print_best(Genotype &target, const size_t pheno_index, const Commander &commander){

	std::string pheno_name =
	        (pheno_info.name.size() > 1) ? pheno_info.name[pheno_index] : "";
	std::string output_prefix = commander.out();
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    std::string out_best = output_prefix + ".best";
    std::fstream best_out(out_best.c_str(), std::fstream::out | std::fstream::in
            | std::fstream::ate);
	auto&& best_info = m_prs_results[m_best_index];
	int best_snp_size = best_info.num_snp;
	if (best_snp_size == 0) {
		fprintf(stderr, "Error: Best R2 obtained when no SNPs were included\n");
		fprintf(stderr, "       Cannot output the best PRS score\n");
	}
	else
	{
		for (size_t sample = 0; sample < target.num_sample(); ++sample) {
			// samples that are extracted are ignored
			// sample excluded will not be output here
			std::string has_pheno =
					target.sample_in_regression(sample) ? "Yes" : "No";

            size_t loc = m_best_file.header_length+sample*m_best_file.line_width+
            		m_best_file.skip_column_length+m_best_file.processed_threshold+m_best_file.processed_threshold*m_numeric_width;
            best_out.seekp(loc);
            best_out << std::setprecision(m_precision) <<  m_best_sample_score[sample];
		}
	}
	best_out.close();
	m_best_file.processed_threshold++;
}
void PRSice::regress_score(Genotype& target, const double threshold,
                           size_t thread, const size_t pheno_index,
                           const size_t iter_threshold)
{
    double r2 = 0.0, r2_adjust = 0.0, p_value = 0.0, coefficient = 0.0,
           se = 0.0;
    const size_t num_include_samples = target.num_sample();
    if (m_num_snp_included == 0
        || (m_num_snp_included == m_prs_results[iter_threshold].num_snp))
    {
        return; // didn't got extra SNPs to process
    }

    for (size_t sample_id = 0; sample_id < num_include_samples; ++sample_id) {
        std::string sample = target.sample_id(sample_id);
        if (m_sample_with_phenotypes.find(sample)
            != m_sample_with_phenotypes.end())
        {
            g_independent_variables(m_sample_with_phenotypes.at(sample), 1) =
                target.calculate_score(m_score, sample_id);
        }
    }

    if (m_target_binary[pheno_index]) {
        try
        {
            Regression::glm(m_phenotype, g_independent_variables, p_value, r2,
                            coefficient, se, 25, thread, true);
        }
        catch (const std::runtime_error& error)
        {
            // This should only happen when the glm doesn't converge.
            // Let's hope that won't happen...
            fprintf(stderr, "Error: GLM model did not converge!\n");
            fprintf(stderr, "       Please send me the DEBUG files\n");
            std::ofstream debug;
            debug.open("DEBUG");
            debug << g_independent_variables << std::endl;
            debug.close();
            debug.open("DEBUG.y");
            debug << m_phenotype << std::endl;
            debug.close();
            fprintf(stderr, "Error: %s\n", error.what());
        }
    }
    else
    {
        Regression::linear_regression(m_phenotype, g_independent_variables,
                                      p_value, r2, r2_adjust, coefficient, se,
                                      thread, true);
    }

    // If this is the best r2, then we will add it
    int best_index = m_best_index;
    if (iter_threshold == 0 || best_index < 0
        || m_prs_results[best_index].r2 < r2)
    {
        m_best_index = iter_threshold;
        // can't remember why I can't just copy the whole vector
        for (size_t s = 0; s < num_include_samples; ++s) {
            m_best_sample_score[s] = target.calculate_score(m_score, s);
        }
    }
    // This should be thread safe as each thread will only mind their own
    // region now add the PRS result to the vectors (hopefully won't be out
    // off scope
    prsice_result cur_result;
    cur_result.threshold = threshold;
    cur_result.r2 = r2;
    cur_result.r2_adj = r2_adjust;
    cur_result.coefficient = coefficient;
    cur_result.p = p_value;
    cur_result.emp_p = -1.0;
    cur_result.num_snp = m_num_snp_included;
    cur_result.se = se;
    m_prs_results[iter_threshold] = cur_result;
}


void PRSice::process_permutations()
{
    // can't generate an empirical p-value if there is no observed p-value
    if (m_best_index == -1) return;
    double best_p = m_prs_results[m_best_index].p;
    size_t num_better = 0;
    for (auto&& p : g_perm_result) num_better += (p <= best_p);
    m_prs_results[m_best_index].emp_p =
        (double) (num_better + 1.0) / (double) (m_num_perm + 1.0);
}

void PRSice::permutation(Genotype& target, const size_t n_thread,
                         bool logit_perm)
{

    int num_iter = m_num_perm / m_perm_per_slice;

    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm_matrix(
        m_phenotype.rows());
    size_t num_include_samples = target.num_sample();
    Eigen::setNbThreads(n_thread);
    for (size_t sample_id = 0; sample_id < num_include_samples; ++sample_id) {
        std::string sample = target.sample_id(sample_id);

        if (m_sample_with_phenotypes.find(sample)
            != m_sample_with_phenotypes.end())
        {
            g_independent_variables(m_sample_with_phenotypes[sample], 1) =
                target.calculate_score(m_score, sample_id);
        }
    }
    int rank = 0;
    // logit_perm can only be true if it is binary trait and user used the
    // --logit-perm flag
    if (!logit_perm) {
        g_perm_pre_decomposed.compute(g_independent_variables);
        rank = g_perm_pre_decomposed.rank();
        Eigen::MatrixXd R = g_perm_pre_decomposed.matrixR()
                                .topLeftCorner(rank, rank)
                                .triangularView<Eigen::Upper>();
        g_pre_se_calulated = (R.transpose() * R).inverse().diagonal();
    }
    int cur_remain = m_remain_slice;
    // we reseed the random number generator in each iteration so that
    // we will always get the same pheontype permutation
    std::mt19937 rand_gen{m_seed};
    // need to do the permutation of phenotype without threading
    // so that we can reserve the sequence
    size_t processed = 0;
    // +1 so that we will also procdess the last slice
    for (int iter = 0; iter < num_iter + 1; ++iter) {
        double* perm_pheno_ptr = g_permuted_pheno.data();
        Eigen::setNbThreads(n_thread);
        size_t cur_perm = m_perm_per_slice;
        cur_perm += (cur_remain > 0) ? 1 : 0;
        cur_remain--;
        if (cur_perm + processed > m_num_perm) {
            cur_perm = m_num_perm - processed;
        }
        if (cur_perm < 1) break;
        // g_permuted_pheno.resize(cur_perm);
        for (size_t p = 0; p < cur_perm; ++p) {
            Eigen::Map<Eigen::VectorXd> perm_vec(perm_pheno_ptr,
                                                 m_phenotype.rows());
            perm_matrix.setIdentity();
            std::shuffle(perm_matrix.indices().data(),
                         perm_matrix.indices().data()
                             + perm_matrix.indices().size(),
                         rand_gen);
            // key point here: g_permuted_pheno is a vector
            perm_vec = perm_matrix * m_phenotype; // permute columns
            perm_pheno_ptr = &(perm_pheno_ptr[m_phenotype.rows()]);
        }
        Eigen::setNbThreads(1);
        // if want to use windows, we need to ditch the use of std::thread
        // now multithread it and get the corresponding p-values

        std::vector<pthread_t> pthread_store(n_thread);
        int job_size = cur_perm / n_thread;
        int remain = cur_perm % n_thread;
        size_t start = 0;
        g_perm_range.clear();

        for (uintptr_t ulii = 0; ulii < n_thread; ulii++) {
            size_t ending = start + job_size + (remain > 0);
            ending = (ending > cur_perm) ? cur_perm : ending;
            perm_info cur_info;
            cur_info.start = start;
            cur_info.end = ending;
            cur_info.rank = rank;
            cur_info.processed = processed;
            g_perm_range[ulii] = cur_info;
            start = ending;
            remain--;
            try
            {
#ifdef _WIN32
                pthread_store[ulii] = (HANDLE) _beginthreadex(
                    nullptr, 4096, thread_perm, (void*) ulii, 0, nullptr);
                if (!pthread_store[ulii]) {
                    join_all_threads(pthread_store.data(), ulii);
                    throw std::runtime_error(
                        "Error: Cannot create thread for permutation!");
                }
#else
                int error_code =
                    pthread_create(&(pthread_store[ulii]), nullptr,
                                   &PRSice::thread_perm, (void*) ulii);
                if (error_code) {

                    std::string error_message =
                        "Error: Cannot create thread for permutation! ("
                        + std::string(strerror(error_code)) + ")";
                    // join_threads(pthread_store.data(), ulii);
                    throw std::runtime_error(error_message);
                }
#endif
            }
            catch (std::exception& ex)
            {
                std::string error_message =
                    "Error: Cannot create thread for permutation!\n";
                error_message.append(ex.what());
                throw std::runtime_error(error_message);
            }
        }
        join_all_threads(pthread_store.data(), n_thread);


        processed += cur_perm;
    }
}

THREAD_RET_TYPE PRSice::thread_perm(void* id)
{
    size_t i_thread = (size_t) id;
    perm_info pi = g_perm_range[i_thread];
    size_t start = pi.start;
    size_t end = pi.end;
    size_t processed = pi.processed;
    size_t rank = pi.rank;
    bool intercept = true;
    size_t n = g_independent_variables.rows();
    std::vector<double> temp_store;
    temp_store.reserve(end - start);

    for (size_t i = start; i < end; ++i) {
        double* perm_pheno_ptr = g_permuted_pheno.data();
        perm_pheno_ptr = &(perm_pheno_ptr[i * g_independent_variables.rows()]);
        Eigen::Map<Eigen::VectorXd> perm_pheno(perm_pheno_ptr,
                                               g_independent_variables.rows());
        double obs_p = 2.0; // for safety reason, make sure it is out bound
        if (g_logit_perm) {
            double r2, coefficient, se;
            Regression::glm(perm_pheno, g_independent_variables, obs_p, r2,
                            coefficient, se, 25, 1, true);
        }
        else
        {
            Eigen::VectorXd beta = g_perm_pre_decomposed.solve(perm_pheno);
            Eigen::MatrixXd fitted = g_independent_variables * beta;

            Eigen::VectorXd residual = perm_pheno - fitted;
            int rdf = n - rank;
            double rss = 0.0;
            for (size_t r = 0; r < n; ++r) {
                rss += residual(r) * residual(r);
            }
            size_t se_index = intercept;
            for (size_t ind = 0; ind < (size_t) beta.rows(); ++ind) {
                if (g_perm_pre_decomposed.colsPermutation().indices()(ind)
                    == intercept)
                {
                    se_index = ind;
                    break;
                }
            }
            double resvar = rss / (double) rdf;
            Eigen::VectorXd se = (g_pre_se_calulated * resvar).array().sqrt();
            double tval = beta(intercept) / se(se_index);
            obs_p = misc::calc_tprob(tval, n);
        }
        // store the best p_value for the processed+i permutaiton
        // this is thread safe as we will never actually touch any overlapped
        // area
        temp_store.push_back(obs_p);
    }
    int index = 0;
    // this might seems odd, but we put it here to minimize false sharing (best
    // if mutex)
    for (size_t i = start; i < end; ++i) {
        double obs_p = temp_store[index++];
        double ori_p = g_perm_result[processed + i];
        g_perm_result[processed + i] = (ori_p > obs_p) ? obs_p : ori_p;
    }
    THREAD_RETURN;
}


void PRSice::prep_output(const Commander& c_commander, Genotype &target, std::vector<std::string> region_name, const size_t pheno_index){
	// As R has a default precision of 7, we will go a bit
	// higher to ensure we use up all precision
	std::string pheno_name =
	        (pheno_info.name.size() > 1) ? pheno_info.name[pheno_index] : "";
	std::string output_prefix = c_commander.out();
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    const bool perm = (c_commander.permutation() != 0);
    std::string output_name = output_prefix;
    std::string out_prsice = output_name + ".prsice";
    std::string out_all = output_name + ".all.score";
    std::string out_best = output_name + ".best";
    std::ofstream prsice_out, best_out, all_out;

    // .prsice output
    prsice_out.open(out_prsice.c_str());
    if (!prsice_out.is_open()) {
        std::string error_message =
            "Error: Cannot open file: " + out_prsice + " to write";
        throw std::runtime_error(error_message);
    }
    prsice_out << "Set\tThreshold\tR2\tP\tCoefficient\tStandard.Error\tNum_SNP";
    if (m_prset) prsice_out << "\tCompetitive_P";
    if (perm) prsice_out << "\tEmpirical_P";
    prsice_out << std::endl;
    prsice_out.close();

    // .best output
    best_out.open(out_best.c_str());
    if (!best_out.is_open()) {
    	std::string error_message =
    			"Error: Cannot open file: " + out_best + " to write";
    	throw std::runtime_error(error_message);
    }
    std::string header_line = "FID IID In_Regression";
    // if not preset, then it is PRS,otherwise, it will be the
    if(!m_prset) header_line.append(" PRS");
    else{
    	for(size_t i = 0; i < region_name.size()-1; ++i){
    		header_line.append(" "+region_name[i]);
    	}
    }
    best_out << header_line << std::endl;
    m_best_file.header_length = header_line.length()+1;
    m_best_file.processed_threshold = 0;
    // each numeric output took 12 spaces, then for each output, there is one space next to each
    m_best_file.line_width = region_name.size()*m_numeric_width + region_name.size()
    		+ m_max_fid_length +1+ m_max_iid_length+1;
    m_best_file.skip_column_length = m_max_fid_length+m_max_iid_length+3+3;


    // also handle all score here
    const bool all_scores = c_commander.all_scores();
    if(all_scores){
    	all_out.open(out_all.c_str());
    	if (!all_out.is_open()) {
    		std::string error_message =
    				"Cannot open file " + out_all + " for write";
    		throw std::runtime_error(error_message);
    	}
    	std::vector<double> avail_thresholds = target.get_thresholds();
    	std::sort(avail_thresholds.begin(), avail_thresholds.end());
    	size_t num_thresholds = avail_thresholds.size();
    	header_line = "FID IID";
    	if(!m_prset){
            for (auto& thres : avail_thresholds) {
            	header_line.append(" " + std::to_string(thres));
            }
    	}else{
        	for(size_t i = 0; i < region_name.size()-1; ++i){
        		for (auto& thres : avail_thresholds) {
            		header_line.append(" " + region_name[i]+"_"+std::to_string(thres));
            	}
            }
    	}
    	m_all_file.header_length = header_line.length()+1;
    	m_all_file.processed_threshold = 0;
    	m_all_file.line_width = num_thresholds*region_name.size()*m_numeric_width+
    			num_thresholds*region_name.size()+
				m_max_fid_length+1+m_max_iid_length+1;
    	m_all_file.skip_column_length = m_max_fid_length+m_max_iid_length+2;
    	all_out << header_line << std::endl;
    }

    // output sample IDs
    size_t num_samples_included = target.num_sample();
    for (size_t i_sample = 0; i_sample < num_samples_included; ++i_sample) {
    	std::string name =
    			target.fid(i_sample) + " " + target.iid(i_sample);
    	std::string best_line = name+" "+((target.sample_in_regression(i_sample))?"Yes":"No");
    	best_out << std::setfill(' ') << std::setw(m_best_file.line_width)
                       << std::left << best_line << std::endl;
    	if(all_scores){
        	all_out << std::setfill(' ') << std::setw(m_all_file.line_width)
                           << std::left << name << std::endl;
    	}
    }
    m_all_file.line_width++;
    m_best_file.line_width++;// now account for new line
    best_out.close();
    if(all_out.is_open()) all_out.close();

}
void PRSice::output(const Commander& c_commander, const Region& region,
                    const size_t pheno_index, const size_t region_index,
                    Genotype& target)
{
    std::vector<double> prev = c_commander.prevalence();
    bool has_prevalence = (prev.size() != 0);
    has_prevalence = has_prevalence && c_commander.is_binary(pheno_index);
    double top = 1.0, bottom = 1.0, prevalence = -1;
    if (has_prevalence) {
        size_t num_binary = 0;
        for (size_t i = 0; i < pheno_index; ++i) {
            if (c_commander.is_binary(i))
                num_binary++; // this is the number of previous binary traits
        }
        int num_case = 0, num_control = 0;
        for (size_t i = 0; i < (size_t) m_phenotype.rows(); ++i) {
            if (m_phenotype(i) == 0)
                num_control++;
            else if (m_phenotype(i) == 1)
                num_case++;
        }
        double case_ratio =
            (double) (num_case) / (double) (num_case + num_control);
        prevalence = prev[num_binary];
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

    std::string pheno_name =
        (pheno_info.name.size() > 1) ? pheno_info.name[pheno_index] : "";
    std::string output_prefix = c_commander.out();
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);

    const bool perm = (c_commander.permutation() != 0);

    bool valid = m_best_index != -1;
    if (!valid
        || region.get_count(region_index)
               == 0) // we know regions with 0 SNP will not have valid PRS
    {
        if (region.get_count(region_index) != 0) {
            fprintf(stderr, "Error: No valid PRS ");
            if (m_prset)
                fprintf(stderr, "for %s",
                        region.get_name(region_index).c_str());
            fprintf(stderr, "!\n");
        }
        return;
    }
    std::string out_prsice = output_prefix + ".prsice";
    std::string out_snp = output_prefix + ".snps";
    // std::string out_summary = output_name + ".summary";
    std::ofstream prsice_out, snp_out;
    prsice_out.open(out_prsice.c_str(), std::fstream::app);
    if (!prsice_out.is_open()) {
        std::string error_message =
            "Error: Cannot open file: " + out_prsice + " to write";
        throw std::runtime_error(error_message);
    }

    for (size_t i = 0; i < m_prs_results.size(); ++i) {
        if (m_prs_results[i].threshold < 0 || m_prs_results[i].p < 0) continue;
        double full = m_prs_results[i].r2;
        double null = m_null_r2;
        if (has_prevalence) {
            full = top * full / (1 + bottom * full);
            null = top * null / (1 + bottom * null);
        }
        double r2 = full - null;
        prsice_out << region.get_name(region_index) << "\t" << m_prs_results[i].threshold << "\t" << r2 << "\t"
                   << m_prs_results[i].p << "\t" << m_prs_results[i].coefficient
                   << "\t" << m_prs_results[i].se << "\t"
                   << m_prs_results[i].num_snp;
        if (perm)
            prsice_out << "\t"
                       << ((m_prs_results[i].emp_p >= 0.0)
                               ? std::to_string(m_prs_results[i].emp_p)
                               : "-");
        prsice_out << std::endl;
    }
    prsice_out.close();
    auto&& best_info = m_prs_results[m_best_index];


    prsice_summary prs_sum;
    prs_sum.pheno = pheno_name;
    prs_sum.set = region.get_name(region_index);
    prs_sum.result = best_info;
    prs_sum.result.r2 = best_info.r2;
    prs_sum.r2_null = m_null_r2;
    prs_sum.top = top;
    prs_sum.bottom = bottom;
    prs_sum.prevalence = prevalence;

    m_prs_summary.push_back(prs_sum);
    if (best_info.p > 0.1)
        m_significant_store[0]++;
    else if (best_info.p > 1e-5)
        m_significant_store[1]++;
    else
        m_significant_store[2]++;

    if (c_commander.print_snp()) {
        target.print_snp(out_snp, m_prs_results[m_best_index].threshold,
                         region_index);
    }
}

void PRSice::summarize(const Commander& commander, Reporter& reporter)
{
    bool prev_out = false;

    const bool perm = (commander.permutation() != 0);
    std::string message = "There are ";
    if (m_significant_store[0] != 0) {
        message.append(std::to_string(m_significant_store[0])
                       + " region(s) with p-value > 0.1 (\033[1;31mnot "
                         "significant\033[0m);");
        prev_out = true;
    }
    if (m_significant_store[1] != 0) {
        if (m_significant_store[2] == 0 && prev_out) {
            message.append(" and ");
        }
        message.append(
            std::to_string(m_significant_store[1])
            + " region(s) with p-value between "
              "0.1 and 1e-5 (\033[1;31mmay not be significant\033[0m);");
        prev_out = true;
    }
    if (m_significant_store[2] != 0) {
        if (prev_out) message.append(" and ");
        message.append(std::to_string(m_significant_store[2])
                       + " region(s) with p-value less than 1e-5.");
    }
    if (!perm) {
        message.append(
            " Please note that these results are inflated due to the "
            "overfitting inherent in finding the best-fit "
            "PRS (but it's still best to find the best-fit PRS!).\n"
            "You can use the --perm option (see manual) to calculate "
            "an empirical P-value.");
    }
    reporter.report(message);
    std::string out_name = commander.out() + ".summary";
    std::ofstream out;
    out.open(out_name.c_str());
    if (!out.is_open()) {
        std::string error_message =
            "Error: Cannot open file: " + out_name + " to write";
        throw std::runtime_error(error_message);
    }
    out << "Phenotype\tSet\tThreshold\tPRS.R2\tFull.R2\tNull."
           "R2\tPrevalence\tCoefficient\tStandard.Error\tP\tNum_SNP";
    if (perm) out << "\tEmpirical-P";
    out << std::endl;
    for (auto&& sum : m_prs_summary) {
        out << ((sum.pheno.empty()) ? "-" : sum.pheno) << "\t" << sum.set
            << "\t" << sum.result.threshold;
        if (sum.prevalence > 0) {
            double full = sum.result.r2;
            double null = sum.r2_null;
            full = sum.top * full / (1 + sum.bottom * full);
            null = sum.top * null / (1 + sum.bottom * null);
            out << "\t" << full - null << "\t" << full << "\t" << null << "\t"
                << sum.prevalence;
        }
        else
        {
            out << "\t" << sum.result.r2 - sum.r2_null << "\t" << sum.result.r2
                << "\t" << sum.r2_null << "\t-";
        }
        out << "\t" << sum.result.coefficient << "\t" << sum.result.se << "\t"
            << sum.result.p << "\t" << sum.result.num_snp;
        if (perm) out << "\t" << sum.result.emp_p;
        out << std::endl;
    }
    out.close();
}

PRSice::~PRSice()
{
    // dtor
}

void PRSice::gen_perm_memory(const Commander &commander,const size_t sample_ct, Reporter& reporter)
{
    intptr_t min_memory_byte = 8 * sample_ct;
    intptr_t max_req_memory = min_memory_byte * m_num_perm;
#ifdef __APPLE__
    int32_t mib[2];
    size_t sztmp;
#endif
    unsigned char* bigstack_ua = nullptr; // ua = unaligned
    int64_t llxx;
    intptr_t default_alloc_mb;
    intptr_t malloc_size_mb = 0;
#ifdef __APPLE__
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    llxx = 0;

    sztmp = sizeof(int64_t);
    sysctl(mib, 2, &llxx, &sztmp, nullptr, 0);
    llxx /= 1048576;
#else
#ifdef _WIN32
    MEMORYSTATUSEX memstatus;
    memstatus.dwLength = sizeof(memstatus);
    GlobalMemoryStatusEx(&memstatus);
    llxx = memstatus.ullTotalPhys / 1048576;
#else
    llxx = ((uint64_t) sysconf(_SC_PHYS_PAGES))
           * ((size_t) sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
    if (!llxx) {
        default_alloc_mb = BIGSTACK_DEFAULT_MB;
    }
    else if (llxx < (BIGSTACK_MIN_MB * 2))
    {
        default_alloc_mb = BIGSTACK_MIN_MB;
    }
    else
    {
        default_alloc_mb = llxx / 2;
    }
    if (!malloc_size_mb) {
        malloc_size_mb = default_alloc_mb;
    }
    else if (malloc_size_mb < BIGSTACK_MIN_MB)
    {
        malloc_size_mb = BIGSTACK_MIN_MB;
    }
    std::string message = "";
#ifndef __LP64__
    if (malloc_size_mb > 2047) {
        malloc_size_mb = 2047;
    }
#endif
    bigstack_ua =
        (unsigned char*) malloc(malloc_size_mb * 1048576 * sizeof(char));
    // if fail, return nullptr which will then get into the while loop
    while (!bigstack_ua) {
        malloc_size_mb = (malloc_size_mb * 3) / 4;
        if (malloc_size_mb < BIGSTACK_MIN_MB) {
            malloc_size_mb = BIGSTACK_MIN_MB;
        }
        bigstack_ua =
            (unsigned char*) malloc(malloc_size_mb * 1048576 * sizeof(char));
        if (bigstack_ua) {
        }
        else if (malloc_size_mb == BIGSTACK_MIN_MB)
        {
            throw std::runtime_error(
                "Failed to allocate required memory for permutation storage");
        }
    }
    if (malloc_size_mb * 1048576 < min_memory_byte) {
        throw std::runtime_error(
            "Failed to allocate required memory for permutation storage");
    }
    delete[] bigstack_ua;
    bigstack_ua = nullptr;
	// * 0.5 to provide room of error
    size_t valid_memory = commander.max_memory(malloc_size_mb*1048576);
    intptr_t final_mb = (valid_memory-misc::current_ram_usage())*0.5;
    // start update here
    if(final_mb < 0){
    	throw std::runtime_error("Error: Insufficient memory for permutation!");
    }
    if (final_mb < min_memory_byte) {
        m_perm_per_slice = 1;
    }
    else if (final_mb > max_req_memory)
    {
        m_perm_per_slice = m_num_perm;
    }
    else
    {
        m_perm_per_slice = final_mb  / min_memory_byte;
    }
    message = std::to_string(((final_mb>max_req_memory)?max_req_memory:final_mb)/1048576.0) + " MB RAM reserved for permutation\n";
    reporter.report(message);
    // wanna use double vector here as the sample size here might not be
    // the one used in the permutation. This might then lead to problem
    // in the permutation (segmentation fault, etc)
    g_permuted_pheno.resize(sample_ct * m_perm_per_slice);
    // g_num_snps.resize(g_max_threshold_store * m_sample_names.size(), 0);
    // g_prs_storage.resize(g_max_threshold_store * m_sample_names.size(), 0.0);
}
