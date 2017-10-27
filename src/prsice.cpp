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

std::mutex PRSice::score_mutex;


void PRSice::pheno_check(const Commander& c_commander)
{
    std::vector<std::string> pheno_header = c_commander.pheno_col();
    std::string pheno_file = c_commander.pheno_file();
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
                "ERROR: Not enough column in Phenotype file."
                "Have you use the --ignore-fid option");
        }
        std::string sample_id = col[0];
        if (!m_ignore_fid && col.size() > 1) sample_id.append("+" + col[1]);
        std::ofstream log_file_stream;
        log_file_stream.open(m_log_file.c_str(), std::ofstream::app);
        if (!log_file_stream.is_open()) {
            std::string error_message =
                "ERROR: Cannot open log file: " + m_log_file;
            throw std::runtime_error(error_message);
        }
        log_file_stream << "Check Phenotype file: " << pheno_file << std::endl;
        log_file_stream << "Column Name of Sample ID: " << sample_id
                        << std::endl;

        bool found = false;
        std::unordered_map<std::string, bool> dup_col;
        if (pheno_header.size() == 0) {
            pheno_info.use_pheno = true;
            pheno_info.col.push_back(1 + !m_ignore_fid);

            pheno_info.name.push_back("");
            pheno_info.order.push_back(0);
            pheno_info.binary.push_back(c_commander.is_binary(0));

            log_file_stream << "Phenotype Name: " << col[pheno_info.col.back()]
                            << std::endl;
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
                        fprintf(
                            stderr,
                            "Phenotype: %s cannot be found in phenotype file\n",
                            pheno_header[i_pheno].c_str());
                        log_file_stream
                            << "Phenotype: " << pheno_header[i_pheno]
                            << " cannot be found in phenotype file"
                            << std::endl;
                    }
                }
            }
        }
        log_file_stream << std::endl;
        log_file_stream.close();
    }
    size_t num_pheno = (pheno_info.use_pheno) ? pheno_info.col.size() : 1;
    fprintf(stderr, "There are a total of %zu phenotype to process\n",
            num_pheno);
}

void PRSice::update_sample_included()
{
    m_max_fid_length = 3;
    m_max_iid_length = 3;
    m_sample_included.clear();
    m_sample_index.clear();
    // anyone that's included in the study are considered
    // therefore, it should work even for multiple different
    // phenotypes
    for (size_t i_sample = 0; i_sample < m_sample_names.size(); ++i_sample) {
        auto&& sample = m_sample_names[i_sample];
        if (!sample.included) continue;
        m_max_fid_length = (m_max_fid_length > sample.FID.length())
                               ? m_max_fid_length
                               : sample.FID.length();
        m_max_iid_length = (m_max_iid_length > sample.IID.length())
                               ? m_max_iid_length
                               : sample.IID.length();

        std::string id =
            (m_ignore_fid) ? sample.IID : sample.FID + "_" + sample.IID;
        m_sample_included.push_back(id);
        m_sample_index.push_back(i_sample);
    }
}
void PRSice::init_matrix(const Commander& c_commander, const size_t pheno_index,
                         Genotype& target, const bool prslice)
{
    m_null_r2 = 0.0;
    m_phenotype = Eigen::VectorXd::Zero(0);
    m_independent_variables.resize(0, 0);
    // m_sample_names.clear();
    m_sample_with_phenotypes.clear();


    const bool no_regress = c_commander.no_regress();
    const std::string pheno_file = c_commander.pheno_file();
    const std::string output_name = c_commander.out();

    // this includes all samples
    if (m_sample_names.empty()) m_sample_names = target.sample_names();
    gen_pheno_vec(pheno_file, pheno_index, !no_regress);
    if (!no_regress) {
        std::vector<std::string> cov_header = c_commander.get_cov_header();
        gen_cov_matrix(c_commander.get_cov_file(), cov_header);
    }
    // now inform PRSice which samples should be included
    update_sample_included();

    // get the null r2
    double null_r2_adjust = 0.0, null_p = 0.0, null_coeff = 0.0;
    int n_thread = c_commander.thread();
    if (m_independent_variables.cols() > 2 && !no_regress) {
        assert(m_independent_variables.rows() == m_phenotype.rows());
        if (c_commander.is_binary(pheno_index)) {
            // ignore the first column
            // this is ok as both the first column (intercept) and the
            // second column (PRS) is currently 1
            Regression::glm(m_phenotype,
                            m_independent_variables.topRightCorner(
                                m_independent_variables.rows(),
                                m_independent_variables.cols() - 1),
                            null_p, m_null_r2, null_coeff, 25, n_thread, true);
        }
        else
        {
            // ignore the first column
            Regression::linear_regression(
                m_phenotype,
                m_independent_variables.topRightCorner(
                    m_independent_variables.rows(),
                    m_independent_variables.cols() - 1),
                null_p, m_null_r2, null_r2_adjust, null_coeff, n_thread, true);
        }
    }
}

void PRSice::gen_pheno_vec(const std::string& pheno_file_name,
                           const int pheno_index, bool regress)
{
    std::vector<double> pheno_store;
    // reserve the maximum size (All samples)
    pheno_store.reserve(m_sample_names.size());
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
        std::getline(pheno_file, line); // remove header line
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
                    + " columns\n"
                      "Have you use the --ignore-fid option?";
                throw std::runtime_error(error_message);
            }
            std::string id =
                (m_ignore_fid) ? token[0] : token[0] + "_" + token[1];
            phenotype_info[id] = token[pheno_col_index];
        }
        pheno_file.close();

        // Add phenotype information to each sample
        for (auto&& sample : m_sample_names) {
            std::string id =
                (m_ignore_fid) ? sample.IID : sample.FID + "_" + sample.IID;
            if (sample.included) num_included++;
            if (phenotype_info.find(id) != phenotype_info.end()
                && sample.included && phenotype_info[id].compare("NA") != 0)
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
                            throw std::runtime_error("");
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
                    sample.has_pheno = true;
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
        for (auto&& sample : m_sample_names) {
            if (sample.included) num_included++;
            if (sample.pheno.compare("NA") == 0 || !sample.included) {
                // it is ok to skip NA as default = sample.has_pheno = false
                continue;
            }
            try
            {
                if (binary) {
                    int temp = misc::convert<int>(sample.pheno);
                    if (temp >= 0 && temp <= 2) {
                        pheno_store.push_back(temp);
                        max_num = (temp > max_num) ? temp : max_num;
                        num_case += (temp == 1);
                        num_control += (temp == 0);
                    }
                    else
                    {
                        throw std::runtime_error("");
                    }
                }
                else
                {
                    pheno_store.push_back(misc::convert<double>(sample.pheno));
                    if (input_sanity_check.size() < 2) {
                        input_sanity_check.insert(pheno_store.back());
                    }
                }
                m_sample_with_phenotypes[m_ignore_fid
                                             ? sample.IID
                                             : sample.FID + "_" + sample.IID] =
                    sample_index_ct++;
                sample.has_pheno = true;
            }
            catch (const std::runtime_error& error)
            {
                invalid_pheno++;
            }
        }
    }

    std::ofstream log_file_stream;
    log_file_stream.open(m_log_file.c_str(), std::ofstream::app);
    if (!log_file_stream.is_open()) {
        std::string error_message =
            "ERROR: Cannot open log file: " + m_log_file;
        throw std::runtime_error(error_message);
    }


    if (num_not_found != 0) {
        log_file_stream << num_not_found << " sample(s) without phenotype"
                        << std::endl;
        fprintf(stderr, "Number of missing samples: %zu\n", num_not_found);
    }
    if (invalid_pheno != 0) {
        fprintf(stderr, "Number of invalid phenotyps: %zu\n", invalid_pheno);
        log_file_stream << invalid_pheno << " sample(s) with invalid phenotype"
                        << std::endl;
    }
    if (num_not_found == num_included) {

        log_file_stream
            << "None of the target samples were found in the phenotype file"
            << std::endl;
        fprintf(
            stderr,
            "None of the target samples were found in the phenotype file\n");
        if (m_ignore_fid) {
            fprintf(
                stderr,
                "Maybe the first column of your phenotype file is the FID?\n");
            log_file_stream
                << "Maybe the first column of your phenotype file is the FID?"
                << std::endl;
        }
        else
        {
            fprintf(stderr,
                    "Maybe your phenotype file doesn not contain the FID?\n");
            log_file_stream
                << "Maybe your phenotype file doesn not contain the FID?"
                << std::endl;
            fprintf(stderr, "Might consider using --ignore-fid\n");
            log_file_stream << "Might consider using --ignore-fid" << std::endl;
        }
        log_file_stream << std::endl;
        log_file_stream.close();
        throw std::runtime_error("ERROR: No sample left");
    }
    if (invalid_pheno == num_included) {
        log_file_stream << "ERROR: No sample left" << std::endl;
        log_file_stream << std::endl;
        log_file_stream.close();
        throw std::runtime_error("ERROR: No sample left");
    }
    if (input_sanity_check.size() < 2 && !binary) {
        fprintf(stderr, "Only one phenotype value detected\n");
        log_file_stream << "Only one phenotype value detected" << std::endl;
        auto itr = input_sanity_check.begin();
        if ((*itr) == -9) {
            fprintf(stderr, "and they are all -9\n");
            log_file_stream << "and they are all -9" << std::endl;
        }
        log_file_stream << "Not enough valid phenotype" << std::endl;
        log_file_stream << std::endl;
        log_file_stream.close();
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
        log_file_stream << "Mixed encoding! Both 0/1 and 1/2 encoding found!"
                        << std::endl;
        log_file_stream << std::endl;
        log_file_stream.close();
        throw std::runtime_error(
            "Mixed encoding! Both 0/1 and 1/2 encoding found!");
    }
    if (pheno_store.size() == 0 && regress) {
        log_file_stream << "No phenotype presented" << std::endl;
        log_file_stream << std::endl;
        log_file_stream.close();
        throw std::runtime_error("No phenotype presented");
    }
    // now store the vector into the m_phenotype vector
    m_phenotype =
        Eigen::Map<Eigen::VectorXd>(pheno_store.data(), pheno_store.size());


    if (binary) {
        log_file_stream << num_control << " control(s)" << std::endl;
        log_file_stream << num_case << " case(s)" << std::endl;
        fprintf(stderr, "Number of controls : %i\n", num_control);
        fprintf(stderr, "Number of cases : %i\n", num_case);
        if (regress) {
            if (num_control == 0)
                throw std::runtime_error("There are no control samples");
            if (num_case == 0) throw std::runtime_error("There are no cases");
        }
    }
    else
    {
        log_file_stream << m_phenotype.rows() << " sample(s) with phenotype"
                        << std::endl;
        fprintf(stderr, "Number of sample(s) with phenotype  : %zu\n",
                m_phenotype.rows());
    }
    log_file_stream << std::endl;
    log_file_stream.close();
}


std::vector<size_t> PRSice::get_cov_index(const std::string& c_cov_file,
                                          std::vector<std::string>& cov_header)
{
    std::vector<size_t> cov_index;
    std::ifstream cov;
    cov.open(c_cov_file.c_str());
    if (!cov.is_open()) {
        std::string error_message =
            "ERROR: Cannot open covariate file: " + c_cov_file;
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
        // same, +1 when fid is include
        for (size_t i_header = 1 + !m_ignore_fid; i_header < token.size();
             ++i_header)
        {
            if (included.find(token[i_header]) != included.end()) {
                cov_index.push_back(i_header);
            }
        }
    }
    std::ofstream log_file_stream;
    log_file_stream.open(m_log_file.c_str(), std::ofstream::app);
    if (!log_file_stream.is_open()) {
        std::string error_message =
            "ERROR: Cannot open log file: " + m_log_file;
        throw std::runtime_error(error_message);
    }

    std::sort(cov_index.begin(), cov_index.end());
    if (cov_index.size() == 0) {
        log_file_stream << "ERROR: No valid covariates!" << std::endl;
        log_file_stream << std::endl;
        log_file_stream.close();
        throw std::runtime_error("ERROR: No valid covariates!");
    }
    else
    {
        if (cov_index.size() == 1) {
            log_file_stream << "1 valid covariate included" << std::endl;
            fprintf(stderr, "1 valid covariate included\n");
        }
        else
            fprintf(stderr, "%zu valid covariates included\n",
                    cov_index.size());
        log_file_stream << cov_index.size() << " valid covariates included"
                        << std::endl;
    }
    log_file_stream << std::endl;
    log_file_stream.close();
    cov_header = token;
    return cov_index;
}


void PRSice::check_factor_cov(
    const std::string& c_cov_file, const std::vector<std::string>& c_cov_header,
    const std::vector<size_t>& cov_index,
    std::vector<std::unordered_map<std::string, int>>& factor_levels)
{
    std::ifstream cov;
    cov.open(c_cov_file.c_str());
    if (!cov.is_open()) {
        std::string error_message =
            "ERROR: Cannot open covariate file: " + c_cov_file;
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
                "ERROR: Malformed covariate file, should contain at least "
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
            "ERROR: Cannot open log file: " + m_log_file;
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
                            std::vector<std::string>& cov_header)
{
    // The size of the map should be informative of the number of sample
    size_t num_sample = m_sample_with_phenotypes.size();
    if (c_cov_file.empty()) {
        // if no covariates, just return a matrix of 1
        m_independent_variables = Eigen::MatrixXd::Ones(num_sample, 2);
        return;
    }
    // obtain the index of each covariate

    std::vector<size_t> cov_index = get_cov_index(c_cov_file, cov_header);

    std::ofstream log_file_stream;
    log_file_stream.open(m_log_file.c_str(), std::ofstream::app);
    if (!log_file_stream.is_open()) {
        std::string error_message =
            "ERROR: Cannot open log file: " + m_log_file;
        throw std::runtime_error(error_message);
    }
    log_file_stream << "Processing the covariate file: " << c_cov_file
                    << std::endl;
    // log_flie_stream.close(); // if check factor covariates
    fprintf(stderr, "\nStart processing the covariates\n");
    fprintf(stderr, "==============================\n");
    std::vector<std::pair<std::string, size_t>> valid_sample_index;
    // Initialize the independent variables matrix with 1s
    // might be worth while to check the covariates
    // not live yet. Wait till we have time to debug it
    // std::vector<std::unordered_map<std::string, int>> factor_levels;
    // check_factor_cov( c_cov_file,c_cov_header, cov_index, factor_levels);

    m_independent_variables =
        Eigen::MatrixXd::Ones(num_sample, cov_index.size() + 2);
    bool valid = true;
    size_t num_valid = 0;
    std::ifstream cov;
    cov.open(c_cov_file.c_str());
    if (!cov.is_open()) {
        std::string error_message =
            "ERROR: Cannot open covariate file: " + c_cov_file;
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
            log_file_stream.close();
            std::string error_message =
                "ERROR: Malformed covariate file, should contain at least "
                + std::to_string(max_index) + " column!";
            throw std::runtime_error(error_message);
        }
        std::string id = (m_ignore_fid) ? token[0] : token[0] + "_" + token[1];
        if (m_sample_with_phenotypes.find(id) != m_sample_with_phenotypes.end())
        {
            // sample is found in the phenotype vector
            int index = m_sample_with_phenotypes[id]; // index on vector
            for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov) {
                try
                {
                    double temp =
                        misc::convert<double>(token[cov_index[i_cov]]);
                    m_independent_variables(index, i_cov + 2) =
                        temp; // + 2 because first line = intercept, second line
                              // = PRS
                }
                catch (const std::runtime_error& error)
                {
                    valid = false;
                    // place holder as 0, will remove it later
                    m_independent_variables(index, i_cov + 2) = 0;
                    missing_count[i_cov]++;
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
        fprintf(stderr, "Number of samples with invalid covariate: %d\n",
                removed);
        log_file_stream << removed << " sample(s) with invalid covariate"
                        << std::endl;
        log_file_stream << "Covariate\tNumber of Missing Samples" << std::endl;
        for (size_t miss = 0; miss < missing_count.size(); ++miss) {
            log_file_stream << cov_header[cov_index[miss]] << "\t"
                            << missing_count[miss] << std::endl;
        }
        double portion = (double) removed / (double) num_sample;
        if (valid_sample_index.size() == 0) {
            // if all samples are removed
            for (size_t miss = 0; miss < missing_count.size(); ++miss) {
                if (missing_count[miss] == num_sample) {
                    // we sorted the column index so we can't tell what the
                    // column name is useless we also store the head of the file
                    // (too troublesome)
                    fprintf(stderr,
                            "Column %zu is invalid, please check it is of the "
                            "correct format\n",
                            miss);
                }
            }

            log_file_stream
                << "All samples removed due to missingness in covariate file!"
                << std::endl;
            log_file_stream << std::endl;
            log_file_stream.close();
            throw std::runtime_error(
                "All samples removed due to missingness in covariate file!");
        }
        if (portion > 0.05) {
            fprintf(
                stderr,
                "WARNING: More than %03.2f%% of the samples were removed!\n",
                portion * 100);

            fprintf(stderr,
                    "         Do check if your covariate file is correct\n");
            log_file_stream << "WARNING: More than " << portion * 100
                            << "% of the samples were removed!" << std::endl;
            log_file_stream
                << "         Do check if your covariate file is correct"
                << std::endl;
        }

        // sort the sample index
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
        for (size_t cur_index = 0; cur_index < valid_sample_index.size();
             ++cur_index)
        {
            std::string name = std::get<0>(valid_sample_index[cur_index]);
            m_sample_with_phenotypes[name] = cur_index;
            size_t original_index = std::get<1>(valid_sample_index[cur_index]);
            if (original_index != cur_index) {
                m_phenotype(cur_index, 0) = m_phenotype(original_index, 0);
                for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov) {
                    m_independent_variables(cur_index, i_cov + 2) =
                        m_independent_variables(original_index, i_cov + 2);
                }
            }
        }
        m_independent_variables.conservativeResize(
            valid_sample_index.size(), m_independent_variables.cols());
        m_phenotype.conservativeResize(valid_sample_index.size(), 1);

        fprintf(stderr, "\nFinal number of samples: %zu\n\n",
                valid_sample_index.size());
    }
    else
    {
        fprintf(stderr, "\nFinal number of samples: %zu\n\n",
                valid_sample_index.size());
    }

    log_file_stream << "After reading covariate file:" << std::endl;
    log_file_stream << valid_sample_index.size()
                    << " sample(s) included in the analysis" << std::endl;
    log_file_stream << std::endl;
    log_file_stream.close();
}


void PRSice::regress_score(const double threshold, size_t thread,
                           const size_t pheno_index,
                           const size_t iter_threshold)
{
    double r2 = 0.0, r2_adjust = 0.0, p_value = 0.0, coefficient = 0.0;
    size_t num_include_samples = m_current_sample_score.size();
    if (m_num_snp_included == 0
        || (m_num_snp_included == m_prs_results[iter_threshold].num_snp))
    {
        return; // didn't got extra SNPs to process
    }

    for (size_t sample_id = 0; sample_id < num_include_samples; ++sample_id) {
        std::string sample = m_sample_included[sample_id];
        if (m_sample_with_phenotypes.find(sample)
            != m_sample_with_phenotypes.end())
        {
            m_independent_variables(m_sample_with_phenotypes.at(sample), 1) =
                (m_average_score)
                    ? ((m_current_sample_score[sample_id].num_snp == 0)
                           ? 0
                           : m_current_sample_score[sample_id].prs
                                 / (double) m_current_sample_score[sample_id]
                                       .num_snp)
                    : m_current_sample_score[sample_id].prs;
        }
    }


    if (m_target_binary[pheno_index]) {
        try
        {
            Regression::glm(m_phenotype, m_independent_variables, p_value, r2,
                            coefficient, 25, thread, true);
        }
        catch (const std::runtime_error& error)
        {
            // This should only happen when the glm doesn't converge.
            // Let's hope that won't happen...
            fprintf(stderr, "ERROR: GLM model did not converge!\n");
            fprintf(stderr, "       Please send me the DEBUG files\n");
            std::ofstream debug;
            debug.open("DEBUG");
            debug << m_independent_variables << std::endl;
            debug.close();
            debug.open("DEBUG.y");
            debug << m_phenotype << std::endl;
            debug.close();
            fprintf(stderr, "ERROR: %s\n", error.what());
        }
    }
    else
    {
        Regression::linear_regression(m_phenotype, m_independent_variables,
                                      p_value, r2, r2_adjust, coefficient,
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
            m_best_sample_score[s] = m_current_sample_score[s];
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
    m_prs_results[iter_threshold] = cur_result;
}


void PRSice::run_prsice(const Commander& c_commander,
                        const std::string region_name, const size_t pheno_index,
                        const size_t region_index, Genotype& target)
{
    // prslice can easily be implemented using PRSet functionality
    // so maybe remove prslice from this function
    const bool no_regress = c_commander.no_regress();
    const bool all = c_commander.all();
    const int num_thread = c_commander.thread();
    Eigen::initParallel();
    Eigen::setNbThreads(num_thread);
    m_best_index = -1;
    m_num_snp_included = 0;
    m_perm_result.resize(m_num_perm, 2);
    m_prs_results.clear();
    m_current_sample_score.clear();
    m_best_sample_score.clear();
    m_prs_results.resize(target.num_threshold());
    // set to -1 to indicate not don
    for (auto&& p : m_prs_results) p.threshold = -1;
    const bool multi = pheno_info.name.size() > 1;
    // initialize score vector
    m_current_sample_score.resize(m_sample_included.size());
    // now initialize them
    for (size_t i_sample = 0; i_sample < m_sample_included.size(); ++i_sample) {
        m_current_sample_score[i_sample].has_pheno =
            m_sample_names[m_sample_index[i_sample]].has_pheno;
    }
    // directly copy it;
    m_best_sample_score = m_current_sample_score;

    // now prepare all score
    // in theory, we only need to calulate it once for every phenotype + sets
    // but it is easier to do it this way
    std::fstream all_out;
    size_t width_of_line = 0;
    size_t num_thresholds = 0;
    size_t header_length = 0;
    if (all) {
        std::vector<double> avail_thresholds = target.get_thresholds();
        std::sort(avail_thresholds.begin(), avail_thresholds.end());
        num_thresholds = avail_thresholds.size();
        // Most of the length below are hard coded. Not sure if
        // they will mess up in some crazy machine
        // i.e if those machine output numbers larger than
        // 12 digits...
        width_of_line = num_thresholds + num_thresholds * 12 + 1
                        + m_max_fid_length + m_max_iid_length;
        std::string header = "FID IID";
        for (auto& thres : avail_thresholds) {
            header.append(" " + std::to_string(thres));
        }
        header_length = header.length() + 1; // +1 for newline
        std::string all_out_name = c_commander.out();
        if (multi) {
            all_out_name.append("." + pheno_info.name[pheno_index]);
        }
        if (m_prset) all_out_name.append("." + region_name);
        all_out_name.append(".all.score");
        all_out.open(all_out_name.c_str(), std::fstream::out | std::fstream::in
                                               | std::fstream::trunc);

        if (!all_out.is_open()) {
            std::string error_message =
                "Cannot open file " + all_out_name + " for write";
            throw std::runtime_error(error_message);
        }
        all_out << header << std::endl;
        for (auto&& sample : m_sample_index) {
            std::string name =
                m_sample_names[sample].FID + " " + m_sample_names[sample].IID;
            all_out << std::setfill(' ') << std::setw(width_of_line)
                    << std::left << name << std::endl;
        }
        width_of_line++; // to account for the new line
    }


    // current threshold iteration
    size_t iter_threshold = 0;
    // +1 such that only 100% when finished
    size_t max_category = target.max_category() + 1;
    int cur_category = 0, cur_index = -1;
    double cur_threshold = 0.0, prev_progress = 0.0;
    while (target.get_score(m_current_sample_score, cur_index, cur_category,
                            cur_threshold, m_num_snp_included, region_index))
    {
        double progress =
            (double) cur_category / (double) (max_category) *100.0;
        if (progress - prev_progress > 0.01 && !m_prset) {
            fprintf(stderr, "\rProcessing %03.2f%%", progress);
            prev_progress = progress;
        }

        if (all) {
            for (size_t sample = 0; sample < m_sample_included.size(); ++sample)
            {
                double score =
                    (m_average_score)
                        ? ((m_current_sample_score[sample].num_snp == 0)
                               ? 0
                               : m_current_sample_score[sample].prs
                                     / (double) m_current_sample_score[sample]
                                           .num_snp)
                        : m_current_sample_score[sample].prs;
                size_t loc = header_length + sample * width_of_line
                             + m_max_fid_length + 1 + m_max_iid_length + 1
                             + iter_threshold + iter_threshold * 12;
                all_out.seekp(loc);
                all_out << score;
            }
        }
        if (no_regress) {
            iter_threshold++;
            continue;
        }
        regress_score(cur_threshold, num_thread, pheno_index, iter_threshold);

        if (c_commander.permute()) {
            permutation(num_thread, c_commander.logit_perm()
                                        && m_target_binary[pheno_index]);
        }
        iter_threshold++;
    }
    if (all_out.is_open()) all_out.close();
    if (!m_prset) fprintf(stderr, "\rProcessing %03.2f%%\n", 100.0);
    process_permutations();
}

void PRSice::process_permutations()
{
    // can't generate an empirical p-value if there is no observed p-value
    if (m_best_index == -1) return;
    double best_p = m_prs_results[m_best_index].p;
    size_t num_better = 0;
    for (auto&& p : m_perm_result) num_better += (p <= best_p);
    m_prs_results[m_best_index].emp_p =
        (double) (num_better + 1.0) / (double) (m_num_perm + 1.0);
}

void PRSice::permutation(const size_t n_thread, bool logit_perm)
{
    int num_iter = m_num_perm / m_perm_per_slice;
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(
        m_phenotype.rows());
    size_t num_include_samples = m_current_sample_score.size();
    Eigen::setNbThreads(n_thread);
    for (size_t sample_id = 0; sample_id < num_include_samples; ++sample_id) {
        std::string sample = (m_ignore_fid)
                                 ? m_sample_names[sample_id].IID
                                 : m_sample_names[sample_id].FID + "_"
                                       + m_sample_names[sample_id].IID;

        if (m_sample_with_phenotypes.find(sample)
            != m_sample_with_phenotypes.end())
        {
            m_independent_variables(m_sample_with_phenotypes.at(sample), 1) =
                (m_average_score)
                    ? ((m_current_sample_score[sample_id].num_snp == 0)
                           ? 0
                           : m_current_sample_score[sample_id].prs
                                 / (double) m_current_sample_score[sample_id]
                                       .num_snp)
                    : m_current_sample_score[sample_id].prs;
        }
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomposed;
    int rank = 0;
    Eigen::VectorXd pre_se;
    if (!logit_perm) {
        decomposed.compute(m_independent_variables);
        rank = decomposed.rank();
        Eigen::MatrixXd R = decomposed.matrixR()
                                .topLeftCorner(rank, rank)
                                .triangularView<Eigen::Upper>();
        pre_se = (R.transpose() * R).inverse().diagonal();
    }
    Eigen::setNbThreads(1);
    int cur_remain = m_remain_slice;
    // we reseed the random number generator in each iteration so that
    // we will always get the same pheontype permutation
    std::mt19937 rand_gen{m_seed};
    // need to do the permutation of phenotype without threading
    // so that we can reserve the sequence
    size_t processed = 0;
    for (int iter = 0; iter < num_iter + 1; ++iter) {
        size_t cur_perm = m_perm_per_slice;
        cur_perm += (cur_remain > 0) ? 1 : 0;
        if (cur_perm + processed > m_num_perm) {
            cur_perm = m_num_perm - processed;
        }
        cur_remain--;
        std::vector<Eigen::MatrixXd> perm_pheno(cur_perm);
        for (size_t p = 0; p < cur_perm; ++p) {
            perm.setIdentity();
            std::shuffle(perm.indices().data(),
                         perm.indices().data() + perm.indices().size(),
                         rand_gen);
            perm_pheno[p] = perm * m_phenotype; // permute columns
        }
        // now multithread it and get the corresponding p-values
        std::vector<std::thread> thread_store;
        int job_size = cur_perm / n_thread;
        int remain = cur_perm % n_thread;
        size_t start = 0;
        for (size_t i_thread = 0; i_thread < n_thread; ++i_thread) {
            size_t ending = start + job_size + (remain > 0);
            ending = (ending > cur_perm) ? cur_perm : ending;
            thread_store.push_back(
                std::thread(&PRSice::thread_perm, this, std::ref(decomposed),
                            std::ref(perm_pheno), start, ending, rank,
                            std::cref(pre_se), processed, logit_perm));
            start = ending;
            remain--;
        }
        for (auto&& thread : thread_store) thread.join();
        processed += cur_perm;
    }
}

void PRSice::thread_perm(
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& decomposed,
    std::vector<Eigen::MatrixXd>& pheno_perm, size_t start, size_t end,
    int rank, const Eigen::VectorXd& pre_se, size_t processed, bool logit_perm)
{

    bool intercept = true;
    size_t n = m_independent_variables.rows();
    for (size_t i = start; i < end; ++i) {
        double ori_p = m_perm_result[processed + i];
        double obs_p = 2.0; // for safety reason, make sure it is out bound
        if (logit_perm) {
            double r2, coefficient;
            Regression::glm(pheno_perm[i], m_independent_variables, obs_p, r2,
                            coefficient, 25, 1, true);
        }
        else
        {
            Eigen::VectorXd beta = decomposed.solve(pheno_perm[i]);
            Eigen::MatrixXd fitted = m_independent_variables * beta;
            Eigen::VectorXd residual = pheno_perm[i] - fitted;
            int rdf = n - rank;
            double rss = 0.0;
            for (size_t r = 0; r < n; ++r) {
                rss += residual(r) * residual(r);
            }
            size_t se_index = intercept;
            for (size_t ind = 0; ind < (size_t) beta.rows(); ++ind) {
                if (decomposed.colsPermutation().indices()(ind) == intercept) {
                    se_index = ind;
                    break;
                }
            }
            double resvar = rss / (double) rdf;
            Eigen::VectorXd se = (pre_se * resvar).array().sqrt();
            double tval = beta(intercept) / se(se_index);
            boost::math::students_t dist(rdf);
            obs_p =
                2 * boost::math::cdf(boost::math::complement(dist, fabs(tval)));
        }
        // store the best p_value for the processed+i permutaiton
        m_perm_result[processed + i] = (ori_p > obs_p) ? obs_p : ori_p;
    }
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

    const bool perm = c_commander.permute();
    std::string output_name = output_prefix;

    bool valid = m_best_index != -1;
    if (!valid
        || region.get_count(region_index)
               == 0) // we know regions with 0 SNP will not have valid PRS
    {
        if (region.get_count(region_index) != 0) {
            fprintf(stderr, "ERROR: No valid PRS ");
            if (m_prset)
                fprintf(stderr, "for %s",
                        region.get_name(region_index).c_str());
            fprintf(stderr, "!\n");
        }
        return;
    }
    if (m_prset)
        output_name = output_prefix + "." + region.get_name(region_index);
    std::string out_best = output_name + ".best";
    std::string out_prsice = output_name + ".prsice";
    std::string out_snp = output_name + ".snps";
    // std::string out_summary = output_name + ".summary";
    std::ofstream best_out, prsice_out, snp_out, summary_out;
    prsice_out.open(out_prsice.c_str());
    if (!prsice_out.is_open()) {
        std::string error_message =
            "ERROR: Cannot open file: " + out_prsice + " to write";
        throw std::runtime_error(error_message);
    }
    prsice_out << "Threshold\tR2\tP\tCoefficient\tNum_SNP";
    if (perm) prsice_out << "\tEmpirical_P";
    prsice_out << std::endl;
    for (size_t i = 0; i < m_prs_results.size(); ++i) {
        if (m_prs_results[i].threshold < 0 || m_prs_results[i].p < 0) continue;
        double full = m_prs_results[i].r2;
        double null = m_null_r2;
        if (has_prevalence) {
            full = top * full / (1 + bottom * full);
            null = top * null / (1 + bottom * null);
        }
        double r2 = full - null;
        prsice_out << m_prs_results[i].threshold << "\t" << r2 << "\t"
                   << m_prs_results[i].p << "\t" << m_prs_results[i].coefficient
                   << "\t" << m_prs_results[i].num_snp;
        if (perm)
            prsice_out << "\t"
                       << ((m_prs_results[i].emp_p >= 0.0)
                               ? std::to_string(m_prs_results[i].emp_p)
                               : "-");
        prsice_out << std::endl;
    }
    prsice_out.close();

    best_out.open(out_best.c_str());
    if (!best_out.is_open()) {
        std::string error_message =
            "ERROR: Cannot open file: " + out_best + " to write";
        throw std::runtime_error(error_message);
    }
    /*
    summary_out.open(out_summary.c_str());
    if (!summary_out.is_open()) {
        std::string error_message =
            "ERROR: Cannot open file: " + out_summary + " to write";
        throw std::runtime_error(error_message);
    }
    */
    auto&& best_info = m_prs_results[m_best_index];
    // summary_out << "Best Threshold:   " << best_info.threshold << std::endl;


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
    /*
    summary_out << "R2 of PRS only:   " << r2 << std::endl;
    summary_out << "R2 of full model: " << full << std::endl;
    summary_out << "Null R2:          " << null << std::endl;
    summary_out << "P-value:          " << best_info.p << std::endl;
    if (perm)
        summary_out << "Empirical P:      " << best_info.emp_p << std::endl;
    summary_out << "Coefficient:      " << best_info.coefficient << std::endl;
    summary_out << "Number of SNPs:   " << best_info.num_snp << std::endl;
    summary_out.close();
    */
    best_out << "FID\tIID\tPRS\tHas_Phenotype" << std::endl;
    int best_snp_size = best_info.num_snp;
    if (best_snp_size == 0) {
        fprintf(stderr, "ERROR: Best R2 obtained when no SNPs were included\n");
        fprintf(stderr, "       Cannot output the best PRS score\n");
    }
    else
    {
        for (size_t sample = 0; sample < m_sample_index.size(); ++sample) {
            // samples that are extracted are ignored
            // sample excluded will not be output here
            if (!m_sample_names[m_sample_index[sample]].included) continue;
            std::string has_pheno =
                (m_sample_names[m_sample_index[sample]].has_pheno) ? "Yes"
                                                                   : "No";
            best_out << m_sample_names[m_sample_index[sample]].FID << "\t"
                     << m_sample_names[m_sample_index[sample]].IID << "\t"
                     << m_best_sample_score[sample].prs / (double) best_snp_size
                     << "\t" << has_pheno << std::endl;
        }
    }
    best_out.close();

    if (c_commander.print_snp()) {
        target.print_snp(out_snp, m_prs_results[m_best_index].threshold,
                         region_index);
    }
}

void PRSice::summarize(const Commander& commander)
{
    bool prev_out;

    const bool perm = commander.permute();

    fprintf(stderr, "There are ");
    if (m_significant_store[0] != 0) {
        fprintf(stderr,
                "%zu region(s) with p-value > 0.1 (\033[1;31mnot "
                "significant\033[0m);\n",
                m_significant_store[0]);
        prev_out = true;
    }
    if (m_significant_store[1] != 0) {
        if (m_significant_store[2] == 0 && prev_out) {
            fprintf(stderr, "and ");
        }
        fprintf(stderr,
                "%zu region(s) with p-value between \n"
                "0.1 and 1e-5  (\033[1;31mmay not be significant\033[0m);\n ",
                m_significant_store[1]);
        prev_out = true;
    }
    if (m_significant_store[2] != 0) {
        if (prev_out) fprintf(stderr, " and ");
        fprintf(stderr, "%zu region(s) with p-value less than 1e-5\n",
                m_significant_store[2]);
    }
    if (!perm) {
        fprintf(stderr,
                "Please note that these results are inflated due to the\n");
        fprintf(stderr, "overfitting inherent in finding the best-fit\n");
        fprintf(stderr,
                "PRS (but it's still best to find the best-fit PRS!).\n\n");
        fprintf(stderr,
                "You can use the --perm option (see manual) to calculate\n");
        fprintf(stderr, "an empirical P-value.\n");
    }
    std::string out_name = commander.out() + ".summary";
    std::ofstream out;
    out.open(out_name.c_str());
    if (!out.is_open()) {
        std::string error_message =
            "ERROR: Cannot open file: " + out_name + " to write";
        throw std::runtime_error(error_message);
    }
    out << "Phenotype\tSet\tThreshold\tPRS.R2\tFull.R2\tNull."
           "R2\tPrevalence\tCoefficient\tP\tNum_SNP";
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
                << sum.prevalence << "\t" << sum.result.coefficient << "\t"
                << sum.result.p << "\t" << sum.result.num_snp << std::endl;
        }
        else
        {
            out << "\t" << sum.result.r2 - sum.r2_null << "\t" << sum.result.r2
                << "\t" << sum.r2_null << "\t-";
        }
        out << "\t" << sum.result.coefficient << "\t" << sum.result.p << "\t"
            << sum.result.num_snp;
        if (perm) out << "\t" << sum.result.emp_p;
        out << std::endl;
    }
    out.close();
}

PRSice::~PRSice()
{
    // dtor
}
