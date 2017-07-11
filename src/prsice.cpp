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

void PRSice::set_lee(double prevalence, double case_ratio, double &top, double &bottom) const
{
    double x = misc::qnorm(1-prevalence);
    double z = misc::dnorm(misc::qnorm(1-prevalence));
    double i = z / prevalence;
    double cc = prevalence*(1-prevalence)*prevalence*(1-prevalence)/(z*z*case_ratio*(1-case_ratio));
    double theta =i*((case_ratio - prevalence)/(1 - prevalence))*(i*((case_ratio-prevalence)/
            (1-prevalence))-x);
    double e=1-pow(case_ratio,(2*case_ratio))*pow((1 -case_ratio),(2 * (1 -case_ratio)));
    top = cc*e;
    bottom = cc*e*theta;

    // in theory, this can be better because
    // most variables can be reuse
}



void PRSice::pheno_check(const Commander &c_commander) 
{
    std::vector < std::string > pheno_header = c_commander.pheno_col();
    std::string pheno_file = c_commander.pheno_file();
    if (pheno_header.size() != 0 && pheno_file.empty()) 
    {
        throw std::runtime_error( "You must provide a phenotype file for multiple phenotype analysis");
    }
    if (pheno_file.empty()) 
    {
        pheno_info.use_pheno = false;
        pheno_info.binary.push_back(c_commander.is_binary(0));
    } 
    else 
    {
        std::ifstream pheno;
        pheno.open(pheno_file.c_str());
        if (!pheno.is_open()) 
        {
            std::string error_message = "Cannot open phenotype file: " + pheno_file;
            throw std::runtime_error(error_message);
        }
        std::string line;
        std::getline(pheno, line);
        if (line.empty()) 
        {
            throw std::runtime_error( "Cannot have empty header line for phenotype file!");
        }
        pheno.close();
        misc::trim(line);
        std::vector < std::string > col = misc::split(line);
        bool found = false;
        std::unordered_map<std::string, bool> dup_col;
        if (pheno_header.size() == 0)
        {
            // use the second column from the pheno file
            pheno_info.use_pheno = true;
            pheno_info.col.push_back(1+!m_ignore_fid);
            pheno_info.name.push_back("");
            pheno_info.order.push_back(0);
            pheno_info.binary.push_back(c_commander.is_binary(0));
        }
        else
        {
            for (size_t i_pheno = 0; i_pheno < pheno_header.size(); ++i_pheno) 
            {
                if (dup_col.find(pheno_header[i_pheno]) == dup_col.end()) 
                {
                    found = false;
                    dup_col[pheno_header[i_pheno]] = true;
                    // start from 1+!m_ignore_fid to skip the iid and fid part
                    for (size_t i_column = 1+!m_ignore_fid; i_column < col.size(); ++i_column) 
                    {
                        if (col[i_column].compare(pheno_header[i_pheno]) == 0) 
                        {
                            found = true;
                            pheno_info.use_pheno=true;
                            pheno_info.col.push_back(i_column);
                            pheno_info.name.push_back(pheno_header[i_pheno]);
                            pheno_info.order.push_back(i_pheno);
                            pheno_info.binary.push_back(c_commander.is_binary(i_pheno));
                            break;
                        }
                    }
                    if (!found) 
                    {
                        fprintf(stderr, "Phenotype: %s cannot be found in phenotype file\n",
                            pheno_header[i_pheno].c_str());
                    }
                }
            }
        }

    }
    size_t num_pheno= (pheno_info.use_pheno)? pheno_info.col.size() : 1;
    fprintf(stderr, "There are a total of %zu phenotype to process\n",num_pheno);
}

void PRSice::init_matrix(const Commander &c_commander, const size_t pheno_index, Genotype &target,
    const bool prslice)
{
    m_null_r2 = 0.0;
    m_phenotype = Eigen::VectorXd::Zero(0);
    m_independent_variables.resize(0,0);
    bool no_regress = c_commander.no_regress();
    bool all = c_commander.all();
    bool transpose = c_commander.transpose();
    std::string pheno_file = c_commander.pheno_file();
    std::string output_name = c_commander.out();

    std::ofstream all_out;
    bool multi = pheno_info.col.size()>1;
    if(all && !prslice)
    {
        std::string all_out_name = output_name;
        if(multi)
        {
            all_out_name.append("."+pheno_info.name[pheno_index]);
        }
        all_out_name.append(".all.score");
        all_out.open(all_out_name.c_str());
        if(!all_out.is_open())
        {
            std::string error_message = "Cannot open file: "+all_out_name+" for write";
            throw std::runtime_error(error_message);
        }
    }
    m_sample_names = target.sample_names();
    gen_pheno_vec(pheno_file, pheno_index, !no_regress);
    if (!no_regress)
    {
        gen_cov_matrix(c_commander.get_cov_file(), c_commander.get_cov_header());
    }


    double null_r2_adjust = 0.0, null_p = 0.0, null_coeff = 0.0;
    // calculate the null r2
    int n_thread = c_commander.thread();
    if (m_independent_variables.cols() > 2 && !no_regress)
    {
        Eigen::MatrixXd covariates_only;
        covariates_only = m_independent_variables;
        covariates_only.block(0, 1, covariates_only.rows(), covariates_only.cols() - 2) =
            covariates_only.topRightCorner(covariates_only.rows(), covariates_only.cols() - 2);
        covariates_only.conservativeResize(covariates_only.rows(),covariates_only.cols() - 1);
        if (c_commander.is_binary(pheno_index)) 
        {
            Regression::glm(m_phenotype, covariates_only, null_p, m_null_r2,
                null_coeff, 25, n_thread, true);
        } 
        else 
        {
            Regression::linear_regression(m_phenotype, covariates_only, null_p,
                m_null_r2, null_r2_adjust, null_coeff, n_thread, true);
        }
    }
    target.update_include(m_sample_names);
}

void PRSice::gen_pheno_vec(const std::string &pheno_file_name, const int pheno_index, bool regress)
{
    std::vector<double> pheno_store;
    bool binary = pheno_info.binary[pheno_index];

    int max_num = 0;
    int num_case =0;
    int num_control =0;
    size_t invalid_pheno = 0;
    size_t num_not_found = 0;
    std::string line;
    if(pheno_info.use_pheno) // use phenotype file
    {
        int pheno_col_index = pheno_info.col[pheno_index];
        std::ifstream pheno_file;
        pheno_file.open(pheno_file_name.c_str());
        if(!pheno_file.is_open())
        {
            std::string error_message = "Cannot open phenotype file: " + pheno_file_name;
            throw std::runtime_error(error_message);
        }

        std::unordered_map<std::string, std::string> phenotype_info;
        while (std::getline(pheno_file, line))
        {

            misc::trim(line);
            if (line.empty()) continue;
            std::vector < std::string > token = misc::split(line);
            if (token.size() <= pheno_index + 1 + !m_ignore_fid) // need to check the range
            {
                std::string error_message = "Malformed pheno file, should contain at least "
                        + std::to_string(pheno_index + 2 + !m_ignore_fid) + " columns\n"
                                "Have you use the --ignore-fid option?";
                throw std::runtime_error(error_message);
            }
            std::string id =(m_ignore_fid)? token[0]:token[0]+"_"+token[1];
            phenotype_info[id] = token[pheno_col_index];
        }
        pheno_file.close();
        // now go through the sample information
        for(auto &&sample : m_sample_names)
        {
            std::string id = (m_ignore_fid)? sample.IID : sample.FID+"_"+sample.IID;
            if(phenotype_info.find(id)!=phenotype_info.end() && sample.included)
            {
                try{
                    if(binary)
                    {
                        int temp = misc::convert<int>(phenotype_info[id]);
                        if(temp >=0 && temp <= 2)
                        {
                            m_sample_with_phenotypes[id]=pheno_store.size();
                            pheno_store.push_back( temp);
                            max_num = (temp>max_num)?temp:max_num;
                            num_case+=(temp==1);
                            num_control+=(temp==0);
                        }
                        else
                        {
                            invalid_pheno++;
                            if(regress) sample.included=false;
                            // we don't care about invalid phenotype when we are not performing regression
                        }
                    }
                    else
                    {
                        m_sample_with_phenotypes[id]=pheno_store.size();
                        pheno_store.push_back(misc::convert<double>(phenotype_info[id]));
                    }
                }catch(const std::runtime_error &error){
                    invalid_pheno++;
                    if(regress) sample.included=false;

                }
            }
            else
            {
                if(regress) sample.included=false;
                num_not_found++;
            }
        }
    }
    else
    {
        // directly extract it from the sample_name stuff
        size_t cur_index =0;
        for(auto &&sample: m_sample_names)
        {
            if(sample.pheno.compare("NA")==0 || !sample.included){
                if(regress) sample.included=false;
                continue;
            }
            try{
                if(binary)
                {
                    int temp = misc::convert<int>(sample.pheno);
                    if(temp >=0 && temp <= 2)
                    {
                        m_sample_with_phenotypes[m_ignore_fid?sample.IID:sample.FID+"_"+sample.IID]=pheno_store.size();
                        pheno_store.push_back(temp);
                        max_num = (temp>max_num)?temp:max_num;
                        num_case+=(temp==1);
                        num_control+=(temp==0);
                    }
                    else
                    {
                        invalid_pheno++;
                        if(regress) sample.included=false;
                    }
                }
                else
                {
                    m_sample_with_phenotypes[m_ignore_fid?sample.IID:sample.FID+"_"+sample.IID]=pheno_store.size();
                    pheno_store.push_back(misc::convert<double>(sample.pheno));
                }
            }catch(const std::runtime_error &error){
                invalid_pheno++;
                if(regress) sample.included=false;
            }
        }
    }
    if (num_not_found != 0)
    {
        fprintf(stderr, "Number of missing samples: %zu\n", num_not_found);
    }
    if(invalid_pheno!=0)
    {
        fprintf(stderr, "Number of invalid phenotyps: %zu\n", invalid_pheno);
    }
    bool error = false;
    if(max_num > 1 && binary)
    {
        num_case = 0;
        num_control = 0;
        size_t check = 0;
        for(auto &&pheno : pheno_store)
        {
            pheno--;
            if(pheno < 0){
            	error = true;
            }
            else (pheno==1)? num_case++: num_control++;
            check++;
        }
    }
    if(error && regress)
    {
        throw std::runtime_error("Mixed encoding! Both 0/1 and 1/2 encoding found!");
    }
    if (pheno_store.size() == 0 && regress) throw std::runtime_error("No phenotype presented");
    m_phenotype = Eigen::Map<Eigen::VectorXd>(pheno_store.data(), pheno_store.size());
    if (binary) 
    {
        fprintf(stderr, "Number of controls : %i\n", num_control);
        fprintf(stderr, "Number of cases : %i\n", num_case);
        if (regress) 
        {
            if (num_control == 0) throw std::runtime_error("There are no control samples");
            if (num_case == 0) throw std::runtime_error("There are no cases");
        }
    } 
    else 
    {
        fprintf(stderr, "Number of sample(s) with phenotype  : %zu\n", m_phenotype.rows());
    }
}


std::vector<size_t> PRSice::get_cov_index(const std::string &c_cov_file,
        const std::vector<std::string> &c_cov_header)
{
    std::vector<size_t> cov_index;
    std::ifstream cov;
    cov.open(c_cov_file.c_str());
    if (!cov.is_open())
    {
        std::string error_message = "ERROR: Cannot open covariate file: " + c_cov_file;
        throw std::runtime_error(error_message);
    }
    std::string line;
    size_t num_valid = 0;
    std::getline(cov, line);
    // obtain the header information of the covariate file
    if(line.empty()) throw std::runtime_error("First line of covariate file is empty!");
    std::vector < std::string > token = misc::split(line);
    if (c_cov_header.size() == 0)
    {
        // if no header is provided, we will use all the covariates included
        for (size_t i = 1+!m_ignore_fid; i < token.size(); ++i) cov_index.push_back(i); //FID, therefore+1
    }
    else
    {
        std::unordered_set<std::string> included;
        for (auto cov: c_cov_header)
        {
            if(cov.empty()) continue;
            if(included.find(cov) == included.end())// to avoid duplicated covariance headers
            {
                // got annoyed with the input of PC.1 PC.2 PC.3, do this automatic thingy to substitute them
                if(cov.at(0)=='@')
                {
                    cov.erase(0,1);
                    std::vector<std::string> open = misc::split(cov, "[");
                    std::vector<std::string> info;
                    std::vector<bool> list;
                    for(auto o : open)
                    {
                        if(o.find("]")!=std::string::npos)
                        {
                            std::vector<std::string> close =misc::split(o, "]");
                            // the first one will always be the list
                            info.push_back(close[0]);
                            list.push_back(true);
                            for(size_t cl = 1; cl< close.size(); ++cl)
                            {
                                info.push_back(close[cl]);
                                list.push_back(false);
                            }
                        }
                        else
                        {
                            info.push_back(o);
                            list.push_back(false);
                        }
                    }
                    std::vector<std::string> final_covariates;
                    for(size_t c=0; c < info.size(); ++c)
                    {
                        if(list[c])
                        {
                            std::vector<std::string> individual =misc::split(info[c], ".");
                            std::vector<int> numeric;
                            for(auto &&ind : individual)
                            {
                                if(ind.find("-")!=std::string::npos)
                                {
                                    std::vector<std::string> range = misc::split(ind, "-");
                                    if(range.size() != 2)
                                    {
                                        throw std::runtime_error("ERROR: Invalid range format, range must be in the form of start-end");
                                    }
                                    try{
                                        int start = misc::convert<int>(range[0]);
                                        int end = misc::convert<int>(range[1]);
                                        if(start > end){
                                            int temp = end;
                                            end = start;
                                            start = temp;
                                        }
                                        for(size_t s = start; s<=end; ++s)
                                        {
                                            numeric.push_back(s);
                                        }
                                    } catch (const std::runtime_error &error){
                                        std::string error_message = "ERROR: Invalid parameter: "+range[0]+" or "+range[1]+", only allow integer!";
                                        throw std::runtime_error(error_message);
                                    }
                                }
                                else
                                {
                                    try {
                                        int temp = misc::convert<int>( ind);
                                        numeric.push_back(temp);
                                    } catch (const std::runtime_error &error){
                                        std::string error_message = "ERROR: Invalid parameter: "+ind+", only allow integer!";
                                        throw std::runtime_error(error_message);
                                    }
                                }
                            }

                            // Now we have all the numeric parameters
                            if(final_covariates.empty())
                            {
                                for(auto n: numeric)
                                {
                                    final_covariates.push_back(std::to_string(n));
                                }
                            }
                            else
                            {
                                size_t cur_size = final_covariates.size();
                                for(size_t final=0; final < cur_size; ++final){
                                    std::string cur = final_covariates[final];
                                    final_covariates[final].append(std::to_string(numeric.front()));
                                    for(size_t s = 1; s < numeric.size(); ++s)
                                    {
                                        final_covariates.push_back(cur+std::to_string(numeric[s]));
                                    }
                                }
                            }
                        }
                        else
                        {
                            for(size_t final=0; final < final_covariates.size(); ++final)
                            {
                                final_covariates[final].append(info[c]);
                            }
                            if(final_covariates.empty()) final_covariates.push_back(info[c]);
                        }
                    }
                    for(auto res : final_covariates)
                    {
                        if(included.find(res)==included.end())
                        {
                            included.insert(res);
                        }
                    }
                }
                else included.insert(cov);
            }
        }
        // same, +1 when fid is include
        for (size_t i_header = 1+!m_ignore_fid; i_header < token.size(); ++i_header)
        {
            if (included.find(token[i_header]) != included.end())
            {
                cov_index.push_back(i_header);
            }
        }
    }
    std::sort(cov_index.begin(), cov_index.end());
    if(cov_index.size()==0)
    {
        throw std::runtime_error("ERROR: No valid covariates!");
    }
    return cov_index;
}

void PRSice::gen_cov_matrix(const std::string &c_cov_file,
        const std::vector<std::string> &c_cov_header)
{
    size_t num_sample = m_sample_with_phenotypes.size();
	if(c_cov_file.empty())
	{
		m_independent_variables = Eigen::MatrixXd::Ones(num_sample, 2);
		return;
	}
    std::vector<size_t> cov_index = get_cov_index(c_cov_file, c_cov_header);
    fprintf(stderr, "\nStart processing the covariates\n");
    fprintf(stderr, "==============================\n");
    std::vector < std::pair<std::string, size_t> > valid_sample_index;
    m_independent_variables = Eigen::MatrixXd::Ones(num_sample, cov_index.size() + 2);
    bool valid=true;
    size_t num_valid=0;
    std::ifstream cov;
    cov.open(c_cov_file.c_str());
    if (!cov.is_open())
    {
        std::string error_message = "ERROR: Cannot open covariate file: " + c_cov_file;
        throw std::runtime_error(error_message);
    }
    std::string line;
    std::getline(cov, line); // remove header
    int max_index = cov_index.back()+1;
    std::vector<int> missing_count(cov_index.size(),0);
    int total_sample =0;
    while(std::getline(cov, line))
    {
        misc::trim(line);
        if(line.empty()) continue;
        valid = true;
        std::vector<std::string> token = misc::split(line);
        if(token.size() < max_index)
        {
            std::string error_message = "ERROR: Malformed covariate file, should contain at least "
                    + std::to_string(max_index) + " column!";
            throw std::runtime_error(error_message);
        }
        std::string id =(m_ignore_fid)? token[0]:token[0]+"_"+token[1];
        if (m_sample_with_phenotypes.find(id) != m_sample_with_phenotypes.end())
        { // sample is found in the phenotype vector
            int index = m_sample_with_phenotypes[id];
            for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov)
            {
                try {
                    double temp = misc::convert<double>( token[cov_index[i_cov]]);
                    m_independent_variables(index, i_cov + 2) = temp; // + 2 because first line = intercept, second line = PRS
                } catch (const std::runtime_error &error) {
                    valid = false;
                    m_independent_variables(index, i_cov + 2) = 0;
                    missing_count[i_cov]++;
                }
            }
            if (valid)
            {
                valid_sample_index.push_back( std::pair<std::string, size_t>(id, index));
                num_valid++;
            }
        }
    }
    // now we need to handle the situation where there are a different number of samples

    if (valid_sample_index.size() != num_sample && num_sample != 0) {
        int removed = num_sample - valid_sample_index.size();
        fprintf(stderr, "Number of samples with invalid covariate: %d\n", removed);
        double portion = (double) removed / (double) num_sample;
        if(valid_sample_index.size() == 0)
        {
            // check which column is the main problem
            for(size_t miss = 0; miss < missing_count.size(); ++miss)
            {
                if(missing_count[miss] == num_sample)
                {
                    // we sorted the column index so we can't tell what the column name is
                    // useless we also store the head of the file (too troublesome)
                    fprintf(stderr, "Column %zu is invalid, please check it is of the correct format\n", miss);
                }
            }
            throw std::runtime_error("All samples removed due to missingness in covariate file!");
        }
        if (portion > 0.05)
        {
            fprintf(stderr, "WARNING! More than %03.2f%% of the samples were removed!\n", portion * 100);
            fprintf(stderr, "         Do check if your covariate file is correct\n");
        }
        std::sort(begin(valid_sample_index), end(valid_sample_index),
            [](std::pair<std::string, size_t> const &t1, std::pair<std::string, size_t> const &t2)
            {
                if(std::get<1>(t1)==std::get<1>(t2)) return std::get<0>(t1).compare(std::get<0>(t2)) < 0;
                else return std::get<1>(t1) < std::get<1>(t2);
            });

        // update the m_phenotype and m_independent
        m_sample_with_phenotypes.clear();
        for (size_t cur_index = 0; cur_index < valid_sample_index.size(); ++cur_index) 
        {
            std::string name = std::get < 0 > (valid_sample_index[cur_index]);
            m_sample_with_phenotypes[name]  = cur_index;
            size_t update_index = std::get < 1 > (valid_sample_index[cur_index]);
            if (update_index != cur_index) 
            {
                m_phenotype(cur_index, 0) = m_phenotype(update_index, 0);
                for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov) 
                {
                    m_independent_variables(cur_index, i_cov + 2) = m_independent_variables(update_index, i_cov + 2);
                }
            }
        }
        m_independent_variables.conservativeResize(valid_sample_index.size(),
                m_independent_variables.cols());
        m_phenotype.conservativeResize(valid_sample_index.size(), 1);

        fprintf(stderr, "\nFinal number of samples: %zu\n\n", valid_sample_index.size());
    }
    else
    {
        fprintf(stderr, "\nFinal number of samples: %zu\n\n", valid_sample_index.size());
    }

}


void PRSice::prsice(const Commander &c_commander, const std::vector<std::string> &region_name,
        const size_t c_pheno_index, Genotype &target, bool prslice)
{
    // Let the Genotype class lead the way
    bool no_regress = c_commander.no_regress() && !prslice;
    bool all = c_commander.all() && !prslice;
    bool transpose = c_commander.transpose();
    bool multi = pheno_info.name.size()>1;
    std::ofstream all_out;
    if(all)
    {
        std::string all_out_name = c_commander.out();
        if(multi)
        {
            all_out_name.append("."+pheno_info.name[c_pheno_index]);
        }
        all_out_name.append(".all.score");
        all_out.open(all_out_name.c_str(), std::ofstream::app);
        if (!all_out.is_open())
        {
            std::string error_message = "Cannot open file " + all_out_name + " for write";
            throw std::runtime_error(error_message);
        }
    }
    Eigen::initParallel();
    std::vector < std::thread > thread_store;
    size_t n_thread = c_commander.thread();
    m_best_index.clear();
    m_best_index.resize(m_region_size);
    m_num_snp_included.resize(m_region_size, 0);

    // previous +1 is because of the cur_category not initiailized correctly
    // this happens when the lowest threshold has no SNPs in it
    m_prs_results =  misc::vec2d<prsice_result>(m_region_size, target.num_threshold());
    for(size_t i_region = 0; i_region < m_region_size; ++i_region)
    {
        for(size_t i = 0; i < m_prs_results.cols(); ++i)
        {
            m_prs_results(i_region, i).threshold = -1;
        }
    }

    /** REMEMBER, WE WANT ALL PRS FOR ALL SAMPLES **/
    /** 1/7 CHANGE BEHAVIOUR, WE ONLY RETAIN THE SELECTED SAMPLES
     *  (THIS IS FOR SCORE READING, i.e MAF calculation)
     **/

    for(size_t i_sample=0; i_sample < m_sample_names.size(); ++i_sample)
    {
        auto &&sample = m_sample_names[i_sample];
        if(sample.included)
        {
            std::string id = (m_ignore_fid)? sample.IID: sample.FID+"_"+sample.IID;
            m_sample_included.push_back(id);
            m_sample_index.push_back(i_sample);
        }
    }
    if (all && !prslice && !transpose)
    {
        //we skip it when it is transposed so that the lines are always regularish?
        all_out << "Threshold\tRegion";
        for (auto &&sample : m_sample_index)
        {
            if(m_ignore_fid) all_out << "\t" << m_sample_names[sample].IID;
            else all_out << "\t" << m_sample_names[sample].FID << "_" << m_sample_names[sample].IID;
        }
        all_out << std::endl;
    }
    size_t num_included_samples = m_sample_included.size();
    // These are lite version. We can ignore the FID and IID because we
    // know they will always follow the sequence in m_sample_names
    // this will help us saving some memory spaces
    // by default, m_current_sample_score only contains samples that are included
    // so doesn't need the included field
    m_current_sample_score = misc::vec2d<Sample_lite>(m_region_size, num_included_samples);
    m_best_sample_score = misc::vec2d<Sample_lite>(m_region_size, num_included_samples);
    // now let Genotype class do the work
    size_t max_category = target.max_category()+1; // so that it won't be 100% until the very end
    int cur_category=0, cur_index =-1;
    double cur_threshold =0.0;
    unsigned int seed = std::random_device()(); // might need to comment out this for valgrind cerr
    if(c_commander.seeded()) seed = c_commander.seed();
    // seed need to be outside the loop so each iteration will return the same sequence
    // therefore the same permutation
    // we also want to know how many samples we can hold within 1gb ram
    int perm_per_slice = 0;
    int remain_slice = 0;
    if(c_commander.permute())
    {
        m_region_perm_result = misc::vec2d<double>(m_region_size, c_commander.num_permutation(), 2.0);
        // first check for ridiculously large sample size
        if(CHAR_BIT*num_included_samples >1000000000)
        {
            perm_per_slice = 1;
        }
        else{
            // in theory, most of the time, perm_per_slice should be
            // equal to c_commander.num_permutation();
            int sample_memory = CHAR_BIT*num_included_samples;
            perm_per_slice = 1000000000/sample_memory;
            perm_per_slice = (perm_per_slice > c_commander.num_permutation())?
                    c_commander.num_permutation() :
                    perm_per_slice;
            // Additional slice to keep
            remain_slice = c_commander.num_permutation()%perm_per_slice;

        }
        if(m_target_binary[c_pheno_index])
        {
            fprintf(stderr, "\nWARNING: To speed up the permutation, we perform\n");
            fprintf(stderr, "         linear regression instead of logistic\n");
            fprintf(stderr, "         regression within the permutation and uses\n");
            fprintf(stderr, "         the p-value to rank the thresholds. Our assumptions\n");
            fprintf(stderr, "         are as follow:\n");
            fprintf(stderr, "         1. Linear Regression & Logistic Regression produce\n");
            fprintf(stderr, "            similar p-values\n");
            if(!c_commander.logit_perm())
            {
                fprintf(stderr, "         2. P-value is correlated with R2\n\n");
                fprintf(stderr, "         If you must, you can run logistic regression instead\n");
                fprintf(stderr, "         by setting the --logit-perm flag\n\n");
            }
        }
    }

    size_t iter_threshold =0;
    size_t cur_process = 0;
    while(target.get_score(m_current_sample_score, cur_index, cur_category, cur_threshold, m_num_snp_included))
    {

        if (!prslice)
            fprintf(stderr, "\rProcessing %03.2f%%", (double) cur_category / (double) (max_category) * 100.0);

        if (all && all_out.is_open()) {
            for (size_t i_region = 0; i_region < m_region_size; ++i_region)
            {
                if(transpose)
                {
                    all_out << cur_threshold << "\t" << i_region;
                }
                else all_out << cur_threshold << "\t" << region_name.at(i_region);
                for (size_t sample = 0; sample < num_included_samples; ++sample)
                {
                    double score = (m_current_sample_score(i_region,sample).num_snp==0)? 0 :m_current_sample_score(i_region,sample).prs / (double) m_current_sample_score(i_region,sample).num_snp;
                    all_out << "\t" << score;
                }
                all_out << std::endl;
            }
        }
        if(no_regress) continue;

        if (n_thread == 1 || m_region_size == 1)
        {
            thread_score(0, m_region_size, cur_threshold, n_thread, c_pheno_index, iter_threshold);
        }
        else
        {
            if (m_region_size < n_thread)
            {
                for (size_t i_region = 0; i_region < m_region_size; ++i_region)
                {
                    thread_store.push_back( std::thread(&PRSice::thread_score, this,
                        i_region, i_region + 1, cur_threshold,
                        1, c_pheno_index, iter_threshold));
                }
            }
            else
            {
                int job_size = m_region_size / n_thread;
                int remain = m_region_size % n_thread;
                size_t start = 0;
                for (size_t i_thread = 0; i_thread < n_thread; ++i_thread)
                {
                    size_t ending = start + job_size + (remain > 0);
                    ending = (ending > m_region_size) ? m_region_size : ending;
                    thread_store.push_back( std::thread(&PRSice::thread_score, this, start,
                        ending, cur_threshold, 1,
                        c_pheno_index, iter_threshold));
                    start = ending;
                    remain--;
                }
            }
            // joining the threads
            for (auto &&thread : thread_store) thread.join();
            thread_store.clear();
        }
        if(c_commander.permute())
        {
            permutation(seed, perm_per_slice, remain_slice,
                    c_commander.num_permutation(), n_thread, c_commander.logit_perm());
        }
        iter_threshold++;
    }
    if (all_out.is_open()) all_out.close();
    if (!prslice) fprintf(stderr, "\rProcessing %03.2f%%\n", 100.0);
    process_permutations();
}

void PRSice::thread_score(size_t region_start, size_t region_end,
        double threshold, size_t thread, const size_t c_pheno_index,
        const size_t iter_threshold)
{

    Eigen::MatrixXd X;
    bool thread_safe = false;
    if (region_start == 0 && region_end == m_region_size)
        thread_safe = true;
    else
        X = m_independent_variables;
    double r2 = 0.0, r2_adjust = 0.0, p_value = 0.0, coefficient = 0.0;
    size_t num_include_samples = m_current_sample_score.cols();
    for (size_t iter = region_start; iter < region_end; ++iter)
    {
        // The m_prs size check is just so that the back will be valid
        // m_prs will only be empty for the first run
		if (m_num_snp_included[iter] == 0 ||
            (m_num_snp_included[iter] == m_prs_results(iter, iter_threshold).num_snp)
        )  continue; // don't bother when there is no additional SNPs added
		// Problem is, if we do sample selection, it is possible for that SNP to
		// have MAF of 0 because of the small resulting sample size

        double total = 0.0;
        for (size_t sample_id = 0; sample_id < num_include_samples; ++sample_id)
        {
            std::string sample =m_sample_included[sample_id];
            // The reason why we need to update the m_sample_with_phenotypes matrix
            if (m_sample_with_phenotypes.find(sample) != m_sample_with_phenotypes.end())
            {
                double score = (m_current_sample_score(iter,sample_id).num_snp==0)? 0.0 :
                        m_current_sample_score(iter,sample_id).prs / (double) m_current_sample_score(iter,sample_id).num_snp ;
                total+= score;
                if(thread_safe)
                {
                    m_independent_variables(m_sample_with_phenotypes.at(sample), 1) = score;

                }
                else
                {
                    X(m_sample_with_phenotypes.at(sample), 1) =score;
                }
            }

        }
        if(total==0.0)
        {
            continue;
        }
        if (m_target_binary[c_pheno_index])
        {
            try {
                if (thread_safe)
                    Regression::glm(m_phenotype, m_independent_variables, p_value, r2, coefficient, 25, thread, true);
                else
                    Regression::glm(m_phenotype, X, p_value, r2, coefficient, 25, thread, true);
            } catch (const std::runtime_error &error) {
                // This should only happen when the glm doesn't converge.
                // Let's hope that won't happen...
                fprintf(stderr, "ERROR: GLM model did not converge!\n");
                fprintf(stderr, "       Please send me the DEBUG files\n");
                std::ofstream debug;
                debug.open("DEBUG");
                if (thread_safe)
                    debug << m_independent_variables << std::endl;
                else
                    debug << X << std::endl;
                debug.close();
                debug.open("DEBUG.y");
                debug << m_phenotype << std::endl;
                debug.close();
                fprintf(stderr, "ERROR: %s\n", error.what());
                exit(-1);
            }
        } 
        else 
        {
            if (thread_safe)
                Regression::linear_regression(m_phenotype, m_independent_variables, p_value, r2, r2_adjust,
                        coefficient, thread, true);
            else
                Regression::linear_regression(m_phenotype, X, p_value, r2, r2_adjust, coefficient, thread, true);
        }


        // If this is the best r2, then we will add it
        size_t best_index = m_best_index[iter];
        if (iter_threshold==0 || m_prs_results(iter, best_index).r2 < r2)
        {
            m_best_index[iter] = iter_threshold;
            for(size_t s =0; s<m_current_sample_score.cols(); ++s)
            {
                m_best_sample_score(iter, s) = m_current_sample_score(iter, s);
            }
        }
        // This should be thread safe as each thread will only mind their own region
        // now add the PRS result to the vectors (hopefully won't be out off scope
        prsice_result cur_result;
        cur_result.threshold=threshold;
        cur_result.r2 = r2;
        cur_result.r2_adj = r2_adjust;
        cur_result.coefficient = coefficient;
        cur_result.p = p_value;
        cur_result.emp_p = -1.0;
        cur_result.num_snp = m_num_snp_included[iter];
        m_prs_results(iter, iter_threshold)=cur_result;
    }
}






void PRSice::process_permutations()
{
    for(size_t i_region=0; i_region < m_region_size; ++i_region)
    {
        double best_p = m_prs_results(i_region, m_best_index[i_region]).p;
        size_t num_perm = m_region_perm_result.cols();
        size_t num_better = 0;
        for(size_t i = 0; i < num_perm; ++i)
        {
            num_better += (m_region_perm_result(i_region, i) <= best_p);
        }
        m_prs_results(i_region, m_best_index[i_region]).emp_p = (double)(num_better+1.0)/(double)(num_perm+1.0);
    }
}

void PRSice::permutation(unsigned int seed, int perm_per_slice, int remain_slice,
        int total_permutation, int n_thread, bool logit_perm)
{


    int num_iter = total_permutation/perm_per_slice;
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm( m_phenotype.rows());
    // the sequence should be as follow
    // 1. For each region, perform decomposition (this ensures we have different permutation for all N perm,
    //    otherwise, the sequence will be the same for each "slice"
    // 2. For each slice, build all permutations
    // 3. For each slice, perform multi-thread calculation
    //
    size_t num_include_samples = m_current_sample_score.cols();
    for(size_t i_region=0; i_region < m_region_size; ++i_region)
    {
        Eigen::setNbThreads(n_thread);
        // get the decomposition
        for (size_t sample_id = 0; sample_id < num_include_samples; ++sample_id)
        {
            std::string sample = m_sample_included[sample_id];
            // The reason why we need to update the m_sample_with_phenotypes matrix
            if (m_sample_with_phenotypes.find(sample) != m_sample_with_phenotypes.end())
            {
                m_independent_variables(m_sample_with_phenotypes.at(sample), 1) =
                        (m_current_sample_score(i_region,sample_id).num_snp==0)? 0.0 :
                                m_current_sample_score(i_region,sample_id).prs / (double) m_current_sample_score(i_region,sample_id).num_snp ;
            }
        }
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomposed;
        int rank=0;
        Eigen::VectorXd pre_se;
        if(!logit_perm)
        {
            decomposed.compute(m_independent_variables);
            rank = decomposed.rank();
            Eigen::MatrixXd R = decomposed.matrixR().topLeftCorner(rank, rank).triangularView<Eigen::Upper>();
            pre_se = (R.transpose()*R).inverse().diagonal();
        }
        Eigen::setNbThreads(1);
        int cur_remain = remain_slice;
        // we reseed the random number generator in each iteration so that
        // we will always get the same pheontype permutation
        std::mt19937 rand_gen{seed};
        // need to do the permutation of phenotype without threading
        // so that we can reserve the sequence
        size_t processed = 0;
        for(int iter = 0; iter < num_iter+1; ++iter)
        {
            int cur_perm = perm_per_slice;
            cur_perm += (cur_remain>0)? 1:0;
            if(cur_perm+processed > total_permutation)
            {
                cur_perm =total_permutation - processed;
            }
            cur_remain--;
            std::vector<Eigen::MatrixXd> perm_pheno(cur_perm);
            for(size_t p = 0; p < cur_perm; ++p)
            {
                perm.setIdentity();
                std::shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size(),
                        rand_gen);
                perm_pheno[p] = perm * m_phenotype; // permute columns
            }
            // now multithread it and get the corresponding p-values
            std::vector<std::thread> thread_store;
            int job_size = cur_perm / n_thread;
            int remain = cur_perm % n_thread;
            size_t start = 0;
            for (size_t i_thread = 0; i_thread < n_thread; ++i_thread)
            {
                size_t ending = start + job_size + (remain > 0);
                ending = (ending > cur_perm) ? cur_perm : ending;
                thread_store.push_back( std::thread(&PRSice::thread_perm, this,
                        std::ref(decomposed), std::ref(perm_pheno), start,
                        ending, i_region, rank, std::cref(pre_se), processed,
                        logit_perm));
                start = ending;
                remain--;
            }
            for(auto &&thread : thread_store) thread.join();
            processed+=cur_perm;
        }
    }
}


void PRSice::thread_perm(Eigen::ColPivHouseholderQR<Eigen::MatrixXd> &decomposed,
        std::vector<Eigen::MatrixXd> &pheno_perm, size_t start, size_t end,
        size_t i_region, int rank, const Eigen::VectorXd &pre_se, size_t processed,
        bool logit_perm)
{

    bool intercept =true;
    int n = m_independent_variables.rows();
    for(size_t i = start; i < end; ++i)
    {
        double ori_p = m_region_perm_result(i_region, processed+i);
        double obs_p = 2.0; // for safety reason, make sure it is out bound
        if(logit_perm)
        {
            double r2, coefficient;
            Regression::glm(pheno_perm[i], m_independent_variables, obs_p, r2, coefficient, 25, 1, true);
        }
        else
        {
            Eigen::VectorXd beta = decomposed.solve(pheno_perm[i]);
            Eigen::MatrixXd fitted = m_independent_variables*beta;
            Eigen::VectorXd residual = pheno_perm[i]-fitted;
            int rdf = n-rank;
            double rss = 0.0;
            for(size_t r = 0; r < n; ++r){
                rss+=residual(r)*residual(r);
            }
            size_t se_index = intercept;
            for(size_t ind=0;ind < beta.rows(); ++ind)
            {
                if(decomposed.colsPermutation().indices()(ind) == intercept)
                {
                    se_index = ind;
                    break;
                }
            }
            double resvar = rss/(double)rdf;
            Eigen::VectorXd se = (pre_se*resvar).array().sqrt();
            double tval = beta(intercept)/se(se_index);
            boost::math::students_t dist(rdf);
            obs_p = 2*boost::math::cdf(boost::math::complement(dist, fabs(tval)));
        }
        m_region_perm_result(i_region, processed+i) = (ori_p>obs_p)? obs_p : ori_p;
    }
}

void PRSice::output(const Commander &c_commander, const Region &c_region,
        size_t pheno_index, Genotype &target) const
{
    // this is ugly, need to make it better
    // check if it needs to be adjusted
    std::vector<double> prev = c_commander.prevalence();
    bool has_prevalence = (prev.size()!=0);
    int num_binary = 0;
    for(size_t i = 0; i < pheno_index; ++i)
    {
        if(c_commander.is_binary(i)) num_binary++;
    }
    double top=0.0, bottom = 0.0;
    if(has_prevalence && c_commander.is_binary(pheno_index))
    {
        int num_case = 0, num_control=0;
        for(size_t i = 0; i < m_phenotype.rows(); ++i)
        {
            if(m_phenotype(i)==0) num_control++;
            else if(m_phenotype(i)==1) num_case++;
        }
        set_lee(prev[num_binary], (double)(num_case)/(double)(num_case+num_control),
                top, bottom);
        // try to get the number of case and control
    }
    has_prevalence = has_prevalence&&c_commander.is_binary(pheno_index);
    std::string pheno_name = (pheno_info.name.size()>1)?pheno_info.name[pheno_index]:"";
    std::string output_prefix = c_commander.out();
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    bool perm = c_commander.permute();
    std::string output_name = output_prefix;
    for(size_t i_region = 0; i_region < m_region_size; ++i_region)
    {
        // check number of valid results
        bool valid = false;
        for(size_t i = 0; i < m_prs_results.cols(); ++i)
        {
            if(m_prs_results(i_region, i).threshold>=0)
            {
                valid = true;
                break;
            }
        }
        if(!valid)
        {
            fprintf(stderr, "ERROR: No valid PRS!\n");
            continue;
        }
        if(m_region_size > 1) output_name = output_prefix+"."+c_region.get_name(i_region);
        std::string out_best = output_name + ".best";
        std::string out_prsice = output_name + ".prsice";
        std::string out_snp = output_name +".snps";
        std::ofstream best_out, prsice_out, snp_out;
        prsice_out.open(out_prsice.c_str());
        if (!prsice_out.is_open())
        {
            std::string error_message = "ERROR: Cannot open file: " + out_prsice + " to write";
            throw std::runtime_error(error_message);
        }
        prsice_out << "Threshold\tR2\tP\tCoefficient\tNum_SNP";
        if (perm) prsice_out << "\tEmpirical_P";
        prsice_out << std::endl;
        for(size_t i = 0; i < m_prs_results.cols(); ++i)
        {
            if(m_prs_results(i_region,i).threshold < 0)  continue;
            double r2 = m_prs_results(i_region, i).r2 - m_null_r2;
            r2 = ((has_prevalence)? lee_adjust(r2, top, bottom):r2 );
            prsice_out <<m_prs_results(i_region, i).threshold << "\t"
                    << r2 << "\t"
                    << m_prs_results(i_region, i).p << "\t"
                    << m_prs_results(i_region, i).coefficient << "\t"
                    << m_prs_results(i_region, i).num_snp;
            if (perm) prsice_out << "\t" << ((m_prs_results(i_region, i).emp_p>=0.0)? std::to_string(m_prs_results(i_region, i).emp_p) : "-");
            prsice_out << std::endl;
        }
        prsice_out.close();

        best_out.open(out_best.c_str());
        if (!best_out.is_open())
        {
            std::string error_message = "ERROR: Cannot open file: " + out_best + " to write";
            throw std::runtime_error(error_message);
        }
        best_out << "FID\tIID\tprs_" << m_prs_results(i_region, m_best_index[i_region]).threshold<< std::endl;
        int best_snp_size = m_prs_results(i_region, m_best_index[i_region]).num_snp;
        if (best_snp_size == 0) 
        {
            fprintf(stderr, "ERROR: Best R2 obtained when no SNPs were included\n");
            fprintf(stderr, "       Cannot output the best PRS score\n");
        } else 
        {
            for(size_t sample=0; sample<m_sample_included.size(); ++sample)
            {
                best_out << m_sample_names[m_sample_index[sample]].FID << "\t"
                        << m_sample_names[m_sample_index[sample]].IID << "\t"
                        << m_best_sample_score(i_region,sample).prs/(double) best_snp_size
                        << std::endl;
            }
        }
        best_out.close();

        if(c_commander.print_snp())
        {
            target.print_snp(out_snp, m_prs_results(i_region, m_best_index[i_region]).threshold);
        }
        double best_p = m_prs_results(i_region, m_best_index[i_region]).p;
        if(!c_commander.permute())
        {
            fprintf(stderr, "\n");
            if(m_region_size>1)
                fprintf(stderr, "For %s\n", c_region.get_name(i_region).c_str());
            if(best_p > 0.1)
            {
                fprintf(stderr, "Your best-fit PRS has a P-Value > 0.1 which is \033[1;31mnot significant\033[0m.\n\n");
                fprintf(stderr, "This result is inflated due to the overfitting \n");
                fprintf(stderr, "inherent in finding the best-fit PRS (but it's still\n");
                fprintf(stderr, "best to find the best-fit PRS!).\n\n");
                fprintf(stderr, "You can use the --perm option (see manual) to calculate\n");
                fprintf(stderr, "an empirical P-value.\n");
            }
            else if(best_p > 1e-5)
            {
                fprintf(stderr, "Your best-fit PRS has a P-value between 0.1 and 1e-5\n");
                fprintf(stderr, "which \033[1;31mmay not be significant\033[0m - after accounting for\n");
                fprintf(stderr, "inflation due to the overfitting inherent in finding\n");
                fprintf(stderr, "the best-fit PRS (but it's still best to find the\n");
                fprintf(stderr, "best-fit PRS!).\n\n");
                fprintf(stderr, "While we suggest a P-value significance threshold of\n");
                fprintf(stderr, "P < 1e-4 for high resolution scoring (Euesden et al. 15),\n");
                fprintf(stderr, "we highly recommend repeating the analysis using the\n");
                fprintf(stderr, "--perm option (see manual) to calculate an empirical\n");
                fprintf(stderr, "P-value here.\n");
            }
            else
            {
                fprintf(stderr, "Your best-fit PRS has a P-value <= 1e-5 which is \033[1;31msignificant\033[0m\n");
                fprintf(stderr, "- according to our suggested P-value significance\n");
                fprintf(stderr, "threshold of P < 1e-4 for high resolution scoring\n");
                fprintf(stderr, "(Euesden et al. 15).\n\n");
                fprintf(stderr, "However, this result is inflated due to the\n");
                fprintf(stderr, "overfitting inherent in finding the best-fit\n");
                fprintf(stderr, "PRS (but it's still best to find the best-fit PRS!).\n\n");
                fprintf(stderr, "You can use the --perm option (see manual) to calculate\n");
                fprintf(stderr, "an empirical P-value.\n");
            }
        }
        //if(!c_commander.print_all()) break;
    }

    if(m_region_size > 1)
    {
        // now print the group information
        std::string out_region = output_prefix + ".prset";
        std::ofstream region_out;
        region_out.open(out_region.c_str());
        region_out << "Region\tThreshold\tR2\tCoefficient\tP\tNum_SNP";
        if (perm) region_out << "\tEmpirical_P";
        region_out    << std::endl;
        size_t i_region = 0;
        for (auto bi : m_best_index)
        {
            double r2 =m_prs_results(i_region, bi).r2 - m_null_r2;
            r2 = ((has_prevalence)? lee_adjust(r2, top, bottom):r2 );
            region_out << c_region.get_name(i_region) << "\t" <<
                    m_prs_results(i_region, bi).threshold << "\t"
                    << r2 << "\t"
                    << m_prs_results(i_region, bi).coefficient<< "\t"
                    << m_prs_results(i_region, bi).p << "\t"
                    << m_prs_results(i_region, bi).num_snp;
            if(perm) region_out << "\t" << ((m_prs_results(i_region, bi).emp_p>=0.0)?
                    std::to_string((double) (m_prs_results(i_region, bi).emp_p)) : "-");
            region_out << std::endl;
            i_region++;
        }
        region_out.close();
    }
}

void PRSice::transpose_all(const Commander &c_commander, const Region &c_region, size_t pheno_index) const
{
    // man... this will be so slow...
    // the whole reason why we don't like the transposed feature is that it will
    // be extremely slow and we thought it won't be helpful especially in the case
    // where no fastscore is used.
    // but nontheless we will let this function fly
    // the output is regular:
    // for each thresold, go through each regions
    fprintf(stderr, "\nTransposing all score file. Might take ages.\n");
    bool multi = pheno_info.name.size()>1;
    size_t num_samples = m_sample_included.size();
    // for any line, there will be
    std::string header = "FID\tIID";
    std::string output_name = c_commander.out();
    std::string all_out_name = output_name;
    if(multi)
    {
        all_out_name.append("."+pheno_info.name[pheno_index]);
    }
    all_out_name.append(".all.score");
    std::ifstream all_out;
    all_out.open(all_out_name.c_str());
    if(!all_out.is_open())
    {
        std::string error_message = "Cannot open file: "+all_out_name+" for write";
        throw std::runtime_error(error_message);
    }

    std::string line;
    std::string prev = "";
    while(std::getline(all_out, line)) // we have removed the header for this output
    {
        misc::trim(line);
        if(line.empty()) continue;
        std::string thres  = misc::get_column(line, 1);
        if(thres.empty() || thres.compare(prev)==0) continue;
        header.append("\t"+thres);
    }
    all_out.clear();
    all_out.seekg(0, std::ios::beg);
    // now generate the files
    std::vector<std::string> file_names;
    if(m_region_size==1)
    {
        std::string file_name = all_out_name;
        file_name.append(".transposed");
        std::ofstream out;
        out.open(file_name.c_str());
        if(!out.is_open())
        {
            std::string error_message = "Cannot open file: "+file_name+" for write";
            throw std::runtime_error(error_message);
        }
        out << header << std::endl;
        out.close();
        file_names.push_back(file_name);
    }
    else
    {
        std::string temp_name = output_name;
        if(multi)
        {
            temp_name.append("."+pheno_info.name[pheno_index]);
        }
        for(size_t i_region=0; i_region < m_region_size; ++i_region)
        {
            std::string file_name = temp_name+"."+c_region.get_name(i_region)+".all.score.transposed";
            std::ofstream out;
            out.open(file_name.c_str());
            if(!out.is_open())
            {
                std::string error_message = "Cannot open file: "+file_name+" for write";
                throw std::runtime_error(error_message);
            }
            out << header << std::endl;
            out.close();
            file_names.push_back(file_name);
        }
    }
    for(size_t i_sample=0; i_sample < m_sample_included.size(); ++i_sample)
    {
        std::vector<std::string> region_lines(m_region_size); // we append to the lines first then output at once
        while(std::getline(all_out, line)) // we have removed the header for this output
        {
            misc::trim(line);
            if(line.empty()) continue;
            std::string prs  = misc::get_column(line, 3+i_sample);
            int region  = misc::convert<int>(misc::get_column(line, 2));
            if(!prs.empty())
            {
                region_lines[region].append("\t"+prs);
            }
        }
        all_out.clear();
        all_out.seekg(0, std::ios::beg);
        // now output the lines to the file
        for(size_t f =0; f < file_names.size(); ++f)
        {
            std::ofstream fo;
            fo.open(file_names[f].c_str(), std::ofstream::app);
            if(!fo.is_open())
            {
                std::string error_message = "Cannot open file: "+file_names[f]+" for write";
                throw std::runtime_error(error_message);
            }
            fo << m_sample_names[m_sample_index[i_sample]].FID << "\t" <<
                    m_sample_names[m_sample_index[i_sample]].IID << region_lines[f] << std::endl;
            fo.close();
        }
    }
    all_out.close();
    std::remove( all_out_name.c_str());
}

PRSice::~PRSice() {
    //dtor
}


