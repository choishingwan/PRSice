#include "prsice.hpp"

std::mutex PRSice::score_mutex;


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

    if (all && !prslice) 
    {
        all_out << "Threshold\tRegion";
        for (auto &&sample : m_sample_names)
            all_out << "\t" << sample.FID<< ":" << sample.IID;
        all_out << std::endl;
        all_out.close();
    }
    double null_r2_adjust = 0.0, null_p = 0.0, null_coeff = 0.0;
    // calculate the null r2
    int n_thread = c_commander.thread();
    if (m_independent_variables.cols() > 2) 
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
}

void PRSice::gen_pheno_vec(const std::string &pheno_file_name, const int pheno_index, bool regress)
{
    std::vector<double> pheno_store;
    bool binary = pheno_info.binary[pheno_index];
    int pheno_col_index = pheno_info.col[pheno_index];
    int max_num = 0;
    int num_case =0;
    int num_control =0;
    size_t invalid_pheno = 0;
    size_t num_not_found = 0;
    std::string line;

    if(pheno_info.use_pheno) // use phenotype file
    {
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
            if (token.size() < pheno_index + 1)
            {
                std::string error_message = "Malformed pheno file, should contain at least "
                        + std::to_string(pheno_index + 1) + " columns";
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
                            sample.included = false;
                        }
                    }
                    else
                    {
                        m_sample_with_phenotypes[id]=pheno_store.size();
                        pheno_store.push_back(misc::convert<double>(phenotype_info[id]));
                    }
                }catch(const std::runtime_error &error){
                    invalid_pheno++;
                    sample.included=false;
                }
            }
            else
            {
                sample.included =false;
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
                sample.included=false;
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
                        sample.included = false;
                    }
                }
                else
                {
                    m_sample_with_phenotypes[m_ignore_fid?sample.IID:sample.FID+"_"+sample.IID]=pheno_store.size();
                    pheno_store.push_back(misc::convert<double>(sample.pheno));
                }
            }catch(const std::runtime_error &error){
                invalid_pheno++;
                sample.included = false;
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
    if(error)
    {
        throw std::runtime_error("Mixed encoding! Both 0/1 and 1/2 encoding found!");
    }
    if (pheno_store.size() == 0) throw std::runtime_error("No phenotype presented");
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
        for (auto &&cov: c_cov_header)
        {
            if(included.find(cov) == included.end())// to avoid duplicated covariance headers
            {
                included.insert(cov);
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
}


void PRSice::prsice(const Commander &c_commander, const std::vector<std::string> &region_name,
        const size_t c_pheno_index, Genotype &target, bool prslice)
{
    // Let the Genotype class lead the way
    bool no_regress = c_commander.no_regress() && !prslice;
    bool all = c_commander.all() && !prslice;
    bool multi = pheno_info.name.size()>0;
    std::ofstream all_out;
    if(all)
    {
        std::string all_out_name = c_commander.out() +".";
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
    m_prs_results.clear();
    m_prs_results.resize(m_region_size);
    /** REMEMBER, WE WANT ALL PRS FOR ALL SAMPLES **/
    size_t total_sample_size = m_sample_names.size();
    // These are lite version. We can ignore the FID and IID because we
    // know they will always follow the sequence in m_sample_names
    // this will help us saving some memory spaces
    m_best_score.resize(m_region_size);
    m_current_score.resize(m_region_size);
    for(size_t i_region = 0; i_region < m_region_size; ++i_region)
    {
        m_best_score[i_region].resize(total_sample_size);
        m_current_score[i_region].resize(total_sample_size);
    }

    // now let Genotype class do the work
    size_t max_category = target.max_category()+1; // so that it won't be 100% until the very end
    int cur_category=0, cur_index =0;
    double cur_threshold =0.0;
    while(target.get_score(m_current_score, cur_index, cur_category, cur_threshold, m_num_snp_included))
    {
        if (!prslice)
            fprintf(stderr, "\rProcessing %03.2f%%", (double) cur_category / (double) (max_category) * 100.0);
        if (all && all_out.is_open()) {
            for (size_t i_region = 0; i_region < m_region_size; ++i_region)
            {
                all_out << cur_threshold << "\t" << region_name.at(i_region);
                for (size_t sample = 0; sample < total_sample_size; ++sample)
                {
                    all_out << "\t" << m_current_score[i_region][sample].prs / (double) m_current_score[i_region][sample].num_snp;
                }
                all_out << std::endl;
            }
        }
        if(no_regress) continue;

        if (n_thread == 1 || m_region_size == 1)
        {
            thread_score(0, m_region_size, cur_threshold, n_thread, c_pheno_index);
        }
        else
        {
            if (m_region_size < n_thread)
            {
                for (size_t i_region = 0; i_region < m_region_size; ++i_region)
                {
                    thread_store.push_back( std::thread(&PRSice::thread_score, this,
                        i_region, i_region + 1, cur_threshold,
                        1, c_pheno_index));
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
                        c_pheno_index));
                    start = ending;
                    remain--;
                }
            }
            // joining the threads
            for (auto &&thread : thread_store) thread.join();
            thread_store.clear();
        }
    }
    if (all_out.is_open()) all_out.close();
    if (!prslice) fprintf(stderr, "\rProcessing %03.2f%%\n", 100.0);

    if(c_commander.permute())
    {
        unsigned int seed = std::random_device()();
        if(c_commander.seeded()) seed = c_commander.seed();
        std::mt19937 rand_gen{seed};
        fprintf(stderr, "\nStart performing permutation with seed %u\n\n", seed);
        // need to perform permutation for the best score
        if (n_thread == 1 || m_region_size == 1)
        {
            thread_perm(c_commander.num_permutation(), 0, m_region_size, n_thread, c_pheno_index, rand_gen);
        }
        else
        {
            if (m_region_size < n_thread)
            {
                for (size_t i_region = 0; i_region < m_region_size; ++i_region)
                {
                    thread_store.push_back( std::thread(&PRSice::thread_perm, this,
                            c_commander.num_permutation(), i_region, i_region + 1, 1,
                            c_pheno_index, rand_gen));
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
                    thread_store.push_back( std::thread(&PRSice::thread_perm, this,
                            c_commander.num_permutation(), start, ending, 1,
                            c_pheno_index, rand_gen));
                    start = ending;
                    remain--;
                }
            }
            // joining the threads
            for (auto &&thread : thread_store) thread.join();
            thread_store.clear();
        }
    }
    fprintf(stderr, "Completed!\n");
}



void PRSice::thread_score(size_t region_start, size_t region_end,
        double threshold, size_t thread, const size_t c_pheno_index)
{

    Eigen::MatrixXd X;
    bool thread_safe = false;
    if (region_start == 0 && region_end == m_current_score.size())
        thread_safe = true;
    else
        X = m_independent_variables;
    double r2 = 0.0, r2_adjust = 0.0, p_value = 0.0, coefficient = 0.0;
    for (size_t iter = region_start; iter < region_end; ++iter)
    {
        // The m_prs size check is just so that the back will be valid
        // m_prs will only be empty for the first run
        if (m_num_snp_included[iter] == 0 ||
            (m_prs_results[iter].size() != 0 &&
                m_num_snp_included[iter] == m_prs_results[iter].back().num_snp)
        )  continue; // don't bother when there is no additional SNPs added

        for (size_t sample_id = 0; sample_id < m_current_score[iter].size(); ++sample_id)
        {
            std::string sample = (m_ignore_fid)? m_sample_names[sample_id].IID:
                m_sample_names[sample_id].FID+"_"+m_sample_names[sample_id].IID;
            // The reason why we need to update the m_sample_with_phenotypes matrix
            if (m_sample_with_phenotypes.find(sample) != m_sample_with_phenotypes.end())
            {
                if(thread_safe)
                {
                    m_independent_variables(m_sample_with_phenotypes.at(sample), 1) =
                        (m_current_score[iter][sample_id].num_snp==0)? 0.0 :
                        m_current_score[iter][sample_id].prs / (double) m_current_score[iter][sample_id].num_snp ;
                }
                else
                {
                    X(m_sample_with_phenotypes.at(sample), 1) =
                        (m_current_score[iter][sample_id].num_snp==0)? 0.0 :
                        m_current_score[iter][sample_id].prs / (double) m_current_score[iter][sample_id].num_snp ;
                }
            }
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


        // It this is the best r2, then we will add it
        size_t best_index = m_best_index[iter];
        if (m_prs_results[iter].empty() || m_prs_results[iter][best_index].r2 < r2)
        {
            m_best_index[iter] = m_prs_results[iter].size();
            m_best_score[iter] = m_current_score[iter];
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
        m_prs_results[iter].push_back(cur_result);
    }
}





void PRSice::thread_perm(const size_t num_perm,  size_t region_start, size_t region_end,
        size_t thread, const size_t c_pheno_index, std::mt19937 rand_gen)
{
    Eigen::MatrixXd X;
    bool thread_safe = false;
    if (region_start == 0 && region_end == m_current_score.size())
        thread_safe = true;
    else
        X = m_independent_variables;
    double r2 = 0.0, r2_adjust = 0.0, p_value = 0.0, coefficient = 0.0;
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm( m_phenotype.rows());
    Eigen::MatrixXd A_perm;
    for (size_t iter = region_start; iter < region_end; ++iter)
    {
        // first build the independent variable matrix
        for (size_t sample_id = 0; sample_id < m_current_score[iter].size(); ++sample_id)
        {
            std::string sample = (m_ignore_fid)? m_sample_names[sample_id].IID:
                    m_sample_names[sample_id].FID+"_"+m_sample_names[sample_id].IID;
            // The reason why we need to update the m_sample_with_phenotypes matrix
            if (m_sample_with_phenotypes.find(sample) != m_sample_with_phenotypes.end())
            {
                if(thread_safe)
                {
                    m_independent_variables(m_sample_with_phenotypes.at(sample), 1) =
                            (m_best_score[iter][sample_id].num_snp==0)? 0.0 :
                                    m_best_score[iter][sample_id].prs / (double) m_best_score[iter][sample_id].num_snp ;
                }
                else
                {
                    X(m_sample_with_phenotypes.at(sample), 1) =
                            (m_best_score[iter][sample_id].num_snp==0)? 0.0 :
                                    m_best_score[iter][sample_id].prs / (double) m_best_score[iter][sample_id].num_snp ;
                }
            }
        }
        size_t best = m_best_index[iter];
        double p_value = m_prs_results[iter][best].p;
        size_t num_better = 0;
        double perm_p, perm_r2,perm_r2_adj, perm_coefficient;
        for (size_t i_perm = 0; i_perm < num_perm; ++i_perm)
        {
            perm.setIdentity();
            std::shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size(),
                    rand_gen);
            A_perm = perm * m_phenotype; // permute columns
            if (m_target_binary[c_pheno_index])
            {
                if (thread_safe)
                    Regression::glm(A_perm, m_independent_variables, perm_p, perm_r2, perm_coefficient, 25, thread, true);
                else
                    Regression::glm(A_perm, X, perm_p, perm_r2, perm_coefficient, 25, thread, true);
            }
            else{
                if (thread_safe)
                    Regression::linear_regression(A_perm, m_independent_variables, perm_p, perm_r2,
                            perm_r2_adj, perm_coefficient, thread, true);
                else
                    Regression::linear_regression(A_perm, X, perm_p, perm_r2, perm_r2_adj, perm_coefficient, thread, true);
           }
            if (perm_p < p_value) num_better++;
        }
        m_prs_results[iter][best].emp_p = (double)(num_better+1.0) /(double)(num_perm+1.0);
    }
}

void PRSice::output(const Commander &c_commander, const Region &c_region,
        size_t pheno_index, Genotype &target) const
{
    // this is ugly, need to make it better
    std::string pheno_name = (pheno_info.name.size()>1)?pheno_info.name[pheno_index]:"";
    std::string output_prefix = c_commander.out();
    if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
    bool perm = c_commander.permute();
    std::string output_name = output_prefix;
    for(size_t i_region = 0; i_region < m_region_size; ++i_region)
    {
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
        for (auto &&prs : m_prs_results[i_region]) 
        {
            prsice_out <<prs.threshold << "\t"
                    << prs.r2 - m_null_r2 << "\t"
                    << prs.p << "\t"
                    << prs.coefficient << "\t"
                    << prs.num_snp;
            if (perm) prsice_out << "\t" << ((prs.emp_p>=0.0)? std::to_string(prs.emp_p) : "_");
            prsice_out << std::endl;
        }
        prsice_out.close();

        best_out.open(out_best.c_str());
        if (!best_out.is_open())
        {
            std::string error_message = "ERROR: Cannot open file: " + out_best + " to write";
            throw std::runtime_error(error_message);
        }
        best_out << "FID\tIID\tIncluded\tprs_" << m_prs_results[i_region][m_best_index[i_region]].threshold<< std::endl;
        int best_snp_size = m_prs_results[i_region][m_best_index[i_region]].num_snp;
        if (best_snp_size == 0) 
        {
            fprintf(stderr, "ERROR: Best R2 obtained when no SNPs were included\n");
            fprintf(stderr, "       Cannot output the best PRS score\n");
        } else 
        {
            for(size_t sample=0; sample<m_sample_names.size(); ++sample)
            {
                best_out << m_sample_names[sample].FID << "\t"
                        << m_sample_names[sample].IID << "\t"
                        << ((m_sample_names[sample].included)? "Y":"N") <<"\t"
                        << m_best_score[i_region][sample].prs/(double) best_snp_size
                        << std::endl;
            }
        }
        best_out.close();

        if(c_commander.print_snp())
        {
            target.print_snp(out_snp, m_prs_results[i_region][m_best_index[i_region]].threshold);
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
            region_out << c_region.get_name(i_region) << "\t" <<
                    m_prs_results[i_region][bi].threshold << "\t"
                    << m_prs_results[i_region][bi].r2 - m_null_r2 << "\t"
                    << m_prs_results[i_region][bi].coefficient<< "\t"
                    << m_prs_results[i_region][bi].p << "\t"
                    << m_prs_results[i_region][bi].num_snp;
            if(perm) region_out << "\t" << ((m_prs_results[i_region][bi].emp_p>=0.0)?
                    std::to_string((double) (m_prs_results[i_region][bi].emp_p)) : "-");
            region_out << std::endl;
            i_region++;
        }
        region_out.close();
    }
}

PRSice::~PRSice() {
    //dtor
}


