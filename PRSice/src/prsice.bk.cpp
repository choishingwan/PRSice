#include "prsice.hpp"

std::mutex PRSice::score_mutex;


void PRSice::individual_pheno_prs(const Commander &c_commander, const Region &c_region, const bool pre_run,
                                  pheno_storage &pheno_info, const bool multi)
{
    bool fastscore = c_commander.fastscore();
    bool no_regress =c_commander.no_regression();
    bool all = c_commander.all();
    bool target_binary = c_commander.target_is_binary();
    double bound_start =  (fastscore)? c_commander.get_bar_lower(): c_commander.get_lower();
    double bound_end = (fastscore)? c_commander.get_bar_upper():c_commander.get_upper();
    double bound_inter = c_commander.get_inter();
    std::string target = c_commander.get_target();
    std::string pheno_file = c_commander.get_pheno();
    std::string output_name = c_commander.get_out();
    std::ofstream all_out;

    if(all)
    {
        std::string all_out_name = output_name+"."+m_current_base;
        if(multi)
        {
            all_out_name.append("."+std::get<pheno_store::NAME>(pheno_info));
        }
        all_out_name.append(".all.score");
        all_out.open(all_out_name.c_str());
        if(!all_out.is_open())
        {
            std::string error_message = "Cannot open file "+all_out_name+" for write";
            throw std::runtime_error(error_message);
        }
    }
    std::unordered_map<std::string,size_t> sample_index;
    std::vector<prs_score> sample_prs;
    Eigen::VectorXd phenotype;
    Eigen::MatrixXd independent_variables;
    if(!no_regress)
    {
        gen_pheno_vec(phenotype, target, std::get<pheno_store::FILE_NAME>(pheno_info),
                      std::get<pheno_store::INDEX>(pheno_info), target_binary, sample_index, sample_prs);
        gen_cov_matrix(independent_variables, c_commander.get_cov_file(), c_commander.get_cov_header(),
                       sample_index);
    }
    else
    {
        std::string fam_name = target+".fam";
        if(fam_name.find("#")!=std::string::npos)
        {
            misc::replace_substring(fam_name, "#", m_chr_list.front());
        }
        std::ifstream fam;
        fam.open(fam_name.c_str());
        if(!fam.is_open())
        {
            std::string error_message = "ERROR: Cannot open fam file: " + fam_name;
            throw std::runtime_error(error_message);
        }
        std::string line;
        while(std::getline(fam, line))
        {
            misc::trim(line);
            if(!line.empty())
            {
                std::vector<std::string> token = misc::split(line);
                if(token.size() < 6) throw std::runtime_error("Malformed fam file, should contain at least 6 columns");
                sample_prs.push_back(prs_score(token[+FAM::IID], 0.0));
            }
        }
        fam.close();
    }
    for(size_t i_region=0; i_region < c_region.size(); ++i_region)
    {
        m_current_prs.push_back(sample_prs);
    }
    m_num_snp_included = std::vector<size_t>(c_region.size());
    size_t cur_start_index = 0;
    if(pre_run) get_prs_score(target, cur_start_index);
    m_best_score = m_current_prs;
    if(all)
    {
        all_out << "Threshold\tRegion";
        for(auto &&sample : sample_prs) all_out << "\t" << std::get<+PRS::IID>(sample); // Get the sample names
        if(pre_run)
        {
            for(size_t i_region=0; i_region < m_current_prs.size(); ++i_region)
            {
                for(auto &&prs : m_current_prs[i_region])
                {
                    all_out << "\t" << std::get<+PRS::PRS>(prs)/(double)m_num_snp_included[i_region];
                }
                all_out << std::endl;
            }
        }
        all_out.close();
    }
    double null_r2 = 0.0, null_r2_adjust=0.0, null_p=0.0;
    int n_thread= c_commander.get_thread();
    if(independent_variables.cols()>2)
    {
        Eigen::MatrixXd covariates_only;
        covariates_only = independent_variables;
        covariates_only.block(0,1,covariates_only.rows(),covariates_only.cols()-2) = covariates_only.topRightCorner(covariates_only.rows(),covariates_only.cols()-2);
        covariates_only.conservativeResize(covariates_only.rows(),covariates_only.cols()-1);
        if(target_binary)
        {
            Regression::glm(phenotype, covariates_only, null_p, null_r2, 25, n_thread, true);
        }
        else
        {
            Regression::linear_regression(phenotype, covariates_only, null_p, null_r2, null_r2_adjust, n_thread, true);
        }
    }
    calculate_scores(c_commander, c_region, cur_start_index,
                     independent_variables, phenotype, std::get<pheno_store::NAME>(pheno_info), sample_index);
    if(!no_regress)
    {
        prs_output(c_commander, c_region, null_r2, std::get<pheno_store::NAME>(pheno_info));
    }

}


void PRSice::pheno_prslice(const Commander &c_commander, const Region &c_region,  pheno_storage &pheno_index,
                           const bool multi)
{
    // get the null first
    m_prs_results.clear();
    m_best_threshold.clear();
    m_best_score.clear();
    m_num_snp_included.clear();
    m_current_prs.clear();
    bool target_binary = c_commander.target_is_binary();
    std::string target = c_commander.get_target();
    std::string pheno_file = c_commander.get_pheno();
    std::unordered_map<std::string,size_t> sample_index;
    std::unordered_map<size_t, std::string> bin_names;
    std::vector<prs_score> sample_prs;
    Eigen::VectorXd phenotype;
    Eigen::MatrixXd independent_variables;
    gen_pheno_vec(phenotype, target, std::get<pheno_store::FILE_NAME>(pheno_index),
                  std::get<pheno_store::INDEX>(pheno_index), target_binary, sample_index, sample_prs);
    gen_cov_matrix(independent_variables, c_commander.get_cov_file(), c_commander.get_cov_header(), sample_index);
    for(size_t i_region=0; i_region < c_region.size(); ++i_region)
    {
        m_current_prs.push_back(sample_prs);
    }
    m_num_snp_included = std::vector<size_t>(c_region.size());
    double null_r2 = 0.0, null_r2_adjust=0.0, null_p=0.0;
    int n_thread= c_commander.get_thread();
    if(independent_variables.cols()>2)
    {
        Eigen::MatrixXd covariates_only;
        covariates_only = independent_variables;
        covariates_only.block(0,1,covariates_only.rows(),covariates_only.cols()-2) = covariates_only.topRightCorner(covariates_only.rows(),covariates_only.cols()-2);
        covariates_only.conservativeResize(covariates_only.rows(),covariates_only.cols()-1);
        if(target_binary)
        {
            Regression::glm(phenotype, covariates_only, null_p, null_r2, 25, n_thread, true);
        }
        else
        {
            Regression::linear_regression(phenotype, covariates_only, null_p, null_r2, null_r2_adjust, n_thread, true);
        }
    }






    int window = c_commander.prslice();
    std::unordered_map<std::string, std::string> file_name_reference;
    std::vector<std::string> file_name;
    if(c_commander.get_target().find("#")!=std::string::npos)
    {
        for(auto &&chr: m_chr_list)
        {
            std::string name = c_commander.get_target();
            misc::replace_substring(name, "#", chr);
            file_name_reference[chr] = name;
            file_name.push_back(name);
        }
    }
    else
    {
        for(auto &&chr : m_chr_list)
        {
            file_name_reference[chr] = c_commander.get_target();
        }
        file_name.push_back(c_commander.get_target());
    }
    // Although in a sense, the partition should be the same for all PRSlice Phenotype run,
    // we will do it separately just to simplify the code (after all, I don't think it is
    // very efficient to run all phenotype together anyway)
    bool fastscore = c_commander.fastscore();
    double bound_start =  (fastscore)? c_commander.get_bar_lower(): c_commander.get_lower();
    double bound_end = (fastscore)? c_commander.get_bar_upper():c_commander.get_upper();
    double bound_inter = c_commander.get_inter();
    bool full_model = c_commander.full();
    std::string prev_chr = "";
    size_t prev_loc = 0;
    std::vector<double> best_r2;
    std::vector<std::vector<p_partition> > best_snp_index;
    std::unordered_map<std::string, size_t> part_ref;
    bool pre_run = false;
    std::string cur_name="";
    std::string last_chr="";
    size_t last_loc=0;
    size_t cur_bin_count = 0;
    for(auto &&snp : m_snp_list) // for each SNP
    {
        if(m_include_snp.find(snp.get_rs_id())!=m_include_snp.end())
        {
            std::string cur_chr = snp.get_chr();
            size_t cur_loc = snp.get_loc();
            last_chr = cur_chr;
            last_loc = cur_loc;
            if((prev_chr.empty() || cur_chr.compare(prev_chr)!=0) || (cur_loc-prev_loc) > window)
            {
                bin_names[cur_bin_count-1].append("-"+cur_loc);
                for(auto file : file_name)
                {
                    size_t cur_line = 0;
                    std::ifstream bim;
                    bim.open(file.c_str());
                    if(!bim.is_open())
                    {
                        std::string error_message = "Cannot open bim file: " + file;
                        throw std::runtime_error(error_message);
                    }
                    std::string line;
                    while(std::getline(bim, line))
                    {
                        misc::trim(line);
                        if(!line.empty())
                        {
                            std::vector<std::string> token = misc::split(line);
                            if(token.size() < 6) throw std::runtime_error("Malformed bim file. Should contain at least 6 column");
                            std::string rs = token[+BIM::RS];
                            if(part_ref.find(rs)!=part_ref.end())
                            {
                                std::get<+PRS::LINE>(m_partition[part_ref[rs]]) = cur_line;
                            }
                        }
                        cur_line++;
                    }
                }
                // Run PRS here
                std::vector<p_partition> cur_best_snp_index;
                double cur_r2 = calculate_prslice_prs(c_commander, c_region, cur_best_snp_index, pre_run, pheno_index, sample_prs, sample_index, phenotype, independent_variables);
                best_r2.push_back(cur_r2);
                best_snp_index.push_back(cur_best_snp_index);
                // Reset stuff here
                m_partition.clear();
                prev_chr = cur_chr;
                prev_loc = cur_loc;
                cur_name = cur_chr+":"+cur_loc;
                bin_names[bin_count] = cur_name;
                bin_count++;
                pre_run = false;
            }
            else
            {
                // add this to the m_partition
                p_partition part;
                std::get<+PRS::RS>(part) = snp.get_rs_id();
                std::get<+PRS::LINE>(part) = 0; // we don't know what line it is
                std::get<+PRS::INDEX>(part) = m_include_snp[snp.get_rs_id()];
                std::get<+PRS::FILENAME>(part) = file_name_reference[cur_chr];
                double p = snp.get_p_value();
                if(p< bound_start)
                {
                    std::get<+PRS::CATEGORY>(part) = -1;
                    m_partition.push_back(part);
                    part_ref[snp.get_rs_id()] = m_partition.size()-1;
                }
                else if(p<bound_end)
                {
                    int category = -1;
                    if(fastscore)
                    {
                        category = c_commander.get_category(p);
                        if(category ==-2)
                        {
                            throw std::runtime_error("Undefined category!");
                        }
                    }
                    else category = (int)((p-bound_start)/bound_inter);
                    if(category <0) pre_run=true;
                    std::get<+PRS::CATEGORY>(part) = category;
                    m_partition.push_back(part);
                    part_ref[snp.get_rs_id()] = m_partition.size()-1;
                }
                else if(full_model)
                {
                    std::get<+PRS::CATEGORY>(part) = (int)(1-bound_start)/bound_inter;
                    m_partition.push_back(part);
                    part_ref[snp.get_rs_id()] = m_partition.size()-1;
                }
            }
        }
    }
    // Finish the left overs
    if(m_partition.size() > 0)
    {
        bin_names[cur_bin_count-1].append("-"+last_loc);
        for(auto file : file_name)
        {
            size_t cur_line = 0;
            std::ifstream bim;
            bim.open(file.c_str());
            if(!bim.is_open())
            {
                std::string error_message = "Cannot open bim file: " + file;
                throw std::runtime_error(error_message);
            }
            std::string line;
            while(std::getline(bim, line))
            {
                misc::trim(line);
                if(!line.empty())
                {
                    std::vector<std::string> token = misc::split(line);
                    if(token.size() < 6) throw std::runtime_error("Malformed bim file. Should contain at least 6 column");
                    std::string rs = token[+BIM::RS];
                    if(part_ref.find(rs)!=part_ref.end())
                    {
                        std::get<+PRS::LINE>(m_partition[part_ref[rs]]) = cur_line;
                    }
                }
                cur_line++;
            }
        }
        std::vector<p_partition> cur_best_snp_index;
        double cur_r2=calculate_prslice_prs(c_commander, c_region, cur_best_snp_index, pre_run, pheno_index, sample_prs, sample_index, phenotype, independent_variables);
        best_r2.push_back(cur_r2);
        best_snp_index.push_back(cur_best_snp_index);
    }

    // Now completed the PRSlice part, need to again
    std::vector<size_t> best_order = sort_by_r2(best_r2);
    int total = best_order.size();
    int processed = 0;
    // need to figure out how to present the data


    size_t cur_start_index=0;
    size_t bin_count = 0;
    for(auto best : best_order)
    {
        fprintf(stderr, "\rProcessing %d out of %d regions\r", processed, total);
        if(best_snp_index[best].size()==0)
        {
            throw std::runtime_error("Best Snp should always contain at least 1 SNP!");
        }
        m_partition.insert(m_partition.end(), best_snp_index[best].begin(), best_snp_index[best].end());
        )
        // just read everything
        PLINK prs(target, m_chr_list);
        prs.initialize();
        prs.get_score(m_partition, m_snp_list, m_current_prs, cur_start_index, m_partition.size());
        // now we've got the PRS, proceed to calculation
        cur_start_index = m_partition.size();
        if(n_thread == 1 || m_current_prs.size()==1)
        {
            thread_score(independent_variables, phenotype, sample_index, 0, m_current_prs.size(), target_binary, (int)bin_count,n_thread);
        }
        else
        {
            // perform multi threading
            if(c_region.size() < n_thread)
            {
                for(size_t i_region = 0; i_region < c_region.size(); ++i_region)
                {
                    thread_store.push_back(std::thread(&PRSice::thread_score, this, std::ref(independent_variables), std::cref(phenotype), std::cref(sample_index), i_region, i_region+1, target_binary, (int)bin_count,1));
                }
            }
            else
            {
                int job_size = c_region.size()/n_thread;
                int remain = c_region.size()%n_thread;
                size_t start =0;
                for(size_t i_thread = 0; i_thread < n_thread; ++i_thread)
                {
                    size_t ending = start+job_size+(remain>0);
                    ending = (ending>c_region.size())? c_region.size(): ending;
                    thread_store.push_back(std::thread(&PRSice::thread_score, this,
                                                       std::ref(independent_variables), std::cref(phenotype), std::cref(sample_index), start, ending, target_binary, (int)bin_count,1));
                    start=ending;
                    remain--;
                }
            }
            // joining the threads
            for(auto &&thread : thread_store) thread.join();
        }
        bin_count++;
    }
    // now the m_best_threhsold should have the info and the prs_result should contain the
    fprintf(stderr, "\nCompleted\n");
}

double PRSice::calculate_prslice_prs(const Commander &commander, const Region &c_region, std::vector<p_partition> &best_snp_index, bool pre_run,  pheno_storage &pheno_index, std::vector<prs_score> &sample_prs,
                                     std::unordered_map<std::string,size_t> &sample_index, Eigen::VectorXd &phenotype, Eigen::MatrixXd &independent_variables)
{
    if(m_partition.size()==0) return 0.0;
    m_prs_results.clear();
    m_best_threshold.clear();
    m_best_score.clear();
    m_num_snp_included.clear();
    m_current_prs.clear();
    std::string target = c_commander.get_target();
    for(size_t i_region=0; i_region < c_region.size(); ++i_region)
    {
        m_current_prs.push_back(sample_prs);
    }
    m_num_snp_included = std::vector<size_t>(c_region.size());
    size_t cur_start_index = 0;
    if(pre_run) get_prs_score(target, cur_start_index);
    m_best_score = m_current_prs;
    calculate_scores(c_commander, c_region, cur_start_index,
                     independent_variables, phenotype, std::get<pheno_store::NAME>(pheno_index), sample_index);
    // now need to get the best threshold information here

    size_t n_snp = std::get<+PRS::NSNP>(m_best_threshold[i_region]);
    for(size_t i_snp = 0; i_snp < m_partition.size(); ++i_snp)
    {
        best_snp_index.push_back(m_partition[i_snp]);
    }
    return std::get<+PRS::R2>(m_best_threshold[i_region]);
}


std::vector<size_t> sort_by_r2(const std::vector<double> &r2) const
{
    std::vector<size_t> idx(r2.size());
    std::iota(idx.begin(), idx.end(),0);
    std::sort(idx.begin(), idx.end(), [&r2](size_t i1, size_t i2)
    {
        return r2[i1] < r2[i2];
    });
    return idx;
}

// basically update the score vector to contain the new polygenic score
PRSice::PRSice()
{
    m_base_index = 0;
}

PRSice::~PRSice()
{
    //dtor
}
