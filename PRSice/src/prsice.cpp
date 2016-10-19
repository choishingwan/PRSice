#include "prsice.hpp"

std::mutex PRSice::score_mutex;


// Seems to be working alright (hope I know how to write test cases...)
void PRSice::get_snp(const Commander &c_commander, Region &region, const double &threshold)
{
    const std::string input = c_commander.get_base(m_base_index);
    const bool beta = c_commander.get_base_binary(m_base_index);
    // just issue the warning. would not terminate
    if(beta && c_commander.statistic().compare("OR")==0)
    {
        fprintf(stderr, "WARNING: OR detected but user suggest the input is beta!\n");
    }
    // First, we need to obtain the index of different columns from the base files
    // NOTE: -1 means missing and index is hard coded such that each index will represent
    //       specific header
    std::vector<int> index = SNP::get_index(c_commander, input);
    if((c_commander.get_target().find("#")!=std::string::npos && index[+SNP_Index::CHR]<0)||
            (c_commander.ld_prefix().find("#")!=std::string::npos && index[+SNP_Index::CHR]<0 ) )
    {
        std::string error_message = "To use chromosome separated PLINK input, you must provide"
                                    " the CHR header as we use the CHR information form the base file to substitute #";
        throw std::runtime_error(error_message);
    }
    else if(region.size() > 1 && (index[+SNP_Index::CHR] <0 || index[+SNP_Index::BP] <0 ))
    {
        std::string error_message = "To perform PRSet, you must provide the CHR and LOC header such"
                                    " that we can determine the set membership";
        throw std::runtime_error(error_message);
    }
    else if(c_commander.prslice() > 0.0 && (index[+SNP_Index::CHR] < 0 || index[+SNP_Index::BP]<0))
    {
        std::string error_message = "To perform PRSlice, you must provide the CHR and LOC header such"
                                    " that we can perform the slicing";
        throw std::runtime_error(error_message);
    }
    // Open the file
    std::ifstream snp_file;
    snp_file.open(input.c_str());
    if(!snp_file.is_open())
    {
        std::string error_message = "Cannot open base file: "+input;
        throw std::runtime_error(error_message);
    }
    // Some QC counts
    int num_duplicated = 0;
    size_t num_stat_not_convertible = 0;
    size_t num_p_not_convertible = 0;
    size_t num_indel = 0;
    size_t num_se_not_convertible=0;
    size_t num_exclude = 0;
    int max_index = index[+SNP_Index::MAX];
    std::string line;
    // remove header if index is not provided
    if(!c_commander.index()) std::getline(snp_file, line);
    bool read_error = false;
    bool not_converted = false;
    bool exclude=false;
    std::vector<std::string> token;
    // Actual reading the file, will do a bunch of QC
    while(std::getline(snp_file, line))
    {
        misc::trim(line);
        if(!line.empty())
        {
            not_converted=false;
            exclude = false;
            token = misc::split(line);
            if(token.size() <= max_index) throw std::runtime_error("More index than column in data");
            else
            {
                std::string rs_id = token[index[+SNP_Index::RS]];
                std::string chr = "";
                if(index[0] >= 0)
                {
                    chr = token[index[+SNP_Index::CHR]];
                }
                std::string ref_allele = (index[+SNP_Index::REF] >= 0)? token[index[+SNP_Index::REF]]:"";
                std::string alt_allele = (index[+SNP_Index::ALT] >= 0)? token[index[+SNP_Index::ALT]]:"";
                double pvalue = 0.0;
                if(index[+SNP_Index::P] >= 0)
                {
                    try
                    {
                        pvalue = misc::convert<double>(token[index[+SNP_Index::P]]);
                        if(pvalue < 0.0 || pvalue > 1.0)
                        {
                            read_error =true;
                            fprintf(stderr, "ERROR: %s's p-value is %f\n", rs_id.c_str(), pvalue);
                        }
                        else if(pvalue > threshold)
                        {
                            exclude=true;
                            num_exclude++;
                        }
                    }
                    catch(const std::runtime_error &error)
                    {
                        num_p_not_convertible++;
                        not_converted = true;
                    }
                }
                double stat = 0.0;
                if(index[+SNP_Index::STAT] >= 0)
                {
                    //Check if it is double
                    try
                    {
                        stat = misc::convert<double>(token[index[+SNP_Index::STAT]]);
                        if(!beta) stat = log(stat);
                    }
                    catch(const std::runtime_error& error)  //we know only runtime error is throw
                    {
                        num_stat_not_convertible++;
                        not_converted = true;
                    }
                }
                double se = 0.0;
                if(index[+SNP_Index::SE] >= 0)
                {
                    try
                    {
                        se = misc::convert<double>(token[index[+SNP_Index::SE]]);
                    }
                    catch(const std::runtime_error &error)
                    {
                        num_se_not_convertible++;
                    }
                }
                int loc = -1;
                if(index[+SNP_Index::BP]>=0)
                {
                    try
                    {
                        int temp = misc::convert<int>(token[index[+SNP_Index::BP]].c_str());
                        if(temp <0)
                        {
                            read_error=true;
                            fprintf(stderr, "ERROR: %s has negative loci\n", rs_id.c_str());
                        }
                        else loc = temp;
                    }
                    catch(const std::runtime_error &error)
                    {

                    }
                }
                if(ref_allele.compare("-")==0 || ref_allele.compare("I") == 0 || ref_allele.compare("D")==0 ||
                        ref_allele.size()>1)
                {
                    num_indel++;
                }
                else if(!alt_allele.empty() &&
                        (alt_allele.compare("-")==0 || alt_allele.compare("I") == 0 ||
                         alt_allele.compare("D")==0 || alt_allele.size()>1))
                {
                    num_indel++;
                }
                else if(! not_converted && ! exclude)
                {
                    m_snp_list.push_back(new SNP(rs_id, chr, loc, ref_allele, alt_allele, stat, se, pvalue, region.empty_flag(), region.size()));
                }
                else
                {
//		    			We skip any SNPs with non-convertible stat and p-value as we don't know how to
//		    			handle them. Most likely those will be NA, which should be ignored anyway
                }
            }
        }
        if(read_error) throw std::runtime_error("Please check if you have the correct input");
    }
    snp_file.close();
    // This may seems unnecessary and time consuming
    // but hopefully this should speed up the region discovery step
    // By performing sorting, we can ensure that the region and SNP
    // is in the same order
    // The only exception is when loc and chr are not provided. In
    // that case, SNPs will be sorted by the rsid and PRSet/PRSlice
    // cannot be run

    m_snp_list.sort();
    // Interestingly, because we already need the sorting of the SNP list, the map method of
    // finding duplication is actually slower (can be as much as 2x in mock test)
    // That is mainly due to the number of SNPs in the map
    // The larger the map, the slower it gets
    size_t before = m_snp_list.size();
    m_snp_list.erase(std::unique(m_snp_list.begin(), m_snp_list.end()), m_snp_list.end());
    size_t after = m_snp_list.size();
    // now write in the index
    // As we don't expect too many chromosome, map will be rather efficient here
    std::map<std::string, bool> unique_chr;
    for(size_t i_snp = 0; i_snp < m_snp_list.size(); ++i_snp)
    {
        m_snp_index[m_snp_list[i_snp].get_rs_id()]=i_snp;
        if(unique_chr.find(m_snp_list[i_snp].get_chr())==unique_chr.end())
        {
            unique_chr[m_snp_list[i_snp].get_chr()] = true;
            // The chr_list should be sorted in the same order as the SNPs' order
            m_chr_list.push_back(m_snp_list[i_snp].get_chr());
        }
        if(index[0] >=0 && index[5]>=0)
        {
            m_snp_list[i_snp].set_flag(region.check(m_snp_list[i_snp].get_chr(), m_snp_list[i_snp].get_loc()));
        }
    }

    // Now output the statistics. Might want to improve the outputs
    num_duplicated = (int)before-(int)after;
    fprintf(stderr, "Number of SNPs from base  : %zu\n", m_snp_list.size());
    if(num_indel!=0) fprintf(stderr, "Number of Indels          : %zu\n", num_indel);
    if(num_duplicated!=0) fprintf(stderr, "Number of duplicated SNPs : %d\n", num_duplicated);
    if(num_stat_not_convertible!=0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_convertible);
    if(num_p_not_convertible!=0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_convertible);
    if(num_se_not_convertible!=0) fprintf(stderr, "Failed to convert %zu SE\n", num_se_not_convertible);
}

void PRSice::clump(const Commander &c_commander, const Region &region)
{
    std::string target = c_commander.get_target();
    std::unordered_map<std::string, size_t> inclusion;
    bool has_ld = !c_commander.ld_prefix().empty();
    std::string ld_file = (has_ld)? c_commander.ld_prefix(): target;
    PLINK clump(ld_file, m_chr_list, c_commander.get_thread());
    // Because we will go through the target anyway, first go through the target for inclusion
    size_t num_ambig=0, not_found=0, num_duplicate=0;
    std::string target_bim_name = target+".bim";
    if(target_bim_name.find("#")!=std::string::npos)
    {
        for(auto chr: m_chr_list)
        {
            std::string target_chr_bim_name = target_bim_name;
            misc::replace_substring(target_chr_bim_name, "#", chr);
            check_inclusion(inclusion, target_chr_bim_name, num_ambig, not_found, num_duplicate);
        }
    }
    else check_inclusion(inclusion, target_bim_name, num_ambig, not_found, num_duplicate);
    fprintf(stderr, "\nIn Target File\n");
    fprintf(stderr,"==============================\n");
    if(num_ambig != 0)	fprintf(stderr, "Number of ambiguous SNPs  : %zu\n", num_ambig);
    if(num_duplicate != 0) fprintf(stderr, "Number of duplicated SNPs : %zu\n", num_duplicate);
    if(not_found != 0) fprintf(stderr, "Number of SNPs not found  : %zu\n", not_found);
    fprintf(stderr, "Number of SNPs in target  : %zu\n", inclusion.size());
    if(has_ld)
    {
        fprintf(stderr, "\nIn LD Reference %s\n", ld_file.c_str());
        fprintf(stderr,"==============================\n");
        clump.clump_initialize(inclusion, m_snp_list, m_snp_index);
    }
    else clump.clump_initialize(inclusion);
    fprintf(stderr,"\nStart performing clumping\n");

    // This will perform clumping. When completed, inclusion will contain the index SNPs
    // However, this should be changed in later version such that we can also
    // handle different regions
    clump.start_clumping(inclusion, m_snp_list, m_snp_index, c_commander.get_clump_p(),
                         c_commander.get_clump_r2(), c_commander.get_clump_kb(), c_commander.get_proxy());
    // Because it is now call from main, clump will run out of scope, which is what we want

}



// Seems alright now
// This function should give us a nice platform to know where and when to get each SNPs
void PRSice::check_inclusion(std::unordered_map<std::string, size_t> &inclusion, const std::string &c_target_bim_name,
                             size_t &num_ambig, size_t &not_found, size_t &num_duplicate)
{
    std::ifstream target_file;
    target_file.open(c_target_bim_name.c_str());
    if(!target_file.is_open())
    {
        std::string error_message = "Cannot open target bim file: "+c_target_bim_name;
        throw std::runtime_error(error_message);
    }
    std::string line;
    while(std::getline(target_file, line))
    {
        misc::trim(line);
        if(!line.empty())
        {
            std::vector<std::string> token = misc::split(line);
            if(token.size() < 6) throw std::runtime_error("Malformed bim file. Should contain at least 6 column");
            std::string chr = token[+BIM::CHR];
            std::string rsid = token[+BIM::RS];
            int loc = -1;
            int temp = 0;
            try
            {
                temp =misc::convert<int>(token[+BIM::BP]);
                if(temp < 0)
                {
                    std::string error_message = "Negative coordinate of SNP in "+c_target_bim_name;
                    throw std::runtime_error(error_message);
                }
                loc = temp;
            }
            catch(std::runtime_error &error)
            {
                std::string error_message = "Non-numeric coordinate of SNP in "+c_target_bim_name;
                throw std::runtime_error(error_message);
            }
            std::string ref_allele = token[+BIM::A1];
            std::string alt_allele = token[+BIM::A2];

            if(inclusion.find(rsid) == inclusion.end() && m_snp_index.find(rsid)!=m_snp_index.end())
            {
                // will do some soft checking, will issue warning if there are any problem
                // first check if ambiguous
                if( (ref_allele.compare("A")==0 && alt_allele.compare("T")==0) ||
                        (ref_allele.compare("a")==0 && alt_allele.compare("t")==0) ||
                        (ref_allele.compare("T")==0 && alt_allele.compare("A")==0) ||
                        (ref_allele.compare("t")==0 && alt_allele.compare("a")==0) ||
                        (ref_allele.compare("G")==0 && alt_allele.compare("C")==0) ||
                        (ref_allele.compare("g")==0 && alt_allele.compare("c")==0) ||
                        (ref_allele.compare("C")==0 && alt_allele.compare("G")==0) ||
                        (ref_allele.compare("c")==0 && alt_allele.compare("g")==0))
                {
                    num_ambig++;
                }
                else
                {
                    // not ambiguous, now do soft checking
                    size_t index = m_snp_index.at(rsid);
                    if(loc!=-1 && m_snp_list[index].get_loc()!=-1)
                    {
                        bool same = m_snp_list[index].check_loc(chr, loc, ref_allele, alt_allele);
                        if(!same)
                        {
                            fprintf(stderr, "WARNING: %s differ between target and base file\n", rsid.c_str());
                            fprintf(stderr, "         It is advised that you check the files are \n");
                            fprintf(stderr, "         From the same genome build\n");
                        }
                    }
                    inclusion[rsid]=index;
                }
            }
            else if(inclusion.find(rsid)!=inclusion.end())
            {
                num_duplicate++;
            }
            else
            {
                not_found++;
            }
        }
    }
    target_file.close();
}




void PRSice::run_prs(const Commander &c_commander, const std::map<std::string, size_t> &inclusion,
                     const Region &c_region)
{
    // Might want to add additional parameter for the region output
    // keep it for now
    // First, get the phenotype and covariate matrix
    // fam_index is used for storing the index of each individual on the matrix
    // it is used because there might be missing sample which we will exclude from
    // the matrix
    size_t num_pheno = 1;
    std::vector<std::string> pheno_col = c_commander.get_pheno_col();
    std::string target = c_commander.get_target();
    std::string pheno_file = c_commander.get_pheno();
    std::vector<int> pheno_index;
//  Try to obtain the phenotype information (if any)
//	This is to handle multiple phenotype input
    if(pheno_col.size() > num_pheno && !pheno_file.empty())
    {
        std::ifstream pheno;
        pheno.open(pheno_file.c_str());
        if(!pheno.is_open())
        {
            std::string error_message = "Cannot open phenotype file: "+pheno_file;
            throw std::runtime_error(error_message);
        }
        std::string line;
        std::getline(pheno, line);
        pheno.close();
        misc::trim(line);
        std::vector<std::string> col = misc::split(line);
        num_pheno=0;
        std::vector<std::string> temp;
        for(auto i:pheno_col)
        {
            bool found = false;
            for(size_t i_col = 1; i_col < col.size(); ++i_col)
            {
                if(i.compare(col[i_col])==0)
                {
                    num_pheno++;
                    pheno_index.push_back(i_col);
                    temp.push_back(i);
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                fprintf(stderr, "Phenotype: %s cannot be found in phenotype file\n", i.c_str());
            }
        }
        if(temp.size()!= pheno_col.size())
        {
            pheno_col.clear();
            pheno_col=temp;
        }
    }
    m_partition.clear();
    bool pre_run = false;

    categorize(c_commander, inclusion, pre_run);

    bool fastscore = c_commander.fastscore();
    double bound_start =  (fastscore)? c_commander.get_bar_lower(): c_commander.get_lower();
    double bound_end = (fastscore)? c_commander.get_bar_upper():c_commander.get_upper();
    double bound_inter = c_commander.get_inter();
    bool no_regress =c_commander.no_regression();
    bool all = c_commander.all();
    bool target_binary = c_commander.target_is_binary();
    std::string output_name = c_commander.get_out();
    for(size_t i_pheno=0; i_pheno < num_pheno; ++i_pheno)
    {
        std::ofstream all_out;
        if(all)
        {
            std::string all_out_name = output_name+"."+m_current_base;
            if(num_pheno > 1) all_out_name.append("."+pheno_col[i_pheno]+".all.score");
            all_out.open(all_out_name.c_str());
            if(!all_out.is_open())
            {
                std::string error_message = "Cannot open file "+all_out_name+" for write";
                throw std::runtime_error(error_message);
            }
        }
        std::map<std::string,size_t> fam_index;
        // Should contain all samples including those that doesn't have a phenotype
        std::vector<prs_score> prs_fam;
        Eigen::VectorXd phenotype;
        Eigen::MatrixXd independent_variables;
        if(!no_regress)
        {
            // only get the phenotype and covariate files when we need to run the regression
            phenotype = gen_pheno_vec(target, pheno_file, (pheno_index.size()>0)? pheno_index[i_pheno]:i_pheno,
                                      target_binary, fam_index, prs_fam);
            independent_variables = gen_cov_matrix(c_commander.get_cov_file(), c_commander.get_cov_header(), fam_index);
        }
        else
        {
            // initialize the prs_fam vector
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
                    prs_fam.push_back(prs_score(token[1], 0.0));
                }
            }
            fam.close();
        }

        // this contain the information of the best cutoff
        //r2, threshold at best r2 and size at best r2
        for(size_t i_region=0; i_region < c_region.size(); ++i_region)
        {
            m_current_prs.push_back(prs_fam);
        }
        // cur_start_index indiciates how much of quick_ref has been processed
        size_t cur_start_index = 0;
        // num_snp_included should contain the number of SNPs included for each individual
        // regions
        std::vector<size_t> num_snp_included(c_region.size());
        // first read everything that are smaller than the lowest threshold
        if(pre_run) get_prs_score(target, cur_start_index);
        //	This will update the best score to equal to the base. When there is no valid threshold
        //	this should be used to provide the PRS score output
        m_best_score = m_current_prs;

        if(all)
        {
            all_out << "Threshold\tRegion";
            for(size_t i = 0; i < prs_fam.size(); ++i) all_out << "\t" << std::get<0>(prs_fam[i]);
            if(pre_run)
            {
                for(size_t i_region=0; i_region < m_current_prs.size(); ++i_region)
                {
                    for(size_t i_prs_score=0; i_prs_score < m_current_prs[i_region].size(); ++i_prs_score)
                    {
                        all_out << "\t" << std::get<1>(m_current_prs[i_region][i_prs_score])/(double)num_snp_included[i_region];
                    }
                    all_out << std::endl;
                }
            }
            all_out.close();
        }
        //    	Calculate the null here if required. This should speed things up a lot when there are a large amount
        //    	of covariates
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

        // now call function to perform the regression stuff

        calculate_scores(c_commander, c_region, cur_start_index,
                         independent_variables, phenotype, (pheno_col.size()==0)?"":pheno_col[i_pheno], fam_index);
        if(!no_regress)
        {
            // now perform the output for all the best scores and PRSice results
            prs_output(c_commander, c_region, null_r2, (pheno_col.size()==0)?"":pheno_col[i_pheno]);

        }
    }
    fprintf(stderr, "Completed\n");
}

void PRSice::categorize(const Commander &c_commander, const std::map<std::string, size_t> &inclusion, bool &pre_run)
{
    bool fastscore = c_commander.fastscore();
    double bound_start =  (fastscore)? c_commander.get_bar_lower(): c_commander.get_lower();
    double bound_end = (fastscore)? c_commander.get_bar_upper():c_commander.get_upper();
    double bound_inter = c_commander.get_inter();
    bool full_model = c_commander.full();
    std::vector<std::string> file_names;
    if(c_commander.get_target().find("#")!=std::string::npos)
    {
        for(auto chr: m_chr_list)
        {
            std::string name = c_commander.get_target();
            misc::replace_substring(name, "#", chr);
            file_names.push_back(name);
        }
    }
    else
    {
        std::string name = c_commander.get_target();
        file_names.push_back(name);
    }
    for(auto name: file_names)
    {
        std::ifstream bim;
        std::string bim_name = name+".bim";
        bim.open(bim_name.c_str());
        if(!bim.is_open())
        {
            std::string error_message = "Cannot open bim file: " +bim_name;
            throw std::runtime_error(error_message);
        }
        std::string line;
        size_t cur_line = 0;
        while(getline(bim, line))
        {
            misc::trim(line);
            if(!line.empty())
            {
                std::vector<std::string> token = misc::split(line);
                if(token.size() < 6) throw std::runtime_error("Malformed bim file, should contain at least 6 columns");
                if(inclusion.find(token[1])!=inclusion.end())
                {
                    double p = m_snp_list.at(inclusion.at(token[1])).get_p_value();
                    if(p< bound_start)
                    {
                        pre_run=true;
                        m_partition.push_back(p_partition(token[1], cur_line, -1,inclusion.at(token[1])));
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
                        m_partition.push_back(p_partition(token[1], cur_line, (category<0)?-1:category,inclusion.at(token[1]), name));
                    }
                    else if(full_model)
                    {
                        m_partition.push_back(p_partition(token[1], cur_line, (int)(1-bound_start)/bound_inter,inclusion.at(token[1]), name));
                    }
                }
                cur_line++;
            }
        }
        bim.close();
    }
    if(m_partition.size() ==0)
    {
        fprintf(stderr, "None of the SNPs met the threshold\n");
        return;
    }
// now we sort quick_ref such that the smaller p_value categories are always in the front
    std::sort(begin(m_partition), end(m_partition),
              [](p_partition const &t1, p_partition const &t2)
    {
        if(std::get<CATEGORY>(t1)==std::get<CATEGORY>(t2))
        {
            if(std::get<FILENAME>(t1).compare(std::get<FILENAME>(t2))==0)
            {
                return std::get<LINE>(t1)<std::get<LINE>(t2);
            }
            else return std::get<FILENAME>(t1).compare(std::get<FILENAME>(t2))<0;
        }
        else return std::get<CATEGORY>(t1)<std::get<CATEGORY>(t2);
    }
             );
}

Eigen::VectorXd PRSice::gen_pheno_vec(const std::string &c_target, const std::string c_pheno,
                                      const int pheno_index, bool target_binary,
                                      std::map<std::string, size_t> &fam_index,
                                      std::vector<prs_score > &prs_score)
{
    std::vector<double> phenotype_store;
    std::ifstream pheno_file;
    std::string fam_name = c_target+".fam";
    if(fam_name.find("#")!=std::string::npos)
    {
        // we can do this because the fam file should be the same for all chromosome
        misc::replace_substring(fam_name, "#", m_chr_list.front());
    }
    std::string line;
    size_t cur_index = 0;
    // add this information just because I have made the mistake regarding the case
    // control label before
    size_t num_case = 0, num_control = 0;
    if(c_pheno.empty())
    {
        if(pheno_index > 0)
        {
            throw std::runtime_error("When no phenotype file is provided, pheno_index should be 0");
        }
        // Use fam file
        pheno_file.open(fam_name.c_str());
        while(std::getline(pheno_file, line))
        {
            misc::trim(line);
            if(!line.empty())
            {
                std::vector<std::string> token = misc::split(line);
                if(token.size() < 6) throw std::runtime_error("Malformed fam file, should contain at least 6 columns");
                prs_score.push_back(std::pair<std::string, double>(token[1], 0.0));
                if(token[5].compare("NA")!=0)
                {
                    try
                    {
                        if(target_binary)
                        {
                            double temp = misc::convert<int>(token[5]);
                            if(temp-1>=0 && temp-1<2)
                            {
                                fam_index[token[1]]=cur_index;
                                phenotype_store.push_back(temp-1);
                                cur_index++;
                                if(temp==2) num_case++;
                                if(temp==1) num_control++;
                            }
                            // anything other than 1 or 2 will be treated as missing
                        }
                        else
                        {
                            double temp = misc::convert<double>(token[5]);
                            fam_index[token[1]]=cur_index;
                            phenotype_store.push_back(temp);
                            cur_index++;
                        }
                    }
                    catch(const std::runtime_error &error)
                    {
                        // if we can't handle it, it is missing
                    }
                }
            }
        }
        pheno_file.close();
    }
    else
    {
        // Use pheno file
        std::ifstream fam;
        fam.open(fam_name.c_str());
        pheno_file.open(c_pheno.c_str());
        // First form the map using the fam
        std::map<std::string, std::string> pheno_info;
        while(std::getline(pheno_file, line))
        {
            misc::trim(line);
            if(!line.empty())
            {
                std::vector<std::string> token = misc::split(line);
                if(token.size() < pheno_index+1)
                {
                    std::string error_message = "Malformed pheno file, should contain at least "+std::to_string(pheno_index+1)+" columns";
                    throw std::runtime_error(error_message);
                }
                pheno_info[token[0]] = token[pheno_index];
            }
        }
        pheno_file.close();
        while(std::getline(fam, line))
        {
            misc::trim(line);
            if(!line.empty())
            {
                std::vector<std::string> token = misc::split(line);
                if(token.size() < 6) std::runtime_error("Malformed fam file, should contain at least 6 columns");
                prs_score.push_back(std::pair<std::string, double>(token[1], 0.0));
                if(pheno_info.find(token[1])!= pheno_info.end())
                {
                    std::string p = pheno_info[token[1]];
                    if(p.compare("NA")!=0)
                    {
                        try
                        {
                            if(target_binary)
                            {
                                double temp = misc::convert<int>(p);
                                if(temp-1>=0 && temp-1<=2)
                                {
                                    fam_index[token[1]]=cur_index;
                                    phenotype_store.push_back(temp-1);
                                    cur_index++;
                                    if(temp==2) num_case++;
                                    if(temp==1) num_control++;
                                }
                                // anything other than 1 or 2 will be treated as missing
                            }
                            else
                            {
                                double temp = misc::convert<double>(p);
                                fam_index[token[1]]=cur_index;
                                phenotype_store.push_back(temp);
                                cur_index++;
                            }
                        }
                        catch(const std::runtime_error &error) { }
                    }
                }
            }
        }
        fam.close();
    }
    if(phenotype_store.size()==0) throw std::runtime_error("No phenotypes present");
    Eigen::Map<Eigen::VectorXd> res(phenotype_store.data(), phenotype_store.size());
    if(target_binary)
    {
        if(num_control==0) throw std::runtime_error("There are no control samples");
        if(num_case==0) throw std::runtime_error("There are no cases");
        fprintf(stderr,"Number of controls : %zu\n", num_control);
        fprintf(stderr,"Number of cases : %zu\n", num_case);
    }
    else
    {
        fprintf(stderr,"Number of sample(s) with phenotype  : %zu\n", res.rows());
    }
    return res;
}

Eigen::MatrixXd PRSice::gen_cov_matrix(const std::string &c_cov_file, const std::vector<std::string> &c_cov_header, std::map<std::string, size_t> &fam_index)
{
    size_t num_sample = fam_index.size();
    if(c_cov_file.empty())
    {
        // changed this to 2 to avoid needing to add the intercept for every regression
        return Eigen::MatrixXd::Ones(num_sample,2);
    }
    else
    {
        // First read in the header
        std::ifstream cov;
        cov.open(c_cov_file.c_str());
        if(!cov.is_open())
        {
            std::string error_message = "ERROR: Cannot open covariate file: "+c_cov_file;
            throw std::runtime_error(error_message);
        }
        std::string line;
        std::vector<size_t> cov_index;
        int max_index = 0;
        std::getline(cov, line);
        if(!line.empty())
        {
            std::vector<std::string> token = misc::split(line);
            if(c_cov_header.size() == 0)
            {
                // if no header is provided, we will use all the covariates included
                for(size_t i = 1; i < token.size(); ++i) cov_index.push_back(i);
                max_index = cov_index.size()-1;
            }
            else
            {
                std::map<std::string, bool> include;
                // if specific headers are provided, we should only include them
                for(size_t i = 0; i < c_cov_header.size(); ++i) include[c_cov_header[i]]= true;
                for(size_t i = 1; i < token.size(); ++i)
                {
                    if(include.find(token[i])!=include.end())
                    {
                        cov_index.push_back(i);
                        if(i > max_index) max_index=i;
                    }
                }
            }
        }
        else throw std::runtime_error("First line of covariate file is empty!");
        // now we know how much we are working with
        // need to include the space for the polygenic risk score
        // +2 because the first row is for intercept and the second row is for the PRS
        Eigen::MatrixXd result = Eigen::MatrixXd::Ones(num_sample, cov_index.size()+2);
        while(std::getline(cov, line))
        {
            misc::trim(line);
            if(!line.empty())
            {
                std::vector<std::string> token = misc::split(line);
                if(token.size() <= max_index)
                {
                    std::string error_message = "ERROR: Malformed covariate file, should contain at least "+std::to_string(max_index+1)+" column!";
                    throw std::runtime_error(error_message);
                }
                if(fam_index.find(token[1])!= fam_index.end())
                {
                    int index = fam_index[token[1]];
                    for(size_t i = 0; i < cov_index.size(); ++i)
                    {
                        try
                        {
                            double temp = misc::convert<double>(token[cov_index[i]]);
                            result(index, i+2) = temp;
                        }
                        catch(const std::runtime_error &error)
                        {
                            // All missing values are treated as 0
                            result(index, i+2) = 0;
                        }
                    }
                }
            }
        }
        return result;
    }
}

bool PRSice::get_prs_score(const std::string &target, size_t &cur_index)
{
    // Here is the actual calculation of the PRS
    if(m_partition.size()==0) return false; // nothing to do
    int prev_index = std::get<CATEGORY>(m_partition[cur_index]);
    int end_index = 0;
    bool ended =false;
    for(size_t i = cur_index; i < m_partition.size(); ++i)
    {
        if(std::get<CATEGORY>(m_partition[i]) != prev_index && std::get<CATEGORY>(m_partition[i])>=0 )
        {
            end_index = i;
            ended=true;
            break;
        }
        else if(std::get<CATEGORY>(m_partition[i])!=prev_index) prev_index=std::get<CATEGORY>(m_partition[i]); // only when the category is still negative
        // Use as part of the output
        for(size_t i_region=0; i_region< m_num_snp_included.size(); ++i_region)
        {
            if(m_snp_list.at(std::get<INDEX>(m_partition[i])).in(i_region)) m_num_snp_included[i_region]++;
        }
    }
    if(!ended) end_index = m_partition.size();
    PLINK prs(target, m_chr_list);
    prs.initialize();
    prs.get_score(m_partition, m_snp_list, m_current_prs, cur_index, end_index);

    cur_index = end_index;
    return true;
}



void PRSice::calculate_scores(const Commander &c_commander,
                              const Region &c_region, size_t cur_start_index, Eigen::MatrixXd &independent_variables,
                              Eigen::VectorXd &phenotype, const std::string &pheno_name,
                              const std::map<std::string, size_t> &fam_index)
{

    // Now we prepare for the PRS analysis
    bool target_binary = c_commander.target_is_binary();
    std::string target = c_commander.get_target();
    double current_upper =0.0;
    // need to initialize this to use multithreading with EIGEN library
    Eigen::initParallel();
    // This is the storage for therad
    std::vector<std::thread> thread_store;
    // getting the number of thread
    size_t n_thread = c_commander.get_thread();
    // some initialization. TBH, don't think this should slow down too much even if we put it in the loop
    double p_value=0.0, r2=0.0, r2_adjust = 0.0;
    // try to initialize them w.r.t regions such that we don't need to worry about the
    // index in the later analysi
    for(size_t i = 0; i < m_best_threshold.size(); ++i)
    {
        m_best_threshold[i] = PRSice_best(0,0,0);
        m_prs_results.push_back(std::vector<PRSice_result>(0));
    }
    // now start going through all the thresholds
    // when cur_start_index == quick_ref.size(), all SNPs should have processed
    if(cur_start_index == m_partition.size())
    {
        fprintf(stderr, "There are no valid threshold to test\n");
        fprintf(stderr, "All SNPs have p-value lower than --lower and higher than --upper\n");
        fprintf(stderr, "We will output the PRS score for you just to be nice\n");
    }

    bool fastscore = c_commander.fastscore();
    double bound_start =  (fastscore)? c_commander.get_bar_lower(): c_commander.get_lower();
    double bound_end = (fastscore)? c_commander.get_bar_upper():c_commander.get_upper();
    double bound_inter = c_commander.get_inter();
    bool no_regress =c_commander.no_regression();
    bool all = c_commander.all();

    std::ofstream all_out;
    if(all)
    {
        std::string all_out_name = c_commander.get_out()+"."+m_current_base;
        if(!pheno_name.empty()) all_out_name.append("."+pheno_name+".all.score");
        all_out.open(all_out_name.c_str(), std::ofstream::app);
        if(!all_out.is_open())
        {
            std::string error_message = "Cannot open file "+all_out_name+" for write";
            throw std::runtime_error(error_message);
        }
    }
    while(cur_start_index != m_partition.size())
    {
        // getting the current cutoff
        current_upper = std::min((std::get<CATEGORY>(m_partition[cur_start_index])+1)*bound_inter+bound_start, bound_end);
        fprintf(stderr, "\rProcessing %f", current_upper);

        // now calculate the PRS for each region
        bool reg = get_prs_score(target, cur_start_index);
        // if regression is not required, we will simply output the score
        if(all)
        {
            for(size_t i_region=0; i_region < m_current_prs.size(); ++i_region)
            {
                all_out << current_upper << "\t" << c_region.get_name(i_region);
                for(size_t i_reg_score=0; i_reg_score < m_current_prs[i_region].size(); ++i_reg_score)
                {
                    all_out << "\t" << std::get<1>(m_current_prs[i_region][i_reg_score])/(double)m_num_snp_included[i_region];
                }
                all_out << std::endl;
            }
        }
        // update the boolean, basically, if we have added new SNPs AND require regression
        // we will perform the regression analysis
        reg=reg&&!no_regress;
        if(reg)
        {
            // here we do multithreading
            if(n_thread == 1 || m_current_prs.size()==1)
            {
                thread_score(independent_variables, phenotype, fam_index, 0, m_current_prs.size(),
                             target_binary, current_upper,n_thread);
            }
            else
            {
                // perform multi threading
                if(c_region.size() < n_thread)
                {
                    for(size_t i_region = 0; i_region < c_region.size(); ++i_region)
                    {
                        thread_store.push_back(std::thread(&PRSice::thread_score, this,
                                                           std::ref(independent_variables), std::cref(phenotype), std::cref(fam_index),
                                                           i_region, i_region+1, target_binary, current_upper,1));
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
                                                           std::ref(independent_variables), std::cref(phenotype), std::cref(fam_index),
                                                           start, ending, target_binary, current_upper,1));
                        start=ending;
                        remain--;
                    }
                }
                // joining the threads
                for(size_t i_thread = 0; i_thread < thread_store.size(); ++i_thread)
                {
                    thread_store[i_thread].join();
                }
            }
        }
    }
    fprintf(stderr, "\n");
    if(all) all_out.close();
}


void PRSice::thread_score( Eigen::MatrixXd &independent_variables, const Eigen::VectorXd &c_pheno,
                           const std::map<std::string, size_t> &c_fam_index, size_t region_start, size_t region_end,
                           bool target_binary, double threshold, size_t thread)
{
    Eigen::MatrixXd X;
    bool thread_safe=false;
    // so we will only copy the matrix when it is not thread safe to do so
    if(region_start==0 && region_end ==m_current_prs.size()) thread_safe = true;
    else X = independent_variables;
    // c_prs_region_score = prs score of the region at the current threshold
    // c_num_snp_include = num of SNP in region at current threshold
    // c_fam_index = index for each individual on the matrix (mainly to deal with missing data)
    // prs_best_info = threshold info with best score. p here is the threshold, not the p-value for PRSice
    // prs_best_score = store the best score
    // region_result = PRSice result, containing all the required information for output
    //std::vector<PRSice::PRSice_result> temp_region_result;
    double r2 = 0.0, r2_adjust=0.0, p_value = 0.0;
    for(size_t iter = region_start; iter < region_end; ++iter)
    {
        std::vector<prs_score>::const_iterator prs_iter = m_current_prs.at(iter).begin();
        std::vector<prs_score>::const_iterator prs_end = m_current_prs.at(iter).end();
        for(; prs_iter != prs_end; ++prs_iter)
        {
            std::string sample = std::get<IID>(*prs_iter);
            if(c_fam_index.find(sample)!=c_fam_index.end())
            {
                if(thread_safe) independent_variables(c_fam_index.at(sample), 1) = std::get<PRS>(*prs_iter)/(double)m_num_snp_included[iter];
                else X(c_fam_index.at(sample), 1) = std::get<PRS>(*prs_iter)/(double)m_num_snp_included[iter];
            }
        }
        if(target_binary)
        {
            try
            {
                if(thread_safe) Regression::glm(c_pheno, independent_variables, p_value, r2, 25, thread, true);
                else Regression::glm(c_pheno, X, p_value, r2, 25, thread, true);
            }
            catch(const std::runtime_error &error)
            {
                // This should only happen when the glm doesn't converge.
                // Let's hope that won't happen...
                fprintf(stderr, "ERROR: GLM model did not converge!\n");
                fprintf(stderr, "       Please send me the DEBUG files\n");
                std::ofstream debug;
                debug.open("DEBUG");
                if(thread_safe) debug << independent_variables << std::endl;
                else 	debug << X<< std::endl;
                debug.close();
                debug.open("DEBUG.y");
                debug << c_pheno << std::endl;
                debug.close();
                std::cerr << "ERROR: " << error.what() << std::endl;
                exit(-1);
            }
        }
        else
        {
            if(thread_safe) Regression::linear_regression(c_pheno, independent_variables, p_value, r2, r2_adjust, thread, true);
            else Regression::linear_regression(c_pheno, X, p_value, r2, r2_adjust, thread, true);
        }
        // This should be thread safe as each thread will only mind their own region
        // now add the PRS result to the vectors (hopefully won't be out off scope
        m_prs_results.at(iter).push_back(PRSice_result(threshold, r2, m_num_snp_included.at(iter), r2_adjust, p_value));
        // It this is the best r2, then we will add it
        if(std::get<R2>(m_best_threshold[iter]) < r2)
        {
            m_best_threshold[iter] = PRSice_best(threshold, r2, m_num_snp_included.at(iter));
            m_best_score[iter] = m_current_prs.at(iter);
        }
    }
}



// This function should be responsible for the output
void PRSice::prs_output(const Commander &c_commander, const Region &c_region, const double null_r2,
                        const std::string pheno_name) const
{

    std::string output_prefix = c_commander.get_out()+"."+m_current_base;
    if(!pheno_name.empty()) output_prefix.append("."+pheno_name+".");
    for(size_t i_region = 0; i_region< m_prs_results.size(); ++i_region)
    {
        std::string output_name = output_prefix+"."+c_region.get_name(i_region);
        std::string out_best = output_name+".best";
        std::string out_prsice = output_name+".prsice";
        std::ofstream best_out, prsice_out;
        best_out.open(out_best.c_str());
        prsice_out.open(out_prsice.c_str());
        if(!best_out.is_open())
        {
            std::string error_message = "ERROR: Cannot open file: " +out_best+" to write";
            throw std::runtime_error(error_message);
        }
        if(!prsice_out.is_open())
        {
            std::string error_message = "ERROR: Cannot open file: " +out_prsice+" to write";
            throw std::runtime_error(error_message);
        }
        best_out << "IID\tprs_"<<std::get<THRESHOLD>(m_best_threshold[i_region]) << std::endl;
        prsice_out << "Threshold\tR2\tP\tNum_SNP" << std::endl;
        // We want to skip the intercept for now
        for(size_t i_prsice=0; i_prsice< m_prs_results[i_region].size(); ++i_prsice)
        {
            prsice_out << std::get<THRESHOLD>(m_prs_results[i_region][i_prsice]) << "\t" <<
                       std::get<R2>(m_prs_results[i_region][i_prsice])-null_r2 << "\t" <<
                       std::get<P>(m_prs_results[i_region][i_prsice])<< "\t" <<
                       std::get<NSNP>(m_prs_results[i_region][i_prsice]) << std::endl;
        }
        int best_snp_size = std::get<NSNP>(m_best_threshold[i_region]);
        if(best_snp_size==0)
        {
            fprintf(stderr, "ERROR: Best R2 obtained when no SNPs were included\n");
            fprintf(stderr, "       Cannot output the best PRS score\n");
        }
        else
        {
            for(size_t i_score=0; i_score < m_best_score[i_region].size(); ++i_score)
            {
                best_out << std::get<IID>(m_best_score[i_region][i_score]) << "\t" <<
                         std::get<PRS>(m_best_score[i_region][i_score])/(double)best_snp_size<< std::endl;
            }
        }
        prsice_out.close();
        best_out.close();
    }
}

// basically update the score vector to contain the new polygenic score
PRSice::PRSice()
{
    //ctor
}

PRSice::~PRSice()
{
    //dtor
}
