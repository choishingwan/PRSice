#include "prsice.hpp"

std::mutex PRSice::score_mutex;


// Seems to be working alright (hope I know how to write test cases...)
void PRSice::get_snp(const Commander &c_commander, Region &region, const double &c_threshold)
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
                        else if(pvalue > c_threshold)
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

    m_snp_list.sort();
    size_t before = m_snp_list.size();
    m_snp_list.erase(std::unique(m_snp_list.begin(), m_snp_list.end()), m_snp_list.end());
    size_t after = m_snp_list.size();
    std::map<std::string, bool> unique_chr;
    for(size_t i_snp = 0; i_snp < m_snp_list.size(); ++i_snp)
    {
        m_snp_index[m_snp_list[i_snp].get_rs_id()]=i_snp;
        if(unique_chr.find(m_snp_list[i_snp].get_chr())==unique_chr.end())
        {
            unique_chr[m_snp_list[i_snp].get_chr()] = true;
            m_chr_list.push_back(m_snp_list[i_snp].get_chr());
        }
        if(index[0] >=0 && index[5]>=0)
        {
            m_snp_list[i_snp].set_flag(region.check(m_snp_list[i_snp].get_chr(), m_snp_list[i_snp].get_loc()));
        }
    }


    num_duplicated = (int)before-(int)after;
    fprintf(stderr, "Number of SNPs from base  : %zu\n", m_snp_list.size());
    if(num_indel!=0) fprintf(stderr, "Number of Indels          : %zu\n", num_indel);
    if(num_duplicated!=0) fprintf(stderr, "Number of duplicated SNPs : %d\n", num_duplicated);
    if(num_stat_not_convertible!=0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_convertible);
    if(num_p_not_convertible!=0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_convertible);
    if(num_se_not_convertible!=0) fprintf(stderr, "Failed to convert %zu SE\n", num_se_not_convertible);
}

void PRSice::clump(const Commander &c_commander)
{
    std::string target = c_commander.get_target();
    bool has_ld = !c_commander.ld_prefix().empty();
    std::string ld_file = (has_ld)? c_commander.ld_prefix(): target;
    PLINK clump(ld_file, m_chr_list, c_commander.get_thread());
    // Because we will go through the target anyway, first go through the target for inclusion
    size_t num_ambig=0, not_found=0, num_duplicate=0;
    std::string target_bim_name = target+".bim";
    if(target_bim_name.find("#")!=std::string::npos)
    {
        for(auto &&chr: m_chr_list)
        {
            std::string target_chr_bim_name = target_bim_name;
            misc::replace_substring(target_chr_bim_name, "#", chr);
            check_inclusion(target_chr_bim_name, num_ambig, not_found, num_duplicate);
        }
    }
    else check_inclusion(target_bim_name, num_ambig, not_found, num_duplicate);
    fprintf(stderr, "\nIn Target File\n");
    fprintf(stderr,"==============================\n");
    if(num_ambig != 0)	fprintf(stderr, "Number of ambiguous SNPs  : %zu\n", num_ambig);
    if(num_duplicate != 0) fprintf(stderr, "Number of duplicated SNPs : %zu\n", num_duplicate);
    if(not_found != 0) fprintf(stderr, "Number of SNPs not found  : %zu\n", not_found);
    fprintf(stderr, "Number of SNPs in target  : %zu\n", m_include_snp.size());
    if(has_ld)
    {
        fprintf(stderr, "\nIn LD Reference %s\n", ld_file.c_str());
        fprintf(stderr,"==============================\n");
        clump.clump_initialize(m_include_snp, m_snp_list, m_snp_index);
    }
    else clump.clump_initialize(m_include_snp);
    fprintf(stderr,"\nStart performing clumping\n");
    clump.start_clumping(m_include_snp, m_snp_list, c_commander.get_clump_p(),
                         c_commander.get_clump_r2(), c_commander.get_clump_kb(), c_commander.get_proxy());
}

void PRSice::check_inclusion(const std::string &c_target_bim_name, size_t &num_ambig, size_t &not_found,
		size_t &num_duplicate)
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

            if(m_include_snp.find(rsid) == m_include_snp.end() && m_snp_index.find(rsid)!=m_snp_index.end())
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
                    m_include_snp[rsid]=index;
                }
            }
            else if(m_include_snp.find(rsid)!=m_include_snp.end())
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


void PRSice::run_prs(const Commander &c_commander, const Region &c_region)
{
	std::vector<std::string> pheno_header = c_commander.get_pheno_col();
	std::vector<pheno_storage> pheno_index;
	std::string target = c_commander.get_target();
    std::string pheno_file = c_commander.get_pheno();
    std::string output_name = c_commander.get_out();

    bool pre_run=false;
    categorize(c_commander, pre_run);

    if(pheno_header.size() !=0&& pheno_file.empty())
    {
    		throw std::runtime_error("You must provide a phenotype file for multiple phenotype analysis");
    }
    if(pheno_file.empty())
    {
    		std::string fam = target+".fam";
    		pheno_storage temp;
    		std::get<pheno_store::FILE_NAME>(temp) = fam;
    		std::get<pheno_store::INDEX>(temp) = +FAM::PHENOTYPE;
    		std::get<pheno_store::NAME>(temp) = "";
    		pheno_index.push_back(std::tuple<std::string, size_t, std::string>());
    }
    else{
    		std::ifstream pheno;
    		pheno.open(pheno_file.c_str());
    		if(!pheno.is_open())
    		{
    			std::string error_message = "Cannot open phenotype file: "+pheno_file;
    			throw std::runtime_error(error_message);
    		}
    		std::string line;
    		std::getline(pheno, line);
    		if(line.empty())
    		{
    			throw std::runtime_error("Cannot have empty header line for phenotype file!");
    		}
    		pheno.close();
    		misc::trim(line);
    		std::vector<std::string> col = misc::split(line);
    		bool found = false;
    		std::sort(pheno_header.begin(), pheno_header.end());
    		pheno_header.erase(std::unique(pheno_header.begin(), pheno_header.end()), pheno_header.end());
    		for(auto &&pheno_info : pheno_header)
    		{
    			found = false;
    			for(size_t i_column =0; i_column < col.size(); ++i_column)
    			{
    				if(col[i_column].compare(pheno_info)==0)
    				{
    					found = true;
    					pheno_storage temp;
		    			std::get<pheno_store::FILE_NAME>(temp) = pheno_file;
		    			std::get<pheno_store::INDEX>(temp) = i_column;
		    			std::get<pheno_store::NAME>(temp) = pheno_info;
    		    			pheno_index.push_back(temp);
    					break;
    				}
    			}
    			if(!found)
    			{
    				fprintf(stderr, "Phenotype: %s cannot be found in phenotype file\n", pheno_info.c_str());
    			}
    		}
    }

    //Start processing each individual phenotypes
    for(auto &&pheno_info : pheno_index)
    {
        m_prs_results.clear();
        m_best_threshold.clear();
        m_best_score.clear();
        m_num_snp_included.clear();
        m_current_prs.clear();
    		individual_pheno_prs(c_commander, c_region, pre_run, pheno_info, pheno_index.size()> 1);
    }

}

void PRSice::categorize(const Commander &c_commander, bool &pre_run)
{
    bool fastscore = c_commander.fastscore();
    double bound_start =  (fastscore)? c_commander.get_bar_lower(): c_commander.get_lower();
    double bound_end = (fastscore)? c_commander.get_bar_upper():c_commander.get_upper();
    double bound_inter = c_commander.get_inter();
    bool full_model = c_commander.full();
    std::vector<std::string> file_names;
    if(c_commander.get_target().find("#")!=std::string::npos)
    {
        for(auto &&chr: m_chr_list)
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
    for(auto &&name: file_names)
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
                if(m_include_snp.find(token[+BIM::RS])!=m_include_snp.end())
                {
                    double p = m_snp_list.at(m_include_snp.at(token[+BIM::RS])).get_p_value();
                    if(p< bound_start)
                    {
                        pre_run=true;
                        p_partition part;
                        std::get<+PRS::RS>(part) = token[+BIM::RS];
                        std::get<+PRS::LINE>(part) = cur_line;
                        std::get<+PRS::CATEGORY>(part) = -1;
                        std::get<+PRS::INDEX>(part) = m_include_snp.at(token[+BIM::RS]);
                        std::get<+PRS::FILENAME>(part) = name;
                        m_partition.push_back(part);
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
                        p_partition part;
                        std::get<+PRS::RS>(part) = token[+BIM::RS];
                        std::get<+PRS::LINE>(part) = cur_line;
                        std::get<+PRS::CATEGORY>(part) = (category<0)?-1:category;
                        std::get<+PRS::INDEX>(part) = m_include_snp.at(token[+BIM::RS]);
                        std::get<+PRS::FILENAME>(part) = name;
                        m_partition.push_back(part);
                    }
                    else if(full_model)
                    {
                        p_partition part;
                        std::get<+PRS::RS>(part) = token[+BIM::RS];
                        std::get<+PRS::LINE>(part) = cur_line;
                        std::get<+PRS::CATEGORY>(part) = (int)(1-bound_start)/bound_inter;
                        std::get<+PRS::INDEX>(part) = m_include_snp.at(token[+BIM::RS]);
                        std::get<+PRS::FILENAME>(part) = name;
                        m_partition.push_back(part);
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
    std::sort(begin(m_partition), end(m_partition),
              [](p_partition const &t1, p_partition const &t2)
				{
					if(std::get<+PRS::CATEGORY>(t1)==std::get<+PRS::CATEGORY>(t2))
					{
						if(std::get<+PRS::FILENAME>(t1).compare(std::get<+PRS::FILENAME>(t2))==0)
						{
							return std::get<+PRS::LINE>(t1)<std::get<+PRS::LINE>(t2);
						}
						else return std::get<+PRS::FILENAME>(t1).compare(std::get<+PRS::FILENAME>(t2))<0;
					}
					else return std::get<+PRS::CATEGORY>(t1)<std::get<+PRS::CATEGORY>(t2);
				}
    );
}
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
		calculate_scores(c_commander, c_region, cur_start_index,
				independent_variables, phenotype, std::get<pheno_store::NAME>(pheno_info), sample_index);
		if(!no_regress)
		{
			prs_output(c_commander, c_region, null_r2, std::get<pheno_store::NAME>(pheno_info));
		}
	}

}

void PRSice::gen_pheno_vec(Eigen::VectorXd &phenotype, const std::string &c_target, const std::string c_pheno,
    		const int pheno_index, bool target_binary, std::unordered_map<std::string, size_t> &sample_index,
			std::vector<prs_score> &sample_prs)
{

	std::vector<double> phenotype_store;
	std::ifstream pheno_file;
	std::string fam_name = c_target+".fam";
	if(fam_name.find("#")!=std::string::npos)
	{
		misc::replace_substring(fam_name, "#", m_chr_list.front());
	}
	std::string line;
	size_t cur_index;
	size_t num_case=0, num_control=0;
	if(fam_name.compare(c_pheno)==0)
	{
		pheno_file.open(fam_name.c_str());
		if(!pheno_file.is_open())
		{
			std::string error_message = "Cannot open phenotype file: " + fam_name;
			throw std::runtime_error(error_message);
		}
		while(std::getline(pheno_file, line))
		{
			misc::trim(line);
			if(!line.empty())
			{
				std::vector<std::string> token = misc::split(line);
				if(token.size() < 6) throw std::runtime_error("Malformed fam file, should contain at least 6 columns");
				sample_prs.push_back(prs_score(token[+FAM::IID], 0.0));
				if(token[+FAM::PHENOTYPE].compare("NA")!=0)
				{
					try
					{
						if(target_binary)
						{
							double temp = misc::convert<int>(token[+FAM::PHENOTYPE]);
							if(temp-1>=0 && temp-1<2)
							{
								sample_index[token[+FAM::IID]]=cur_index;
								phenotype_store.push_back(temp-1);
								cur_index++;
								if(temp==2) num_case++;
								if(temp==1) num_control++;
							}
						}
						else
						{
							double temp = misc::convert<double>(token[+FAM::PHENOTYPE]);
							sample_index[token[+FAM::IID]]=cur_index;
							phenotype_store.push_back(temp);
							cur_index++;
						}
					}
					catch(const std::runtime_error &error){}
				}
			}
		}
		pheno_file.close();
	}
	else
	{
		std::ifstream fam;
		fam.open(fam_name.c_str());
		pheno_file.open(c_pheno.c_str());
		if(!fam.is_open())
		{
			std::string error_message = "Cannot open fam file: " + fam_name;
			throw std::runtime_error(error_message);
		}
		if(!pheno_file.is_open())
		{
			std::string error_message = "Cannot open phenotype file: " + fam_name;
			throw std::runtime_error(error_message);
		}
		// Main problem: we want the order following the fam file
		std::unordered_map<std::string, std::string> phenotype_info;
		cur_index = 0;
		while(std::getline(pheno_file, line))
		{
			misc::trim(line);
			if(!line.empty()){
				std::vector<std::string> token = misc::split(line);
				if(token.size() < pheno_index+1)
				{
                    std::string error_message = "Malformed pheno file, should contain at least "+std::to_string(pheno_index+1)+" columns";
                    throw std::runtime_error(error_message);
                }
				phenotype_info[token[0]] = token[pheno_index];
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
				sample_prs.push_back(prs_score(token[+FAM::IID], 0.0));
				if(phenotype_info.find(token[+FAM::IID])!= phenotype_info.end())
				{
					std::string p = phenotype_info[token[+FAM::IID]];
					if(p.compare("NA")!=0)
					{
						try
						{
							if(target_binary)
							{
								double temp = misc::convert<int>(p);
								if(temp-1>=0 && temp-1<=2)
								{
									sample_index[token[+FAM::IID]]=cur_index;
									phenotype_store.push_back(temp-1);
									cur_index++;
									if(temp==2) num_case++;
									if(temp==1) num_control++;
								}
							}
							else
							{
								double temp = misc::convert<double>(p);
								sample_index[token[+FAM::IID]]=cur_index;
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
	if(phenotype_store.size() == 0) throw std::runtime_error("No phenotype presented");
	phenotype = Eigen::Map<Eigen::VectorXd>(phenotype_store.data(), phenotype_store.size());
	if(target_binary)
	{
		if(num_control==0) throw std::runtime_error("There are no control samples");
		if(num_case==0) throw std::runtime_error("There are no cases");
		fprintf(stderr,"Number of controls : %zu\n", num_control);
		fprintf(stderr,"Number of cases : %zu\n", num_case);
	}
	else
	{
		fprintf(stderr,"Number of sample(s) with phenotype  : %zu\n", phenotype.rows());
	}
}

void PRSice::gen_cov_matrix(Eigen::MatrixXd &independent_variables, const std::string &c_cov_file,
		const std::vector<std::string> &c_cov_header, std::unordered_map<std::string, size_t> &sample_index)
{
    size_t num_sample = sample_index.size();
    if(c_cov_file.empty())
    {
        independent_variables = Eigen::MatrixXd::Ones(num_sample,2);
        return;
    }
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
    			std::unordered_map<std::string, bool> included;
    			// if specific headers are provided, we should only include them
    			for(size_t i = 0; i < c_cov_header.size(); ++i) included[c_cov_header[i]]= true;
    			for(size_t i_header = 1; i_header < token.size(); ++i_header)
    			{
    				if(included.find(token[i_header])!=included.end())
    				{
    					cov_index.push_back(i_header);
    					if(i_header > max_index) max_index=i_header;
    				}
    			}
    		}
    }
    else throw std::runtime_error("First line of covariate file is empty!");
    independent_variables = Eigen::MatrixXd::Ones(num_sample, cov_index.size()+2);
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
    			if(sample_index.find(token[0])!= sample_index.end())
    			{
    				int index = sample_index[token[0]];
    				for(size_t i = 0; i < cov_index.size(); ++i)
    				{
    					try
    					{
    						double temp = misc::convert<double>(token[cov_index[i]]);
    						independent_variables(index, i+2) = temp; // + 2 because first line = intercept, second line = PRS
    					}
    					catch(const std::runtime_error &error)
    					{
    						independent_variables(index, i+2) = 0;
    					}
    				}
    			}
    		}
    }
}

bool PRSice::get_prs_score(const std::string &target, size_t &cur_index)
{
    if(m_partition.size()==0) return false; // nothing to do
    int prev_index = std::get<+PRS::CATEGORY>(m_partition[cur_index]);
    int end_index = 0;
    bool ended =false;
    for(size_t i = cur_index; i < m_partition.size(); ++i)
    {
        if(std::get<+PRS::CATEGORY>(m_partition[i]) != prev_index && std::get<+PRS::CATEGORY>(m_partition[i])>=0 )
        {
            end_index = i;
            ended=true;
            break;
        }
        else if(std::get<+PRS::CATEGORY>(m_partition[i])!=prev_index) prev_index=std::get<+PRS::CATEGORY>(m_partition[i]); // only when the category is still negative
        // Use as part of the output
        for(size_t i_region=0; i_region< m_num_snp_included.size(); ++i_region)
        {
            if(m_snp_list[std::get<+PRS::INDEX>(m_partition[i])].in(i_region)) m_num_snp_included[i_region]++;
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
                              const std::unordered_map<std::string, size_t> &sample_index)
{
    bool target_binary = c_commander.target_is_binary();
    std::string target = c_commander.get_target();
    double current_upper =0.0;
    Eigen::initParallel();
    std::vector<std::thread> thread_store;
    size_t n_thread = c_commander.get_thread();
    double p_value=0.0, r2=0.0, r2_adjust = 0.0;
    for(size_t i = 0; i < m_best_threshold.size(); ++i)
    {
        m_best_threshold[i] = PRSice_best(0,0,0);
        m_prs_results.push_back(std::vector<PRSice_result>(0));
    }
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
        current_upper = std::min((std::get<+PRS::CATEGORY>(m_partition[cur_start_index])+1)*bound_inter+bound_start, bound_end);
        fprintf(stderr, "\rProcessing %f", current_upper);

        // now calculate the PRS for each region
        bool reg = get_prs_score(target, cur_start_index);
        // if regression is not required, we will simply output the score
        if(all)
        {
            for(size_t i_region=0; i_region < m_current_prs.size(); ++i_region)
            {
                all_out << current_upper << "\t" << c_region.get_name(i_region);
                for(auto &&prs : m_current_prs[i_region])
                {
                		all_out << "\t" << std::get<+PRS::PRS>(prs)/(double)m_num_snp_included[i_region];
                }
                all_out << std::endl;
            }
        }
        reg=reg&&!no_regress;
        if(reg)
        {
        		if(n_thread == 1 || m_current_prs.size()==1)
            {
        			thread_score(independent_variables, phenotype, sample_index, 0, m_current_prs.size(),
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
                        		std::ref(independent_variables), std::cref(phenotype), std::cref(sample_index),
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
                        		std::ref(independent_variables), std::cref(phenotype), std::cref(sample_index),
								start, ending, target_binary, current_upper,1));
                        start=ending;
                        remain--;
                    }
                }
                // joining the threads
                for(auto &&thread : thread_store) thread.join();
            }
        }
    }
    fprintf(stderr, "\n");
    if(all) all_out.close();
}


void PRSice::thread_score( Eigen::MatrixXd &independent_variables, const Eigen::VectorXd &c_pheno,
                           const std::unordered_map<std::string, size_t> &c_sample_index, size_t region_start, size_t region_end,
                           bool target_binary, double threshold, size_t thread)
{
    Eigen::MatrixXd X;
    bool thread_safe=false;
    if(region_start==0 && region_end ==m_current_prs.size()) thread_safe = true;
    else X = independent_variables;
    double r2 = 0.0, r2_adjust=0.0, p_value = 0.0;
    for(size_t iter = region_start; iter < region_end; ++iter)
    {
    		for(auto prs : m_current_prs.at(iter))
    		{
    			std::string sample = std::get<+PRS::IID>(prs);
    			if(c_sample_index.find(sample)!=c_sample_index.end())
    			{
    				if(thread_safe) independent_variables(c_sample_index.at(sample), 1) = std::get<+PRS::PRS>(prs)/(double)m_num_snp_included[iter];
    				else X(c_sample_index.at(sample), 1) = std::get<+PRS::PRS>(prs)/(double)m_num_snp_included[iter];
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
        PRSice_result res;
        std::get<+PRS::THRESHOLD>(res) = threshold;
        std::get<+PRS::R2>(res) = r2;
        std::get<+PRS::NSNP>(res) =  m_num_snp_included[iter];
        std::get<+PRS::R2ADJ>(res) =  r2_adjust;
        std::get<+PRS::P>(res) =  p_value;
        m_prs_results[iter].push_back(res);
        // It this is the best r2, then we will add it
        if(std::get<+PRS::R2>(m_best_threshold[iter]) < r2)
        {
        		PRSice_best best;
        		std::get<+PRS::THRESHOLD>(best) = threshold;
        		std::get<+PRS::R2>(best) = r2;
        		std::get<+PRS::NSNP>(best) =  m_num_snp_included[iter];
            m_best_threshold[iter] = best;
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
        best_out << "IID\tprs_"<<std::get<+PRS::THRESHOLD>(m_best_threshold[i_region]) << std::endl;
        prsice_out << "Threshold\tR2\tP\tNum_SNP" << std::endl;
        // We want to skip the intercept for now
        for(auto &&prs : m_prs_results[i_region])
        {
        		prsice_out << std::get<+PRS::THRESHOLD>(prs) << "\t" <<
        				std::get<+PRS::R2>(prs)-null_r2 << "\t" <<
						std::get<+PRS::P>(prs)<< "\t" <<
						std::get<+PRS::NSNP>(prs) << std::endl;
        }
        int best_snp_size = std::get<+PRS::NSNP>(m_best_threshold[i_region]);
        if(best_snp_size==0)
        {
            fprintf(stderr, "ERROR: Best R2 obtained when no SNPs were included\n");
            fprintf(stderr, "       Cannot output the best PRS score\n");
        }
        else
        {
            for(auto &&prs : m_best_score[i_region])
            {
                best_out << std::get<+PRS::IID>(prs) << "\t" <<
                         std::get<+PRS::PRS>(prs)/(double)best_snp_size<< std::endl;
            }
        }
        prsice_out.close();
        best_out.close();
    }
}


void PRSice::run_prslice(const Commander &c_commander){
	int window_size = c_commander.prslice();
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
