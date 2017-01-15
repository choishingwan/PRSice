#include "prsice.hpp"

void PRSice::prslice_windows(const Commander &c_commander, const Region &c_region)
{
	std::string prev_chr = "";
	size_t prev_loc = 0;
	size_t wind_size = c_commander.prslice();
	m_partition.clear();
	m_best_snps.clear();
	bool fastscore = c_commander.fastscore();
	double bound_start =  (fastscore)? c_commander.get_bar_lower(): c_commander.get_lower();
	double bound_end = (fastscore)? c_commander.get_bar_upper():c_commander.get_upper();
	double bound_inter = c_commander.get_inter();
	std::unordered_map<std::string, std::string> chr_files;
	bool full_model = c_commander.full();
	bool multi_chr = c_commander.get_target().find("#")!=std::string::npos;
	for(auto &&chr : m_chr_list)
	{
		if(multi_chr)
		{
			std::string name = m_target;
			misc::replace_substring(name, "#", chr);
			chr_files[chr] = name;
		}
		else
		{
			chr_files[chr] = m_target;
		}
	}

	std::unordered_map<std::string, size_t> partition_index;
	for(auto snp : m_snp_list)
	{
		if(m_include_snp.find(snp.get_rs_id())!=m_include_snp.end())
		{
			std::string cur_chr = snp.get_chr();
			size_t cur_loc = snp.get_loc();
			if(m_partition.size()!=0 && (prev_chr.empty() || cur_chr != prev_chr || cur_loc - prev_loc > wind_size))
			{
				//get line number, must be same file as same chromosome
				update_line(partition_index);
				// run prsice
				prsice(c_commander, c_region, true);
				std::string window_name = cur_chr+":"+std::to_string(m_snp_list[std::get<+PRS::INDEX>(m_partition.front())].get_loc()) +
						"-"+std::to_string(m_snp_list[std::get<+PRS::INDEX>(m_partition.back())].get_loc());
				size_t num_snps = std::get<+PRS::NSNP>(m_best_threshold[0]);
				double best_r2 = std::get<+PRS::R2>(m_best_threshold[0]);
				std::vector<p_partition> best_snp;
				for(size_t i_snp=0; i_snp < num_snps; ++i_snp)
				{
					best_snp.push_back(m_partition[i_snp]);
				}
				windows wind;
				std::get<prslice_wind::WIND>(wind) = window_name;
				std::get<prslice_wind::R2>(wind) = best_r2;
				std::get<prslice_wind::SNPS>(wind) = best_snp;
				std::get<prslice_wind::P>(wind) = std::get<+PRS::P>(m_best_threshold[0]);
				std::get<prslice_wind::NSNP>(wind) = std::get<+PRS::NSNP>(m_best_threshold[0]);
				std::get<prslice_wind::COEFF>(wind) = std::get<+PRS::COEFF>(m_best_threshold[0]);
				m_best_snps.push_back(wind);
				m_best_threshold.clear();
				m_partition.clear();
				partition_index.clear();
			}
			else
			{
				double p = snp.get_p_value();
				p_partition part;
				std::get<+PRS::RS>(part) = snp.get_rs_id();
				std::get<+PRS::LINE>(part) = 0;
				std::get<+PRS::INDEX>(part) = m_include_snp[snp.get_rs_id()];
				std::get<+PRS::FILENAME>(part) = chr_files[cur_chr];
				if(p<bound_end)
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
					std::get<+PRS::CATEGORY>(part) = (category<0)?0:category;
					m_partition.push_back(part);
					partition_index[snp.get_rs_id()] = m_partition.size()-1;
				}
				else if(full_model)
				{
					std::get<+PRS::CATEGORY>(part) = (int)((bound_end+0.1)-bound_start)/bound_inter; // This ensure they all fall into the same category
					m_partition.push_back(part);
					partition_index[snp.get_rs_id()] = m_partition.size()-1;
				}
			}
		}
	}
	if(m_partition.size()!=0)
	{
		size_t front_index = std::get<+PRS::INDEX>(m_partition.front());
		update_line(partition_index);
		std::string window_name = m_snp_list[front_index].get_chr()+":"+std::to_string(m_snp_list[front_index].get_loc()) +
				"-"+std::to_string(m_snp_list[std::get<+PRS::INDEX>(m_partition.back())].get_loc());
		size_t num_snps = std::get<+PRS::NSNP>(m_best_threshold[0]);
		double best_r2 = std::get<+PRS::R2>(m_best_threshold[0]);
		std::vector<p_partition> best_snp;
		for(size_t i_snp=0; i_snp < num_snps; ++i_snp)
		{
			best_snp.push_back(m_partition[i_snp]);
		}
		windows wind;
		std::get<prslice_wind::WIND>(wind) = window_name;
		std::get<prslice_wind::R2>(wind) = best_r2;
		std::get<prslice_wind::SNPS>(wind) = best_snp;
		std::get<prslice_wind::P>(wind) = std::get<+PRS::P>(m_best_threshold[0]);
		std::get<prslice_wind::NSNP>(wind) = std::get<+PRS::NSNP>(m_best_threshold[0]);
		std::get<prslice_wind::COEFF>(wind) = std::get<+PRS::COEFF>(m_best_threshold[0]);
		m_best_snps.push_back(wind);
	}


    		std::sort(begin(m_best_snps), end(m_best_snps),
    	              [](windows const &t1, windows const &t2)
    	    {
    			if(std::get<prslice_wind::R2>(t1)==std::get<prslice_wind::R2>(t2))
    			{
    				return std::get<prslice_wind::WIND>(t1).compare(std::get<prslice_wind::WIND>(t2))<0;
    			}
    			else return std::get<prslice_wind::R2>(t1) <std::get<prslice_wind::R2>(t2);
    	    }
    	             );
}

void PRSice::prslice(const Commander &c_commander, const Region &c_region, const size_t c_pheno_index)
{
	size_t n_thread = c_commander.get_thread();
	m_partition.clear();
	m_best_threshold.clear();
	m_num_snp_included.clear();
	m_best_score.clear(); //
	m_current_prs.clear(); //
	m_prs_results.clear(); //
	m_current_prs.push_back(m_sample_names);
	m_best_threshold.push_back(PRSice_best(0,0,0,0,0,0));
	m_prs_results.push_back(std::vector<PRSice_result>(0));
	m_best_score = m_current_prs;
	m_num_snp_included.push_back(0);
	size_t cur_index=0;
	size_t bin_count=0;
	for(auto &&snps : m_best_snps)
	{
		m_partition.insert(m_partition.end(), std::get<prslice_wind::SNPS>(snps).begin(), std::get<prslice_wind::SNPS>(snps).end());
		PLINK prs(m_target, m_chr_list);
		prs.initialize();
		// after this, the m_current_prs[i_region] will be updated
		prs.get_score(m_partition, m_snp_list, m_current_prs, cur_index, m_partition.size(), 0);
		cur_index=m_partition.size();
		m_num_snp_included.front() += m_partition.size();
		thread_score(0, 1, bin_count,n_thread, c_pheno_index);
		bin_count++;
	}
}

void PRSice::output(const Commander &c_commander, size_t pheno_index) const
{
	std::string pheno_name = std::get<pheno_store::NAME>(m_pheno_names[pheno_index]);
	std::string output_prefix = c_commander.get_out()+"."+m_current_base;
	if(!pheno_name.empty()) output_prefix.append("."+pheno_name+".");
	std::string out_best = output_prefix+".best";
	std::string out_prslice = output_prefix+".prslice";
	std::string out_wind = output_prefix+".windows";
	std::ofstream best_prs, prslice, window;
	best_prs.open(out_best.c_str());
	prslice.open(out_prslice.c_str());
	window.open(out_wind.c_str());
	if(!best_prs.is_open())
	{
		std::string error_message = "Cannot open file " +out_best+" for write";
		throw std::runtime_error(error_message);
	}
	if(!prslice.is_open())
	{
		std::string error_message = "Cannot open file " +out_prslice+" for write";
		throw std::runtime_error(error_message);
	}
	if(!window.is_open())
	{
		std::string error_message = "Cannot open file " +out_wind+" for write";
		throw std::runtime_error(error_message);
	}
	window << "Window\tR2\tCoefficient\tP\tNum_SNP" << std::endl;
	best_prs << "IID\tprs_"<<std::get<+PRS::THRESHOLD>(m_best_threshold[0]) << std::endl;
	prslice << "Num_Windows\tR2\tCoefficient\tP\tNum_SNP" << std::endl;
	for(auto &&wind : m_best_snps)
	{
		window << std::get<prslice_wind::WIND>(wind) << "\t"
				<< std::get<prslice_wind::R2>(wind)-m_null_r2 << "\t"
				<< std::get<prslice_wind::COEFF>(wind) << "\t"
				<< std::get<prslice_wind::P>(wind) << "\t"
				<< std::get<prslice_wind::NSNP>(wind) << "\t"
				<< std::endl;
	}
	window.close();
	int best_snp_size = std::get<+PRS::NSNP>(m_best_threshold[0]);
	if(best_snp_size==0)
	{
		fprintf(stderr, "ERROR: Best R2 obtained when no SNPs were included\n");
		fprintf(stderr, "       Cannot output the best PRS score\n");
	}
	else
	{
		for(auto &&prs : m_best_score[0])
		{
			best_prs << std::get<+PRS::IID>(prs) << "\t" <<
					std::get<+PRS::PRS>(prs)/(double)best_snp_size<< std::endl;
		}
	}
	best_prs.close();
	for(auto &&prs : m_prs_results[0])
	{
		prslice << std::get<+PRS::THRESHOLD>(prs) << "\t" <<
				std::get<+PRS::R2>(prs)-m_null_r2 << "\t" <<
				std::get<+PRS::COEFF>(prs)<< "\t" <<
				std::get<+PRS::P>(prs)<< "\t" <<
				std::get<+PRS::NSNP>(prs) << std::endl;
	}
	prslice.close();

}


void PRSice::update_line(std::unordered_map<std::string, size_t> &partition_index)
{
	std::ifstream bim;
	std::string bim_name = std::get<+PRS::FILENAME>(m_partition.front())+".bim";
	bim.open(bim_name.c_str());
	if(!bim.is_open())
	{
		std::string error_message = "Cannot open bim file " + bim_name+" for PRSlice";
		throw std::runtime_error(error_message);
	}
	std::string line;
	size_t cur_line = 0;
	while(std::getline(bim, line))
	{
		misc::trim(line);
		if(!line.empty())
		{
			std::vector<std::string> token = misc::split(line);
			if(token.size() < 6) throw std::runtime_error("Malformed bim file, should contain at least 6 columns");
			if(partition_index.find(token[+BIM::RS])!=partition_index.end())
			{
				std::get<+PRS::LINE>(m_partition[partition_index[token[+BIM::RS]]]) =cur_line;
			}
			cur_line++;
		}
	}
	bim.close();
}

