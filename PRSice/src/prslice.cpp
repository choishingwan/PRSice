#include "prsice.hpp"

void prslice_windows(const Commander &c_commander, const Region &c_region)
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
		if(m_include_snp.find(snp.get_rs_id()))
		{
			std::string cur_chr = snp.get_chr();
			size_t cur_loc = snp.get_loc();
			if(prev_chr.empty() || cur_chr != prev_chr || cur_loc - prev_loc > wind_size)
			{
				//get line number, must be same file as same chromosome
				std::string window_name = cur_chr+":"+std::to_string(m_snp_list[std::get<+PRS::INDEX>(m_partition.front())].get_loc()) +
						"-"+std::to_string(m_snp_list[std::get<+PRS::INDEX>(m_partition.back())].get_loc());
				m_window_names.push_back(window_name);
				update_line(partition_index);
				// run prsice
				prsice(c_commander, c_region, true);
				for(auto &&threshold : m_best_threshold)
				{
					size_t best_snp_num = std::get<+PRS::NSNP>(threshold);
					double best_r2 = std::get<+PRS::R2>(threshold);
					for(size_t i_part=0; i_part <  best_snp_num; ++i_part)
					{
						//WORK HERE
					}
				}
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
		std::string window_name = m_snp_list[front_index].get_chr()+":"+std::to_string(m_snp_list[front_index].get_loc()) +
				"-"+std::to_string(m_snp_list[std::get<+PRS::INDEX>(m_partition.back())].get_loc());
		m_window_names.push_back(window_name);
		update_line(partition_index);
	}
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
			if(partition_index.find(token[+BIM::RS])!=partition_index.size())
			{
				std::get<+PRS::LINE>(m_partition[partition_index[token[+BIM::RS]]]]) =cur_line;
			}
			cur_line++;
		}
	}
	bim.close();
}
