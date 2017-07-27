// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly
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


#include "region.hpp"

Region::Region(std::vector<std::string> feature, const std::unordered_map<std::string, int> &chr_order)
: m_chr_order(chr_order)
{
    m_bit_size = sizeof(long_type)*CHAR_BIT;
    // Make the base region which includes everything
    m_duplicated_names.insert("Base");
    m_region_name.push_back("Base");
    m_gtf_feature = feature;
    m_region_list.push_back(std::vector<region_bound>(1));
}

void Region::run(const std::string &gtf, const std::string &msigdb, const std::vector<std::string> &bed,
		const std::string &out)
{
    process_bed(bed);


    std::unordered_map<std::string, std::set<std::string> > id_to_name;
    if(!gtf.empty())  // without the gtf file, we will not process the msigdb file
    {
    	fprintf(stderr, "Processing the GTF file\n");
        std::unordered_map<std::string, region_bound > gtf_boundary;
        try
        {
        	gtf_boundary=process_gtf(gtf, id_to_name, out);
        }
        catch(const std::runtime_error &error)
        {
        	gtf_boundary.clear();
            fprintf(stderr, "ERROR: Cannot process GTF file: %s\n", error.what());
            fprintf(stderr, "       Will not process any of the msigdb items\n");
        }
        fprintf(stderr, "A total of %lu genes found in the GTF file\n", gtf_boundary.size());
        if(gtf_boundary.size() != 0)
        {
            process_msigdb(msigdb, gtf_boundary, id_to_name);
        }
    }
    m_snp_check_index = std::vector<size_t>(m_region_name.size());
    m_region_snp_count = std::vector<int>(m_region_name.size());

    m_duplicated_names.clear();
}

void Region::process_bed(const std::vector<std::string> &bed)
{
	for(auto &b : bed)
	{
		fprintf(stderr, "Reading: %s\n", b.c_str());
		std::ifstream bed_file;
		bool error = false;
		bed_file.open(b.c_str());
		if(!bed_file.is_open())
		{
			fprintf(stderr, "WARNING: %s cannot be open. It will be ignored\n", b.c_str());
			continue;
		}
		if(m_duplicated_names.find(b)!=m_duplicated_names.end())
		{
			fprintf(stderr, "%s is duplicated, it will be ignored\n", b.c_str());
			continue;
		}
		std::vector<region_bound> current_region;
		std::string line;
		size_t num_line = 0;
		while(std::getline(bed_file, line))
		{
			num_line++;
			misc::trim(line);
			if(line.empty()) continue;
			std::vector<std::string> token = misc::split(line);
			if(token.size() < 3)
			{
				fprintf(stderr, "ERROR: %s contain less than 3 column\n", b.c_str());
				fprintf(stderr, "       This file will be ignored\n");
				error=true;
				break;
			}
			int temp = 0;
			size_t start=0, end=0;
			try{
				temp = misc::convert<int>(token[1]);
				if(temp >= 0) start = temp+1;//That's because bed is 0 based
				else
				{
					fprintf(stderr, "ERROR: Negative Start Coordinate at line %zu!\n", num_line);
					fprintf(stderr, "       This file will be ignored\n");
					error=true;
					break;
				}
			}
			catch(const std::runtime_error &er){
				fprintf(stderr, "ERROR: Cannot convert start coordinate! (line: %zu)\n", num_line);
				fprintf(stderr, "       This file will be ignored\n");
				error=true;
				break;
			}
			try{
				temp = misc::convert<int>(token[2]);
				if(temp >= 0) end = temp+1;//That's because bed is 0 based
				else
				{
					fprintf(stderr, "ERROR: Negative End Coordinate at line %zu!\n", num_line);
					fprintf(stderr, "       This file will be ignored\n");
					error=true;
					break;
				}
			}
			catch(const std::runtime_error &er){
				fprintf(stderr, "ERROR: Cannot convert end coordinate! (line: %zu)\n", num_line);
				fprintf(stderr, "       This file will be ignored\n");
				error=true;
				break;
			}
			if(m_chr_order.find(token[0])!= m_chr_order.end())
			{
				region_bound cur_bound;
				cur_bound.chr = m_chr_order[token[0]];
				cur_bound.start = start;
				cur_bound.end = end;
				current_region.push_back(cur_bound);
			}
		}

		if(!error)
		{
			std::sort(begin(current_region), end(current_region),
					[](region_bound const &t1, region_bound const &t2)
					{
				if(t1.chr==t2.chr)
				{
					if(t1.start==t2.start)
						return t1.end < t2.end;
					return t1.start < t2.end;
				}
				else return t1.chr < t2.chr;
					}
			);
			m_region_list.push_back(current_region);
			m_region_name.push_back(b);
			m_duplicated_names.insert(b);
		}
		bed_file.close();
	}
}


std::unordered_map<std::string, Region::region_bound > Region::process_gtf(const std::string &gtf,
		std::unordered_map<std::string, std::set<std::string> > &id_to_name, const std::string &out_prefix)
{
    std::unordered_map<std::string, Region::region_bound > result_boundary;
    if(gtf.empty()) return result_boundary; // basically return an empty map

    std::ifstream gtf_file;
    gtf_file.open(gtf.c_str());
    if(!gtf_file.is_open())
    {
        std::string error_message = "Cannot open gtf file: "+gtf;
        throw std::runtime_error(error_message);
    }
    std::string line;
    size_t num_line = 0;
    while(std::getline(gtf_file, line))
    {
    	num_line++;
        misc::trim(line);
        if(line.empty() || line[0]=='#') continue ;
        std::vector<std::string> token = misc::split(line, "\t");
        std::string chr = token[+GTF::CHR];
        if(in_feature(token[+GTF::FEATURE]) && m_chr_order.find(chr)!=m_chr_order.end())
        {
        	int temp=0;
        	size_t start=0, end=0;
        	try
        	{
        		temp = misc::convert<int>(token[+GTF::START]);
        		if(temp < 0)
        		{
        			fprintf(stderr, "ERROR: Negative Start Coordinate! (line: %zu)\n", num_line);
        			fprintf(stderr, "       Will ignore the gtf file\n");
        			result_boundary.clear();
        			id_to_name.clear();
        			return result_boundary;
        			break;
        		}
        		start=temp;
        	}
        	catch(const std::runtime_error &er)
        	{
        		fprintf(stderr, "ERROR: Cannot convert the start coordinate! (line: %zu)\n", num_line);
        		fprintf(stderr, "       Will ignore the gtf file\n");
        		result_boundary.clear();
        		id_to_name.clear();
        		return result_boundary;
        	}
        	try
        	{
        		temp = misc::convert<int>(token[+GTF::END]);
        		if(temp < 0)
        		{
        			fprintf(stderr, "ERROR: Negative End Coordinate! (line: %zu)\n", num_line);
        			fprintf(stderr, "       Will ignore the gtf file\n");
        			result_boundary.clear();
        			id_to_name.clear();
        			return result_boundary;
        			break;
        		}
        		end=temp;
        	}
        	catch(const std::runtime_error &er)
        	{
        		fprintf(stderr, "ERROR: Cannot convert the end coordinate! (line: %zu)\n", num_line);
        		fprintf(stderr, "       Will ignore the gtf file\n");
        		result_boundary.clear();
        		id_to_name.clear();
        		return result_boundary;
        	}
                //Now extract the name
        	std::vector<std::string> attribute = misc::split(token[+GTF::ATTRIBUTE], ";");
        	std::string name="", id="";
        	for(auto &info : attribute)
        	{
        		if(info.find("gene_id")!=std::string::npos)
        		{
        			std::vector<std::string> extract = misc::split(info);
        			if(extract.size() > 1)
        			{
// WARNING: HARD CODING HERE cerr
        				extract[1].erase(std::remove(extract[1].begin(), extract[1].end(), '\"'), extract[1].end());
        				id = extract[1];
        			}
        		}
        		else if(info.find("gene_name")!=std::string::npos)
        		{
        			std::vector<std::string> extract = misc::split(info);
        			if(extract.size() > 1)
        			{
        				extract[1].erase(std::remove(extract[1].begin(), extract[1].end(), '\"'), extract[1].end());
        				name = extract[1];
        			}
        		}
        	}
        	if(!id.empty())
        	{
        		id_to_name[name].insert(id);
        	}
        	//Now add the information to the map using the id
        	if(result_boundary.find(id)!=result_boundary.end())
        	{
        		if(result_boundary[id].chr!=m_chr_order[chr])
        		{
        			fprintf(stderr, "ERROR: Same gene occur on two separate chromosome!\n");
        			fprintf(stderr, "       Will ignore the gtf file\n");
        			result_boundary.clear();
        			id_to_name.clear();
        			return result_boundary;
        		}
        		if(result_boundary[id].start > start)
        			result_boundary[id].start = start;
        		if(result_boundary[id].end < end)
        			result_boundary[id].end = end;
        	}
        	else
        	{
        		region_bound cur_bound;
        		cur_bound.chr = m_chr_order[chr];
        		cur_bound.start = start;
        		cur_bound.end = end;
        		result_boundary[id]=cur_bound;
        	}
        }
    }
    return result_boundary;
}

void Region::process_msigdb(const std::string &msigdb,
                            const std::unordered_map<std::string, Region::region_bound > &gtf_info,
                            const std::unordered_map<std::string, std::set<std::string> > &id_to_name)
{
    if(msigdb.empty() || gtf_info.size()==0) return; // Got nothing to do
    //Assume format = Name URL Gene
    std::ifstream input;
    input.open(msigdb.c_str());
    if(!input.is_open()) fprintf(stderr, "Cannot open %s. Will skip this file\n", msigdb.c_str());
    else
    {
        std::string line;
        while(std::getline(input, line))
        {
        	misc::trim(line);
        	if(line.empty()) continue;
        	std::vector<std::string> token = misc::split(line);
        	if(token.size() < 2)  // Will treat the url as gene just in case
        	{
        		fprintf(stderr, "Each line require at least 2 information\n");
        		fprintf(stderr, "%s\n", line.c_str());
        	}
        	else if(m_duplicated_names.find(token[0])==m_duplicated_names.end())
        	{
        		std::string name = token[0];
        		std::vector<region_bound> current_region;
        		for(auto &gene: token)
        		{
        			if(gtf_info.find(gene)==gtf_info.end())
        			{
        				if(id_to_name.find(gene)!= id_to_name.end())
        				{
        					auto &name = id_to_name.at(gene);
        					for(auto &&translate: name)
        					{
        						if(gtf_info.find(translate) != gtf_info.end())
        						{
        							current_region.push_back(gtf_info.at(translate));
        						}
        					}
        				}
        			}
        			else
        			{
        				current_region.push_back(gtf_info.at(gene));
        			}
        		}
        		std::sort(begin(current_region), end(current_region),
        				[](region_bound const &t1, region_bound const &t2)
						{
        			if(t1.chr == t2.chr)
        			{
        				if(t1.start==t2.start) return t1.end < t2.end;
        				return t1.start < t2.end;
        			}
        			else return t1.chr < t2.chr;
						}
        		);
        		m_region_list.push_back(current_region);
        		m_region_name.push_back(name);
        		m_duplicated_names.insert(name);
        	}
        	else if(m_duplicated_names.find(token[0])!=m_duplicated_names.end())
        	{
        		fprintf(stderr, "Duplicated Set: %s. It will be ignored\n", token[0].c_str());
        	}
        }
        input.close();
    }
}

void Region::print_file(std::string output) const
{
	std::ofstream region_out;
	region_out.open(output.c_str());
	if (!region_out.is_open()) {
		std::string error = "Cannot open region information file to write: " + output;
		throw std::runtime_error(error);
	}
	region_out << "Region\t#SNPs" << std::endl;
	for(size_t i_region = 0; i_region < m_region_name.size(); ++i_region)
	{
		region_out << m_region_name[i_region] << "\t" << m_region_snp_count[i_region] << std::endl;
	}
	region_out.close();
}

Region::~Region() {}

std::vector<long_type> Region::check(int chr, size_t loc)
{
    std::vector<long_type> res = std::vector<long_type>(((m_region_name.size()+1)/m_bit_size)+1);
    res[0]=1; // base region which contains everything
    for(size_t i_region = 0; i_region < m_region_list.size(); ++i_region)
    {
        if(i_region==0)
        {
            res[0] |= ONE;
            m_region_snp_count[0]++;
        }
        else
        {
            size_t current_region_size = m_region_list[i_region].size();
            while(m_snp_check_index[i_region]< current_region_size)
            {
                // do the checking
                auto &&current_bound = m_region_list[i_region][m_snp_check_index[i_region]];
                int region_chr = current_bound.chr;
                size_t region_start = current_bound.start;
                size_t region_end = current_bound.end;
                if(chr != region_chr) m_snp_check_index[i_region]++;
                else  // same chromosome
                {
                    if(region_start <= loc && region_end >=loc)
                    {
                        // This is the region
                        res[i_region/m_bit_size] |= ONE << i_region%m_bit_size;
                        m_region_snp_count[i_region]++;
                        break;
                    }
                    else if(region_start> loc) break;
                    else if(region_end < loc) m_snp_check_index[i_region]++;
                }
            }
        }
    }
    return res;
}

void Region::info() const{
	if (m_region_name.size() == 1)
		fprintf(stderr, "1 region included\n");
	else if (m_region_name.size() > 1)
		fprintf(stderr, "A total of %zu regions are included\n", m_region_name.size());
}
