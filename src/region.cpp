/*
 * region.cpp
 *
 *  Created on: 25 Aug 2016
 *      Author: shingwanchoi
 */

#include "region.hpp"

Region::Region(std::vector<std::string> feature)
{
    m_bit_size = sizeof(long_type)*CHAR_BIT;
    // Make the base region which includes everything
    m_duplicated_names.insert("Base");
    m_region_name.push_back("Base");
    m_gtf_feature =feature;
    m_region_list.push_back(std::vector<boundary>(1));
}

void Region::run(const std::string &gtf, const std::string &msigdb, const std::vector<std::string> &bed,
		const std::string &out)
{
    process_bed(bed);


    std::unordered_map<std::string, std::set<std::string> > id_to_name;
    if(!gtf.empty())  // without the gtf file, we will not process the msigdb file
    {
    		fprintf(stderr, "Processing the GTF file\n");
        std::unordered_map<std::string, boundary > gtf_boundary;
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
}

void Region::process_bed(const std::vector<std::string> &bed)
{
	for(auto &b : bed)
	{
		fprintf(stderr, "Reading: %s\n", b.c_str());
		std::ifstream bed_file;
		bool error = false;
		bed_file.open(b.c_str());
		if(!bed_file.is_open()) fprintf(stderr, "WARNING: %s cannot be open. It will be ignored\n", b.c_str());
		else if(m_duplicated_names.find(b)==m_duplicated_names.end())
		{
			std::vector<boundary> current_region;
			std::string line;
			size_t num_line = 0;
			while(std::getline(bed_file, line))
			{
				num_line++;
				misc::trim(line);
				if(!line.empty())
				{
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
                    		temp = misc::convert<int>(token[+BOUNDARY::START]);
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
                    		temp = misc::convert<int>(token[+BOUNDARY::END]);
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
                    boundary cur_bound;
                    std::get<+BOUNDARY::CHR>(cur_bound) = token[0];
                    std::get<+BOUNDARY::START>(cur_bound) = start;
                    std::get<+BOUNDARY::END>(cur_bound) = end;
                    current_region.push_back(cur_bound);
				}
			}
			if(!error)
			{
				std::sort(begin(current_region), end(current_region),
						[](boundary const &t1, boundary const &t2)
						{
							if(std::get<+BOUNDARY::CHR>(t1).compare(std::get<+BOUNDARY::CHR>(t2))==0)
							{
								if(std::get<+BOUNDARY::START>(t1)==std::get<+BOUNDARY::START>(t2))
									return std::get<+BOUNDARY::END>(t1)<std::get<+BOUNDARY::END>(t2);
								return std::get<+BOUNDARY::START>(t1) < std::get<+BOUNDARY::START>(t2);
							}
							else return std::get<+BOUNDARY::CHR>(t1).compare(std::get<+BOUNDARY::CHR>(t2))<0;
						}
					);
				m_region_list.push_back(current_region);m_region_name.push_back(b);
				m_duplicated_names.insert(b);
			}
            bed_file.close();
		}
		else
		{
			fprintf(stderr, "%s is duplicated, it will be ignored\n", b.c_str());
		}
    }
}

std::unordered_map<std::string, boundary > Region::process_gtf(const std::string &gtf,
		std::unordered_map<std::string, std::set<std::string> > &id_to_name, const std::string &out_prefix)
{
    std::unordered_map<std::string, boundary > result_boundary;
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
        if(!line.empty() && line[0]!='#')
        {
            std::vector<std::string> token = misc::split(line, "\t");
            if(in_feature(token[+GTF::FEATURE]))
            {
                std::string chr = token[+GTF::CHR];
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
//                				WARNING: HARD CODING HERE
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
                		if(std::get<+BOUNDARY::CHR>(result_boundary[id]).compare(chr)!=0)
                		{
                			fprintf(stderr, "ERROR: Same gene occur on two separate chromosome!\n");
                			fprintf(stderr, "       Will ignore the gtf file\n");
                			result_boundary.clear();
                			id_to_name.clear();
                			return result_boundary;
                		}
                		if(std::get<+BOUNDARY::START>(result_boundary[id]) > start)
                				std::get<+BOUNDARY::START>(result_boundary[id]) = start;
                		if(std::get<+BOUNDARY::END>(result_boundary[id]) < end)
                				std::get<+BOUNDARY::END>(result_boundary[id]) = end;
                }
                else
                {
                		boundary cur_bound;
                		std::get<+BOUNDARY::CHR>(cur_bound) = chr;
                		std::get<+BOUNDARY::START>(cur_bound) = start;
                		std::get<+BOUNDARY::END>(cur_bound) = end;
                		result_boundary[id]=boundary(cur_bound);
                }
            }
        }
    }
    return result_boundary;
}

void Region::process_msigdb(const std::string &msigdb,
                            const std::unordered_map<std::string, boundary > &gtf_info,
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
            if(!line.empty())
            {
                std::vector<std::string> token = misc::split(line);
                if(token.size() < 2)  // Will treat the url as gene just in case
                {
                    fprintf(stderr, "Each line require at least 2 information\n");
                    fprintf(stderr, "%s\n", line.c_str());
                }
                else
                {
                    std::string name = token[0];
                    std::vector<boundary> current_region;
                    for(auto &gene: token)
                    {
                    		if(gtf_info.find(gene)==gtf_info.end())
                    		{
                    			if(id_to_name.find(gene)!= id_to_name.end())
                    			{
                    				auto name = id_to_name.at(gene);
                    				for(auto &translate: name)
                    				{
                    					if(gtf_info.find(translate) != gtf_info.end())
                    					{
                    						current_region.push_back(gtf_info.at(translate));
                    					}
                    				}
                    			}
                    			else
                    			{
                            		current_region.push_back(gtf_info.at(gene));
                    			}
                    		}
                    }
                    std::sort(begin(current_region), end(current_region),
                    		[](boundary const &t1, boundary const &t2)
							{
                        			if(std::get<+BOUNDARY::CHR>(t1).compare(std::get<+BOUNDARY::CHR>(t2))==0)
                        			{
                        				if(std::get<+BOUNDARY::START>(t1)==std::get<+BOUNDARY::START>(t2))
                        					return std::get<+BOUNDARY::END>(t1)<std::get<+BOUNDARY::END>(t2);
                        				return std::get<+BOUNDARY::START>(t1) < std::get<+BOUNDARY::START>(t2);
                        			}
                        			else return std::get<+BOUNDARY::CHR>(t1).compare(std::get<+BOUNDARY::CHR>(t2))<0;
							}
                    );
                }
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

std::vector<long_type> Region::check(std::string chr, size_t loc)
{
    std::vector<long_type> res = std::vector<long_type>(((m_region_name.size()+1)/m_bit_size)+1);
    res[0]=1; // base region which contains everything
    for(size_t i_region = 0; i_region < m_region_list.size(); ++i_region)
    {
        if(i_region==0)
        {
            res[0] |= ONE;
        }
        else
        {
            size_t current_region_size = m_region_list[i_region].size();
            while(m_snp_check_index[i_region]< current_region_size)
            {
                // do the checking
                boundary current_bound = m_region_list[i_region][m_snp_check_index[i_region]];
                std::string region_chr = std::get<+BOUNDARY::CHR>(current_bound);
                size_t region_start = std::get<+BOUNDARY::START>(current_bound);
                size_t region_end = std::get<+BOUNDARY::END>(current_bound);
                if(chr.compare(region_chr) != 0) m_snp_check_index[i_region]++;
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
	fprintf(stderr, "\nRegion Information\n");
	fprintf(stderr, "==============================\n");
	if (m_region_name.size() == 1)
		fprintf(stderr, "1 region is included\n");
	else if (m_region_name.size() > 1)
		fprintf(stderr, "A total of %zu regions are included\n", m_region_name.size());
}
