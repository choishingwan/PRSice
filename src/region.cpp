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

Region::Region(std::vector<std::string> feature,
               const std::unordered_map<std::string, int>& chr_order, const int window_5,
			   const int window_3)
    : m_chr_order(chr_order), m_5prime(window_5), m_3prime(window_3)
{
    // Make the base region which includes everything
    m_duplicated_names.insert("Base");
    m_region_name.push_back("Base");
    m_gtf_feature = feature;
    m_region_list.push_back(std::vector<region_bound>(1));
}

void Region::run(const std::string& gtf, const std::string& msigdb,
                 const std::vector<std::string>& bed, const std::string& out, const std::string &background, Reporter &reporter)
{
	if(gtf.empty() && bed.size()==0){
	    m_snp_check_index = std::vector<size_t>(m_region_name.size(), 0);
	    m_region_snp_count = std::vector<int>(m_region_name.size());
		return;
	}
	if((m_5prime >0 ||  m_3prime > 0) && (m_5prime!=m_3prime)){
		std::string message = "Warning: We will assume a positive strand for any features with unspecific strand information e.g. \".\"";
		reporter.report(message);
	}
    process_bed(bed, reporter);
    m_out_prefix = out;
    size_t num_bed_region = m_region_list.size();
    std::unordered_map<std::string, std::set<std::string>> id_to_name;
    std::unordered_map<std::string, region_bound> gtf_boundary;
    if (!gtf.empty()) // without the gtf file, we will not process the msigdb
                      // file
    {
        reporter.report("Processing the GTF file");
        try
        {
            gtf_boundary = process_gtf(gtf, id_to_name, out, reporter);
        }
        catch (const std::runtime_error& error)
        {
            gtf_boundary.clear();
            throw std::runtime_error(error.what());
        }
        std::string message = "A total of "+std::to_string(gtf_boundary.size())+" genes found in the GTF file";
        reporter.report(message);
        if (gtf_boundary.size() != 0) {
            process_msigdb(msigdb, gtf_boundary, id_to_name, reporter);
        }
    }
    if(background.empty()){
    	generate_background(gtf_boundary, num_bed_region, reporter);
    }
    else{
    	// what is a good way to determine the input type for background?
    	read_background(background, gtf_boundary, id_to_name, reporter);
    }
    m_snp_check_index = std::vector<size_t>(m_region_name.size(), 0);
    m_region_snp_count = std::vector<int>(m_region_name.size());

    m_duplicated_names.clear();
}


std::vector<Region::region_bound> Region::solve_overlap(std::vector<Region::region_bound> &current_region){
    std::sort(begin(current_region), end(current_region),
              [](region_bound const& t1, region_bound const& t2) {
                  if (t1.chr == t2.chr) {
                      if (t1.start == t2.start) return t1.end < t2.end;
                      return t1.start < t2.end;
                  }
                  else
                      return t1.chr < t2.chr;
              });
    std::vector<Region::region_bound> result;
    int prev_chr = -1;
    size_t prev_start = 0;
    size_t prev_end = 0;
    for(auto && bound  : current_region){
    	if(prev_chr==-1){
    		prev_chr = bound.chr;
    		prev_start = bound.start;
    		prev_end = bound.end;
    	}else if(prev_chr != bound.chr || bound.start > prev_end){
    		// new region
    		region_bound cur_bound;
    		cur_bound.chr = prev_chr;
    		cur_bound.start = prev_start;
    		cur_bound.end = prev_end;
    		result.push_back(cur_bound);
    		prev_chr = bound.chr;
    		prev_start = bound.start;
    		prev_end = bound.end;
    	}else{
    		prev_end =bound.end;
    	}
    }
    if(prev_chr != -1){
    	region_bound cur_bound ;
    	cur_bound.chr= prev_chr;
    	cur_bound.start = prev_start;
    	cur_bound.end = prev_end;
    	result.push_back(cur_bound);
    }
    return result;
}
void Region::process_bed(const std::vector<std::string>& bed, Reporter &reporter)
{
	// TODO: Allow user define name by modifying their input (e.g. Bed:Name
	bool print_warning = false;
    for (auto& b : bed) {
    	std::string message = "Reading: "+b;
    	reporter.report(message);
        std::ifstream bed_file;
        bool error = false;
        bed_file.open(b.c_str());
        if (!bed_file.is_open()) {
        	message = "Warning: "+b+" cannot be open. It will be ignored";
        	reporter.report(message);
            continue;
        }
        if (m_duplicated_names.find(b) != m_duplicated_names.end()) {
        	message = "Warning: "+b+" is duplicated, it will be ignored";
        	reporter.report(message);
            continue;
        }
        std::vector<region_bound> current_region;
        std::string line;
        size_t num_line = 0;
        while (std::getline(bed_file, line)) {
            num_line++;
            misc::trim(line);
            if (line.empty()) continue;
            std::vector<std::string> token = misc::split(line);
            if (token.size() < 3) {

            	message = "Error: "+b+" contain less than 3 columns, it will be ignored";
            	reporter.report(message);
            	break;
            }
            if(token.size() <= +BED::STRAND && (m_5prime >0 ||  m_3prime > 0) && (m_5prime!=m_3prime) && !print_warning){
            		std::string message = "Warning: You bed file does not contain strand information, we will assume all regions are on the positive strand, e.g. start coordinates always on the 5' end";
            		reporter.report(message);
            		print_warning = true;
            }
            int temp = 0;
            size_t start = 0, end = 0;
            message="";
            try
            {
                temp = misc::convert<int>(token[+BED::START]);
                if (temp >= 0)
                    start = temp + 1; // That's because bed is 0 based
                else
                {

                	message.append("Error: Negative Start Coordinate at line "+std::to_string(num_line)+"!");
                    error = true;
                }
            }
            catch (const std::runtime_error& er)
            {
            	message.append("Error: Cannot convert start coordinate! (line: "+std::to_string(num_line)+")!");
                error = true;
            }
            try
            {
                temp = misc::convert<int>(token[+BED::END]);
                if (temp >= 0)
                    end = temp + 1; // That's because bed is 0 based
                else
                {
                	message.append("Error: Negative End Coordinate at line "+std::to_string(num_line)+"!");
                    error = true;
                }
            }
            catch (const std::runtime_error& er)
            {
            	message.append("Error: Cannot convert end coordinate! (line: "+std::to_string(num_line)+")!");
                error = true;
            }
            if(start > end){
            	error = true;
            	message.append("Error: Start coordinate should be smaller than end coordinate!\n");
            	message.append("start: "+std::to_string(start)+"\n");
            	message.append("end: "+std::to_string(end)+"\n");
            }
            if(error) break;
            // only include regions that falls into the chromosome of interest
            if (m_chr_order.find(token[+BED::CHR]) != m_chr_order.end()) {
            	if(token.size() > +BED::STRAND){
                    if(token[+BED::STRAND].compare("-")==0){
                    	if(start-m_3prime <1){
                    		start = 1;
                    	}
                    	else start -= m_3prime;
                    	end+=m_5prime;
                    }
                    else if(token[+BED::STRAND].compare("+")==0 || token[+BED::STRAND].compare(".")==0){
                    	if(start-m_5prime < 1){
                    		start = 1;
                    	}else start -=m_5prime;
                    	end+=m_3prime;
                    }else{
                    	std::string error = "Error: Undefined strand information. Possibly a malform BED file: "+token[+BED::STRAND];
                        throw std::runtime_error(error);
                    }
            	}
            	else{
					if(start-m_5prime<1){
						start = 1;
					}
					else start -= m_5prime;
					end += m_3prime;
				}
                region_bound cur_bound;
                cur_bound.chr = m_chr_order[token[0]];
                cur_bound.start = start;
                cur_bound.end = end;
                // this should help us to avoid problem
                current_region.push_back(cur_bound);
            }
        }

        if (!error) {
        	// TODO: DEBUG This!! Might go out of scope
            m_region_list.push_back(solve_overlap(current_region));
            m_region_name.push_back(b);
            m_duplicated_names.insert(b);
        }
        else{
        	throw std::runtime_error(message);
        }
        bed_file.close();
    }
}


std::unordered_map<std::string, Region::region_bound> Region::process_gtf(
    const std::string& gtf,
    std::unordered_map<std::string, std::set<std::string>>& id_to_name,
    const std::string& out_prefix,
	Reporter &reporter)
{

    std::unordered_map<std::string, Region::region_bound> result_boundary;
    if (gtf.empty()) return result_boundary; // basically return an empty map

    std::string line;
    size_t num_line = 0, exclude_feature = 0;


    bool gz_input = false;
    GZSTREAM_NAMESPACE::igzstream gz_gtf_file;
    if (gtf.substr(gtf.find_last_of(".") + 1).compare("gz") == 0) {
        gz_gtf_file.open(gtf.c_str());
        if (!gz_gtf_file.good()) {
            std::string error_message =
                "Error: Cannot open GTF (gz) to read!\n";
            throw std::runtime_error(error_message);
        }
        gz_input = true;
    }

    std::ifstream gtf_file;
    if (!gz_input) {
        gtf_file.open(gtf.c_str());
        if (!gtf_file.is_open()) {
            std::string error_message = "Cannot open gtf file: " + gtf;
            throw std::runtime_error(error_message);
        }
    }

    while ((!gz_input && std::getline(gtf_file, line))
           || (gz_input && std::getline(gz_gtf_file, line)))
    {
        num_line++;
        misc::trim(line);
        // skip headers
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> token = misc::split(line, "\t");
        std::string chr = token[+GTF::CHR];
        if (in_feature(token[+GTF::FEATURE])
            && m_chr_order.find(chr) != m_chr_order.end())
        {
            int temp = 0;
            size_t start = 0, end = 0;
            try
            {
                temp = misc::convert<int>(token[+GTF::START]);
                if (temp < 0) {
                	// well, this is a bit extreme. Alternatively we can opt for
                	// just skipping problematic entries
                    std::string error =  "Error: Negative Start Coordinate! (line: "+std::to_string(num_line)+"). Will ignore the gtf file\n";
                    result_boundary.clear();
                    id_to_name.clear();
                    throw std::runtime_error(error);
                }
                start = temp;
            }
            catch (const std::runtime_error& er)
            {
                std::string error = "Error: Cannot convert the start coordinate! (line: "+std::to_string(num_line)+"). Will ignore the gtf file";
                result_boundary.clear();
                id_to_name.clear();
                throw std::runtime_error(error);
            }
            try
            {
                temp = misc::convert<int>(token[+GTF::END]);
                if (temp < 0) {
                	std::string error =  "Error: Negative End Coordinate! (line: "+std::to_string(num_line)+"). Will ignore the gtf file\n";
                    result_boundary.clear();
                    id_to_name.clear();
                    throw std::runtime_error(error);
                }
                end = temp;
            }
            catch (const std::runtime_error& er)
            {
            	std::string error = "Error: Cannot convert the end coordinate! (line: "+std::to_string(num_line)+"). Will ignore the gtf file";
                result_boundary.clear();
                id_to_name.clear();
                throw std::runtime_error(error);
            }
            // Now extract the name
            std::vector<std::string> attribute =
                misc::split(token[+GTF::ATTRIBUTE], ";");
            std::string name = "", id = "";
            for (auto& info : attribute) {
                if (info.find("gene_id") != std::string::npos) {
                    std::vector<std::string> extract = misc::split(info);
                    if (extract.size() > 1) {
                        // TODO: WARNING: HARD CODING HERE
                        extract[1].erase(std::remove(extract[1].begin(),
                                                     extract[1].end(), '\"'),
                                         extract[1].end());
                        id = extract[1];
                    }
                }
                else if (info.find("gene_name") != std::string::npos)
                {
                    std::vector<std::string> extract = misc::split(info);
                    if (extract.size() > 1) {
                        extract[1].erase(std::remove(extract[1].begin(),
                                                     extract[1].end(), '\"'),
                                         extract[1].end());
                        name = extract[1];
                    }
                }
            }
            // the GTF only mandate the ID field, so GTF can miss out the name field.
            if (!id.empty() && !name.empty()) {
                id_to_name[name].insert(id);
            }

            if(start > end){
            	std::string message="Error: Start coordinate should be smaller than end coordinate!\n";
            	message.append("start: "+std::to_string(start)+"\n");
            	message.append("end: "+std::to_string(end)+"\n");
                result_boundary.clear();
                id_to_name.clear();
                throw std::runtime_error(message);

            }
            if(token[+GTF::STRAND].compare("-")==0){
            	if(start-m_3prime <1){
            		start = 1;
            	}
            	else start -= m_3prime;
            	end+=m_5prime;
            }
            else if(token[+GTF::STRAND].compare("+")==0 || token[+GTF::STRAND].compare(".")==0){
            	if(start-m_5prime < 1){
            		start = 1;
            	}else start -=m_5prime;
            	end+=m_3prime;
            }else{
            	std::string error = "Error: Undefined strand information. Possibly a malform GTF file: "+token[+GTF::STRAND];
                result_boundary.clear();
                id_to_name.clear();
                throw std::runtime_error(error);
            }
            // Now add the information to the map using the id
            if (result_boundary.find(id) != result_boundary.end()) {
                if (result_boundary[id].chr != m_chr_order[chr]) {
                	std::string error = "Error: Same gene occur on two separate chromosome!";
                    result_boundary.clear();
                    id_to_name.clear();
                    throw std::runtime_error(error);
                }
                if (result_boundary[id].start > start)
                    result_boundary[id].start = start;
                if (result_boundary[id].end < end)
                    result_boundary[id].end = end;
            }
            else
            {
                region_bound cur_bound;
                cur_bound.chr = m_chr_order[chr];
                cur_bound.start = start;
                cur_bound.end = end;
                result_boundary[id] = cur_bound;

            }
        }
        else
        {
            exclude_feature++;
        }
    }
    std::string message ="";
    if (exclude_feature == 1){
        message.append("A total of "+std::to_string(exclude_feature)+" entry removed due to feature selection");
    	reporter.report(message);
    }
    if (exclude_feature > 1){
    	message.append("A total of "+std::to_string(exclude_feature)+" entries removed due to feature selection");
    	reporter.report(message);
    }

    return result_boundary;
}

void Region::process_msigdb(
    const std::string& msigdb,
    const std::unordered_map<std::string, Region::region_bound>& gtf_info,
    const std::unordered_map<std::string, std::set<std::string>>& id_to_name,
	Reporter &reporter)
{
    if (msigdb.empty() || gtf_info.size() == 0) return; // Got nothing to do
    // Assume format = Name URL Gene
    std::ifstream input;
    // in theory, it should be easy for us to support multiple msigdb file.
    // but ignore that for now TODO
    input.open(msigdb.c_str());
    if (!input.is_open())
    	reporter.report("Cannot open "+msigdb+". Will skip this file");
    else
    {
        std::string line;
        while (std::getline(input, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            std::vector<std::string> token = misc::split(line);
            if (token.size() < 2) // Will treat the url as gene just in case
            {
            	std::string message = "Error: Each line require at least 2 information\n";
            	message.append(line);
            	reporter.report(message);
            }
            else if (m_duplicated_names.find(token[0])
                     == m_duplicated_names.end())
            {
                std::string name = token[0];
                std::vector<region_bound> current_region;
                for (auto& gene : token) {
                    if (gtf_info.find(gene) == gtf_info.end()) {
                    	// we cannot find the gene name in the gtf information (which uses gene ID)
                        if (id_to_name.find(gene) != id_to_name.end()) {
                        	// we found a way to convert the gene name to gene id
                            auto& name = id_to_name.at(gene);
                            // problem is, one gene name can correspond to multiple gene id
                            // in that case, wee will take all of them
                            for (auto&& translate : name) {
                                if (gtf_info.find(translate) != gtf_info.end())
                                {
                                    current_region.push_back(
                                        gtf_info.at(translate));
                                }
                            }
                        }
                    }
                    else
                    {
                        current_region.push_back(gtf_info.at(gene));
                    }
                }
                m_region_list.push_back(solve_overlap(current_region));
                m_region_name.push_back(name);
                m_duplicated_names.insert(name);
            }
            else if (m_duplicated_names.find(token[0])
                     != m_duplicated_names.end())
            {
            	reporter.report("Duplicated Set: "+token[0]+". It will be ignored");
            }
        }
        input.close();
    }
}

void Region::generate_background(const std::unordered_map<std::string, Region::region_bound> &gtf_info, const size_t num_bed_region, Reporter &reporter){
	// this will be very ineffective. But whatever
	std::vector<Region::region_bound> temp_storage;
	for(auto &&gtf: gtf_info){
		temp_storage.push_back(gtf.second);
	}
	for(size_t i = 0; i< num_bed_region; ++i){
		for(size_t j = 0; j < m_region_list[i].size(); ++j){
			temp_storage.push_back(m_region_list[i][j]);
		}
	}
	m_region_list.push_back(solve_overlap(temp_storage));
    m_region_name.push_back("Background");
}

void Region::read_background(const std::string &background,
        const std::unordered_map<std::string, region_bound>& gtf_info,
        const std::unordered_map<std::string, std::set<std::string>>&
            id_to_name, Reporter &reporter){
    std::unordered_map<std::string, int> file_type{{"bed", 1},
                                                   {"range", 0},
                                                   {"gene", 2}};
    std::vector<std::string> background_info = misc::split(background, ":");
    if(background_info.size() != 2){
    	std::string error = "Error: Format of --background should be <File Name>:<File Type>";
    	throw std::runtime_error(error);
    }
    auto type = file_type.find(background_info[1]);
    if(type == file_type.end()){
    	std::string error = "Error: Undefined file type. Supported formats are bed, gene or range";
    	throw std::runtime_error(error);
    }
    std::ifstream input;
    input.open(background_info[0].c_str());
    if(!input.is_open()){
    	std::string error = "Error: Cannot open background file: "+background_info[0];
    	throw std::runtime_error(error);
    }
    bool print_warning = false, error = false;
    std::vector<Region::region_bound> current_bound;
    std::string line;
    if(type==0 || type==1){
    	// range or bed
        size_t num_line = 0;
    	while(std::getline(input, line)){
    		num_line++;
    		misc::trim(line);
    		if (line.empty()) continue;
    		std::vector<std::string> token = misc::split(line);
    		if (token.size() < 3) {
    			std::string message = "Error: "+background_info[0]+" contain less than 3 columns, it will be ignored";
    			throw std::runtime_error(message);
    		}
    		if(token.size() <= +BED::STRAND && (m_5prime >0 ||  m_3prime > 0) && (m_5prime!=m_3prime) && !print_warning){
    			std::string message = "Warning: You bed file does not contain strand information, we will assume all regions are on the positive strand, e.g. start coordinates always on the 5' end";
    			reporter.report(message);
    			print_warning = true;
    		}
    		int temp = 0;
    		size_t start = 0, end = 0;
    		std::string message="";
    		try
    		{
    			temp = misc::convert<int>(token[+BED::START]);
    			if (temp >= 0)
    				start = temp + type; // That's because bed is 0 based and range format is 1 based
    			else
    			{
    				message.append("Error: Negative Start Coordinate at line "+std::to_string(num_line)+"!");
    				error = true;
    			}
    		}
    		catch (const std::runtime_error& er)
    		{
    			message.append("Error: Cannot convert start coordinate! (line: "+std::to_string(num_line)+")!");
    			error = true;
    		}
    		try
    		{
    			temp = misc::convert<int>(token[+BED::END]);
    			if (temp >= 0)
    				end = temp + 1; // That's because bed is 0 based
    			else
    			{
    				message.append("Error: Negative End Coordinate at line "+std::to_string(num_line)+"!");
    				error = true;
    			}
    		}
    		catch (const std::runtime_error& er)
    		{
    			message.append("Error: Cannot convert end coordinate! (line: "+std::to_string(num_line)+")!");
    			error = true;
    		}
    		if(start > end){
    			error = true;
    			message.append("Error: Start coordinate should be smaller than end coordinate!\n");
    			message.append("start: "+std::to_string(start)+"\n");
    			message.append("end: "+std::to_string(end)+"\n");
    		}
    		if(error) break;
    		// only include regions that falls into the chromosome of interest
    		if (m_chr_order.find(token[+BED::CHR]) != m_chr_order.end()) {
    			int strand_index = (type==1)?(+BED::STRAND) : (+BED::END+1);
    			if(token.size() > strand_index){
    				if(token[strand_index].compare("-")==0){
    					if(start-m_3prime <1){
    						start = 1;
    					}
    					else start -= m_3prime;
    					end+=m_5prime;
    				}
    				else if(token[strand_index].compare("+")==0 || token[strand_index].compare(".")==0){
    					if(start-m_5prime < 1){
    						start = 1;
    					}else start -=m_5prime;
    					end+=m_3prime;
    				}else{
    					std::string error = "Error: Undefined strand information. Possibly a malform BED file: "+token[+BED::STRAND];
    					throw std::runtime_error(error);
    				}
    			}
    			else{
    				if(start-m_5prime<1){
    					start = 1;
    				}
    				else start -= m_5prime;
    				end += m_3prime;
    			}
    			region_bound cur_bound;
    			cur_bound.chr = m_chr_order[token[0]];
    			cur_bound.start = start;
    			cur_bound.end = end;
    			// this should help us to avoid problem
    			current_bound.push_back(cur_bound);
    		}
    	}
    }else{
    	// gene list format

        while (std::getline(input, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            if (gtf_info.find(line) == gtf_info.end()) {
            	// we cannot find the gene name in the gtf information (which uses gene ID)
            	if (id_to_name.find(line) != id_to_name.end()) {
            		// we found a way to convert the gene name to gene id
            		auto& name = id_to_name.at(line);
            		// problem is, one gene name can correspond to multiple gene id
            		// in that case, wee will take all of them
            		for (auto&& translate : name) {
            			if (gtf_info.find(translate) != gtf_info.end())
            			{
            				current_bound.push_back(
            						gtf_info.at(translate));
            			}
            		}
            	}
            }
            else
            {
            	current_bound.push_back(gtf_info.at(line));
            }
        }
    }
    input.close();
	m_region_list.push_back(solve_overlap(current_bound));
    m_region_name.push_back("Background");
}


void Region::print_file(std::string output) const
{
    std::ofstream region_out;
    region_out.open(output.c_str());
    if (!region_out.is_open()) {
        std::string error =
            "Cannot open region information file to write: " + output;
        throw std::runtime_error(error);
    }
    region_out << "Region\t#SNPs" << std::endl;
    for (size_t i_region = 0; i_region < m_region_name.size(); ++i_region) {
        region_out << m_region_name[i_region] << "\t"
                   << m_region_snp_count[i_region] << std::endl;
    }
    region_out.close();
}

Region::~Region() {}

void Region::check(std::string chr, size_t loc, std::vector<uintptr_t>& flag)
{
    flag[0] |= ONELU;
    m_region_snp_count[0]++;
    if (m_chr_order.find(chr) == m_chr_order.end())
        return; // chromosome not found
    int chr_index = m_chr_order[chr];
    // note: the chr is actually the order on the m_chr_order instead of the
    // sactual chromosome
    for (size_t i_region = 1; i_region < m_region_name.size(); ++i_region) {
        size_t current_region_size = m_region_list[i_region].size();
        while (m_snp_check_index[i_region] < current_region_size) {
            auto&& current_bound =
                m_region_list[i_region][m_snp_check_index[i_region]];
            int region_chr = current_bound.chr;
            size_t region_start = current_bound.start;
            size_t region_end = current_bound.end;
            if (chr_index
                > region_chr) // only increment if we have passed the chromosome
            {
                m_snp_check_index[i_region]++;
            }
            else if (chr_index == region_chr) // same chromosome
            {
                if (region_start <= loc && region_end >= loc) {
                    // This is the region
                    flag[i_region / BITCT] |= ONELU << ((i_region) % BITCT);
                    m_region_snp_count[i_region]++;
                    break;
                }
                else if (region_start > loc)
                    break;
                else if (region_end < loc)
                {
                    m_snp_check_index[i_region]++;
                }
            }
            else
            {
                // not the same chromosome
                break;
            }
        }
    }
}

void Region::info(Reporter& reporter) const
{
    std::string message = "";
    if (m_region_name.size() == 1) {
        message = "1 region included";
    }
    else if (m_region_name.size() > 1)
    {
        message = "A total of " + std::to_string(m_region_name.size())
                  + " regions are included";
    }
    reporter.report(message);
}
