/*
 * region.cpp
 *
 *  Created on: 25 Aug 2016
 *      Author: shingwanchoi
 */

#include "../inc/region.hpp"

Region::Region() {
	// TODO Auto-generated constructor stub

}

Region::~Region() {
	// TODO Auto-generated destructor stub
}

void Region::process_bed(const std::vector<std::string> &bed){
	for(size_t i = 0; i < bed.size(); ++i){
		// go through each bed file
		std::ifstream bed_file;
		bool error=false;
		bed_file.open(bed[i].c_str());
		if(!bed_file.is_open()) fprintf(stderr, "WARNING: %s cannot be open. It will be ignored\n", bed[i].c_str());
		else{
			m_region_name.push_back(bed[i]);
			std::vector<std::tuple<std::string, size_t, size_t> > current_region;
			std::string line;
			while(std::getline(bed_file, line)){
				misc::trim(line);
				if(!line.empty()){
					std::vector<std::string> token = misc::split(line);
					if(token.size() < 3){
						fprintf(stderr, "ERROR: %s contain less than 3 column\n", bed[i].c_str());
						fprintf(stderr, "       This file will be ignored\n");
						error=true;
						break;
					}
					int temp = 0;
					size_t start=0, end=0;
					try{
						temp =misc::convert<int>(token[1]);
						if(temp >= 0) start = temp+1;//That's because bed is 0 based
						else{
							fprintf(stderr, "ERROR: Negative Coordinates!\n");
							fprintf(stderr, "       This file will be ignored\n");
							error=true;
							break;
						}

					}
					catch(const std::runtime_error &error){
						fprintf(stderr, "ERROR: Cannot convert some of the coordinates!\n");
						fprintf(stderr, "       This file will be ignored\n");
						error=true;
						break;
					}
					try{
						temp =misc::convert<int>(token[2]); //because bed originally is exclusive
						if(temp >= 0) end = temp;
						else{
							fprintf(stderr, "ERROR: Negative Coordinates!\n");
							fprintf(stderr, "       This file will be ignored\n");
							error=true;
							break;
						}

					}
					catch(const std::runtime_error &error){
						fprintf(stderr, "ERROR: Cannot convert some of the coordinates!\n");
						fprintf(stderr, "       This file will be ignored\n");
						error=true;
						break;
					}
					current_region.push_back(std::tuple<std::string, size_t, size_t>(token[0], start, end));
				}
			}
			if(!error){
				m_region_list.push_back(current_region);
			}
			else{
				m_region_name.pop_back();
			}
			bed_file.close();
		}
	}
}

std::map<std::string, std::tuple<std::string, size_t, size_t> > process_gtf(const std::string &gtf){
	std::map<std::string, std::tuple<std::string, size_t, size_t> > res;
	if(gtf.empty()) return res; // basically return an empty map
	std::ifstream gtf_file;
	gtf_file.open(gtf.c_str());
	if(!gtf_file.is_open()){
		std::string error_message = "Cannot open gtf file: "+gtf;
		throw std::runtime_error(error_message);
	}
	std::string line;
	while(std::getline(gtf_file, line)){
		misc::trim(line);

	}
	gtf_file.close();
	return res;
}

void Region::process_msigdb(const std::vector<std::string> &msigdb,
							const std::map<std::string, std::tuple<std::string, size_t, size_t> > &gtf_info){

}
void Region::run(const std::string &gtf, const std::vector<std::string> &msigdb, const std::vector<std::string> &bed){
	process_bed(bed);
	if(!gtf.empty()){ // without the gtf file, we will not process the msigdb file
		std::map<std::string, std::tuple<std::string, size_t, size_t> > gtf_info;
		try{
			gtf_info=process_gtf(gtf);
		}
		catch(const std::runtime_error &error){
			fprintf(stderr, "ERROR: Cannot process gtf file: %s\n", error.what());
			fprintf(stderr, "       Will not process any of the msigdb items\n");
		}
		if(gtf_info.size() != 0){
			process_msigdb(msigdb, gtf_info);
		}
	}

}

