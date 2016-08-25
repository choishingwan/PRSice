/*
 * region.h
 *
 *  Created on: 25 Aug 2016
 *      Author: shingwanchoi
 */

#ifndef PRSICE_INC_REGION_HPP_
#define PRSICE_INC_REGION_HPP_

#include <string>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <map>
#include "misc.hpp"

class Region {
public:
	Region();
	virtual ~Region();
	void run(const std::string &gtf, const std::vector<std::string> &msigdb, const std::vector<std::string> &bed);
private:
	void process_bed(const std::vector<std::string> &bed);
	std::map<std::string, std::tuple<std::string, size_t, size_t> > process_gtf(const std::string &gtf);
	void process_msigdb(const std::vector<std::string> &msigdb,
						const std::map<std::string, std::tuple<std::string, size_t, size_t> > &gtf_info);
	std::vector<std::string> m_region_name;
	std::vector< std::vector<std::tuple<std::string, size_t, size_t> > m_region_list;
};

#endif /* PRSICE_INC_REGION_HPP_ */
