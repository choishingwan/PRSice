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
#include <sys/stat.h>
#include <stdio.h>
#include <map>
#include "misc.hpp"

class Region {
public:
	Region();
	virtual ~Region();
	void run(const std::string &gtf, const std::string &msigdb, const std::vector<std::string> &bed, bool gen_bed);
	void reset(){ m_index = std::vector<size_t>(m_region_name.size());};
#if defined(__LP64__) || defined(_WIN64)
	uint64_t* check(std::string chr, size_t loc);
#else
	uint32_t* check(std::string chr, size_t loc);
#endif
private:
	typedef std::tuple<std::string, size_t, size_t> boundary;
	void process_bed(const std::vector<std::string> &bed);
	std::map<std::string, boundary > process_gtf(const std::string &gtf, std::map<std::string, std::string> &id_to_name, bool gen_bed);
	void process_msigdb(const std::string &msigdb,
						const std::map<std::string, boundary > &gtf_info,
						const std::map<std::string, std::string> &id_to_name);
	std::vector<std::string> m_region_name;
	std::vector< std::vector<boundary> > m_region_list;
	std::vector<size_t> m_index;
};

#endif /* PRSICE_INC_REGION_HPP_ */
