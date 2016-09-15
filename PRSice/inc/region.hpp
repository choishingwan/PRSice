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
#include <tuple>
#include <limits.h>
#include <string.h>
#include "misc.hpp"

class Region {
public:
	Region();
	virtual ~Region();
	void run(const std::string &gtf, const std::string &msigdb, const std::vector<std::string> &bed, const std::string &out, bool gen_bed);
	void reset(){ m_index = std::vector<size_t>(m_region_name.size());};
#if defined(__LP64__) || defined(_WIN64)
	typedef uint64_t long_type;
#else
	typedef uint32_t long_type;
#endif
	long_type* check(std::string chr, size_t loc);
    size_t flag_size() const { return (m_region_name.size()+1)/(m_bit_size)+1; };
    size_t size() const { return m_region_name.size(); };
    std::string get_name(size_t i) const { return m_region_name.at(i); };
private:
	typedef std::tuple<std::string, size_t, size_t> boundary;
	void process_bed(const std::vector<std::string> &bed);
	std::map<std::string, boundary > process_gtf(const std::string &gtf, std::map<std::string, std::string> &id_to_name, const std::string &out_prefix, bool gen_bed);
	void process_msigdb(const std::string &msigdb,
						const std::map<std::string, boundary > &gtf_info,
						const std::map<std::string, std::string> &id_to_name);
	std::vector<std::string> m_region_name;
	std::vector< std::vector<boundary> > m_region_list;
	// This is to indicate the current location onf each region
	// This work because we assume all SNPs are sorted by their coordinates
	// in the same way as the region files.
	std::vector<size_t> m_index;
	size_t m_bit_size;
};

#endif /* PRSICE_INC_REGION_HPP_ */
