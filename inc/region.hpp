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
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <tuple>
#include <limits.h>
#include <string.h>
#include <utility>
#include "misc.hpp"
#include "storage.hpp"
#include <iostream>

class Region
{
public:
    Region(std::vector<std::string> feature);
    virtual ~Region();
    void run(const std::string &gtf, const std::string &msigdb, const std::vector<std::string> &bed, const std::string &out);
    void reset()
    {
        m_snp_check_index = std::vector<size_t>(m_region_name.size());
        m_region_snp_count = std::vector<int>(m_region_name.size());
    };

    std::vector<long_type> check(std::string chr, size_t loc);
    size_t flag_size() const
    {
        return (m_region_name.size()+1)/(m_bit_size)+1;
    };
    size_t size() const
    {
        return m_region_name.size();
    };
    std::string get_name(size_t i) const
    {
        return m_region_name.at(i);
    };

    std::vector<long_type> empty_flag()
    {
        std::vector<long_type> res = std::vector<long_type>(((m_region_name.size()+1)/m_bit_size)+1);
        res[0]=1; // base region which contains everything
        return res;
    }
    int get_count(size_t i) const { return m_region_snp_count.at(i); };
    void info() const;
    void print_file(std::string output) const;
    void prslice()
    {
    		m_region_name.clear();
    		m_region_name.push_back("Base");
    }
private:
    std::unordered_set<std::string> m_duplicated_names;
    std::vector<std::string> m_region_name;
    std::vector<std::string> m_gtf_feature;
    std::vector< std::vector<boundary> > m_region_list;
    // This is to indicate the current location of each region
    // This work because we assume all SNPs are sorted by their coordinates
    // in the same way as the region files.
    std::vector<size_t> m_snp_check_index;
    std::vector<int> m_region_snp_count;

    bool in_feature(std::string in) const
    {
    		for(auto &feature: m_gtf_feature)
    		{
    			if(in.compare(feature)==0) return true;
    		}
    		return false;
    }




    void process_bed(const std::vector<std::string> &bed);

    std::unordered_map<std::string, boundary > process_gtf(const std::string &gtf,
    		std::unordered_map<std::string, std::set<std::string> > &id_to_name, const std::string &out_prefix);

    void process_msigdb(const std::string &msigdb,
                        const std::unordered_map<std::string, boundary > &gtf_info,
                        const std::unordered_map<std::string, std::set<std::string> > &id_to_name);


    size_t m_bit_size;
};

#endif /* PRSICE_INC_REGION_HPP_ */
