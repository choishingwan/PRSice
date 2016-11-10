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

#include <iostream>

class Region
{
public:
    Region();
    virtual ~Region();
    void run(const std::string &gtf, const std::string &msigdb, const std::vector<std::string> &bed, const std::string &out, bool gen_bed);
    void reset()
    {
        m_index = std::vector<size_t>(m_region_name.size());
        m_region_count = std::vector<int>(m_region_name.size());
    };
#if defined(__LP64__) || defined(_WIN64)
    typedef uint64_t long_type;
#define ONE  0x1LLU
#else
    typedef uint32_t long_type;
#define ONE  0x1LU
#endif
    long_type* check(std::string chr, size_t loc);
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

    Region::long_type* empty_flag()
    {
        long_type* res = new long_type[((m_region_name.size()+1)/m_bit_size)+1];
        memset(res, 0x0,(((m_region_name.size()+1)/m_bit_size)+1)*sizeof(long_type));
        res[0]=1; // base region which contains everything
        return res;
    }
    std::vector<std::pair<std::string, double> > get_info() const {
    		return m_processed_regions;
    }
    int get_count(size_t i) const { return m_region_count.at(i); };

private:
    typedef std::tuple<std::string, size_t, size_t> boundary;
    void process_bed(const std::vector<std::string> &bed);
    std::unordered_map<std::string, boundary > process_gtf(const std::string &gtf,
    		std::unordered_map<std::string, std::set<std::string> > &id_to_name, const std::string &out_prefix, bool gen_bed);
    void process_msigdb(const std::string &msigdb,
                        const std::unordered_map<std::string, boundary > &gtf_info,
                        const std::unordered_map<std::string, std::set<std::string> > &id_to_name);
    std::unordered_set<std::string> m_dup_names;
    std::vector<std::string> m_region_name;
    std::vector< std::vector<boundary> > m_region_list;
    std::vector<std::pair<std::string, double> > m_processed_regions;
    std::vector<int> m_region_count;
    // This is to indicate the current location onf each region
    // This work because we assume all SNPs are sorted by their coordinates
    // in the same way as the region files.
    std::vector<size_t> m_index;
    size_t m_bit_size;
};

#endif /* PRSICE_INC_REGION_HPP_ */
