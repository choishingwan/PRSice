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

#ifndef REGION_H
#define REGION_H

#include "gzstream.h"
#include "misc.hpp"
#include "plink_common.hpp"
#include "reporter.hpp"
#include "storage.hpp"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <map>
#include <set>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

class Region
{
public:
	Region(const std::string &exclusion_range, Reporter& reporter);
    Region(std::vector<std::string> feature,
           const int window_5, const int window_3);
    virtual ~Region();
    void run(const std::string& gtf, const std::string& msigdb,
             const std::vector<std::string>& bed, const std::string& out,
             const std::string& background, Reporter& reporter);
    void reset()
    {
        m_snp_check_index = std::vector<size_t>(m_region_name.size());
        m_region_snp_count = std::vector<int>(m_region_name.size());
    };

    void update_flag(std::string chr, size_t loc, std::vector<uintptr_t>& flag);
    size_t size() const { return m_region_name.size(); };
    std::string get_name(size_t i) const { return m_region_name.at(i); };
    std::vector<std::string> names() const { return m_region_name; };
    int get_count(size_t i) const { return m_region_snp_count.at(i); };
    void info(Reporter& reporter) const;
    void print_file(std::string output) const;
    void prslice()
    {
        m_region_name.clear();
        m_region_name.push_back("Base");
    }
    void clean()
    {
        m_region_list = std::vector<std::vector<region_bound>>();
        m_snp_check_index = std::vector<size_t>();
    }
    void post_clump_count(std::vector<int>& count)
    {
        int max = 0;
        m_region_post_clump_count.resize(count.size());
        for (size_t i = 0; i < count.size(); ++i) {
            max = (count[i] > max && i != count.size() - 1 && i != 0) ? count[i]
                                                                      : max;
            m_region_post_clump_count[i] = count[i];
            if (m_region_size_duplicated.find(count[i])
                != m_region_size_duplicated.end())
                m_region_size_duplicated[count[i]] = true;
            else
                m_region_size_duplicated[count[i]] = false;
        }
        // we won't do background with the base group
        if (count.size() > 1 && max > count.back()) {
            throw std::runtime_error("Error: Not enough background SNP for "
                                     "calculation of competitive P-value!");
        }
    }
    size_t num_post_clump_snp(size_t i_region) const
    {
        return m_region_post_clump_count.at(i_region);
    }
    bool duplicated_size(size_t i_region) const
    {
        return m_region_size_duplicated.at(
            m_region_post_clump_count.at(i_region));
    }
    bool check_exclusion(const std::string &chr, const size_t loc);
private:
    struct region_bound
    {
        int chr;
        size_t start;
        size_t end;
    };
    std::string m_out_prefix; // for log file
    // for checking duplicated region
    // use member variable because both bed and msigdb needs this
    // and don't want to pass this around
    std::unordered_set<std::string> m_duplicated_names;
    // the name of the regions
    std::vector<std::string> m_region_name;
    // features that we'd like to capture
    std::vector<std::string> m_gtf_feature;
    // the actual region boundary
    // can't use vec2d because we don't know the size in advance
    std::vector<std::vector<region_bound>> m_region_list;
    // This is to indicate the current location of each region
    // This work because we assume all SNPs are sorted by their coordinates
    // in the same way as the region files.
    std::vector<size_t> m_snp_check_index;
    // the number of SNPs from the base+target that falls into the region
    std::vector<int> m_region_snp_count;
    std::vector<int> m_region_post_clump_count;
    // this is use for informing us if we would bother to store the permutation
    // results
    std::unordered_map<int, bool> m_region_size_duplicated;
    std::unordered_map<int, std::vector<int>> m_chr_index;
    int m_5prime = 0;
    int m_3prime = 0;

    bool in_feature(std::string in) const
    {
        // number of feature should be small enough such that
        // iterating the vector should be alright?
        for (auto& feature : m_gtf_feature) {
            if (in.compare(feature) == 0) return true;
        }
        return false;
    }


    void process_bed(const std::vector<std::string>& bed, Reporter& reporter);

    std::unordered_map<std::string, region_bound> process_gtf(
        const std::string& gtf,
        std::unordered_map<std::string, std::set<std::string>>& id_to_name,
        const std::string& out_prefix, Reporter& reporter);
    std::vector<Region::region_bound>
    solve_overlap(std::vector<Region::region_bound>& current_region);
    void process_msigdb(
        const std::string& msigdb,
        const std::unordered_map<std::string, region_bound>& gtf_info,
        const std::unordered_map<std::string, std::set<std::string>>&
            id_to_name,
        Reporter& reporter);
    void generate_background(
        const std::unordered_map<std::string, region_bound>& gtf_info,
        const size_t num_bed_region, Reporter& reporter);
    void read_background(
        const std::string& background,
        const std::unordered_map<std::string, region_bound>& gtf_info,
        const std::unordered_map<std::string, std::set<std::string>>&
            id_to_name,
        Reporter& reporter);
};

#endif /* PRSICE_INC_REGION_HPP_ */
