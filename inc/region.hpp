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

#include "cgranges.h"
#include "gzstream.h"
#include "misc.hpp"
#include "plink_common.hpp"
#include "reporter.hpp"
#include "storage.hpp"
#include <fstream>
#include <iostream>
#include <iterator>
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

class Genotype;
class Region
{
public:
    /*!
     * \brief Default constructor that do nothing
     */
    Region() {}
    virtual ~Region();

    /*!
     * \brief Function to construct the flag for SNP membership. This is rather
     * fragile as it require the SNP being sorted in the same order as the
     * region. This is achieved by using the same chromosome code.
     * \param chr the chromosome code for the SNP
     * \param rs the RS ID of the SNP, useful for debugging
     * \param loc the coordinate of the SNP
     * \param flag the binary flag of the SNP regional membership
     */
    void update_flag(const intptr_t chr, const std::string& rs, intptr_t loc,
                     std::vector<uintptr_t>& flag);
    /*!
     * \brief This will print a log of the number of region included
     * \param reporter the logger
     */
    void print_region_number(Reporter& reporter) const;

    /*!
     * \brief This function is use to clean out the m_region list and
     * m_snp_check_index to release some memory as they will not be used after
     * the construction of flags
     */
    void clean()
    {
        m_region_list = std::vector<std::vector<region_bound>>();
        m_snp_check_index = std::vector<size_t>();
    }
    /*!
     * \brief Store information fo the number of SNP in each region
     * \param count is the input containing the count for each region
     */
    void post_clump_count(std::vector<int>& count)
    {
        int max = 0;
        // we resize our region count storage
        m_region_post_clump_count.resize(count.size());
        // and calculate the last region
        const size_t last_region_index = count.size() - 1;
        for (size_t i = 0; i < count.size(); ++i) {
            if (i != last_region_index && i != 0) {
                // we want to know what is the maximum number of SNP for all the
                // set involves (except base and background)
                max = std::max(count[i], max);
            }
            // then store in the count
            m_region_post_clump_count[i] = count[i];
        }
        // we won't do background with the base group
        if (count.size() > 1 && max > count.back()
            && m_region_name.back() == "Background")
        {
            // if we have more than one group and the maximum SNP count is
            // larger than the number of SNPs in background (can happen if the
            // backgroudn isn't generated from GTF &/ region input is bed
            throw std::runtime_error("Error: Not enough background SNP for "
                                     "calculation of competitive P-value!");
        }
    }
    /*!
     * \brief Get the number of SNP contain in the region
     * \param i_region is the region index
     * \return the number of SNP
     */
    int num_post_clump_snp(size_t i_region) const
    {
        return m_region_post_clump_count.at(i_region);
    }
    /*!
     * \brief Return the number of regions involved
     * \return the number of regions involved
     */
    size_t size() const { return m_region_name.size(); }
    /*!
     * \brief Return the number of boundaries for the i th set. Use for unit
     * testing
     * \param i is the index of the set
     * \return The number of boundaries for the i th set
     */
    size_t num_bound(size_t i) const { return m_region_list.at(i).size(); }
    /*!
     * \brief Return the name of the i th set
     * \param i is the set index
     * \return  return the name of the i th set
     */

    std::string get_name(size_t i) const { return m_region_name.at(i); }
    /*!
     * \brief return the names of all sets
     * \return a string vector containing the name of all sets
     */
    std::vector<std::string> names() const { return m_region_name; }

    static void generate_exclusion(cgranges_t* cr,
                                   const std::string& exclusion_range);
    static void add_flags(const std::vector<std::string>& feature,
                          const int window_5, const int window_3,
                          const bool genome_wide_background,
                          const std::string& gtf, const std::string& msigdb,
                          const std::vector<std::string>& bed,
                          const std::string& snp_set,
                          const std::string& background, Genotype& target,
                          Reporter& reporter);

private:
    static void load_background(
        const std::string& background, const int window_5, const int window_3,
        const uint32_t max_chr,
        std::unordered_map<std::string, std::vector<int>>& msigdb_list,
        bool printed_warning, cgranges_t* gene_sets, Reporter& reporter);
    static void
    load_msigdb(const std::string& msig,
                std::unordered_map<std::string, std::vector<int>>& msigdb_list,
                std::vector<std::string>& region_names,
                std::unordered_set<std::string> duplicated_sets, int& set_idx,
                Reporter& reporter);
    static void load_gtf(
        const std::string& gtf,
        const std::unordered_map<std::string, std::vector<int>>& msigdb_list,
        const std::vector<std::string>& features, const uint32_t max_chr,
        const int window_5, const int window_3, cgranges_t* gene_sets,
        const bool genome_wide_background, Reporter& reporter);
    static bool load_bed_regions(
        const std::string& bed_file, cgranges_t* gene_sets, const int window_5,
        const int window_3, bool& print_warning, const int set_idx,
        const uint32_t max_chr, std::vector<std::string>& region_names,
        std::unordered_set<std::string> duplicated_sets, Reporter& reporter);
    static void load_snp_sets(
        std::string snp_file,
        std::unordered_map<std::string, std::vector<int>>& snp_in_sets,
        std::vector<std::string>& region_names,
        std::unordered_set<std::string>& duplicated_sets, int& set_idx,
        Reporter& reporter);
    static size_t num_snp_set(const std::string& name)
    {
        std::ifstream input;
        std::vector<std::string> file_name = misc::split(name, ":");
        input.open(file_name[0].c_str());
        if (!input.is_open()) {
            std::string error =
                "Error: Cannot open file: " + file_name[0] + "\n";
            throw std::runtime_error(error);
        }
        std::string line;
        std::vector<std::string> token;
        size_t num_set = 0;
        while (std::getline(input, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            if (token.size() == 1) return 1;
            num_set++;
        }
        return num_set;
    }


    std::string m_out_prefix; // for log file
    // for checking duplicated region
    // use member variable because both bed and msigdb needs this
    // and don't want to pass this around


    static bool in_feature(const std::string& in,
                           const std::vector<std::string>& feature)
    {
        return std::find(feature.begin(), feature.end(), in) != feature.end();
    }


    /*!
     * \brief This function will take in the gtf information and generate the
     *        background region containing all the feature
     * \param gtf_info is the boundaries generated from the gtf file
     */
    void generate_background(
        const std::unordered_multimap<std::string, region_bound>& gtf_info);
    /*!
     * \brief This function will read in a background file and generate the
     *        background region based on the input
     * \param background is the string of format name:type
     * \param gtf_info is the gtf multimap generated by process_gtf
     * \param id_to_name is the gene id to gene name translator
     * \param reporter is the logger
     */
    void read_background(
        const std::string& background,
        const std::unordered_multimap<std::string, region_bound>& gtf_info,
        const std::unordered_map<std::string, std::set<std::string>>&
            id_to_name,
        Reporter& reporter);
};

#endif /* PRSICE_INC_REGION_HPP_ */
