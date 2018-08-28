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

class Genotype;
class Region
{
public:
    /*!
     * \brief Contructor for the exclusion region. Use for excluding SNPs
     * \param exclusion_range is the user input string containing either a file
     *        or a string indicating the region to be excluded
     * \param reporter is the logger
     */
    Region(const std::string& exclusion_range, Reporter& reporter);
    /*!
     * \brief The constructor for the region object. Use for defining gene sets
     * \param feature is the vector string containing features to be included
     *        from the gtf file
     * \param window_5 is the length to be extended to the 5' end
     * \param window_3 is the length to be extended to the 3' end
     * \param genome_wide_background indicate if we want the whole genome to act
     * as background
     */
    Region(std::vector<std::string> feature, const int window_5,
           const int window_3, const bool run_perm,
           const bool genome_wide_background);
    virtual ~Region();
    /*!
     * \brief To generate the region boundaries
     * \param gtf is the gtf file input
     * \param msigdb is the GMT file from msigdb
     * \param bed is a list of bed file
     * \param snp_set is a file containing a single snp set
     * \param multi_snp_sets is a filt contain multiple gene set similar to GMT
     * \param out is the output prefix
     * \param background is the file containing the background information
     * \param genome_wide_background is a boolean indicate if we want to use
     *        whole genome as background
     * \param target is the target genotype
     * \param reporter is the logger
     */
    void generate_regions(const std::string& gtf, const std::string& msigdb,
                          const std::vector<std::string>& bed,
                          const std::string& snp_set,
                          const std::string& multi_snp_sets,
                          const std::string& background, const Genotype& target,
                          Reporter& reporter);
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
     * \brief This function will take the coordinate of a SNP and check if it
     * falls within the exclusion region
     * \param chr is the chromosome string
     * \param loc is the coordinate
     * \return true if the snp falls within an exclusion region
     */
    bool check_exclusion(const std::string& chr, const intptr_t loc);
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

private:
    // IMPORTANT: The end is non-inclusive
    struct region_bound
    {
        intptr_t chr;
        intptr_t start;
        intptr_t end;
    };
    std::string m_out_prefix; // for log file
    // for checking duplicated region
    // use member variable because both bed and msigdb needs this
    // and don't want to pass this around
    std::unordered_set<std::string> m_duplicated_names;
    std::vector<std::unordered_set<std::string>> m_snp_sets;
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

    intptr_t m_5prime = 0;
    intptr_t m_3prime = 0;
    bool m_has_background = false;
    bool m_run_perm = false;
    bool m_genome_wide_background = false;
    /*!
     * \brief A simple function to check if the string equals to one of the
     * feature
     * \param in is the input string
     * \return true if the input string equals to one of the feature
     */
    bool in_feature(std::string in) const
    {
        return std::find(m_gtf_feature.begin(), m_gtf_feature.end(), in)
               != m_gtf_feature.end();
    }
    /*!
     * \brief Function to process SNP set input (both single and multiple),
     * might want to separate that out for clarity
     * \param single_snp_set file contain a single snp set, contain one column
     * \param multi_snp_set file containing multiple snp sets. Similar to GMT
     * \param target is the target genotype. Use to get SNP coordinate
     * \param reporter is the logger
     */
    void process_snp_sets(const std::string& single_snp_set,
                          const std::string& multi_snp_set,
                          const Genotype& target, Reporter& reporter);
    /*!
     * \brief This function will read in a list of bed file names and construct
     * the region for each bed file
     * \param bed is the vector containing all bed file name
     * \param reporter is the logger
     */
    void process_bed(const std::vector<std::string>& bed, Reporter& reporter);
    /*!
     * \brief This function will process the gtf file and return a multimap
     *        containing all regions related to the gene id
     * \param gtf is the file name
     * \param id_to_name is the gene id to gene name translator
     * \param max_chr is the maximum possible chromosome code. Allow us to skip
     * some chromosomes
     * \param reporter is the logger
     * \return the multimap containing all regions of a gene
     */
    std::unordered_multimap<std::string, region_bound> process_gtf(
        const std::string& gtf,
        std::unordered_map<std::string, std::set<std::string>>& id_to_name,
        const uint32_t max_chr, Reporter& reporter);
    /*!
     * \brief Remove overlapped region within a list of region
     * \param current_region is the region input
     * \return a region vector contain non-overlapping regions
     */
    std::vector<Region::region_bound>
    solve_overlap(std::vector<Region::region_bound>& current_region);
    /*!
     * \brief Given the gtf information, this funciton will read in the GMT file
     *        from msigdb and transform all gene sets into regions
     * \param msigdb a string containing comma separated list of msigdb file
     * \param gtf_info is the gtf multimap generated by process_gtf
     * \param id_to_name is the gene id to gene name translator
     * \param reporter is the reporter
     */
    void process_msigdb(
        const std::string& msigdb,
        const std::unordered_multimap<std::string, region_bound>& gtf_info,
        const std::unordered_map<std::string, std::set<std::string>>&
            id_to_name,
        Reporter& reporter);
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
