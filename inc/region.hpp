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
#include "snp.hpp"
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
    static void generate_exclusion(cgranges_t* cr,
                                   const std::string& exclusion_range);
    static size_t add_flags(std::vector<std::string>& region_names,
                            const std::vector<std::string>& feature,
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

    static bool in_feature(const std::string& in,
                           const std::vector<std::string>& feature)
    {
        return std::find(feature.begin(), feature.end(), in) != feature.end();
    }
};

#endif /* PRSICE_INC_REGION_HPP_ */
