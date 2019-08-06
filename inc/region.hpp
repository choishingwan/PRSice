// This file is part of PRSice-2, copyright (C) 2016-2019
// Shing Wan Choi, Paul F. O’Reilly
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

#include "IITree.h"
#include "cgranges.h"
#include "genotype.hpp"
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

class Region
{
public:
    /*!
     * \brief Default constructor that do nothing
     */
    Region() {}
    virtual ~Region();
    static void generate_exclusion(std::vector<IITree<size_t, size_t>>& cr,
                                   const std::string& exclusion_range);
    static size_t generate_regions(
        std::vector<IITree<size_t, size_t>>& gene_sets,
        std::vector<std::string>& region_names,
        std::unordered_map<std::string, std::vector<size_t>>& snp_in_sets,
        const std::vector<std::string>& feature, const size_t window_5,
        const size_t window_3, const bool genome_wide_background,
        const std::string& gtf, const std::string& msigdb,
        const std::vector<std::string>& bed, const std::string& snp_set,
        const std::string& background, const size_t max_chr,
        Reporter& reporter);

protected:
    static void load_background(
        const std::string& background, const size_t window_5,
        const size_t window_3, const size_t max_chr,
        std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        std::unordered_map<std::string, std::vector<size_t>>& snp_in_sets,
        bool printed_warning, std::vector<IITree<size_t, size_t>>& gene_sets,
        Reporter& reporter);
    static void load_msigdb(
        const std::string& msig,
        std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        std::vector<std::string>& region_names,
        std::unordered_set<std::string>& duplicated_sets, size_t& set_idx,
        Reporter& reporter);
    static void load_gtf(
        const std::string& gtf,
        const std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        const std::vector<std::string>& features, const size_t max_chr,
        const size_t window_5, const size_t window_3,
        std::vector<IITree<size_t, size_t>>& gene_sets,
        const bool genome_wide_background, const bool provided_background,
        Reporter& reporter);
    static bool load_bed_regions(
        const std::string& bed_file,
        std::vector<IITree<size_t, size_t>>& gene_sets, const size_t window_5,
        const size_t window_3, bool& printed_warning, const size_t set_idx,
        const size_t max_chr, std::vector<std::string>& region_names,
        std::unordered_set<std::string>& duplicated_sets, Reporter& reporter);
    static void load_snp_sets(
        std::string snp_file,
        std::unordered_map<std::string, std::vector<size_t>>& snp_in_sets,
        std::vector<std::string>& region_names,
        std::unordered_set<std::string>& duplicated_sets, size_t& set_idx,
        Reporter& reporter);
    static void is_bed_line(const std::vector<std::string>& bed_line,
                            size_t& column_size, bool& is_header)
    {
        if (bed_line.front() == "track" || bed_line.front() == "browser")
        {
            is_header = true;
            return;
        }

        if (column_size == 0) { column_size = bed_line.size(); }
        else if (column_size != bed_line.size())
        {
            std::string message =
                "Error: Inconsistent number of column in BED file. "
                "Current line has "
                + std::to_string(bed_line.size())
                + " "
                  "column and previous line has "
                + std::to_string(column_size) + "\n";
            throw std::runtime_error(message);
        }
        // don't bother to check the coordinate as those are kinda check later
        // on and allow for better error report (though we can also pull in
        // reporter here)
        std::string chr = bed_line.front();
        std::transform(chr.begin(), chr.end(), chr.begin(), ::toupper);
        if (chr.rfind("CHR") != 0)
        {
            if (!misc::isNumeric(bed_line.front()))
            {
                std::string message =
                    "Error: Invalid BED format. First field "
                    "of BED file should be chromosomal information\n";
                throw std::runtime_error(message);
            }
        }
        if (bed_line.size() > 5)
        {
            if (bed_line[5] != "." && bed_line[5] != "+" && bed_line[5] != "-")
            {
                std::string message = "Error: Undefined strand information: "
                                      + bed_line[5] + "\n";
                message.append(
                    "Valid strand characters are '.', '+' and '-'\n");
                throw std::runtime_error(message);
            }
        }
    }

    static bool in_feature(const std::string& in,
                           const std::vector<std::string>& feature)
    {
        return std::find(feature.begin(), feature.end(), in) != feature.end();
    }

    static bool parse_attribute(const std::string& attribute_str,
                                std::string& gene_id, std::string& gene_name)
    {
        assert(!attribute_str.empty());
        gene_id = "";
        gene_name = "";
        std::vector<std::string> attributes = misc::split(attribute_str, ";");
        std::vector<std::string> token;
        bool found_id = false, found_name = false;
        for (auto&& a : attributes)
        {
            // remove space before the attribute name
            misc::trim(a);
            if (a.rfind("gene_id", 0) == 0)
            {
                token = misc::split(a, " ");
                if (token.size() != 2)
                {
                    throw std::runtime_error(
                        "Error: Malformed attribute value: " + a);
                }
                gene_id = token.back();
                gene_id.erase(std::remove(gene_id.begin(), gene_id.end(), '\"'),
                              gene_id.end());
                if (found_name) return true;
                found_id = true;
            }
            else if (a.rfind("gene_name", 0) == 0)
            {
                token = misc::split(a, " ");
                if (token.size() != 2)
                {
                    throw std::runtime_error(
                        "Error: Malformed attribute value: " + a);
                }
                gene_name = token.back();
                gene_name.erase(
                    std::remove(gene_name.begin(), gene_name.end(), '\"'),
                    gene_name.end());
                if (found_id) return true;
                found_name = true;
            }
        }
        return false;
    }
};

#endif /* PRSICE_INC_REGION_HPP_ */
