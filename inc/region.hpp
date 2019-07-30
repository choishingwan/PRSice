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
    static void generate_exclusion(std::vector<IITree<int, int>>& cr,
                                   const std::string& exclusion_range);
    static size_t generate_regions(
        std::vector<IITree<int, int>>& gene_sets,
        std::vector<std::string>& region_names,
        std::unordered_map<std::string, std::vector<int>>& snp_in_sets,
        const std::vector<std::string>& feature, const int window_5,
        const int window_3, const bool genome_wide_background,
        const std::string& gtf, const std::string& msigdb,
        const std::vector<std::string>& bed, const std::string& snp_set,
        const std::string& background, const uint32_t max_chr,
        Reporter& reporter);

protected:
    static void load_background(
        const std::string& background, const int window_5, const int window_3,
        const uint32_t max_chr,
        std::unordered_map<std::string, std::vector<int>>& msigdb_list,
        bool printed_warning, std::vector<IITree<int, int>>& gene_sets,
        Reporter& reporter);
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
        const int window_5, const int window_3,
        std::vector<IITree<int, int>>& gene_sets,
        const bool genome_wide_background, const bool provided_background,
        Reporter& reporter);
    static bool load_bed_regions(
        const std::string& bed_file, std::vector<IITree<int, int>>& gene_sets,
        const int window_5, const int window_3, bool& print_warning,
        const int set_idx, const uint32_t max_chr,
        std::vector<std::string>& region_names,
        std::unordered_set<std::string> duplicated_sets, Reporter& reporter);
    static void load_snp_sets(
        std::string snp_file,
        std::unordered_map<std::string, std::vector<int>>& snp_in_sets,
        std::vector<std::string>& region_names,
        std::unordered_set<std::string>& duplicated_sets, int& set_idx,
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
        if (bed_line.front().rfind("CHR") != 0)
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


    static std::vector<std::string> get_attribute(const std::string& attribute)
    {
        std::size_t prev = 0, pos;
        std::string temp;
        std::vector<std::string> result(2, "");
        bool add_next = false;
        size_t add_id = 3;
        size_t num_added = 0;
        while ((pos = attribute.find_first_of(" ;", prev)) != std::string::npos
               && num_added < 2)
        {
            if (pos > prev)
            {
                temp = attribute.substr(prev, pos - prev);
                if (temp.size() == 7 && temp.substr(5) == "id")
                {
                    add_next = true;
                    add_id = 0;
                }
                else if (temp.size() == 9 && temp.substr(5) == "name")
                {
                    add_next = true;
                    add_id = 1;
                }
                else if (add_next)
                {
                    add_next = false;
                    temp.erase(std::remove(temp.begin(), temp.end(), '\"'),
                               temp.end());
                    result[add_id] = temp;
                    num_added++;
                }
                assert(add_id < 2);
            }
            prev = pos + 1;
        }
        if (add_next)
            result[add_id] = (attribute.substr(prev, std::string::npos));
        return result;
    }
};

#endif /* PRSICE_INC_REGION_HPP_ */
