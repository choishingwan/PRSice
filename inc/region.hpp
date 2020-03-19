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
#include "commander.hpp"
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
#include <string_view>
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
    Region(const GeneSets& set, Reporter* reporter)
        : m_bed(set.bed)
        , m_feature(set.feature)
        , m_msigdb(set.msigdb)
        , m_snp_set(set.snp)
        , m_background(set.background)
        , m_gtf(set.gtf)
        , m_window_5(set.wind_5)
        , m_window_3(set.wind_3)
        , m_genome_wide_background(set.full_as_background)
        , m_reporter(reporter)
    {
    }
    virtual ~Region();
    static void generate_exclusion(std::vector<IITree<size_t, size_t>>& cr,
                                   const std::string& exclusion_range);
    size_t generate_regions(const size_t max_chr);
    std::vector<std::string> get_names() const { return m_region_name; }

    const std::vector<IITree<size_t, size_t>>& get_gene_sets() const
    {
        return m_gene_sets;
    }
    const std::unordered_map<std::string, std::vector<size_t>>&
    get_snp_sets() const
    {
        return m_snp_in_sets;
    }

protected:
    void load_background(
        const size_t max_chr,
        std::unordered_map<std::string, std::vector<size_t>>& msigdb_list);

    void extend_region(size_t& start, size_t& end, std::string_view strand)
    {
        if (strand == "-")
        {
            if (m_window_3 > start)
                start = 1;
            else
                start -= m_window_3;
            end += m_window_5;
        }
        else if (strand == "+" || strand == ".")
        {
            if (m_window_5 > start)
                start = 1;
            else
                start -= m_window_5;
            end += m_window_3;
        }
        else
        {
            throw std::runtime_error("Error: Undefined strand "
                                     "information. Possibly a malform "
                                     "file: "
                                     + std::string(strand));
        }
    }
    static void start_end(std::string_view start_str, std::string_view end_str,
                          const size_t pad, size_t& start, size_t& end)
    {
        start_end(std::string(start_str), std::string(end_str), pad, start,
                  end);
    }

    static void start_end(const std::string& start_str,
                          const std::string& end_str, const size_t pad,
                          size_t& start, size_t& end)
    {
        try
        {
            // That's because bed is 0 based and range format is 1
            // based. and type for bed is 1 and type for range is 0
            start = misc::Convertor::convert<size_t>(start_str) + pad;
        }
        catch (...)
        {
            throw std::runtime_error("start");
        }
        try
        {
            end = misc::Convertor::convert<size_t>(end_str);
        }
        catch (...)
        {
            throw std::runtime_error("end ");
        }
        if (start > end)
        {
            // don't check if it's already error out
            std::string message =
                "Error: Start coordinate should be smaller than "
                "end coordinate!\n";
            message.append("start: " + std::to_string(start) + "\n");
            message.append("end: " + std::to_string(end) + "\n");
            throw std::logic_error(message);
        }
    }
    void load_msigdb(
        const std::string& msig,
        std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        size_t& set_idx);
    void load_gtf(
        const std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        const size_t max_chr);
    bool duplicated_set(std::string_view set_name)
    {
        return duplicated_set(std::string(set_name));
    }
    bool duplicated_set(const std::string& set_name)
    {
        if (m_processed_sets.find(set_name) != m_processed_sets.end())
        {
            m_reporter->report("Warning: Set name of " + set_name
                               + " is duplicated, it will be ignored");
            return true;
        }
        m_processed_sets.insert(set_name);
        m_region_name.push_back(set_name);
        return false;
    }
    bool load_bed_regions(const std::string& bed_file, const size_t set_idx,
                          const size_t max_chr);
    void load_snp_sets(const std::string& snp_file, size_t& set_idx);
    bool get_set_name(const std::string& input, std::string& file_name,
                      std::string& set_name)
    {
        std::vector<std::string> token = misc::split(input, ":");
        if (token.size() > 2)
        {
            throw std::runtime_error("Error: Undefine file input format: "
                                     + input);
        }
        file_name = token.front();
        set_name = token.back();
        return token.size() == 2;
    }
    static void read_bed(const std::string& bed);
    static void read_bed(std::string_view bed,
                         std::vector<IITree<size_t, size_t>>& m_gene_sets);
    static bool valid_chr(const std::string& input)
    {
        std::string chr = input;
        misc::to_upper(chr);
        if (chr.rfind("CHR") != 0)
        {
            if (!misc::isNumeric(chr)) { return false; }
        }
        return true;
    }

    static bool valid_strand(std::string_view input)
    {
        return !(input != "." && input != "+" && input != "-");
    }
    static void is_bed_line(const std::vector<std::string_view>& bed_line,
                            size_t& column_size, bool& is_header)
    {
        if (bed_line.front() == "track" || bed_line.front() == "browser")
        {
            is_header = true;
            return;
        }
        if (bed_line.size() < 3)
        {
            throw std::runtime_error("Error: Malformed BED file. BED file "
                                     "should contain at least 3 column "
                                     "for all rows!\n");
        }
        if (column_size == 0) { column_size = bed_line.size(); }
        else if (column_size != bed_line.size())
        {
            throw std::runtime_error(
                "Error: Inconsistent number of column in file. "
                "Current line has "
                + std::to_string(bed_line.size())
                + " "
                  "column and previous line has "
                + std::to_string(column_size) + "\n");
        }
        if (!valid_chr(std::string(bed_line.front())))
        {
            throw std::runtime_error(
                "Error: Invalid file format. First field "
                "of file should be chromosomal information\n");
        }
        if (bed_line.size() > 5)
        {
            if (!valid_strand(bed_line[+BED::STRAND]))
            {
                throw std::runtime_error(
                    "Error: Undefined strand information: "
                    + std::string(bed_line[+BED::STRAND])
                    + "\nValid strand characters are '.', '+' and '-'\n");
            }
        }
    }

    bool in_feature(const std::string_view& in,
                    const std::vector<std::string>& feature)
    {
        return std::any_of(feature.begin(), feature.end(),
                           [&](const std::string& elem) { return elem == in; });
    }
    bool in_feature(const std::string& in,
                    const std::vector<std::string>& feature)
    {
        return std::find(feature.begin(), feature.end(), in) != feature.end();
    }

    bool find_gene_info(std::string_view substr, std::string& gene_id,
                        std::string& gene_name, bool& found_id,
                        bool& found_name)
    {
        if (substr.rfind("gene_id", 0) == 0)
        {
            auto token = misc::tokenize(substr, " ");
            if (token.size() != 2)
            {
                throw std::runtime_error("Error: Malformed attribute value: "
                                         + std::string(substr));
            }
            gene_id = std::string(token.back());
            gene_id.erase(std::remove(gene_id.begin(), gene_id.end(), '\"'),
                          gene_id.end());
            if (found_name) return true;
            found_id = true;
        }
        else if (substr.rfind("gene_name", 0) == 0)
        {
            auto token = misc::tokenize(substr, " ");
            if (token.size() != 2)
            {
                throw std::runtime_error("Error: Malformed attribute value: "
                                         + std::string(substr));
            }
            gene_name = std::string(token.back());
            gene_name.erase(
                std::remove(gene_name.begin(), gene_name.end(), '\"'),
                gene_name.end());
            if (found_id) return true;
            found_name = true;
        }
        return false;
    }

    bool parse_attribute(std::string_view attribute_str, std::string& gene_id,
                         std::string& gene_name)
    {
        assert(!attribute_str.empty());
        gene_id = "";
        gene_name = "";
        std::vector<std::string_view> token;
        bool found_id = false, found_name = false;
        token = misc::tokenize(attribute_str, ";");
        for (auto&& item : token)
        {
            // only return true when both gene_name and gene_id are found
            if (find_gene_info(item, gene_id, gene_name, found_id, found_name))
            { return true; }
        }
        return false;
    }

    std::vector<IITree<size_t, size_t>> m_gene_sets;
    std::unordered_map<std::string, std::vector<size_t>> m_snp_in_sets;
    std::vector<std::string> m_bed;
    std::vector<std::string> m_feature;
    std::vector<std::string> m_msigdb;
    std::vector<std::string> m_snp_set;
    std::vector<std::string> m_region_name;
    std::unordered_set<std::string> m_processed_sets;
    std::string m_background;
    std::string m_gtf;
    size_t m_window_5 = 0;
    size_t m_window_3 = 0;
    bool m_genome_wide_background;
    bool m_printed_bed_strand_warning = false;
    Reporter* m_reporter;
};

#endif /* PRSICE_INC_REGION_HPP_ */
