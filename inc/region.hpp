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
    size_t generate_regions(
        const std::unordered_map<std::string, size_t>& included_snp_idx,
        const std::vector<std::unique_ptr<SNP>>& included_snps,
        const size_t max_chr);

    const std::vector<std::string>& get_names() const { return m_region_name; }

    const std::vector<IITree<size_t, size_t>>& get_gene_sets() const
    {
        return m_gene_sets;
    }

protected:
    void load_background(
        const std::unordered_map<std::string, size_t>& snp_list_idx,
        const std::vector<std::unique_ptr<SNP>>& snp_list, const size_t max_chr,
        std::unordered_map<std::string, std::vector<size_t>>& msigdb_list);

    static void extend_region(std::string_view strand, const size_t wind_5,
                              const size_t wind_3, size_t& start, size_t& end)
    {
        if (strand == "-")
        {
            if (wind_3 > start)
                start = 1;
            else
                start -= wind_3;
            end += wind_5;
        }
        else if (strand == "+" || strand == ".")
        {
            if (wind_5 > start)
                start = 1;
            else
                start -= wind_5;
            end += wind_3;
        }
        else
        {
            throw std::runtime_error("Error: Undefined strand "
                                     "information. Possibly a malform "
                                     "file: "
                                     + std::string(strand));
        }
    }
    static std::tuple<size_t, size_t>
    start_end(const std::string_view& start_str,
              const std::string_view& end_str, const bool zero_based)
    {
        return start_end(std::string(start_str), std::string(end_str),
                         zero_based);
    }
    static std::tuple<size_t, size_t> start_end(const std::string& start_str,
                                                const std::string& end_str,
                                                const bool zero_based)
    {
        size_t start, end;
        try
        {
            // range format is one based, we need to add one if input is zero
            // based
            start = misc::convert<size_t>(start_str) + zero_based;
        }
        catch (...)
        {
            throw std::runtime_error(
                "Error: Invalid start coordinate: " + start_str + "\n");
        }
        try
        {
            end = misc::convert<size_t>(end_str);
        }
        catch (...)
        {
            throw std::runtime_error(
                "Error: Invalid end coordinate: " + start_str + "\n");
        }
        if (start > end)
        {
            throw std::runtime_error(
                "Error: Start coordinate should be smaller "
                "than end coordinate!\nstart:"
                + start_str + "\nend: " + end_str + "\n");
        }
        return {start, end};
    }
    bool add_gene_region_from_gtf(
        const std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        const std::string id, const size_t chr, const size_t start,
        const size_t end)
    {
        if (id.empty()) return false;
        auto&& id_search = msigdb_list.find(id);
        if (id_search != msigdb_list.end())
        {
            for (auto&& idx : id_search->second)
            { m_gene_sets[chr].add(start, end, idx); }
            return true;
        }
        return false;
    }
    void load_msigdb(
        std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        std::unique_ptr<std::istream> input, size_t& set_idx);
    std::tuple<size_t, size_t, size_t> transverse_gtf(
        const std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        const std::streampos file_length, const size_t max_chr,
        const bool gz_input, std::unique_ptr<std::istream> gtf_stream);
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
    void transverse_snp_file(
        const std::unordered_map<std::string, size_t>& snp_list_idx,
        const std::vector<std::unique_ptr<SNP>>& snp_list,
        const bool is_set_file, std::unique_ptr<std::istream> input,
        size_t& set_idx);
    void
    load_snp_sets(const std::unordered_map<std::string, size_t>& snp_list_idx,
                  const std::vector<std::unique_ptr<SNP>>& snp_list,
                  const std::string& snp_file, size_t& set_idx);
    std::tuple<std::string, std::string, bool>
    get_set_name(const std::string& input)
    {
        std::vector<std::string> token = misc::split(input, ":");
        if (token.size() > 2)
        {
            throw std::runtime_error("Error: Undefine file input format: "
                                     + input);
        }
        // when used didn't provide a set name, then it will automatically use
        // the file name as the set name (as token size = 1, back == front)
        return {token.front(), token.back(), token.size() == 2};
    }
    static void read_bed(std::unique_ptr<std::istream> bed,
                         std::vector<IITree<size_t, size_t>>& cr,
                         bool& print_bed_strand_warning,
                         const size_t wind_5 = 0, const size_t wind_3 = 0,
                         const size_t max_chr = MAX_POSSIBLE_CHROM,
                         const size_t set_idx = 0,
                         const bool ZERO_BASED = true);

    static bool valid_strand(std::string_view input)
    {
        return !(input != "." && input != "+" && input != "-");
    }
    static bool is_bed_header(const std::vector<std::string_view>& bed_line,
                              size_t& column_size)
    {
        if (bed_line.front() == "track" || bed_line.front() == "browser"
            || bed_line.front().at(0) == '#')
        { return true; }
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
        if (bed_line.size() > 5 && !valid_strand(bed_line[+BED::STRAND]))
        {
            throw std::runtime_error(
                "Error: Undefined strand information: "
                + std::string(bed_line[+BED::STRAND])
                + "\nValid strand characters are '.', '+' and '-'\n");
        }
        return false;
    }

    bool in_feature(const std::string_view& in,
                    const std::vector<std::string>& feature)
    {
        return std::any_of(feature.begin(), feature.end(),
                           [&](const std::string& elem) { return elem == in; });
    }

    std::tuple<std::string, bool> parse_gene_id(std::string_view substr,
                                                std::string name)
    {
        std::string gene_id;
        if (substr.rfind(name, 0) == 0 || substr.rfind(" " + name, 0) == 0)
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
            return {gene_id, true};
        }
        return {"", false};
    }
    bool find_gene_info(std::string_view substr, std::string& gene_id,
                        std::string& gene_name, bool& found_id,
                        bool& found_name)
    {
        if (!found_id)
        {
            std::tie(gene_id, found_id) = parse_gene_id(substr, "gene_id");
            if (found_name && found_id) return true;
        }
        if (!found_name)
        {
            std::tie(gene_name, found_name) =
                parse_gene_id(substr, "gene_name");
            if (found_name && found_id) return true;
        }
        return false;
    }

    void parse_attribute(std::string_view attribute_str, std::string& gene_id,
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
            { return; }
        }
    }

    std::vector<IITree<size_t, size_t>> m_gene_sets;
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
    Reporter* m_reporter;
};

#endif /* PRSICE_INC_REGION_HPP_ */
