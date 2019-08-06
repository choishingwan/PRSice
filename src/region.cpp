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


#include "region.hpp"
#include "genotype.hpp"

// This is for exclusion region
// end boundary is inclusive
// e.g. bound with chr1:1-10 will remove any SNPs on chr1 with coordinate
// from 1 to 10

void Region::generate_exclusion(std::vector<IITree<int, int>>& cr,
                                const std::string& exclusion_range)
{
    // do nothing when no exclusion is required.
    if (exclusion_range.empty()) return;
    std::vector<std::string> exclude_regions =
        misc::split(exclusion_range, ",");
    // file input is represented by the lack of :
    std::vector<std::string> range, boundary;
    for (auto&& region : exclude_regions)
    {
        range = misc::split(region, ":");
        if (range.size() == 1)
        {
            // bed file input
            std::ifstream input_file;
            input_file.open(region);
            if (!input_file.is_open())
            {
                std::string message = "Error: " + region
                                      + " cannot be open. Please check you "
                                        "have the correct input";
                throw std::runtime_error(message);
            }
            // BED starts from 0 and end is exclusive
            std::string line;
            size_t column_size = 0;
            bool is_header = false;
            while (std::getline(input_file, line))
            {
                misc::trim(line);
                is_header = false;
                boundary = misc::split(line);
                if (boundary.size() < 3)
                {
                    std::string message = "Error: Malformed BED file. BED file "
                                          "should contain at least 3 column "
                                          "for all rows!\n";
                    throw std::runtime_error(message);
                }
                try
                {
                    is_bed_line(boundary, column_size, is_header);
                }
                catch (const std::exception& ex)
                {
                    std::string message = ex.what();
                    message.append("Problematic line: " + line);
                    throw std::runtime_error(message);
                }
                if (is_header) continue; // skip header
                int chr = get_chrom_code_raw(boundary[0].c_str());
                int low_bound = misc::string_to_int(boundary[1].c_str()) + 1;
                // +1 because start at 0
                int upper_bound = misc::string_to_int(boundary[2].c_str());
                // Do nothing because while BED end is exclusive, it is 0 base.
                // For an inclusive end bound, we will need to do -1 (make it
                // inclusive) and +1 (transform to 1 base)
                if (low_bound > upper_bound || low_bound < 1 || upper_bound < 1)
                {
                    std::string message =
                        "Error: Invalid exclusion coordinate. "
                        "Coordinate must be larger than 1 and the end "
                        "coordinate "
                        "must be larger than or equal to the start "
                        "coordinate\n";
                    throw std::runtime_error(message);
                }
                if (chr >= 0)
                {
                    // ignore any chromosome that we failed to parse, which
                    // will have chr < 0
                    size_t chr_num = static_cast<size_t>(chr);
                    if (cr.size() < chr_num + 1) { cr.resize(chr_num + 1); }
                    cr[chr_num].add(low_bound, upper_bound, 0);
                }
            }
            input_file.close();
        }
        else if (range.size() == 2)
        {
            int chr = get_chrom_code_raw(range[0].c_str());
            boundary = misc::split(range[1], "-");
            // we take this as absolute position.
            // E.g. remove SNP with coordinate == 10 if
            // boundary = 10-20
            // boundary = 1-10
            // boundary = 10
            if (boundary.size() > 2)
            {
                std::string message =
                    "Error: Invalid exclusion range format. "
                    "Should be chr:start, chr:start-end or a bed file\n";
                throw std::runtime_error(message);
            }
            int low_bound = misc::string_to_int(boundary.front().c_str());
            int upper_bound = misc::string_to_int(boundary.back().c_str());
            // the library find overlap, which for SNP at 10
            // the boundary should be defined as 9-11 when we read in the SNP
            // here we do nothing but sainity check of the input
            if (low_bound > upper_bound || low_bound < 1 || upper_bound < 1)
            {
                std::string message =
                    "Error: Invalid exclusion coordinate. "
                    "Coordinate must be larger than 1 and the end coordinate "
                    "must be larger than or equal to the start coordinate\n";
                throw std::runtime_error(message);
            }
            // the range is completely inclusive []
            if (chr >= 0)
            {
                // ignore any chromosome that we failed to parse, which
                // will have chr < 0
                size_t chr_num = static_cast<size_t>(chr);
                if (cr.size() < chr_num + 1) { cr.resize(chr_num + 1); }
                cr[chr_num].add(low_bound, upper_bound, 0);
            }
        }
        else
        {
            std::string message =
                "Error: Invalid exclusion range format. "
                "Should be chr:start, chr:start-end or a bed file\n";
            throw std::runtime_error(message);
        }
    }
    // also index the tree
    for (auto&& tree : cr) { tree.index(); }
}

size_t Region::generate_regions(
    std::vector<IITree<int, int>>& gene_sets,
    std::vector<std::string>& region_names,
    std::unordered_map<std::string, std::vector<int>>& snp_in_sets,
    const std::vector<std::string>& feature, const int window_5,
    const int window_3, const bool genome_wide_background,
    const std::string& gtf, const std::string& msigdb,
    const std::vector<std::string>& bed, const std::string& snp_set,
    const std::string& background, const uint32_t max_chr, Reporter& reporter)
{
    // should be a fresh start each time
    gene_sets.clear();
    region_names.clear();
    snp_in_sets.clear();
    std::string message = "Start processing gene set information\n";
    message.append("==================================================");
    reporter.report(message);
    // we can now utilize the last field of cgranges as the index of gene
    // set of interest
    std::unordered_set<std::string> duplicated_sets;
    region_names.push_back("Base");
    region_names.push_back("Background");
    duplicated_sets.insert("Base");
    duplicated_sets.insert("Background");
    // 0 reserved for base
    // 1 reserved for background
    int set_idx = 2;
    bool printed_warning = false;
    for (auto&& bed_file : bed)
    {
        message = "Loading " + bed_file + " (BED)";
        reporter.report(message);
        set_idx += load_bed_regions(bed_file, gene_sets, window_5, window_3,
                                    printed_warning, set_idx, max_chr,
                                    region_names, duplicated_sets, reporter);
    }
    // SNP sets
    // we allow ambiguous snp set input
    // It all depends on the file content
    // (only one column = SNP list, mutliple column = SNP sets)
    std::vector<std::string> snp_sets = misc::split(snp_set, ",");
    for (auto&& s : snp_sets)
    {
        message = "Loading " + s + " (SNP sets)";
        reporter.report(message);
        load_snp_sets(s, snp_in_sets, region_names, duplicated_sets, set_idx,
                      reporter);
    }
    // the SNP sets information are stored within snp_in_sets
    if (!msigdb.empty() && gtf.empty())
    {
        std::string error_message =
            "Error: MSigDB input requires a complementary"
            "GTF file!\n";
        throw std::runtime_error(error_message);
    }
    std::unordered_map<std::string, std::vector<int>> msigdb_list;
    std::vector<std::string> msigdb_name = misc::split(msigdb, ",");
    for (auto&& msig : msigdb_name)
    {
        message = "Loading " + msig + " (MSigDB)";
        reporter.report(message);
        load_msigdb(msig, msigdb_list, region_names, duplicated_sets, set_idx,
                    reporter);
    }
    // now process the gtf file and add the regions
    if (!background.empty() && !genome_wide_background)
    {
        // we need to read in a background file, but ignore the
        // case where we use the gtf as background as we should've
        // already done that
        message = "Loading background set from " + background;
        reporter.report(message);
        load_background(background, window_5, window_3, max_chr, msigdb_list,
                        printed_warning, gene_sets, reporter);
    }
    if (!gtf.empty() && (!msigdb.empty() || !genome_wide_background))
    {
        // generate the region information based on the msigdb info
        // or generate for the background
        message = "Loading GTF file: " + gtf;
        reporter.report(message);
        load_gtf(gtf, msigdb_list, feature, max_chr, window_5, window_3,
                 gene_sets, genome_wide_background, !background.empty(),
                 reporter);
    }
    // index gene list
    for (auto&& tree : gene_sets) tree.index();
    // now we can go through each SNP and update their flag
    const size_t num_sets = region_names.size();
    if (num_sets == 2)
    {
        // because we will always have base and background
        message = "1 region included";
    }
    else
    {
        // -1 to remove the background count, as we are not going to print
        // the background anyway
        message = "A total of " + misc::to_string(num_sets - 2)
                  + " regions plus the base region are included";
    }
    reporter.report(message);
    return num_sets;
}

void Region::load_background(
    const std::string& background, const int window_5, const int window_3,
    const uint32_t max_chr,
    std::unordered_map<std::string, std::vector<int>>& msigdb_list,
    bool printed_warning, std::vector<IITree<int, int>>& gene_sets,
    Reporter& reporter)
{
    const std::unordered_map<std::string, int> file_type {
        {"bed", 1}, {"range", 0}, {"gene", 2}};
    // format of the background string should be name:format
    const std::vector<std::string> background_info =
        misc::split(background, ":");
    if (background_info.size() != 2)
    {
        std::string error =
            "Error: Format of --background should be <File Name>:<File Type>";
        throw std::runtime_error(error);
    }
    // check if we know the format
    auto&& type = file_type.find(background_info[1]);
    if (type == file_type.end())
    {
        std::string error = "Error: Undefined file type. Supported formats are "
                            "bed, gene or range";
        throw std::runtime_error(error);
    }
    // open the file
    std::ifstream input;
    input.open(background_info[0].c_str());
    if (!input.is_open())
    {
        std::string error_message = "Error: Cannot open background file: "
                                    + background_info[0] + " to read\n";
        throw std::runtime_error(error_message);
    }
    std::string line;
    bool error = false;
    if (type->second == 0 || type->second == 1)
    {
        // this is either a range format or a bed file
        // the only difference is range is 1 based and bed is 0 based
        size_t num_line = 0;
        std::vector<std::string> token;
        while (std::getline(input, line))
        {
            num_line++;
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            if (token.size() < 3)
            {
                std::string message =
                    "Error: " + background_info[0]
                    + " contain less than 3 columns, it will be ignored";
                throw std::runtime_error(message);
            }

            int32_t chr_code = get_chrom_code_raw(token[0].c_str());
            if (chr_code > static_cast<int32_t>(max_chr) || chr_code < 0)
                continue;

            if (token.size() <= +BED::STRAND && (window_5 > 0 || window_3 > 0)
                && (window_5 != window_3) && !printed_warning)
            {
                std::string message = "Warning: You bed file does not contain "
                                      "strand information, we will assume all "
                                      "regions are on the positive strand, "
                                      "e.g. start coordinates always on the 5' "
                                      "end";
                reporter.report(message);
                printed_warning = true;
            }
            int start = 0, end = 0;
            std::string message = "";
            try
            {
                start = misc::string_to_int(token[+BED::START].c_str());
                if (start >= 0)
                    // That's because bed is 0 based and range format is 1
                    // based. and type for bed is 1 and type for range is 0
                    start += type->second;
                else
                {
                    message.append("Error: Negative Start Coordinate at line "
                                   + misc::to_string(num_line) + "!");
                    reporter.report(message);
                    error = true;
                }
            }
            catch (...)
            {
                message.append("Error: Cannot convert start coordinate! (line: "
                               + misc::to_string(num_line) + ")!");
                reporter.report(message);
                error = true;
            }
            try
            {
                end = misc::string_to_int(token[+BED::END].c_str());
                if (end < 0)
                {
                    message.append("Error: Negative End Coordinate at line "
                                   + misc::to_string(num_line) + "!");
                    reporter.report(message);
                    error = true;
                }
            }
            catch (...)
            {
                message.append("Error: Cannot convert end coordinate! (line: "
                               + std::to_string(num_line) + ")!");
                reporter.report(message);
                error = true;
            }
            if (!error && start > end)
            {
                // don't check if it's already error out
                error = true;
                message.append("Error: Start coordinate should be smaller than "
                               "end coordinate!\n");
                message.append("start: " + std::to_string(start) + "\n");
                message.append("end: " + std::to_string(end) + "\n");
                throw std::runtime_error(message);
            }
            else if (error)
            {
                throw std::runtime_error("");
            }
            // the strand location is different depending on the type
            // if it is bed, then we use the STRAND index
            // if not, then we assume the format is CHR START END STRAND
            std::vector<std::string>::size_type strand_index =
                (type->second == 1) ? (+BED::STRAND) : (+BED::END + 1);

            if (token.size() > strand_index)
            {
                if (token[strand_index] == "-")
                {
                    start -= window_3;
                    if (start < 1) start = window_3;
                    end += window_5;
                }
                else if (token[strand_index].compare("+") == 0
                         || token[strand_index].compare(".") == 0)
                {
                    start -= window_5;
                    if (start < 1) start = 1;
                    end += window_3;
                }
                else
                {
                    std::string error = "Error: Undefined strand "
                                        "information. Possibly a malform "
                                        "BED file: "
                                        + token[+BED::STRAND];
                    throw std::runtime_error(error);
                }
            }
            else
            {
                start -= window_5;
                if (start < 1) start = 1;
                end += window_3;
            }
            // 1 because that is the index reserved for background
            if (chr_code < 0) continue;
            if (gene_sets.size() <= static_cast<size_t>(chr_code) + 1)
            { gene_sets.resize(static_cast<size_t>(chr_code) + 1); }
            gene_sets[static_cast<size_t>(chr_code)].add(start, end, 1);
        }
    }
    else
    {
        // gene list format
        std::vector<std::string> token;
        while (std::getline(input, line))
        {
            // read in the gene file
            misc::trim(line);
            if (line.empty()) continue;
            // this allow flexibility for both Gene Gene Gene and
            // Gene\nGene\n format
            token = misc::split(line);
            for (auto&& g : token) { msigdb_list[g].push_back(1); }
        }
    }
    input.close();
}


void Region::load_gtf(
    const std::string& gtf,
    const std::unordered_map<std::string, std::vector<int>>& msigdb_list,
    const std::vector<std::string>& feature, const uint32_t max_chr,
    const int window_5, const int window_3,
    std::vector<IITree<int, int>>& gene_sets, const bool genome_wide_background,
    const bool provided_background, Reporter& reporter)
{
    // don't bother if there's no msigdb genes and we are using genome wide
    // background

    if (msigdb_list.empty() && genome_wide_background) return;
    bool gz_input = false;
    // we want to allow gz file input (as GTF file can be big)
    GZSTREAM_NAMESPACE::igzstream gz_gtf_file;
    if (gtf.substr(gtf.find_last_of(".") + 1).compare("gz") == 0)
    {
        gz_gtf_file.open(gtf.c_str());
        if (!gz_gtf_file.good())
        {
            std::string error_message =
                "Error: Cannot open GTF (gz) to read!\n";
            throw std::runtime_error(error_message);
        }
        gz_input = true;
    }
    std::ifstream gtf_file;
    if (!gz_input)
    {
        gtf_file.open(gtf.c_str());
        if (!gtf_file.is_open())
        {
            std::string error_message = "Cannot open gtf file: " + gtf;
            throw std::runtime_error(error_message);
        }
    }
    std::istream* stream;
    if (gz_input) { stream = &(gz_gtf_file); }
    else
    {
        stream = &(gtf_file);
    }
    std::vector<std::string> token, attribute, extract;
    std::string chr, name, id, line;
    int chr_code, start, end;
    size_t num_line = 0;
    size_t exclude_feature = 0, chr_exclude = 0;

    // this should ensure we will be reading either from the gz stream or
    // ifstream
    while (std::getline(*stream, line))
    {
        misc::trim(line);
        // skip headers
        if (line.empty() || line[0] == '#') continue;
        // we don't count the header lines
        num_line++;
        misc::split(token, line, "\t");
        // convert chr string into consistent chr_coding
        chr_code = get_chrom_code_raw(token[+GTF::CHR].c_str());
        if (in_feature(token[+GTF::FEATURE], feature) && chr_code >= 0
            && chr_code <= static_cast<int>(max_chr) && token.size() == 9)
        {
            start = 0;
            end = 0;
            try
            {
                start = misc::string_to_int(token[+GTF::START].c_str());
                if (start < 0)
                {
                    // we opt for extreme stringency. Will definitely want the
                    // whole gtf file to contain valid entries
                    std::string error =
                        "Error: Negative Start Coordinate! (line: "
                        + std::to_string(num_line)
                        + "). Will ignore the gtf file\n";
                    throw std::runtime_error(error);
                }
            }
            catch (...)
            {
                std::string error =
                    "Error: Cannot convert the start coordinate! (line: "
                    + misc::to_string(num_line) + "). Will ignore the gtf file";
                throw std::runtime_error(error);
            }
            try
            {
                end = misc::string_to_int(token[+GTF::END].c_str());
                if (end < 0)
                {
                    std::string error =
                        "Error: Negative End Coordinate! (line: "
                        + misc::to_string(num_line)
                        + "). Will ignore the gtf file\n";
                    throw std::runtime_error(error);
                }
                // GTF end is inclusive
            }
            catch (...)
            {
                std::string error =
                    "Error: Cannot convert the end coordinate! (line: "
                    + std::to_string(num_line) + "). Will ignore the gtf file";
                throw std::runtime_error(error);
            }
            // Now extract the name
            bool found = parse_attribute(token[+GTF::ATTRIBUTE], id, name);
            if (!found)
            {
                // lack both
                std::string message =
                    "Error: GTF file should contain the "
                    "gene_id field. This GTF file does not contain either the "
                    "gene_id field or gene_name field. Please check if you "
                    "have "
                    "the correct file\n";
                throw std::runtime_error(message);
            }
            else if (id.empty())
            {
                // lack ID, but mandate
                std::string message = "Error: GTF file should contain the "
                                      "gene_id field. Please check if you have "
                                      "the correct file\n";
                throw std::runtime_error(message);
            }

            // it is ok to not check for previous error as they all result in
            // throw
            if (start > end)
            {
                std::string message = "Error: Start coordinate should be "
                                      "smaller than end coordinate!\n";
                message.append("start: " + misc::to_string(start) + "\n");
                message.append("end: " + misc::to_string(end) + "\n");
                throw std::runtime_error(message);
            }
            // add padding
            if (token[+GTF::STRAND].compare("-") == 0)
            {
                start -= window_3;
                if (start < 1) start = 1;
                end += window_5;
            }
            else if (token[+GTF::STRAND].compare("+") == 0
                     || token[+GTF::STRAND].compare(".") == 0)
            {
                start -= window_5;
                if (start < 1) start = 1;
                end += window_3;
            }
            else
            {
                // GTF by definition should contain the strand information. If
                // there isn't, then this is a malformed GTF file
                std::string error = "Error: Undefined strand information. "
                                    "Possibly a malform GTF file: "
                                    + token[+GTF::STRAND];
                throw std::runtime_error(error);
            }
            // now check if we can find it in the MSigDB entry
            auto&& id_search = msigdb_list.find(id);
            auto&& name_search = msigdb_list.find(name);
            if (chr_code < 0) continue;
            if (gene_sets.size() <= static_cast<size_t>(chr_code) + 1)
            { gene_sets.resize(static_cast<size_t>(chr_code) + 1); }
            if (id_search != msigdb_list.end())
            {
                for (auto&& idx : id_search->second)
                {
                    gene_sets[static_cast<size_t>(chr_code)].add(start, end,
                                                                 idx);
                }
            }
            if (name_search != msigdb_list.end())
            {
                for (auto&& idx : name_search->second)
                {
                    gene_sets[static_cast<size_t>(chr_code)].add(start, end,
                                                                 idx);
                }
            }
            if (!genome_wide_background & !provided_background)
            {
                // we want to generate the background from GTF
                // but not when a background file is provided
                // background index is 1
                gene_sets[static_cast<size_t>(chr_code)].add(start, end, 1);
            }
        }
        else if (chr_code >= 0 && chr_code <= static_cast<int>(max_chr)
                 && token.size() == 9)
        {
            exclude_feature++;
        }
        else if (chr_code >= 0 && chr_code <= static_cast<int>(max_chr))
        {
            // this is a malformed file
            std::string error_message = "Error: Malformed GTF file! GTF should "
                                        "always contain 9 columns!\n";
            throw std::runtime_error(error_message);
        }
        else
        {
            chr_exclude++;
        }
    }
    std::string message = "";
    if (!num_line)
    { throw std::runtime_error("Error: Empty GTF file detected!\n"); }
    if (exclude_feature + chr_exclude >= num_line)
    {
        throw std::runtime_error("Error: No GTF entry remain after filter by "
                                 "feature and chromosome!\n");
    }
    if (exclude_feature == 1)
    {
        message.append("A total of " + std::to_string(exclude_feature)
                       + " entry removed due to feature selection\n");
    }
    else if (exclude_feature > 1)
    {
        message.append("A total of " + std::to_string(exclude_feature)
                       + " entries removed due to feature selection\n");
    }
    if (chr_exclude > 1)
    {
        message.append(
            "A total of " + std::to_string(chr_exclude)
            + " entries removed as they are not on autosomal chromosome\n");
    }
    else if (chr_exclude == 1)
    {
        message.append(
            "A total of " + std::to_string(chr_exclude)
            + " entry removed as they are not on autosomal chromosome\n");
    }
    reporter.report(message);
}


void Region::load_snp_sets(
    std::string snp_file,
    std::unordered_map<std::string, std::vector<int>>& snp_in_sets,
    std::vector<std::string>& region_names,
    std::unordered_set<std::string>& duplicated_sets, int& set_idx,
    Reporter& reporter)
{
    std::string file_name, set_name, line, message;
    std::vector<std::string> token = misc::split(snp_file, ":");
    if (token.size() == 2)
    {
        file_name = token[0];
        set_name = token[1];
    }
    else if (token.size() == 1)
    {
        set_name = file_name = snp_file;
    }
    else
    {
        std::string error_message =
            "Error: Undefine SNP set file input format: " + snp_file;
        throw std::runtime_error(error_message);
    }
    // first check if it is a set file
    std::ifstream input;
    input.open(file_name.c_str());
    if (!input.is_open())
    {
        std::string error =
            "Error: Cannot open SNP set file: " + file_name + "\n";
        throw std::runtime_error(error);
    }
    bool is_set_file = false;
    while (std::getline(input, line))
    {
        misc::trim(line);
        token = misc::split(line);
        is_set_file = (token.size() > 1);
        break;
    }
    // we want to only use the file name
    set_name = misc::base_name(set_name);
    if (!is_set_file)
    {
        if (duplicated_sets.find(set_name) != duplicated_sets.end())
        {
            std::string message =
                "Warning: Set name of " + set_name
                + " is duplicated, this set will be ignored\n";
            reporter.report(message);
            return;
        }
        duplicated_sets.insert(set_name);
        region_names.push_back(set_name);
    }
    input.clear();
    input.seekg(0, input.beg);
    while (std::getline(input, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (is_set_file)
        {
            if (duplicated_sets.find(token[0]) != duplicated_sets.end())
            {
                std::string message = "Warning: Set name of " + token[0]
                                      + " is duplicated, it will be ignored";
                reporter.report(message);
            }
            else
            {
                duplicated_sets.insert(token[0]);
                region_names.push_back(token[0]);
                for (size_t i = 1; i < token.size(); ++i)
                { snp_in_sets[token[i]].push_back(set_idx); }
                set_idx++;
            }
        }
        else
        {
            snp_in_sets[token[0]].push_back(set_idx);
        }
    }
    input.close();
    if (!is_set_file) { set_idx++; }
}

bool Region::load_bed_regions(const std::string& bed_file,
                              std::vector<IITree<int, int>>& gene_sets,
                              const int window_5, const int window_3,
                              bool& printed_warning, const int set_idx,
                              const uint32_t max_chr,
                              std::vector<std::string>& region_names,
                              std::unordered_set<std::string> duplicated_sets,
                              Reporter& reporter)
{
    /* If we have
     *
     *  chr1   | T | A | C | C | G |
     *         | | | | | | | | | | |
     * 1-base  | 1 | 2 | 3 | 4 | 5 |
     * 0-base  0   1   2   3   4   5
     *
     * For bed file with input of
     * chr1 1 5
     * we will want the region_bound variable to be
     * 1 2 6
     * Such that any SNP with coordinate of chr1:2-5 on the 1
     * base coordinate system will be included in the analysis
     */
    // check if we have the name
    std::string file_name, set_name, line, message;
    std::vector<std::string> token = misc::split(bed_file, ":");
    int chr_code;
    if (token.size() == 2)
    {
        file_name = token[0];
        set_name = token[1];
    }
    else if (token.size() == 1)
    {
        set_name = file_name = bed_file;
    }
    else
    {
        std::string error_message =
            "Error: Undefine bed file input format: " + bed_file;
        throw std::runtime_error(error_message);
    }
    set_name = misc::base_name(set_name);
    if (duplicated_sets.find(set_name) != duplicated_sets.end())
    {
        std::string message = "Warning: Set name of " + set_name
                              + " is duplicated, it will be ignored";
        reporter.report(message);
        return false;
    }
    duplicated_sets.insert(set_name);
    region_names.push_back(set_name);
    std::ifstream input;
    input.open(file_name.c_str());
    if (!input.is_open())
    {
        std::string error = "Error: Cannot open bed file: " + file_name + "\n";
        throw std::runtime_error(error);
    }

    // now read in the file
    bool error = false, is_header = false;
    size_t num_line = 0, column_size = 0;
    while (std::getline(input, line))
    {
        is_header = false;
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        try
        {
            is_bed_line(token, column_size, is_header);
        }
        catch (const std::exception& ex)
        {
            std::string message = ex.what();
            message.append("Problematic line: " + line);
            throw std::runtime_error(message);
        }
        if (is_header) continue; // skip header

        if (token.size() < 3)
        {
            message = "Error: " + file_name
                      + " contain less than 3 columns, please check your "
                        "bed files in the correct format";
            throw std::runtime_error(message);
        }

        // skip all check later
        chr_code = get_chrom_code_raw(token[0].c_str());
        if (chr_code > static_cast<int32_t>(max_chr) || chr_code < 0) continue;
        if (token.size() <= +BED::STRAND && (window_5 > 0 || window_3 > 0)
            && (window_5 != window_3) && !printed_warning)
        {
            // as bed file does not necessarily contain the strand
            // information, we will issue a warning and assume positive
            // strand. Though for this situation, it might be better for
            // user to use the same padding for 3' and 5'
            std::string message = "Warning: You bed file does not contain "
                                  "strand information, we will assume all "
                                  "regions are on the positive strand, "
                                  "e.g. start coordinates always on the 5' "
                                  "end";
            reporter.report(message);
            printed_warning = true;
        }
        int start = 0, end = 0;
        message = "";
        try
        {
            start = misc::string_to_int(token[+BED::START].c_str());
            if (start >= 0)
                ++start; // That's because bed is 0 based
            else
            {
                message.append("Error: Negative Start Coordinate at line "
                               + misc::to_string(num_line) + "!");
                error = true;
            }
        }
        catch (...)
        {
            message.append("Error: Cannot convert start coordinate! (line: "
                           + misc::to_string(num_line) + ")!");
            error = true;
        }

        try
        {
            end = misc::string_to_int(token[+BED::END].c_str());
            if (end < 0)
            {
                message.append("Error: Negative End Coordinate at line "
                               + misc::to_string(num_line) + "!");
                error = true;
            }
        }
        catch (...)
        {
            message.append("Error: Cannot convert end coordinate! (line: "
                           + misc::to_string(num_line) + ")!");
            error = true;
        }
        if (!error && start > end)
        {
            // only do this check if there isn't an error before
            error = true;
            message.append("Error: Start coordinate should be smaller than "
                           "end coordinate!\n");
            message.append("start: " + misc::to_string(start) + "\n");
            message.append("end: " + misc::to_string(end) + "\n");
        }
        if (error)
        {
            message.append("Please check your input is correct");
            throw std::runtime_error(message);
        }
        if (token.size() > +BED::STRAND)
        {
            // we have the strand information, therefore can add the padding
            // accordingly
            if (token[+BED::STRAND] == "-")
            {
                // negative strand, so add 3' to start and 5' to end
                start -= window_3;
                if (start < 1) start = 1;
                end += window_5;
            }
            else if (token[+BED::STRAND] == "+" || token[+BED::STRAND] == ".")
            {
                // positive or unknown strand, add 5' to start and 3' to end
                start -= window_5;
                if (start < 1) start = 1;
                end += window_3;
            }
            else
            {
                std::string error = "Error: Undefined strand "
                                    "information. Possibly a malform "
                                    "BED file: "
                                    + token[+BED::STRAND];
                throw std::runtime_error(error);
            }
        }
        else
        {
            start -= window_5;
            if (start < 1) start = 1;
            end += window_3;
        }
        // now add the boundary to the region
        // will screw up if the number of set is higher than int32_t.
        // But PRSice should crash before that, right?
        if (chr_code < 0) continue;
        if (gene_sets.size() <= static_cast<size_t>(chr_code) + 1)
        { gene_sets.resize(static_cast<size_t>(chr_code) + 1); }
        gene_sets[static_cast<size_t>(chr_code)].add(start, end, set_idx);
    }
    return true;
}

void Region::load_msigdb(
    const std::string& msig,
    std::unordered_map<std::string, std::vector<int>>& msigdb_list,
    std::vector<std::string>& region_names,
    std::unordered_set<std::string> duplicated_sets, int& set_idx,
    Reporter& reporter)
{
    std::ifstream input;
    input.open(msig.c_str());
    if (!input.is_open())
    {
        std::string error = "Error: Cannot open MSigDB file: " + msig + "\n";
        throw std::runtime_error(error);
    }
    std::string line;
    std::vector<std::string> token;
    while (std::getline(input, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() < 2)
        {
            std::string message =
                "Error: Each line of MSigDB require at least 2 information: "
                + msig + "\n";
            message.append(line);
            throw std::runtime_error(message);
        }
        if (duplicated_sets.find(token[0]) != duplicated_sets.end())
        {
            std::string message = "Warning: Set name of " + token[0]
                                  + " is duplicated, it will be ignored";
            reporter.report(message);
        }
        else
        {
            duplicated_sets.insert(token[0]);
            region_names.push_back(token[0]);
            for (size_t i = 1; i < token.size(); ++i)
            { msigdb_list[token[i]].push_back(set_idx); }
            set_idx++;
        }
    }
    input.close();
}

Region::~Region() {}
