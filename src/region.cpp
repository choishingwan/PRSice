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


#include "region.hpp"
#include "genotype.hpp"

Region::Region(const std::string& exclusion_range, Reporter& reporter)
{
    // First check if it is a simple range in format of
    // chr:start-end,chr:start-end
    if (exclusion_range.empty()) return;
    std::vector<std::string> region_range = misc::split(exclusion_range, ",");
    bool file_input = false;
    if (region_range.size() == 1) {
        // this can either be: 1 range, or an input bed file
        std::vector<std::string> token = misc::split(region_range[0], ":");
        if (token.size() == 1) {
            file_input = true;
        }
    }
    if (file_input) {
        std::vector<std::string> bed;
        bed.push_back(exclusion_range);
        process_bed(bed, reporter);
    }
    else
    {
        // we will manually insert the regions
        std::vector<region_bound> current_region;
        for (auto&& range : region_range) {
            std::vector<std::string> token = misc::split(range, ":");
            if (token.size() >= 2) {
                int chr = get_chrom_code_raw(token[0].c_str());
                int temp, start, end;
                try
                {
                    std::vector<std::string> coordinates =
                        misc::split(token[1], "-");
                    temp = misc::convert<int>(coordinates[0].c_str());
                    if (temp < 0) {
                        std::string error = "Error: Negative Start Coordinate "
                                            "for exclusion range!";
                        throw std::runtime_error(error);
                    }
                    start = temp;
                    temp = misc::convert<int>(
                        coordinates[(coordinates.size() > 1)].c_str());
                    if (temp < 0) {
                        std::string error = "Error: Negative End Coordinate "
                                            "for exclusion range!";
                        throw std::runtime_error(error);
                    }
                    end = temp + (!(coordinates.size() > 1));
                }
                catch (const std::runtime_error& error)
                {
                    fprintf(stderr, "Error: Non-numeric coordinate(s)!\n");
                    throw std::runtime_error(error.what());
                }
                region_bound cur_bound;
                cur_bound.chr = chr;
                cur_bound.start = start;
                cur_bound.end = end;
                current_region.push_back(cur_bound);
            }
            else
            {
                std::string message = "Error: Invalid exclusion range format. "
                                      "Should be chr:start or chr:start-end";
                throw std::runtime_error(message);
            }
        }
        size_t i_region = m_region_list.size();
        m_region_list.push_back(solve_overlap(current_region, i_region));
    }
    m_snp_check_index = std::vector<size_t>(1, 0);
}

Region::Region(std::vector<std::string> feature, const int window_5,
               const int window_3)
    : m_5prime(window_5), m_3prime(window_3)
{
    // Make the base region which includes everything
    m_duplicated_names.insert("Base");
    m_region_name.push_back("Base");
    m_gtf_feature = feature;
    m_region_list.push_back(std::vector<region_bound>(1));
}

void Region::run(const std::string& gtf, const std::string& msigdb,
                 const std::vector<std::string>& bed,
                 const std::string& snp_set, const std::string& multi_snp_sets,
                 const Genotype& target, const std::string& out,
                 const std::string& background, Reporter& reporter)
{

    if (gtf.empty() && bed.size() == 0 && snp_set.empty()
        && multi_snp_sets.empty())
    {
        m_snp_check_index = std::vector<size_t>(m_region_name.size(), 0);
        m_region_snp_count = std::vector<int>(m_region_name.size());
        return;
    }
    if ((m_5prime > 0 || m_3prime > 0) && (m_5prime != m_3prime)) {
        std::string message = "Warning: We will assume a positive strand for "
                              "any features with unspecific strand information "
                              "e.g. \".\"";
        reporter.report(message);
    }
    process_bed(bed, reporter);
    m_out_prefix = out;
    size_t num_bed_region = m_region_list.size();
    std::unordered_map<std::string, std::set<std::string>> id_to_name;
    std::unordered_map<std::string, region_bound> gtf_boundary;
    // without the gtf file, we will not process the msigdb file
    if (!gtf.empty()) {
        reporter.report("Processing the GTF file");
        try
        {
            uint32_t max_chr = target.max_chr();
            gtf_boundary = process_gtf(gtf, id_to_name, out, max_chr, reporter);
        }
        catch (const std::runtime_error& error)
        {
            gtf_boundary.clear();
            throw std::runtime_error(error.what());
        }
        std::string message = "A total of "
                              + std::to_string(gtf_boundary.size())
                              + " genes found in the GTF file";
        reporter.report(message);
        if (gtf_boundary.size() != 0) {
            process_msigdb(msigdb, gtf_boundary, id_to_name, reporter);
        }
    }
    process_snp_sets(snp_set, multi_snp_sets, target, reporter);
    if (background.empty()) {
        generate_background(gtf_boundary, num_bed_region, reporter);
    }
    else
    {
        read_background(background, gtf_boundary, id_to_name, reporter);
    }
    m_snp_check_index = std::vector<size_t>(m_region_name.size(), 0);
    m_region_snp_count = std::vector<int>(m_region_name.size());
    m_duplicated_names.clear();
}
void Region::process_snp_sets(const std::string& single_snp_set,
                              const std::string& multi_snp_set,
                              const Genotype& target, Reporter& reporter)
{
    if (single_snp_set.empty() && multi_snp_set.empty()) return;
    std::string message = "";
    size_t i_region = 0;
    if (!single_snp_set.empty()) {
        std::ifstream input;
        input.open(single_snp_set.c_str());
        if (!input.is_open()) {
            message = "Error: " + single_snp_set + " cannot be open!";
            throw std::runtime_error(message);
        }
        if (m_duplicated_names.find(single_snp_set) != m_duplicated_names.end())
        {
            message = "Warning: Set name of " + single_snp_set
                      + " is duplicated, it will be ignored";
            reporter.report(message);
        }
        else
        {
            std::vector<region_bound> current_region;
            std::string line;
            while (std::getline(input, line)) {
                misc::trim(line);
                int chr, loc;
                if (target.get_snp_loc(line, chr, loc)) {
                    region_bound cur_bound;
                    cur_bound.chr = chr;
                    cur_bound.start = loc;
                    cur_bound.end = loc;
                    current_region.push_back(cur_bound);
                }
            }
            input.close();
            if (current_region.size() > 0) {
                i_region = m_region_list.size();
                m_region_list.push_back(
                    solve_overlap(current_region, i_region));
                m_region_name.push_back(single_snp_set);
                m_duplicated_names.insert(single_snp_set);
            }
        }
    }
    if (!multi_snp_set.empty()) {
        std::ifstream input;
        input.open(multi_snp_set.c_str());
        if (!input.is_open()) {
            message = "Error: " + multi_snp_set + " cannot be open!";
            throw std::runtime_error(message);
        }
        std::string line;
        while (std::getline(input, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            std::vector<std::string> token = misc::split(line);
            if (token.size() <= 1) {
                message = "Error: Multi-set file should contain at least 2 "
                          "columns. Did you want to use --snp-set instead?";
                // can be less stringent
                throw std::runtime_error(message);
            }
            if (m_duplicated_names.find(token[0]) != m_duplicated_names.end()) {
                message = "Warning: Set name of " + token[0]
                          + " is duplicated, it will be ignored";
                reporter.report(message);
                continue;
            }
            std::vector<region_bound> current_region;
            for (auto&& snp : token) {
                int chr, loc;
                if (target.get_snp_loc(snp, chr, loc)) {
                    region_bound cur_bound;
                    cur_bound.chr = chr;
                    cur_bound.start = loc;
                    cur_bound.end = loc;
                    current_region.push_back(cur_bound);
                }
            }
            if (current_region.size() > 0) {
                i_region = m_region_list.size();
                m_region_list.push_back(
                    solve_overlap(current_region, i_region));
                m_region_name.push_back(token[0]);
                m_duplicated_names.insert(token[0]);
            }
        }
        input.close();
    }
}

void Region::process_bed(const std::vector<std::string>& bed,
                         Reporter& reporter)
{
    // TODO: Allow user define name by modifying their input (e.g. Bed:Name
    bool print_warning = false;
    size_t i_region = 0;
    for (auto& b : bed) {
        std::string message = "Reading: " + b;
        reporter.report(message);
        std::ifstream bed_file;
        bool error = false;
        bed_file.open(b.c_str());
        if (!bed_file.is_open()) {
            message = "Warning: " + b + " cannot be open. It will be ignored";
            reporter.report(message);
            continue;
        }
        if (m_duplicated_names.find(b) != m_duplicated_names.end()) {
            message =
                "Warning: " + b + " is duplicated, it will only be read once";
            reporter.report(message);
            continue;
        }
        std::vector<region_bound> current_region;
        std::string line;
        size_t num_line = 0;
        while (std::getline(bed_file, line)) {
            num_line++;
            misc::trim(line);
            if (line.empty()) continue;
            std::vector<std::string> token = misc::split(line);
            if (token.size() < 3) {

                message = "Error: " + b
                          + " contain less than 3 columns, it will be ignored";
                reporter.report(message);
                break;
            }
            if (token.size() <= +BED::STRAND && (m_5prime > 0 || m_3prime > 0)
                && (m_5prime != m_3prime) && !print_warning)
            {
                std::string message = "Warning: You bed file does not contain "
                                      "strand information, we will assume all "
                                      "regions are on the positive strand, "
                                      "e.g. start coordinates always on the 5' "
                                      "end";
                reporter.report(message);
                print_warning = true;
            }
            int temp = 0;
            size_t start = 0, end = 0;
            message = "";
            try
            {
                temp = misc::convert<int>(token[+BED::START]);
                if (temp >= 0)
                    start = temp + 1; // That's because bed is 0 based
                else
                {

                    message.append("Error: Negative Start Coordinate at line "
                                   + std::to_string(num_line) + "!");
                    error = true;
                }
            }
            catch (const std::runtime_error& er)
            {
                message.append("Error: Cannot convert start coordinate! (line: "
                               + std::to_string(num_line) + ")!");
                error = true;
            }
            try
            {
                temp = misc::convert<int>(token[+BED::END]);
                if (temp >= 0)
                    end = temp;
                else
                {
                    message.append("Error: Negative End Coordinate at line "
                                   + std::to_string(num_line) + "!");
                    error = true;
                }
            }
            catch (const std::runtime_error& er)
            {
                message.append("Error: Cannot convert end coordinate! (line: "
                               + std::to_string(num_line) + ")!");
                error = true;
            }
            if (start > end) {
                error = true;
                message.append("Error: Start coordinate should be smaller than "
                               "end coordinate!\n");
                message.append("start: " + std::to_string(start) + "\n");
                message.append("end: " + std::to_string(end) + "\n");
            }
            if (error) {
                message.append("We will ignore this file");
                break;
            }
            // only include regions that falls into the chromosome of interest
            // ditch m_chr_order as it is not really useful in this new way of
            // handling things
            if (token.size() > +BED::STRAND) {
                if (token[+BED::STRAND].compare("-") == 0) {
                    if (start - m_3prime < 1) {
                        start = 1;
                    }
                    else
                        start -= m_3prime;
                    end += m_5prime;
                }
                else if (token[+BED::STRAND].compare("+") == 0
                         || token[+BED::STRAND].compare(".") == 0)
                {
                    if (start - m_5prime < 1) {
                        start = 1;
                    }
                    else
                        start -= m_5prime;
                    end += m_3prime;
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
                if (start - m_5prime < 1) {
                    start = 1;
                }
                else
                    start -= m_5prime;
                end += m_3prime;
            }
            region_bound cur_bound;
            cur_bound.chr = get_chrom_code_raw(token[0].c_str());
            cur_bound.start = start;
            cur_bound.end = end;
            // this should help us to avoid problem
            current_region.push_back(cur_bound);
        }

        if (!error) {
            // TODO: DEBUG This!! Might go out of scope
            i_region = m_region_list.size();
            m_region_list.push_back(solve_overlap(current_region, i_region));
            m_region_name.push_back(b);
            m_duplicated_names.insert(b);
        }
        else
        {
            throw std::runtime_error(message);
        }
        bed_file.close();
    }
}


// it return an unordered_map where the key is the gene ID, and the value is the
// boundary
std::unordered_map<std::string, Region::region_bound> Region::process_gtf(
    const std::string& gtf,
    std::unordered_map<std::string, std::set<std::string>>& id_to_name,
    const std::string& out_prefix, const uint32_t max_chr, Reporter& reporter)
{

    std::unordered_map<std::string, Region::region_bound> result_boundary;
    if (gtf.empty()) return result_boundary; // basically return an empty map

    std::string line;
    size_t num_line = 0, exclude_feature = 0;


    bool gz_input = false;
    GZSTREAM_NAMESPACE::igzstream gz_gtf_file;
    if (gtf.substr(gtf.find_last_of(".") + 1).compare("gz") == 0) {
        gz_gtf_file.open(gtf.c_str());
        if (!gz_gtf_file.good()) {
            std::string error_message =
                "Error: Cannot open GTF (gz) to read!\n";
            throw std::runtime_error(error_message);
        }
        gz_input = true;
    }

    std::ifstream gtf_file;
    if (!gz_input) {
        gtf_file.open(gtf.c_str());
        if (!gtf_file.is_open()) {
            std::string error_message = "Cannot open gtf file: " + gtf;
            throw std::runtime_error(error_message);
        }
    }
    std::vector<std::string> token, attribute, extract;
    std::string chr, name, id;
    int chr_code, temp;
    uint32_t start, end;
    while ((!gz_input && std::getline(gtf_file, line))
           || (gz_input && std::getline(gz_gtf_file, line)))
    {
        num_line++;
        misc::trim(line);
        // skip headers
        if (line.empty() || line[0] == '#') continue;
        token = misc::split(line, "\t");
        chr = token[+GTF::CHR];
        chr_code = get_chrom_code_raw(chr.c_str());
        if (in_feature(token[+GTF::FEATURE]) && chr_code <= max_chr) {
            temp = 0;
            start = 0;
            end = 0;
            try
            {
                temp = misc::convert<int>(token[+GTF::START]);
                if (temp < 0) {
                    // well, this is a bit extreme. Alternatively we can opt for
                    // just skipping problematic entries
                    std::string error =
                        "Error: Negative Start Coordinate! (line: "
                        + std::to_string(num_line)
                        + "). Will ignore the gtf file\n";
                    result_boundary.clear();
                    id_to_name.clear();
                    throw std::runtime_error(error);
                }
                start = temp;
            }
            catch (const std::runtime_error& er)
            {
                std::string error =
                    "Error: Cannot convert the start coordinate! (line: "
                    + std::to_string(num_line) + "). Will ignore the gtf file";
                result_boundary.clear();
                id_to_name.clear();
                throw std::runtime_error(error);
            }
            try
            {
                temp = misc::convert<int>(token[+GTF::END]);
                if (temp < 0) {
                    std::string error =
                        "Error: Negative End Coordinate! (line: "
                        + std::to_string(num_line)
                        + "). Will ignore the gtf file\n";
                    result_boundary.clear();
                    id_to_name.clear();
                    throw std::runtime_error(error);
                }
                end = temp;
            }
            catch (const std::runtime_error& er)
            {
                std::string error =
                    "Error: Cannot convert the end coordinate! (line: "
                    + std::to_string(num_line) + "). Will ignore the gtf file";
                result_boundary.clear();
                id_to_name.clear();
                throw std::runtime_error(error);
            }
            // Now extract the name
            attribute = misc::split(token[+GTF::ATTRIBUTE], ";");
            name = "";
            id = "";
            for (auto& info : attribute) {
                if (info.find("gene_id") != std::string::npos) {
                    extract = misc::split(info);
                    if (extract.size() > 1) {
                        // TODO: WARNING: HARD CODING HERE
                        extract[1].erase(std::remove(extract[1].begin(),
                                                     extract[1].end(), '\"'),
                                         extract[1].end());
                        id = extract[1];
                    }
                }
                else if (info.find("gene_name") != std::string::npos)
                {
                    extract = misc::split(info);
                    if (extract.size() > 1) {
                        extract[1].erase(std::remove(extract[1].begin(),
                                                     extract[1].end(), '\"'),
                                         extract[1].end());
                        name = extract[1];
                    }
                }
            }
            // the GTF only mandate the ID field, so GTF can miss out the name
            // field.
            if (!id.empty() && !name.empty()) {
                id_to_name[name].insert(id);
            }

            if (start > end) {
                std::string message = "Error: Start coordinate should be "
                                      "smaller than end coordinate!\n";
                message.append("start: " + std::to_string(start) + "\n");
                message.append("end: " + std::to_string(end) + "\n");
                result_boundary.clear();
                id_to_name.clear();
                throw std::runtime_error(message);
            }
            if (token[+GTF::STRAND].compare("-") == 0) {
                if (start - m_3prime < 1) {
                    start = 1;
                }
                else
                    start -= m_3prime;
                end += m_5prime;
            }
            else if (token[+GTF::STRAND].compare("+") == 0
                     || token[+GTF::STRAND].compare(".") == 0)
            {
                if (start - m_5prime < 1) {
                    start = 1;
                }
                else
                    start -= m_5prime;
                end += m_3prime;
            }
            else
            {
                std::string error = "Error: Undefined strand information. "
                                    "Possibly a malform GTF file: "
                                    + token[+GTF::STRAND];
                result_boundary.clear();
                id_to_name.clear();
                throw std::runtime_error(error);
            }
            // Now add the information to the map using the id
            auto&& result_search = result_boundary.find(id);
            if (result_search != result_boundary.end()) {
                if (result_search->second.chr != chr_code) {
                    std::string error =
                        "Error: Same gene occur on two separate chromosome!";
                    result_boundary.clear();
                    id_to_name.clear();
                    throw std::runtime_error(error);
                }
                result_search->second.start =
                    std::min(result_search->second.start, start);
                result_search->second.end =
                    std::max(result_search->second.end, end);
            }
            else
            {
                region_bound cur_bound;
                cur_bound.chr = chr_code;
                cur_bound.start = start;
                cur_bound.end = end;
                result_boundary[id] = cur_bound;
            }
        }
        else
        {
            exclude_feature++;
        }
    }
    std::string message = "";
    if (exclude_feature == 1) {
        message.append("A total of " + std::to_string(exclude_feature)
                       + " entry removed due to feature selection");
        reporter.report(message);
    }
    if (exclude_feature > 1) {
        message.append("A total of " + std::to_string(exclude_feature)
                       + " entries removed due to feature selection");
        reporter.report(message);
    }
    return result_boundary;
}

void Region::process_msigdb(
    const std::string& msigdb,
    const std::unordered_map<std::string, Region::region_bound>& gtf_info,
    const std::unordered_map<std::string, std::set<std::string>>& id_to_name,
    Reporter& reporter)
{
    if (msigdb.empty() || gtf_info.size() == 0) return; // Got nothing to do
    // Assume format = Name URL Gene
    std::ifstream input;
    size_t i_region = 0;
    // in theory, it should be easy for us to support multiple msigdb file.
    // but ignore that for now TODO
    input.open(msigdb.c_str());
    if (!input.is_open())
        reporter.report("Cannot open " + msigdb + ". Will skip this file");
    else
    {
        std::string line, name;
        std::vector<std::string> token;
        std::vector<region_bound> current_region;
        while (std::getline(input, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            if (token.size() < 2) // Will treat the url as gene just in case
            {
                std::string message =
                    "Error: Each line require at least 2 information\n";
                message.append(line);
                reporter.report(message);
            }
            else if (m_duplicated_names.find(token[0])
                     == m_duplicated_names.end())
            {
                name = token[0];
                current_region.clear();
                current_region.shrink_to_fit();
                current_region.reserve(token.size());
                for (auto& gene : token) {
                    auto&& gtf_search = gtf_info.find(gene);
                    if (gtf_search == gtf_info.end()) {
                        // we cannot find the gene name in the gtf information
                        // (which uses gene ID)
                        auto&& gene_name_iter = id_to_name.find(gene);
                        if (gene_name_iter != id_to_name.end()) {
                            // we found a way to convert the gene name to gene
                            // id
                            auto& gene_name = gene_name_iter->second;
                            // problem is, one gene name can correspond to
                            // multiple gene id in that case, wee will take all
                            // of them
                            for (auto&& translate : gene_name) {
                                auto&& gene_gtf = gtf_info.find(translate);
                                if (gene_gtf != gtf_info.end()) {
                                    current_region.push_back(gene_gtf->second);
                                }
                            }
                        }
                    }
                    else
                    {
                        current_region.push_back(gtf_search->second);
                    }
                }
                i_region = m_region_list.size();
                m_region_list.push_back(
                    solve_overlap(current_region, i_region));
                m_region_name.push_back(name);
                m_duplicated_names.insert(name);
            }
            else
            {
                reporter.report("Duplicated Set: " + token[0]
                                + ". It will be ignored");
            }
        }
        input.close();
    }
}

void Region::read_background(
    const std::string& background,
    const std::unordered_map<std::string, region_bound>& gtf_info,
    const std::unordered_map<std::string, std::set<std::string>>& id_to_name,
    Reporter& reporter)
{
    std::unordered_map<std::string, int> file_type{
        {"bed", 1}, {"range", 0}, {"gene", 2}};
    std::vector<std::string> background_info = misc::split(background, ":");
    if (background_info.size() != 2) {
        std::string error =
            "Error: Format of --background should be <File Name>:<File Type>";
        throw std::runtime_error(error);
    }
    auto&& type = file_type.find(background_info[1]);
    if (type == file_type.end()) {
        std::string error = "Error: Undefined file type. Supported formats are "
                            "bed, gene or range";
        throw std::runtime_error(error);
    }
    std::ifstream input;
    input.open(background_info[0].c_str());
    if (!input.is_open()) {
        std::string error =
            "Error: Cannot open background file: " + background_info[0];
        throw std::runtime_error(error);
    }
    bool print_warning = false, error = false;
    std::vector<Region::region_bound> current_bound;
    std::string line;
    if (type->second == 0 || type->second == 1) {
        // range or bed
        size_t num_line = 0;
        std::vector<std::string> token;
        while (std::getline(input, line)) {
            num_line++;
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            if (token.size() < 3) {
                std::string message =
                    "Error: " + background_info[0]
                    + " contain less than 3 columns, it will be ignored";
                throw std::runtime_error(message);
            }
            if (token.size() <= +BED::STRAND && (m_5prime > 0 || m_3prime > 0)
                && (m_5prime != m_3prime) && !print_warning)
            {
                std::string message = "Warning: You bed file does not contain "
                                      "strand information, we will assume all "
                                      "regions are on the positive strand, "
                                      "e.g. start coordinates always on the 5' "
                                      "end";
                reporter.report(message);
                print_warning = true;
            }
            int temp = 0;
            size_t start = 0, end = 0;
            std::string message = "";
            try
            {
                temp = misc::convert<int>(token[+BED::START]);
                if (temp >= 0)
                    start = temp + type->second; // That's because bed is 0
                                                 // based and range format is 1
                                                 // based
                else
                {
                    message.append("Error: Negative Start Coordinate at line "
                                   + std::to_string(num_line) + "!");
                    error = true;
                }
            }
            catch (const std::runtime_error& er)
            {
                message.append("Error: Cannot convert start coordinate! (line: "
                               + std::to_string(num_line) + ")!");
                error = true;
            }
            try
            {
                temp = misc::convert<int>(token[+BED::END]);
                if (temp >= 0)
                    end = temp + 1; // That's because bed is 0 based
                else
                {
                    message.append("Error: Negative End Coordinate at line "
                                   + std::to_string(num_line) + "!");
                    error = true;
                }
            }
            catch (const std::runtime_error& er)
            {
                message.append("Error: Cannot convert end coordinate! (line: "
                               + std::to_string(num_line) + ")!");
                error = true;
            }
            if (start > end) {
                error = true;
                message.append("Error: Start coordinate should be smaller than "
                               "end coordinate!\n");
                message.append("start: " + std::to_string(start) + "\n");
                message.append("end: " + std::to_string(end) + "\n");
            }
            if (error) break;
            // only include regions that falls into the chromosome of interest
            int strand_index =
                (type->second == 1) ? (+BED::STRAND) : (+BED::END + 1);
            if (token.size() > strand_index) {
                if (token[strand_index].compare("-") == 0) {
                    if (start - m_3prime < 1) {
                        start = 1;
                    }
                    else
                        start -= m_3prime;
                    end += m_5prime;
                }
                else if (token[strand_index].compare("+") == 0
                         || token[strand_index].compare(".") == 0)
                {
                    if (start - m_5prime < 1) {
                        start = 1;
                    }
                    else
                        start -= m_5prime;
                    end += m_3prime;
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
                if (start - m_5prime < 1) {
                    start = 1;
                }
                else
                    start -= m_5prime;
                end += m_3prime;
            }
            region_bound cur_bound;
            cur_bound.chr = get_chrom_code_raw(token[0].c_str());
            cur_bound.start = start;
            cur_bound.end = end;
            // this should help us to avoid problem
            current_bound.push_back(cur_bound);
        }
    }
    else
    {
        // gene list format

        while (std::getline(input, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            auto&& gtf_search = gtf_info.find(line);
            if (gtf_search == gtf_info.end()) {
                // we cannot find the gene name in the gtf information (which
                // uses gene ID)
                auto&& id_search = id_to_name.find(line);
                if (id_search != id_to_name.end()) {
                    // we found a way to convert the gene name to gene id
                    auto& name = id_search->second;
                    // problem is, one gene name can correspond to multiple gene
                    // id in that case, wee will take all of them
                    for (auto&& translate : name) {
                        auto&& gtf_name_search = gtf_info.find(translate);
                        if (gtf_name_search != gtf_info.end()) {
                            current_bound.push_back(gtf_name_search->second);
                        }
                    }
                }
            }
            else
            {
                current_bound.push_back(gtf_search->second);
            }
        }
    }
    input.close();
    size_t i_region = m_region_list.size();
    m_region_list.push_back(solve_overlap(current_bound, i_region));
    m_region_name.push_back("Background");
}


void Region::generate_background(
    const std::unordered_map<std::string, Region::region_bound>& gtf_info,
    const size_t num_bed_region, Reporter& reporter)
{
    // this will be very ineffective. But whatever
    std::vector<Region::region_bound> temp_storage;
    temp_storage.reserve(gtf_info.size());
    for (auto&& gtf : gtf_info) {
        temp_storage.push_back(gtf.second);
    }
    /* don't include the bed as background as they are not
     * necessary defined as genic regions
     * If users want, they should provide the bed background
     * as a separate input
     */
    /*
    for (size_t i = 0; i < num_bed_region; ++i) {
        for (size_t j = 0; j < m_region_list[i].size(); ++j) {
            temp_storage.push_back(m_region_list[i][j]);
        }
    }
    */
    size_t i_region = m_region_list.size();
    m_region_list.push_back(solve_overlap(temp_storage, i_region));
    m_region_name.push_back("Background");
}

std::vector<Region::region_bound>
Region::solve_overlap(std::vector<Region::region_bound>& current_region,
                      size_t i_region)
{
    std::sort(begin(current_region), end(current_region),
              [](region_bound const& t1, region_bound const& t2) {
                  if (t1.chr == t2.chr) {
                      if (t1.start == t2.start) return t1.end < t2.end;
                      return t1.start < t2.end;
                  }
                  else
                      return t1.chr < t2.chr;
              });
    std::vector<Region::region_bound> result;
    result.reserve(current_region.size());
    int prev_chr = -1;
    size_t prev_start = 0;
    size_t prev_end = 0;
    // optimize for human for now
    for (auto&& bound : current_region) {
        if (prev_chr == -1) {
            prev_chr = bound.chr;
            prev_start = bound.start;
            prev_end = bound.end;
        }
        else if (prev_chr != bound.chr || bound.start > prev_end)
        {
            // new region
            region_bound cur_bound;
            cur_bound.chr = prev_chr;
            cur_bound.start = prev_start;
            cur_bound.end = prev_end;
            result.push_back(cur_bound);
            prev_chr = bound.chr;
            prev_start = bound.start;
            prev_end = bound.end;
        }
        else
        {
            prev_end = bound.end;
        }
    }
    if (prev_chr != -1) {
        region_bound cur_bound;
        cur_bound.chr = prev_chr;
        cur_bound.start = prev_start;
        cur_bound.end = prev_end;
        result.push_back(cur_bound);
    }
    result.shrink_to_fit();
    return result;
}
void Region::print_file(std::string output) const
{
    std::ofstream region_out;
    region_out.open(output.c_str());
    if (!region_out.is_open()) {
        std::string error =
            "Cannot open region information file to write: " + output;
        throw std::runtime_error(error);
    }
    region_out << "Region\t#SNPs" << std::endl;
    for (size_t i_region = 0; i_region < m_region_name.size(); ++i_region) {
        region_out << m_region_name[i_region] << "\t"
                   << m_region_snp_count[i_region] << std::endl;
    }
    region_out.close();
}

Region::~Region() {}

bool Region::check_exclusion(const std::string& chr, const size_t loc)
{
    int cur_chr = get_chrom_code_raw(chr.c_str());
    if (m_region_list.empty() || m_region_list.front().empty()) return false;
    // there is only one region, so we can ignore the for loop
    size_t i_region = 0;
    size_t cur_region_size = m_region_list.front().size();
    auto&& snp_check_index = m_snp_check_index[i_region];
    while (snp_check_index < cur_region_size) {
        auto&& current_bound = m_region_list[i_region][snp_check_index];
        int region_chr = current_bound.chr;
        size_t region_start = current_bound.start;
        size_t region_end = current_bound.end;
        if (cur_chr != region_chr || region_end < loc) {
            snp_check_index++;
        }
        else if (region_start <= loc && region_end >= loc)
        {
            return true;
        }
        else if (region_start > loc)
            break;
    }
    return false;
}

void Region::update_flag(const int chr, const std::string& rs, size_t loc,
                         std::vector<uintptr_t>& flag)
{
    // we first make sure the base region is true for all SNPs
    flag[0] |= ONELU;
    m_region_snp_count[0]++;
    // note: the chr is actually the order on the m_chr_order instead of the
    // actual chromosome
    // boolean to check if we have switched to next chromosome, we only allow
    // one chromosome switches
    size_t region_size = m_region_name.size();
    size_t region_start, region_end;
    int region_chr;
    // go through each region (including backgroudn)
    for (size_t i_region = 1; i_region < region_size; ++i_region) {
        // find out how many boundaries are there within the current region
        size_t current_region_size = m_region_list[i_region].size();
        // while we still have boundary left
        // allow to go through all boundaries
        // as we assume SNPs are read in the correct order (sorted)
        // we use m_snp_check_index to track the last SNP's boundary index
        // and start our search from there so that we can skip un-necessary
        // comparison
        while (m_snp_check_index[i_region] < current_region_size) {
            // obtain the current boundary as defined by m_snp_check_index
            // with YF's test data set, we need to use rs3748592
            auto&& current_bound =
                m_region_list[i_region][m_snp_check_index[i_region]];
            region_chr = current_bound.chr;
            region_start = current_bound.start;
            region_end = current_bound.end;
            if (chr != region_chr || region_end < loc) {
                // not the same chromosome, so jump to the correct location
                m_snp_check_index[i_region]++;
            }
            else if (region_start <= loc && region_end >= loc)
            {
                // This is the region
                flag[i_region / BITCT] |= ONELU << ((i_region) % BITCT);
                m_region_snp_count[i_region]++;
                break;
            }
            else if (region_start > loc)
                break;
        }
    }
}

void Region::info(Reporter& reporter) const
{
    std::string message = "";
    if (m_region_name.size() == 1) {
        message = "1 region included";
    }
    else if (m_region_name.size() > 1)
    {
        // -1 to remove the background count, as we are not going to print the
        // background anyway
        message = "A total of " + std::to_string(m_region_name.size() - 1)
                  + " regions are included";
    }
    reporter.report(message);
}
