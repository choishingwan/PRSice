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

void Region::read_bed(const std::string& bed,
                      std::vector<IITree<size_t, size_t>>& cr)
{
    std::ifstream input_file;
    input_file.open(bed);
    std::vector<std::string> range, boundary;
    if (!input_file.is_open())
    {
        throw std::runtime_error("Error: " + bed
                                 + " cannot be open. Please check you "
                                   "have the correct input");
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

        size_t low_bound, upper_bound;
        start_end(boundary[+BED::START], boundary[+BED::END], true, low_bound,
                  upper_bound);
        ++upper_bound;
        int chr = get_chrom_code_raw(boundary[0].c_str());
        if (chr >= 0 && chr < MAX_POSSIBLE_CHROM)
        {
            // ignore any chromosome that we failed to parse (chr < 0)
            // and sex chromosomes (chr > MAX_POSSIBLE_CHROM)
            size_t chr_num = static_cast<size_t>(chr);
            if (cr.size() < chr_num + 1) { cr.resize(chr_num + 1); }
            cr[chr_num].add(low_bound, upper_bound, 0);
        }
    }
    input_file.close();
}

void Region::generate_exclusion(std::vector<IITree<size_t, size_t>>& cr,
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
            read_bed(region, cr);
        }
        else if (range.size() == 2)
        {
            int chr = get_chrom_code_raw(range[0].c_str());
            boundary = misc::split(range[1], "-");
            if (boundary.size() > 2)
            {
                throw std::runtime_error(
                    "Error: Invalid exclusion range format. "
                    "Should be chr:start, chr:start-end or a bed file\n");
            }
            size_t low_bound, upper_bound;
            start_end(boundary.front(), boundary.back(), false, low_bound,
                      upper_bound);
            // the library find overlap, which for SNP at 10
            // the boundary should be defined as 10-11 when we read in the SNP
            // here we do nothing but sainity check of the input
            ++upper_bound;
            if (chr >= 0 && chr < MAX_POSSIBLE_CHROM)
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


size_t Region::generate_regions(const size_t max_chr)
{
    // should be a fresh start each time
    m_gene_sets.clear();
    m_region_name.clear();
    m_snp_in_sets.clear();

    // we can now utilize the last field of cgranges as the index of gene
    // set of interest
    m_region_name.push_back("Base");
    m_region_name.push_back("Background");
    m_processed_sets.insert("Base");
    m_processed_sets.insert("Background");

    // don't want to output gene set info if we are not working on gene sets
    // technically, if only gtf is provided, we should also have an early
    // return. But I don't bother remodelling the unit test, so we will not do
    // an early return when only gtf is provided
    if (m_snp_set.empty() && m_bed.empty() && m_msigdb.empty() && m_gtf.empty())
    { return m_region_name.size(); }
    std::string message = "Start processing gene set information\n";
    message.append("==================================================");
    m_reporter->report(message);

    // 0 reserved for base
    // 1 reserved for background
    size_t set_idx = 2;
    for (auto&& bed_file : m_bed)
    {
        message = "Loading " + bed_file + " (BED)";
        m_reporter->report(message);
        set_idx += load_bed_regions(bed_file, set_idx, max_chr);
    }
    // SNP sets
    // we allow ambiguous snp set input
    // It all depends on the file content
    // (only one column = SNP list, mutliple column = SNP sets)
    for (auto&& s : m_snp_set)
    {
        message = "Loading " + s + " (SNP sets)";
        m_reporter->report(message);
        load_snp_sets(s, set_idx);
    }
    // the SNP sets information are stored within snp_in_sets
    if (!m_msigdb.empty() && m_gtf.empty())
    {
        std::string error_message =
            "Error: MSigDB input requires a complementary"
            "GTF file!\n";
        throw std::runtime_error(error_message);
    }
    std::unordered_map<std::string, std::vector<size_t>> msigdb_list;
    for (auto&& msig : m_msigdb)
    {
        message = "Loading " + msig + " (MSigDB)";
        m_reporter->report(message);
        try
        {
            load_msigdb(msig, msigdb_list, set_idx);
        }
        catch (const std::runtime_error& e)
        {
            throw std::runtime_error(e.what());
        }
    }
    // now process the gtf file and add the regions
    if (!m_background.empty() && !m_genome_wide_background)
    {
        // we need to read in a background file, but ignore the
        // case where we use the gtf as background as we should've
        // already done that
        message = "Loading background set from " + m_background;
        m_reporter->report(message);
        load_background(max_chr, msigdb_list);
    }
    if (!m_gtf.empty() && (!m_msigdb.empty() || !m_genome_wide_background))
    {
        // generate the region information based on the msigdb info
        // or generate for the background
        message = "Loading GTF file: " + m_gtf;
        m_reporter->report(message);
        load_gtf(msigdb_list, max_chr);
    }
    // index gene list
    for (auto&& tree : m_gene_sets) tree.index();
    // now we can go through each SNP and update their flag
    const size_t num_sets = m_region_name.size();
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
    m_reporter->report(message);
    return num_sets;
}

void Region::load_background(
    const size_t max_chr,
    std::unordered_map<std::string, std::vector<size_t>>& msigdb_list)
{
    const std::unordered_map<std::string, size_t> file_type {
        {"bed", 1}, {"range", 0}, {"gene", 2}, {"snp", 3}};
    // format of the background string should be name:format
    std::string user_type, file_name;
    if (!get_set_name(m_background, file_name, user_type))
    {
        std::string error =
            "Error: Format of --background should be <File Name>:<File Type>";
        throw std::runtime_error(error);
    }
    // check if we know the format
    auto&& type = file_type.find(user_type);
    if (type == file_type.end())
    {
        std::string error = "Error: Undefined file type. Supported formats are "
                            "bed, gene, range or snp";
        throw std::runtime_error(error);
    }
    // open the file
    std::ifstream input;
    input.open(file_name.c_str());
    if (!input.is_open())
    {
        std::string error_message =
            "Error: Cannot open background file: " + file_name + " to read\n";
        throw std::runtime_error(error_message);
    }
    std::string line;
    std::vector<std::string> token;
    bool is_header = false;
    size_t column_size = 0;
    if (type->second == 0 || type->second == 1)
    {
        // this is either a range format or a bed file
        // the only difference is range is 1 based and bed is 0 based
        size_t num_line = 0;
        while (std::getline(input, line))
        {
            num_line++;
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            is_bed_line(token, column_size, is_header);
            if (is_header) continue;

            int32_t chr_code = get_chrom_code_raw(token[0].c_str());
            if (chr_code < 0 || chr_code >= MAX_POSSIBLE_CHROM) continue;
            size_t chr = static_cast<size_t>(chr_code);
            if (chr > max_chr) continue;
            if (token.size() <= +BED::STRAND
                && (m_window_5 > 0 || m_window_3 > 0)
                && (m_window_5 != m_window_3) && !m_printed_bed_strand_warning)
            {
                std::string message = "Warning: You bed file does not contain "
                                      "strand information, we will assume all "
                                      "regions are on the positive strand, "
                                      "e.g. start coordinates always on the 5' "
                                      "end";
                m_reporter->report(message);
                m_printed_bed_strand_warning = true;
            }
            size_t start = 0, end = 0;
            std::string message = "";
            try
            {
                start_end(token[+BED::START], token[+BED::END], type->second,
                          start, end);
            }
            catch (const std::runtime_error& e)
            {
                message.append("Error: Invalid " + std::string(e.what())
                               + " coordinate! (line: "
                               + std::to_string(num_line) + ")!");
                throw std::runtime_error(message);
            }
            catch (const std::logic_error& e)
            {
                throw std::runtime_error(e.what());
            }
            // the strand location is different depending on the type
            // if it is bed, then we use the STRAND index
            // if not, then we assume the format is CHR START END STRAND
            size_t strand_index =
                (type->second == 1) ? (+BED::STRAND) : (+BED::END + 1);

            if (token.size() > strand_index)
            { extend_region(start, end, token[strand_index]); }
            else
            {
                extend_region(start, end, ".");
            }
            // 1 because that is the index reserved for background
            if (m_gene_sets.size() <= chr + 1) { m_gene_sets.resize(chr + 1); }
            m_gene_sets[chr].add(start, end, 1);
        }
    }
    else if (type->second == 2)
    {
        // gene list format
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
    else
    {
        // this is a file contain all the SNPs to be included for background
        // doesn't matter what format it is
        while (std::getline(input, line))
        {
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            for (auto& s : token) { m_snp_in_sets[s].push_back(1); }
        }
    }
    input.close();
}


void Region::load_gtf(
    const std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
    const size_t max_chr)
{
    // don't bother if there's no msigdb genes and we are using genome wide
    // background
    const bool provided_background = !m_background.empty();
    if (msigdb_list.empty() && m_genome_wide_background) return;
    bool gz_input = misc::is_gz_file(m_gtf);
    // we want to allow gz file input (as GTF file can be big)
    GZSTREAM_NAMESPACE::igzstream gz_gtf_file;
    std::ifstream gtf_file;
    if (gz_input)
    {
        gz_gtf_file.open(m_gtf.c_str());
        if (!gz_gtf_file.good())
        {
            std::string error_message =
                "Error: Cannot open GTF (gz) to read!\n";
            throw std::runtime_error(error_message);
        }
    }
    else
    {
        gtf_file.open(m_gtf.c_str());
        if (!gtf_file.is_open())
        {
            std::string error_message = "Cannot open gtf file: " + m_gtf;
            throw std::runtime_error(error_message);
        }
    }
    std::istream* stream;
    if (gz_input) { stream = &(gz_gtf_file); }
    else
    {
        stream = &(gtf_file);
    }
    std::vector<std::string> token(9), attribute, extract;
    std::string chr_str, name, id, line;
    int chr_code;
    size_t chr, start, end;
    size_t num_line = 0;
    size_t exclude_feature = 0, chr_exclude = 0;

    // this should ensure we will be reading either from the gz stream or
    // ifstream
    while (std::getline(*stream, line))
    {
        misc::trim(line);
        // skip headers
        if (line.empty() || line[0] == '#') continue;
        ++num_line;
        misc::split(token, line, "\t");
        if (token.size() != +GTF::MAX)
        {
            std::string error_message = "Error: Malformed GTF file! GTF should "
                                        "always contain 9 columns!\n";
            throw std::runtime_error(error_message);
        }
        else if (in_feature(token[+GTF::FEATURE], m_feature))
        {
            // convert chr string into consistent chr_coding
            chr_code = get_chrom_code_raw(token[+GTF::CHR].c_str());
            chr = static_cast<size_t>(chr_code);
            if (chr_code < 0 || chr_code >= MAX_POSSIBLE_CHROM || chr > max_chr)
            {
                ++chr_exclude;
                continue;
            }
            start = 0;
            end = 0;
            try
            {
                start_end(token[+GTF::START], token[+GTF::END], 0, start, end);
            }
            catch (const std::runtime_error& e)
            {
                std::string error =
                    "Error: Cannot Invalid " + std::string(e.what())
                    + " coordinate! (line: " + misc::to_string(num_line)
                    + "). Will ignore the gtf file";
                throw std::runtime_error(error);
            }
            catch (const std::logic_error& e)
            {
                throw std::runtime_error(e.what());
            }

            // Now extract the name
            parse_attribute(token[+GTF::ATTRIBUTE], id, name);
            if (id.empty())
            {
                // lack ID, but mandate
                std::string message = "Error: GTF file should contain the "
                                      "gene_id field. Please check if you have "
                                      "the correct file\n";
                throw std::runtime_error(message);
            }
            // add padding
            extend_region(start, end, token[+GTF::STRAND]);
            // now check if we can find it in the MSigDB entry
            auto&& id_search = msigdb_list.find(id);
            if (m_gene_sets.size() <= chr + 1) { m_gene_sets.resize(chr + 1); }
            if (id_search != msigdb_list.end())
            {
                for (auto&& idx : id_search->second)
                { m_gene_sets[chr].add(start, end, idx); }
            }
            if (!name.empty())
            {
                auto&& name_search = msigdb_list.find(name);
                if (name_search != msigdb_list.end())
                {
                    for (auto&& idx : name_search->second)
                    { m_gene_sets[chr].add(start, end, idx); }
                }
            }
            if (!m_genome_wide_background && !provided_background)
            {
                // we want to generate the background from GTF
                // but not when a background file is provided
                // background index is 1
                m_gene_sets[chr].add(start, end, 1);
            }
        }
        else
        {
            ++exclude_feature;
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
    if (exclude_feature > 0)
    {
        std::string entry = (exclude_feature == 1) ? "entry" : "entries";
        message.append("A total of " + std::to_string(exclude_feature) + " "
                       + entry + " removed due to feature selection\n");
    }
    if (chr_exclude > 0)
    {
        std::string entry = (chr_exclude == 1) ? "entry" : "entries";
        message.append("A total of " + std::to_string(chr_exclude) + " " + entry
                       + " removed as they are not on autosomal chromosome\n");
    }
    m_reporter->report(message);
}


void Region::load_snp_sets(const std::string& snp_file, size_t& set_idx)
{
    std::string file_name, set_name, line, message;
    bool user_set_input = get_set_name(snp_file, file_name, set_name);
    // first check if it is a set file
    std::ifstream input;
    input.open(file_name.c_str());
    if (!input.is_open())
    {
        throw std::runtime_error("Error: Cannot open SNP set file: " + file_name
                                 + "\n");
    }
    bool is_set_file = false;
    std::vector<std::string> token;
    while (std::getline(input, line))
    {
        misc::trim(line);
        token = misc::split(line);
        is_set_file = (token.size() > 1);
        break;
    }
    if (!is_set_file)
    {
        // we want to only use the file name
        if (!user_set_input) set_name = misc::base_name(set_name);
        duplicated_set(set_name);
    }
    else if (user_set_input)
    {
        throw std::runtime_error(
            "Error: Undefine SNP set file input format: " + snp_file
            + ". Set name is not allowed for multi-set SNP file");
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
            if (!duplicated_set(token[0]))
            {
                for (auto&& snp : token)
                { m_snp_in_sets[snp].push_back(set_idx); }
                ++set_idx;
            }
        }
        else
        {
            m_snp_in_sets[token[0]].push_back(set_idx);
        }
    }
    input.close();
    if (!is_set_file) { ++set_idx; }
}

bool Region::load_bed_regions(const std::string& bed_file, const size_t set_idx,
                              const size_t max_chr)
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
    int chr_code;
    bool user_set_name = get_set_name(bed_file, file_name, set_name);
    if (!user_set_name) set_name = misc::base_name(set_name);
    if (duplicated_set(set_name)) { return false; }
    std::ifstream input;
    input.open(file_name.c_str());
    if (!input.is_open())
    {
        std::string error = "Error: Cannot open bed file: " + file_name + "\n";
        throw std::runtime_error(error);
    }

    // now read in the file
    std::vector<std::string> token;
    bool is_header = false;
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

        // skip all check later
        chr_code = get_chrom_code_raw(token[0].c_str());
        if (chr_code < 0 || chr_code >= MAX_POSSIBLE_CHROM) continue;
        size_t chr = static_cast<size_t>(chr_code);
        if (chr > max_chr) continue;

        if (token.size() <= +BED::STRAND && (m_window_5 > 0 || m_window_3 > 0)
            && (m_window_5 != m_window_3) && !m_printed_bed_strand_warning)
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
            m_reporter->report(message);
            m_printed_bed_strand_warning = true;
        }
        size_t start = 0, end = 0;
        message = "";
        try
        {
            start_end(token[+BED::START], token[+BED::END], 1, start, end);
        }
        catch (const std::runtime_error& e)
        {
            message.append("Error: Invalid " + std::string(e.what())
                           + " coordinate at line "
                           + misc::to_string(num_line));
            message.append("Please check your input is correct");
            throw std::runtime_error(message);
        }
        catch (const std::logic_error& e)
        {
            message = e.what();
            message.append("Please check your input is correct");
            throw std::runtime_error(message);
        }

        if (token.size() > +BED::STRAND)
        { extend_region(start, end, token[+BED::STRAND]); }
        else
        {
            extend_region(start, end, ".");
        }
        if (m_gene_sets.size() <= chr + 1) { m_gene_sets.resize(chr + 1); }
        m_gene_sets[chr].add(start, end, set_idx);
    }
    return true;
}

void Region::load_msigdb(
    const std::string& msig,
    std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
    size_t& set_idx)
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
        if (!duplicated_set(token[0]))
        {
            for (size_t i = 1; i < token.size(); ++i)
            { msigdb_list[token[i]].push_back(set_idx); }
            ++set_idx;
        }
    }
    input.close();
}

Region::~Region() {}
