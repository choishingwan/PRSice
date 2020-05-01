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

// This is for exclusion region
// end boundary is inclusive
// e.g. bound with chr1:1-10 will remove any SNPs on chr1 with coordinate
// from 1 to 10

void Region::read_bed(std::unique_ptr<std::istream> bed,
                      std::vector<IITree<size_t, size_t>>& cr,
                      bool& print_bed_strand_warning, const size_t wind_5,
                      const size_t wind_3, const size_t max_chr,
                      const size_t set_idx, const bool ZERO_BASED)
{
    print_bed_strand_warning = false;
    std::string line;
    size_t column_size = 0, low_bound, upper_bound;
    std::vector<std::string_view> range, boundary;
    size_t num_line = 0;
    while (std::getline(*bed, line))
    {
        ++num_line;
        misc::trim(line);
        boundary = misc::tokenize(line);
        try
        {
            // skip header
            if (is_bed_header(boundary, column_size)) continue;
        }
        catch (const std::exception& ex)
        {
            throw std::runtime_error(std::string(ex.what())
                                     + "Problematic line: " + line);
        }

        if (boundary.size() <= +BED::STRAND && (wind_5 > 0 || wind_3 > 0)
            && (wind_5 != wind_3))
        {
            // as bed file does not necessarily contain the strand
            // information, we will issue a warning and assume positive
            // strand. Though for this situation, it might be better for
            // user to use the same padding for 3' and 5'
            print_bed_strand_warning = true;
        }
        try
        {
            std::tie(low_bound, upper_bound) = start_end(
                boundary[+BED::START], boundary[+BED::END], ZERO_BASED);
        }
        catch (const std::runtime_error& e)
        {
            throw std::runtime_error(std::string(e.what())
                                     + "Please check your input is correct");
        }
        if (boundary.size() > +BED::STRAND)
        {
            extend_region(boundary[+BED::STRAND], wind_5, wind_3, low_bound,
                          upper_bound);
        }
        else
        {
            extend_region(".", wind_5, wind_3, low_bound, upper_bound);
        }
        auto chr = Genotype::get_chrom_code(boundary[+BED::CHR]);
        if (chr >= 0 && chr < MAX_POSSIBLE_CHROM
            && chr <= static_cast<int32_t>(max_chr))
        {
            // ignore any chromosome that we failed to parse (chr < 0)
            // and sex chromosomes (chr > MAX_POSSIBLE_CHROM)
            size_t chr_num = static_cast<size_t>(chr);
            if (cr.size() < chr_num + 1) { cr.resize(chr_num + 1); }
            cr[chr_num].add(low_bound, upper_bound, set_idx);
        }
    }
    bed.reset();
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
    bool dummy;
    const bool ZERO_BASED = true;
    for (auto&& region : exclude_regions)
    {
        range = misc::split(region, ":");
        if (range.size() == 1)
        {
            // bed file input
            auto bed_file = misc::load_stream(region);
            read_bed(std::move(bed_file), cr, dummy);
        }
        else if (range.size() == 2)
        {
            boundary = misc::split(range[1], "-");
            if (boundary.size() > 2)
            {
                throw std::runtime_error(
                    "Error: Invalid exclusion range format. "
                    "Should be chr:start, chr:start-end or a bed file\n");
            }
            auto [low_bound, upper_bound] =
                start_end(boundary.front(), boundary.back(), !ZERO_BASED);
            int chr = Genotype::get_chrom_code(range[0]);
            if (chr < 0)
            {
                throw std::runtime_error(
                    "Error: Invalid chromosome --x-range input\n");
            }
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
            throw std::runtime_error(
                "Error: Invalid exclusion range format. "
                "Should be chr:start, chr:start-end or a bed file\n");
        }
    }
    // also index the tree
    for (auto&& tree : cr) { tree.index(); }
}


size_t Region::generate_regions(
    const std::unordered_map<std::string, size_t>& included_snp_idx,
    const std::vector<SNP>& included_snps, const size_t max_chr)
{
    // should be a fresh start each time
    m_gene_sets.clear();
    m_region_name.clear();

    // we can now utilize the last field of cgranges as the index of gene
    // set of interest
    m_region_name = {"Base", "Background"};
    m_processed_sets = {"Base", "Background"};
    if (m_snp_set.empty() && m_bed.empty() && m_msigdb.empty())
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
        load_snp_sets(included_snp_idx, included_snps, s, set_idx);
    }
    // the SNP sets information are stored within snp_in_sets
    if (!m_msigdb.empty() && m_gtf.empty())
    {
        throw std::runtime_error("Error: MSigDB input requires a complementary"
                                 "GTF file!\n");
    }
    std::unordered_map<std::string, std::vector<size_t>> msigdb_list;
    for (auto&& msig : m_msigdb)
    {
        message = "Loading " + msig + " (MSigDB)";
        m_reporter->report(message);
        try
        {
            auto msig_file = misc::load_stream(msig);
            load_msigdb(msigdb_list, std::move(msig_file), set_idx);
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
        load_background(included_snp_idx, included_snps, max_chr, msigdb_list);
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
    const std::unordered_map<std::string, size_t>& snp_list_idx,
    const std::vector<SNP>& snp_list, const size_t max_chr,
    std::unordered_map<std::string, std::vector<size_t>>& msigdb_list)
{
    const std::unordered_map<std::string, size_t> file_type {
        {"bed", 1}, {"range", 0}, {"gene", 2}, {"snp", 3}};
    // format of the background string should be name:format
    auto [file_name, user_type, valid] = get_set_name(m_background);
    if (!valid)
    {
        throw std::runtime_error(
            "Error: Format of --background should be <File Name>:<File Type>");
    }
    // check if we know the format
    misc::to_lower(user_type);
    auto&& type = file_type.find(user_type);
    if (type == file_type.end())
    {
        throw std::runtime_error(
            "Error: Undefined file type. Supported formats are "
            "bed, gene, range or snp");
    }
    // open the file
    auto input = misc::load_stream(file_name);
    std::string line;
    std::vector<std::string> token;
    const size_t BACKGROUND_IDX = 1;
    if (type->second == 0 || type->second == 1)
    {
        bool print_bed_warning;
        // this is either a range format or a bed file
        // the only difference is range is 1 based and bed is 0 based
        read_bed(std::move(input), m_gene_sets, print_bed_warning, m_window_5,
                 m_window_3, max_chr, BACKGROUND_IDX, type->second);
        if (print_bed_warning)
        {
            m_reporter->report("Warning: " + file_name
                               + " for background does not contain "
                                 "strand information, we will assume all "
                                 "regions are on the positive strand, "
                                 "e.g. start coordinates always on the 5' "
                                 "end");
        }
    }
    else if (type->second == 2)
    {
        // gene list format
        while (std::getline(*input, line))
        {
            // read in the gene file
            misc::trim(line);
            if (line.empty()) continue;
            // this allow flexibility for both Gene Gene Gene and
            // Gene\nGene\n format
            token = misc::split(line);
            for (auto&& g : token) { msigdb_list[g].push_back(BACKGROUND_IDX); }
        }
        input.reset();
    }
    else
    {
        // this is a file contain all the SNPs to be included for background
        // doesn't matter what format it is
        size_t chr_num, start, end;
        while (std::getline(*input, line))
        {
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            for (auto& s : token)
            {
                auto snp_idx = snp_list_idx.find(s);
                if (snp_idx != snp_list_idx.end())
                {
                    auto&& cur_snp = snp_list[snp_idx->second];
                    chr_num = cur_snp.chr();
                    start = cur_snp.loc();
                    end = cur_snp.loc() + 1;
                    if (m_gene_sets.size() < chr_num + 1)
                    { m_gene_sets.resize(chr_num + 1); }
                    m_gene_sets[chr_num].add(start, end, BACKGROUND_IDX);
                }
            }
        }
        input.reset();
    }
}

std::tuple<size_t, size_t, size_t> Region::transverse_gtf(
    const std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
    const std::streampos file_length, const size_t max_chr, const bool gz_input,
    std::unique_ptr<std::istream> gtf_stream)
{
    const bool add_background =
        !m_genome_wide_background && m_background.empty();
    std::vector<std::string_view> token(9), attribute, extract;
    std::string chr_str, name, id, line;
    int chr_code;
    size_t chr, start, end;
    size_t num_line = 0;
    size_t exclude_feature = 0, chr_exclude = 0;
    const bool ZERO_BASED = true;
    const size_t BACKGROUND_IDX = 1;
    // this should ensure we will be reading either from the gz stream or
    // ifstream
    double progress, prev_progress = 0.0;
    while (std::getline(*gtf_stream, line))
    {
        if (!gz_input)
        {
            progress = static_cast<double>(gtf_stream->tellg())
                       / static_cast<double>(file_length) * 100;
            if (!m_reporter->unit_testing() && progress - prev_progress > 0.01)
            {
                fprintf(stderr, "\rReading %03.2f%%", progress);
                prev_progress = progress;
            }
        }
        misc::trim(line);
        // skip headers
        if (line.empty() || line[0] == '#') continue;
        ++num_line;
        token = misc::tokenize(line, "\t");
        if (token.size() != +GTF::MAX)
        {
            throw std::runtime_error("Error: Malformed GTF file! GTF should "
                                     "always contain 9 columns!\n");
        }
        if (!in_feature(token[+GTF::FEATURE], m_feature))
        {
            ++exclude_feature;
            continue;
        }
        // convert chr string into consistent chr_coding
        chr_code = Genotype::get_chrom_code(token[+GTF::CHR]);
        chr = static_cast<size_t>(chr_code);
        if (chr_code < 0 || chr_code >= MAX_POSSIBLE_CHROM || chr > max_chr)
        {
            ++chr_exclude;
            continue;
        }
        try
        {
            std::tie(start, end) =
                start_end(token[+GTF::START], token[+GTF::END], !ZERO_BASED);
        }
        catch (const std::runtime_error& e)
        {
            throw std::runtime_error(std::string(e.what())
                                     + "Will ignore the gtf file");
        }
        extend_region(token[+GTF::STRAND], m_window_5, m_window_3, start, end);

        // Now extract the name
        parse_attribute(token[+GTF::ATTRIBUTE], id, name);
        if (id.empty())
        {
            // lack ID, but mandate
            throw std::runtime_error("Error: GTF file should contain the "
                                     "gene_id field. Please check if you have "
                                     "the correct file\n");
        }
        // now check if we can find it in the MSigDB entry
        if (m_gene_sets.size() <= chr + 1) { m_gene_sets.resize(chr + 1); }
        if (!add_gene_region_from_gtf(msigdb_list, id, chr, start, end))
        { add_gene_region_from_gtf(msigdb_list, name, chr, start, end); }
        if (add_background)
        { m_gene_sets[chr].add(start, end, BACKGROUND_IDX); }
    }
    return {num_line, exclude_feature, chr_exclude};
}

void Region::load_gtf(
    const std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
    const size_t max_chr)
{
    // don't bother if there's no msigdb genes and we are using genome wide
    // background
    if (msigdb_list.empty() && m_genome_wide_background) return;
    bool gz_input = false;
    auto stream = misc::load_stream(m_gtf, gz_input);
    std::streampos file_length = 0;
    if (!gz_input)
    {
        stream->seekg(0, stream->end);
        file_length = stream->tellg();
        stream->clear();
        stream->seekg(0, stream->beg);
    }
    auto [num_line, exclude_feature, chr_exclude] = transverse_gtf(
        msigdb_list, file_length, max_chr, gz_input, std::move(stream));
    if (!m_reporter->unit_testing())
    { fprintf(stderr, "\rReading %03.2f%%\n", 100.0); }
    if (!num_line)
    { throw std::runtime_error("Error: Empty GTF file detected!\n"); }
    if (exclude_feature + chr_exclude >= num_line)
    {
        throw std::runtime_error("Error: No GTF entry remain after filter by "
                                 "feature and chromosome!\n");
    }
    std::string message = "";
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

void Region::transverse_snp_file(
    const std::unordered_map<std::string, size_t>& snp_list_idx,
    const std::vector<SNP>& snp_list, const bool is_set_file,
    std::unique_ptr<std::istream> input, size_t& set_idx)
{
    std::string line;
    std::vector<std::string> token;
    size_t chr_num, low_bound, upper_bound;
    while (std::getline(*input, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (is_set_file)
        {
            if (!duplicated_set(token[0]))
            {
                for (auto&& snp : token)
                {
                    auto snp_idx = snp_list_idx.find(snp);
                    if (snp_idx != snp_list_idx.end())
                    {
                        auto&& cur_snp = snp_list[snp_idx->second];
                        chr_num = cur_snp.chr();
                        low_bound = cur_snp.loc();
                        upper_bound = cur_snp.loc() + 1;
                        if (m_gene_sets.size() < chr_num + 1)
                        { m_gene_sets.resize(chr_num + 1); }
                        m_gene_sets[chr_num].add(low_bound, upper_bound,
                                                 set_idx);
                    }
                }
                ++set_idx;
            }
        }
        else
        {
            auto snp_idx = snp_list_idx.find(token.front());
            if (snp_idx != snp_list_idx.end())
            {
                auto&& cur_snp = snp_list[snp_idx->second];
                chr_num = cur_snp.chr();
                low_bound = cur_snp.loc();
                upper_bound = cur_snp.loc();
                if (m_gene_sets.size() < chr_num + 1)
                { m_gene_sets.resize(chr_num + 1); }
                m_gene_sets[chr_num].add(low_bound, upper_bound, set_idx);
            }
        }
    }
    input.reset();
}
void Region::load_snp_sets(
    const std::unordered_map<std::string, size_t>& snp_list_idx,
    const std::vector<SNP>& snp_list, const std::string& snp_file,
    size_t& set_idx)
{
    std::string line, message;
    auto [file_name, set_name, user_input] = get_set_name(snp_file);
    // first check if it is a set file
    auto input = misc::load_stream(file_name);
    bool is_set_file = false;
    std::vector<std::string_view> token;
    while (std::getline(*input, line))
    {
        misc::trim(line);
        token = misc::tokenize(line);
        is_set_file = (token.size() > 1);
        break;
    }
    if (!is_set_file) { duplicated_set(set_name); }
    else if (user_input)
    {
        throw std::runtime_error(
            "Error: Undefine SNP set file input format: " + snp_file
            + ". Set name is not allowed for multi-set SNP file");
    }

    input->clear();
    input->seekg(0, input->beg);
    transverse_snp_file(snp_list_idx, snp_list, is_set_file, std::move(input),
                        set_idx);
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
    std::string line, message;
    auto [file_name, set_name, user_input] = get_set_name(bed_file);
    if (duplicated_set(set_name)) { return false; }
    auto input = misc::load_stream(file_name);
    // now read in the file
    std::vector<std::string_view> token;
    bool print_warning = false;
    read_bed(std::move(input), m_gene_sets, print_warning, m_window_5,
             m_window_3, max_chr, set_idx);
    if (print_warning)
    {
        m_reporter->report("Warning: " + file_name
                           + " does not contain "
                             "strand information, we will assume all "
                             "regions are on the positive strand, "
                             "e.g. start coordinates always on the 5' "
                             "end");
    }
    return true;
}

void Region::load_msigdb(
    std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
    std::unique_ptr<std::istream> input, size_t& set_idx)
{
    std::string line;
    std::vector<std::string> token;
    while (std::getline(*input, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() < 2)
        {
            throw std::runtime_error("Error: Each line of MSigDB require "
                                     "at least 2 information: "
                                     + line);
        }
        if (!duplicated_set(token[0]))
        {
            for (size_t i = 1; i < token.size(); ++i)
            { msigdb_list[token[i]].push_back(set_idx); }
            ++set_idx;
        }
    }
    input.reset();
}

Region::~Region() {}
