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

#include "genotype.hpp"

void Genotype::build_clump_windows()
{
    // should sort w.r.t reference
    std::sort(begin(m_existed_snps), end(m_existed_snps),
              [](SNP const& t1, SNP const& t2) {
                  if (t1.chr() == t2.chr())
                  {
                      if (t1.loc() < t2.loc())
                      {
                          if (t1.ref_file_name() == t2.ref_file_name())
                          { return t1.ref_byte_pos() < t2.ref_byte_pos(); }
                          return t1.ref_file_name() < t2.ref_file_name();
                      }
                      else
                          return (t1.loc() < t2.loc());
                  }
                  else
                      return (t1.chr() < t2.chr());
              });
    // we do it here such that the m_existed_snps is sorted correctly
    // low_bound is where the current snp should read from and last_snp is where
    // the last_snp in the vector which doesn't have the up_bound set
    size_t low_bound = 0, last_snp = 0, prev_loc = 0, vector_index = 0,
           diff = 0;
    size_t prev_chr = ~size_t(0);
    bool first_snp = true;
    m_max_window_size = 0;
    // now we iterate thorugh all the SNPs to define the clumping window
    for (auto&& cur_snp : m_existed_snps)
    {
        if (first_snp || prev_chr != cur_snp.chr())
        {
            prev_chr = cur_snp.chr();
            prev_loc = cur_snp.loc();
            low_bound = vector_index;
            first_snp = false;
        }
        // we can safely assume current location always bigger than prev_loc
        // as we have sorted the vector
        else if (cur_snp.loc() - prev_loc > m_clump_distance)
        {
            // now the chromosome didn't change, and the distance of our current
            // SNP is further away from the previous SNP than our required
            // threshold
            while (cur_snp.loc() - prev_loc > m_clump_distance
                   && low_bound < vector_index)
            {
                ++low_bound;
                prev_loc = m_existed_snps[low_bound].loc();
            }
        }

        // now low_bound should be the first SNP where the core index SNP need
        // to read from
        cur_snp.set_low_bound(low_bound);
        // set the end of the vector as the default up bound
        cur_snp.set_up_bound(m_existed_snps.size());
        // update all previous SNPs that are out bounud
        while ((m_existed_snps[last_snp].chr() != cur_snp.chr())
               || (cur_snp.loc() - m_existed_snps[last_snp].loc()
                   > m_clump_distance))
        {
            // if the last SNP is on a differenet chromosome or it is to far
            // from the current SNP
            diff = vector_index - m_existed_snps[last_snp].low_bound();
            // we will set the up bound of that SNP to the current SNP
            m_existed_snps[last_snp].set_up_bound(vector_index);
            ++last_snp;
            if (m_max_window_size < diff) { m_max_window_size = diff; }
        }
        // assign the index
        ++vector_index;
    }
    size_t idx = m_existed_snps.size() - 1;
    while (true)
    {
        auto&& cur_snp = m_existed_snps[idx];
        if (cur_snp.up_bound() != m_existed_snps.size()) break;
        if (m_max_window_size < cur_snp.up_bound() - cur_snp.low_bound())
        { m_max_window_size = cur_snp.up_bound() - cur_snp.low_bound(); }
        if (idx == 0) break;
        --idx;
    }
}
// std::mutex Genotype::m_mutex;
std::vector<std::string> Genotype::set_genotype_files(const std::string& prefix)
{
    std::vector<std::string> genotype_files;
    if (prefix.find("#") != std::string::npos)
    {
        for (size_t chr = 1; chr <= m_autosome_ct; ++chr)
        {
            std::string name = prefix;
            misc::replace_substring(name, "#", std::to_string(chr));
            genotype_files.push_back(name);
        }
    }
    else
    {
        genotype_files.push_back(prefix);
    }
    return genotype_files;
}

void Genotype::add_flags(
    const std::vector<IITree<size_t, size_t>>& gene_sets,
    const std::unordered_map<std::string, std::vector<size_t>>& snp_in_sets,
    const size_t num_sets, const bool genome_wide_background)
{
    const size_t num_snps = m_existed_snps.size();
    const size_t required_size = BITCT_TO_WORDCT(num_sets);
    size_t chr, bp;
    std::vector<uintptr_t> flag(required_size, 0);
    std::unordered_map<std::string, std::vector<size_t>>::const_iterator
        snp_idx;
    for (size_t i = 0; i < num_snps; ++i)
    {
        auto&& snp = m_existed_snps[i];
        chr = snp.chr();
        bp = snp.loc();
        // if we want more speed, we can move b and max_b from the construct
        // flag function to avoid re-allocation of memory. Also for the flag
        // structure as that can almost always be reused (we copy when we set
        // flag)
        construct_flag(snp.rs(), gene_sets, snp_in_sets, flag, required_size,
                       chr, bp, genome_wide_background);
        m_existed_snps[i].set_flag(num_sets, flag);
    }
}

void Genotype::read_base(
    const std::string& base_file, const std::vector<size_t>& col_index,
    const std::vector<bool>& has_col, const std::vector<double>& barlevels,
    const double& bound_start, const double& bound_inter,
    const double& bound_end, const std::string& exclude_snps,
    const std::string& extract_snps,
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const double& maf_control, const double& maf_case,
    const double& info_threshold, const bool maf_control_filter,
    const bool maf_case_filter, const bool info_filter, const bool fastscore,
    const bool no_full, const bool is_beta, const bool is_index,
    const bool keep_ambig, Reporter& reporter)
{
    // can assume region is of the same order as m_existed_snp
    // because they use the same chr encoding and similar sorting algorithm
    // (chr, bp)
    // doesn't matter if two SNPs have the same coordinates
    assert(col_index.size() == +BASE_INDEX::MAX + 1);
    assert(!barlevels.empty());
    if (!m_is_ref)
    {
        if (!extract_snps.empty())
        {
            m_exclude_snp = false;
            m_snp_selection_list = load_snp_list(extract_snps, reporter);
        }
        else if (!exclude_snps.empty())
        {
            m_snp_selection_list = load_snp_list(exclude_snps, reporter);
        }
    }
    const size_t max_index = col_index[+BASE_INDEX::MAX];
    GZSTREAM_NAMESPACE::igzstream gz_snp_file;
    std::ifstream snp_file;
    std::ofstream mismatch_snp_record;
    // Read in threshold information
    // barlevel is sorted, last element is the biggest
    const double max_threshold =
        no_full ? (fastscore ? barlevels.back() : bound_end) : 1.0;
    // Start reading the base file. If the base file contain gz as its suffix,
    // we will read it as a gz file
    std::vector<std::string> token;
    std::string line;
    std::string message = "Base file: " + base_file + "\n";
    // Some QC counts
    std::string rs_id;
    std::string ref_allele;
    std::string alt_allele;
    double maf = 1;
    double maf_case_temp = 1;
    double info_score = 1;
    double pvalue = 2.0;
    double stat = 0.0;
    double pthres = 0.0;
    size_t chr;
    size_t loc = 0;
    size_t num_duplicated = 0;
    size_t num_excluded = 0;
    size_t num_selected = 0;
    size_t num_region_exclude = 0;
    size_t num_ambiguous = 0;
    size_t num_haploid = 0;
    size_t num_not_converted = 0; // this is for NA
    size_t num_negative_stat = 0;
    size_t num_line_in_base = 0;
    size_t num_info_filter = 0;
    size_t num_chr_filter = 0;
    size_t num_maf_filter = 0;
    bool maf_filtered = false;
    bool exclude = false;
    bool to_remove = false;
    int category = -1;
    int32_t chr_code;
    std::streampos file_length = 0;
    bool gz_input = false;
    try
    {
        gz_input = misc::is_gz_file(base_file);
    }
    catch (const std::runtime_error& e)
    {
        throw std::runtime_error(e.what());
    }

    if (gz_input)
    {
        gz_snp_file.open(base_file.c_str());
        if (!gz_snp_file.good())
        {
            std::string error_message = "Error: Cannot open base file: "
                                        + base_file + " (gz) to read!\n";
            throw std::runtime_error(error_message);
        }
        if (!is_index)
        {
            // if the input is index, we will keep the header, otherwise, we
            // will remove the header
            std::getline(gz_snp_file, line);
            message.append("GZ file detected. Header of file is:\n");
            message.append(line + "\n\n");
        }
        else
        {
            message.append("GZ file detected.");
        }
        reporter.report("Due to library restrictions, we cannot display "
                        "progress bar for gz");
    }
    else
    {
        snp_file.open(base_file.c_str());
        if (!snp_file.is_open())
        {
            std::string error_message =
                "Error: Cannot open base file: " + base_file;
            throw std::runtime_error(error_message);
        }
        snp_file.seekg(0, snp_file.end);
        file_length = snp_file.tellg();
        snp_file.clear();
        snp_file.seekg(0, snp_file.beg);
        // if the input is index, we will keep the header, otherwise, we
        // will remove the header
        if (!is_index) std::getline(snp_file, line);
    }
    double prev_progress = 0.0;

    std::unordered_set<std::string> dup_index;
    while ((gz_input && std::getline(gz_snp_file, line))
           || (!gz_input && std::getline(snp_file, line)))
    {
        if (!gz_input)
        {
            double progress = static_cast<double>(snp_file.tellg())
                              / static_cast<double>(file_length) * 100;
            if (progress - prev_progress > 0.01)
            {
                fprintf(stderr, "\rReading %03.2f%%", progress);
                prev_progress = progress;
            }
        }
        misc::trim(line);
        if (line.empty()) continue;
        num_line_in_base++;
        exclude = false;
        token = misc::split(line);
        if (token.size() <= max_index)
        {
            std::string error_message = line;
            error_message.append("\nMore index than column in data\n");
            throw std::runtime_error(error_message);
        }

        rs_id = token[col_index[+BASE_INDEX::RS]];
        if (dup_index.find(rs_id) != dup_index.end())
        {
            ++num_duplicated;
            continue;
        }
        // if this is not a duplicated SNP
        auto&& selection = m_snp_selection_list.find(rs_id);
        if ((!m_exclude_snp && selection == m_snp_selection_list.end())
            || (m_exclude_snp && selection != m_snp_selection_list.end()))
        {
            num_selected++;
            continue;
        }
        dup_index.insert(rs_id);
        chr_code = -1;
        if (has_col[+BASE_INDEX::CHR])
        {
            chr_code =
                get_chrom_code_raw(token[col_index[+BASE_INDEX::CHR]].c_str());
            // only two check, if it is haploid, or if it is sex chromosome
            if (chr_code < 0)
            {
                // invalid chromosome
                ++num_chr_filter;
                continue;
            }
            if (chr_code >= MAX_POSSIBLE_CHROM)
            {
                // sex chromosome
                num_haploid++;
                continue;
            }
            if (is_set(m_haploid_mask.data(), static_cast<uint32_t>(chr_code)))
            {
                // this is a haploid chromosome outside of the sex chromosome
                ++num_haploid;
                continue;
            }
        }
        if (chr_code == -1) { chr = ~size_t(0); }
        else
        {
            chr = static_cast<size_t>(chr_code);
        }
        ref_allele = (has_col[+BASE_INDEX::REF])
                         ? token[col_index[+BASE_INDEX::REF]]
                         : "";
        alt_allele = (has_col[+BASE_INDEX::ALT])
                         ? token[col_index[+BASE_INDEX::ALT]]
                         : "";
        std::transform(ref_allele.begin(), ref_allele.end(), ref_allele.begin(),
                       ::toupper);
        std::transform(alt_allele.begin(), alt_allele.end(), alt_allele.begin(),
                       ::toupper);
        if (has_col[+BASE_INDEX::BP])
        {
            // obtain the SNP coordinate
            try
            {
                loc = misc::string_to_size_t(
                    token[col_index[+BASE_INDEX::BP]].c_str());
            }
            catch (...)
            {
                std::string error_message =
                    "Error: Invalid loci for " + rs_id + ": "
                    + token[col_index[+BASE_INDEX::BP]] + "\n";
                throw std::runtime_error(error_message);
            }
        }
        else
        {
            loc = ~size_t(0);
        }
        to_remove = Genotype::within_region(exclusion_regions, chr, loc);
        if (to_remove && !exclude)
        {
            num_region_exclude++;
            continue;
        }
        maf = 1;
        // note, this is for reference
        maf_filtered = false;
        if (maf_control_filter && has_col[+BASE_INDEX::MAF])
        {
            // only read in if we want to perform MAF filtering
            try
            {
                maf = misc::convert<double>(
                    token[col_index[+BASE_INDEX::MAF]].c_str());
            }
            catch (...)
            {
                // exclude because we can't read the MAF, therefore assume
                // this is problematic
                ++num_maf_filter;
                continue;
            }
            if (maf < maf_control)
            {
                ++num_maf_filter;
                continue;
            }
        }
        maf_case_temp = 1;
        if (maf_case_filter && has_col[+BASE_INDEX::MAF_CASE])
        {
            // only read in the case MAf if we want to perform filtering on
            // it
            try
            {
                maf_case_temp = misc::convert<double>(
                    token[col_index[+BASE_INDEX::MAF_CASE]].c_str());
            }
            catch (...)
            {
                // again, any problem in parsing the MAF will lead to SNP
                // being deleted
                // we don't want to double count the MAF filtering, thus we
                // only add one to the maf filter count if we haven't
                // already filtered this SNP based on the control MAf
                ++num_maf_filter;
                continue;
            }
            if (maf_case_temp < maf_case)
            {
                ++num_maf_filter;
                continue;
            }
        }
        info_score = 1;
        if (info_filter && has_col[+BASE_INDEX::INFO])
        {
            // obtain the INFO score
            try
            {
                info_score = misc::convert<double>(
                    token[col_index[+BASE_INDEX::INFO]].c_str());
            }
            catch (...)
            {
                // if no info score, just assume it doesn't pass the QC
                ++num_info_filter;
                continue;
            }
            if (info_score < info_threshold)
            {
                ++num_info_filter;
                continue;
            }
        }
        pvalue = 2.0;
        try
        {
            pvalue = misc::convert<double>(token[col_index[+BASE_INDEX::P]]);
            if (pvalue < 0.0 || pvalue > 1.0)
            {
                std::string error_message =
                    "Error: Invalid p-value for " + rs_id + ": "
                    + token[col_index[+BASE_INDEX::P]] + "!\n";
                throw std::runtime_error(error_message);
            }
            else if (pvalue > max_threshold)
            {
                ++num_excluded;
                continue;
            }
        }
        catch (...)
        {
            ++num_not_converted;
            continue;
        }
        stat = 0.0;
        try
        {
            stat = misc::convert<double>(token[col_index[+BASE_INDEX::STAT]]);
            if (stat < 0 && !is_beta)
            {
                ++num_negative_stat;
                continue;
            }
            else if (misc::logically_equal(stat, 0.0) && !is_beta)
            {
                // we can't calculate log(0), so we will say we can't
                // convert it
                ++num_not_converted;
                continue;
            }
            else if (!is_beta)
                // if it is OR, we will perform natural log to convert it to
                // beta
                stat = log(stat);
        }
        catch (...)
        {
            // non-numeric statistic
            ++num_not_converted;
            continue;
        }
        if (!alt_allele.empty() && ambiguous(ref_allele, alt_allele))
        {
            // only coun the number if the snp is kept
            if (!exclude && keep_ambig) num_ambiguous++;
            if (!exclude) exclude = !keep_ambig;
        }
        if (!exclude)
        {
            category = -1;
            pthres = 0.0;
            if (fastscore)
            { category = calculate_category(pvalue, barlevels, pthres); }
            else
            {
                // calculate the threshold instead
                category = calculate_category(pvalue, bound_start, bound_inter,
                                              bound_end, pthres, no_full);
            }
            // now add SNP
            m_existed_snps.emplace_back(SNP(rs_id, chr, loc, ref_allele,
                                            alt_allele, stat, pvalue, category,
                                            pthres));
            m_existed_snps_index[rs_id] = m_existed_snps.size() - 1;
        }
    }
    if (gz_input)
        gz_snp_file.close();
    else
        snp_file.close();

    fprintf(stderr, "\rReading %03.2f%%\n", 100.0);
    message.append(std::to_string(num_line_in_base)
                   + " variant(s) observed in base file, with:\n");
    // should emit an error and terminate PRSice. But don't bother
    // at the moment. Can implement the kill switch later
    if (num_duplicated)
    {
        message.append(std::to_string(num_duplicated)
                       + " duplicated variant(s)\n");
    }
    if (num_selected)
    {
        message.append(std::to_string(num_selected)
                       + " variant(s) excluded based on user input\n");
    }
    if (num_chr_filter)
    {
        message.append(
            std::to_string(num_chr_filter)
            + " variant(s) excluded as they are on unknown/sex chromosome\n");
    }
    if (num_haploid)
    {
        message.append(std::to_string(num_haploid)
                       + " variant(s) located on haploid chromosome\n");
    }
    if (num_region_exclude)
    {
        message.append(
            std::to_string(num_region_exclude)
            + " variant(s) excluded as they fall within x-range region(s)\n");
    }
    if (num_maf_filter)
    {
        message.append(std::to_string(num_maf_filter)
                       + " variant(s) excluded due to MAF threshold\n");
    }
    if (num_info_filter)
    {
        message.append(std::to_string(num_info_filter)
                       + " variant(s) with INFO score less than "
                       + std::to_string(info_threshold) + "\n");
    }
    if (num_excluded)
    {
        message.append(std::to_string(num_excluded)
                       + " variant(s) excluded due to p-value threshold\n");
    }
    if (num_not_converted)
    {
        message.append(std::to_string(num_not_converted)
                       + " NA stat/p-value observed\n");
    }
    if (num_negative_stat)
    {
        message.append(std::to_string(num_negative_stat)
                       + " negative statistic observed. Maybe you have "
                         "forgotten the --beta flag?\n");
    }
    if (num_ambiguous)
    {
        message.append(std::to_string(num_ambiguous) + " ambiguous variant(s)");
        if (!keep_ambig) { message.append(" excluded"); }
        message.append("\n");
    }
    message.append(std::to_string(m_existed_snps.size())
                   + " total variant(s) included from base file\n\n");
    reporter.report(message);
    if (m_existed_snps.size() == 0)
    { throw std::runtime_error("Error: No valid variant remaining"); }
}


std::vector<std::string>
Genotype::load_genotype_prefix(const std::string& file_name)
{
    std::vector<std::string> genotype_files;
    std::ifstream multi;
    multi.open(file_name.c_str());
    if (!multi.is_open())
    {
        throw std::runtime_error(
            std::string("Error: Cannot open list file: " + file_name));
    }
    std::string line;
    while (std::getline(multi, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        genotype_files.push_back(line);
    }
    multi.close();
    return genotype_files;
}
void Genotype::init_chr(int num_auto, bool no_x, bool no_y, bool no_xy,
                        bool no_mt)
{
    // this initialize haploid mask as the maximum possible number

    if (num_auto < 0)
    {
        num_auto = -num_auto;
        m_autosome_ct = static_cast<uint32_t>(num_auto);
        m_xymt_codes[X_OFFSET] = -1;
        m_xymt_codes[Y_OFFSET] = -1;
        m_xymt_codes[XY_OFFSET] = -1;
        m_xymt_codes[MT_OFFSET] = -1;
        m_max_code = static_cast<uint32_t>(num_auto);
        fill_all_bits((static_cast<uint32_t>(num_auto)) + 1,
                      m_haploid_mask.data());
    }
    else
    {
        uint32_t unsign_num_auto = static_cast<uint32_t>(num_auto);
        m_autosome_ct = unsign_num_auto;
        m_xymt_codes[X_OFFSET] = num_auto + 1;
        m_xymt_codes[Y_OFFSET] = num_auto + 2;
        m_xymt_codes[XY_OFFSET] = num_auto + 3;
        m_xymt_codes[MT_OFFSET] = num_auto + 4;
        set_bit(unsign_num_auto + 1, m_haploid_mask.data());
        set_bit(unsign_num_auto + 2, m_haploid_mask.data());
        if (no_x)
        {
            m_xymt_codes[X_OFFSET] = -1;
            clear_bit(unsign_num_auto + 1, m_haploid_mask.data());
        }
        if (no_y)
        {
            m_xymt_codes[Y_OFFSET] = -1;
            clear_bit(unsign_num_auto + 2, m_haploid_mask.data());
        }
        if (no_xy) { m_xymt_codes[XY_OFFSET] = -1; }
        if (no_mt) { m_xymt_codes[MT_OFFSET] = -1; }
        if (m_xymt_codes[MT_OFFSET] != -1) { m_max_code = unsign_num_auto + 4; }
        else if (m_xymt_codes[XY_OFFSET] != -1)
        {
            m_max_code = unsign_num_auto + 3;
        }
        else if (m_xymt_codes[Y_OFFSET] != -1)
        {
            m_max_code = unsign_num_auto + 2;
        }
        else if (m_xymt_codes[X_OFFSET] != -1)
        {
            m_max_code = unsign_num_auto + 1;
        }
        else
        {
            m_max_code = unsign_num_auto;
        }
    }
    /*
    fill_all_bits(m_autosome_ct + 1, m_chrom_mask.data());
    for (uint32_t xymt_idx = 0; xymt_idx < XYMT_OFFSET_CT; ++xymt_idx) {
        int32_t cur_code = m_xymt_codes[xymt_idx];
        if (cur_code != -1) {
            set_bit(m_xymt_codes[xymt_idx], m_chrom_mask.data());
        }
    }
    m_chrom_start.resize(m_max_code); // 1 extra for the info
    */
}

std::unordered_set<std::string>
Genotype::load_snp_list(const std::string& input, Reporter& reporter)
{
    std::ifstream in;
    // first, we read in the file
    in.open(input.c_str());
    if (!in.is_open())
    {
        std::string error_message =
            "Error: Cannot open extract / exclude file: " + input;
        throw std::runtime_error(error_message);
    }
    std::string line;
    // we will return the "result" variable
    std::unordered_set<std::string> result;
    std::string message;
    // allow the input to be slightly more flexible. This flexibility relies on
    // the header or the column name of the file
    std::getline(in, line);
    misc::trim(line);
    std::vector<std::string> token = misc::split(line);
    // rs_index store the location of the RS ID
    size_t rs_index = 0;
    if (token.size() != 1)
    {
        // if we have more than one column, we will try and see if it is some
        // format we can recognized
        bool has_snp_colname = false;
        for (auto&& name : token)
        {
            // we will go through the header and identify any of the followings
            // (case insensitive): 1) SNP 2) RS 3) RS_ID 4) RS.ID 5) RSID

            std::transform(name.begin(), name.end(), name.begin(), ::toupper);
            if (name == "SNP" || name == "RS" || name == "RS_ID"
                || name == "RS.ID" || name == "RSID" || name == "VARIANT.ID"
                || name == "VARIANT_ID")
            {
                /// we will assume this column to contain the SNP ID
                has_snp_colname = true;
                message = name + " assume to be column containing SNP ID";
                reporter.report(message);
                break;
            }
            // we will continue to iterate until we reaches the end or found a
            // column with the specific names
            ++rs_index;
        }
        if (!has_snp_colname)
        {
            // if we reaches to the end of the column name list and didn't find
            // the required header, we will resort to familiar file formats
            if (token.size() == 6)
            {
                // with 6 column, we will assume this to be a bim file, where
                // the SNP ID is at the second column
                message = "SNP extraction/exclusion list contains 6 columns, "
                          "will assume this is a bim file, with the "
                          "second column contains the SNP ID";
                reporter.report(message);
                // we set the rs index to 1 (use the second column as the RS ID)
                rs_index = 1;
            }
            else
            {
                // otherwise, we will assume the first column contain the SNP ID
                message = "SNP extraction/exclusion list contains "
                          + misc::to_string(token.size())
                          + " columns, "
                            "will assume first column contains the SNP ID";
                reporter.report(message);
                // we set the rs index to 0 (use the first column as the RS ID)
                rs_index = 0;
            }
        }
    }
    else
    {
        // with only one column, it is easy, just use that column as the SNP ID
        // rs_index will be 0
        message =
            "Only one column detected, will assume only SNP ID is provided";
        reporter.report(message);
    }
    in.clear();
    in.seekg(0, std::ios::beg);
    // now that we have identify the format of the file we can read in the file
    while (std::getline(in, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        // don't think the parsing . into chr:bp will be helpful in the context
        // of an extraction / exclusion list
        // as this is a set, we don't need to worry about duplicates
        result.insert(token[rs_index]);
    }
    return result;
}

std::unordered_set<std::string> Genotype::load_ref(const std::string& input,
                                                   const std::string& delim,
                                                   bool ignore_fid)
{
    std::ifstream in;
    in.open(input.c_str());
    if (!in.is_open())
    {
        std::string error_message =
            "Error: Cannot open keep / remove file: " + input;
        throw std::runtime_error(error_message);
    }
    // read in the file
    std::string line;
    // now go through the sample file. We require the FID (if any) and IID  must
    // be the first 1/2 column of the file
    std::unordered_set<std::string> result;
    while (std::getline(in, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (ignore_fid) { result.insert(token[0]); }
        else
        {
            if (token.size() < 2)
                throw std::runtime_error(
                    "Error: Require FID and IID for extraction. "
                    "You can ignore the FID by using the --ignore-fid flag");
            result.insert(token[0] + delim + token[1]);
        }
    }
    in.close();
    return result;
}

// return true if  we need to work on it
bool Genotype::chr_code_check(int32_t chr_code, bool& sex_error,
                              bool& chr_error, std::string& error_message)
{
    if (chr_code < 0)
    {
        // invalid chr code
        chr_error = true;
        error_message = "Error: Invalid chromosome number for SNP";
        return true;
    }
    if (chr_code > MAX_POSSIBLE_CHROM
        || is_set(m_haploid_mask.data(), static_cast<uint32_t>(chr_code)))
    {
        // this is sex / mt chromosome
        sex_error = true;
        error_message = "Warning: Currently not support "
                        "haploid chromosome and sex "
                        "chromosomes\n";
        return true;
    }
    return false;
}

void Genotype::load_samples(const std::string& keep_file,
                            const std::string& remove_file,
                            const std::string& delim, bool verbose,
                            Reporter& reporter)
{
    if (!remove_file.empty())
    { m_sample_selection_list = load_ref(remove_file, delim, m_ignore_fid); }
    else if (!keep_file.empty())
    {
        m_remove_sample = false;
        m_sample_selection_list = load_ref(keep_file, delim, m_ignore_fid);
    }
    if (!m_is_ref)
    {
        // m_sample_names = gen_sample_vector();
        m_sample_id = gen_sample_vector(delim);
    }
    else
    {
        // don't bother loading up the sample vector as it should
        // never be used for reference panel (except for the founder_info)
        gen_sample_vector(delim);
    }
    std::string message = misc::to_string(m_unfiltered_sample_ct) + " people ("
                          + misc::to_string(m_num_male) + " male(s), "
                          + misc::to_string(m_num_female)
                          + " female(s)) observed\n";
    message.append(misc::to_string(m_founder_ct) + " founder(s) included");
    if (verbose) reporter.report(message);
    m_sample_selection_list.clear();
}

void Genotype::calc_freqs_and_intermediate(
    const double& maf_threshold, const double& geno_threshold,
    const double& info_threshold, const bool maf_filter, const bool geno_filter,
    const bool info_filter, const bool hard_coded, bool verbose,
    Reporter& reporter, Genotype* target)
{
    std::string message = "";
    m_num_geno_filter = 0;
    m_num_maf_filter = 0;
    m_num_info_filter = 0;
    calc_freq_gen_inter(maf_threshold, geno_threshold, info_threshold,
                        maf_filter, geno_filter, info_filter, hard_coded,
                        target);
    if (m_num_geno_filter != 0)
    {
        message.append(
            std::to_string(m_num_geno_filter)
            + " variant(s) excluded based on genotype missingness threshold\n");
    }
    if (m_num_maf_filter != 0)
    {
        message.append(std::to_string(m_num_maf_filter)
                       + " variant(s) excluded based on MAF threshold\n");
    }
    if (m_num_info_filter != 0)
    {
        message.append(
            std::to_string(m_num_info_filter)
            + " variant(s) excluded based on INFO score threshold\n");
    }
    if (!m_is_ref)
    { message.append(std::to_string(m_marker_ct) + " variant(s) included"); }
    else
    {
        message.append(std::to_string(target->m_existed_snps.size())
                       + " variant(s) remained");
    }

    if (verbose) reporter.report(message);
}

void Genotype::load_snps(
    const std::string& out,
    const std::vector<IITree<size_t, size_t>>& exclusion_regions, bool verbose,
    Reporter& reporter, Genotype* target)
{
    m_base_missed = 0;
    m_num_ambig = 0;
    m_num_xrange = 0;
    gen_snp_vector(exclusion_regions, out, target);
    m_marker_ct = m_existed_snps.size();
    std::string message = "";
    if (m_base_missed != 0)
    {
        message.append(std::to_string(m_base_missed)
                       + " variant(s) not found in previous data\n");
    }
    if (m_num_ambig != 0 && !m_keep_ambig)
    {
        message.append(std::to_string(m_num_ambig)
                       + " ambiguous variant(s) excluded\n");
    }
    else if (m_num_ambig != 0)
    {
        message.append(std::to_string(m_num_ambig)
                       + " ambiguous variant(s) kept\n");
    }
    if (m_num_xrange != 0)
    {
        message.append(std::to_string(m_num_xrange)
                       + " variant(s) removed as they fall within the "
                         "--x-range region(s)\n");
    }
    if (!m_is_ref)
    {
        message.append(std::to_string(m_marker_ct) + " variant(s) included\n");
    }
    else
    {
        message.append(std::to_string(target->m_existed_snps.size())
                       + " variant(s) remained\n");
    }

    if (verbose) reporter.report(message);
    m_snp_selection_list.clear();
    if (!m_is_ref)
    {
        if (m_marker_ct == 0)
        {
            message = "Error: No vairant remained!\n";
            throw std::runtime_error(message);
        }
    }
    else
    {
        if (target->m_existed_snps.size() == 0)
        {
            message =
                "Error: No vairant remained after matching with reference!\n";
            throw std::runtime_error(message);
        }
    }
}


Genotype::~Genotype() {}


void Genotype::efficient_clumping(Genotype& reference, Reporter& reporter)
{
    // the m_existed_snp must be sorted before coming into this equation
    reporter.report("Start performing clumping");
    // we want to initialize the vectors with size correspond to the sample in
    // the reference panel. Need to use the unfiltered sample size
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(reference.m_unfiltered_sample_ct);
    // again, we want to initialize the vector containing the founder
    // membership, which require us to know the number of founders in the
    // reference panel
    const uint32_t founder_ctv3 =
        BITCT_TO_ALIGNED_WORDCT(static_cast<uint32_t>(reference.m_founder_ct));
    const uint32_t founder_ctsplit = 3 * founder_ctv3;
    const uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(reference.m_founder_ct);
    const uintptr_t founder_ctv2 =
        QUATERCT_TO_ALIGNED_WORDCT(reference.m_founder_ct);

    // We only want to perform clumping if our R2 is higher than a minimum
    // threshold. Depending on whethre we do proxy clumping or not, the minimum
    // threshold can be the clump_p or clump_proxy parameter
    const double min_r2 =
        (m_use_proxy) ? std::min(m_clump_proxy, m_clump_r2) : m_clump_r2;
    // is_x is used in PLINK to indicate if the genotype is from the X
    // chromsome, as PRSice ignore any sex chromosome, we can set it as a
    // constant false
    const bool is_x = false;
    // when we read in the genotype, the genotype is stored in byte represented
    // by uintptr_t. Then we will use the PLINK 2 functions to obtain the R2
    // As this vector can be big, we will only initialize it once (memory
    // allocation can be expensive)
    std::vector<uintptr_t> genotype_vector(unfiltered_sample_ctl * 2);
    // We need a vector to indicate which SNPs are remaining after clumping so
    // that we can remove any clumped SNPs away. This is achieved by using this
    // boolean vector
    std::vector<bool> remain_core(m_existed_snps.size(), false);
    // next few parameters are used to store the intermediate output from PLINK
    // function, corresponding to the count of different genotypes
    double freq11;
    double freq11_expected;
    double freq1x;
    double freq2x;
    double freqx1;
    double freqx2;
    double dxx;
    // and this is the storage to result R2
    double r2 = -1.0;
    // The following two vectors are used for storing the intermediate output.
    // Again, put memory allocation at the beginning
    std::vector<uintptr_t> index_data(3 * founder_ctsplit + founder_ctv3);
    std::vector<uint32_t> index_tots(6);
    // pre-allocate the memory without bothering the memory pool stuff
    std::vector<uint32_t> ld_missing_count(
        static_cast<std::vector<uint32_t>::size_type>(m_max_window_size));
    // This is a data mask used by PLINK in the calculation of R2. Preallocate
    // to speed up
    std::vector<uintptr_t> founder_include2(founder_ctv2, 0);
    fill_quatervec_55(static_cast<uint32_t>(reference.m_founder_ct),
                      founder_include2.data());


// one way to speed things up as in PLINK 2 is to pre-allocate the memory space
// for what we need to do next. The following code did precisely that (borrow
// from PLINK2)
#ifdef __APPLE__
    int32_t mib[2];
    size_t sztmp;
#endif
    unsigned char* bigstack_ua = nullptr; // ua = unaligned
    unsigned char* bigstack_initial_base;
    int64_t llxx;
    intptr_t default_alloc_mb;
    intptr_t malloc_size_mb = 0;
#ifdef __APPLE__
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    llxx = 0;

    sztmp = sizeof(int64_t);
    sysctl(mib, 2, &llxx, &sztmp, nullptr, 0);
    llxx /= 1048576;
#else
#ifdef _WIN32
    MEMORYSTATUSEX memstatus;
    memstatus.dwLength = sizeof(memstatus);
    GlobalMemoryStatusEx(&memstatus);
    llxx = memstatus.ullTotalPhys / 1048576;
#else
    llxx = ((uint64_t) sysconf(_SC_PHYS_PAGES))
           * ((size_t) sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
    if (!llxx) { default_alloc_mb = BIGSTACK_DEFAULT_MB; }
    else if (llxx < (BIGSTACK_MIN_MB * 2))
    {
        default_alloc_mb = BIGSTACK_MIN_MB;
    }
    else
    {
        default_alloc_mb = llxx / 2;
    }
    if (!malloc_size_mb) { malloc_size_mb = default_alloc_mb; }
    else if (malloc_size_mb < BIGSTACK_MIN_MB)
    {
        malloc_size_mb = BIGSTACK_MIN_MB;
    }
    std::string message = "";
#ifndef __LP64__
    if (malloc_size_mb > 2047) { malloc_size_mb = 2047; }
#endif

    // m_max_window_size represent the maximum number of SNPs required for any
    // one window

    // calculate the maximum amount of memory required (the thing is, we can't
    // do the analysis if we don't have this amount as the largest region will
    // fail. So should I call this maximum amount of memory required or minimum
    // amount of memory required?
    malloc_size_mb = (static_cast<uintptr_t>(m_max_window_size) + 1)
                         * founder_ctv2 * sizeof(intptr_t) / 1048576
                     + 1;
    if (llxx)
    {
        // we have detected the memory, but need to check if that's enough
        if (malloc_size_mb > llxx)
        {
            std::string error_message =
                "Error: Insufficient memory for clumping! Require "
                + misc::to_string(malloc_size_mb) + " MB but detected only "
                + misc::to_string(llxx) + " MB";
            throw std::runtime_error(error_message);
        }
        else
        {
            message = misc::to_string(llxx) + " MB RAM detected; reserving "
                      + misc::to_string(malloc_size_mb) + " MB for clumping\n";
        }
    }
    else
    {
        message = "Failed to calculate system memory. Attemping to reserve"
                  + misc::to_string(malloc_size_mb) + " MB for clumping\n";
    }
    // now allocate the memory into a pointer
    bigstack_ua = reinterpret_cast<unsigned char*>(malloc(
        static_cast<uintptr_t>(malloc_size_mb) * 1048576 * sizeof(char)));
    // if fail, return nullptr which will then get into the while loop
    if (!bigstack_ua)
    { throw std::runtime_error("Failed to allocate required memory"); }
    else
    {
        message.append("Allocated " + misc::to_string(malloc_size_mb)
                       + " MB successfully");
    }
    reporter.report(message);
    // force 64-byte align to make cache line sensitivity work (from PLINK, not
    // familiar with computer programming to know this...)
    // will stay with the old style cast to avoid trouble
    bigstack_initial_base =
        (unsigned char*) round_up_pow2((uintptr_t) bigstack_ua, CACHELINE);

    // window data is the pointer walking through the allocated memory
    uintptr_t* window_data = nullptr;
    // now we can point the window pointer to the start of the allocated memory
    window_data = (uintptr_t*) bigstack_initial_base;
    uintptr_t* window_data_ptr = nullptr;
    unsigned char* g_bigstack_end = &(
        bigstack_initial_base[(static_cast<uintptr_t>(malloc_size_mb) * 1048576
                               - static_cast<uintptr_t>(bigstack_initial_base
                                                        - bigstack_ua))
                              & (~(CACHELINE - ONELU))]);

    uintptr_t num_core_snps = 0;
    // and max_window_size is the number of windows we can handle in one round
    // given the memory that we have. 1 window = 1 SNP
    uintptr_t max_window_size =
        (((uintptr_t) g_bigstack_end) - ((uintptr_t) bigstack_initial_base))
        / (founder_ctv2 * sizeof(intptr_t));

    g_bigstack_end = nullptr;
    // a counter to count how many windows have we read
    uintptr_t cur_window_size = 0;
    // if max_window size is 0, it means we don't have enough memory for
    // clumping
    if (!max_window_size)
    { throw std::runtime_error("Error: Not enough memory for clumping!"); }

    double prev_progress = -1.0;
    const auto num_snp = m_existed_snps.size();
    for (size_t i_snp = 0; i_snp < num_snp; ++i_snp)
    {
        // now start iterate through each SNP
        double progress =
            static_cast<double>(i_snp) / static_cast<double>(num_snp) * 100;
        if (progress - prev_progress > 0.01)
        {
            fprintf(stderr, "\rClumping Progress: %03.2f%%", progress);
            prev_progress = progress;
        }
        // get the index on target.m_existed_snp to the i_th most significant
        // SNP
        auto&& cur_snp_index = m_sort_by_p_index[i_snp];
        // read in the current SNP
        auto&& cur_target_snp = m_existed_snps[cur_snp_index];
        if (cur_target_snp.clumped() || cur_target_snp.p_value() > m_clump_p)
            // ignore any SNP that are clumped or that has a p-value higher than
            // the clump-p threshold
            continue;
        // with the new change, all SNPs are stored in target
        // so we can safely ignore finding it in the reference
        // though we still need to get the sample size information
        // from the reference panel
        // Any SNP with p-value less than clump-p will be ignored
        // because they can never be an index SNP and thus are not of our
        // interested

        // this is the first SNP we should read from
        std::vector<size_t>::size_type start =
            static_cast<std::vector<size_t>::size_type>(
                cur_target_snp.low_bound());
        // this is the first SNP we should ignore
        std::vector<size_t>::size_type end =
            static_cast<std::vector<size_t>::size_type>(
                cur_target_snp.up_bound());
        // reset our pointer to the start of the memory stack as we are working
        // on a new core SNP
        window_data_ptr = window_data;
        cur_window_size = 0;
        // transversing on TARGET
        // now we will transverse any SNP that comes before the index SNP in the
        // file
        for (size_t i_pair = start; i_pair < cur_snp_index; i_pair++)
        {
            // the start and end correspond to index on m_existed_snps instead
            // of m_sort_by_p, so we can skip reading from m_sort_by_p
            // the current SNP
            auto&& pair_target_snp = m_existed_snps[i_pair];
            if (pair_target_snp.clumped()
                || pair_target_snp.p_value() > m_clump_p)
                // ignore SNP that are clumped or that has higher p-value than
                // threshold
                continue;
            // Something PLINK does. I suspect this is to reset the content of
            // the pointer to 0
            window_data_ptr[founder_ctv2 - 2] = 0;
            window_data_ptr[founder_ctv2 - 1] = 0;
            if (++cur_window_size == max_window_size)
            { throw std::runtime_error("Error: Out of memory!"); }
            // read in the genotype data from the memory
            // this depends on the type of the reference.
            // Most important information is the ref_byte_pos (reference byte
            // position for reading) and ref_file_name (which reference file
            // should we read from)
            reference.read_genotype(window_data_ptr,
                                    pair_target_snp.ref_byte_pos(),
                                    pair_target_snp.ref_file_name());
            // we then move the pointer forward to the next space in the memory
            window_data_ptr = &(window_data_ptr[founder_ctv2]);
        }

        if (++cur_window_size == max_window_size)
        { throw std::runtime_error("Error: Out of memory!"); }
        // now we want to read in the index / core SNP
        // reset the content of the pointer again
        window_data_ptr[founder_ctv2 - 2] = 0;
        window_data_ptr[founder_ctv2 - 1] = 0;
        // reset the index_data information
        std::fill(index_data.begin(), index_data.end(), 0);
        // then we can read in the genotype from the reference panel
        // note the use of cur_target_snp
        reference.read_genotype(window_data_ptr, cur_target_snp.ref_byte_pos(),
                                cur_target_snp.ref_file_name());
        // generate the required data mask
        // Disclaimer: For the next few lines, they are from PLINK and I don't
        // fully understand what they are doing
        vec_datamask(reference.m_founder_ct, 0, window_data_ptr,
                     founder_include2.data(), index_data.data());
        // then populate the index_tots
        index_tots[0] = popcount2_longs(index_data.data(), founder_ctl2);
        vec_datamask(reference.m_founder_ct, 2, window_data_ptr,
                     founder_include2.data(), &(index_data[founder_ctv2]));
        index_tots[1] =
            popcount2_longs(&(index_data[founder_ctv2]), founder_ctl2);
        vec_datamask(reference.m_founder_ct, 3, window_data_ptr,
                     founder_include2.data(), &(index_data[2 * founder_ctv2]));
        index_tots[2] =
            popcount2_longs(&(index_data[2 * founder_ctv2]), founder_ctl2);
        // we have finished reading the index and stored the necessary
        // inforamtion, we can now calculate the R2 between the index and
        // previous SNPs
        // move back to the front of the memory (we don't need to read the SNP
        // again as the info is stored in the memory)
        window_data_ptr = window_data;
        for (size_t i_pair = start; i_pair < cur_snp_index; i_pair++)
        {
            auto&& pair_target_snp = m_existed_snps[i_pair];
            if (pair_target_snp.clumped()
                || pair_target_snp.p_value() > m_clump_p)
                // Again, ignore unwanted SNP
                continue;
            r2 = -1;
            uint32_t counts[18];
            // calculate the counts
            // these counts are then used for calculation of R2. However, I
            // don't fully understand the algorithm here (copy from PLINK2)
            genovec_3freq(window_data_ptr, index_data.data(), founder_ctl2,
                          &(counts[0]), &(counts[1]), &(counts[2]));
            counts[0] = index_tots[0] - counts[0] - counts[1] - counts[2];
            genovec_3freq(window_data_ptr, &(index_data[founder_ctv2]),
                          founder_ctl2, &(counts[3]), &(counts[4]),
                          &(counts[5]));
            counts[3] = index_tots[1] - counts[3] - counts[4] - counts[5];
            genovec_3freq(window_data_ptr, &(index_data[2 * founder_ctv2]),
                          founder_ctl2, &(counts[6]), &(counts[7]),
                          &(counts[8]));
            counts[6] = index_tots[2] - counts[6] - counts[7] - counts[8];
            if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x,
                                        &freqx1, &freqx2, &freq11))
            {
                // if the calculation is sucessful, we can then calculate the R2
                freq11_expected = freqx1 * freq1x;
                dxx = freq11 - freq11_expected;
                // message from PLINK:
                // if r^2 threshold is 0, let everything else through but
                // exclude the apparent zeroes.  Zeroes *are* included if
                // r2_thresh is negative,
                // though (only nans are rejected then).
                if (fabs(dxx) < SMALL_EPSILON
                    || fabs(freq11_expected * freq2x * freqx2) < SMALL_EPSILON)
                { r2 = 0.0; }
                else
                {
                    r2 = dxx * dxx / (freq11_expected * freq2x * freqx2);
                }
            }
            if (r2 >= min_r2)
            {
                // if the R2 between two SNP is higher than the minim threshold,
                // we will perform clumping
                // use the core SNP to clump the pair_target_snp
                cur_target_snp.clump(pair_target_snp, r2, m_use_proxy,
                                     m_clump_proxy);
            }
            // travel to the next snp
            window_data_ptr = &(window_data_ptr[founder_ctv2]);
        }
        // now we can read the SNPs that come after the index SNP in the file

        for (size_t i_pair = cur_snp_index + 1; i_pair < end; ++i_pair)
        {
            // we don't need to store the SNP information, as we can process
            // each SNP immediately. Always reset the pointer to the beginning
            // of the stack
            window_data_ptr = window_data;
            // read in the SNP information from teh target
            auto&& pair_target_snp = m_existed_snps[i_pair];
            if (pair_target_snp.clumped()
                || pair_target_snp.p_value() > m_clump_p)
                // skip if not required
                continue;
            // reset data
            window_data_ptr[founder_ctv2 - 2] = 0;
            window_data_ptr[founder_ctv2 - 1] = 0;
            // read in the genotype information
            reference.read_genotype(window_data_ptr,
                                    pair_target_snp.ref_byte_pos(),
                                    pair_target_snp.ref_file_name());

            r2 = -1;
            // obtain count using PLINK magic
            uint32_t counts[18];
            genovec_3freq(window_data_ptr, index_data.data(), founder_ctl2,
                          &(counts[0]), &(counts[1]), &(counts[2]));
            counts[0] = index_tots[0] - counts[0] - counts[1] - counts[2];
            genovec_3freq(window_data_ptr, &(index_data[founder_ctv2]),
                          founder_ctl2, &(counts[3]), &(counts[4]),
                          &(counts[5]));
            counts[3] = index_tots[1] - counts[3] - counts[4] - counts[5];
            genovec_3freq(window_data_ptr, &(index_data[2 * founder_ctv2]),
                          founder_ctl2, &(counts[6]), &(counts[7]),
                          &(counts[8]));
            counts[6] = index_tots[2] - counts[6] - counts[7] - counts[8];
            if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x,
                                        &freqx1, &freqx2, &freq11))
            {
                freq11_expected = freqx1 * freq1x;
                dxx = freq11 - freq11_expected;
                // Message in PLINK:
                // if r^2 threshold is 0, let everything else through but
                // exclude the apparent zeroes.  Zeroes *are* included if
                // r2_thresh is negative,
                // though (only nans are rejected then).
                if (fabs(dxx) < SMALL_EPSILON
                    || fabs(freq11_expected * freq2x * freqx2) < SMALL_EPSILON)
                { r2 = 0.0; }
                else
                {
                    r2 = dxx * dxx / (freq11_expected * freq2x * freqx2);
                }
            }
            // now perform clumping if required
            if (r2 >= min_r2)
            {
                cur_target_snp.clump(pair_target_snp, r2, m_use_proxy,
                                     m_clump_proxy);
            }
        }
        // we set the core SNP to be "clumped" so that it will no longer be
        // considered by other SNP
        cur_target_snp.set_clumped();
        // we set the remain_core to true so that we will keep it at the end
        remain_core[cur_snp_index] = true;
        num_core_snps++;
        // we also get the p-value threshold of the SNP so that we know what
        // thresholds are included in our analysis. This information is vital
        // for the generation of all score file
    }
    fprintf(stderr, "\rClumping Progress: %03.2f%%\n\n", 100.0);
    // now we release the memory stack
    window_data = nullptr;
    window_data_ptr = nullptr;
    free(bigstack_ua);
    bigstack_ua = nullptr;
    bigstack_initial_base = nullptr;
    if (num_core_snps != m_existed_snps.size())
    {
        //  if there are clumped SNPs, we would like to remove them from the
        //  vector
        // the remain_core SNP has the same order as the m_existed_snps
        // (core_snp_index is the index correspond to the m_existed_snp of
        // target) so it is ok for us to simply perform the remove if

        // because of the algorithm, any SNP with p-value above the clump-p
        // threshold will never be process, thus the remain_core for those SNPs
        // will always be 0, and will be removed (though they will also
        // misleadingly be considered as "clumped")
        // TODO: Issue a warning if the highest p-value threshold is higher than
        // the clump-p threshold
        m_existed_snps.erase(
            std::remove_if(
                m_existed_snps.begin(), m_existed_snps.end(),
                [&remain_core, this](const SNP& s) {
                    return !remain_core[&s - &*begin(m_existed_snps)];
                }),
            m_existed_snps.end());
        m_existed_snps.shrink_to_fit();
    }
    // we no longer require the index. might as well clear it (and hope it will
    // release the memory)
    m_existed_snps_index.clear();

    message = "Number of variant(s) after clumping : "
              + misc::to_string(m_existed_snps.size());
    reporter.report(message);
}


bool Genotype::prepare_prsice()
{
    if (m_existed_snps.size() == 0) return false;
    std::sort(begin(m_existed_snps), end(m_existed_snps),
              [](SNP const& t1, SNP const& t2) {
                  if (t1.category() == t2.category())
                  {
                      if (t1.file_name() == t2.file_name())
                      { return t1.byte_pos() < t2.byte_pos(); }
                      else
                          return t1.file_name() < t2.file_name();
                  }
                  else
                      return t1.category() < t2.category();
              });
    return true;
}
void Genotype::build_membership_matrix(
    std::vector<size_t>& region_membership,
    std::vector<size_t>& region_start_idx, const size_t num_sets,
    const std::string& out, const std::vector<std::string>& region_name,
    const bool print_snps)
{
    std::vector<std::vector<size_t>> temporary_storage(num_sets);
    std::vector<size_t> idx;
    std::unordered_set<double> threshold;
    std::ofstream snp_out;
    const std::string snp_name = out + ".snp";
    const bool is_prset = (num_sets != 2);
    if (print_snps)
    {
        snp_out.open(snp_name.c_str());
        if (!snp_out.is_open())
        {
            std::string error_message =
                "Error: Cannot open file: " + snp_name + " to write!\n";
            throw std::runtime_error(error_message);
        }
        snp_out << "CHR\tSNP\tBP\tP";
        for (size_t i = 0; i < region_name.size() - !is_prset; ++i)
        { snp_out << "\t" << region_name[i]; }
        snp_out << "\n";
    }
    size_t prev_idx;
    bool has_snp = false;
    if (is_prset)
    {
        for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp)
        {
            auto&& snp = m_existed_snps[i_snp];
            prev_idx = 0;
            if (threshold.find(snp.get_threshold()) == threshold.end())
            {
                ++m_num_thresholds;
                threshold.insert(snp.get_threshold());
                m_thresholds.push_back(snp.get_threshold());
            }
            if (print_snps)
            {
                snp_out << snp.chr() << "\t" << snp.rs() << "\t" << snp.loc()
                        << "\t" << snp.p_value();
                idx = snp.get_set_idx(num_sets);
                for (auto&& index : idx)
                {
                    assert(index >= prev_idx);
                    for (; prev_idx < index; ++prev_idx) { snp_out << "\t0"; }
                    snp_out << "\t1";
                    if (index > 1) has_snp = true;
                    temporary_storage[index].push_back(i_snp);
                    prev_idx = index + 1;
                }
                for (; prev_idx < num_sets; ++prev_idx) { snp_out << "\t0"; }
                snp_out << "\n";
            }
            else
            {
                idx = snp.get_set_idx(num_sets);
                for (auto&& index : idx)
                {
                    temporary_storage[index].push_back(i_snp);
                    if (index > 1) has_snp = true;
                }
            }
        }
        size_t cur_idx = 0;
        for (size_t i = 0; i < num_sets; ++i)
        {
            region_start_idx.push_back(cur_idx);
            if (!temporary_storage[i].empty())
            {
                region_membership.insert(region_membership.end(),
                                         temporary_storage[i].begin(),
                                         temporary_storage[i].end());
                cur_idx = region_membership.size();
            }
        }
        if (!has_snp)
        {
            std::string error_message =
                "Error: None of the gene sets contain any SNP(s) after "
                "clumping. Have you provided the correct input? E.g. GMT file "
                "containing Entrez ID with GTF files that uses the Ensembl "
                "gene ID?\n";
            throw std::runtime_error(error_message);
        }
    }
    else
    {
        // not PRSet. We know the membership is 1:num_snp
        region_start_idx.push_back(0);
        region_start_idx.push_back(m_existed_snps.size());
        if (print_snps)
        {
            for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp)
            {
                auto&& snp = m_existed_snps[i_snp];
                if (threshold.find(snp.get_threshold()) == threshold.end())
                {
                    threshold.insert(snp.get_threshold());
                    m_num_thresholds++;
                    m_thresholds.push_back(snp.get_threshold());
                }
                snp_out << snp.chr() << "\t" << snp.rs() << "\t" << snp.loc()
                        << "\t" << snp.p_value() << "\t1\n";
                region_membership.push_back(i_snp);
            }
        }
        else
        {
            // directly initialize the vector
            region_membership.resize(m_existed_snps.size());
            for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp)
            {
                auto&& snp = m_existed_snps[i_snp];
                if (threshold.find(snp.get_threshold()) == threshold.end())
                {
                    threshold.insert(snp.get_threshold());
                    m_num_thresholds++;
                    m_thresholds.push_back(snp.get_threshold());
                }
                region_membership[i_snp] = i_snp;
            }
        }
    }
}

void Genotype::get_null_score(const size_t& set_size, const size_t& prev_size,
                              std::vector<size_t>& background_list,
                              const bool first_run,
                              const bool require_statistic,
                              const bool use_ref_maf)
{

    if (m_existed_snps.empty() || set_size >= m_existed_snps.size()) return;
    // we will initailize a selected_snp_index containing the index of SNPs that
    // we'd like to add / assign to our PRS in the current round.
    // we will get anything from (prev_size , set_size]
    assert(prev_size < set_size);
    std::vector<size_t>::iterator select_start = background_list.begin();
    std::advance(select_start, prev_size);
    std::vector<size_t>::iterator select_end = background_list.begin();
    std::advance(select_end, set_size);
    std::sort(select_start, select_end);
    read_score(select_start, select_end, first_run, use_ref_maf);
    if (require_statistic)
    {
        misc::RunningStat rs;
        size_t num_prs = m_prs_info.size();
        for (size_t i = 0; i < num_prs; ++i)
        {
            if (!IS_SET(m_sample_include, i)) continue;
            if (m_prs_info[i].num_snp == 0) { rs.push(0.0); }
            else
            {
                rs.push(m_prs_info[i].prs
                        / static_cast<double>(m_prs_info[i].num_snp));
            }
        }
        m_mean_score = rs.mean();
        m_score_sd = rs.sd();
    }
}

bool Genotype::get_score(std::vector<size_t>::const_iterator& start_index,
                         const std::vector<size_t>::const_iterator& end_index,
                         double& cur_threshold, uint32_t& num_snp_included,
                         const bool non_cumulate, const bool require_statistic,
                         const bool first_run, const bool use_ref_maf)
{
    // if there are no SNPs or we are at the end
    if (m_existed_snps.size() == 0 || start_index == end_index
        || (*start_index) == m_existed_snps.size())
        return false;
    // reset number of SNPs if we don't need cumulative PRS
    if (non_cumulate) num_snp_included = 0;
    int cur_category = m_existed_snps[(*start_index)].category();
    cur_threshold = m_existed_snps[(*start_index)].get_threshold();
    std::vector<size_t>::const_iterator region_end = start_index;
    for (; region_end != end_index; ++region_end)
    {
        if (m_existed_snps[(*region_end)].category() != cur_category)
        {
            cur_category = m_existed_snps[(*region_end)].category();
            break;
        }
        num_snp_included++;
    }
    read_score(start_index, region_end, (non_cumulate || first_run),
               use_ref_maf);
    // update the current index
    start_index = region_end;
    // if ((*start_index) == 0) return -1;
    if (require_statistic)
    {
        misc::RunningStat rs;
        size_t num_prs = m_prs_info.size();
        for (size_t i = 0; i < num_prs; ++i)
        {
            if (!IS_SET(m_sample_include, i)) continue;
            if (m_prs_info[i].num_snp == 0) { rs.push(0.0); }
            else
            {
                rs.push(m_prs_info[i].prs
                        / static_cast<double>(m_prs_info[i].num_snp));
            }
        }
        m_mean_score = rs.mean();
        m_score_sd = rs.sd();
    }
    return true;
}

/**
 * DON'T TOUCH AREA
 *
 */

double Genotype::calc_lnlike(double known11, double known12, double known21,
                             double known22, double center_ct_d, double freq11,
                             double freq12, double freq21, double freq22,
                             double half_hethet_share, double freq11_incr)
{
    double lnlike;
    freq11 += freq11_incr;
    freq22 += freq11_incr;
    freq12 += half_hethet_share - freq11_incr;
    freq21 += half_hethet_share - freq11_incr;
    lnlike = center_ct_d * log(freq11 * freq22 + freq12 * freq21);
    if (known11 != 0.0) { lnlike += known11 * log(freq11); }
    if (known12 != 0.0) { lnlike += known12 * log(freq12); }
    if (known21 != 0.0) { lnlike += known21 * log(freq21); }
    if (known22 != 0.0) { lnlike += known22 * log(freq22); }
    return lnlike;
}

// This is where the magic happens
uint32_t Genotype::em_phase_hethet(double known11, double known12,
                                   double known21, double known22,
                                   uint32_t center_ct, double* freq1x_ptr,
                                   double* freq2x_ptr, double* freqx1_ptr,
                                   double* freqx2_ptr, double* freq11_ptr,
                                   uint32_t* onside_sol_ct_ptr)
{
    // Returns 1 if at least one SNP is monomorphic over all valid
    // observations; returns 0 otherwise, and fills all frequencies
    // using the maximum likelihood solution to the cubic equation.
    // (We're discontinuing most use of EM phasing since better
    // algorithms have been developed, but the two marker case is
    // mathematically clean and fast enough that it'll probably remain
    // useful as an input for some of those better algorithms...)
    double center_ct_d = static_cast<int32_t>(center_ct);
    double twice_tot = known11 + known12 + known21 + known22 + 2 * center_ct_d;
    uint32_t sol_start_idx = 0;
    uint32_t sol_end_idx = 1;
    double solutions[3];
    double twice_tot_recip;
    double half_hethet_share;
    double freq11;
    double freq12;
    double freq21;
    double freq22;
    double prod_1122;
    double prod_1221;
    double incr_1122;
    double best_sol;
    double best_lnlike;
    double cur_lnlike;
    double freq1x;
    double freq2x;
    double freqx1;
    double freqx2;
    double lbound;
    double dxx;
    uint32_t cur_sol_idx;
    // shouldn't have to worry about subtractive cancellation problems
    // here
    if (twice_tot == 0.0) { return 1; }
    twice_tot_recip = 1.0 / twice_tot;
    freq11 = known11 * twice_tot_recip;
    freq12 = known12 * twice_tot_recip;
    freq21 = known21 * twice_tot_recip;
    freq22 = known22 * twice_tot_recip;
    prod_1122 = freq11 * freq22;
    prod_1221 = freq12 * freq21;
    half_hethet_share = center_ct_d * twice_tot_recip;
    // the following four values should all be guaranteed nonzero except
    // in the NAN case
    freq1x = freq11 + freq12 + half_hethet_share;
    freq2x = 1.0 - freq1x;
    freqx1 = freq11 + freq21 + half_hethet_share;
    freqx2 = 1.0 - freqx1;
    if (center_ct)
    {
        if ((prod_1122 != 0.0) || (prod_1221 != 0.0))
        {
            sol_end_idx = cubic_real_roots(
                0.5
                    * (freq11 + freq22 - freq12 - freq21
                       - 3 * half_hethet_share),
                0.5
                    * (prod_1122 + prod_1221
                       + half_hethet_share
                             * (freq12 + freq21 - freq11 - freq22
                                + half_hethet_share)),
                -0.5 * half_hethet_share * prod_1122, solutions);
            while (sol_end_idx
                   && (solutions[sol_end_idx - 1]
                       > half_hethet_share + SMALLISH_EPSILON))
            {
                sol_end_idx--;
                // assert(sol_end_idx && sol_end_idx-1 >= 0);
            }
            while ((sol_start_idx < sol_end_idx)
                   && (solutions[sol_start_idx] < -SMALLISH_EPSILON))
            {
                sol_start_idx++;
                // assert((sol_start_idx < sol_end_idx) &&sol_start_idx
                // < 3);
            }
            if (sol_start_idx == sol_end_idx)
            {
                // Lost a planet Master Obi-Wan has.  How
                // embarrassing... lost root must be a double root at
                // one of the boundary points, just check their
                // likelihoods
                sol_start_idx = 0;
                sol_end_idx = 2;
                solutions[0] = 0;
                solutions[1] = half_hethet_share;
            }
            else
            {
                if (solutions[sol_start_idx] < 0)
                { solutions[sol_start_idx] = 0; } // checking here
                if (solutions[sol_end_idx - 1] > half_hethet_share)
                { solutions[sol_end_idx - 1] = half_hethet_share; }
            }
        }
        else
        {
            solutions[0] = 0;
            if ((freq22 + SMALLISH_EPSILON < half_hethet_share + freq21)
                && (freq21 + SMALLISH_EPSILON < half_hethet_share + freq22))
            {
                sol_end_idx = 3;
                solutions[1] = (half_hethet_share + freq21 - freq22) * 0.5;
                solutions[2] = half_hethet_share;
            }
            else
            {
                sol_end_idx = 2;
                solutions[1] = half_hethet_share;
            }
        }
        best_sol = solutions[sol_start_idx];
        if (sol_end_idx > sol_start_idx + 1)
        {
            // select largest log likelihood
            best_lnlike = calc_lnlike(known11, known12, known21, known22,
                                      center_ct_d, freq11, freq12, freq21,
                                      freq22, half_hethet_share, best_sol);
            cur_sol_idx = sol_start_idx + 1;
            do
            {
                incr_1122 = solutions[cur_sol_idx];
                cur_lnlike = calc_lnlike(known11, known12, known21, known22,
                                         center_ct_d, freq11, freq12, freq21,
                                         freq22, half_hethet_share, incr_1122);
                if (cur_lnlike > best_lnlike)
                {
                    cur_lnlike = best_lnlike;
                    best_sol = incr_1122;
                }
            } while (++cur_sol_idx < sol_end_idx);
        }
        if (onside_sol_ct_ptr && (sol_end_idx > sol_start_idx + 1))
        {
            if (freqx1 * freq1x >= freq11)
            {
                dxx = freq1x * freqx1 - freq11;
                if (dxx > half_hethet_share) { dxx = half_hethet_share; }
            }
            else
            {
                dxx = 0.0;
            }
            // okay to NOT count suboptimal boundary points because they
            // don't permit direction changes within the main interval
            // this should exactly match haploview_blocks_classify()'s D
            // sign check
            if ((freq11 + best_sol) - freqx1 * freq1x >= 0.0)
            {
                if (best_sol > dxx + SMALLISH_EPSILON)
                { lbound = dxx + SMALLISH_EPSILON; }
                else
                {
                    lbound = dxx;
                }
                if (best_sol < half_hethet_share - SMALLISH_EPSILON)
                { half_hethet_share -= SMALLISH_EPSILON; }
            }
            else
            {
                if (best_sol > SMALLISH_EPSILON) { lbound = SMALLISH_EPSILON; }
                else
                {
                    lbound = 0.0;
                }
                if (best_sol < dxx - SMALLISH_EPSILON)
                { half_hethet_share = dxx - SMALLISH_EPSILON; }
                else
                {
                    half_hethet_share = dxx;
                }
            }
            for (cur_sol_idx = sol_start_idx; cur_sol_idx < sol_end_idx;
                 cur_sol_idx++)
            {
                if (solutions[cur_sol_idx] < lbound) { sol_start_idx++; }
                if (solutions[cur_sol_idx] > half_hethet_share) { break; }
            }
            if (cur_sol_idx >= sol_start_idx + 2)
            { *onside_sol_ct_ptr = cur_sol_idx - sol_start_idx; }
        }
        freq11 += best_sol;
    }
    else if ((prod_1122 == 0.0) && (prod_1221 == 0.0))
    {
        return 1;
    }
    *freq1x_ptr = freq1x;
    *freq2x_ptr = freq2x;
    *freqx1_ptr = freqx1;
    *freqx2_ptr = freqx2;
    *freq11_ptr = freq11;
    return 0;
}

uint32_t Genotype::em_phase_hethet_nobase(uint32_t* counts, uint32_t is_x1,
                                          uint32_t is_x2, double* freq1x_ptr,
                                          double* freq2x_ptr,
                                          double* freqx1_ptr,
                                          double* freqx2_ptr,
                                          double* freq11_ptr)
{
    // if is_x1 and/or is_x2 is set, counts[9]..[17] are male-only
    // counts.
    double known11 = (double) (2 * counts[0] + counts[1] + counts[3]);
    double known12 = (double) (2 * counts[2] + counts[1] + counts[5]);
    double known21 = (double) (2 * counts[6] + counts[3] + counts[7]);
    double known22 = (double) (2 * counts[8] + counts[5] + counts[7]);
    if (is_x1 || is_x2)
    {
        if (is_x1 && is_x2)
        {
            known11 -= (double) ((int32_t) counts[9]);
            known12 -= (double) ((int32_t) counts[11]);
            known21 -= (double) ((int32_t) counts[15]);
            known22 -= (double) ((int32_t) counts[17]);
        }
        else if (is_x1)
        {
            known11 -=
                ((double) (2 * counts[9] + counts[10])) * (1.0 - SQRT_HALF);
            known12 -=
                ((double) (2 * counts[11] + counts[10])) * (1.0 - SQRT_HALF);
            known21 -=
                ((double) (2 * counts[15] + counts[16])) * (1.0 - SQRT_HALF);
            known22 -=
                ((double) (2 * counts[17] + counts[16])) * (1.0 - SQRT_HALF);
        }
        else
        {
            known11 -=
                ((double) (2 * counts[9] + counts[12])) * (1.0 - SQRT_HALF);
            known12 -=
                ((double) (2 * counts[11] + counts[12])) * (1.0 - SQRT_HALF);
            known21 -=
                ((double) (2 * counts[15] + counts[14])) * (1.0 - SQRT_HALF);
            known22 -=
                ((double) (2 * counts[17] + counts[14])) * (1.0 - SQRT_HALF);
        }
    }
    return em_phase_hethet(known11, known12, known21, known22, counts[4],
                           freq1x_ptr, freq2x_ptr, freqx1_ptr, freqx2_ptr,
                           freq11_ptr, nullptr);
}

uint32_t Genotype::load_and_split3(uintptr_t* rawbuf,
                                   uint32_t unfiltered_sample_ct,
                                   uintptr_t* casebuf, uint32_t case_ctv,
                                   uint32_t ctrl_ctv, uint32_t do_reverse,
                                   uint32_t is_case_only,
                                   uintptr_t* nm_info_ptr)
{
    uintptr_t* rawbuf_end = &(rawbuf[unfiltered_sample_ct / BITCT2]);
    uintptr_t* ctrlbuf = &(casebuf[3 * case_ctv]);
    uintptr_t case_words[4];
    uintptr_t ctrl_words[4];
    uint32_t case_rem = 0;
    uint32_t ctrl_rem = 0;
    uint32_t read_shift_max = BITCT2;
    uint32_t sample_uidx = 0;
    uint32_t offset0_case = do_reverse * 2 * case_ctv;
    uint32_t offset2_case = (1 - do_reverse) * 2 * case_ctv;
    uint32_t offset0_ctrl = do_reverse * 2 * ctrl_ctv;
    uint32_t offset2_ctrl = (1 - do_reverse) * 2 * ctrl_ctv;
    uint32_t read_shift;
    uintptr_t read_word;
    uintptr_t ulii;

    case_words[0] = 0;
    case_words[1] = 0;
    case_words[2] = 0;
    case_words[3] = 0;
    ctrl_words[0] = 0;
    ctrl_words[1] = 0;
    ctrl_words[2] = 0;
    ctrl_words[3] = 0;
    while (1)
    {
        while (rawbuf < rawbuf_end)
        {
            read_word = *rawbuf++;
            for (read_shift = 0; read_shift < read_shift_max;
                 sample_uidx++, read_shift++)
            {
                ulii = read_word & 3; // Both is_set is always true, because
                                      // dummy_nm is set
                case_words[ulii] |= ONELU << case_rem;
                if (++case_rem == BITCT)
                {
                    casebuf[offset0_case] = case_words[0];
                    casebuf[case_ctv] = case_words[2];
                    casebuf[offset2_case] = case_words[3];
                    casebuf++;
                    case_words[0] = 0;
                    case_words[2] = 0;
                    case_words[3] = 0;
                    case_rem = 0;
                }
                read_word >>= 2;
            }
        }
        if (sample_uidx == unfiltered_sample_ct)
        {
            if (case_rem)
            {
                casebuf[offset0_case] = case_words[0];
                casebuf[case_ctv] = case_words[2];
                casebuf[offset2_case] = case_words[3];
            }
            if (ctrl_rem)
            {
                ctrlbuf[offset0_ctrl] = ctrl_words[0];
                ctrlbuf[ctrl_ctv] = ctrl_words[2];
                ctrlbuf[offset2_ctrl] = ctrl_words[3];
            }
            ulii = 3;
            if (case_words[1]) { ulii -= 1; }
            if (ctrl_words[1]) { ulii -= 2; }
            *nm_info_ptr = ulii;
            return 0;
        }
        rawbuf_end++;
        read_shift_max = unfiltered_sample_ct % BITCT2;
    }
}

void Genotype::two_locus_count_table(uintptr_t* lptr1, uintptr_t* lptr2,
                                     uint32_t* counts_3x3, uint32_t sample_ctv3,
                                     uint32_t is_zmiss2)
{
#ifdef __LP64__
    uint32_t uii;
    fill_uint_zero(9, counts_3x3);
    if (!is_zmiss2)
    {
        two_locus_3x3_tablev((__m128i*) lptr1, (__m128i*) lptr2, counts_3x3,
                             sample_ctv3 / 2, 3);
    }
    else
    {
        two_locus_3x3_tablev((__m128i*) lptr2, (__m128i*) lptr1, counts_3x3,
                             sample_ctv3 / 2, 2);
        uii = counts_3x3[1];
        counts_3x3[1] = counts_3x3[3];
        counts_3x3[3] = uii;
        counts_3x3[6] = counts_3x3[2];
        counts_3x3[7] = counts_3x3[5];
    }
#else
    counts_3x3[0] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
    counts_3x3[3] =
        popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
    counts_3x3[6] =
        popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
    lptr2 = &(lptr2[sample_ctv3]);
    counts_3x3[1] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
    counts_3x3[4] =
        popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
    counts_3x3[7] =
        popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
    if (!is_zmiss2)
    {
        lptr2 = &(lptr2[sample_ctv3]);
        counts_3x3[2] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
        counts_3x3[5] =
            popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
        counts_3x3[8] = popcount_longs_intersect(
            lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
    }
#endif
}

void Genotype::two_locus_count_table_zmiss1(uintptr_t* lptr1, uintptr_t* lptr2,
                                            uint32_t* counts_3x3,
                                            uint32_t sample_ctv3,
                                            uint32_t is_zmiss2)
{

#ifdef __LP64__
    fill_uint_zero(6, counts_3x3);
    if (is_zmiss2)
    {
        two_locus_3x3_zmiss_tablev((__m128i*) lptr1, (__m128i*) lptr2,
                                   counts_3x3, sample_ctv3 / 2);
    }
    else
    {
        two_locus_3x3_tablev((__m128i*) lptr1, (__m128i*) lptr2, counts_3x3,
                             sample_ctv3 / 2, 2);
    }
#else
    counts_3x3[0] = popcount_longs_intersect(lptr1, lptr2, sample_ctv3);
    counts_3x3[1] =
        popcount_longs_intersect(lptr1, &(lptr2[sample_ctv3]), sample_ctv3);
    if (!is_zmiss2)
    {
        counts_3x3[2] = popcount_longs_intersect(
            lptr1, &(lptr2[2 * sample_ctv3]), sample_ctv3);
        counts_3x3[5] = popcount_longs_intersect(
            &(lptr1[sample_ctv3]), &(lptr2[2 * sample_ctv3]), sample_ctv3);
    }
    lptr1 = &(lptr1[sample_ctv3]);
    counts_3x3[3] = popcount_longs_intersect(lptr1, lptr2, sample_ctv3);
    counts_3x3[4] =
        popcount_longs_intersect(lptr1, &(lptr2[sample_ctv3]), sample_ctv3);
#endif
}

#ifdef __LP64__
void Genotype::two_locus_3x3_tablev(__m128i* vec1, __m128i* vec2,
                                    uint32_t* counts_3x3, uint32_t sample_ctv6,
                                    uint32_t iter_ct)
{
    const __m128i m1 = {FIVEMASK, FIVEMASK};
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    __m128i* vec20;
    __m128i* vec21;
    __m128i* vec22;
    __m128i* vend1;
    __m128i loader1;
    __m128i loader20;
    __m128i loader21;
    __m128i loader22;
    __m128i count10;
    __m128i count11;
    __m128i count12;
    __m128i count20;
    __m128i count21;
    __m128i count22;
    __univec acc0;
    __univec acc1;
    __univec acc2;
    uint32_t ct;
    uint32_t ct2;
    while (iter_ct--)
    {
        ct = sample_ctv6;
        vec20 = vec2;
        vec21 = &(vec20[sample_ctv6]);
        vec22 = &(vec20[2 * sample_ctv6]);
        while (ct >= 30)
        {
            ct -= 30;
            vend1 = &(vec1[30]);
            acc0.vi = _mm_setzero_si128();
            acc1.vi = _mm_setzero_si128();
            acc2.vi = _mm_setzero_si128();
            do
            {
            two_locus_3x3_tablev_outer:
                loader1 = *vec1++;
                loader20 = *vec20++;
                loader21 = *vec21++;
                loader22 = *vec22++;
                count10 = _mm_and_si128(loader1, loader20);
                count11 = _mm_and_si128(loader1, loader21);
                count12 = _mm_and_si128(loader1, loader22);
                count10 = _mm_sub_epi64(
                    count10, _mm_and_si128(_mm_srli_epi64(count10, 1), m1));
                count11 = _mm_sub_epi64(
                    count11, _mm_and_si128(_mm_srli_epi64(count11, 1), m1));
                count12 = _mm_sub_epi64(
                    count12, _mm_and_si128(_mm_srli_epi64(count12, 1), m1));
            two_locus_3x3_tablev_two_left:
                // unlike the zmiss variant, this apparently does not
                // suffer from enough register spill to justify
                // shrinking the inner loop
                loader1 = *vec1++;
                loader20 = *vec20++;
                loader21 = *vec21++;
                loader22 = *vec22++;
                count20 = _mm_and_si128(loader1, loader20);
                count21 = _mm_and_si128(loader1, loader21);
                count22 = _mm_and_si128(loader1, loader22);
                count20 = _mm_sub_epi64(
                    count20, _mm_and_si128(_mm_srli_epi64(count20, 1), m1));
                count21 = _mm_sub_epi64(
                    count21, _mm_and_si128(_mm_srli_epi64(count21, 1), m1));
                count22 = _mm_sub_epi64(
                    count22, _mm_and_si128(_mm_srli_epi64(count22, 1), m1));
            two_locus_3x3_tablev_one_left:
                loader1 = *vec1++;
                loader20 = *vec20++;
                loader21 = _mm_and_si128(loader1, loader20); // half1
                loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1),
                                         m1); // half2
                count10 = _mm_add_epi64(count10, _mm_and_si128(loader21, m1));
                count20 = _mm_add_epi64(count20, loader22);
                loader20 = *vec21++;
                loader21 = _mm_and_si128(loader1, loader20);
                loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
                count11 = _mm_add_epi64(count11, _mm_and_si128(loader21, m1));
                count21 = _mm_add_epi64(count21, loader22);
                loader20 = *vec22++;
                loader21 = _mm_and_si128(loader1, loader20);
                loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
                count12 = _mm_add_epi64(count12, _mm_and_si128(loader21, m1));
                count22 = _mm_add_epi64(count22, loader22);

                count10 = _mm_add_epi64(
                    _mm_and_si128(count10, m2),
                    _mm_and_si128(_mm_srli_epi64(count10, 2), m2));
                count11 = _mm_add_epi64(
                    _mm_and_si128(count11, m2),
                    _mm_and_si128(_mm_srli_epi64(count11, 2), m2));
                count12 = _mm_add_epi64(
                    _mm_and_si128(count12, m2),
                    _mm_and_si128(_mm_srli_epi64(count12, 2), m2));
                count10 = _mm_add_epi64(
                    count10,
                    _mm_add_epi64(
                        _mm_and_si128(count20, m2),
                        _mm_and_si128(_mm_srli_epi64(count20, 2), m2)));
                count11 = _mm_add_epi64(
                    count11,
                    _mm_add_epi64(
                        _mm_and_si128(count21, m2),
                        _mm_and_si128(_mm_srli_epi64(count21, 2), m2)));
                count12 = _mm_add_epi64(
                    count12,
                    _mm_add_epi64(
                        _mm_and_si128(count22, m2),
                        _mm_and_si128(_mm_srli_epi64(count22, 2), m2)));
                acc0.vi = _mm_add_epi64(
                    acc0.vi,
                    _mm_add_epi64(
                        _mm_and_si128(count10, m4),
                        _mm_and_si128(_mm_srli_epi64(count10, 4), m4)));
                acc1.vi = _mm_add_epi64(
                    acc1.vi,
                    _mm_add_epi64(
                        _mm_and_si128(count11, m4),
                        _mm_and_si128(_mm_srli_epi64(count11, 4), m4)));
                acc2.vi = _mm_add_epi64(
                    acc2.vi,
                    _mm_add_epi64(
                        _mm_and_si128(count12, m4),
                        _mm_and_si128(_mm_srli_epi64(count12, 4), m4)));
            } while (vec1 < vend1);
            const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
            acc0.vi =
                _mm_add_epi64(_mm_and_si128(acc0.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc0.vi, 8), m8));
            acc1.vi =
                _mm_add_epi64(_mm_and_si128(acc1.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
            acc2.vi =
                _mm_add_epi64(_mm_and_si128(acc2.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
            counts_3x3[0] +=
                ((acc0.u8[0] + acc0.u8[1]) * 0x1000100010001LLU) >> 48;
            counts_3x3[1] +=
                ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
            counts_3x3[2] +=
                ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
        }
        if (ct)
        {
            vend1 = &(vec1[ct]);
            ct2 = ct % 3;
            acc0.vi = _mm_setzero_si128();
            acc1.vi = _mm_setzero_si128();
            acc2.vi = _mm_setzero_si128();
            ct = 0;
            if (ct2)
            {
                count10 = _mm_setzero_si128();
                count11 = _mm_setzero_si128();
                count12 = _mm_setzero_si128();
                if (ct2 == 2) { goto two_locus_3x3_tablev_two_left; }
                count20 = _mm_setzero_si128();
                count21 = _mm_setzero_si128();
                count22 = _mm_setzero_si128();
                goto two_locus_3x3_tablev_one_left;
            }
            goto two_locus_3x3_tablev_outer;
        }
        counts_3x3 = &(counts_3x3[3]);
    }
}
#endif
