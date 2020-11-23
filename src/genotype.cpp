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

std::string Genotype::print_duplicated_snps(
    const std::unordered_set<std::string>& duplicated_snp,
    const std::string& out_prefix)
{
    // there are duplicated SNPs, we will need to terminate with the
    // information
    std::ofstream log_file_stream;
    std::string dup_name = out_prefix + ".valid";
    log_file_stream.open(dup_name.c_str());
    if (!log_file_stream.is_open())
    { throw std::runtime_error("Error: Cannot open file: " + dup_name); }
    // we should not use m_existed_snps unless it is for reference
    for (auto&& snp : m_existed_snps)
    {
        // we only output the valid SNPs.
        if (duplicated_snp.find(snp.rs()) == duplicated_snp.end()
            && (!m_has_chr_id_formula
                || duplicated_snp.find(get_chr_id(snp))
                       == duplicated_snp.end()))
            log_file_stream << snp.rs() << "\t" << snp.chr() << "\t"
                            << snp.loc() << "\t" << snp.ref() << "\t"
                            << snp.alt() << "\n";
    }
    log_file_stream.close();
    return std::string(
        "Error: A total of " + std::to_string(duplicated_snp.size())
        + " duplicated SNP ID detected out of "
        + misc::to_string(m_existed_snps.size())
        + " input SNPs! Valid SNP ID (post --extract / "
          "--exclude, non-duplicated SNPs) stored at "
        + dup_name + ". You can avoid this error by using --extract "
        + dup_name);
}

void Genotype::build_clump_windows(const unsigned long long& clump_distance)
{
    // should sort w.r.t reference
    std::sort(begin(m_existed_snps), end(m_existed_snps),
              [](SNP const& t1, SNP const& t2) {
                  if (t1.chr() == t2.chr())
                  {
                      if (t1.loc() == t2.loc())
                      {
                          if (t1.get_file_idx(true) == t2.get_file_idx(true))
                          {
                              return t1.get_byte_pos(true)
                                     < t2.get_byte_pos(true);
                          }
                          return t1.get_file_idx(true) < t2.get_file_idx(true);
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
    size_t low_bound = 0, prev_loc = 0, cur_dist;
    size_t prev_chr = ~size_t(0);
    bool first_snp = true;
    m_max_window_size = 0;
    // now we iterate thorugh all the SNPs to define the clumping window
    for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp)
    {
        auto&& cur_snp = m_existed_snps[i_snp];
        if (first_snp || prev_chr != cur_snp.chr())
        {
            if (prev_chr != cur_snp.chr())
            {
                for (size_t i = low_bound; i < i_snp; ++i)
                { m_existed_snps[i].set_up_bound(i_snp); }
            }
            prev_chr = cur_snp.chr();
            prev_loc = cur_snp.loc();
            low_bound = i_snp;
            first_snp = false;
        }
        auto snp_distance = i_snp - low_bound;
        if (m_max_window_size < snp_distance) m_max_window_size = snp_distance;
        // new chr will not go into following loop
        cur_dist = cur_snp.loc() - prev_loc;
        while (cur_dist > clump_distance && low_bound < i_snp)
        {
            snp_distance = i_snp - low_bound;
            m_existed_snps[low_bound].set_up_bound(i_snp);
            // go to next SNP
            ++low_bound;
            prev_loc = m_existed_snps[low_bound].loc();
            cur_dist = cur_snp.loc() - prev_loc;
        }
        if (m_max_window_size < snp_distance) m_max_window_size = snp_distance;
        // now low_bound should be the first SNP where the core index SNP need
        // to read from
        cur_snp.set_low_bound(low_bound);
        // set the end of the vector as the default up bound
        cur_snp.set_up_bound(m_existed_snps.size());
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

void Genotype::add_flags(const std::vector<IITree<size_t, size_t>>& gene_sets,
                         const size_t num_sets,
                         const bool genome_wide_background)
{
    const size_t required_size = BITCT_TO_WORDCT(num_sets);
    for (auto&& snp : m_existed_snps)
    {
        construct_flag(gene_sets, required_size, genome_wide_background, &snp);
    }
}

void Genotype::snp_extraction(const std::string& extract_snps,
                              const std::string& exclude_snps)
{
    if (!extract_snps.empty())
    {
        m_exclude_snp = false;
        auto input = misc::load_stream(extract_snps);
        m_snp_selection_list = load_snp_list(std::move(input));
    }
    else if (!exclude_snps.empty())
    {
        auto input = misc::load_stream(exclude_snps);
        m_snp_selection_list = load_snp_list(std::move(input));
    }
}
std::string
Genotype::get_chr_id_from_base(const BaseFile& base_file,
                               const std::vector<std::string_view>& token) const
{
    // check if we have all the column we need
    assert(!token.empty());
    std::string chr_id = "";
    for (auto col : m_chr_id_column)
    {
        if (col < 0)
        {
            auto idx = (-1 * col - 1);
            chr_id += m_chr_id_symbol[idx];
        }
        else if (!base_file.has_column[col])
        {
            throw std::runtime_error("Error: Required column for chr id "
                                     "construction not found in base file!");
        }
        else
        {
            auto idx = base_file.column_index[col];
            auto str = std::string(token[idx]);
            misc::to_upper(str);
            chr_id += str;
        }
    }
    return chr_id;
}
std::tuple<std::vector<size_t>, std::unordered_set<std::string>>
Genotype::transverse_base_file(
    const BaseFile& base_file, const QCFiltering& base_qc,
    const PThresholding& threshold_info,
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const std::streampos file_length, const bool gz_input,
    std::unique_ptr<std::istream> input)
{
    const unsigned long long max_index =
        base_file.column_index[+BASE_INDEX::MAX];
    const double max_threshold =
        threshold_info.no_full
            ? (threshold_info.fastscore ? threshold_info.bar_levels.back()
                                        : threshold_info.upper)
            : 1.0;
    double progress, prev_progress = 0.0;
    std::vector<std::string_view> token;
    std::string rs_id;
    std::string chr_id;
    std::string ref_allele;
    std::string alt_allele;
    double pvalue = 2.0;
    double stat = 0.0;
    double pthres = 0.0;
    size_t chr = 0;
    size_t loc = 0;
    unsigned long long category = 0;
    std::unordered_set<std::string> processed_rs, dup_rs;
    std::vector<size_t> filter_count(+FILTER_COUNT::MAX, 0);
    std::string line;
    while (std::getline(*input, line))
    {
        if (!gz_input)
        {
            progress = static_cast<double>(input->tellg())
                       / static_cast<double>(file_length) * 100;
            if (!m_reporter->unit_testing() && progress - prev_progress > 0.01)
            {
                fprintf(stderr, "\rReading %03.2f%%", progress);
                prev_progress = progress;
            }
        }
        misc::trim(line);
        if (line.empty()) continue;

        ++filter_count[+FILTER_COUNT::NUM_LINE];
        token = misc::tokenize(line);
        for (auto&& t : token) { misc::trim(t); }
        if (token.size() <= max_index)
        {
            throw std::runtime_error(line
                                     + "\nMore index than column in data\n");
        }

        /*
                if (SNP::use_chr_id()
                    && !parse_chr_id(token, base_file, processed_rs, dup_rs,
                                     filter_count, rs_id))
                {}
                else */
        if (!parse_rs_id(token, base_file, processed_rs, dup_rs, filter_count,
                         rs_id, chr_id))
        { continue; }
        if (!parse_chr(token, base_file, filter_count, chr)) { continue; }
        parse_allele(token, base_file, +BASE_INDEX::EFFECT, ref_allele);
        parse_allele(token, base_file, +BASE_INDEX::NONEFFECT, alt_allele);
        if (!parse_loc(token, base_file, loc))
        {
            throw std::runtime_error(
                "Error: Invalid loci for " + rs_id + ": "
                + std::string(token[base_file.column_index[+BASE_INDEX::BP]])
                + "\n");
        }
        if (base_file.has_column[+BASE_INDEX::BP]
            && base_file.has_column[+BASE_INDEX::CHR]
            && Genotype::within_region(exclusion_regions, chr, loc))
        {
            ++filter_count[+FILTER_COUNT::REGION];
            continue;
        }

        if (!base_filter_by_value(token, base_file, base_qc.maf, filter_count,
                                  +FILTER_COUNT::MAF, +BASE_INDEX::MAF)
            || !base_filter_by_value(token, base_file, base_qc.maf_case,
                                     filter_count, +FILTER_COUNT::MAF,
                                     +BASE_INDEX::MAF_CASE))
        { continue; }
        if (!base_filter_by_value(token, base_file, base_qc.info_score,
                                  filter_count, +FILTER_COUNT::INFO,
                                  +BASE_INDEX::INFO))
        { continue; }
        if (!parse_pvalue(token[base_file.column_index[+BASE_INDEX::P]],
                          max_threshold, filter_count, pvalue))
        { continue; }
        if (!parse_stat(token[base_file.column_index[+BASE_INDEX::STAT]],
                        base_file.is_or, filter_count, stat))
        { continue; }
        if (!alt_allele.empty() && ambiguous(ref_allele, alt_allele))
        {
            ++filter_count[+FILTER_COUNT::AMBIG];
            if (!m_keep_ambig) continue;
        }
        category = 0;
        pthres = 0.0;
        if (threshold_info.fastscore)
        {
            category =
                cal_bar_category(pvalue, threshold_info.bar_levels, pthres);
        }
        else
        {
            try
            {
                category = calculate_category(threshold_info, pvalue, pthres);
            }
            catch (const std::runtime_error&)
            {
                m_very_small_thresholds = true;
                category = 0;
            }
        }
        m_existed_snps_index[rs_id] = m_existed_snps.size();
        // we should also load the chr_id
        if (!chr_id.empty())
            m_existed_snps_index[chr_id] = m_existed_snps.size();
        m_existed_snps.emplace_back(SNP(rs_id, chr, loc, ref_allele, alt_allele,
                                        stat, pvalue, category, pthres));
    }
    if (!m_reporter->unit_testing())
    { fprintf(stderr, "\rReading %03.2f%%\n", 100.0); }
    input.reset();
    return {filter_count, dup_rs};
}
std::tuple<std::vector<size_t>, std::unordered_set<std::string>>
Genotype::read_base(
    const BaseFile& base_file, const QCFiltering& base_qc,
    const PThresholding& threshold_info,
    const std::vector<IITree<size_t, size_t>>& exclusion_regions)
{
    std::string line;
    std::string message = "Base file: " + base_file.file_name + "\n";
    std::streampos file_length = 0;
    bool gz_input;
    auto stream = misc::load_stream(base_file.file_name, gz_input);
    if (!gz_input)
    {
        stream->seekg(0, stream->end);
        file_length = stream->tellg();
        stream->clear();
        stream->seekg(0, stream->beg);
    }
    else
    {
        message.append("GZ file detected. ");
    }
    if (!base_file.is_index)
    {
        std::getline(*stream, line);
        message.append("Header of file is:\n" + line + "\n\n");
    }
    m_reporter->report(message);
    message.clear();
    return transverse_base_file(base_file, base_qc, threshold_info,
                                exclusion_regions, file_length, gz_input,
                                std::move(stream));
}


void Genotype::print_base_stat(const std::vector<size_t>& filter_count,
                               const std::unordered_set<std::string>& dup_index,
                               const std::string& out, const double info_score)
{
    std::string message = std::to_string(filter_count[+FILTER_COUNT::NUM_LINE])
                          + " variant(s) observed in base file, with:\n";
    if (filter_count[+FILTER_COUNT::SELECT])
    {
        message.append(std::to_string(filter_count[+FILTER_COUNT::SELECT])
                       + " variant(s) excluded based on user input\n");
    }
    if (filter_count[+FILTER_COUNT::CHR])
    {
        message.append(
            std::to_string(filter_count[+FILTER_COUNT::CHR])
            + " variant(s) excluded as they are on unknown/sex chromosome\n");
    }
    if (filter_count[+FILTER_COUNT::HAPLOID])
    {
        message.append(std::to_string(filter_count[+FILTER_COUNT::HAPLOID])
                       + " variant(s) located on haploid chromosome\n");
    }
    if (filter_count[+FILTER_COUNT::REGION])
    {
        message.append(
            std::to_string(filter_count[+FILTER_COUNT::REGION])
            + " variant(s) excluded as they fall within x-range region(s)\n");
    }
    if (filter_count[+FILTER_COUNT::MAF])
    {
        message.append(std::to_string(filter_count[+FILTER_COUNT::MAF])
                       + " variant(s) excluded due to MAF threshold\n");
    }
    if (filter_count[+FILTER_COUNT::INFO])
    {
        message.append(std::to_string(filter_count[+FILTER_COUNT::INFO])
                       + " variant(s) with INFO score less than "
                       + std::to_string(info_score) + "\n");
    }
    if (filter_count[+FILTER_COUNT::P_EXCLUDED])
    {
        message.append(std::to_string(filter_count[+FILTER_COUNT::P_EXCLUDED])
                       + " variant(s) excluded due to p-value threshold\n");
    }
    if (filter_count[+FILTER_COUNT::NOT_CONVERT])
    {
        message.append(std::to_string(filter_count[+FILTER_COUNT::NOT_CONVERT])
                       + " NA stat/p-value observed\n");
    }
    if (filter_count[+FILTER_COUNT::NEGATIVE])
    {
        message.append(std::to_string(filter_count[+FILTER_COUNT::NEGATIVE])
                       + " negative statistic observed. Maybe you have "
                         "forgotten the --beta flag?\n");
    }
    if (filter_count[+FILTER_COUNT::AMBIG])
    {
        message.append(std::to_string(filter_count[+FILTER_COUNT::AMBIG])
                       + " ambiguous variant(s)");
        if (!m_keep_ambig) { message.append(" excluded"); }
        message.append("\n");
    }
    if (filter_count[+FILTER_COUNT::DUP_SNP])
    {
        if (!message.empty()) { m_reporter->report(message); }
        throw std::runtime_error(print_duplicated_snps(dup_index, out));
    }
    message.append(std::to_string(m_existed_snps.size())
                   + " total variant(s) included from base file\n\n");
    m_reporter->report(message);
    if (m_existed_snps.size() == 0)
    { throw std::runtime_error("Error: No valid variant remaining"); }
}

bool Genotype::has_parent(const std::unordered_set<std::string>& founder_info,
                          const std::vector<std::string>& token,
                          const std::string& fid, const size_t idx)
{

    if (idx == ~size_t(0)) return false;
    auto found = founder_info.find(fid + m_delim + token.at(idx));
    return found != founder_info.end();
}


bool Genotype::check_chr(const std::string& chr_str, std::string& prev_chr,
                         size_t& chr_num, bool& chr_error, bool& sex_error)
{
    if (chr_str != prev_chr)
    {
        auto chr_code = get_chrom_code(chr_str);
        if (chr_code < 0)
        {
            if (!chr_error)
                m_reporter->report("Error: Invalid chromosome number for SNP");
            chr_error = true;
            return false;
        }
        if (chr_code > static_cast<int>(m_autosome_ct)
            || is_set(m_haploid_mask.data(), static_cast<uint32_t>(chr_code)))
        {
            // this is sex / mt chromosome
            if (!sex_error)
                m_reporter->report("Warning: Currently not support "
                                   "haploid chromosome and sex "
                                   "chromosomes\n");
            sex_error = true;
            return false;
        }
        chr_num = static_cast<size_t>(chr_code);
        prev_chr = chr_str;
    }
    return true;
}
bool Genotype::check_rs(const std::string& snpid, const std::string& chrid,
                        std::string& rsid,
                        std::unordered_set<std::string>& processed_snps,
                        std::unordered_set<std::string>& duplicated_snps,
                        Genotype* genotype)
{
    if ((snpid.empty() || snpid == ".") && (rsid.empty() || rsid == "."))
    {
        ++m_base_missed;
        return false;
    }
    auto&& find_rs = genotype->m_existed_snps_index.find(rsid);
    if (find_rs == genotype->m_existed_snps_index.end())
    {
        if (snpid.empty()
            || genotype->m_existed_snps_index.find(snpid)
                   == genotype->m_existed_snps_index.end())
        {
            if (chrid.empty()
                || genotype->m_existed_snps_index.find(chrid)
                       == genotype->m_existed_snps_index.end())
            {
                ++m_base_missed;
                return false;
            }
            rsid = chrid;
        }
        else
        {
            rsid = snpid;
        }
    }
    if (processed_snps.find(rsid) != processed_snps.end())
    {
        // no need to add m_base_missed as this will completley error out
        duplicated_snps.insert(rsid);
        return false;
    }
    return true;
}

bool Genotype::check_ambig(const std::string& a1, const std::string& a2,
                           const std::string& ref, bool& flipping)
{
    bool ambig = ambiguous(a1, a2);
    if (ambig)
    {
        ++m_num_ambig;
        if (!m_keep_ambig) return false;
        // we assume allele matching has been done
        // therefore ref = a1 or complementary ref = a1
        // or alt = a1 or complementary(alt = a1)
        // but as this is ambiguous, the complemenary check will always return
        // true AT match to TA without need of flipping (as comp TA = AT)
        // but we want to flip this (we want AT to match to AT)
        // so we will check if a1 == ref and only flip when that is not true
        flipping = (a1 != ref);
    }
    return true;
}
bool Genotype::not_in_xregion(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const SNP& base, const SNP& target)
{
    if (base.chr() != ~size_t(0) && base.loc() != ~size_t(0)) return true;
    if (Genotype::within_region(exclusion_regions, target.chr(), target.loc()))
    {
        ++m_num_xrange;
        return false;
    }
    return true;
}
std::string Genotype::chr_id_from_genotype(const SNP& snp) const
{
    std::string chr_id = "";
    if (!m_has_chr_id_formula) return chr_id;

    for (auto&& col : m_chr_id_column)
    {
        if (col < 0)
        {
            auto idx = -1 * col - 1;
            chr_id += m_chr_id_symbol[idx];
        }
        else
        {
            switch (col)
            {
            case +BASE_INDEX::CHR: chr_id += std::to_string(snp.chr()); break;
            case +BASE_INDEX::BP: chr_id += std::to_string(snp.loc()); break;
            case +BASE_INDEX::EFFECT: chr_id += snp.ref(); break;
            case +BASE_INDEX::NONEFFECT: chr_id += snp.alt(); break;
            default: throw std::logic_error("Error: This should never happen");
            }
        }
    }
    return chr_id;
}
bool Genotype::process_snp(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const std::string& mismatch_snp_record_name,
    const std::string& mismatch_source, const std::string& snpid, SNP& snp,
    std::unordered_set<std::string>& processed_snps,
    std::unordered_set<std::string>& duplicated_snps,
    std::vector<bool>& retain_snp, Genotype* genotype)
{
    misc::to_upper(snp.ref());
    misc::to_upper(snp.alt());
    auto chr_id = chr_id_from_genotype(snp);
    if (!check_rs(snpid, chr_id, snp.rs(), processed_snps, duplicated_snps,
                  genotype))
        return false;
    auto snp_idx = genotype->m_existed_snps_index[snp.rs()];
    auto&& target_snp = genotype->m_existed_snps[snp_idx];

    bool flipping = false;
    if (!target_snp.matching(snp, flipping))
    {
        genotype->print_mismatch(mismatch_snp_record_name, mismatch_source,
                                 target_snp, snp);
        ++m_num_ref_target_mismatch;
        return false;
    }
    if (!check_ambig(snp.ref(), snp.alt(), target_snp.ref(), flipping))
        return false;
    // only do region test if we know we haven't done it during read_base
    // we will do it in read_base if we have chr and loc info.
    if (!not_in_xregion(exclusion_regions, target_snp, snp)) { return false; }
    //  only add valid SNPs
    processed_snps.insert(snp.rs());
    target_snp.add_snp_info(snp, flipping, m_is_ref);
    retain_snp[snp_idx] = true;
    return true;
}

void Genotype::gen_sample(const size_t fid_idx, const size_t iid_idx,
                          const size_t sex_idx, const size_t dad_idx,
                          const size_t mum_idx, const size_t cur_idx,
                          const std::unordered_set<std::string>& founder_info,
                          const std::string& pheno,
                          std::vector<std::string>& token,
                          std::vector<Sample_ID>& sample_storage,
                          std::unordered_set<std::string>& sample_in_file,
                          std::vector<std::string>& duplicated_sample_id)
{
    assert(m_vector_initialized);
    for (size_t i = 0; i < token.size(); ++i) { misc::trim(token[i]); }
    // we have already checked for malformed file
    const std::string fid = (m_ignore_fid) ? "-" : token[fid_idx];
    const std::string id =
        (m_ignore_fid) ? token[iid_idx] : fid + m_delim + token[iid_idx];
    if (m_max_fid_length < token[fid_idx].length())
        m_max_fid_length = token[fid_idx].length();
    if (m_max_iid_length < token[iid_idx].length())
        m_max_iid_length = token[iid_idx].length();

    // end immediately if duplicated samples are found
    if (sample_in_file.find(id) != sample_in_file.end())
    {
        duplicated_sample_id.push_back(id);
        return;
    }
    auto&& find_id =
        m_sample_selection_list.find(id) != m_sample_selection_list.end();
    bool inclusion = m_remove_sample ^ find_id;
    bool in_regression = false;
    // we can't check founder if there isn't fid

    const bool has_father = has_parent(founder_info, token, fid, dad_idx);
    const bool has_mother = has_parent(founder_info, token, fid, mum_idx);
    if (inclusion)
    {
        if (!m_ignore_fid && (has_father || has_mother))
        {
            // we still calculate PRS for this sample
            SET_BIT(cur_idx, m_calculate_prs.data());
            ++m_num_non_founder;
            // but will only include it in the regression model if users asked
            // to include non-founders
            in_regression = m_keep_nonfounder;
        }
        else
        {
            ++m_founder_ct;
            SET_BIT(cur_idx, m_sample_for_ld.data());
            SET_BIT(cur_idx, m_calculate_prs.data());
            in_regression = true;
        }
    }
    m_sample_ct += inclusion;
    // TODO: Better sex parsing? Can also be 0, 1 or F and M
    if (sex_idx != ~size_t(0) && token[sex_idx] == "1") { ++m_num_male; }
    else if (sex_idx != ~size_t(0) && token[sex_idx] == "2")
    {
        ++m_num_female;
    }
    else
    {
        ++m_num_ambig_sex;
    }
    // this must be incremented within each loop
    if (inclusion && !m_is_ref)
    {
        sample_storage.emplace_back(
            Sample_ID(token[fid_idx], token[iid_idx], pheno, in_regression));
    }
    sample_in_file.insert(id);
}

std::vector<std::string>
Genotype::load_genotype_prefix(std::unique_ptr<std::istream> in)
{
    std::vector<std::string> genotype_files;
    std::string line;
    while (std::getline(*in, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        genotype_files.push_back(line);
    }
    // we no longer need in
    in.reset();
    return genotype_files;
}

void Genotype::init_chr(int num_auto, bool no_x, bool no_y, bool no_xy,
                        bool no_mt)
{
    // this initialize haploid mask as the maximum possible number
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    m_haploid_mask.resize(CHROM_MASK_WORDS, 0);
    if (num_auto < 0)
    {
        // not required at the moment
        /*
        num_auto = -num_auto;
        m_autosome_ct = static_cast<uint32_t>(num_auto);
        m_xymt_codes[X_OFFSET] = -1;
        m_xymt_codes[Y_OFFSET] = -1;
        m_xymt_codes[XY_OFFSET] = -1;
        m_xymt_codes[MT_OFFSET] = -1;
        m_max_code = static_cast<uint32_t>(num_auto);
        fill_all_bits((static_cast<uint32_t>(num_auto)) + 1,
                      m_haploid_mask.data());*/
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

size_t Genotype::get_rs_column(const std::string& input)
{
    std::vector<std::string> token = misc::split(input);
    // rs_index store the location of the RS ID
    if (token.size() != 1)
    {
        size_t rs_index = 0;
        for (auto&& name : token)
        {
            misc::to_upper(name);
            if (name == "SNP" || name == "RS" || name == "RS_ID"
                || name == "RS.ID" || name == "RSID" || name == "VARIANT.ID"
                || name == "VARIANT_ID" || name == "SNP_ID" || name == "SNP.ID")
            {
                m_reporter->report(name
                                   + " assume to be column containing SNP ID");
                return rs_index;
            }
            ++rs_index;
        }
        // if user accidentally add an empty column at the back of their bim
        // file
        // now it is possible that the file has an empty column for header, but
        // that will just be an edge case and can be easily solved by adding the
        // appropriate header
        if (token.size() == 6 || (token.size() == 7 && token.back().empty()))
        {
            m_reporter->report(
                "SNP extraction/exclusion list contains 6 columns, "
                "will assume this is a bim file, with the "
                "second column contains the SNP ID");
            return 1;
        }
        else
        {
            m_reporter->report(
                "SNP extraction/exclusion list contains "
                + misc::to_string(token.size())
                + " columns, "
                  "will assume first column contains the SNP ID");
            return 0;
        }
    }
    else
    {
        m_reporter->report(
            "Only one column detected, will assume only SNP ID is provided");
        return 0;
    }
}

std::unordered_set<std::string>
Genotype::load_snp_list(std::unique_ptr<std::istream> input)
{
    std::string line;
    std::getline(*input, line);
    input->clear();
    input->seekg(0, std::ios::beg);
    misc::trim(line);
    size_t rs_index = get_rs_column(line);
    std::vector<std::string_view> token;
    std::unordered_set<std::string> result;
    while (std::getline(*input, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::tokenize(line);
        result.insert(std::string(token[rs_index]));
    }
    input.reset();
    return result;
}

std::unordered_set<std::string>
Genotype::load_ref(std::unique_ptr<std::istream> input, bool ignore_fid)
{
    std::string line;
    // now go through the sample file. We require the FID (if any) and IID  must
    // be the first 1/2 column of the file
    std::vector<std::string> token;
    std::unordered_set<std::string> result;
    while (std::getline(*input, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (ignore_fid) { result.insert(token[0]); }
        else
        {
            if (token.size() < 2)
                throw std::runtime_error(
                    "Error: Require FID and IID for extraction. "
                    "You can ignore the FID by using the --ignore-fid flag");
            result.insert(token[0] + m_delim + token[1]);
        }
    }
    input.reset();
    return result;
}

void Genotype::load_samples(bool verbose)
{
    if (!m_remove_file.empty())
    {
        auto input = misc::load_stream(m_remove_file);
        m_sample_selection_list = load_ref(std::move(input), m_ignore_fid);
    }
    else if (!m_keep_file.empty())
    {
        m_remove_sample = false;
        auto input = misc::load_stream(m_keep_file);
        m_sample_selection_list = load_ref(std::move(input), m_ignore_fid);
    }
    if (!m_is_ref) { m_sample_id = gen_sample_vector(); }
    else
    {
        // don't bother loading up the sample vector as it should
        // never be used for reference panel (except for the founder_info)
        gen_sample_vector();
    }
    std::string message = misc::to_string(m_unfiltered_sample_ct) + " people ("
                          + misc::to_string(m_num_male) + " male(s), "
                          + misc::to_string(m_num_female)
                          + " female(s)) observed\n";
    message.append(misc::to_string(m_founder_ct) + " founder(s) included");

    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    m_exclude_from_std.resize(unfiltered_sample_ctl, 0);
    if (verbose) m_reporter->report(message);
}

bool Genotype::perform_freqs_and_inter(const QCFiltering& filter_info,
                                       const std::string& prefix,
                                       Genotype* target)
{
    if (!m_intermediate
        && (misc::logically_equal(filter_info.geno, 1.0)
            || filter_info.geno > 1.0)
        && (misc::logically_equal(filter_info.maf, 0.0)
            || filter_info.maf > 1.0)
        && (misc::logically_equal(filter_info.info_score, 0.0)
            || filter_info.info_score < 0.0))
    { return false; }
    const std::string print_target = (m_is_ref) ? "reference" : "target";
    m_reporter->report("Calculate MAF and perform filtering on " + print_target
                       + " SNPs\n"
                         "==================================================");
    auto&& genotype = (m_is_ref) ? target : this;
    std::sort(
        begin(genotype->m_existed_snps), end(genotype->m_existed_snps),
        [this](SNP const& t1, SNP const& t2) {
            if (t1.get_file_idx(m_is_ref) == t2.get_file_idx(m_is_ref))
            { return t1.get_byte_pos(m_is_ref) < t2.get_byte_pos(m_is_ref); }
            else
                return t1.get_file_idx(m_is_ref) == t2.get_file_idx(m_is_ref);
        });
    return calc_freq_gen_inter(filter_info, prefix, genotype);
}

void Genotype::calc_freqs_and_intermediate(const QCFiltering& filter_info,
                                           const std::string& prefix,
                                           bool verbose, Genotype* target)
{
    std::string message = "";
    m_num_geno_filter = 0;
    m_num_maf_filter = 0;
    m_num_info_filter = 0;
    // only print the filtering message if filtering was performed
    if (perform_freqs_and_inter(filter_info, prefix, target))
    {
        m_marker_ct = m_existed_snps.size();
        if (m_num_geno_filter != 0)
        {
            message.append(std::to_string(m_num_geno_filter)
                           + " variant(s) excluded based on genotype "
                             "missingness threshold\n");
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
        if (m_num_miss_filter != 0)
        {
            message.append(std::to_string(m_num_miss_filter)
                           + " variant(s) excluded as they are completely "
                             "missed on all founder samples\n");
        }
        if (!m_is_ref)
        {
            message.append(std::to_string(m_marker_ct)
                           + " variant(s) included");
        }
        else
        {
            message.append(std::to_string(target->m_existed_snps.size())
                           + " variant(s) remained");
        }
        if (verbose) m_reporter->report(message);
    }
}

void Genotype::print_mismatch(const std::string& out, const std::string& type,
                              const SNP& target, const SNP& new_snp)
{
    // mismatch found between base target and reference
    if (!m_mismatch_snp_record.is_open())
    {
        m_mismatch_snp_record.open(out.c_str());
        if (!m_mismatch_snp_record.is_open())
        {
            throw std::runtime_error(std::string("Cannot open mismatch file to "
                                                 "write: "
                                                 + out));
        }
        m_mismatch_snp_record << "File_Type\tRS_ID\tCHR_Target\tCHR_"
                                 "File\tBP_Target\tBP_File\tA1_"
                                 "Target\tA1_File\tA2_Target\tA2_"
                                 "File\n";
    }
    m_mismatch_snp_record << type << "\t" << new_snp.rs() << "\t"
                          << new_snp.chr() << "\t";
    if (target.chr() == ~size_t(0)) { m_mismatch_snp_record << "-\t"; }
    else
    {
        m_mismatch_snp_record << target.chr() << "\t";
    }
    m_mismatch_snp_record << new_snp.loc() << "\t";
    if (target.loc() == ~size_t(0)) { m_mismatch_snp_record << "-\t"; }
    else
    {
        m_mismatch_snp_record << target.loc() << "\t";
    }
    m_mismatch_snp_record << new_snp.ref() << "\t" << target.ref() << "\t"
                          << new_snp.alt() << "\t" << target.alt() << std::endl;
}

void Genotype::load_snps(
    const std::string& out,
    const std::vector<IITree<size_t, size_t>>& exclusion_regions, bool verbose,
    Genotype* target)
{
    m_base_missed = 0;
    m_num_ambig = 0;
    m_num_xrange = 0;
    gen_snp_vector(exclusion_regions, out, target);
    auto&& snp_store_location = m_is_ref ? target : this;
    snp_store_location->m_marker_ct = snp_store_location->m_existed_snps.size();
    std::string message = "";
    if (m_base_missed != 0)
    {
        message.append(std::to_string(m_base_missed)
                       + " variant(s) not found in previous data\n");
    }
    std::string action = "excluded";
    if (m_keep_ambig) action = "kept";
    if (m_num_ref_target_mismatch != 0)
    {
        message.append(std::to_string(m_num_ref_target_mismatch)
                       + " variant(s) with mismatch information\n");
    }
    if (m_num_ambig != 0)
    {
        message.append(std::to_string(m_num_ambig) + " ambiguous variant(s) "
                       + action + "\n");
    }
    if (m_num_xrange != 0)
    {
        message.append(std::to_string(m_num_xrange)
                       + " variant(s) removed as they fall within the "
                         "--x-range region(s)\n");
    }
    message.append(std::to_string(snp_store_location->m_marker_ct)
                   + " variant(s) included\n");

    if (verbose) m_reporter->report(message);
    m_snp_selection_list.clear();
    if (snp_store_location->m_marker_ct == 0)
    {
        message = "Error: No vairant remained!\n";
        throw std::runtime_error(message);
    }
    // only remove selection list here, as we might need this for bgen file
    // check
    m_sample_selection_list.clear();
}


Genotype::~Genotype() {}
intptr_t Genotype::cal_avail_memory(const uintptr_t founder_ctv2)
{
#ifdef __APPLE__
    int32_t mib[2];
    size_t sztmp;
#endif
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
    std::string message = "";
    if (llxx)
    {
        // we have detected the memory, but need to check if that's enough
        if (malloc_size_mb > llxx)
        {
            throw std::runtime_error(
                "Error: Insufficient memory for clumping! Require "
                + misc::to_string(malloc_size_mb) + " MB but detected only "
                + misc::to_string(llxx) + " MB");
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
    m_reporter->report(message);
    return malloc_size_mb;
}
std::vector<std::pair<size_t, size_t>> Genotype::get_chrom_boundary()
{
    std::vector<std::pair<size_t, size_t>> chrom_bound;
    std::uint32_t prev_chr = 0;
    auto total = m_existed_snps.size();
    for (size_t i = 0; i < total; ++i)
    {
        auto&& cur_snp = m_existed_snps[m_sort_by_p_index[i]];
        if (cur_snp.chr() != prev_chr)
        {
            if (!chrom_bound.empty())
            {
                // new chr
                std::get<1>(chrom_bound.back()) = i;
            }
            chrom_bound.push_back(std::pair<size_t, size_t>(i, total));
            prev_chr = cur_snp.chr();
        }
    }

    return chrom_bound;
}
void Genotype::clumping(const Clumping& clump_info, Genotype& reference,
                        size_t threads)
{
    m_reporter->report("Start performing clumping");
    std::vector<std::atomic<bool>> remain_snps(m_existed_snps.size());
    for (auto&& s : remain_snps) { s = false; }
    std::atomic<size_t> num_core = 0;
    using range = std::pair<size_t, size_t>;
    // get boundaries
    std::vector<range> snp_range = get_chrom_boundary();
    if (threads == 1)
    {
        dummy_reporter progress_reporter(m_existed_snps.size(),
                                         !m_reporter->unit_testing());
        threaded_clumping(snp_range, clump_info, progress_reporter, remain_snps,
                          num_core, reference);
    }
    else
    {
        if (threads > m_autosome_ct) { threads = m_autosome_ct; }
        Thread_Queue<size_t> progress_observer;
        std::thread observer(&Genotype::clump_progress_observer, this,
                             std::ref(progress_observer), m_existed_snps.size(),
                             threads, !m_reporter->unit_testing());
        std::vector<std::thread> subjects;
        size_t job_per_thread = snp_range.size() / static_cast<size_t>(threads);
        int remain = static_cast<int>(static_cast<size_t>(snp_range.size())
                                      % static_cast<size_t>(threads));
        size_t job_start = 0;
        for (size_t i_thread = 0; i_thread < threads; ++i_thread)
        {
            std::vector<range> job_sets(snp_range.begin() + job_start,
                                        snp_range.begin() + job_start
                                            + job_per_thread + (remain > 0));
            subjects.push_back(
                std::thread(&Genotype::threaded_clumping<Thread_Queue<size_t>>,
                            this, job_sets, std::cref(clump_info),
                            std::ref(progress_observer), std::ref(remain_snps),
                            std::ref(num_core), std::ref(reference)));
            job_start += job_per_thread + (remain > 0);
            remain--;
        }
        observer.join();
        for (auto&& thread : subjects) thread.join();
    }
    if (!m_reporter->unit_testing())
    { fprintf(stderr, "\rClumping Progress: %03.2f%%\n", 100.0); }
    /*
        std::vector<bool> non_atomic_remain(m_existed_snps.size());
        for (size_t i = 0; i < remain_snps.size(); ++i)
        { non_atomic_remain[i] = remain_snps[i]; }*/
    if (num_core != m_existed_snps.size()) { shrink_snp_vector(remain_snps); }
    m_existed_snps_index.clear();
    m_reporter->report("Number of variant(s) after clumping : "
                       + misc::to_string(m_existed_snps.size()));
}

void Genotype::clump_progress_observer(Thread_Queue<size_t>& progress_observer,
                                       size_t total_snp, size_t num_thread,
                                       bool verbose)
{
    size_t progress = 0, total_run = 0;
    double cur_progress = 0, prev_progress = 0;
    while (!progress_observer.pop(progress, num_thread))
    {
        total_run += progress;
        cur_progress = static_cast<double>(total_run)
                       / static_cast<double>(total_snp) * 100;
        if (verbose && cur_progress - prev_progress > 0.01)
        {
            fprintf(stderr, "\rClumping Progress: %03.2f%%", cur_progress);
            prev_progress = cur_progress;
        }
    }
}

template <typename T>
void Genotype::threaded_clumping(
    const std::vector<std::pair<size_t, size_t>> snp_range,
    const Clumping& clump_info, T& progress_observer,
    std::vector<std::atomic<bool>>& remain_snps, std::atomic<size_t>& num_core,
    Genotype& reference)
{
    const double min_r2 = clump_info.use_proxy
                              ? std::min(clump_info.proxy, clump_info.r2)
                              : clump_info.r2;
    const uint32_t founder_ctv3 =
        BITCT_TO_ALIGNED_WORDCT(static_cast<uint32_t>(reference.m_founder_ct));
    const uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(reference.m_founder_ct);
    const uint32_t founder_ctsplit = 3 * founder_ctv3;
    const uintptr_t founder_ctv2 =
        QUATERCT_TO_ALIGNED_WORDCT(reference.m_founder_ct);
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(reference.m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    std::vector<uintptr_t> index_data(3 * founder_ctsplit + founder_ctv3);
    std::vector<uintptr_t> index_tots(6);

    std::vector<uintptr_t> founder_include2(founder_ctv2, 0);
    fill_quatervec_55(static_cast<uint32_t>(reference.m_founder_ct),
                      founder_include2.data());

    // available memory size
    // can set number to way higher with ultra (e.g. all SNPs)
    size_t num_snp_in_chr = 0;
    size_t max_snp_in_chr = 0;
    // only reserves the max amount of SNPs required for all chromosome involved
    for (auto&& range : snp_range)
    {
        num_snp_in_chr = std::get<1>(range) - std::get<0>(range);
        if (num_snp_in_chr > max_snp_in_chr) max_snp_in_chr = num_snp_in_chr;
    }
    // define the max size as the maximum window size if it is bigger than 20%
    // of the size of current chromosome. Otherwise, define the max size as 20%
    // of the chromosome.
    const auto max_size = m_max_window_size > std::floor(max_snp_in_chr * 0.2)
                              ? m_max_window_size
                              : std::floor(max_snp_in_chr * 0.2);
    GenotypePool genotype_pool(max_size + 1, unfiltered_sample_ctv2);
    auto tmp_genotype = genotype_pool.alloc();
    double r2 = -1;
    FileRead genotype_file;
    size_t num_processed = 0, prev_processed = 0;
    double local_progress = 0.0, prev_progress = 0.0;
    size_t local_num_core = 0;
    auto&& sample_for_ld = reference.m_sample_for_ld.data();
    for (auto&& range : snp_range)
    {
        for (size_t i_snp = std::get<0>(range); i_snp < std::get<1>(range);
             ++i_snp)
        {
            local_progress = static_cast<double>(num_processed - num_snp_in_chr)
                             / static_cast<double>(num_snp_in_chr);
            if (local_progress - prev_progress > 0.01)
            {
                progress_observer.emplace(num_processed - prev_processed);
                prev_progress = num_processed;
                prev_processed = num_processed;
            }

            // first implement it without the progress bar
            auto&& core_snp_idx = m_sort_by_p_index[i_snp];
            auto&& core_snp = m_existed_snps[core_snp_idx];
            if (core_snp.clumped() || core_snp.p_value() > clump_info.pvalue)
            { continue; }
            const size_t clump_start_idx = core_snp.low_bound();
            const size_t clump_end_idx = core_snp.up_bound();
            // the reason this is a two part process is so that we can reduce
            // the number of fseek
            for (size_t clump_idx = clump_start_idx; clump_idx < core_snp_idx;
                 ++clump_idx)
            {
                auto&& clump_snp = m_existed_snps[clump_idx];
                if (clump_snp.clumped()
                    || clump_snp.p_value() > clump_info.pvalue)
                { continue; }
                if (clump_snp.current_genotype() == nullptr)
                {
                    // store clump SNP's genotype into our genotype pool
                    clump_snp.set_genotype_storage(genotype_pool.alloc());
                    reference.read_genotype(
                        clump_snp, reference.m_founder_ct, genotype_file,
                        tmp_genotype->get_geno(), clump_snp.current_genotype(),
                        sample_for_ld, true);
                }
            }
            if (core_snp.current_genotype() == nullptr)
            {
                // store core SNP's genotype into our genotype pool
                core_snp.set_genotype_storage(genotype_pool.alloc());
                reference.read_genotype(core_snp, reference.m_founder_ct,
                                        genotype_file, tmp_genotype->get_geno(),
                                        core_snp.current_genotype(),
                                        sample_for_ld, true);
            }
            update_index_tot(founder_ctl2, founder_ctv2, reference.m_founder_ct,
                             index_data, index_tots, founder_include2,
                             core_snp.current_genotype());
            // free core SNP's genotype form the genotype pool as we will no
            // longer need it. (it will never be clumped by another SNP)
            core_snp.freed_geno_storage(genotype_pool);

            for (size_t clump_idx = clump_start_idx; clump_idx < core_snp_idx;
                 ++clump_idx)
            {
                auto&& clump_snp = m_existed_snps[clump_idx];
                if (clump_snp.clumped()
                    || clump_snp.p_value() > clump_info.pvalue)
                { continue; }
                r2 = get_r2(founder_ctl2, founder_ctv2,
                            clump_snp.current_genotype(), index_data,
                            index_tots);

                if (r2 >= min_r2)
                {
                    core_snp.clump(clump_snp, r2, clump_info.use_proxy,
                                   clump_info.proxy);
                    // remove SNP's genotype data from the genotype pool if it
                    // is clumped out
                    if (clump_snp.clumped())
                    { clump_snp.freed_geno_storage(genotype_pool); }
                }
            }
            // now we can read the SNPs that come after the index SNP in the
            // file
            for (size_t clump_idx = core_snp_idx + 1; clump_idx < clump_end_idx;
                 ++clump_idx)
            {
                auto&& clump_snp = m_existed_snps[clump_idx];
                if (clump_snp.clumped()
                    || clump_snp.p_value() > clump_info.pvalue)
                    continue;
                if (clump_snp.current_genotype() == nullptr)
                {
                    // store clump SNP's genotype into our genotype pool
                    clump_snp.set_genotype_storage(genotype_pool.alloc());
                    reference.read_genotype(
                        clump_snp, reference.m_founder_ct, genotype_file,
                        tmp_genotype->get_geno(), clump_snp.current_genotype(),
                        sample_for_ld, true);
                }
                r2 = get_r2(founder_ctl2, founder_ctv2,
                            clump_snp.current_genotype(), index_data,
                            index_tots);

                if (r2 >= min_r2)
                {
                    // remove SNP's genotype data from the genotype pool if it
                    // is clumped out
                    core_snp.clump(clump_snp, r2, clump_info.use_proxy,
                                   clump_info.proxy);
                    if (clump_snp.clumped())
                    { clump_snp.freed_geno_storage(genotype_pool); }
                }
            }
            core_snp.set_clumped();
            // we set the remain_core to true so that we will keep it at the end
            remain_snps[core_snp_idx] = true;
            ++num_processed;
            ++local_num_core;
        }
        // in theory, by the time we reached here, genotype_pool should be empty
        // as all SNPs should either be clumped out or are index
    }

    progress_observer.completed();
    num_core += local_num_core;
    genotype_pool.free(tmp_genotype);
}
void Genotype::recalculate_categories(const PThresholding& p_info)
{ // need to loop through the SNPs to check
    std::sort(begin(m_existed_snps), end(m_existed_snps),
              [](SNP const& t1, SNP const& t2) {
                  if (misc::logically_equal(t1.p_value(), t2.p_value()))
                  {
                      if (t1.chr() == t2.chr())
                      {
                          if (t1.loc() == t2.loc())
                          { return t1.rs() < t2.rs(); }
                          else
                              return t1.loc() < t2.loc();
                      }
                      else
                          return t1.chr() < t2.chr();
                  }
                  else
                      return t1.p_value() < t2.p_value();
              });
    unsigned long long cur_category = 0;
    double prev_p = p_info.lower;
    bool has_warned = false, cur_warn;
    for (auto&& snp : m_existed_snps)
    {
        snp.set_category(cur_category, prev_p, p_info.upper, p_info.inter,
                         cur_warn);
        if (cur_warn && !has_warned)
        {
            has_warned = true;
            m_reporter->report(
                "Warning: P-value threshold steps are too small. While we "
                "can still try and generate appropriate thresholds, it is "
                "likely to suffer from numeric instability, i.e. the "
                "number representing the p-value threshold of the current "
                "SNP will only have precision up to 15 digits");
        }
    }
}
bool Genotype::prepare_prsice()
{
    if (m_existed_snps.size() == 0) return false;
    if (m_very_small_thresholds)
    {
        // simply run it SNP by SNP
        std::sort(begin(m_existed_snps), end(m_existed_snps),
                  [](SNP const& t1, SNP const& t2) {
                      if (misc::logically_equal(t1.p_value(), t2.p_value()))
                      {
                          if (t1.get_file_idx() == t2.get_file_idx())
                          { return t1.get_byte_pos() < t2.get_byte_pos(); }
                          else
                              return t1.get_file_idx() < t2.get_file_idx();
                      }
                      else
                          return t1.p_value() < t2.p_value();
                  });
        unsigned long long idx = 0;
        for (auto&& snp : m_existed_snps)
        {
            snp.set_category(idx, snp.p_value());
            ++idx;
        }
    }
    else
    {
        std::sort(begin(m_existed_snps), end(m_existed_snps),
                  [](SNP const& t1, SNP const& t2) {
                      if (t1.category() == t2.category())
                      {
                          if (t1.get_file_idx() == t2.get_file_idx())
                          { return t1.get_byte_pos() < t2.get_byte_pos(); }
                          else
                              return t1.get_file_idx() < t2.get_file_idx();
                      }
                      else
                          return t1.category() < t2.category();
                  });
    }
    return true;
}
void Genotype::parse_chr_id_formula(const std::string& chr_id_formula)
{
    if (chr_id_formula.empty()) return;
    m_has_chr_id_formula = true;
    auto length = chr_id_formula.length();
    std::unordered_map<char, int> symbol_idx;
    // c = chr
    // l = loc
    // a = a1
    // b = a2
    // start at 1 so that we know and number < 0 is for symbol
    int num_symbol = 1;
    for (size_t i = 0; i < length; ++i)
    {
        auto cur_char = chr_id_formula.at(i);
        switch (cur_char)
        {
        case 'c':
        case 'C': m_chr_id_column.push_back(+BASE_INDEX::CHR); break;
        case 'l':
        case 'L': m_chr_id_column.push_back(+BASE_INDEX::BP); break;
        case 'a':
        case 'A': m_chr_id_column.push_back(+BASE_INDEX::EFFECT); break;
        case 'b':
        case 'B': m_chr_id_column.push_back(+BASE_INDEX::NONEFFECT); break;
        default:
        {
            if (symbol_idx.find(cur_char) == symbol_idx.end())
            {
                symbol_idx[cur_char] = -num_symbol;
                m_chr_id_symbol.push_back(cur_char);
                ++num_symbol;
            }
            m_chr_id_column.push_back(symbol_idx[cur_char]);
        }
        }
    }
}
// TODO: This function is likely to have bug (2019-12-09 It does)
std::vector<std::vector<size_t>>
Genotype::build_membership_matrix(const size_t num_sets,
                                  const std::vector<std::string>& region_name,
                                  const bool print_snps, std::ostream& out)
{
    // set_thresholds contain the thresholds in each set.
    if (print_snps && !out.good())
    { throw std::runtime_error("Error: Cannot open snp file to write"); }
    if (num_sets != region_name.size())
    {
        throw std::runtime_error("Error: Number of set(s) does not match the "
                                 "number of region name!");
    }
    // Structure = std::vector<std::set<double>>
    m_set_thresholds.resize(num_sets);

    const bool is_prset = (num_sets != 2);
    if (print_snps)
    {
        out << "CHR\tSNP\tBP\tP";
        for (size_t i = 0; i < region_name.size() - !is_prset; ++i)
        { out << "\t" << region_name[i]; }
        out << "\n";
    }
    bool has_snp = false;
    // temporary storage is a 2D vector. For each Set, what are the SNP idx in
    // this set
    std::vector<std::vector<size_t>> region_membership(num_sets);
    for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp)
    {
        auto&& snp = m_existed_snps[i_snp];
        auto&& flags = snp.get_flag();

        if (print_snps)
        {
            out << snp.chr() << "\t" << snp.rs() << "\t" << snp.loc() << "\t"
                << snp.p_value();
        }
        for (size_t s = 0; s < num_sets; ++s)
        {
            if (print_snps && (is_prset || s != 1))
            { out << "\t" << IS_SET(flags.data(), s); }
            if (IS_SET(flags.data(), s))
            {
                m_set_thresholds[s].insert(snp.get_threshold());
                has_snp = true;
                region_membership[s].push_back(i_snp);
            }
        }
        if (print_snps) out << "\n";
    }
    if (!has_snp)
    {
        throw std::runtime_error(
            "Error: None of the gene sets contain any SNP(s) after "
            "clumping. Have you provided the correct input? E.g. GMT "
            "file "
            "containing Entrez ID with GTF files that uses the Ensembl "
            "gene ID?\n");
    }
    return region_membership;
}

void Genotype::standardize_prs()
{
    misc::RunningStat rs;
    const size_t num_prs = m_prs_info.size();
    for (size_t i = 0; i < num_prs; ++i)
    {
        // only standardize using samples that are selected and have valid pheno
        if (!IS_SET(m_calculate_prs, i) || !m_sample_id[i].in_regression
            || IS_SET(m_exclude_from_std, i))
            continue;
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

void Genotype::get_null_score(std::vector<PRS>& prs_list,
                              const size_t& set_size, const size_t& prev_size,
                              std::vector<size_t>& background_list,
                              const bool first_run)
{

    if (m_existed_snps.empty() || set_size >= m_existed_snps.size()) return;
    // we will initailize a selected_snp_index containing the index of SNPs
    // that we'd like to add / assign to our PRS in the current round. we
    // will get anything from (prev_size , set_size]
    assert(prev_size < set_size);
    std::vector<size_t>::iterator select_start = background_list.begin();
    std::advance(select_start, static_cast<long>(prev_size));
    std::vector<size_t>::iterator select_end = background_list.begin();
    std::advance(select_end, static_cast<long>(set_size));
    std::sort(select_start, select_end);
    read_score(prs_list, select_start, select_end, first_run);
    if (m_prs_calculation.scoring_method == SCORING::STANDARDIZE
        || m_prs_calculation.scoring_method == SCORING::CONTROL_STD)
    { standardize_prs(); }
}
void Genotype::load_genotype_to_memory()
{
    // don't reserve memory if we don't need to run hard coding
    if (!m_hard_coded) { return; }
    m_genotype_stored = true;
    // this is use for initialize the array sizes
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    std::streampos cur_line;
    m_genotype_pool =
        GenotypePool(m_existed_snps.size(), unfiltered_sample_ctv2);
    std::sort(begin(m_existed_snps), end(m_existed_snps),
              [](SNP const& t1, SNP const& t2) {
                  if (t1.get_file_idx() == t2.get_file_idx())
                  { return t1.get_byte_pos() < t2.get_byte_pos(); }
                  else
                      return t1.get_file_idx() < t2.get_file_idx();
              });
    for (auto&& snp : m_existed_snps)
    {
        snp.set_genotype_storage(m_genotype_pool.alloc());
        this->count_and_read_genotype(snp);
    }
}


bool Genotype::get_score(std::vector<size_t>::const_iterator& start_index,
                         const std::vector<size_t>::const_iterator& end_index,
                         double& cur_threshold, uint32_t& num_snp_included,
                         const bool first_run)
{
    // if there are no SNPs or we are at the end
    if (m_existed_snps.size() == 0 || start_index == end_index
        || (*start_index) == m_existed_snps.size())
        return false;
    // reset number of SNPs if we don't need cumulative PRS
    if (m_prs_calculation.non_cumulate) num_snp_included = 0;
    std::vector<size_t>::const_iterator region_end = start_index;
    if (!m_very_small_thresholds)
    {
        unsigned long long cur_category =
            m_existed_snps[(*start_index)].category();
        cur_threshold = m_existed_snps[(*start_index)].get_threshold();
        for (; region_end != end_index; ++region_end)
        {
            if (m_existed_snps[(*region_end)].category() != cur_category)
            { break; }
            ++num_snp_included;
        }
    }
    else
    {
        // when we have very small thresholds, we use the p-value as the
        // indicator
        auto cur_pvalue = m_existed_snps[(*start_index)].p_value();
        cur_threshold = cur_pvalue;
        for (; region_end != end_index; ++region_end)
        {
            if (!misc::logically_equal(m_existed_snps[(*region_end)].p_value(),
                                       cur_pvalue))
            { break; }
            ++num_snp_included;
        }
    }
    read_score(start_index, region_end,
               (m_prs_calculation.non_cumulate || first_run));
    // update the current index
    start_index = region_end;
    // if ((*start_index) == 0) return -1;
    if (m_prs_calculation.scoring_method == SCORING::STANDARDIZE
        || m_prs_calculation.scoring_method == SCORING::CONTROL_STD)
    { standardize_prs(); }
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
    double known11 = static_cast<double>(2 * counts[0] + counts[1] + counts[3]);
    double known12 = static_cast<double>(2 * counts[2] + counts[1] + counts[5]);
    double known21 = static_cast<double>(2 * counts[6] + counts[3] + counts[7]);
    double known22 = static_cast<double>(2 * counts[8] + counts[5] + counts[7]);
    if (is_x1 || is_x2)
    {
        if (is_x1 && is_x2)
        {
            known11 -= static_cast<double>(static_cast<int32_t>(counts[9]));
            known12 -= static_cast<double>(static_cast<int32_t>(counts[11]));
            known21 -= static_cast<double>(static_cast<int32_t>(counts[15]));
            known22 -= static_cast<double>(static_cast<int32_t>(counts[17]));
        }
        else if (is_x1)
        {
            known11 -= (static_cast<double>(2 * counts[9] + counts[10]))
                       * (1.0 - SQRT_HALF);
            known12 -= (static_cast<double>(2 * counts[11] + counts[10]))
                       * (1.0 - SQRT_HALF);
            known21 -= (static_cast<double>(2 * counts[15] + counts[16]))
                       * (1.0 - SQRT_HALF);
            known22 -= (static_cast<double>(2 * counts[17] + counts[16]))
                       * (1.0 - SQRT_HALF);
        }
        else
        {
            known11 -= (static_cast<double>(2 * counts[9] + counts[12]))
                       * (1.0 - SQRT_HALF);
            known12 -= (static_cast<double>(2 * counts[11] + counts[12]))
                       * (1.0 - SQRT_HALF);
            known21 -= (static_cast<double>(2 * counts[15] + counts[14]))
                       * (1.0 - SQRT_HALF);
            known22 -= (static_cast<double>(2 * counts[17] + counts[14]))
                       * (1.0 - SQRT_HALF);
        }
    }
    return em_phase_hethet(known11, known12, known21, known22, counts[4],
                           freq1x_ptr, freq2x_ptr, freqx1_ptr, freqx2_ptr,
                           freq11_ptr, nullptr);
}
