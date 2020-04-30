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
        if (duplicated_snp.find(snp.rs()) == duplicated_snp.end())
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
        else if (cur_snp.loc() - prev_loc > clump_distance)
        {
            // now the chromosome didn't change, and the distance of our current
            // SNP is further away from the previous SNP than our required
            // threshold
            while (cur_snp.loc() - prev_loc > clump_distance
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
                   > clump_distance))
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

void Genotype::add_flags(const std::vector<IITree<size_t, size_t>>& gene_sets,
                         const size_t num_sets,
                         const bool genome_wide_background)
{
    const size_t num_snps = m_existed_snps.size();
    const size_t required_size = BITCT_TO_WORDCT(num_sets);
    size_t chr, bp;
    std::vector<uintptr_t> flag(required_size, 0);
    for (size_t i = 0; i < num_snps; ++i)
    {
        auto&& snp = m_existed_snps[i];
        chr = snp.chr();
        bp = snp.loc();
        construct_flag(snp.rs(), gene_sets, flag, required_size, chr, bp,
                       genome_wide_background);
        m_existed_snps[i].set_flag(num_sets, flag);
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
                         rs_id))
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
bool Genotype::process_snp(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const std::string& mismatch_snp_record_name,
    const std::string& mismatch_source, const std::string& snpid, SNP& snp,
    std::unordered_set<std::string>& processed_snps,
    std::unordered_set<std::string>& duplicated_snps,
    std::vector<bool>& retain_snp, Genotype* genotype)
{
    // TODO: allow  chr_id in future
    // std::string chr_id = SNP::chr_id(chr_num, bp);
    if (!check_rs(snpid, "", snp.rs(), processed_snps, duplicated_snps,
                  genotype))
        return false;
    auto snp_idx = genotype->m_existed_snps_index[snp.rs()];
    auto&& target_snp = genotype->m_existed_snps[snp_idx];

    misc::to_upper(snp.ref());
    misc::to_upper(snp.alt());
    bool flipping;
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
    const std::string fid = (m_ignore_fid) ? "" : token[fid_idx];
    const std::string id =
        (m_ignore_fid) ? token[iid_idx] : fid + m_delim + token[iid_idx];
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
            Sample_ID(fid, token[iid_idx], pheno, in_regression));
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
    m_sample_selection_list.clear();
}

bool Genotype::perform_freqs_and_inter(const QCFiltering& filter_info,
                                       const std::string& prefix,
                                       Genotype* target)
{
    if (!m_intermediate
        && (misc::logically_equal(filter_info.geno, 1.0)
            || filter_info.geno > 1.0)
        && (misc::logically_equal(filter_info.maf, 0.0)
            || filter_info.maf > 0.0)
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

void Genotype::efficient_clumping(const Clumping& clump_info,
                                  Genotype& reference)
{
    m_reporter->report("Start performing clumping");
    const double min_r2 = clump_info.use_proxy
                              ? std::min(clump_info.proxy, clump_info.r2)
                              : clump_info.r2;
    std::vector<bool> remain_snps(m_existed_snps.size(), false);
    // intermediates
    const uint32_t founder_ctv3 =
        BITCT_TO_ALIGNED_WORDCT(static_cast<uint32_t>(reference.m_founder_ct));
    const uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(reference.m_founder_ct);
    const uint32_t founder_ctsplit = 3 * founder_ctv3;
    const uintptr_t founder_ctv2 =
        QUATERCT_TO_ALIGNED_WORDCT(reference.m_founder_ct);
    std::vector<uintptr_t> index_data(3 * founder_ctsplit + founder_ctv3);
    std::vector<uintptr_t> index_tots(6);
    std::vector<uintptr_t> founder_include2(founder_ctv2, 0);
    fill_quatervec_55(static_cast<uint32_t>(reference.m_founder_ct),
                      founder_include2.data());
    // available memory size
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(reference.m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    std::vector<uintptr_t> geno_storage;
    try
    {
        geno_storage.resize(static_cast<size_t>((m_max_window_size + 1)
                                                * unfiltered_sample_ctv2));
        std::fill(geno_storage.begin(), geno_storage.end(), 0);
    }
    catch (...)
    {
        throw std::runtime_error("Error: Not enough memory for clumping");
    };
    double prev_progress = -1.0;
    const auto num_snp = m_existed_snps.size();
    double r2 = -1.0;
    uintptr_t* clump_geno_idx = nullptr;
    size_t num_core_snps = 0;
    for (size_t i_snp = 0; i_snp < num_snp; ++i_snp)
    {
        double progress =
            static_cast<double>(i_snp) / static_cast<double>(num_snp) * 100;
        if (progress - prev_progress > 0.01)
        {
            fprintf(stderr, "\rClumping Progress: %03.2f%%", progress);
            prev_progress = progress;
        }
        auto&& core_snp_idx = m_sort_by_p_index[i_snp];
        auto&& core_snp = m_existed_snps[core_snp_idx];
        if (core_snp.clumped() || core_snp.p_value() > clump_info.pvalue)
        { continue; }
        const size_t clump_start_idx = core_snp.low_bound();
        const size_t clump_end_idx = core_snp.up_bound();
        // we are reusing the genotype storage. So the idx starts at 1
        // the reason this is a two part process is so that we can reduce the
        // number of fseek
        clump_geno_idx = geno_storage.data();
        for (size_t clump_idx = clump_start_idx; clump_idx < core_snp_idx;
             ++clump_idx)
        {
            auto&& clump_snp = m_existed_snps[clump_idx];
            if (clump_snp.clumped() || clump_snp.p_value() > clump_info.pvalue)
            { continue; }
            reference.read_genotype(clump_geno_idx,
                                    clump_snp.get_byte_pos(true),
                                    clump_snp.get_file_idx(true));
            //++clump_geno_idx;
            clump_geno_idx = &(clump_geno_idx[unfiltered_sample_ctv2]);
        }
        // now read in core snp
        reference.read_genotype(clump_geno_idx, core_snp.get_byte_pos(true),
                                core_snp.get_file_idx(true));

        update_index_tot(founder_ctl2, founder_ctv2, reference.m_founder_ct,
                         index_data, index_tots, founder_include2,
                         clump_geno_idx);

        clump_geno_idx = geno_storage.data();
        for (size_t clump_idx = clump_start_idx; clump_idx < core_snp_idx;
             ++clump_idx)
        {
            auto&& clump_snp = m_existed_snps[clump_idx];
            if (clump_snp.clumped() || clump_snp.p_value() > clump_info.pvalue)
                // Again, ignore unwanted SNP
                continue;
            r2 = get_r2(founder_ctl2, founder_ctv2, clump_geno_idx, index_data,
                        index_tots);
            if (r2 >= min_r2)
            {
                // if the R2 between two SNP is higher than the minim threshold,
                // we will perform clumping
                // use the core SNP to clump the pair_target_snp
                core_snp.clump(clump_snp, r2, clump_info.use_proxy,
                               clump_info.proxy);
            }
            // travel to the next snp
            clump_geno_idx = &(clump_geno_idx[unfiltered_sample_ctv2]);
        }
        // now we can read the SNPs that come after the index SNP in the file

        for (size_t clump_idx = core_snp_idx + 1; clump_idx < clump_end_idx;
             ++clump_idx)
        {
            clump_geno_idx = geno_storage.data();
            auto&& clump_snp = m_existed_snps[clump_idx];
            if (clump_snp.clumped() || clump_snp.p_value() > clump_info.pvalue)
                continue;
            // read in the genotype information
            reference.read_genotype(clump_geno_idx,
                                    clump_snp.get_byte_pos(true),
                                    clump_snp.get_file_idx(true));
            r2 = get_r2(founder_ctl2, founder_ctv2, clump_geno_idx, index_data,
                        index_tots);
            // now perform clumping if required
            if (r2 >= min_r2)
            {
                core_snp.clump(clump_snp, r2, clump_info.use_proxy,
                               clump_info.proxy);
            }
        }
        core_snp.set_clumped();
        // we set the remain_core to true so that we will keep it at the end
        remain_snps[core_snp_idx] = true;
        ++num_core_snps;
    }
    fprintf(stderr, "\rClumping Progress: %03.2f%%\n\n", 100.0);
    if (num_core_snps != m_existed_snps.size())
    { shrink_snp_vector(remain_snps); }

    // we no longer require the index. might as well clear it (and hope it will
    // release the memory)
    m_existed_snps_index.clear();
    m_reporter->report("Number of variant(s) after clumping : "
                       + misc::to_string(m_existed_snps.size()));
}


void Genotype::plink_clumping(const Clumping& clump_info, Genotype& reference)
{
    // the m_existed_snp must be sorted before coming into this equation
    m_reporter->report("Start performing clumping");
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
    const double min_r2 = (clump_info.use_proxy)
                              ? std::min(clump_info.proxy, clump_info.r2)
                              : clump_info.r2;

    // when we read in the genotype, the genotype is stored in byte represented
    // by uintptr_t. Then we will use the PLINK 2 functions to obtain the R2
    // As this vector can be big, we will only initialize it once (memory
    // allocation can be expensive)
    std::vector<uintptr_t> genotype_vector(unfiltered_sample_ctl * 2);
    // We need a vector to indicate which SNPs are remaining after clumping to
    // remove any clumped SNPs. This is achieved by using this boolean vector
    std::vector<bool> remain_core(m_existed_snps.size(), false);
    // next few parameters are used to store the intermediate output from PLINK
    // function, corresponding to the count of different genotypes

    // and this is the storage to result R2
    double r2 = -1.0;
    // The following two vectors are used for storing the intermediate output.
    // Again, put memory allocation at the beginning
    std::vector<uintptr_t> index_data(3 * founder_ctsplit + founder_ctv3);
    std::vector<uintptr_t> index_tots(6);
    // This is a data mask used by PLINK in the calculation of R2. Preallocate
    // to speed up
    std::vector<uintptr_t> founder_include2(founder_ctv2, 0);
    fill_quatervec_55(static_cast<uint32_t>(reference.m_founder_ct),
                      founder_include2.data());
    // one way to speed things up as in PLINK 2 is to pre-allocate the memory
    // space for what we need to do next. The following code did precisely that
    // (borrow from PLINK2)

    intptr_t malloc_size_mb = cal_avail_memory(founder_ctv2);
    // now allocate the memory into a pointer

    // window data is the pointer walking through the allocated memory
    size_t max_window_size, num_core_snps = 0;
    unsigned char* bigstack_ua = nullptr; // ua = unaligned
    unsigned char* bigstack_initial_base;
    bigstack_ua = reinterpret_cast<unsigned char*>(malloc(
        static_cast<uintptr_t>(malloc_size_mb) * 1048576 * sizeof(char)));
    // if fail, return nullptr which will then get into the while loop
    if (!bigstack_ua)
    { throw std::runtime_error("Failed to allocate required memory"); }
    else
    {
        m_reporter->report("Allocated " + misc::to_string(malloc_size_mb)
                           + " MB successfully");
    }
    // force 64-byte align to make cache line sensitivity work (from PLINK, not
    // familiar with computer programming to know this...)
    bigstack_initial_base = reinterpret_cast<unsigned char*>(
        round_up_pow2(reinterpret_cast<uintptr_t>(bigstack_ua), CACHELINE));

    // window data is the pointer walking through the allocated memory
    unsigned char* g_bigstack_end = &(
        bigstack_initial_base[(static_cast<uintptr_t>(malloc_size_mb) * 1048576
                               - static_cast<uintptr_t>(bigstack_initial_base
                                                        - bigstack_ua))
                              & (~(CACHELINE - ONELU))]);

    // and max_window_size is the number of windows we can handle in one round
    // given the memory that we have. 1 window = 1 SNP
    max_window_size = ((reinterpret_cast<uintptr_t>(g_bigstack_end))
                       - (reinterpret_cast<uintptr_t>(bigstack_initial_base)))
                      / (founder_ctv2 * sizeof(intptr_t));

    g_bigstack_end = nullptr;
    uintptr_t* window_data =
        reinterpret_cast<uintptr_t*>(bigstack_initial_base);
    if (!max_window_size)
    { throw std::runtime_error("Error: Not enough memory for clumping!"); }
    uintptr_t* window_data_ptr = nullptr;

    // a counter to count how many windows have we read
    uintptr_t cur_window_size = 0;
    // if max_window size is 0, it means we don't have enough memory for
    // clumping

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
        if (cur_target_snp.clumped()
            || cur_target_snp.p_value() > clump_info.pvalue)
        {
            // ignore any SNP that are clumped or that has a p-value higher than
            // the clump-p threshold
            continue;
        }
        // Any SNP with p-value less than clump-p will be ignored
        // because they can never be an index SNP and thus are not of our
        // interested

        // this is the first SNP we should read from
        const size_t start = cur_target_snp.low_bound();
        // this is the first SNP we should ignore
        const size_t end = cur_target_snp.up_bound();
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
                || pair_target_snp.p_value() > clump_info.pvalue)
            {
                // ignore SNP that are clumped or that has higher p-value than
                // threshold
                continue;
            }
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
                                    pair_target_snp.get_byte_pos(true),
                                    pair_target_snp.get_file_idx(true));
            // we then move the pointer forward to the next space in the memory
            window_data_ptr = &(window_data_ptr[founder_ctv2]);
        }

        if (++cur_window_size == max_window_size)
        { throw std::runtime_error("Error: Out of memory!"); }
        // now we want to read in the index / core SNP
        // reset the content of the pointer again
        window_data_ptr[founder_ctv2 - 2] = 0;
        window_data_ptr[founder_ctv2 - 1] = 0;
        // then we can read in the genotype from the reference panel
        // note the use of cur_target_snp
        reference.read_genotype(window_data_ptr,
                                cur_target_snp.get_byte_pos(true),
                                cur_target_snp.get_file_idx(true));
        // reset the index_data information
        std::fill(index_data.begin(), index_data.end(), 0);
        // generate the required data mask
        // Disclaimer: For the next few lines, they are from PLINK and I don't
        // fully understand what they are doing
        // then populate the index_tots
        update_index_tot(founder_ctl2, founder_ctv2, reference.m_founder_ct,
                         index_data, index_tots, founder_include2,
                         window_data_ptr);
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
                || pair_target_snp.p_value() > clump_info.pvalue)
                // Again, ignore unwanted SNP
                continue;
            r2 = get_r2(founder_ctl2, founder_ctv2, window_data_ptr, index_data,
                        index_tots);
            if (r2 >= min_r2)
            {
                // if the R2 between two SNP is higher than the minim threshold,
                // we will perform clumping
                // use the core SNP to clump the pair_target_snp
                cur_target_snp.clump(pair_target_snp, r2, clump_info.use_proxy,
                                     clump_info.proxy);
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
                || pair_target_snp.p_value() > clump_info.pvalue)
                // skip if not required
                continue;
            // reset data
            window_data_ptr[founder_ctv2 - 2] = 0;
            window_data_ptr[founder_ctv2 - 1] = 0;
            // read in the genotype information
            reference.read_genotype(window_data_ptr,
                                    pair_target_snp.get_byte_pos(true),
                                    pair_target_snp.get_file_idx(true));
            r2 = get_r2(founder_ctl2, founder_ctv2, window_data_ptr, index_data,
                        index_tots);
            // now perform clumping if required
            if (r2 >= min_r2)
            {
                cur_target_snp.clump(pair_target_snp, r2, clump_info.use_proxy,
                                     clump_info.proxy);
            }
        }
        // we set the core SNP to be "clumped" so that it will no longer be
        // considered by other SNP
        cur_target_snp.set_clumped();
        // we set the remain_core to true so that we will keep it at the end
        remain_core[cur_snp_index] = true;
        ++num_core_snps;
        // we also get the p-value threshold of the SNP so that we know what
        // thresholds are included in our analysis. This information is vital
        // for the generation of all score file
    }
    fprintf(stderr, "\rClumping Progress: %03.2f%%\n\n", 100.0);
    // now we release the memory stack
    free(bigstack_ua);
    window_data = nullptr;
    window_data_ptr = nullptr;
    bigstack_initial_base = nullptr;
    bigstack_ua = nullptr;
    /*
    if (num_core_snps != m_existed_snps.size())
    { shrink_snp_vector(remain_core); }
    */
    // we no longer require the index. might as well clear it (and hope it will
    // release the memory)
    m_existed_snps_index.clear();
    m_reporter->report("Number of variant(s) after clumping : "
                       + misc::to_string(m_existed_snps.size()));
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
bool Genotype::prepare_prsice(const PThresholding& p_info)
{
    if (m_existed_snps.size() == 0) return false;
    if (m_very_small_thresholds) { recalculate_categories(p_info); }
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
    return true;
}

// TODO: This function is likely to have bug (2019-12-09 It does)
void Genotype::build_membership_matrix(
    std::vector<std::vector<size_t>>& region_membership, const size_t num_sets,
    const std::string& out, const std::vector<std::string>& region_name,
    const bool print_snps)
{
    // set_thresholds contain the thresholds in each set.
    // Structure = std::vector<std::set<double>>
    m_set_thresholds.resize(num_sets);
    std::vector<size_t> idx;
    std::ofstream snp_out;
    const std::string snp_name = out + ".snp";
    const bool is_prset = (num_sets != 2);
    if (print_snps)
    {
        snp_out.open(snp_name.c_str());
        if (!snp_out.is_open())
        {
            throw std::runtime_error("Error: Cannot open file: " + snp_name
                                     + " to write!\n");
        }
        snp_out << "CHR\tSNP\tBP\tP";
        for (size_t i = 0; i < region_name.size() - !is_prset; ++i)
        { snp_out << "\t" << region_name[i]; }
        snp_out << "\n";
    }
    bool has_snp = false;
    std::vector<bool> membership(num_sets);
    // temporary storage is a 2D vector. For each Set, what are the SNP idx in
    // this set
    region_membership.resize(num_sets);
    for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp)
    {
        auto&& snp = m_existed_snps[i_snp];
        std::fill(membership.begin(), membership.end(), false);
        for (size_t s = 0; s < num_sets; ++s)
        {
            if (snp.in(s))
            {
                m_set_thresholds[s].insert(snp.get_threshold());
                membership[s] = true;
                has_snp = true;
                region_membership[s].push_back(i_snp);
            }
        }
        if (print_snps)
        {
            snp_out << snp.chr() << "\t" << snp.rs() << "\t" << snp.loc()
                    << "\t" << snp.p_value();
            for (auto m : membership) { snp_out << "\t" << m; }
            snp_out << "\n";
        }
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
}
void Genotype::standardize_prs()
{
    misc::RunningStat rs;
    const size_t num_prs = m_prs_info.size();
    for (size_t i = 0; i < num_prs; ++i)
    {
        if (!IS_SET(m_calculate_prs, i) || IS_SET(m_exclude_from_std, i))
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
    // first, sort the SNPs for easy reading
    std::sort(begin(m_existed_snps), end(m_existed_snps),
              [](SNP const& t1, SNP const& t2) {
                  if (t1.get_file_idx() == t2.get_file_idx())
                  { return t1.get_byte_pos() < t2.get_byte_pos(); }
                  else
                      return t1.get_file_idx() < t2.get_file_idx();
              });
    // now iterate and read each SNP one by one, get the counts too
    std::vector<size_t> idx(m_existed_snps.size());
    std::iota(std::begin(idx), std::end(idx), 0);
    read_score(idx.begin(), idx.end(), true, true);
    m_genotype_stored = true;
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
    unsigned long long cur_category = m_existed_snps[(*start_index)].category();
    cur_threshold = m_existed_snps[(*start_index)].get_threshold();
    std::vector<size_t>::const_iterator region_end = start_index;
    for (; region_end != end_index; ++region_end)
    {
        if (m_existed_snps[(*region_end)].category() != cur_category)
        {
            cur_category = m_existed_snps[(*region_end)].category();
            break;
        }
        ++num_snp_included;
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
