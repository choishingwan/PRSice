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

#include "binarygen.hpp"

BinaryGen::BinaryGen(const GenoFile& geno, const Phenotype& pheno,
                     const std::string& delim, Reporter* reporter)
{
    m_sample_file = "";
    m_hard_coded = geno.hard_coded;
    const std::string message =
        initialize(geno, pheno, delim, "bgen", reporter);
    if (m_sample_file.empty() && pheno.pheno_file.empty())
    {
        throw std::runtime_error("Error: You must provide a phenotype "
                                 "file for bgen format!\n");
    }
    else if (m_sample_file.empty())
    {
        m_sample_file = pheno.pheno_file;
    }
    // initialize the context size
    m_context_map.resize(m_genotype_file_names.size());
    m_reporter->report(message);
}

size_t BinaryGen::get_sex_col(const std::string& header,
                              const std::string& format_line)
{
    std::vector<std::string> header_names = misc::split(header);
    size_t sex_col = ~size_t(0);
    for (size_t i = 3; i < header_names.size(); ++i)
    {
        misc::to_upper(header_names[i]);
        if (header_names[i].compare("SEX") == 0)
        {
            sex_col = i;
            break;
        }
    }
    // now we read in the second line
    if (sex_col != ~size_t(0))
    {
        // double check if the format is alright
        header_names = misc::split(format_line);
        if (header_names[sex_col] != "D")
        {
            m_reporter->report("Warning: Sex must be coded as "
                               "\"D\" in bgen sample file!\n"
                               "We will ignore the sex information.");
            sex_col = ~size_t(0);
        }
    }
    return sex_col;
}

void BinaryGen::handle_pheno_header(std::unique_ptr<std::istream>& sample)
{
    std::string line;
    bool have_header = false;
    size_t num_line = 0;
    while (std::getline(*sample, line))
    {
        misc::trim(line);
        if (!line.empty()) ++num_line;
    }
    if (num_line == m_unfiltered_sample_ct + 1) { have_header = true; }
    else if (num_line != m_unfiltered_sample_ct)
    {
        throw std::runtime_error(
            "Error: Number of sample in phenotype file does not match "
            "number of samples specified in bgen file. Please check "
            "you "
            "have the correct phenotype file input. Note: Phenotype "
            "file "
            "should have the same number of samples as the bgen file "
            "and "
            "they should appear in the same order");
    }
    (*sample).clear();
    (*sample).seekg(0);
    if (have_header)
    {
        std::getline(*sample, line);
        m_reporter->report("Assume phenotype file has header line: " + line);
    }
}
std::vector<Sample_ID> BinaryGen::gen_sample_vector()
{
    // this is the first time we do something w.r.t bgen file
    // first initialize the context map
    for (size_t i = 0; i < m_genotype_file_names.size(); ++i)
    { m_context_map[i] = get_context(i); }
    std::vector<Sample_ID> sample_name;
    // we always know the sample size from context
    m_unfiltered_sample_ct = m_context_map[0].number_of_samples;
    init_sample_vectors();
    if (m_is_ref && m_sample_selection_list.empty())
    {
        for (size_t i = 0; i < m_unfiltered_sample_ct; ++i)
        {
            ++m_sample_ct;
            SET_BIT(i, m_calculate_prs.data());
            // we assume all bgen samples to be founder
            SET_BIT(i, m_sample_for_ld.data());
        }
    }
    else if (!m_is_ref || !m_sample_selection_list.empty())
    {
        // this is the target, where the m_sample_file must be correct, or
        // this is the reference, which we asked for --keep or --remove and
        // an external sample file was provided (that's why we don't get
        // into the runtime_error)
        if (m_is_ref && m_sample_file.empty())
        {
            throw std::runtime_error("Error: Cannot perform sample "
                                     "filtering on the LD reference "
                                     "file without the sample file!");
        }
        auto sample = misc::load_stream(m_sample_file);
        const bool is_sample_format = check_is_sample_format(sample);
        std::string line;
        size_t sex_col = ~size_t(0);
        // now check if there's a sex information
        if (is_sample_format)
        {
            // only do this if the file is sample format
            std::getline(*sample, line);
            std::string format;
            std::getline(*sample, format);
            sex_col = get_sex_col(line, format);
        }
        // now start reading the file
        size_t line_id = 0;
        std::unordered_set<std::string> sample_in_file;
        std::vector<std::string> duplicated_sample_id;
        std::vector<std::string> token;
        const size_t required_column =
            ((sex_col != ~size_t(0)) ? (sex_col + 1) : (1 + !m_ignore_fid));
        const size_t iid_idx = (is_sample_format || !m_ignore_fid) ? 1 : 0;
        const size_t fid_idx = 0;
        if (!is_sample_format) { handle_pheno_header(sample); }
        while (std::getline(*sample, line))
        {
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            // if it is not the sample file, check if this has a header
            // not the best way, but will do it
            if (token.size() < required_column)
            {
                throw std::runtime_error(
                    "Error: Line " + std::to_string(line_id + 1)
                    + " must have at least " + std::to_string(required_column)
                    + " columns! Number of column="
                    + misc::to_string(token.size()));
            }
            gen_sample(fid_idx, iid_idx, sex_col, ~size_t(0), ~size_t(0),
                       line_id, std::unordered_set<std::string> {}, "", token,
                       sample_name, sample_in_file, duplicated_sample_id);
            ++line_id;
        }
        if (!duplicated_sample_id.empty())
        {
            // TODO: Produce a file containing id of all valid samples
            throw std::runtime_error(
                "Error: A total of "
                + misc::to_string(duplicated_sample_id.size())
                + " duplicated samples detected! Please ensure all samples "
                  "have an "
                  "unique identifier");
        }
        sample.reset();
    }
    post_sample_read_init();
    return sample_name;
}

bool BinaryGen::check_is_sample_format(std::unique_ptr<std::istream>& input)
{
    // read the sample file
    // might want to change it according to the new sample file,
    // which only mandate the first column
    // get the first two line of input
    std::string first_line, second_line;
    std::getline(*input, first_line);
    // we must have at least 2 row for a sample file
    if (!std::getline(*input, second_line)) { return false; }
    (*input).clear();
    (*input).seekg(0);
    // split the first two lines
    const std::vector<std::string_view> first_row = misc::tokenize(first_line);
    const std::vector<std::string_view> second_row =
        misc::tokenize(second_line);
    // each row should have the same number of column
    if (first_row.size() != second_row.size() || first_row.size() < 3)
    { return false; }
    // first 3 must be 0
    for (size_t i = 0; i < 3; ++i)
    {
        if (second_row[i] != "0") return false;
    }
    // Allowed character codes are: DCPB
    for (size_t i = 4; i < second_row.size(); ++i)
    {
        if (second_row[i].length() > 1) return false;
        // then any remaining column must contain either D, C, P or B
        switch (second_row[i].at(0))
        {
        case 'D':
        case 'C':
        case 'P':
        case 'B': break;
        default: return false;
        }
    }
    return true;
}

genfile::bgen::Context BinaryGen::get_context(const size_t& idx)
{
    // get the context information for the input bgen file
    // most of these codes are copy from the bgen library
    const std::string prefix = m_genotype_file_names.at(idx);
    const std::string bgen_name = prefix + ".bgen";
    std::ifstream bgen_file(bgen_name.c_str(), std::ifstream::binary);
    if (!bgen_file.is_open())
    { throw std::runtime_error("Error: Cannot open bgen file " + bgen_name); }
    genfile::bgen::Context context;
    genfile::bgen::read_header_block(bgen_file, &context);
    return context;
}

bool BinaryGen::check_sample_consistent(const std::string& bgen_name,
                                        const genfile::bgen::Context& context)
{
    // only do this if our gen file contains the sample information
    if (context.flags & genfile::bgen::e_SampleIdentifiers)
    {
        std::ifstream bgen_file(bgen_name.c_str(), std::ifstream::binary);
        uint32_t tmp_offset;
        genfile::bgen::Context tmp_context;
        genfile::bgen::read_offset(bgen_file, &tmp_offset);
        genfile::bgen::read_header_block(bgen_file, &tmp_context);
        uint32_t sample_block_size = 0;
        uint32_t actual_number_of_samples = 0;
        uint16_t identifier_size;
        std::string identifier;
        std::size_t bytes_read = 0;
        std::unordered_set<std::string> dup_check;
        genfile::bgen::read_little_endian_integer(bgen_file,
                                                  &sample_block_size);
        genfile::bgen::read_little_endian_integer(bgen_file,
                                                  &actual_number_of_samples);
        // +8 here to account for the sample_block_size and
        // actual_number_of_samples read above. Direct Copy and paste from
        // BGEN lib function read_sample_identifier_block
        bytes_read += 8;
        if (actual_number_of_samples != context.number_of_samples)
        {
            throw std::runtime_error("Error: Number of sample from your "
                                     ".sample/ phenotype file does not match "
                                     "the number of sample included in the "
                                     "bgen file. Maybe check if you have used "
                                     "a filtered sample or phenotype file?");
        }
        // we don't need to check the sample name when we are doing the LD
        // as we don't store any sample ID
        // though might be good practice to actually check the names matched
        // between the sample name and the ld bgen file to avoid mismatch
        // That can be done later on
        if (!m_is_ref)
        {
            const bool has_fid = !m_sample_id.front().FID.empty();
            size_t sample_vector_idx = 0;
            for (size_t i = 0; i < actual_number_of_samples; ++i)
            {
                genfile::bgen::read_length_followed_by_data(
                    bgen_file, &identifier_size, &identifier);
                if (!bgen_file)
                    throw std::runtime_error(
                        "Error: Problem reading bgen file!");
                bytes_read += sizeof(identifier_size) + identifier_size;
                // Need to double check. BGEN format might differ depends
                // if FID is provided. When FID is provided, then the ID
                // should be FID + delimitor + IID; otherwise it'd be IID
                if (IS_SET(m_calculate_prs.data(), i))
                {
                    if (m_sample_id[sample_vector_idx].IID != identifier
                        && (m_sample_id[sample_vector_idx].FID + m_delim
                            + m_sample_id[sample_vector_idx].IID)
                               != identifier)
                    {
                        std::string error_message =
                            "Error: Sample mismatch "
                            "between bgen and phenotype file! Name in BGEN "
                            "file is "
                            ":"
                            + identifier + " and in phentoype file is: ";
                        if (has_fid)
                            error_message.append(
                                m_sample_id[sample_vector_idx].FID + m_delim
                                + m_sample_id[sample_vector_idx].IID);
                        else
                            error_message.append(
                                m_sample_id[sample_vector_idx].IID);
                        error_message.append(
                            ". Please note that PRSice require the bgen file "
                            "and "
                            "the .sample (or phenotype file if sample file is "
                            "not provided) to have sample in the same order. "
                            "(We "
                            "might be able to losen this requirement in future "
                            "when we have more time)");
                        throw std::runtime_error(error_message);
                    }
                    ++sample_vector_idx;
                }
            }
        }
        assert(bytes_read == sample_block_size);
    }
    return true;
}


void BinaryGen::gen_snp_vector(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const std::string& out_prefix, Genotype* target)
{
    const std::string mismatch_snp_record_name = out_prefix + ".mismatch";
    const std::string mismatch_source = m_is_ref ? "Reference" : "Base";
    std::unordered_set<std::string> duplicated_snps;
    std::unordered_set<std::string> processed_snps;
    auto&& genotype = (m_is_ref) ? target : this;
    std::vector<bool> retain_snp(genotype->m_existed_snps.size(), false);
    std::ifstream bgen_file;
    std::string bgen_name;
    std::string allele;
    std::string SNPID;
    std::string RSID;
    std::string cur_id;
    std::string chromosome;
    std::string prev_chr = "";
    std::string file_name;
    std::string error_message = "";
    std::string A1, A2, prefix;
    std::streampos byte_pos;
    size_t total_unfiltered_snps = 0;
    size_t ref_target_match = 0;
    size_t num_snp;
    size_t chr_num = 0;
    uint32_t SNP_position = 0;
    uint32_t offset;
    bool sex_error = false;
    bool chr_error = false;
    for (size_t i = 0; i < m_genotype_file_names.size(); ++i)
    {
        // get the total unfiltered snp size so that we can initalize the vector
        // TODO: Check if there are too many SNPs, therefore cause integer
        // overflow
        total_unfiltered_snps += m_context_map[i].number_of_variants;
    }
    if (!m_is_ref)
    {
        check_sample_consistent(
            std::string(m_genotype_file_names.front() + ".bgen"),
            m_context_map[0]);
    }
    for (size_t file_idx = 0; file_idx < m_genotype_file_names.size();
         ++file_idx)
    {
        // now start reading each bgen file
        prefix = m_genotype_file_names[file_idx];
        bgen_name = prefix + ".bgen";
        if (bgen_file.is_open()) bgen_file.close();
        bgen_file.clear();
        bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
        if (!bgen_file.is_open())
        {
            std::string error_message =
                "Error: Cannot open bgen file " + bgen_name;
            throw std::runtime_error(error_message);
        }
        // read in the offset
        auto&& context = m_context_map[file_idx];
        offset = context.offset;
        // skip the offest
        bgen_file.seekg(offset + 4);
        num_snp = context.number_of_variants;
        // obtain the context information. We don't check out of bound as that
        // is unlikely to happen
        for (size_t i_snp = 0; i_snp < num_snp; ++i_snp)
        {
            // go through each SNP in the file
            if (i_snp % 1000 == 0)
            {
                fprintf(stderr, "\r%zuK SNPs processed in %s   \r",
                        i_snp / 1000, bgen_name.c_str());
            }
            else if (i_snp < 1000)
            {
                fprintf(stderr, "\r%zu SNPs processed in %s\r", i_snp,
                        bgen_name.c_str());
            }
            ++m_unfiltered_marker_ct;
            // directly use the library without decompressing the genotype
            read_snp_identifying_data(bgen_file, context, &SNPID, &RSID,
                                      &chromosome, &SNP_position, &A1, &A2);
            // get the current location of bgen file, this will be used to skip
            // to current location later on
            byte_pos = bgen_file.tellg();
            if (check_chr(chromosome, prev_chr, chr_num, chr_error, sex_error))
            {
                SNP cur_snp(RSID, chr_num, SNP_position, A1, A2, file_idx,
                            byte_pos);
                if (process_snp(exclusion_regions, mismatch_snp_record_name,
                                mismatch_source, SNPID, cur_snp, processed_snps,
                                duplicated_snps, retain_snp, genotype))
                { ++ref_target_match; }
            }

            // read in the genotype data block so that we advance the ifstream
            // pointer to the next SNP entry
            read_genotype_data_block(bgen_file, context, &m_buffer1);
        }
        if (num_snp % 1000 == 0)
        {
            fprintf(stderr, "\r%zuK SNPs processed in %s   \r", num_snp / 1000,
                    bgen_name.c_str());
        }
        else if (num_snp < 1000)
        {
            fprintf(stderr, "\r%zu SNPs processed in %s\r", num_snp,
                    bgen_name.c_str());
        }
        bgen_file.close();
        fprintf(stderr, "\n");
    }
    if (ref_target_match != genotype->m_existed_snps.size())
    {
        // there are mismatch, so we need to update the snp vector
        genotype->shrink_snp_vector(retain_snp);
        // now update the SNP vector index
        genotype->update_snp_index();
    }
    if (duplicated_snps.size() != 0)
    {
        throw std::runtime_error(
            genotype->print_duplicated_snps(duplicated_snps, out_prefix));
    }
}


bool BinaryGen::calc_freq_gen_inter(const QCFiltering& filter_info,
                                    const std::string& prefix, Genotype* target,
                                    bool force_cal)
{
    if (!m_intermediate
        && (misc::logically_equal(filter_info.geno, 1.0)
            || filter_info.geno > 1.0)
        && (misc::logically_equal(filter_info.maf, 0.0)
            || filter_info.maf < 0.0)
        && (misc::logically_equal(filter_info.info_score, 0.0)
            || filter_info.maf < 0.0)
        && !force_cal)
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
    const double sample_ct_recip = 1.0 / (static_cast<double>(m_sample_ct));
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    std::vector<bool> retain_snps(genotype->m_existed_snps.size(), false);

    std::string bgen_name = "";
    const std::string intermediate_name = prefix + ".inter";
    std::ifstream bgen_file;
    double cur_maf, cur_geno;
    std::streampos byte_pos, tmp_byte_pos;
    size_t processed_count = 0;
    size_t cur_file_idx = 0;
    size_t retained = 0;
    size_t ll_ct = 0;
    size_t lh_ct = 0;
    size_t hh_ct = 0;
    size_t uii = 0;
    size_t missing = 0;
    // initialize the sample inclusion mask
    std::vector<uintptr_t> sample_include2(unfiltered_sample_ctv2);
    std::vector<uintptr_t> founder_include2(unfiltered_sample_ctv2);
    // fill it with the required mask (copy from PLINK2)
    init_quaterarr_from_bitarr(m_calculate_prs.data(), m_unfiltered_sample_ct,
                               sample_include2.data());
    init_quaterarr_from_bitarr(m_sample_for_ld.data(), m_unfiltered_sample_ct,
                               founder_include2.data());
    m_tmp_genotype.resize(unfiltered_sample_ctv2, 0);
    // we initialize the plink converter with the sample inclusion vector and
    // also the tempory genotype vector list. We also provide the hard coding
    // threshold
    PLINK_generator setter(&m_calculate_prs, m_tmp_genotype.data(),
                           m_hard_threshold, m_dose_threshold);
    // now consider if we are generating the intermediate file
    std::ofstream inter_out;
    if (m_intermediate)
    {
        // allow generation of intermediate file
        if (m_target_plink && m_is_ref)
        {
            // target already generated some intermediate, now append for
            // reference
            inter_out.open(intermediate_name.c_str(),
                           std::ios::binary | std::ios::app);
        }
        else
        {
            // a new intermediate file
            inter_out.open(intermediate_name.c_str(), std::ios::binary);
        }
    }
    // now start processing the bgen file
    double progress = 0, prev_progress = -1.0;
    const size_t total_snp = genotype->m_existed_snps.size();
    genfile::bgen::Context context;
    for (auto&& snp : genotype->m_existed_snps)
    {
        progress = static_cast<double>(processed_count)
                   / static_cast<double>(total_snp) * 100;
        if (progress - prev_progress > 0.01)
        {
            fprintf(stderr, "\rCalculating allele frequencies: %03.2f%%",
                    progress);
            prev_progress = progress;
        }
        snp.get_file_info(cur_file_idx, byte_pos, m_is_ref);
        context = m_context_map[cur_file_idx];
        ++processed_count;
        // now read in the genotype information
        genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
            m_genotype_file, m_genotype_file_names[cur_file_idx] + ".bgen",
            context, setter, &m_buffer1, &m_buffer2, byte_pos);
        // no founder, much easier
        setter.get_count(ll_ct, lh_ct, hh_ct, missing);
        uii = ll_ct + lh_ct + hh_ct;
        cur_geno = 1.0 - ((static_cast<int32_t>(uii)) * sample_ct_recip);
        uii = 2 * (ll_ct + lh_ct + hh_ct);
        // tmp_total = (ll_ctf + lh_ctf + hh_ctf);
        // assert(m_founder_ct >= tmp_total);
        // missing = m_founder_ct - tmp_total;
        if (!uii) { cur_maf = 0.5; }
        else
        {
            cur_maf = (static_cast<double>(2 * hh_ct + lh_ct))
                      / (static_cast<double>(uii));
            cur_maf = (cur_maf > 0.5) ? 1 - cur_maf : cur_maf;
        }
        // filter by genotype missingness
        if (filter_info.geno < cur_geno)
        {
            ++m_num_geno_filter;
            continue;
        }
        // filter by MAF
        // do not flip the MAF for now, so that we
        // are not confuse later on
        // remove SNP if maf lower than threshold
        if (cur_maf < filter_info.maf)
        {
            ++m_num_maf_filter;
            continue;
        }
        else if (ll_ct == m_sample_ct || hh_ct == m_sample_ct)
        {
            // none of the sample contain this SNP
            // still count as MAF filtering (for now)
            ++m_num_maf_filter;
            continue;
        }

        if (setter.info_score() < filter_info.info_score)
        {
            ++m_num_info_filter;
            continue;
        }
        // if we can reach here, it is not removed
        snp.set_counts(ll_ct, lh_ct, hh_ct, missing, m_is_ref);
        if (m_is_ref) { snp.set_ref_expected(setter.expected()); }
        else
        {
            snp.set_expected(setter.expected());
        }
        ++retained;
        // we need to -1 because we put processed_count ++ forward
        // to avoid continue skipping out the addition
        retain_snps[processed_count - 1] = true;
        if (m_intermediate
            && (m_is_ref || !m_expect_reference || (!m_is_ref && m_hard_coded)))
        {
            // we will only generate the intermediate file if
            // the following happen:
            // 1. User want to generate the intermediate file
            // 2. We are dealing with reference file format
            // 3. We are dealing with target file and there is
            // no reference file
            // 4. We are dealing with target file and we are
            // expected to use hard_coding
            tmp_byte_pos = inter_out.tellp();
            inter_out.write(reinterpret_cast<char*>(&m_tmp_genotype[0]),
                            m_tmp_genotype.size() * sizeof(m_tmp_genotype[0]));
            if (!m_is_ref)
            {
                // target file
                if (m_hard_coded)
                {
                    m_target_plink = true;
                    snp.update_file(m_genotype_file_names.size(), tmp_byte_pos,
                                    false);
                }
                if (!m_expect_reference)
                {
                    // we don't have reference
                    m_ref_plink = true;
                    snp.update_file(m_genotype_file_names.size(), tmp_byte_pos,
                                    true);
                }
            }
            else
            {
                // this is the reference file
                m_ref_plink = true;
                snp.update_file(m_genotype_file_names.size(), tmp_byte_pos,
                                true);
            }
        }
    }
    if (m_intermediate
        && (m_is_ref || (!m_is_ref && m_hard_coded) || !m_expect_reference))
    { // update our genotype file
        inter_out.close();
        m_genotype_file_names.push_back(intermediate_name);
    }
    fprintf(stderr, "\rCalculating allele frequencies: %03.2f%%\n", 100.0);
    // now update the vector
    if (retained != genotype->m_existed_snps.size())
    { genotype->shrink_snp_vector(retain_snps); }
    return true;
}

BinaryGen::~BinaryGen()
{
    if (m_target_plink || m_ref_plink)
    {
        // if we have constructed the intermediate file, we should remove it
        // to save space (plus that file isn't of any useful format and
        // can't be used by any other problem nor can it be reused)
        std::remove(m_genotype_file_names.back().c_str());
    }
}

void BinaryGen::dosage_score(
    std::vector<PRS>& prs_list,
    const std::vector<size_t>::const_iterator& start_idx,
    const std::vector<size_t>::const_iterator& end_idx, bool reset_zero)
{
    // currently, use_ref_maf doesn't work on bgen dosage file
    // main reason is we need expected value instead of
    // the MAF
    bool not_first = !reset_zero;
    // we initialize the PRS interpretor with the required information.
    // m_prs_info is where we store the PRS information
    // and m_sample_include let us know if the sample is required.
    // m_missing_score will inform us as to how to handle the missingness
    PRS_Interpreter setter(&prs_list, &m_calculate_prs,
                           m_prs_calculation.missing_score);
    std::vector<size_t>::const_iterator cur_idx = start_idx;
    size_t file_idx;
    std::streampos byte_pos;
    for (; cur_idx != end_idx; ++cur_idx)
    {
        auto&& snp = m_existed_snps[(*cur_idx)];
        snp.get_file_info(file_idx, byte_pos, m_is_ref);
        // if the file name differ, or the file isn't open, we will open it
        auto&& context = m_context_map[file_idx];
        setter.set_stat(snp.stat(), m_homcom_weight, m_het_weight,
                        m_homrar_weight, snp.is_flipped(), not_first);

        // start performing the parsing
        genfile::bgen::read_and_parse_genotype_data_block<PRS_Interpreter>(
            m_genotype_file, m_genotype_file_names[file_idx] + ".bgen", context,
            setter, &m_buffer1, &m_buffer2, byte_pos);
        // check if this SNP has some non-missing sample, if not, invalidate
        // it
        // after reading in this SNP, we no longer need to reset the PRS
        not_first = true;
    }
}


void BinaryGen::hard_code_score(
    std::vector<PRS>& prs_list,
    const std::vector<size_t>::const_iterator& start_idx,
    const std::vector<size_t>::const_iterator& end_idx, bool reset_zero,
    bool read_only)
{
    // we need to calculate the size of possible vectors
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    // genotype counts
    size_t homrar_ct = 0;
    size_t missing_ct = 0;
    size_t het_ct = 0;
    size_t homcom_ct = 0;
    // weight of each genotype
    double homcom_weight = m_homcom_weight;
    double het_weight = m_het_weight;
    double homrar_weight = m_homrar_weight;
    // this is needed if we want to calculate the MAF of the sample
    const size_t ploidy = 2;
    const size_t miss_count =
        (m_prs_calculation.missing_score != MISSING_SCORE::SET_ZERO) * ploidy;
    const bool is_centre =
        (m_prs_calculation.missing_score == MISSING_SCORE::CENTER);
    const bool mean_impute =
        (m_prs_calculation.missing_score == MISSING_SCORE::MEAN_IMPUTE);
    // check if we need to reset the sample's PRS
    bool not_first = !reset_zero;
    double stat, maf, adj_score, miss_score;
    std::streampos byte_pos;
    // initialize the data structure for storing the genotype
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    genfile::bgen::Context context;
    PLINK_generator setter(&m_calculate_prs, genotype.data(), m_hard_threshold,
                           m_dose_threshold);
    std::vector<size_t>::const_iterator cur_idx = start_idx;
    size_t idx;
    for (; cur_idx != end_idx; ++cur_idx)
    {
        auto&& cur_snp = m_existed_snps[(*cur_idx)];
        // read in the genotype using the modified load_and_collapse_incl
        // function. m_target_plink will inform the function wheter there's
        // an intermediate file
        // if it has the intermediate file, then we should have already
        // calculated the counts
        if (!cur_snp.stored_genotype())
        {
            cur_snp.get_file_info(idx, byte_pos, m_is_ref);

            auto&& file_name = m_genotype_file_names[idx];
            if (m_intermediate
                && cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct,
                                      m_prs_calculation.use_ref_maf))
            {
                // Have intermediate file and have the counts
                // read in the genotype information to the genotype vector
                m_genotype_file.read(file_name, byte_pos, unfiltered_sample_ct4,
                                     reinterpret_cast<char*>(genotype.data()));
            }
            else if (m_intermediate)
            {
                // Have intermediate file but not have the counts
                throw std::logic_error("Error: Sam has a logic error in bgen");
            }
            else
            {
                // now read in the genotype information
                context = m_context_map[idx];
                // start performing the parsing
                genfile::bgen::read_and_parse_genotype_data_block<
                    PLINK_generator>(m_genotype_file, file_name + ".bgen",
                                     context, setter, &m_buffer1, &m_buffer2,
                                     byte_pos);
                setter.get_count(homcom_ct, het_ct, homrar_ct, missing_ct);
            }
            if (read_only) cur_snp.assign_genotype(genotype);
        }
        else
        {
            genotype = cur_snp.get_genotype();
            cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct,
                               m_prs_calculation.use_ref_maf);
        }
        // TODO: if we haven't got the count from the genotype matrix, we will
        // need to calculate that, we might not need to do the
        // counting as that is already done when we convert the dosages into
        // the binary genotypes
        homcom_weight = m_homcom_weight;
        het_weight = m_het_weight;
        homrar_weight = m_homrar_weight;
        maf = 1.0
              - static_cast<double>(homcom_weight * homcom_ct
                                    + het_ct * het_weight
                                    + homrar_weight * homrar_ct)
                    / (static_cast<double>(homcom_ct + het_ct + homrar_ct)
                       * ploidy);
        if (cur_snp.is_flipped())
        {
            // change the mean to reflect flipping
            maf = 1.0 - maf;
            // swap the weighting
            std::swap(homcom_weight, homrar_weight);
        }
        stat = cur_snp.stat(); // Multiply by ploidy
        // only set these value to the imputed value if we require them
        adj_score = 0;
        if (is_centre) { adj_score = ploidy * stat * maf; }
        miss_score = 0;
        if (mean_impute)
        {
            // again, mean_impute is stable, branch prediction should be ok
            miss_score = ploidy * stat * maf;
        }
        // start reading the genotype
        if (!read_only)
        {
            read_prs(genotype, prs_list, ploidy, stat, adj_score, miss_score,
                     miss_count, homcom_weight, het_weight, homrar_weight,
                     not_first);
        }
        // we've finish processing the first SNP no longer need to reset the
        // PRS
        not_first = true;
    }
}


void BinaryGen::read_score(std::vector<PRS>& prs_list,
                           const std::vector<size_t>::const_iterator& start_idx,
                           const std::vector<size_t>::const_iterator& end_idx,
                           bool reset_zero, bool ultra)
{
    // because I don't want to touch the code in dosage_score, we will reset
    // the sample here reset_sample_prs();
    if (m_hard_coded)
    {
        // for hard coded, we need to check if intermediate file is used
        // instead
        hard_code_score(prs_list, start_idx, end_idx, reset_zero, ultra);
    }
    else
    {
        dosage_score(prs_list, start_idx, end_idx, reset_zero);
    }
}
