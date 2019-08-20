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


#include "binaryplink.hpp"


BinaryPlink::BinaryPlink(const std::string& file_list, const std::string& file,
                         size_t thread, const bool ignore_fid,
                         const bool keep_nonfounder, const bool keep_ambig,
                         const bool is_ref, Reporter* reporter)
{
    m_thread = thread;
    m_ignore_fid = ignore_fid;
    m_keep_nonfounder = keep_nonfounder;
    m_keep_ambig = keep_ambig;
    m_is_ref = is_ref;
    m_reporter = reporter;
    // set the chromosome information
    // will need to add more script here if we want to support something other
    // than human
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    m_haploid_mask.resize(CHROM_MASK_WORDS, 0);
    // main use of following function is to set the max code
    init_chr();
    std::string message = "Initializing Genotype ";
    std::string listed_input = file_list;
    std::string input = file;
    std::vector<std::string> token;
    if (!listed_input.empty())
        token = misc::split(listed_input, ",");
    else
        token = misc::split(input, ",");
    bool external_sample = (token.size() == 2);
    if (external_sample) { m_sample_file = token[1]; }

    if (!listed_input.empty())
    {
        // has listed input
        // check if there is an external sample file
        listed_input = token.front();
        m_genotype_file_names = load_genotype_prefix(listed_input);
        message.append("info from file: " + listed_input + " (bed)\n");
    }
    else
    {
        // single file input, check for # and replace it with 1-22
        input = token.front();
        m_genotype_file_names = set_genotype_files(input);
        message.append("file: " + input + " (bed)\n");
    }
    if (external_sample)
    { message.append("With external fam file: " + m_sample_file + "\n"); }
    else
    {
        m_sample_file = m_genotype_file_names.front() + ".fam";
    }
    m_reporter->report(message);
}

std::vector<Sample_ID> BinaryPlink::gen_sample_vector(const std::string& delim)
{
    assert(m_genotype_file_names.size() > 0);
    std::ifstream famfile;
    famfile.open(m_sample_file.c_str());
    if (!famfile.is_open())
    {
        std::string error_message =
            "Error: Cannot open fam file: " + m_sample_file;
        throw std::runtime_error(error_message);
    }
    // number of unfiltered samples
    // this must be correct as this value is use for all subsequent size
    // intiailization
    m_unfiltered_sample_ct = 0;

    std::string line;
    // capture all founder name and check if they exists within the file
    std::unordered_set<std::string> founder_info;
    // first pass to get the number of samples and also get the founder ID
    std::vector<std::string> token;
    while (std::getline(famfile, line))
    {
        misc::trim(line);
        if (!line.empty())
        {
            token = misc::split(line);
            if (token.size() < 6)
            {
                std::string message =
                    "Error: Malformed fam file. Less than 6 column on "
                    "line: "
                    + std::to_string(m_unfiltered_sample_ct + 1) + "\n";
                throw std::runtime_error(message);
            }
            founder_info.insert(token[+FAM::FID] + delim + token[+FAM::IID]);
            m_unfiltered_sample_ct++;
        }
    }
    // now reset the fam file to the start
    famfile.clear();
    famfile.seekg(0);
    // the unfiltered_sampel_ct is used to define the size of all vector
    // used within the program
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);

    // Currently ignore sex information
    m_founder_info.resize(unfiltered_sample_ctl, 0);
    m_sample_include.resize(unfiltered_sample_ctl, 0);

    m_num_male = 0;
    m_num_female = 0;
    m_num_ambig_sex = 0;
    m_num_non_founder = 0;
    std::vector<Sample_ID> sample_name;
    std::unordered_set<std::string> duplicated_samples;
    std::vector<std::string> duplicated_sample_id;
    uintptr_t sample_index = 0; // this is just for error message
    bool inclusion = false;
    bool founder = false;
    while (std::getline(famfile, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        // we have already checked for malformed file
        std::string id = (m_ignore_fid)
                             ? token[+FAM::IID]
                             : token[+FAM::FID] + delim + token[+FAM::IID];
        if (!m_remove_sample)
        {
            inclusion = (m_sample_selection_list.find(id)
                         != m_sample_selection_list.end());
        }
        else
        {
            inclusion = (m_sample_selection_list.find(id)
                         == m_sample_selection_list.end());
        }

        if (founder_info.find(token[+FAM::FID] + delim + token[+FAM::FATHER])
                == founder_info.end()
            && founder_info.find(token[+FAM::FID] + delim + token[+FAM::MOTHER])
                   == founder_info.end()
            && inclusion)
        {
            // this is a founder (with no dad / mum)
            ++m_founder_ct;
            // m_founder_info is a subset of m_sample_include
            // and is used for LD calculation
            SET_BIT(sample_index, m_founder_info.data());
            SET_BIT(sample_index, m_sample_include.data());
            founder = true;
        }
        else if (inclusion)
        {
            // we want to include this sample, but this is not a founder
            // As we also want to generate the PRS for any non-founders, we
            // will include them, but later in regression, stat that they should
            // not be included in the regression unless otherwise
            SET_BIT(sample_index, m_sample_include.data());
            ++m_num_non_founder;
            // if we want to keep it as found, we will store it as founder in
            // the sample ID, which is only use for PRS calculation, not for LD
            founder = m_keep_nonfounder;
        }
        m_sample_ct += inclusion;
        if (token[+FAM::SEX] == "1") { ++m_num_male; }
        else if (token[+FAM::SEX] == "2")
        {
            ++m_num_female;
        }
        else
        {
            ++m_num_ambig_sex;
        }
        // this must be incremented within each loop
        ++sample_index;
        if (duplicated_samples.find(id) != duplicated_samples.end())
            duplicated_sample_id.push_back(id);
        // only store samples that we need, and use the m_sample_include and
        // m_founder_info to indicate if sample is needed for subsequent
        // operations
        if (inclusion && !m_is_ref)
        {
            sample_name.emplace_back(
                Sample_ID(token[+FAM::FID], token[+FAM::IID],
                          token[+FAM::PHENOTYPE], founder));
        }
        duplicated_samples.insert(id);
    }

    if (!duplicated_sample_id.empty())
    {
        // TODO: Produce a file containing id of all valid samples
        std::string error_message =
            "Error: A total of " + misc::to_string(duplicated_sample_id.size())
            + " duplicated samples detected!\n";
        error_message.append(
            "Please ensure all samples have an unique identifier");
        throw std::runtime_error(error_message);
    }

    famfile.close();
    // initialize the m_tmp_genotype vector
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    m_tmp_genotype.resize(unfiltered_sample_ctv2, 0);
    // m_prs_info.reserve(m_sample_ct);
    // now we add the prs information. For some reason, we can't do a simple
    // reserve
    for (size_t i = 0; i < m_sample_ct; ++i) { m_prs_info.emplace_back(PRS()); }
    // also resize the in_regression flag
    m_in_regression.resize(m_sample_include.size(), 0);
    // initialize the sample_include2 and founder_include2 which are
    // both needed in cal_maf or read_score for MAF calculation
    m_sample_include2.resize(unfiltered_sample_ctv2);
    m_founder_include2.resize(unfiltered_sample_ctv2);
    // fill it with the required mask (copy from PLINK2)
    init_quaterarr_from_bitarr(m_sample_include.data(), m_unfiltered_sample_ct,
                               m_sample_include2.data());
    init_quaterarr_from_bitarr(m_founder_info.data(), m_unfiltered_sample_ct,
                               m_founder_include2.data());
    return sample_name;
}

void BinaryPlink::calc_freq_gen_inter(
    const double& maf_threshold, const double& geno_threshold,
    const double& /*info_threshold*/, const bool maf_filter,
    const bool geno_filter, const bool /*info_filter*/,
    const bool /*hard_coded*/, Genotype* target)
{
    // we will go through all the SNPs
    auto&& genotype = (m_is_ref) ? target : this;
    // sort SNPs by the read order to minimize skipping
    if (!m_is_ref)
    {
        // not target
        std::sort(begin(genotype->m_existed_snps),
                  end(genotype->m_existed_snps),
                  [](SNP const& t1, SNP const& t2) {
                      if (t1.file_index() == t2.file_index())
                      { return t1.byte_pos() < t2.byte_pos(); }
                      else
                          return t1.file_index() < t2.file_index();
                  });
    }
    else
    {
        // reference file
        std::sort(begin(genotype->m_existed_snps),
                  end(genotype->m_existed_snps),
                  [](SNP const& t1, SNP const& t2) {
                      if (t1.ref_file_index() == t2.ref_file_index())
                      { return t1.ref_byte_pos() < t2.ref_byte_pos(); }
                      else
                          return t1.ref_file_index() < t2.ref_file_index();
                  });
    }

    // now process the SNPs
    const double sample_ct_recip = 1.0 / (static_cast<double>(m_sample_ct));
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    const size_t total_snp = genotype->m_existed_snps.size();
    std::vector<bool> retain_snps(genotype->m_existed_snps.size(), false);
    std::ifstream bed_file;
    std::string bed_name;
    std::string prev_file = "";
    double progress = 0.0, prev_progress = -1.0;
    double cur_maf, cur_geno;
    unsigned long long byte_pos = 1;
    size_t processed_count = 0;
    size_t retained = 0;
    size_t cur_file_idx = 0;
    uint32_t ll_ct = 0;
    uint32_t lh_ct = 0;
    uint32_t hh_ct = 0;
    uint32_t ll_ctf = 0;
    uint32_t lh_ctf = 0;
    uint32_t hh_ctf = 0;
    uint32_t uii = 0;
    uint32_t missing = 0;
    uint32_t tmp_total = 0;
    // initialize the sample inclusion mask
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
        ++processed_count;
        if (m_is_ref)
        {
            cur_file_idx = snp.ref_file_index();
            byte_pos = snp.ref_byte_pos();
        }
        else
        {
            cur_file_idx = snp.file_index();
            byte_pos = snp.byte_pos();
        }

        // read in the genotype information to the genotype vector
        char* geno = reinterpret_cast<char*>(m_tmp_genotype.data());
        for (uintptr_t i = 0; i < unfiltered_sample_ct4; ++i)
        {
            *geno = m_genotype_files[cur_file_idx][byte_pos + i];
            ++geno;
        }
        // calculate the MAF using PLINK2 function
        single_marker_freqs_and_hwe(
            unfiltered_sample_ctv2, m_tmp_genotype.data(),
            m_sample_include2.data(), m_founder_include2.data(), m_sample_ct,
            &ll_ct, &lh_ct, &hh_ct, m_founder_ct, &ll_ctf, &lh_ctf, &hh_ctf);
        uii = ll_ct + lh_ct + hh_ct;
        cur_geno = 1.0 - (static_cast<int32_t>(uii)) * sample_ct_recip;
        uii = 2 * (ll_ctf + lh_ctf + hh_ctf);
        tmp_total = (ll_ctf + lh_ctf + hh_ctf);
        assert(m_founder_ct >= tmp_total);
        missing = static_cast<uint32_t>(m_founder_ct) - tmp_total;
        if (!uii) { cur_maf = 0.5; }
        else
        {
            cur_maf = (static_cast<double>(2 * hh_ctf + lh_ctf))
                      / (static_cast<double>(uii));

            cur_maf = (cur_maf > 0.5) ? 1 - cur_maf : cur_maf;
        }
        if (misc::logically_equal(cur_maf, 0.0)
            || misc::logically_equal(cur_maf, 1.0))
        {
            // none of the sample contain this SNP
            // still count as MAF filtering (for now)
            ++m_num_maf_filter;
            continue;
        }
        // filter by genotype missingness
        if (geno_filter && geno_threshold < cur_geno)
        {
            ++m_num_geno_filter;
            continue;
        }
        // filter by MAF
        // do not flip the MAF for now, so that we
        // are not confuse later on
        // remove SNP if maf lower than threshold
        if (maf_filter && cur_maf < maf_threshold)
        {
            ++m_num_maf_filter;
            continue;
        }
        // if we can reach here, it is not removed
        if (m_is_ref) { snp.set_ref_counts(ll_ctf, lh_ctf, hh_ctf, missing); }
        else
        {
            snp.set_counts(ll_ctf, lh_ctf, hh_ctf, missing);
        }
        ++retained;
        // we need to -1 because we put processed_count ++ forward
        // to avoid continue skipping out the addition
        retain_snps[processed_count - 1] = true;
    }

    fprintf(stderr, "\rCalculating allele frequencies: %03.2f%%\n", 100.0);
    // now update the vector
    if (retained != genotype->m_existed_snps.size())
    { genotype->shrink_snp_vector(retain_snps); }
}

void BinaryPlink::gen_snp_vector(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const std::string& out_prefix, Genotype* target)
{
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    std::unordered_set<std::string> processed_snps;
    std::unordered_set<std::string> duplicated_snp;
    std::vector<std::string> bim_token;
    std::vector<bool> retain_snp;
    auto&& genotype = (m_is_ref) ? target : this;
    retain_snp.resize(genotype->m_existed_snps.size(), false);
    std::ifstream bim;
    std::ofstream mismatch_snp_record;
    std::string bim_name, bed_name, chr, line;
    std::string prev_chr = "", error_message = "";
    std::string mismatch_snp_record_name = out_prefix + ".mismatch";
    std::string prefix;
    // temp holder for region match
    uintptr_t bed_offset;
    size_t num_retained = 0;
    size_t chr_num = 0;
    size_t num_snp_read = 0;
    unsigned long long byte_pos;
    int chr_code = 0;
    bool chr_error = false, chr_sex_error = false, prev_chr_sex_error = false,
         prev_chr_error = false, flipping = false, to_remove = false;
    for (size_t idx = 0; idx < m_genotype_file_names.size(); ++idx)
    {
        // go through each genotype file
        prefix = m_genotype_file_names[idx];
        bim_name = prefix + ".bim";
        bed_name = prefix + ".bed";
        // make sure we reset the flag of the ifstream by closing it before use
        if (bim.is_open()) bim.close();
        bim.clear();
        bim.open(bim_name.c_str());
        if (!bim.is_open())
        {
            std::string error_message =
                "Error: Cannot open bim file: " + bim_name;
            throw std::runtime_error(error_message);
        }
        // First pass, get the number of marker in bed & bim
        // as we want the number, num_snp_read will start at 0
        num_snp_read = 0;
        prev_chr = "";
        while (std::getline(bim, line))
        {
            misc::trim(line);
            if (line.empty()) continue;
            // don't bother to check, if the
            // bim file is malformed i.e. with less than 6
            // column. It will likely be captured when we do
            // size check and when we do the reading later on
            ++num_snp_read;
        }
        bim.clear();
        bim.seekg(0, bim.beg);
        // check if the bed file is valid
        check_bed(bed_name, num_snp_read, bed_offset);
        // now go through the bim file and perform filtering
        num_snp_read = 0;
        while (std::getline(bim, line))
        {
            misc::trim(line);
            if (line.empty()) continue;
            // we need to remember the actual number read is num_snp_read+1
            ++num_snp_read;
            bim_token = misc::split(line);
            if (bim_token.size() < 6)
            {
                std::string error_message =
                    "Error: Malformed bim file. Less than 6 column on "
                    "line: "
                    + misc::to_string(num_snp_read) + "\n";
                throw std::runtime_error(error_message);
            }
            if (genotype->m_existed_snps_index.find(bim_token[+BIM::RS])
                == genotype->m_existed_snps_index.end())
            {
                ++m_base_missed;
                continue;
            }
            // ensure all alleles are capitalized for easy matching
            std::transform(bim_token[+BIM::A1].begin(),
                           bim_token[+BIM::A1].end(),
                           bim_token[+BIM::A1].begin(), ::toupper);
            std::transform(bim_token[+BIM::A2].begin(),
                           bim_token[+BIM::A2].end(),
                           bim_token[+BIM::A2].begin(), ::toupper);

            // read in the chromosome string
            chr = bim_token[+BIM::CHR];
            // check if this is a new chromosome. If this is a new chromosome,
            // check if we want to remove it
            if (chr != prev_chr)
            {
                // get the chromosome code using PLINK 2 function
                chr_code = get_chrom_code_raw(chr.c_str());
                // check if we want to skip this chromosome
                if (chr_code_check(chr_code, chr_sex_error, chr_error,
                                   error_message))
                {
                    // only print chr error message if we haven't already
                    if (chr_error && !prev_chr_error)
                    {
                        std::cerr << error_message << "\n";
                        prev_chr_sex_error = chr_error;
                    }
                    // only print sex chr error message if we haven't already
                    else if (chr_sex_error && !prev_chr_sex_error)
                    {
                        std::cerr << error_message << "\n";
                        prev_chr_sex_error = chr_sex_error;
                    }
                    continue;
                }
                // only update the prev_chr after we have done the checking
                // this will help us to continue to skip all SNPs that are
                // supposed to be removed instead of the first entry
                prev_chr = chr;
                if (chr_code < 0) chr_num = ~size_t(0);
                chr_num = static_cast<size_t>(chr_code);
            }

            // now read in the coordinate
            size_t loc = ~size_t(0);
            try
            {
                loc = misc::string_to_size_t(bim_token[+BIM::BP].c_str());
            }
            catch (...)
            {
                std::string error_message =
                    "Error: Invalid SNP coordinate: " + bim_token[+BIM::RS]
                    + ":" + bim_token[+BIM::BP] + "\n";
                error_message.append("Please check you have the correct input");
                throw std::runtime_error(error_message);
            }
            to_remove =
                Genotype::within_region(exclusion_regions, chr_num, loc);
            if (to_remove)
            {
                m_num_xrange++;
                continue;
            }
            // check if this is a duplicated SNP
            if (processed_snps.find(bim_token[+BIM::RS])
                != processed_snps.end())
            { duplicated_snp.insert(bim_token[+BIM::RS]); }
            else if (!ambiguous(bim_token[+BIM::A1], bim_token[+BIM::A2])
                     || m_keep_ambig)
            {
                // if the SNP is not ambiguous (or if we want to keep ambiguous
                // SNPs), we will start processing the bed file (if required)

                // now read in the binary information and determine if we
                // want to keep this SNP only do the filtering if we need to
                // as my current implementation isn't as efficient as PLINK
                m_num_ambig +=
                    ambiguous(bim_token[+BIM::A1], bim_token[+BIM::A2]);
                auto&& ref_index =
                    genotype->m_existed_snps_index[bim_token[+BIM::RS]];
                if (!genotype->m_existed_snps[ref_index].matching(
                        chr_num, loc, bim_token[+BIM::A1], bim_token[+BIM::A2],
                        flipping))
                {
                    // mismatch found between base target and reference
                    if (!mismatch_snp_record.is_open())
                    {
                        // open the file accordingly
                        if (m_mismatch_file_output)
                        {
                            mismatch_snp_record.open(
                                mismatch_snp_record_name.c_str(),
                                std::ofstream::app);
                            if (!mismatch_snp_record.is_open())
                            {
                                throw std::runtime_error(
                                    std::string("Cannot open mismatch file to "
                                                "write: "
                                                + mismatch_snp_record_name));
                            }
                        }
                        else
                        {
                            mismatch_snp_record.open(
                                mismatch_snp_record_name.c_str());
                            if (!mismatch_snp_record.is_open())
                            {
                                throw std::runtime_error(
                                    std::string("Cannot open mismatch file to "
                                                "write: "
                                                + mismatch_snp_record_name));
                            }
                            mismatch_snp_record
                                << "File_Type\tRS_ID\tCHR_Target\tCHR_"
                                   "File\tBP_Target\tBP_File\tA1_"
                                   "Target\tA1_File\tA2_Target\tA2_"
                                   "File\n";
                        }
                    }
                    m_mismatch_file_output = true;
                    if (m_is_ref)
                    {
                        mismatch_snp_record
                            << "Reference\t" << bim_token[+BIM::RS] << "\t"
                            << target->m_existed_snps[ref_index].chr() << "\t"
                            << chr_num << "\t"
                            << target->m_existed_snps[ref_index].loc() << "\t"
                            << loc << "\t"
                            << target->m_existed_snps[ref_index].ref() << "\t"
                            << bim_token[+BIM::A1] << "\t"
                            << target->m_existed_snps[ref_index].alt() << "\t"
                            << bim_token[+BIM::A2] << "\n";
                        m_num_ref_target_mismatch++;
                    }
                    else
                    {
                        // base can have missing chr and loc info, which will be
                        // parse to ~size_t(0)
                        mismatch_snp_record << "Base\t" << bim_token[+BIM::RS]
                                            << "\t" << chr_num << "\t";
                        if (m_existed_snps[ref_index].chr() == ~size_t(0))
                        { mismatch_snp_record << "-\t"; }
                        else
                        {
                            mismatch_snp_record
                                << m_existed_snps[ref_index].chr() << "\t";
                        }
                        mismatch_snp_record << loc << "\t";
                        if (m_existed_snps[ref_index].loc() == ~size_t(0))
                        { mismatch_snp_record << "-\t"; }
                        else
                        {
                            mismatch_snp_record
                                << m_existed_snps[ref_index].loc() << "\t";
                        }
                        mismatch_snp_record
                            << bim_token[+BIM::A1] << "\t"
                            << m_existed_snps[ref_index].ref() << "\t"
                            << bim_token[+BIM::A2] << "\t"
                            << m_existed_snps[ref_index].alt() << "\n";
                        m_num_ref_target_mismatch++;
                    }
                }
                else
                {
                    byte_pos = bed_offset
                               + ((num_snp_read - 1) * (unfiltered_sample_ct4));
                    if (m_is_ref)
                        // here, the flipping is relative to target
                        target->m_existed_snps[ref_index].add_reference(
                            idx, byte_pos, flipping);
                    else
                        m_existed_snps[ref_index].add_target(
                            idx, byte_pos, chr_num, loc, bim_token[+BIM::A1],
                            bim_token[+BIM::A2], flipping);
                    processed_snps.insert(bim_token[+BIM::RS]);
                    retain_snp[ref_index] = true;
                    num_retained++;
                }
            }
            else
            {
                ++m_num_ambig;
            }
        }
        bim.close();
    }
    // try to release memory
    if (num_retained != genotype->m_existed_snps.size())
    {
        // remove any SNP that is not retained, the ref_retain vector should
        // have the same ID as the target->m_existed_snps so we can use it
        // directly for SNP removal
        genotype->shrink_snp_vector(retain_snp);
        // need to update index search after we updated the vector
        genotype->update_snp_index();
    }
    if (duplicated_snp.size() != 0)
    {
        // there are duplicated SNPs, we will need to terminate with the
        // information
        std::ofstream log_file_stream;
        std::string dup_name = out_prefix + ".valid";
        log_file_stream.open(dup_name.c_str());
        if (!log_file_stream.is_open())
        {
            std::string error_message = "Error: Cannot open file: " + dup_name;
            throw std::runtime_error(error_message);
        }
        // we should not use m_existed_snps unless it is for reference
        for (auto&& snp : genotype->m_existed_snps)
        {
            // we only output the valid SNPs.
            if (duplicated_snp.find(snp.rs()) == duplicated_snp.end())
                log_file_stream << snp.rs() << "\t" << snp.chr() << "\t"
                                << snp.loc() << "\t" << snp.ref() << "\t"
                                << snp.alt() << "\n";
        }
        log_file_stream.close();
        std::string error_message =
            "Error: A total of " + std::to_string(duplicated_snp.size())
            + " duplicated SNP ID detected out of "
            + misc::to_string(m_existed_snps.size())
            + " input SNPs! Valid SNP ID (post --extract / "
              "--exclude, non-duplicated SNPs) stored at "
            + dup_name + ". You can avoid this error by using --extract "
            + dup_name;
        throw std::runtime_error(error_message);
    }
}


void BinaryPlink::check_bed(const std::string& bed_name, size_t num_marker,
                            uintptr_t& bed_offset)
{
    bed_offset = 3;
    uint32_t uii = 0;
    int64_t llxx = 0;
    int64_t llyy = 0;
    int64_t llzz = 0;
    uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    std::ifstream bed(bed_name.c_str(), std::ios::binary);
    if (!bed.is_open())
    {
        std::string error_message = "Cannot read bed file: " + bed_name;
        throw std::runtime_error(error_message);
    }
    bed.seekg(0, bed.end);
    llxx = bed.tellg();
    if (!llxx)
    {
        std::string error_message = "Error: Empty .bed file: " + bed_name;
        throw std::runtime_error(error_message);
    }
    bed.seekg(0, bed.beg);
    bed.clear();
    char version_check[3];
    bed.read(version_check, 3);
    uii = static_cast<uint32_t>(bed.gcount());
    llyy = static_cast<int64_t>((static_cast<uint64_t>(unfiltered_sample_ct4))
                                * num_marker);
    llzz = static_cast<int64_t>(static_cast<uint64_t>(m_unfiltered_sample_ct)
                                * ((num_marker + 3) / 4));
    bool sample_major = false;
    // compare only the first 3 bytes
    if ((uii == 3) && (!memcmp(version_check, "l\x1b\x01", 3))) { llyy += 3; }
    else if ((uii == 3) && (!memcmp(version_check, "l\x1b", 3)))
    {
        // v1.00 sample-major
        sample_major = true;
        llyy = llzz + 3;
        bed_offset = 2;
    }
    else if (uii && (*version_check == '\x01'))
    {
        // v0.99 SNP-major
        llyy += 1;
        bed_offset = 1;
    }
    else if (uii && (!(*version_check)))
    {
        // v0.99 sample-major
        sample_major = true;
        llyy = llzz + 1;
        bed_offset = 2;
    }
    else
    {
        // pre-v0.99, sample-major, no header bytes
        sample_major = true;
        if (llxx != llzz)
        {
            // probably not PLINK-format at all, so give this error instead
            // of "invalid file size"

            std::string error_message =
                "Error: Invalid header bytes in .bed file: " + bed_name;
            throw std::runtime_error(error_message);
        }
        llyy = llzz;
        bed_offset = 2;
    }
    if (llxx != llyy)
    {
        if ((*version_check == '#')
            || ((uii == 3) && (!memcmp(version_check, "chr", 3))))
        {
            std::string error_message = "Error: Invalid header bytes in PLINK "
                                        "1 .bed file: "
                                        + bed_name
                                        + "  (Is this a UCSC "
                                          "Genome\nBrowser BED file instead?)";
            throw std::runtime_error(error_message);
        }
        else
        {
            std::string error_message =
                "Error: Invalid .bed file size for " + bed_name;
            throw std::runtime_error(error_message);
        }
    }
    if (sample_major)
    {
        throw std::runtime_error(
            "Error: Currently do not support sample major format");
    }
    bed.close();
}

BinaryPlink::~BinaryPlink() {}

void BinaryPlink::read_score(
    const std::vector<size_t>::const_iterator& start_idx,
    const std::vector<size_t>::const_iterator& end_idx, bool reset_zero,
    const bool use_ref_maf)
{
    // for removing unwanted bytes from the end of the genotype vector
    const uintptr_t final_mask =
        get_final_mask(static_cast<uint32_t>(m_sample_ct));
    // this is use for initialize the array sizes
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    uint32_t ll_ct, lh_ct, hh_ct;
    uint32_t ll_ctf, lh_ctf, hh_ctf;
    // for storing the count of each observation
    size_t homrar_ct = 0;
    size_t missing_ct = 0;
    size_t het_ct = 0;
    size_t homcom_ct = 0;
    size_t tmp_total = 0;
    const size_t ploidy = 2;
    // those are the weight (0,1,2) for each genotype observation
    double homcom_weight = m_homcom_weight;
    double het_weight = m_het_weight;
    double homrar_weight = m_homrar_weight;
    // this is required if we want to calculate the MAF from the genotype (for
    // imputation of missing genotype)
    // if we want to set the missing score to zero, miss_count will equal to 0,
    // 1 otherwise
    const size_t miss_count =
        (m_missing_score != MISSING_SCORE::SET_ZERO) * ploidy;
    // this indicate if we want the mean of the genotype to be 0 (missingness =
    // 0)
    const bool is_centre = (m_missing_score == MISSING_SCORE::CENTER);
    // this indicate if we want to impute the missing genotypes using the
    // population mean
    const bool mean_impute = (m_missing_score == MISSING_SCORE::MEAN_IMPUTE);
    // check if it is not the frist run, if it is the first run, we will reset
    // the PRS to zero instead of addint it up
    bool not_first = !reset_zero;
    double stat, maf, adj_score, miss_score;
    // m_cur_file = ""; // just close it
    // if (m_bed_file.is_open()) { m_bed_file.close(); }
    // initialize the genotype vector to store the binary genotypes
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    std::vector<size_t>::const_iterator cur_idx = start_idx;
    unsigned long long cur_line;
    size_t cur_file_idx;
    for (; cur_idx != end_idx; ++cur_idx)
    {
        auto&& cur_snp = m_existed_snps[(*cur_idx)];
        cur_file_idx = cur_snp.file_index();
        auto&& cur_map = m_genotype_files[cur_file_idx];
        /*
                if (m_cur_file != cur_snp.file_name())
                { update_bed(cur_snp.file_name()); }
        */
        // current location of the snp in the bed file
        // allow for quick jumping
        // very useful for read score as most SNPs might not
        // be next to each other
        cur_line = cur_snp.byte_pos();

        const size_t max_file_size = cur_map.mapped_length();
        // we now read the genotype from the file by calling
        // load_and_collapse_incl
        // important point to note here is the use of m_sample_include and
        // m_sample_ct instead of using the m_founder m_founder_info as the
        // founder vector is for LD calculation whereas the sample_include is
        // for PRS

        if (cur_line + unfiltered_sample_ct4 >= max_file_size)
        {
            std::string error_message =
                "Erorr: Reading out of bound: " + misc::to_string(cur_line)
                + " " + misc::to_string(unfiltered_sample_ct4) + " "
                + misc::to_string(max_file_size);
            throw std::runtime_error(error_message);
        }
        char* geno = reinterpret_cast<char*>(m_tmp_genotype.data());
        for (size_t i = 0; i < unfiltered_sample_ct4; ++i)
        {
            *geno = cur_map[cur_line + i];
            ++geno;
        }
        if (!cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct,
                                use_ref_maf))
        {
            // we need to calculate the MA
            single_marker_freqs_and_hwe(
                unfiltered_sample_ctv2, m_tmp_genotype.data(),
                m_sample_include2.data(), m_founder_include2.data(),
                m_sample_ct, &ll_ct, &lh_ct, &hh_ct, m_founder_ct, &ll_ctf,
                &lh_ctf, &hh_ctf);
            homcom_ct = ll_ctf;
            het_ct = lh_ctf;
            homrar_ct = hh_ctf;
            tmp_total = (homcom_ct + het_ct + homrar_ct);
            assert(m_founder_ct >= tmp_total);
            missing_ct = m_founder_ct - tmp_total;
            cur_snp.set_counts(homcom_ct, het_ct, homrar_ct, missing_ct);
        }
        if (m_unfiltered_sample_ct != m_sample_ct)
        {
            copy_quaterarr_nonempty_subset(
                m_tmp_genotype.data(), m_sample_include.data(),
                static_cast<uint32_t>(m_unfiltered_sample_ct),
                static_cast<uint32_t>(m_sample_ct), genotype.data());
        }
        else
        {
            genotype = m_tmp_genotype;
            genotype[(m_unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
        }
        // directly read in the current location
        if (m_founder_ct == missing_ct)
        {
            // problematic snp
            cur_snp.invalid();
            continue;
        }
        homcom_weight = m_homcom_weight;
        het_weight = m_het_weight;
        homrar_weight = m_homrar_weight;
        if (cur_snp.is_flipped()) { std::swap(homcom_weight, homrar_weight); }
        maf =
            static_cast<double>(homcom_weight * homcom_ct + het_ct * het_weight
                                + homrar_weight * homrar_ct)
            / (static_cast<double>((homcom_ct + het_ct + homrar_ct) * ploidy));
        stat = cur_snp.stat();
        adj_score = 0;
        if (is_centre) { adj_score = ploidy * stat * maf; }
        miss_score = 0;
        if (mean_impute) { miss_score = ploidy * stat * maf; }
        // now we go through the SNP vector
        read_prs(genotype, ploidy, stat, adj_score, miss_score, miss_count,
                 homcom_weight, het_weight, homrar_weight, not_first);
        // indicate that we've already read in the first SNP and no longer need
        // to reset the PRS
        not_first = true;
    }
}
