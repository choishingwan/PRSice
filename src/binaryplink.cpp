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

#include "binaryplink.hpp"


BinaryPlink::BinaryPlink(const Commander& commander, Reporter& reporter,
                         const bool is_ref)
{
    m_thread = static_cast<uint32_t>(commander.thread());
    m_ignore_fid = commander.ignore_fid();
    m_keep_nonfounder = commander.nonfounders();
    m_keep_ambig = commander.keep_ambig();
    m_is_ref = is_ref;
    // set the chromosome information
    // will need to add more script here if we want to support something other
    // than human
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    m_haploid_mask.resize(CHROM_MASK_WORDS, 0);
    // main use of following function is to set the max code
    init_chr();
    std::string message = "Loading Genotype ";
    if (is_ref) {
        std::string reference_name;
        if (commander.ref_list(reference_name)) {
            // has listed input
            // check if there is an external sample file
            std::vector<std::string> token = misc::split(reference_name, ",");
            bool external_sample = false;
            if (token.size() == 2) {
                m_sample_file = token[1];
                reference_name = token[0];
                external_sample = true;
            }
            message.append(" info from file " + reference_name + " (bed)\n");
            if (external_sample) {
                message.append("With external fam file: " + m_sample_file
                               + "\n");
            }
            m_genotype_files = load_genotype_prefix(reference_name);
            if (!external_sample) {
                m_sample_file = m_genotype_files.front() + ".fam";
            }
        }
        else
        {
            // single file input, check for # and replace it with 1-22
            commander.ref_name(reference_name);
            std::vector<std::string> token = misc::split(reference_name, ",");
            bool external_sample = false;
            if (token.size() == 2) {
                m_sample_file = token[1];
                reference_name = token[0];
                external_sample = true;
            }
            message.append(" file: " + reference_name + " (bed)\n");
            if (external_sample) {
                message.append("With external fam file: " + m_sample_file
                               + "\n");
            }
            m_genotype_files = set_genotype_files(reference_name);
            if (!external_sample) {
                m_sample_file = m_genotype_files.front() + ".fam";
            }
        }
    }
    else
    {
        std::string target_name;
        if (commander.target_list(target_name)) {
            // has listed input
            // check if there is an external sample file
            std::vector<std::string> token = misc::split(target_name, ",");
            bool external_sample = false;
            if (token.size() == 2) {
                m_sample_file = token[1];
                target_name = token[0];
                external_sample = true;
            }
            message.append(" info from file " + target_name + " (bed)\n");
            if (external_sample) {
                message.append("With external fam file: " + m_sample_file
                               + "\n");
            }
            m_genotype_files = load_genotype_prefix(target_name);
            if (!external_sample) {
                m_sample_file = m_genotype_files.front() + ".fam";
            }
        }
        else
        {
            // single file input, check for # and replace it with 1-22
            target_name = commander.target_name();
            std::vector<std::string> token = misc::split(target_name, ",");
            bool external_sample = false;
            if (token.size() == 2) {
                m_sample_file = token[1];
                target_name = token[0];
                external_sample = true;
            }
            message.append(" file: " + target_name + " (bed)\n");
            if (external_sample) {
                message.append("With external fam file: " + m_sample_file
                               + "\n");
            }
            m_genotype_files = set_genotype_files(target_name);
            if (!external_sample) {
                m_sample_file = m_genotype_files.front() + ".fam";
            }
        }
    }
    reporter.report(message);
}

std::vector<Sample_ID> BinaryPlink::gen_sample_vector()
{
    assert(m_genotype_files.size() > 0);
    std::ifstream famfile;
    famfile.open(m_sample_file.c_str());
    if (!famfile.is_open()) {
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
    while (std::getline(famfile, line)) {
        misc::trim(line);
        if (!line.empty()) {
            std::vector<std::string> token = misc::split(line);
            if (token.size() < 6) {
                std::string message =
                    "Error: Malformed fam file. Less than 6 column on "
                    "line: "
                    + std::to_string(m_unfiltered_sample_ct + 1) + "\n";
                throw std::runtime_error(message);
            }
            // fam file will always have both the FID and IID
            // though this _ concatination might cause a bug if
            // there is a situation where a sample with FID A_B IID A and
            // a sample with FID A and IID B_A, then this will be
            // in-distinguishable
            founder_info.insert(token[+FAM::FID] + "_" + token[+FAM::IID]);
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
    while (std::getline(famfile, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (token.size() < 6) {
            std::string error_message =
                "Error: Malformed fam file. Less than 6 column on line: "
                + std::to_string(sample_index + 1);
            throw std::runtime_error(error_message);
        }
        std::string id = (m_ignore_fid)
                             ? token[+FAM::IID]
                             : token[+FAM::FID] + "_" + token[+FAM::IID];
        if (!m_remove_sample) {
            // we don't want to include this sample if it is not found in the
            // selection_list
            inclusion = (m_sample_selection_list.find(id)
                         != m_sample_selection_list.end());
        }
        else
        {
            // we don't want to include this sample if it is found in the
            // selection_list
            inclusion = (m_sample_selection_list.find(id)
                         == m_sample_selection_list.end());
        }

        if (founder_info.find(token[+FAM::FATHER]) == founder_info.end()
            && founder_info.find(token[+FAM::MOTHER]) == founder_info.end()
            && inclusion)
        {
            // this is a founder (with no dad / mum)
            m_founder_ct++;
            // so m_founder_info is a subset of m_sample_include
            // and is used for LD calculation as we force LD calculation
            // to only be performed on founders
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
            m_num_non_founder++;
        }
        m_sample_ct += inclusion;
        if (token[+FAM::SEX] == "1") {
            m_num_male++;
        }
        else if (token[+FAM::SEX] == "2")
        {
            m_num_female++;
        }
        else
        {
            m_num_ambig_sex++;
        }
        sample_index++;
        if (duplicated_samples.find(id) != duplicated_samples.end())
            duplicated_sample_id.push_back(id);
        // only store samples that we need, and use the m_sample_include and
        // m_founder_info to indicate if sample is needed for subsequent
        // operations
        if (inclusion && !m_is_ref) {
            sample_name.emplace_back(
                Sample_ID(token[+FAM::FID], token[+FAM::IID],
                          token[+FAM::PHENOTYPE], founder));
        }
        duplicated_samples.insert(id);
    }

    if (!duplicated_sample_id.empty()) {
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
    m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    // m_prs_info.reserve(m_sample_ct);
    // now we add the prs information. For some reason, we can't do a simple
    // reserve
    for (size_t i = 0; i < m_sample_ct; ++i) {
        m_prs_info.emplace_back(PRS());
    }
    // also resize the in_regression flag
    m_in_regression.resize(m_sample_include.size(), 0);
    return sample_name;
}


std::vector<SNP> BinaryPlink::gen_snp_vector(const Commander& commander,
                                             Region& exclusion,
                                             Genotype* target)
{
    const std::string out_prefix = commander.out();
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    double maf_threshold = 0.0;
    const bool maf_filter = m_is_ref ? commander.ref_maf(maf_threshold)
                                     : commander.target_maf(maf_threshold);
    double geno_threshold = 0.0;
    const bool geno_filter = m_is_ref ? commander.ref_geno(geno_threshold)
                                      : commander.target_geno(geno_threshold);
    const uintptr_t final_mask =
        get_final_mask(static_cast<uint32_t>(m_sample_ct));
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    const uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_sample_ct);
    std::unordered_set<std::string> duplicated_snp;
    std::vector<SNP> snp_info;
    std::vector<std::string> bim_token;
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    std::vector<bool> ref_retain;
    if (m_is_ref) ref_retain.resize(target->m_existed_snps.size(), false);
    std::ifstream bim, bed;
    std::ofstream mismatch_snp_record;
    std::string bim_name, bed_name, chr, line;
    std::string prev_chr = "", error_message = "";
    std::string mismatch_snp_record_name = out_prefix + ".mismatch";
    double cur_maf;
    std::streampos byte_pos;

    intptr_t nanal = 0;
    uint32_t homrar_ct = 0, missing_ct = 0, het_ct = 0, homcom_ct = 0,
             num_ref_target_match = 0;
    int chr_code = 0;
    int prev_snp_processed = 0, num_snp_read = -1;
    bool chr_error = false, chr_sex_error = false, prev_chr_sex_error = false,
         prev_chr_error = false, has_count = false, dummy;
    // initialize the sample inclusion mask
    m_sample_mask.resize(pheno_nm_ctv2);
    // fill it with the required mask (copy from PLINK2)
    fill_quatervec_55(static_cast<uint32_t>(m_sample_ct), m_sample_mask.data());
    uintptr_t bed_offset;
    for (auto prefix : m_genotype_files) {
        // go through each genotype file
        bim_name = prefix + ".bim";
        bed_name = prefix + ".bed";
        // make sure we reset the flag of the ifstream by closing it before use
        if (bim.is_open()) bim.close();
        if (bed.is_open()) bed.close();
        bim.clear();
        bed.clear();
        bim.open(bim_name.c_str());
        if (!bim.is_open()) {
            std::string error_message =
                "Error: Cannot open bim file: " + bim_name;
            throw std::runtime_error(error_message);
        }
        // First pass, get the number of marker in bed & bim
        // as we want the number, num_snp_read will start at 0
        num_snp_read = 0;
        prev_chr = "";
        while (std::getline(bim, line)) {
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
        check_bed(bed_name, static_cast<size_t>(num_snp_read), bed_offset);
        bed.open(bed_name.c_str());
        if (!bed.is_open()) {
            std::string error_message =
                "Error: Cannot open bed file: " + bed_name;
            throw std::runtime_error(error_message);
        }

        // now go through the bim & bed file and perform filtering
        // reset # SNP read to -1 such that we can avoid troublesome -1 lateron
        num_snp_read = -1;

        // ensure prev_snp_processed = -2 so that we will always perform
        // seek for the first SNP (which will account for the bed offset)
        prev_snp_processed = -2;
        while (std::getline(bim, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            // we need to remember the actual number read is num_snp_read+1
            ++num_snp_read;
            bim_token = misc::split(line);
            if (bim_token.size() < 6) {
                std::string error_message =
                    "Error: Malformed bim file. Less than 6 column on "
                    "line: "
                    + misc::to_string(num_snp_read + 1) + "\n";
                throw std::runtime_error(error_message);
            }
            if (m_is_ref) {
                // for the reference panel
                if (target->m_existed_snps_index.find(bim_token[+BIM::RS])
                    == target->m_existed_snps_index.end())
                {
                    // Skip SNPs not found in the target file
                    continue;
                }
            }
            // ensure all alleles are capitalized for easy matching
            std::transform(bim_token[+BIM::A1].begin(),
                           bim_token[+BIM::A1].end(),
                           bim_token[+BIM::A1].begin(), ::toupper);
            std::transform(bim_token[+BIM::A2].begin(),
                           bim_token[+BIM::A2].end(),
                           bim_token[+BIM::A2].begin(), ::toupper);

            // exclude SNPs that are not required
            if (!m_is_ref) {
                // don't bother doing it when reading reference genome
                // as all SNPs should have been removed in target
                if (!m_exclude_snp
                    && m_snp_selection_list.find(bim_token[+BIM::RS])
                           == m_snp_selection_list.end())
                {
                    continue;
                }
                else if (m_exclude_snp
                         && m_snp_selection_list.find(bim_token[+BIM::RS])
                                != m_snp_selection_list.end())
                {
                    continue;
                }
            }

            // read in the chromosome string
            chr = bim_token[+BIM::CHR];
            // check if this is a new chromosome. If this is a new chromosome,
            // check if we want to remove it
            if (chr != prev_chr) {
                // get the chromosome code using PLINK 2 function
                chr_code = get_chrom_code_raw(chr.c_str());
                // check if we want to skip this chromosome
                if (chr_code_check(chr_code, chr_sex_error, chr_error,
                                   error_message))
                {
                    // only print chr error message if we haven't already
                    if (chr_error && !prev_chr_error) {
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
            }

            // now read in the coordinate
            int loc = -1;
            try
            {
                loc = misc::convert<int>(bim_token[+BIM::BP]);
                if (loc < 0) {
                    // coordinate must >= 0
                    std::string error_message =
                        "Error: SNP with negative corrdinate: "
                        + bim_token[+BIM::RS] + ":" + bim_token[+BIM::BP]
                        + "\n";
                    error_message.append(
                        "Please check you have the correct input");
                    throw std::runtime_error(error_message);
                }
            }
            catch (...)
            {
                std::string error_message =
                    "Error: SNP with non-numeric corrdinate: "
                    + bim_token[+BIM::RS] + ":" + bim_token[+BIM::BP] + "\n";
                error_message.append("Please check you have the correct input");
                throw std::runtime_error(error_message);
            }
            // check if we want to exclude this SNP because this fall within the
            // exclusion region(s)
            if (exclusion.check_exclusion(chr, loc)) {
                continue;
            }
            // check if this is a duplicated SNP
            if (m_existed_snps_index.find(bim_token[+BIM::RS])
                != m_existed_snps_index.end())
            {
                duplicated_snp.insert(bim_token[+BIM::RS]);
            }
            else if (!ambiguous(bim_token[+BIM::A1], bim_token[+BIM::A2])
                     || m_keep_ambig)
            {
                // if the SNP is not ambiguous (or if we want to keep ambiguous
                // SNPs), we will start processing the bed file (if required)

                // now read in the binary information and determine if we
                // want to keep this SNP only do the filtering if we need to
                // as my current implementation isn't as efficient as PLINK
                byte_pos = bed_offset
                           + (num_snp_read
                              * (static_cast<uint64_t>(unfiltered_sample_ct4)));
                if (geno_filter || maf_filter) {
                    // indicate we've already read the maf count
                    has_count = true;
                    if (num_snp_read - prev_snp_processed > 1) {
                        // only skip line if we are not reading sequentially
                        if (!bed.seekg(byte_pos, std::ios_base::beg)) {
                            std::string error_message =
                                "Error: Cannot read the bed file(seek): "
                                + bed_name;
                            throw std::runtime_error(error_message);
                        }
                    }
                    prev_snp_processed = num_snp_read;
                    // read in the genotype information to the genotype vector
                    if (load_and_collapse_incl(
                            static_cast<uint32_t>(m_unfiltered_sample_ct),
                            static_cast<uint32_t>(m_sample_ct),
                            m_sample_include.data(), final_mask, false, bed,
                            m_tmp_genotype.data(), genotype.data()))
                    {
                        std::string error_message =
                            "Error: Cannot read the bed file(read): "
                            + bed_name;
                        throw std::runtime_error(error_message);
                    }
                    // calculate the MAF using PLINK2 function
                    genovec_3freq(genotype.data(), m_sample_mask.data(),
                                  pheno_nm_ctv2, &missing_ct, &het_ct,
                                  &homcom_ct);
                    // calculate the remaining samples
                    nanal = static_cast<uint32_t>(m_sample_ct) - missing_ct;
                    // calculate the hom rare count
                    homrar_ct =
                        static_cast<uint32_t>(nanal) - het_ct - homcom_ct;
                    if (nanal == 0) {
                        // none of the sample contain this SNP
                        // still count as MAF filtering (for now)
                        m_num_maf_filter++;
                        continue;
                    }
                    // filter by genotype missingness
                    if (geno_filter
                        && static_cast<double>(missing_ct)
                                   / static_cast<double>(m_sample_ct)
                               > geno_threshold)
                    {
                        m_num_geno_filter++;
                        continue;
                    }
                    // filter by MAF
                    cur_maf = (static_cast<double>(het_ct + homrar_ct * 2)
                               / (static_cast<double>(nanal) * 2.0));
                    if (cur_maf > 0.5) cur_maf = 1.0 - cur_maf;
                    // remove SNP if maf lower than threshold
                    if (maf_filter && cur_maf < maf_threshold) {
                        m_num_maf_filter++;
                        continue;
                    }
                }
                // we have now completed the geno / maf filtering
                m_num_ambig +=
                    ambiguous(bim_token[+BIM::A1], bim_token[+BIM::A2]);
                if (!m_is_ref) {
                    // only push in the SNP if this is not the reference panel.
                    // For reference panel, we just add the coordinate to the
                    // target to save memory usage
                    m_existed_snps_index[bim_token[+BIM::RS]] = snp_info.size();
                    if (has_count)
                        snp_info.emplace_back(SNP(bim_token[+BIM::RS], chr_code,
                                                  loc, bim_token[+BIM::A1],
                                                  bim_token[+BIM::A2], prefix,
                                                  byte_pos, homcom_ct, het_ct,
                                                  homrar_ct, missing_ct));
                    else
                        snp_info.emplace_back(SNP(bim_token[+BIM::RS], chr_code,
                                                  loc, bim_token[+BIM::A1],
                                                  bim_token[+BIM::A2], prefix,
                                                  byte_pos));
                }
                else
                {

                    auto&& target_index =
                        target->m_existed_snps_index[bim_token[+BIM::RS]];
                    if (!target->m_existed_snps[target_index].matching(
                            chr_code, loc, bim_token[+BIM::A1],
                            bim_token[+BIM::A2], dummy))
                    {
                        // For reference panel, check if the reference and
                        // target matches.
                        // here, this is a mismatch, so we will output the
                        // information to .mismatch file
                        if (!mismatch_snp_record.is_open()) {
                            // open the file accordingly
                            if (m_mismatch_file_output) {
                                mismatch_snp_record.open(
                                    mismatch_snp_record_name.c_str(),
                                    std::ofstream::app);
                                if (!mismatch_snp_record.is_open()) {
                                    throw std::runtime_error(std::string(
                                        "Cannot open mismatch file to "
                                        "write: "
                                        + mismatch_snp_record_name));
                                }
                            }
                            else
                            {
                                mismatch_snp_record.open(
                                    mismatch_snp_record_name.c_str());
                                if (!mismatch_snp_record.is_open()) {
                                    throw std::runtime_error(std::string(
                                        "Cannot open mismatch file to "
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
                        mismatch_snp_record
                            << "Reference\t" << bim_token[+BIM::RS] << "\t"
                            << target->m_existed_snps[target_index].chr()
                            << "\t" << chr_code << "\t"
                            << target->m_existed_snps[target_index].loc()
                            << "\t" << loc << "\t"
                            << target->m_existed_snps[target_index].ref()
                            << "\t" << bim_token[+BIM::A1] << "\t"
                            << target->m_existed_snps[target_index].alt()
                            << "\t" << bim_token[+BIM::A2] << "\n";
                        m_num_ref_target_mismatch++;
                    }
                    else
                    {
                        // this is not a mismatch, we can add it to the
                        // target
                        target->m_existed_snps[target_index].add_reference(
                            prefix, byte_pos);
                        ref_retain[target_index] = true;
                        num_ref_target_match++;
                    }
                }
            }
            else if (!m_keep_ambig)
            {
                // directly skip the ambiguous SNP if we don't want to keep it
                m_num_ambig++;
            }
        }
        bim.close();
        if (bed.is_open()) bed.close();
    }
    // try to release memory
    snp_info.shrink_to_fit();
    if (m_is_ref && num_ref_target_match != target->m_existed_snps.size()) {
        // remove any SNP that is not retained, the ref_retain vector should
        // have the same ID as the target->m_existed_snps so we can use it
        // directly for SNP removal
        target->m_existed_snps.erase(
            std::remove_if(
                target->m_existed_snps.begin(), target->m_existed_snps.end(),
                [&ref_retain, &target](const SNP& s) {
                    return !ref_retain[(&s - &*begin(target->m_existed_snps))];
                }),
            target->m_existed_snps.end());
        target->m_existed_snps.shrink_to_fit();

        // When reading the reference file, we actually update the
        // SNP list in target file. This lead to the index search in
        // target to have the wrong index. To avoid that, we need to
        // update the SNP index accordingly
        target->update_snp_index();
    }
    if (duplicated_snp.size() != 0) {
        // there are duplicated SNPs, we will need to terminate with the
        // information
        std::ofstream log_file_stream;
        std::string dup_name = out_prefix + ".valid";
        log_file_stream.open(dup_name.c_str());
        if (!log_file_stream.is_open()) {
            std::string error_message = "Error: Cannot open file: " + dup_name;
            throw std::runtime_error(error_message);
        }
        for (auto&& snp : (m_is_ref ? target->m_existed_snps : m_existed_snps))
        {
            if (duplicated_snp.find(snp.rs()) != duplicated_snp.end())
                log_file_stream << snp.rs() << "\t" << snp.chr() << "\t"
                                << snp.loc() << "\t" << snp.ref() << "\t"
                                << snp.alt() << "\n";
        }
        log_file_stream.close();
        std::string error_message =
            "Error: A total of " + std::to_string(duplicated_snp.size())
            + " duplicated SNP ID detected out of "
            + misc::to_string(snp_info.size())
            + " input SNPs! Valid SNP ID (post --extract / "
              "--exclude, non-duplicated SNPs) stored at "
            + dup_name + ". You can avoid this error by using --extract "
            + dup_name;
        throw std::runtime_error(error_message);
    }
    // mismatch_snp_record will close itself when it goes out of scope
    return snp_info;
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
    if (!bed.is_open()) {
        std::string error_message = "Cannot read bed file: " + bed_name;
        throw std::runtime_error(error_message);
    }
    bed.seekg(0, bed.end);
    llxx = bed.tellg();
    if (!llxx) {
        throw std::runtime_error("Error: Empty .bed file.");
    }
    bed.seekg(0, bed.beg);
    char version_check[3];
    bed.read(version_check, 3);
    uii = static_cast<uint32_t>(bed.gcount());
    llyy = static_cast<int64_t>((static_cast<uint64_t>(unfiltered_sample_ct4))
                                * num_marker);
    llzz = static_cast<int64_t>(static_cast<uint64_t>(m_unfiltered_sample_ct)
                                * ((num_marker + 3) / 4));
    bool sample_major = false;
    // compare only the first 3 bytes
    if ((uii == 3) && (!memcmp(version_check, "l\x1b\x01", 3))) {
        llyy += 3;
    }
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
        if (llxx != llzz) {
            // probably not PLINK-format at all, so give this error instead
            // of "invalid file size"
            throw std::runtime_error(
                "Error: Invalid header bytes in .bed file.");
        }
        llyy = llzz;
        bed_offset = 2;
    }
    if (llxx != llyy) {
        if ((*version_check == '#')
            || ((uii == 3) && (!memcmp(version_check, "chr", 3))))
        {
            throw std::runtime_error("Error: Invalid header bytes in PLINK "
                                     "1 .bed file.  (Is this a UCSC "
                                     "Genome\nBrowser BED file instead?)");
        }
        else
        {
            throw std::runtime_error("Error: Invalid .bed file size.");
        }
    }
    if (sample_major) {
        throw std::runtime_error(
            "Error: Currently do not support sample major format");
    }
    bed.close();
}

BinaryPlink::~BinaryPlink() {}


void BinaryPlink::read_score(std::vector<size_t>& index, bool reset_zero)
{
    // region_index should be the background index
    switch (m_model)
    {
    case MODEL::HETEROZYGOUS: read_score(index, 0, 1, 0, reset_zero); break;
    case MODEL::DOMINANT: read_score(index, 0, 1, 1, reset_zero); break;
    case MODEL::RECESSIVE: read_score(index, 0, 0, 1, reset_zero); break;
    default: read_score(index, 0, 1, 2, reset_zero); break;
    }
}


void BinaryPlink::read_score(std::vector<size_t>& index_bound,
                             uint32_t homcom_wt, uint32_t het_wt,
                             uint32_t homrar_wt, bool reset_zero)
{
    double stat, maf, adj_score, miss_score;
    const uintptr_t final_mask = get_final_mask(m_sample_ct);
    // for array size
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    uintptr_t* lbptr;
    uintptr_t ulii;
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uint32_t homrar_ct;
    uint32_t missing_ct = 0;
    uint32_t het_ct = 0;
    uint32_t homcom_ct = 0;
    uint32_t homcom_weight = homcom_wt;
    uint32_t het_weight = het_wt;
    uint32_t homrar_weight = homrar_wt;
    intptr_t nanal;
    // For set zero, miss_count will become 0
    const uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_sample_ct);
    const uint32_t miss_count = (m_missing_score != MISSING_SCORE::SET_ZERO);
    const bool is_centre = (m_missing_score == MISSING_SCORE::CENTER);
    const bool mean_impute = (m_missing_score == MISSING_SCORE::MEAN_IMPUTE);
    bool not_first = !reset_zero;
    m_cur_file = ""; // just close it
    if (m_bed_file.is_open()) {
        m_bed_file.close();
    }
    // index is w.r.t. partition, which contain all the information
    for (auto&& i_snp : index_bound) {
        // for each SNP
        auto&& cur_snp = m_existed_snps[i_snp];
        if (m_cur_file.empty() || m_cur_file.compare(cur_snp.file_name()) != 0)
        {
            // If we are processing a new file
            if (m_bed_file.is_open()) {
                m_bed_file.close();
            }
            m_cur_file = cur_snp.file_name();
            std::string bedname = m_cur_file + ".bed";
            m_bed_file.open(bedname.c_str(), std::ios::binary);
            if (!m_bed_file.is_open()) {
                std::string error_message =
                    "Error: Cannot open bed file: " + bedname;
                throw std::runtime_error(error_message);
            }
            m_prev_loc = 0;
        }
        // current location of the snp in the bed file
        // allow for quick jumping
        // very useful for read score as most SNPs might not
        // be next to each other
        std::streampos cur_line = cur_snp.byte_pos();
        if (m_prev_loc != cur_line
            && !m_bed_file.seekg(cur_line, std::ios_base::beg))
        {
            throw std::runtime_error("Error: Cannot read the bed file!");
        }
        m_prev_loc = cur_line + (std::streampos) unfiltered_sample_ct4;
        // loadbuf_raw is the temporary
        // loadbuff is where the genotype will be located
        if (load_and_collapse_incl(m_unfiltered_sample_ct, m_sample_ct,
                                   m_sample_include.data(), final_mask, false,
                                   m_bed_file, m_tmp_genotype.data(),
                                   genotype.data()))
        {
            throw std::runtime_error("Error: Cannot read the bed file!");
        }
        // try to calculate MAF here
        /*
        genovec_3freq(genotype.data(), sample_include2.data(),
        pheno_nm_ctv2, &missing_ct, &het_ct, &homcom_ct);
                      */
        cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct);
        if (homcom_ct + het_ct + homrar_ct + missing_ct == 0) {
            genovec_3freq(genotype.data(), m_sample_mask.data(), pheno_nm_ctv2,
                          &missing_ct, &het_ct, &homcom_ct);
            cur_snp.set_counts(homcom_ct, het_ct, homrar_ct, missing_ct);
        }
        nanal = m_sample_ct - missing_ct;
        if (nanal == 0) {
            cur_snp.invalidate();
            continue;
        }
        homcom_weight = homcom_wt;
        het_weight = het_wt;
        homrar_weight = homrar_wt;
        maf = (double) (het_ct * het_weight + homrar_weight * homrar_ct)
              / (double) (nanal * 2.0);
        if (cur_snp.is_flipped()) {
            // change the mean to reflect flipping
            maf = 1.0 - maf;
            // swap the weighting
            std::swap(homcom_weight, homrar_weight);
        }
        stat = cur_snp.stat() * 2; // Multiply by ploidy
        // we don't allow the use of center and mean impute together
        // if centre, missing = 0 anyway (kinda like mean imputed)
        // centre is 0 if false
        adj_score = stat * maf * is_centre;
        miss_score = stat * maf * mean_impute;


        // now we go through the SNP vector

        lbptr = genotype.data();
        uii = 0;
        ulii = 0;
        do
        {
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct) {
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            ujj = 0;
            while (ujj < BITCT) {
                ukk = (ulii >> ujj) & 3;
                auto&& sample_prs = m_prs_info[uii + (ujj / 2)];
                // now we will get all genotypes (0, 1, 2, 3)
                switch (ukk)
                {
                default:
                    sample_prs.num_snp = sample_prs.num_snp * not_first + 1;
                    sample_prs.prs = sample_prs.prs * not_first
                                     + homcom_weight * stat * 0.5 - adj_score;
                    break;
                case 1:
                    sample_prs.num_snp = sample_prs.num_snp * not_first + 1;
                    sample_prs.prs = sample_prs.prs * not_first
                                     + het_weight * stat * 0.5 - adj_score;
                    break;
                case 3:
                    sample_prs.num_snp = sample_prs.num_snp * not_first + 1;
                    sample_prs.prs = sample_prs.prs * not_first
                                     + homrar_weight * stat * 0.5 - adj_score;
                    break;
                case 2:
                    sample_prs.prs = sample_prs.prs * not_first + miss_score;
                    sample_prs.num_snp =
                        sample_prs.num_snp * not_first + miss_count;
                    break;
                }
                ulii &= ~((3 * ONELU) << ujj);
                ujj += 2;
            }
            uii += BITCT2;
        } while (uii < m_sample_ct);
        not_first = true;
    }
}


void BinaryPlink::read_score(size_t start_index, size_t end_bound,
                             uint32_t homcom_wt, uint32_t het_wt,
                             uint32_t homrar_wt, const size_t region_index,
                             bool set_zero)
{
    const uintptr_t final_mask = get_final_mask(m_sample_ct);
    // for array size
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    uintptr_t* lbptr;
    uintptr_t ulii;
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uint32_t homrar_ct = 0;
    uint32_t missing_ct = 0;
    uint32_t het_ct = 0;
    uint32_t homcom_ct = 0;
    uint32_t homcom_weight = homcom_wt;
    uint32_t het_weight = het_wt;
    uint32_t homrar_weight = homrar_wt;
    uint32_t temp_weight = 0;
    // For set zero, miss_count will become 0
    const uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_sample_ct);
    const uint32_t miss_count = (m_missing_score != MISSING_SCORE::SET_ZERO);
    const bool is_centre = (m_missing_score == MISSING_SCORE::CENTER);
    const bool mean_impute = (m_missing_score == MISSING_SCORE::MEAN_IMPUTE);
    bool not_first = !set_zero;
    intptr_t nanal;
    double stat, maf, adj_score, miss_score;
    m_cur_file = ""; // just close it
    if (m_bed_file.is_open()) {
        m_bed_file.close();
    }
    // index is w.r.t. partition, which contain all the information
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    for (size_t i_snp = start_index; i_snp < end_bound; ++i_snp) {
        // for each SNP
        auto&& cur_snp = m_existed_snps[i_snp];
        if (m_cur_file.empty() || m_cur_file.compare(cur_snp.file_name()) != 0)
        {
            // If we are processing a new file
            if (m_bed_file.is_open()) {
                m_bed_file.close();
            }
            m_cur_file = cur_snp.file_name();
            std::string bedname = m_cur_file + ".bed";
            m_bed_file.open(bedname.c_str(), std::ios::binary);
            if (!m_bed_file.is_open()) {
                std::string error_message =
                    "Error: Cannot open bed file: " + bedname;
                throw std::runtime_error(error_message);
            }
            m_prev_loc = 0;
        }
        // only read this SNP if it falls within our region of interest
        if (!cur_snp.in(region_index)) continue;
        // current location of the snp in the bed file
        // allow for quick jumping
        // very useful for read score as most SNPs might not
        // be next to each other
        std::streampos cur_line = cur_snp.byte_pos();
        if (m_prev_loc != cur_line
            && !m_bed_file.seekg(cur_line, std::ios_base::beg))
        {
            throw std::runtime_error("Error: Cannot read the bed file!");
        }
        m_prev_loc = cur_line + (std::streampos) unfiltered_sample_ct4;
        // loadbuf_raw is the temporary
        // loadbuff is where the genotype will be located
        if (load_and_collapse_incl(m_unfiltered_sample_ct, m_sample_ct,
                                   m_sample_include.data(), final_mask, false,
                                   m_bed_file, m_tmp_genotype.data(),
                                   genotype.data()))
        {
            throw std::runtime_error("Error: Cannot read the bed file!");
        }
        // try to calculate MAF here
        /*
        genovec_3freq(genotype.data(), sample_include2.data(),
        pheno_nm_ctv2, &missing_ct, &het_ct, &homcom_ct);
                      */
        bool has_count =
            cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct);

        if (!has_count) {
            genovec_3freq(genotype.data(), m_sample_mask.data(), pheno_nm_ctv2,
                          &missing_ct, &het_ct, &homcom_ct);
            cur_snp.set_counts(homcom_ct, het_ct, homrar_ct, missing_ct);
        }
        nanal = m_sample_ct - missing_ct;
        if (nanal == 0) {
            cur_snp.invalidate();
            continue;
        }
        homcom_weight = homcom_wt;
        het_weight = het_wt;
        homrar_weight = homrar_wt;
        maf = (double) (het_ct * het_weight + homrar_weight * homrar_ct)
              / (double) (nanal * 2.0);
        if (cur_snp.is_flipped()) {
            // change the mean to reflect flipping
            maf = 1.0 - maf;
            // swap the weighting
            temp_weight = homcom_weight;
            homcom_weight = homrar_weight;
            homrar_weight = temp_weight;
        }
        stat = cur_snp.stat() * 2; // Multiply by ploidy
        // we don't allow the use of center and mean impute together
        // if centre, missing = 0 anyway (kinda like mean imputed)
        // centre is 0 if false
        adj_score = stat * maf * is_centre;
        miss_score = stat * maf * mean_impute;


        // now we go through the SNP vector

        lbptr = genotype.data();
        uii = 0;
        ulii = 0;
        do
        {
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct) {
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            ujj = 0;
            while (ujj < BITCT) {
                ukk = (ulii >> ujj) & 3;
                auto&& sample_prs = m_prs_info[uii + (ujj / 2)];
                // now we will get all genotypes (0, 1, 2, 3)
                switch (ukk)
                {
                default:
                    sample_prs.num_snp = sample_prs.num_snp * not_first + 1;
                    sample_prs.prs = sample_prs.prs * not_first
                                     + homcom_weight * stat * 0.5 - adj_score;
                    break;
                case 1:
                    sample_prs.num_snp = sample_prs.num_snp * not_first + 1;
                    sample_prs.prs = sample_prs.prs * not_first
                                     + het_weight * stat * 0.5 - adj_score;
                    break;
                case 3:
                    sample_prs.num_snp = sample_prs.num_snp * not_first + 1;
                    sample_prs.prs = sample_prs.prs * not_first
                                     + homrar_weight * stat * 0.5 - adj_score;
                    break;
                case 2:
                    sample_prs.prs = sample_prs.prs * not_first + miss_score;
                    sample_prs.num_snp =
                        sample_prs.num_snp * not_first + miss_count;
                    break;
                }
                ulii &= ~((3 * ONELU) << ujj);
                ujj += 2;
            }
            uii += BITCT2;
        } while (uii < m_sample_ct);
        not_first = true;
    }
}
void BinaryPlink::read_score(size_t start_index, size_t end_bound,
                             const size_t region_index, bool set_zero)
{
    //	std::vector<size_t> index_bound(end_bound-start_index);
    //  std::iota(index_bound.begin(), index_bound.end(), start_index);

    switch (m_model)
    {
    case MODEL::HETEROZYGOUS:
        read_score(start_index, end_bound, 0, 1, 0, region_index, set_zero);
        break;
    case MODEL::DOMINANT:
        read_score(start_index, end_bound, 0, 1, 1, region_index, set_zero);
        break;
    case MODEL::RECESSIVE:
        read_score(start_index, end_bound, 0, 0, 1, region_index, set_zero);
        break;
    default:
        read_score(start_index, end_bound, 0, 1, 2, region_index, set_zero);
        break;
    }
}
