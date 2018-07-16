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

BinaryPlink::BinaryPlink(const std::string& prefix,
                         const std::string& sample_file,
                         const std::string& multi_input, const size_t thread,
                         const bool ignore_fid, const bool keep_nonfounder,
                         const bool keep_ambig, const bool is_ref)
    : Genotype(thread, ignore_fid, keep_nonfounder, keep_ambig, is_ref)
{
    // place holder. Currently set default to human.
    /** setting the chromosome information **/
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    // we are not using the following script for now as we only support human
    m_haploid_mask.resize(CHROM_MASK_WORDS, 0);
    // m_chrom_mask.resize(CHROM_MASK_WORDS, 0);
    init_chr();
    // get the bed file names
    if (multi_input.empty())
        m_genotype_files = set_genotype_files(prefix);
    else
        m_genotype_files = load_genotype_prefix(multi_input);
    m_sample_file =
        sample_file.empty() ? m_genotype_files.front() + ".fam" : sample_file;
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
    // the unfiltered_sampel_ct is used to define the size of all vector used
    // within the program
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);

    // Currently ignore sex information
    m_founder_info.resize(unfiltered_sample_ctl, 0);
    m_sample_include.resize(unfiltered_sample_ctl, 0);

    m_num_male = 0, m_num_female = 0, m_num_ambig_sex = 0,
    m_num_non_founder = 0;
    std::vector<Sample_ID> sample_name;
    std::unordered_set<std::string> duplicated_samples;
    std::vector<std::string> duplicated_sample_id;
    uintptr_t sample_index = 0; // this is just for error message
    bool inclusion = false;
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
        Sample_ID cur_sample;
        cur_sample.FID = token[+FAM::FID];
        cur_sample.IID = token[+FAM::IID];
        std::string id = (m_ignore_fid)
                             ? token[+FAM::IID]
                             : token[+FAM::FID] + "_" + token[+FAM::IID];
        cur_sample.pheno = token[+FAM::PHENOTYPE];
        // cur_sample.in_regression = false;
        // false as we have not check if the pheno information is valid
        if (!m_remove_sample) {
            inclusion = (m_sample_selection_list.find(id)
                         != m_sample_selection_list.end());
        }
        else
        {
            inclusion = (m_sample_selection_list.find(id)
                         == m_sample_selection_list.end());
        }

        if (founder_info.find(token[+FAM::FATHER]) == founder_info.end()
            && founder_info.find(token[+FAM::MOTHER]) == founder_info.end()
            && inclusion)
        {
            // only set this if no parents were found in the fam file
            m_founder_ct++;
            // so m_founder_info is a subset of m_sample_include
            SET_BIT(sample_index, m_founder_info.data());
            SET_BIT(sample_index, m_sample_include.data());
        }
        else if (inclusion)
        {
            SET_BIT(sample_index, m_sample_include.data());
            m_num_non_founder++;
            // cur_sample.founder = m_keep_nonfounder;
        }
        m_sample_ct += inclusion;
        if (token[+FAM::SEX].compare("1") == 0) {
            m_num_male++;
        }
        else if (token[+FAM::SEX].compare("2") == 0)
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
        // m_founder_info to indicate if sample is needed
        if (inclusion && !m_is_ref) {
            sample_name.push_back(cur_sample);
        }
        duplicated_samples.insert(id);
    }

    if (!duplicated_sample_id.empty()) {
        // TODO: Produce a file containing id of all valid samples
        std::string error_message =
            "Error: A total of " + std::to_string(duplicated_sample_id.size())
            + " duplicated samples detected!\n";
        error_message.append(
            "Please ensure all samples have an unique identifier");
        throw std::runtime_error(error_message);
    }

    famfile.close();
    m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    return sample_name;
}


std::vector<SNP>
BinaryPlink::gen_snp_vector(const double geno, const double maf,
                            const double info, const double hard_threshold,
                            const bool hard_coded, Region& exclusion,
                            const std::string& out_prefix, Genotype* target)
{
    std::unordered_set<std::string> duplicated_snp;
    std::vector<SNP> snp_info;


    std::string line;
    const uintptr_t final_mask = get_final_mask(m_sample_ct);
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    int chr_index = 0;
    int chr_code = 0;
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);


    // initialize the mask for read score
    const uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_sample_ct);
    uint32_t homrar_ct;
    uint32_t missing_ct;
    uint32_t het_ct;
    uint32_t homcom_ct;
    uint32_t num_ref_target_match = 0;
    intptr_t nanal;
    double cur_maf;
    bool chr_error = false, chr_sex_error = false,
         dummy; // to limit error report
    int num_snp_read = 0, prev_snp_processed = 0;
    std::string bim_name, bed_name, chr;
    std::ifstream bim, bed;
    std::string prev_chr = "";
    std::streampos byte_pos;
    std::vector<std::string> bim_info;
    m_sample_mask.resize(pheno_nm_ctv2);
    fill_quatervec_55(m_sample_ct, m_sample_mask.data());
    for (auto prefix : m_genotype_files) {
        bim_name = prefix + ".bim";
        bed_name = prefix + ".bed";
        bim.open(bim_name.c_str());
        if (!bim.is_open()) {
            std::string error_message =
                "Error: Cannot open bim file: " + bim_name;
            throw std::runtime_error(error_message);
        }
        // First pass, get the number of marker in bed & bim
        num_snp_read = 0;
        prev_chr = "";
        // first, we need the number of SNPs included so that
        // we can check the bit of the bim
        while (std::getline(bim, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            bim_info = misc::split(line);
            if (bim_info.size() < 6) {
                std::string error_message =
                    "Error: Malformed bim file. Less than 6 column on "
                    "line: "
                    + std::to_string(num_snp_read) + "\n";
                throw std::runtime_error(error_message);
            }
            num_snp_read++;
        }
        bim.clear();
        bim.seekg(0, bim.beg);
        // check if the bed file is valid
        check_bed(bed_name, num_snp_read);

        bed.open(bed_name.c_str());
        if (!bed.is_open()) {
            std::string error_message =
                "Error: Cannot open bed file: " + bed_name;
            throw std::runtime_error(error_message);
        }
        // bed.seekg(m_bed_offset, std::ios_base::beg);
        // now go through the bim & bed file and perform filtering
        num_snp_read = 0;
        prev_snp_processed = 0;
        while (std::getline(bim, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            num_snp_read++;
            bim_info = misc::split(line);
            // doesn't need to do the format check as we have already done it in
            // the previous pass change them to upper case to avoid match
            // problems
            if (m_is_ref) {
                // SNP not found in the target file
                if (target->m_existed_snps_index.find(bim_info[+BIM::RS])
                    == target->m_existed_snps_index.end())
                {
                    continue;
                }
            }
            std::transform(bim_info[+BIM::A1].begin(), bim_info[+BIM::A1].end(),
                           bim_info[+BIM::A1].begin(), ::toupper);
            std::transform(bim_info[+BIM::A2].begin(), bim_info[+BIM::A2].end(),
                           bim_info[+BIM::A2].begin(), ::toupper);
            chr = bim_info[+BIM::CHR];
            // exclude SNPs that are not required
            if (!m_is_ref) {
                if (!m_exclude_snp
                    && m_snp_selection_list.find(bim_info[+BIM::RS])
                           == m_snp_selection_list.end())
                {
                    continue;
                }
                else if (m_exclude_snp
                         && m_snp_selection_list.find(bim_info[+BIM::RS])
                                != m_snp_selection_list.end())
                {
                    continue;
                }
            }
            /** check if this is from a new chromosome **/
            if (chr.compare(prev_chr) != 0) {
                prev_chr = chr;
                if (m_chr_order.find(chr) != m_chr_order.end()) {
                    throw std::runtime_error("Error: SNPs on the same "
                                             "chromosome must be clustered "
                                             "together!");
                }
                // m_chr_order is to provide consistent sorting for the
                // region class
                m_chr_order[chr] = chr_index++;
                // get the chromosome codes
                chr_code = get_chrom_code_raw(chr.c_str());
                if (((const uint32_t) chr_code) > m_max_code) {
                    // bigger than the maximum code, ignore it
                    if (!chr_error) {
                        // only print this if an error isn't previously
                        // given
                        std::string error_message =
                            "Warning: SNPs with chromosome number larger "
                            "than "
                            + std::to_string(m_max_code) + "."
                            + " They will be ignored!\n";
                        std::cerr << error_message << "\n";
                        // currently avoid passing in reporter here so that  I
                        // don't need to pass the reporter as a parameter
                        chr_error = true;
                        continue;
                    }
                    else if (!chr_sex_error
                             && (is_set(m_haploid_mask.data(), chr_code)
                                 || chr_code == m_xymt_codes[X_OFFSET]
                                 || chr_code == m_xymt_codes[Y_OFFSET]))
                    {
                        // we ignore Sex chromosomes and haploid chromosome

                        fprintf(stderr, "Warning: Currently not support "
                                        "haploid chromosome and sex "
                                        "chromosomes\n");
                        chr_sex_error = true;
                        continue;
                    }
                }
            }

            // now read in the coordinate
            int loc = -1;
            try
            {
                loc = misc::convert<int>(bim_info[+BIM::BP]);
                if (loc < 0) {
                    std::string error_message =
                        "Error: SNP with negative corrdinate: "
                        + bim_info[+BIM::RS] + ":" + bim_info[+BIM::BP] + "\n";
                    error_message.append(
                        "Please check you have the correct input");
                    throw std::runtime_error(error_message);
                }
            }
            catch (const std::runtime_error& er)
            {

                std::string error_message =
                    "Error: SNP with non-numeric corrdinate: "
                    + bim_info[+BIM::RS] + ":" + bim_info[+BIM::BP] + "\n";
                error_message.append("Please check you have the correct input");
                throw std::runtime_error(error_message);
            }

            if (exclusion.check_exclusion(chr, loc)) {
                continue;
            }

            if (m_existed_snps_index.find(bim_info[+BIM::RS])
                != m_existed_snps_index.end())
            {
                duplicated_snp.insert(bim_info[+BIM::RS]);
            }
            else if (!ambiguous(bim_info[+BIM::A1], bim_info[+BIM::A2])
                     || m_keep_ambig)
            {
                // now read in the binary information and determine if we want
                // to keep this SNP
                // only do the filtering if we need to as my current
                // implementation isn't as efficient as PLINK
                byte_pos =
                    m_bed_offset
                    + ((num_snp_read - 1) * ((uint64_t) unfiltered_sample_ct4));
                if (maf > 0 || geno < 1) {
                    if (num_snp_read - prev_snp_processed > 1) {
                        // skip unread lines
                        if (!bed.seekg(byte_pos, std::ios_base::beg)) {
                            std::string error_message =
                                "Error: Cannot read the bed file(seek): "
                                + bed_name;
                            throw std::runtime_error(error_message);
                        }
                    }
                    prev_snp_processed = (num_snp_read - 1);
                    if (load_and_collapse_incl(
                            m_unfiltered_sample_ct, m_sample_ct,
                            m_sample_include.data(), final_mask, false, bed,
                            m_tmp_genotype.data(), genotype.data()))
                    {
                        std::string error_message =
                            "Error: Cannot read the bed file(read): "
                            + bed_name;
                        throw std::runtime_error(error_message);
                    }
                    genovec_3freq(genotype.data(), m_sample_mask.data(),
                                  pheno_nm_ctv2, &missing_ct, &het_ct,
                                  &homcom_ct);
                    nanal = m_sample_ct - missing_ct;
                    homrar_ct = nanal - het_ct - homcom_ct;
                    if (nanal == 0) {
                        // still count as MAF filtering (for now)
                        m_num_maf_filter++;
                        continue;
                    }

                    if ((double) missing_ct / (double) m_sample_ct > geno) {
                        m_num_geno_filter++;
                        continue;
                    }

                    cur_maf = ((double) (het_ct + homrar_ct * 2)
                               / ((double) nanal * 2.0));
                    if (cur_maf > 0.5) cur_maf = 1.0 - cur_maf;
                    // remove SNP if maf lower than threshold
                    if (cur_maf < maf) {
                        m_num_maf_filter++;
                        continue;
                    }
                }
                m_num_ambig +=
                    ambiguous(bim_info[+BIM::A1], bim_info[+BIM::A2]);
                if (!m_is_ref) {
                    m_existed_snps_index[bim_info[+BIM::RS]] = snp_info.size();
                    // TODO: When working with SNP class, we need to add in the
                    // aA AA aa variable to avoid re-calculating the mean
                    snp_info.push_back(SNP(
                        bim_info[+BIM::RS], chr_code, loc, bim_info[+BIM::A1],
                        bim_info[+BIM::A2], prefix, byte_pos, homcom_ct, het_ct,
                        homrar_ct, missing_ct));
                }
                else
                {
                    auto&& target_index =
                        target->m_existed_snps_index[bim_info[+BIM::RS]];
                    if (!target->m_existed_snps[target_index].matching(
                            chr_code, loc, bim_info[+BIM::A1],
                            bim_info[+BIM::A2], dummy))
                    {
                        m_num_ref_target_mismatch++;
                    }
                    else
                    {
                        target->m_existed_snps[target_index].add_reference(
                            prefix, byte_pos);
                        num_ref_target_match++;
                    }
                }
            }
            else if (!m_keep_ambig)
            {
                m_num_ambig++;
            }
        }
    }
    snp_info.shrink_to_fit();
    if (m_is_ref && num_ref_target_match != target->m_existed_snps.size()) {
        target->update_snps();
    }
    if (duplicated_snp.size() != 0) {
        std::ofstream log_file_stream;
        std::string dup_name = out_prefix + ".valid";
        log_file_stream.open(dup_name.c_str());
        if (!log_file_stream.is_open()) {
            std::string error_message = "Error: Cannot open file: " + dup_name;
            throw std::runtime_error(error_message);
        }
        for (auto&& snp : snp_info) {
            if (duplicated_snp.find(snp.rs()) != duplicated_snp.end()) continue;
            log_file_stream << snp.rs() << "\n";
        }
        log_file_stream.close();
        std::string error_message =
            "Error: A total of " + std::to_string(duplicated_snp.size())
            + " duplicated SNP ID detected out of "
            + std::to_string(snp_info.size())
            + " input SNPs! Valid SNP ID (post --extract / "
              "--exclude, non-duplicated SNPs) stored at "
            + dup_name + ". You can avoid this error by using --extract "
            + dup_name;
        throw std::runtime_error(error_message);
    }

    return snp_info;
}

void BinaryPlink::check_bed(const std::string& bed_name, size_t num_marker)
{
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
    uii = bed.gcount();
    llyy = ((uint64_t) unfiltered_sample_ct4) * num_marker;
    llzz = ((uint64_t) m_unfiltered_sample_ct) * ((num_marker + 3) / 4);
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
        m_bed_offset = 2;
    }
    else if (uii && (*version_check == '\x01'))
    {
        // v0.99 SNP-major
        llyy += 1;
        m_bed_offset = 1;
    }
    else if (uii && (!(*version_check)))
    {
        // v0.99 sample-major
        sample_major = true;
        llyy = llzz + 1;
        m_bed_offset = 2;
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
        m_bed_offset = 2;
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
        genovec_3freq(genotype.data(), sample_include2.data(), pheno_nm_ctv2,
                      &missing_ct, &het_ct, &homcom_ct);
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
        genovec_3freq(genotype.data(), sample_include2.data(), pheno_nm_ctv2,
                      &missing_ct, &het_ct, &homcom_ct);
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
