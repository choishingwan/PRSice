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
                         const std::string& sample_file, const size_t thread,
                         const bool ignore_fid, const bool keep_nonfounder,
                         const bool keep_ambig)
    : Genotype(thread, ignore_fid, keep_nonfounder, keep_ambig)
{
    // place holder. Currently set default to human.
    /** setting the chromosome information **/
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    // we are not using the following script for now as we only support human
    m_haploid_mask.resize(CHROM_MASK_WORDS, 0);
    m_chrom_mask.resize(CHROM_MASK_WORDS, 0);
    init_chr();
    // get the bed file names
    m_genotype_files = set_genotype_files(prefix);
    m_sample_file =
        sample_file.empty() ? m_genotype_files.front() + ".fam" : sample_file;
}


std::vector<Sample> BinaryPlink::gen_sample_vector()
{
    assert(m_genotype_files.size() > 0);
    // open the fam file
    std::ifstream famfile;
    famfile.open(m_sample_file.c_str());
    if (!famfile.is_open()) {
        std::string error_message =
            "ERROR: Cannot open fam file: " + m_sample_file;
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
    std::vector<Sample> sample_name;
    std::unordered_set<std::string> duplicated_samples;
    std::vector<std::string> duplicated_sample_id;
    uintptr_t sample_index = 0; // this is just for error message
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
        Sample cur_sample;
        cur_sample.FID = token[+FAM::FID];
        cur_sample.IID = token[+FAM::IID];
        std::string id = (m_ignore_fid)
                             ? token[+FAM::IID]
                             : token[+FAM::FID] + "_" + token[+FAM::IID];
        cur_sample.pheno = token[+FAM::PHENOTYPE];
        // false as we have not check if the pheno information is valid
        cur_sample.has_pheno = false;
        if (!m_remove_sample) {
            cur_sample.included = (m_sample_selection_list.find(id)
                                   != m_sample_selection_list.end());
        }
        else
        {
            cur_sample.included = (m_sample_selection_list.find(id)
                                   == m_sample_selection_list.end());
        }

        if (founder_info.find(token[+FAM::FATHER]) == founder_info.end()
            && founder_info.find(token[+FAM::MOTHER]) == founder_info.end()
            && cur_sample.included)
        {
            // only set this if no parents were found in the fam file
            m_founder_ct++;
            // so m_founder_info is a subset of m_sample_include
            SET_BIT(sample_index, m_founder_info.data());
        }
        else if (cur_sample.included && m_keep_nonfounder)
        {
            // nonfounder but we want to keep it
            SET_BIT(sample_index, m_sample_include.data());
            m_num_non_founder++;
        }
        else
        {
            // nonfounder / unwanted sample
            if (cur_sample.included) {
                // user didn't specify they want the nonfounder
                // so ignore it
                cur_sample.included = false;
                m_num_non_founder++;
            }
        }
        m_sample_ct += cur_sample.included;
        if (cur_sample.included) SET_BIT(sample_index, m_sample_include.data());
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
        if (!cur_sample.included) {
            // try to reduce memory usage...
            cur_sample.FID = "";
            cur_sample.IID = "";
            cur_sample.pheno = "";
        }
        duplicated_samples.insert(id);
        sample_name.push_back(cur_sample);
    }
    if (!duplicated_sample_id.empty()) {
        // TODO: Produce a file containing id of all valid samples
        std::string error_message =
            "ERROR: A total of " + std::to_string(duplicated_sample_id.size())
            + " duplicated samples detected!\n";
        error_message.append(
            "Please ensure all samples have an unique identifier");
        throw std::runtime_error(error_message);
    }

    famfile.close();
    m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    return sample_name;
}


std::vector<SNP> BinaryPlink::gen_snp_vector(const double geno,
                                             const double maf,
                                             const double info,
                                             const double hard_threshold,
                                             const bool hard_coded,
                                             const std::string& out_prefix)
{
    std::unordered_set<std::string> duplicated_snp;
    std::vector<SNP> snp_info;
    std::string line;
    uintptr_t final_mask = get_final_mask(m_sample_ct);
    uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    int chr_index = 0;
    int chr_code = 0;
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    for (auto prefix : m_genotype_files) {
        std::string bim_name = prefix + ".bim";
        std::string bed_name = prefix + ".bed";
        std::ifstream bim(bim_name.c_str());
        if (!bim.is_open()) {
            std::string error_message =
                "ERROR: Cannot open bim file: " + bim_name;
            throw std::runtime_error(error_message);
        }
        // First pass, get the number of marker in bed & bim
        int num_snp_read = 0;
        std::string prev_chr = "";
        bool chr_error = false, chr_sex_error = false; // to limit error report
        while (std::getline(bim, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            std::vector<std::string> bim_info = misc::split(line);
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
        check_bed(bed_name, num_snp_read);

        std::ifstream bed(bed_name.c_str());
        if (!bed.is_open()) {
            std::string error_message =
                "ERROR: Cannot open bed file: " + bed_name;
            throw std::runtime_error(error_message);
        }
        bed.seekg(m_bed_offset, std::ios_base::beg);
        // now go through the bim & bed file and perform filtering
        num_snp_read = 0;
        int prev_snp_processed = 0;
        while (std::getline(bim, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            num_snp_read++;
            std::vector<std::string> bim_info = misc::split(line);
            // doesn't need to do the format check as we have already done it in
            // the previous pass change them to upper case to avoid match
            // problems
            std::transform(bim_info[+BIM::A1].begin(), bim_info[+BIM::A1].end(),
                           bim_info[+BIM::A1].begin(), ::toupper);
            std::transform(bim_info[+BIM::A2].begin(), bim_info[+BIM::A2].end(),
                           bim_info[+BIM::A2].begin(), ::toupper);
            std::string chr = bim_info[+BIM::CHR];
            // exclude SNPs that are not required
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
            /** check if this is from a new chromosome **/
            if (chr.compare(prev_chr) != 0) {
                // only work on this if this is a new chromosome
                prev_chr = chr;
                if (m_chr_order.find(chr) != m_chr_order.end()) {
                    throw std::runtime_error("ERROR: SNPs on the same "
                                             "chromosome must be clustered "
                                             "together!");
                }
                m_chr_order[chr] = chr_index++;
                // get the chromosome codes
                chr_code = get_chrom_code_raw(chr.c_str());
                if (((const uint32_t) chr_code) > m_max_code)
                { // bigger than the maximum code, ignore it
                    if (!chr_error) {
                        // only print this if an error isn't previously given
                        std::string error_message =
                            "WARNING: SNPs with chromosome number larger "
                            "than "
                            + std::to_string(m_max_code) + "."
                            + " They will be ignored!\n";
                        std::cerr
                            << error_message
                            << std::endl; // currently avoid passing in
                                          // reporter here so that I don't need
                        // to pass the reporter as a parameter
                        chr_error = true;
                        continue;
                    }
                    else if (!chr_sex_error
                             && (is_set(m_haploid_mask.data(), chr_code)
                                 || chr_code == m_xymt_codes[X_OFFSET]
                                 || chr_code == m_xymt_codes[Y_OFFSET]))
                    {
                        // we ignore Sex chromosomes and haploid chromosome

                        fprintf(stderr, "WARNING: Currently not support "
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
                        "ERROR: SNP with negative corrdinate: "
                        + bim_info[+BIM::RS] + ":" + bim_info[+BIM::BP] + "\n";
                    error_message.append(
                        "Please check you have the correct input");
                    throw std::runtime_error(error_message);
                }
            }
            catch (const std::runtime_error& er)
            {

                std::string error_message =
                    "ERROR: SNP with non-numeric corrdinate: "
                    + bim_info[+BIM::RS] + ":" + bim_info[+BIM::BP] + "\n";
                error_message.append("Please check you have the correct input");
                throw std::runtime_error(error_message);
            }
            if (m_existed_snps_index.find(bim_info[+BIM::RS])
                != m_existed_snps_index.end())
            {
                duplicated_snp.insert(bim_info[+BIM::RS]);
                // throw std::runtime_error(
                //    "ERROR: Duplicated SNP ID detected!\n");
            }
            else if (!ambiguous(bim_info[+BIM::A1], bim_info[+BIM::A2])
                     || m_keep_ambig)
            {

                // now read in the binary information and determine if we want
                // to keep this SNP
                // only do the filtering if we need to as my current
                // implementation isn't as efficient as PLINK
                if (num_snp_read - prev_snp_processed > 1) {
                    // skip unread lines
                    if (!bed.seekg(m_bed_offset
                                       + ((num_snp_read - 1)
                                          * ((uint64_t) unfiltered_sample_ct4)),
                                   std::ios_base::beg))
                    {
                        std::string error_message =
                            "ERROR: Cannot read the bed file(seek): "
                            + bed_name;
                        throw std::runtime_error(error_message);
                    }
                }
                prev_snp_processed = (num_snp_read - 1);
                // get the location of the SNP in the binary file
                std::streampos byte_pos = bed.tellg();
                if (maf > 0 || geno < 1) {
                    if (load_and_collapse_incl(
                            m_unfiltered_sample_ct, m_sample_ct,
                            m_sample_include.data(), final_mask, false, bed,
                            m_tmp_genotype.data(), genotype.data()))
                    {
                        std::string error_message =
                            "ERROR: Cannot read the bed file(read): "
                            + bed_name;
                        throw std::runtime_error(error_message);
                    }
                    // Now genotype contain the genotype binary vector
                    uintptr_t* lbptr = genotype.data();
                    uint32_t uii = 0;
                    uintptr_t ulii = 0;
                    uint32_t ujj;
                    uint32_t ukk;
                    uint32_t sample_idx = 0;
                    int aa = 0, aA = 0, AA = 0;
                    size_t nmiss = 0;
                    do
                    {
                        ulii = ~(*lbptr++);
                        if (uii + BITCT2 > m_unfiltered_sample_ct) {
                            ulii &= (ONELU
                                     << ((m_unfiltered_sample_ct & (BITCT2 - 1))
                                         * 2))
                                    - ONELU;
                        }
                        while (ulii) {
                            ujj = CTZLU(ulii) & (BITCT - 2);
                            ukk = (ulii >> ujj) & 3;
                            sample_idx = uii + (ujj / 2);
                            if (ukk == 1
                                || ukk == 3) // Because 01 is coded as missing
                            {
                                // 3 is homo alternative
                                // int flipped_geno =
                                // snp_list[snp_index].geno(ukk);
                                if (sample_idx < m_sample_ct) {
                                    int g = (ukk == 3) ? 2 : ukk;
                                    switch (g)
                                    {
                                    case 0: aa++; break;
                                    case 1: aA++; break;
                                    case 2: AA++; break;
                                    }
                                }
                            }
                            else // this should be 2
                            {
                                nmiss++;
                            }
                            ulii &= ~((3 * ONELU) << ujj);
                        }
                        uii += BITCT2;
                    } while (uii < m_sample_ct);
                    // remove SNP if we have higher missingess than specified
                    if ((double) nmiss / (double) m_sample_ct > geno) {
                        m_num_geno_filter++;
                        continue;
                    }
                    double cur_maf = ((double) (aA + AA * 2)
                                      / ((double) (m_sample_ct - nmiss) * 2.0));
                    cur_maf = cur_maf > 0.5 ? 1 - cur_maf : cur_maf;
                    // remove SNP if maf lower than threshold
                    if (cur_maf < maf) {
                        m_num_maf_filter++;
                        continue;
                    }
                }

                m_num_ambig +=
                    ambiguous(bim_info[+BIM::A1], bim_info[+BIM::A2]);
                m_existed_snps_index[bim_info[+BIM::RS]] = snp_info.size();
                // TODO: When working with SNP class, we need to add in the aA
                // AA aa variable to avoid re-calculating the mean
                snp_info.emplace_back(
                    SNP(bim_info[+BIM::RS], chr_code, loc, bim_info[+BIM::A1],
                        bim_info[+BIM::A2], prefix, byte_pos));
            }
            else if (!m_keep_ambig)
            {
                m_num_ambig++;
            }
        }
    }

    if (duplicated_snp.size() != 0) {
        std::ofstream log_file_stream;
        std::string dup_name = out_prefix + ".valid";
        log_file_stream.open(dup_name.c_str());
        if (!log_file_stream.is_open()) {
            std::string error_message = "ERROR: Cannot open file: " + dup_name;
            throw std::runtime_error(error_message);
        }
        for (auto&& snp : m_existed_snps) {
            if (duplicated_snp.find(snp.rs()) != duplicated_snp.end()) continue;
            log_file_stream << snp.rs() << std::endl;
        }
        log_file_stream.close();
        std::string error_message =
            "ERROR: Duplicated SNP ID detected!.Valid SNP ID stored at "
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

void BinaryPlink::read_score(size_t start_index, size_t end_bound,
                             const size_t region_index)
{
    uintptr_t final_mask = get_final_mask(m_sample_ct);
    // for array size
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    size_t num_included_samples = m_sample_names.size();

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
                    "ERROR: Cannot open bed file: " + bedname;
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
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        m_prev_loc = cur_line + (std::streampos) unfiltered_sample_ct4;
        // loadbuf_raw is the temporary
        // loadbuff is where the genotype will be located
        if (load_and_collapse_incl(m_unfiltered_sample_ct, m_sample_ct,
                                   m_sample_include.data(), final_mask, false,
                                   m_bed_file, m_tmp_genotype.data(),
                                   genotype.data()))
        {
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }

        uintptr_t* lbptr = genotype.data();
        uint32_t uii = 0;
        uintptr_t ulii = 0;
        uint32_t ujj;
        uint32_t ukk;
        std::vector<size_t> missing_samples;
        std::vector<double> sample_genotype(num_included_samples);
        double stat = cur_snp.stat() * 2; // Multiply by ploidy
        bool flipped = cur_snp.is_flipped();
        uint32_t sample_idx = 0;

        int aa = 0, aA = 0, AA = 0;
        size_t nmiss = 0;
        do
        {
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct) {
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            while (ulii) {
                ujj = CTZLU(ulii) & (BITCT - 2);
                ukk = (ulii >> ujj) & 3;
                sample_idx = uii + (ujj / 2);
                if (ukk == 1 || ukk == 3) // Because 10 is coded as missing
                {
                    // 3 is homo alternative
                    // int flipped_geno = snp_list[snp_index].geno(ukk);
                    if (sample_idx < num_included_samples) {
                        int g = (ukk == 3) ? 2 : ukk;
                        switch (g)
                        {
                        case 0: aa++; break;
                        case 1: aA++; break;
                        case 2: AA++; break;
                        }
                        sample_genotype[sample_idx] = g;
                    }
                }
                else // this should be 2
                {
                    missing_samples.push_back(sample_idx);
                    nmiss++;
                }
                ulii &= ~((3 * ONELU) << ujj);
            }
            uii += BITCT2;
        } while (uii < num_included_samples);

        if (num_included_samples - nmiss == 0) {
            cur_snp.invalidate();
            continue;
        }
        // due to the way the binary code works, the aa will always be 0
        // added there just for fun tbh
        aa = num_included_samples - nmiss - aA - AA;
        assert(aa >= 0);
        if (flipped) {
            int temp = aa;
            aa = AA;
            AA = temp;
        }
        if (m_model == MODEL::HETEROZYGOUS) {
            // 010
            aa += AA;
            AA = 0;
        }
        else if (m_model == MODEL::DOMINANT)
        {
            // 011;
            aA += AA;
            AA = 0;
        }
        else if (m_model == MODEL::RECESSIVE)
        {
            // 001
            aa += aA;
            aA = AA;
            AA = 0;
        }
        double maf = ((double) (aA + AA * 2)
                      / ((double) (num_included_samples - nmiss)
                         * 2.0)); // MAF does not count missing
        double center_score = stat * maf;
        size_t num_miss = missing_samples.size();
        size_t i_missing = 0;
        // actual index should differ due to PLINK automatically remove samples
        // that are not included
        size_t actual_index = 0;
        for (size_t i_sample = 0; i_sample < num_included_samples; ++i_sample) {
            auto&& sample = m_sample_names[i_sample];
            if (!sample.included) continue;
            if (i_missing < num_miss
                && actual_index == missing_samples[i_missing])
            {
                if (m_missing_score == MISSING_SCORE::MEAN_IMPUTE)
                    sample.prs += center_score;
                // g_prs_storage[vector_pad + i_sample] += center_score;
                if (m_missing_score != MISSING_SCORE::SET_ZERO)
                    sample.num_snp++;
                // g_num_snps[vector_pad + i_sample]++;

                i_missing++;
            }
            else
            { // not missing sample
                if (m_missing_score == MISSING_SCORE::CENTER) {
                    // if centering, we want to keep missing at 0
                    sample.prs -= center_score;
                    // g_prs_storage[vector_pad + i_sample] -= center_score;
                }
                int g = (flipped) ? fabs(sample_genotype[actual_index] - 2)
                                  : sample_genotype[actual_index];
                if (m_model == MODEL::HETEROZYGOUS) {
                    g = (g == 2) ? 0 : g;
                }
                else if (m_model == MODEL::RECESSIVE)
                {
                    g = std::max(0, g - 1);
                }
                else if (m_model == MODEL::DOMINANT)
                {
                    g = (g == 2) ? 1 : g;
                }
                sample.prs += g * stat * 0.5;
                // g_prs_storage[vector_pad + i_sample] += g * stat * 0.5;
                sample.num_snp++;
                // g_num_snps[vector_pad + i_sample]++;
            }
            actual_index++;
        }
    }
}
