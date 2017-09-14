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


BinaryPlink::BinaryPlink(std::string prefix, std::string remove_sample,
                         std::string keep_sample, std::string extract_snp,
                         std::string exclude_snp, std::string fam_name,
                         std::string log_file, bool ignore_fid, bool nonfounder, int num_auto,
                         bool no_x, bool no_y, bool no_xy, bool no_mt,
                         bool keep_ambig, const size_t thread, bool verbose)
{
	m_nonfounder = nonfounder;
    m_fam_name = fam_name;
    /** simple assignments **/
    m_log_file = log_file;
    filter.keep_ambig = keep_ambig;
    m_thread = thread;
    // get the exclusion and extraction list
    if (!remove_sample.empty()) {
        m_sample_selection_list = load_ref(remove_sample, ignore_fid);
    }
    if (!keep_sample.empty()) {
        m_remove_sample = false;
        m_sample_selection_list = load_ref(keep_sample, ignore_fid);
    }
    if (!extract_snp.empty()) {
        m_exclude_snp = false;
        m_snp_selection_list = load_snp_list(extract_snp);
    }
    if (!exclude_snp.empty()) {
        m_snp_selection_list = load_snp_list(exclude_snp);
    }


    /** setting the chromosome information **/
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    // we are not using the following script for now as we only support human
    m_haploid_mask.resize(CHROM_MASK_WORDS, 0);
    m_chrom_mask.resize(CHROM_MASK_WORDS, 0);

    // now initialize the chromosome
    init_chr(num_auto, no_x, no_y, no_xy, no_mt);


    /** now get the chromosome information we've got by replacing the # in the
     * name **/
    // if there are multiple #, they will all be replaced by the same number
    set_genotype_files(prefix);


    /** now read the sample information **/
    m_sample_names = load_samples(ignore_fid);

    /** now read the SNP information **/
    m_existed_snps = load_snps();
    m_marker_ct = m_existed_snps.size();

    if (verbose) {
        std::ofstream log_file_stream;
        log_file_stream.open(log_file.c_str(), std::ofstream::app);
        if (!log_file_stream.is_open()) {
            std::string error_message =
                "ERROR: Cannot open log file: " + log_file;
            throw std::runtime_error(error_message);
        }
        fprintf(stderr, "%zu people (%zu males, %zu females) observed\n",
                m_unfiltered_sample_ct, m_num_male, m_num_female);
        fprintf(stderr, "%zu founder(s) included\n", m_founder_ct);
        log_file_stream << m_unfiltered_sample_ct << " people (" << m_num_male
                        << " male(s), " << m_num_female
                        << " female(s)) observed" << std::endl;
        log_file_stream << m_founder_ct << " founder(s) included" << std::endl;
        if (m_num_ambig != 0 && !keep_ambig) {
            fprintf(stderr, "%u ambiguous variant(s) excluded\n", m_num_ambig);
            log_file_stream << m_num_ambig << " ambiguous variant(s) excluded"
                            << std::endl;
        }
        else if (m_num_ambig != 0)
        {
            fprintf(stderr, "%u ambiguous variants kept\n", m_num_ambig);
            log_file_stream << m_num_ambig << " ambiguous variant(s) kept"
                            << std::endl;
        }
        fprintf(stderr, "%zu variants included\n", m_marker_ct);
        log_file_stream << m_marker_ct << " variant(s) included" << std::endl;
        log_file_stream << std::endl;
        log_file_stream.close();
    }

    check_bed();
    m_cur_file = "";

    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    m_sample_selection_list.clear();
    m_snp_selection_list.clear();
}

BinaryPlink::~BinaryPlink()
{
    if (m_bedfile != nullptr) {
        fclose(m_bedfile);
        m_bedfile = nullptr;
    }
}


std::vector<Sample> BinaryPlink::load_samples(bool ignore_fid)
{
    assert(m_genotype_files.size() > 0);
    // get the name of the first fam file (we only need the first as they should
    // all contain the same information)
    std::string famName = "";
    if (!m_fam_name.empty())
        famName = m_fam_name;
    else
        m_genotype_files.front() + ".fam";
    // open the fam file
    std::ifstream famfile;
    famfile.open(famName.c_str());
    if (!famfile.is_open()) {
        std::string error_message = "ERROR: Cannot open fam file: " + famName;
        throw std::runtime_error(error_message);
    }
    // number of unfiltered samples
    m_unfiltered_sample_ct = 0;
    std::string line;
    std::unordered_set<std::string> founder_info;
    // first pass to get the number of samples and also get the founder ID
    while (std::getline(famfile, line)) {
        misc::trim(line);
        if (!line.empty()) {
            std::vector<std::string> token = misc::split(line);
            if (token.size() < 6) {
                fprintf(stderr,
                        "Error: Malformed fam file. Less than 6 column on "
                        "line: %zu\n",
                        m_unfiltered_sample_ct + 1);
                throw std::runtime_error("");
            }
            founder_info.insert(token[+FAM::FID] + "_" + token[+FAM::IID]);
            m_unfiltered_sample_ct++;
        }
    }
    // now reset the fam file to the start
    famfile.clear();
    famfile.seekg(0);

    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);

    // we don't work with the sex for now, so better ignore them first
    // m_sex_male = new uintptr_t[unfiltered_sample_ctl];
    // std::fill(m_sex_male, m_sex_male+unfiltered_sample_ctl, 0);
    // std::memset(m_sex_male, 0x0, unfiltered_sample_ctl*sizeof(uintptr_t));

    // try to use fill instead of memset for better readability (will be tiny
    // bit slower according to stackoverflow)
    m_founder_info.resize(unfiltered_sample_ctl, 0);

    // Initialize this, but will copy founder into this later on
    m_sample_include.resize(unfiltered_sample_ctl, 0);

    m_num_male = 0, m_num_female = 0, m_num_ambig_sex = 0,
    m_num_non_founder = 0;
    std::vector<Sample> sample_name;
    uintptr_t sample_uidx = 0; // this is just for error message
    while (std::getline(famfile, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (token.size() < 6) {
            std::string error_message =
                "Error: Malformed fam file. Less than 6 column on line: "
                + std::to_string(sample_uidx + 1);
            throw std::runtime_error(error_message);
        }
        Sample cur_sample;
        cur_sample.FID = token[+FAM::FID];
        cur_sample.IID = token[+FAM::IID];
        std::string id = (ignore_fid)
                             ? token[+FAM::IID]
                             : token[+FAM::FID] + "_" + token[+FAM::IID];
        cur_sample.pheno = token[+FAM::PHENOTYPE];
        cur_sample.has_pheno = false; // only true when we have evaluated it to be true
        if (!m_remove_sample) {
            cur_sample.included = (m_sample_selection_list.find(id)
                                   != m_sample_selection_list.end());
        }
        else
        {
            cur_sample.included = (m_sample_selection_list.find(id)
                                   == m_sample_selection_list.end());
        }

        cur_sample.prs = 0;
        cur_sample.num_snp = 0;

        if (founder_info.find(token[+FAM::FATHER]) == founder_info.end()
            && founder_info.find(token[+FAM::MOTHER]) == founder_info.end()
            && cur_sample.included)
        {
            // only set this if no parents were found in the fam file
            m_founder_ct++;
            SET_BIT(sample_uidx, m_founder_info.data()); // essentially, m_founder is a subset of m_sample_include
            SET_BIT(sample_uidx, m_sample_include.data());
        }
        else if (cur_sample.included && m_nonfounder)
        {
            // nonfounder but we want to keep it
        		SET_BIT(sample_uidx, m_sample_include.data());
            m_num_non_founder++;
        }
        else
        {
        		// nonfounder / unwanted sample
        		if(cur_sample.included) m_num_non_founder++;
        }
        if (token[+FAM::SEX].compare("1") == 0) {
            m_num_male++;
            // SET_BIT(sample_uidx, m_sex_male); // if that individual is male,
            // need to set bit
        }
        else if (token[+FAM::SEX].compare("2") == 0)
        {
            m_num_female++;
        }
        else
        {
            m_num_ambig_sex++; // currently ignore as we don't do sex chromosome
            // SET_BIT(sample_uidx, m_sample_exclude); // exclude any samples
            // without sex information
        }
        sample_uidx++;
        sample_name.push_back(cur_sample);
    }
    famfile.close();
    return sample_name;
}

std::vector<SNP> BinaryPlink::load_snps()
{
    assert(m_genotype_files.size() > 0);
    m_unfiltered_marker_ct = 0;
    std::ifstream bimfile;
    std::string prev_chr = "";
    int chr_code = 0;
    int chr_index = 0;
    bool chr_error = false, chr_sex_error = false;
    m_num_ambig = 0;
    std::vector<SNP> snp_info;
    m_num_snp_per_file.resize(m_genotype_files.size());
    size_t cur_file = 0;
    for (auto&& prefix : m_genotype_files) {
        std::string bimname = prefix + ".bim";
        bimfile.open(bimname.c_str());
        if (!bimfile.is_open()) {
            std::string error_message =
                "Error: Cannot open bim file: " + bimname;
            throw std::runtime_error(error_message);
        }
        std::string line;
        int num_line = 0;
        while (std::getline(bimfile, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            std::vector<std::string> token = misc::split(line);
            if (token.size() < 6) {
                fprintf(stderr,
                        "Error: Malformed bim file. Less than 6 column on "
                        "line: %i\n",
                        num_line);
                throw std::runtime_error("");
            }
            // change them to upper case
            std::transform(token[+BIM::A1].begin(), token[+BIM::A1].end(),
                           token[+BIM::A1].begin(), ::toupper);
            std::transform(token[+BIM::A2].begin(), token[+BIM::A2].end(),
                           token[+BIM::A2].begin(), ::toupper);
            std::string chr = token[+BIM::CHR];
            // exclude SNPs that are not required
            if (!m_exclude_snp
                && m_snp_selection_list.find(token[+BIM::RS])
                       == m_snp_selection_list.end())
            {
                m_unfiltered_marker_ct++;
                m_num_snp_per_file[cur_file]++;
                num_line++;
                continue;
            }
            else if (m_exclude_snp
                     && m_snp_selection_list.find(token[+BIM::RS])
                            != m_snp_selection_list.end())
            {
                m_unfiltered_marker_ct++;
                m_num_snp_per_file[cur_file]++;
                num_line++;
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
                        fprintf(stderr,
                                "WARNING: SNPs with chromosome number larger "
                                "than %du\n",
                                m_max_code);
                        fprintf(stderr, "         They will be ignored!\n");
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
            // now get other information of the SNP
            int loc = misc::convert<int>(token[+BIM::BP]);
            if (loc < 0) {
                fprintf(stderr, "ERROR: SNP with negative corrdinate: %s:%s\n",
                        token[+BIM::RS].c_str(), token[+BIM::BP].c_str());
                throw std::runtime_error(
                    "Please check you have the correct input");
            }
            // we really don't like duplicated SNPs right?
            // or a more gentle way will be to exclude the subsequent SNP
            if (m_existed_snps_index.find(token[+BIM::RS])
                != m_existed_snps_index.end())
            {
                throw std::runtime_error(
                    "ERROR: Duplicated SNP ID detected!\n");
            }
            else if (ambiguous(token[+BIM::A1], token[+BIM::A2]))
            {
                m_num_ambig++;
                if (filter.keep_ambig) {
                    // keep it if user want to
                    m_existed_snps_index[token[+BIM::RS]] = snp_info.size();
                    snp_info.push_back(SNP(token[+BIM::RS], chr_code, loc,
                                           token[+BIM::A1], token[+BIM::A2],
                                           prefix, num_line));
                }
            }
            else
            {
                m_existed_snps_index[token[+BIM::RS]] = snp_info.size();
                snp_info.push_back(SNP(token[+BIM::RS], chr_code, loc,
                                       token[+BIM::A1], token[+BIM::A2], prefix,
                                       num_line));
            }
            m_unfiltered_marker_ct++; // add in the checking later on
            m_num_snp_per_file[cur_file]++;
            num_line++;
        }
        bimfile.close();
        cur_file++;
    }

    if (m_unfiltered_marker_ct > 2147483645) {
        throw std::runtime_error(
            "Error: PLINK does not suport more than 2^31 -3 variants. "
            "As we are using PLINK for some of our functions, we might "
            "encounter problem too. "
            "Sorry.");
    }
    return snp_info;
}

void BinaryPlink::check_bed()
{
    uint32_t uii = 0;
    int64_t llxx = 0;
    int64_t llyy = 0;
    int64_t llzz = 0;
    size_t cur_file = 0;

    uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    for (auto&& prefix : m_genotype_files) {
        std::string bedname = prefix + ".bed";
        m_bedfile = fopen(bedname.c_str(), FOPEN_RB);
        if (fseeko(m_bedfile, 0, SEEK_END)) {
            std::string error_message = "Cannot read bed file: " + bedname;
            throw std::runtime_error(error_message);
        }
        llxx = ftello(m_bedfile);
        if (!llxx) {
            throw std::runtime_error("Error: Empty .bed file.");
        }
        rewind(m_bedfile);
        // will let the g_textbuf stay for now
        char version_check[3];
        uii = fread(version_check, 1, 3, m_bedfile);
        size_t marker_ct = m_num_snp_per_file[cur_file];
        llyy = ((uint64_t) unfiltered_sample_ct4) * marker_ct;
        llzz = ((uint64_t) m_unfiltered_sample_ct) * ((marker_ct + 3) / 4);
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
        fclose(m_bedfile);
        m_bedfile = nullptr;
        cur_file++;
    }
}


void BinaryPlink::read_score(misc::vec2d<Sample_lite>& current_prs_score,
                             size_t start_index, size_t end_bound)
{
    uintptr_t final_mask = get_final_mask(m_sample_ct);
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uintptr_t ulii = 0;
    size_t num_included_samples = current_prs_score.cols();
    // m_final_mask

    // a lot of code in PLINK is to handle the sex chromosome
    // which suggest that PRS can be done on sex chromosome
    // that should be something later
    m_cur_file = ""; // just close it
    if (m_bedfile != nullptr) {
        fclose(m_bedfile);
        m_bedfile = nullptr;
    }

    // haven't figured out what this max_reverse does
    // a simple size check
    if (current_prs_score.rows() != m_region_size) {
        throw std::runtime_error(
            "Size of Matrix doesn't match number of region!!");
    }
    std::vector<bool> in_region(m_region_size);
    // index is w.r.t. partition, which contain all the information
    uintptr_t* genotype = new uintptr_t[unfiltered_sample_ctl * 2];
    for (size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
    { // for each SNP
        if (m_cur_file.empty()
            || m_cur_file.compare(m_existed_snps[i_snp].file_name()) != 0)
        {
            if (m_bedfile != nullptr) {
                fclose(m_bedfile);
                m_bedfile = nullptr;
            }
            m_cur_file = m_existed_snps[i_snp].file_name();
            std::string bedname = m_cur_file + ".bed";
            m_bedfile = fopen(
                bedname.c_str(),
                FOPEN_RB); // again, we are assuming that the file is correct
        }
        for (size_t i_region = 0; i_region < m_region_size; ++i_region) {
            in_region[i_region] = m_existed_snps[i_snp].in(i_region);
        }
        size_t cur_line = m_existed_snps[i_snp].snp_id();
        if (fseeko(m_bedfile,
                   m_bed_offset
                       + (cur_line * ((uint64_t) unfiltered_sample_ct4)),
                   SEEK_SET))
        {
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        // loadbuf_raw is the temporary
        // loadbuff is where the genotype will be located
        std::fill(genotype, genotype + unfiltered_sample_ctl * 2, 0);
        std::fill(m_tmp_genotype.begin(), m_tmp_genotype.end(), 0);
        if (load_and_collapse_incl(m_unfiltered_sample_ct, m_sample_ct,
                                   m_sample_include.data(), final_mask, false,
                                   m_bedfile, m_tmp_genotype.data(), genotype))
        {
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        uintptr_t* lbptr = genotype;
        uii = 0;
        std::vector<size_t> missing_samples;
        std::vector<double> sample_genotype(num_included_samples);
        double stat = m_existed_snps[i_snp].stat();
        bool flipped = m_existed_snps[i_snp].is_flipped();
        uint32_t sample_idx = 0;
        size_t total_num = 0;
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
                if (ukk == 1 || ukk == 3) // Because 01 is coded as missing
                {
                    // 3 is homo alternative
                    // int flipped_geno = snp_list[snp_index].geno(ukk);
                    if (sample_idx < num_included_samples) {
                        total_num += (ukk == 3) ? 2 : ukk;
                        sample_genotype[sample_idx] = (ukk == 3) ? 2 : ukk;
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

        double maf = ((double) total_num
                      / ((double) (num_included_samples - nmiss)
                         * 2.0)); // MAF does not count missing
        if (flipped) maf = 1.0 - maf;
        double center_score = stat * maf;
        size_t num_miss = missing_samples.size();
        size_t i_missing = 0;
        for (size_t i_sample = 0; i_sample < num_included_samples; ++i_sample) {
            if (i_missing < num_miss && i_sample == missing_samples[i_missing])
            {
                for (size_t i_region = 0; i_region < m_region_size; ++i_region)
                {
                    if (in_region[i_region]) {
                        if (m_scoring == SCORING::MEAN_IMPUTE)
                            current_prs_score(i_region, i_sample).prs +=
                                center_score;
                        if (m_scoring != SCORING::SET_ZERO)
                            current_prs_score(i_region, i_sample).num_snp++;
                    }
                }
                i_missing++;
            }
            else
            { // not missing sample
                for (size_t i_region = 0; i_region < m_region_size; ++i_region)
                {
                    if (in_region[i_region]) {
                        if (m_scoring == SCORING::CENTER) {
                            // if centering, we want to keep missing at 0
                            current_prs_score(i_region, i_sample).prs -=
                                center_score;
                        }
                        int g = (flipped) ? fabs(sample_genotype[i_sample] - 2)
                                          : sample_genotype[i_sample];
                        current_prs_score(i_region, i_sample).prs +=
                            g * stat * 0.5;
                        current_prs_score(i_region, i_sample).num_snp++;
                    }
                }
            }
        }
    }
    delete[] genotype;
}
