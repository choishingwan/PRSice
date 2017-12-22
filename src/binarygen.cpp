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

#include "binarygen.hpp"

namespace bgenlib = genfile::bgen;

BinaryGen::BinaryGen(const std::string& prefix, const std::string& sample_file,
                     const size_t thread, const bool ignore_fid,
                     const bool keep_nonfounder, const bool keep_ambig)
    : m_thread = thread
, m_ignore_fid = ignore_fid
, m_keep_nonfounder = keep_nonfounder
, m_keep_ambig = keep_ambig
{
    /** setting the chromosome information **/
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    // we are not using the following script for now as we only support human
    m_haploid_mask.resize(CHROM_MASK_WORDS, 0);
    m_chrom_mask.resize(CHROM_MASK_WORDS, 0);
    // place holder. Currently set default to human.
    init_chr();
    // get the bed file names
    m_genotype_files = set_genotype_files(prefix);
    m_sample_file =
        sample_file.empty() ? m_genotype_files.front() : sample_file;
}

std::vector<Sample> BinaryGen::gen_sample_vector(Reporter& reporter)
{
    bool is_sample_format = is_sample_format();
    std::ifstream sample_file(m_sample_file.c_str());
    if (!sample_file.is_open()) {
        std::string error_message =
            "ERROR: Cannot open sample file: " + m_sample_file;
        throw std::runtime_error(error_message)
    }
    std::string line;
    int sex_col = -1;
    if (is_sample_format) {
        std::getline(sample_file, line);
        std::vector<std::string> header_names = misc::split(line);
        std::getline(sample_file, line);
        reporter.report("Detected bgen sample file format");
        for (size_t i = 3; i < header_names.size(); ++i) {
            std::transform(header_names[i].begin(), header_names[i].end(),
                           header_names[i].begin(), ::toupper);
            if (header_names[i].compare("SEX") == 0) {
                sex_col = i;
                break;
            }
        }
        if (sex_col != -1) {
            // double check if the format is alright
            std::vector<std::string> header_format = misc::split(line);
            if (header_format[sex_col].compare("D") != 0) {
                std::string error_message =
                    "ERROR: Sex must be coded as \"D\" in bgen sample file!";
                throw std::runtime_error(error_message);
            }
        }
    }

    int line_id = 0;
    std::vector<Sample> sample_name;
    std::vector<int> sex;
    while (std::getline(sample_file, line)) {
        misc::trim(line);
        if (!line.empty()) {
            std::vector<std::string> token = misc::split(line);
            if (token.size()
                < ((sex_col != -1) ? (sex_col) : (1 + !ignore_fid)))
            {
                std::string error_message =
                    "ERROR: Line " + std::to_string(line_id)
                    + " must have at least "
                    + std::to_string((has_sex) ? (sex_col) : (1 + !ignore_fid))
                    + " columns! Number of column="
                    + std::to_string(token.size());
                throw std::runtime_error(error_message);
            }
            m_unfiltered_sample_ct++;
            Sample cur_sample;
            if (is_sample_format || !m_ignore_fid) {
                cur_sample.FID = token[0];
                cur_sample.IID = token[1];
            }
            else
            {
                cur_sample.FID = "";
                cur_sample.IID = token[1];
            }
            std::string id = (m_ignore_fid)
                                 ? cur_sample.IID
                                 : cur_sample.FID + "_" + cur_sample.IID;
            cur_sample.pheno = "";
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
            if (sex_col != -1) {
                try
                {
                    sex.push_back(misc::convert<int>(token[sex_col]));
                }
                catch (const std::runtime_error& er)
                {
                    throw std::runtime_error("ERROR: Invalid sex coding!\n");
                }
            }
            sample_name.push_back(cur_sample);
        }
    }


    m_unfiltered_sample_ct = sample_name.size();
    m_founder_ct = m_unfiltered_sample_ct;
    // now set all those vectors
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    // don't bother with founder info here as we don't have this information
    m_sample_include.resize(unfiltered_sample_ctl, 0);

    m_num_male = 0, m_num_female = 0, m_num_ambig_sex = 0,
    m_num_non_founder = 0;
    for (size_t i = 0; i < sample_name; i++) {
        if (sample_name[i].included) SET_BIT(i, m_sample_include.data());
        if (sex_col != -1) {
            m_num_male += (sex[i] == 1);
            m_num_female += (sex[i] == 2);
            m_num_ambig_sex += (sex[i] != 1 && sex[i] != 2);
        }
        else
        {
            m_num_ambig_sex++;
        }
    }
    sample_file.close();
    return sample_name;
}

bool BinaryGen::is_sample_format()
{
    std::ifstream sample_file(m_sample_file.c_str());
    if (!sample_file.is_open()) {
        std::string error_message =
            "ERROR: Cannot open sample file: " + m_sample_file;
        throw std::runtime_error(error_message)
    }
    std::string first_line, second_line;
    std::getline(sample_file, first_line);
    std::getline(sample_file, second_line);
    sample_file.close();
    std::vector<std::string> first_row = misc::split(first_line);
    std::vector<std::string> second_row = misc::split(second_line);
    if (first_row.size() != second_row.size() || first_row.size() < 3) {
        return false;
    }
    // first 3 must be 0
    for (size_t i = 0; i < 3; ++i) {
        if (second_row[i].compare("0") != 0) return false;
    }
    // DCPB
    for (size_t i = 4; i < second_row.size(); ++i) {
        if (second_row[i].length() > 1) return false;
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


std::vector<SNP> BinaryGen::gen_snp_vector(const double geno, const double maf,
                                           const double info,
                                           const double hard_threshold,
                                           const bool hard_coded,
                                           const std::string& out_prefix)
{
}

BinaryGen::BinaryGen(std::string prefix, std::string pheno_file, bool header,
                     std::string remove_sample, std::string keep_sample,
                     std::string extract_snp, std::string exclude_snp,
                     const std::string& out_prefix, Reporter& reporter,
                     bool ignore_fid, int num_auto, bool no_x, bool no_y,
                     bool no_xy, bool no_mt, bool keep_ambig,
                     const size_t thread, bool verbose)
{
    m_thread = thread;
    filter.keep_ambig = keep_ambig;
    if (!remove_sample.empty()) {
        m_sample_selection_list = load_ref(remove_sample, ignore_fid);
    }
    if (!keep_sample.empty()) {
        m_remove_sample = false;
        m_sample_selection_list = load_ref(keep_sample, ignore_fid);
    }
    if (!extract_snp.empty()) {
        m_exclude_snp = false;
        m_snp_selection_list = load_snp_list(extract_snp, reporter);
    }
    if (!exclude_snp.empty()) {
        m_snp_selection_list = load_snp_list(exclude_snp, reporter);
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
    set_genotype_files(prefix);

    m_sample_names = preload_samples(pheno_file, reporter, header, ignore_fid);
    // now we update the sample information accordingly
    // unlike binaryPlink, we don't actually want the return value
    // mainly because this is a QC step, checking consistency between
    // the bgen file and the pheno file if the bgen also contain
    // the sample information
    load_samples(ignore_fid);

    /** Load SNP information from the file **/
    m_existed_snps = load_snps(out_prefix);
    m_marker_ct = m_existed_snps.size();

    if (verbose) {
        std::string message = std::to_string(m_unfiltered_sample_ct)
                              + " people (" + std::to_string(m_num_male)
                              + " male(s), " + std::to_string(m_num_female)
                              + " female(s)) observed\n";
        message.append(std::to_string(m_founder_ct) + " founder(s) included\n");
        if (m_num_ambig != 0 && !keep_ambig) {
            message.append(std::to_string(m_num_ambig)
                           + " ambiguous variant(s) excluded\n");
        }
        else if (m_num_ambig != 0)
        {
            message.append(std::to_string(m_num_ambig)
                           + " ambiguous variant(s) kept\n");
        }
        message.append(std::to_string(m_marker_ct) + " variant(s) included\n");
        reporter.report(message);
    }

    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    m_sample_selection_list.clear();
    m_snp_selection_list.clear();
    m_bgen_file.close();
}

BinaryGen::~BinaryGen()
{
    if (m_bgen_file.is_open()) m_bgen_file.close();
}

std::vector<SNP> BinaryGen::load_snps(const std::string& out_prefix)
{
    std::vector<SNP> snp_res;
    bool chr_sex_error = false;
    bool chr_error = false;
    size_t chr_index = 0;
    size_t expected_total = 0;
    m_num_ambig = 0;
    m_unfiltered_marker_ct = 0;
    for (auto&& info : m_bgen_info) {
        expected_total += info.second.number_of_variants;
    }
    std::unordered_set<std::string> dup_list;
    snp_res.reserve(expected_total);
    for (auto&& prefix : m_genotype_files) {
        std::string bgen_name = prefix + ".bgen";
        if (m_bgen_file.is_open()) m_bgen_file.close();
        m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
        if (!m_bgen_file.is_open()) {
            std::string error_message =
                "ERROR: Cannot open bgen file " + bgen_name;
            throw std::runtime_error(error_message);
        }
        uint32_t offset = m_offset_map[prefix];
        m_bgen_file.seekg(offset + 4);
        uint32_t num_snp = m_bgen_info[prefix].number_of_variants;
        std::string prev_chr = "";
        int chr_code = 0;

        auto&& context = m_bgen_info[prefix];
        for (size_t i_snp = 0; i_snp < num_snp; ++i_snp) {
            if (m_unfiltered_marker_ct % 1000 == 0
                && m_unfiltered_marker_ct > 0)
            {
                fprintf(stderr, "\r%zuK SNPs processed\r",
                        m_unfiltered_marker_ct / 1000);
            }


            std::string allele;
            std::string SNPID;
            std::string RSID;
            std::string chromosome;
            uint32_t SNP_position;
            std::vector<std::string> alleles;
            // directly use the library
            genfile::bgen::read_snp_identifying_data(
                m_bgen_file, context, &SNPID, &RSID, &chromosome, &SNP_position,
                [&alleles](std::size_t n) { alleles.resize(n); },
                [&alleles](std::size_t i, std::string const& allele) {
                    alleles.at(i) = allele;
                });
            for (auto&& a : alleles) {
                std::transform(a.begin(), a.end(), a.begin(), ::toupper);
            }

            std::streampos byte_pos = m_bgen_file.tellg();
            bool exclude_snp = false;
            // but we will not process anything
            if (chromosome.compare(prev_chr) != 0) {
                prev_chr = chromosome;
                if (m_chr_order.find(chromosome) != m_chr_order.end()) {
                    throw std::runtime_error("ERROR: SNPs on the same "
                                             "chromosome must be clustered "
                                             "together!");
                }
                m_chr_order[chromosome] = chr_index++;
                chr_code = get_chrom_code_raw(chromosome.c_str());
                if (((const uint32_t) chr_code) > m_max_code)
                { // bigger than the maximum code, ignore it
                    if (!chr_error) {
                        fprintf(stderr,
                                "WARNING: SNPs with chromosome number larger "
                                "than %du\n",
                                m_max_code);
                        fprintf(stderr, "         They will be ignored!\n");
                        chr_error = true;
                        exclude_snp = true;
                    }
                    else if (!chr_sex_error
                             && (is_set(m_haploid_mask.data(), chr_code)
                                 || chr_code == m_xymt_codes[X_OFFSET]
                                 || chr_code == m_xymt_codes[Y_OFFSET]))
                    {
                        fprintf(stderr, "WARNING: Currently not support "
                                        "haploid chromosome and sex "
                                        "chromosomes\n");
                        chr_sex_error = true;
                        exclude_snp = true;
                    }
                }
            }

            if (RSID.compare(".") == 0) // when the rs id isn't available,
                                        // change it to chr:loc coding
            {
                RSID = std::to_string(chr_code) + ":"
                       + std::to_string(SNP_position);
            }
            if ((!m_exclude_snp
                 && m_snp_selection_list.find(RSID)
                        == m_snp_selection_list.end())
                || (m_exclude_snp
                    && m_snp_selection_list.find(RSID)
                           != m_snp_selection_list.end()))
            {
                m_unfiltered_marker_ct++;
                exclude_snp = true;
            }

            if (m_existed_snps_index.find(RSID) != m_existed_snps_index.end()) {
                dup_list.insert(RSID);
            }
            else if (ambiguous(alleles.front(), alleles.back()))
            {
                m_num_ambig++;
                if (!filter.keep_ambig) exclude_snp = true;
            }


            Data probability;
            ProbSetter setter(&probability);
            std::vector<genfile::byte_t> buffer1, buffer2;
            // if we want to exclude this SNP, we will not perform decompression
            genfile::bgen::read_and_parse_genotype_data_block<ProbSetter>(
                m_bgen_file, context, setter, &buffer1, &buffer2, exclude_snp);
            if (!exclude_snp) {
                // if not exclude SNP, then we need to calculate the INFO score
                // and also the MAF
                m_existed_snps_index[RSID] = snp_res.size();
                snp_res.emplace_back(SNP(RSID, chr_code, SNP_position,
                                         alleles.front(), alleles.back(),
                                         prefix, byte_pos));
                m_unfiltered_marker_ct++;
                // directly ignore all others?
            }
        }
    }
    fprintf(stderr, "\n");
    snp_res.shrink_to_fit(); // so that it will be more suitable
    if (m_bgen_file.is_open()) m_bgen_file.close();

    if (dup_list.size() != 0) {
        std::ofstream log_file_stream;
        std::string dup_name = out_prefix + ".valid";
        log_file_stream.open(dup_name.c_str());
        if (!log_file_stream.is_open()) {
            std::string error_message = "ERROR: Cannot open file: " + dup_name;
            throw std::runtime_error(error_message);
        }
        for (auto&& snp : snp_res) {
            if (dup_list.find(snp.rs()) != dup_list.end()) continue;
            log_file_stream << snp.rs() << "\t" << snp.chr() << "\t"
                            << snp.loc() << "\t" << snp.ref() << "\t"
                            << snp.alt() << std::endl;
        }
        log_file_stream.close();
        std::string error_message =
            "ERROR: Duplicated SNP ID detected!. Valid SNP ID stored at "
            + dup_name + ". You can avoid this error by using --extract "
            + dup_name;
        throw std::runtime_error(error_message);
    }

    return snp_res;
}


void BinaryGen::dosage_score(std::vector<Sample_lite>& current_prs_score,
                             size_t start_index, size_t end_bound,
                             const size_t region_index)
{
    m_cur_file = "";
    std::vector<genfile::byte_t> buffer1, buffer2;
    size_t num_included_samples = current_prs_score.size();
    for (size_t i_snp = start_index; i_snp < end_bound; ++i_snp) {
        if (!m_existed_snps[i_snp].in(region_index)) continue;

        auto&& snp = m_existed_snps[i_snp];
        if (m_cur_file.empty() || snp.file_name().compare(m_cur_file) != 0) {
            if (m_bgen_file.is_open()) m_bgen_file.close();
            std::string bgen_name = snp.file_name() + ".bgen";
            m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
            if (!m_bgen_file.is_open()) {
                std::string error_message =
                    "ERROR: Cannot open bgen file: " + snp.file_name();
                throw std::runtime_error(error_message);
            }
            m_cur_file = snp.file_name();
        }
        m_bgen_file.seekg(snp.byte_pos(), std::ios_base::beg);

        Data probability;
        ProbSetter setter(&probability);
        genfile::bgen::read_and_parse_genotype_data_block<ProbSetter>(
            m_bgen_file, m_bgen_info[snp.file_name()], setter, &buffer1,
            &buffer2, false);
        std::vector<size_t> missing_samples;
        double total = 0.0;
        std::vector<double> score(num_included_samples);
        size_t cur_sample = 0;
        for (size_t i_sample = 0; i_sample < probability.size(); ++i_sample) {
            auto&& prob = probability[i_sample];
            if (prob.size() != 3) {
                // this is likely phased
                std::string message = "ERROR: Currently don't support phased "
                                      "data (It is because the lack of "
                                      "development time)\n";
                throw std::runtime_error(message);
            }
            double expected = 0.0;
            if (IS_SET(m_sample_include.data(),
                       i_sample)) // to ignore unwanted samples
            {
                // we want g to be signed so that when -2, it will not cause us
                // troubles
                for (int g = 0; g < (int) prob.size(); ++g) {
                    if (*max_element(prob.begin(), prob.end())
                        < filter.hard_threshold)
                    {
                        missing_samples.push_back(i_sample);
                        break;
                    }
                    else
                    {
                        int geno = (!snp.is_flipped()) ? 2 - g : g;
                        if (m_model == +MODEL::HETEROZYGOUS && geno == 2)
                            geno = 0;
                        else if (m_model == +MODEL::DOMINANT && geno == 2)
                            geno = 1;
                        else if (m_model == +MODEL::RECESSIVE)
                            geno = std::max(geno - 1, 0);
                        expected += prob[g] * geno;
                    }
                }
                score[cur_sample++] = expected;
                total += expected;
            }
        }
        // now process the missing and clean stuff
        // we divide the mean by 2 so that in situation where there is no
        // dosage, it will behave the same as the genotype data.

        size_t num_miss = missing_samples.size();
        if (num_included_samples - num_miss == 0) {
            m_existed_snps[i_snp].invalidate();
            continue;
        }
        double mean =
            total / (((double) num_included_samples - (double) num_miss) * 2);
        size_t i_missing = 0;
        double stat = snp.stat();
        for (size_t i_sample = 0; i_sample < num_included_samples; ++i_sample) {
            if (i_missing < num_miss && i_sample == missing_samples[i_missing])
            {
                if (m_scoring == SCORING::MEAN_IMPUTE)
                    current_prs_score[i_sample].prs += stat * mean;
                if (m_scoring != SCORING::SET_ZERO)
                    current_prs_score[i_sample].num_snp++;

                i_missing++;
            }
            else
            { // not missing sample
                if (m_scoring == SCORING::CENTER) {
                    // if centering, we want to keep missing at 0
                    current_prs_score[i_sample].prs -= stat * mean;
                }
                // again, so that it will generate the same result as
                // genotype file format when we are 100% certain of the
                // genotypes
                current_prs_score[i_sample].prs += score[i_sample] * stat * 0.5;
                current_prs_score[i_sample].num_snp++;
            }
        }
    }
}

void BinaryGen::hard_code_score(std::vector<Sample_lite>& current_prs_score,
                                size_t start_index, size_t end_bound,
                                const size_t region_index)
{

    m_cur_file = "";
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uintptr_t ulii = 0;
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    uintptr_t final_mask = get_final_mask(m_sample_ct);

    size_t num_included_samples = current_prs_score.size();
    // index is w.r.t. partition, which contain all the information
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    //    uintptr_t* genotype = new uintptr_t[unfiltered_sample_ctl * 2];
    for (size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
    { // for each SNP
        if (!m_existed_snps[i_snp].in(region_index)) continue;
        if (load_and_collapse_incl(m_existed_snps[i_snp].byte_pos(),
                                   m_existed_snps[i_snp].file_name(),
                                   m_unfiltered_sample_ct, m_sample_ct,
                                   m_sample_include.data(), final_mask, false,
                                   m_tmp_genotype.data(), genotype.data()))
        {
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        uintptr_t* lbptr = genotype.data();
        uii = 0;
        std::vector<size_t> missing_samples;
        double stat = m_existed_snps[i_snp].stat();
        bool flipped = m_existed_snps[i_snp].is_flipped();
        std::vector<double> genotypes(num_included_samples);

        uint32_t sample_idx = 0;
        int nmiss = 0;
        int aa = 0, aA = 0, AA = 0;
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
                        int g = (ukk == 3) ? 2 : ukk;
                        switch (g)
                        {
                        case 0: aa++; break;
                        case 1: aA++; break;
                        case 2: AA++; break;
                        }
                        genotypes[sample_idx] = g;
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


        size_t i_missing = 0;

        if (num_included_samples - nmiss == 0) {
            m_existed_snps[i_snp].invalidate();
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
        if (m_model == +MODEL::HETEROZYGOUS) {
            // 010
            aa += AA;
            AA = 0;
        }
        else if (m_model == +MODEL::DOMINANT)
        {
            // 011;
            aA += AA;
            AA = 0;
        }
        else if (m_model == +MODEL::RECESSIVE)
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
        for (size_t i_sample = 0; i_sample < num_included_samples; ++i_sample) {
            if (i_missing < num_miss && i_sample == missing_samples[i_missing])
            {

                if (m_scoring == SCORING::MEAN_IMPUTE)
                    current_prs_score[i_sample].prs += center_score;
                if (m_scoring != SCORING::SET_ZERO)
                    current_prs_score[i_sample].num_snp++;
                i_missing++;
            }
            else
            { // not missing sample
                if (m_scoring == SCORING::CENTER) {
                    // if centering, we want to keep missing at 0
                    current_prs_score[i_sample].prs -= center_score;
                }

                int g = (flipped) ? fabs(genotypes[i_sample] - 2)
                                  : genotypes[i_sample];
                if (m_model == +MODEL::HETEROZYGOUS) {
                    g = (g == 2) ? 0 : g;
                }
                else if (m_model == +MODEL::RECESSIVE)
                {
                    g = std::max(0, g - 1);
                }
                else if (m_model == +MODEL::DOMINANT)
                {
                    g = (g == 2) ? 1 : g;
                }

                current_prs_score[i_sample].prs += g * stat * 0.5;
                current_prs_score[i_sample].num_snp++;
            }
        }
    }
}

void BinaryGen::read_score(std::vector<Sample_lite>& current_prs_score,
                           size_t start_index, size_t end_bound,
                           const size_t region_index)
{
    if (filter.use_hard) {
        hard_code_score(current_prs_score, start_index, end_bound,
                        region_index);
        return;
    }
    else
        dosage_score(current_prs_score, start_index, end_bound, region_index);
}


Sample BinaryGen::get_sample(std::vector<std::string>& token, bool ignore_fid,
                             bool has_sex, int sex_col,
                             std::vector<int>& sex_info)
{
    assert(ignore_fid && token.size() > 1);
    std::string id = (ignore_fid) ? token[0] : token[0] + "_" + token[1];
    // this will pose problem when there are duplicated IID names even if they
    // are from different family. However, we don't know how bgen store the
    // sample information (do they contain the FID?) so we will have to work
    // this way
    if (m_sample_index_check.find((ignore_fid) ? token[0] : token[1])
        == m_sample_index_check.end())
    {
        Sample cur_sample;
        cur_sample.FID = token[0];
        cur_sample.IID = (ignore_fid) ? token[0] : token[1];
        cur_sample.pheno = "NA";
        cur_sample.included = false;
        cur_sample.has_pheno = false;
        if (m_remove_sample) {
            cur_sample.included = (m_sample_selection_list.find(id)
                                   == m_sample_selection_list.end());
        }
        else
        {
            cur_sample.included = (m_sample_selection_list.find(id)
                                   != m_sample_selection_list.end());
        }
        // there isn't any founder/nonfounder, so directly using this is ok
        m_sample_ct += cur_sample.included;
        cur_sample.num_snp = 0;
        if (has_sex) {
            try
            {
                sex_info.push_back(misc::convert<int>(token[sex_col]));
            }
            catch (const std::runtime_error& er)
            {
                throw std::runtime_error("ERROR: Invalid sex coding!\n");
            }
        }
        return cur_sample;
    }
    else
    {
        std::string error_message = "ERROR: Duplicated sample: " + id;
        throw std::runtime_error(error_message);
    }
}

std::vector<Sample> BinaryGen::preload_samples(std::string pheno,
                                               Reporter& reporter,
                                               bool has_header, bool ignore_fid)
{
    std::vector<Sample> sample_res;
    std::ifstream pheno_file;
    pheno_file.open(pheno.c_str());
    if (!pheno_file.is_open()) {
        std::string error_message =
            "ERROR: Cannot open phenotype file: " + pheno;
        throw std::runtime_error(error_message);
    }
    std::string line;
    // we will assume we have all
    std::string first_line;
    std::getline(pheno_file, first_line);
    misc::trim(first_line);
    std::vector<std::string> possible_header = misc::split(first_line);
    std::string second_line;
    std::getline(pheno_file, second_line);
    misc::trim(second_line);
    std::vector<std::string> token = misc::split(second_line);
    bool bgen_sample = false;
    if (token.size() > 3) // FID IID and Missing are required, then phenotype
    {
        if (token[0].compare("0") == 0 && token[1].compare("0") == 0
            && token[2].compare("0") == 0)
        {
            bgen_sample = true;
            for (size_t i = 3; i < token.size(); ++i) {
                if (token[i].compare("D") != 0 && token[i].compare("C") != 0
                    && token[i].compare("P") != 0 && token[i].compare("B") != 0)
                {
                    bgen_sample = false;
                    break;
                }
            }
        }
    }
    bool has_sex = false;
    int sex_col = 0;
    std::vector<int> sex_info;
    if (bgen_sample) // then we know the first line is header
    {
        reporter.report("Detected bgen sample file format");
        for (size_t i = 3; i < possible_header.size(); ++i) {
            if (possible_header[i].compare("Sex") == 0
                || possible_header[i].compare("sex") == 0
                || possible_header[i].compare("SEX") == 0)
            {
                has_sex = true;
                sex_col = i;
                break;
            }
        }
        if (has_sex) {
            if (token[sex_col].compare("D") != 0) {
                std::string error_message =
                    "ERROR: Sex must be coded as \"D\" in bgen sample file!";
                throw std::runtime_error(error_message);
            }
        }
    }
    else
    {
        // this is just a normal line
        if (!has_header) {
            // need to check that for normal pheno file, whether there is
            // a header
            if (possible_header.front().compare("FID") == 0
                || (possible_header.size() > 2
                    && possible_header[1].compare("IID") == 0)
                || possible_header.front().compare("IID") == 0)
            {
                // these are simple test. might not be correct though
                // if fit, this is header
            }
            else
            {
                // this isn't a header
                if (possible_header.size() < 2 && ignore_fid) {
                    std::string error_message =
                        "Malformed phenotype file: line:1 only has 1 column.\n";
                    throw std::runtime_error(error_message);
                }
                sample_res.push_back(get_sample(possible_header, ignore_fid,
                                                has_sex, sex_col, sex_info));
                m_sample_index_check[sample_res.back().IID] =
                    sample_res.size() - 1;
            }
        }
        if (token.size() < 2 && ignore_fid) {
            std::string error_message =
                "Malformed phenotype file: line:2 only has 1 column.\n";
            throw std::runtime_error(error_message);
        }
        sample_res.push_back(
            get_sample(token, ignore_fid, has_sex, sex_col, sex_info));
        m_sample_index_check[sample_res.back().IID] = sample_res.size() - 1;
    }
    size_t line_id = 2;
    while (std::getline(pheno_file, line)) {
        misc::trim(line);
        line_id++;
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (token.size() < ((has_sex) ? (sex_col) : (1 + !ignore_fid))) {
            std::string error_message =
                "ERROR: Line " + std::to_string(line_id)
                + " must have at least "
                + std::to_string((has_sex) ? (sex_col) : (1 + !ignore_fid))
                + " columns! Number of column=" + std::to_string(token.size());
            throw std::runtime_error(error_message);
        }
        sample_res.push_back(
            get_sample(token, ignore_fid, has_sex, sex_col, sex_info));
        m_sample_index_check[sample_res.back().IID] = sample_res.size() - 1;
    }
    m_unfiltered_sample_ct = sample_res.size();
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    // uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;

    // m_sex_male = new uintptr_t[m_unfiltered_sample_ctl];
    // std::memset(m_sex_male, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

    // assume all are founder
    m_founder_info.resize(unfiltered_sample_ctl, 0);
    m_founder_ct = 0;

    m_sample_include.resize(unfiltered_sample_ctl, 0);
    for (size_t i = 0; i < sample_res.size(); ++i) {
        if (sample_res[i].included) {
            m_founder_ct++;
            SET_BIT(i, m_founder_info.data());
            SET_BIT(i, m_sample_include.data());
        }
    }

    m_num_male = 0, m_num_female = 0,
    m_num_ambig_sex = (has_sex) ? 0 : m_unfiltered_sample_ct;
    for (size_t i = 0; i < sex_info.size() && has_sex; ++i) {
        if (!sample_res[i].included) continue;
        switch (sex_info[i])
        {
        case 1:
            m_num_male++;
            // SET_BIT(i, m_sex_male);
            break;
        case 2: m_num_female++; break;
        default: m_num_ambig_sex++;
        }
    }
    pheno_file.close();
    return sample_res;
}

std::vector<Sample> BinaryGen::load_samples(bool ignore_fid)
{
    std::unordered_set<std::string> dup_check;
    bool first = true;
    for (auto&& prefix : m_genotype_files) {
        if (m_bgen_file.is_open()) m_bgen_file.close();
        std::string bgen_name = prefix + ".bgen";
        m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
        if (!m_bgen_file.is_open()) {
            std::string error_message =
                "ERROR: Cannot open bgen file: " + bgen_name;
            throw std::runtime_error(error_message);
        }
        uint32_t offset = 0;
        bgenlib::read_little_endian_integer(m_bgen_file, &offset);
        uint32_t header_size = 0, number_of_snp_blocks = 0,
                 number_of_samples = 0, flags = 0;
        char magic[4];
        std::size_t fixed_data_size = 20;
        std::vector<char> free_data;
        bgenlib::read_little_endian_integer(m_bgen_file, &header_size);
        assert(header_size >= fixed_data_size);
        bgenlib::read_little_endian_integer(m_bgen_file, &number_of_snp_blocks);
        bgenlib::read_little_endian_integer(m_bgen_file, &number_of_samples);
        m_bgen_file.read(&magic[0], 4);
        free_data.resize(header_size - fixed_data_size);
        m_bgen_file.read(&free_data[0], free_data.size());
        bgenlib::read_little_endian_integer(m_bgen_file, &flags);
        if ((magic[0] != 'b' || magic[1] != 'g' || magic[2] != 'e'
             || magic[3] != 'n')
            && (magic[0] != 0 || magic[1] != 0 || magic[2] != 0
                || magic[3] != 0))
        {
            throw std::runtime_error("ERROR: Incorrect magic string!\nPlease "
                                     "check you have provided a valid bgen "
                                     "file!");
        }
        if (m_bgen_file) {
            bgenlib::Context current_context;
            current_context.number_of_samples = number_of_samples;
            current_context.number_of_variants = number_of_snp_blocks;
            current_context.magic.assign(&magic[0], &magic[0] + 4);
            // current_context.free_data.assign( free_data.begin(),
            // free_data.end() ) ;
            current_context.flags = flags;
            m_bgen_info[prefix] = current_context;
            m_offset_map[prefix] = offset;
        }
        else
        {
            throw std::runtime_error("ERROR: Problem reading bgen file!");
        }
        if (m_bgen_info[prefix].number_of_samples != m_sample_names.size()) {
            std::string error_message =
                "ERROR: Number of sample in bgen does not match those in "
                "phenotype file! ("
                + std::to_string(m_bgen_info[prefix].number_of_samples) + " vs "
                + std::to_string(m_sample_names.size()) + ")";
            throw std::runtime_error(error_message);
        }

        uint32_t const compressionType =
            m_bgen_info[prefix].flags & bgenlib::e_CompressedSNPBlocks;
        if (compressionType == bgenlib::e_ZstdCompression) {
            throw std::runtime_error(
                "ERROR: zstd compression currently not supported");
        }
        if ((m_bgen_info[prefix].flags & bgenlib::e_SampleIdentifiers) && first)
        { // only read in the sample information for the first bgen
            first = false;
            uint32_t sample_block_size = 0;
            uint32_t actual_number_of_samples = 0;
            uint16_t identifier_size;
            std::string identifier;
            std::size_t bytes_read = 0;
            // the bgen format actually double stored the number of samples and
            // the block size is the sample_block size
            bgenlib::read_little_endian_integer(m_bgen_file,
                                                &sample_block_size);
            bgenlib::read_little_endian_integer(m_bgen_file,
                                                &actual_number_of_samples);
            bytes_read += 8;
            assert(actual_number_of_samples
                   == m_bgen_info[prefix].number_of_samples);
            for (size_t i = 0; i < actual_number_of_samples; ++i) {
                bgenlib::read_length_followed_by_data(
                    m_bgen_file, &identifier_size, &identifier);
                if (m_bgen_file) {
                    bytes_read += sizeof(identifier_size) + identifier_size;
                    if (m_sample_index_check.find(identifier)
                        == m_sample_index_check.end())
                    {
                        throw std::runtime_error("ERROR: Sample mismatch "
                                                 "between bgen and phenotype "
                                                 "file!");
                    }
                    else if (m_sample_index_check[identifier] != i)
                    {
                        throw std::runtime_error("ERROR: Sample sequence "
                                                 "differ between bgen and "
                                                 "phenotype file!");
                    }
                    if (dup_check.find(identifier) == dup_check.end()) {
                        dup_check.insert(identifier);
                        // only for checking
                        if (m_sample_names[i].IID.compare(identifier) != 0) {
                            throw std::runtime_error("ERROR: Sample name "
                                                     "mismatch between bgen "
                                                     "and phenotype file!");
                        }
                    }
                    else
                    {
                        throw std::runtime_error(
                            "ERROR: Duplicated sample within bgen file!");
                    }
                }
                else
                {
                    throw std::runtime_error(
                        "ERROR: Problem reading bgen file!");
                }
            }
            assert(bytes_read == sample_block_size);
        }
    }
    return std::vector<Sample>();
}
/*
void BinaryGen::get_header()
{

}
*/
