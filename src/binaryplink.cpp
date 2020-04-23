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

BinaryPlink::BinaryPlink(const GenoFile& geno, const Phenotype& pheno,
                         const std::string& delim, Reporter* reporter)
{
    const std::string message = initialize(geno, pheno, delim, "bed", reporter);
    if (m_sample_file.empty())
    { m_sample_file = m_genotype_file_names.front() + ".fam"; }
    m_reporter->report(message);
}

std::unordered_set<std::string>
BinaryPlink::get_founder_info(std::unique_ptr<std::istream>& famfile)
{
    std::string line;
    std::vector<std::string> token;
    std::unordered_set<std::string> founder_info;
    while (std::getline(*famfile, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        misc::split(token, line);
        if (token.size() != 6)
        {
            throw std::runtime_error(
                "Error: Malformed fam file. Fam file should have 6 columns."
                "Line: "
                + std::to_string(m_unfiltered_sample_ct + 1) + "\n");
        }
        founder_info.insert(token[+FAM::FID] + m_delim + token[+FAM::IID]);
        ++m_unfiltered_sample_ct;
    }
    (*famfile).clear();
    (*famfile).seekg(0);
    return founder_info;
}
std::vector<Sample_ID> BinaryPlink::gen_sample_vector()
{
    assert(m_genotype_file_names.size() > 0);
    auto famfile = misc::load_stream(m_sample_file);
    m_unfiltered_sample_ct = 0;
    // will also count number of samples here. Which initialize the important
    // m_unfiltered_sample_ct
    std::unordered_set<std::string> founder_info = get_founder_info(famfile);
    init_sample_vectors();
    // we will return the sample_name
    std::vector<Sample_ID> sample_name;
    std::unordered_set<std::string> samples_in_fam;
    std::vector<std::string> duplicated_sample_id;
    // for purpose of output
    uintptr_t sample_index = 0; // this is just for error message
    std::vector<std::string> token;
    std::string line;
    while (std::getline(*famfile, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        misc::split(token, line);
        // we have already checked for malformed file
        gen_sample(+FAM::FID, +FAM::IID, +FAM::SEX, +FAM::FATHER, +FAM::MOTHER,
                   sample_index, founder_info, token[+FAM::PHENOTYPE], token,
                   sample_name, samples_in_fam, duplicated_sample_id);
        ++sample_index;
    }

    if (duplicated_sample_id.size() > 0)
    {
        throw std::runtime_error(
            "Error: A total of " + misc::to_string(duplicated_sample_id.size())
            + " duplicated samples detected!\n"
            + "Please ensure all samples have an unique identifier");
    }
    famfile.reset();
    post_sample_read_init();
    return sample_name;
}

bool BinaryPlink::calc_freq_gen_inter(const QCFiltering& filter_info,
                                      const std::string&, Genotype* genotype)
{
    const double sample_ct_recip = 1.0 / (static_cast<double>(m_sample_ct));
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    const size_t total_snp = genotype->m_existed_snps.size();
    std::vector<bool> retain_snps(total_snp, false);
    double progress = 0.0, prev_progress = -1.0;
    double cur_maf, cur_geno;
    std::streampos byte_pos;
    size_t processed_count = 0;
    size_t retained = 0;
    size_t cur_file_idx = 0;
    uint32_t ref_count = 0;
    uint32_t het_count = 0;
    uint32_t alt_count = 0;
    uint32_t ref_founder_count = 0;
    uint32_t het_founder_count = 0;
    uint32_t alt_founder_count = 0;
    uint32_t total_alleles = 0;
    uint32_t missing = 0;
    uint32_t total_founder_alleles = 0;
    // initialize the sample inclusion mask
    for (auto&& snp : genotype->m_existed_snps)
    {
        progress = static_cast<double>(processed_count)
                   / static_cast<double>(total_snp) * 100;
        if (!m_reporter->unit_testing() && !m_reporter->unit_testing()
            && progress - prev_progress > 0.01)
        {
            fprintf(stderr, "\rCalculating allele frequencies: %03.2f%%",
                    progress);
            prev_progress = progress;
        }
        ++processed_count;

        snp.get_file_info(cur_file_idx, byte_pos, m_is_ref);

        m_genotype_file.read(m_genotype_file_names[cur_file_idx] + ".bed",
                             byte_pos,
                             static_cast<long long>(unfiltered_sample_ct4),
                             reinterpret_cast<char*>(m_tmp_genotype.data()));
        // calculate the MAF using PLINK2 function (take into account of founder
        // status)
        single_marker_freqs_and_hwe(
            unfiltered_sample_ctv2, m_tmp_genotype.data(),
            m_sample_include2.data(), m_founder_include2.data(), m_sample_ct,
            &ref_count, &het_count, &alt_count, m_founder_ct,
            &ref_founder_count, &het_founder_count, &alt_founder_count);
        total_alleles = ref_count + het_count + alt_count;
        cur_geno =
            1.0 - (static_cast<int32_t>(total_alleles)) * sample_ct_recip;
        // filter by genotype missingness
        if (filter_info.geno < cur_geno)
        {
            ++m_num_geno_filter;
            continue;
        }
        // filter by maf
        total_founder_alleles =
            (ref_founder_count + het_founder_count + alt_founder_count);
        assert(m_founder_ct >= total_founder_alleles);
        missing = static_cast<uint32_t>(m_founder_ct) - total_founder_alleles;
        if (missing == m_founder_ct)
        {
            // invalid SNP as all missing. Should just remove it
            ++m_num_miss_filter;
            continue;
        }
        cur_maf =
            (static_cast<double>(2 * alt_founder_count + het_founder_count))
            / (static_cast<double>(2 * total_founder_alleles));
        cur_maf = (cur_maf > 0.5) ? 1 - cur_maf : cur_maf;
        if (misc::logically_equal(cur_maf, 0.0)
            || misc::logically_equal(cur_maf, 1.0))
        {
            // none of the sample contain this SNP
            // still count as MAF filtering (for now)
            ++m_num_maf_filter;
            continue;
        }
        if (cur_maf < filter_info.maf)
        {
            ++m_num_maf_filter;
            continue;
        }
        // if we can reach here, it is not removed
        snp.set_counts(ref_founder_count, het_founder_count, alt_founder_count,
                       missing, m_is_ref);
        ++retained;
        // we need to -1 because we put processed_count ++ forward
        // to avoid continue skipping out the addition
        retain_snps[processed_count - 1] = true;
    }
    if (!m_reporter->unit_testing())
        fprintf(stderr, "\rCalculating allele frequencies: %03.2f%%\n", 100.0);
    // now update the vector
    if (retained != total_snp) { genotype->shrink_snp_vector(retain_snps); }
    return true;
}

size_t BinaryPlink::transverse_bed_for_snp(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const std::string mismatch_snp_record_name, const size_t idx,
    const uintptr_t unfiltered_sample_ct4, const uintptr_t bed_offset,
    std::unique_ptr<std::istream> bim,
    std::unordered_set<std::string>& duplicated_snps,
    std::unordered_set<std::string>& processed_snps,
    std::vector<bool>& retain_snp, bool& chr_error, bool& sex_error,
    Genotype* genotype)
{
    assert(bim->is_open());
    assert(genotype != nullptr);
    const std::string mismatch_source = m_is_ref ? "Reference" : "Base";
    std::vector<std::string> bim_token(6, "");
    std::string line;
    std::string prev_chr = "";
    std::streampos byte_pos;
    size_t num_snp_read = 0;
    size_t num_retained = 0;
    size_t chr_num = 0;
    while (std::getline(*bim, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        ++num_snp_read;
        misc::split(bim_token, line);
        if (bim_token.size() < 6)
        {
            throw std::runtime_error(
                "Error: Malformed bim file. Less than 6 column on "
                "line: "
                + misc::to_string(num_snp_read) + "\n");
        }
        size_t loc = ~size_t(0);
        try
        {
            loc = misc::convert<size_t>(bim_token[+BIM::BP]);
        }
        catch (...)
        {
            throw std::runtime_error(
                "Error: Invalid SNP coordinate: " + bim_token[+BIM::RS] + ":"
                + bim_token[+BIM::BP]
                + "\nPlease check you have the correct input");
        }
        byte_pos = static_cast<std::streampos>(
            bed_offset + ((num_snp_read - 1) * (unfiltered_sample_ct4)));

        if (!check_chr(bim_token[+BIM::CHR], prev_chr, chr_num, chr_error,
                       sex_error))
        { continue; }
        SNP cur_snp(bim_token[+BIM::RS], chr_num, loc, bim_token[+BIM::A1],
                    bim_token[+BIM::A2], idx, byte_pos);
        if (process_snp(exclusion_regions, mismatch_snp_record_name,
                        mismatch_source, "", cur_snp, processed_snps,
                        duplicated_snps, retain_snp, genotype))
        { ++num_retained; }
    }
    bim.reset();
    return num_retained;
}
void BinaryPlink::gen_snp_vector(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const std::string& out_prefix, Genotype* target)
{
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    const std::string mismatch_snp_record_name = out_prefix + ".mismatch";
    const std::string mismatch_print_type = (m_is_ref) ? "Reference" : "Base";
    std::unordered_set<std::string> processed_snps;
    std::unordered_set<std::string> duplicated_snp;
    std::vector<std::string> bim_token;
    auto&& genotype = (m_is_ref) ? target : this;
    std::vector<bool> retain_snp(genotype->m_existed_snps.size(), false);
    std::string bed_name, line;
    std::string prev_chr = "", snpid = "";
    std::string prefix;
    uintptr_t bed_offset;
    size_t num_retained = 0;
    size_t num_snp_read = 0;
    std::streampos byte_pos;
    bool chr_error = false, sex_error = false;
    for (size_t idx = 0; idx < m_genotype_file_names.size(); ++idx)
    {
        // go through each genotype file
        prefix = m_genotype_file_names[idx];
        bed_name = prefix + ".bed";
        // make sure we reset the flag of the ifstream by closing it before use
        auto bim = misc::load_stream(prefix + ".bim");
        // First pass, get the number of marker in bed & bim
        // as we want the number for checking, num_snp_read will start at 0
        num_snp_read = misc::get_num_line(bim);
        check_bed(bed_name, num_snp_read, bed_offset);
        num_retained += transverse_bed_for_snp(
            exclusion_regions, mismatch_snp_record_name, idx,
            unfiltered_sample_ct4, bed_offset, std::move(bim), duplicated_snp,
            processed_snps, retain_snp, chr_error, sex_error, genotype);
    }
    // try to release memory
    if (num_retained != genotype->m_existed_snps.size())
    {
        genotype->shrink_snp_vector(retain_snp);
        // need to update index search after we updated the vector
        genotype->update_snp_index();
    }
    if (duplicated_snp.size() != 0)
    {
        throw std::runtime_error(
            genotype->print_duplicated_snps(duplicated_snp, out_prefix));
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
    { throw std::runtime_error("Cannot read bed file: " + bed_name); }
    bed.seekg(0, bed.end);
    llxx = bed.tellg();
    if (!llxx)
    { throw std::runtime_error("Error: Empty .bed file: " + bed_name); }
    bed.clear();
    bed.seekg(0, bed.beg);
    char version_check[3];
    bed.read(version_check, 3);
    uii = static_cast<uint32_t>(bed.gcount());
    llyy = static_cast<int64_t>(unfiltered_sample_ct4 * num_marker);
    llzz =
        static_cast<int64_t>(m_unfiltered_sample_ct * ((num_marker + 3) / 4));
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
            throw std::runtime_error(
                "Error: Invalid header bytes in .bed file: " + bed_name);
        }
        llyy = llzz;
        bed_offset = 2;
    }
    if (llxx != llyy)
    {
        if ((*version_check == '#')
            || ((uii == 3) && (!memcmp(version_check, "chr", 3))))
        {
            throw std::runtime_error("Error: Invalid header bytes in PLINK "
                                     "1 .bed file: "
                                     + bed_name
                                     + "  (Is this a UCSC "
                                       "Genome\nBrowser BED file instead?)");
        }
        else
        {
            throw std::runtime_error("Error: Invalid .bed file size for "
                                     + bed_name);
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
    std::vector<PRS>& prs_list,
    const std::vector<size_t>::const_iterator& start_idx,
    const std::vector<size_t>::const_iterator& end_idx, bool reset_zero,
    bool read_only)
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
    // currently hard code ploidy to 2. Will keep it this way unil we know how
    // to properly handly non-diploid chromosomes in human and other organisms
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
        (m_prs_calculation.missing_score != MISSING_SCORE::SET_ZERO) * ploidy;
    // this indicate if we want the mean of the genotype to be 0 (missingness =
    // 0)
    const bool is_centre =
        (m_prs_calculation.missing_score == MISSING_SCORE::CENTER);
    // this indicate if we want to impute the missing genotypes using the
    // population mean
    const bool mean_impute =
        (m_prs_calculation.missing_score == MISSING_SCORE::MEAN_IMPUTE);
    // check if it is not the frist run, if it is the first run, we will reset
    // the PRS to zero instead of addint it up
    bool not_first = !reset_zero;
    double stat, maf, adj_score, miss_score;
    // m_cur_file = ""; // just close it
    // if (m_bed_file.is_open()) { m_bed_file.close(); }
    // initialize the genotype vector to store the binary genotypes
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    std::vector<size_t>::const_iterator cur_idx = start_idx;
    std::streampos cur_line;
    std::string file_name;
    size_t file_idx;
    for (; cur_idx != end_idx; ++cur_idx)
    {
        auto&& cur_snp = m_existed_snps[(*cur_idx)];
        if (!cur_snp.stored_genotype())
        {
            cur_snp.get_file_info(file_idx, cur_line, false);
            file_name = m_genotype_file_names[file_idx] + ".bed";
            // we now read the genotype from the file by calling
            // load_and_collapse_incl
            // important point to note here is the use of m_sample_include and
            // m_sample_ct instead of using the m_founder m_founder_info as the
            // founder vector is for LD calculation whereas the sample_include
            // is for PRS
            m_genotype_file.read(
                file_name, cur_line, unfiltered_sample_ct4,
                reinterpret_cast<char*>(m_tmp_genotype.data()));
            if (!cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct,
                                    m_prs_calculation.use_ref_maf))
            {
                // we need to calculate the MA
                // if we want to use reference, we will always have calculated
                // the MAF
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
                cur_snp.set_counts(homcom_ct, het_ct, homrar_ct, missing_ct,
                                   false);
            }
            if (m_unfiltered_sample_ct != m_sample_ct)
            {
                copy_quaterarr_nonempty_subset(
                    m_tmp_genotype.data(), m_calculate_prs.data(),
                    static_cast<uint32_t>(m_unfiltered_sample_ct),
                    static_cast<uint32_t>(m_sample_ct), genotype.data());
            }
            else
            {
                genotype = m_tmp_genotype;
                genotype[(m_unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
            }
            if (read_only) { cur_snp.assign_genotype(genotype); }
        }
        else
        {
            genotype = cur_snp.get_genotype();
            cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct,
                               m_prs_calculation.use_ref_maf);
        }
        if (m_founder_ct == missing_ct)
        {
            // problematic snp
            cur_snp.invalid();
            continue;
        }
        homcom_weight = m_homcom_weight;
        het_weight = m_het_weight;
        homrar_weight = m_homrar_weight;
        maf = 1.0
              - static_cast<double>(homcom_weight * homcom_ct
                                    + het_ct * het_weight
                                    + homrar_weight * homrar_ct)
                    / (static_cast<double>((homcom_ct + het_ct + homrar_ct)
                                           * ploidy));
        if (cur_snp.is_flipped())
        {
            std::swap(homcom_weight, homrar_weight);
            maf = 1.0 - maf;
        }
        stat = cur_snp.stat();
        adj_score = 0;
        if (is_centre) { adj_score = ploidy * stat * maf; }
        miss_score = 0;
        if (mean_impute) { miss_score = ploidy * stat * maf; }
        // now we go through the SNP vector
        if (!read_only)
        {
            read_prs(genotype, prs_list, ploidy, stat, adj_score, miss_score,
                     miss_count, homcom_weight, het_weight, homrar_weight,
                     not_first);
        }
        // indicate that we've already read in the first SNP and no longer need
        // to reset the PRS
        not_first = true;
    }
}
