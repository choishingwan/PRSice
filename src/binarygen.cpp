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
    uint32_t offset;
    genfile::bgen::read_offset(bgen_file, &offset);
    genfile::bgen::read_header_block(bgen_file, &context);
    context.offset = offset;
    return context;
}

void BinaryGen::check_sample_consistent(const genfile::bgen::Context& context,
                                        std::istream& bgen_file)
{
    // only do this if our gen file contains the sample information
    if (context.flags & genfile::bgen::e_SampleIdentifiers)
    {
        uint32_t offset;
        genfile::bgen::read_offset(bgen_file, &offset);
        genfile::bgen::Context tmp_context;
        genfile::bgen::read_header_block(bgen_file, &tmp_context);
        size_t sample_idx = 0;
        genfile::bgen::read_sample_identifier_block(
            bgen_file, tmp_context, [this, &sample_idx](const std::string& id) {
                auto&& find_id = m_sample_selection_list.find(id)
                                 != m_sample_selection_list.end();
                bool inclusion = m_remove_sample ^ find_id;
                if (inclusion)
                {
                    std::string ref_id =
                        (m_ignore_fid ? ""
                                      : m_sample_id[sample_idx].FID + m_delim)
                        + m_sample_id[sample_idx].IID;
                    if (ref_id != id && m_sample_id[sample_idx].IID != id)
                    {
                        throw std::runtime_error(
                            "Error: Sample mismatch "
                            "between bgen and phenotype file! Name in BGEN "
                            "file is "
                            ":"
                            + id + " and in phentoype file is: " + ref_id
                            + ". Please note that PRSice require the bgen file "
                              "and "
                              "the .sample (or phenotype file if sample file "
                              "is "
                              "not provided) to have sample in the same order. "
                              "(We "
                              "might be able to losen this requirement in "
                              "future "
                              "when we have more time)");
                    }
                    ++sample_idx;
                }
            });
    }
}

size_t BinaryGen::transverse_bgen_for_snp(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const std::string mismatch_snp_record_name, const size_t file_idx,
    std::unique_ptr<std::istream> bgen_file,
    std::unordered_set<std::string>& duplicated_snps,
    std::unordered_set<std::string>& processed_snps,
    std::vector<bool>& retain_snp, bool& chr_error, bool& sex_error,
    Genotype* genotype)
{
    // skip the offset (first 4 are used to store the offset)
    // offset contains the byte location of the first variant data block
    assert(m_context_map.size() > file_idx);
    auto&& context = m_context_map[file_idx];
    bgen_file->seekg(context.offset + 4);
    const size_t num_snp = context.number_of_variants;
    const std::string mismatch_source = m_is_ref ? "Reference" : "Base";
    const std::string name = m_genotype_file_names[file_idx] + ".bgen";
    std::string SNPID;
    std::string RSID;
    std::string chromosome;
    std::string prev_chr = "";
    std::string A1, A2;
    size_t chr_num = 0;
    uint32_t SNP_position = 0;

    size_t ref_target_match = 0;
    // obtain the context information. We don't check out of bound as
    // that is unlikely to happen
    for (size_t i_snp = 0; i_snp < num_snp; ++i_snp)
    {
        // go through each SNP in the file
        if (i_snp < 1000 && !m_reporter->unit_testing())
        {
            fprintf(stderr, "\r%zu SNPs processed in %s\r", i_snp,
                    name.c_str());
        }
        else if (i_snp % 1000 == 0 && !m_reporter->unit_testing())
        {
            fprintf(stderr, "\r%zuK SNPs processed in %s   \r", i_snp / 1000,
                    name.c_str());
        }
        ++m_unfiltered_marker_ct;
        // directly use the library without decompressing the genotype
        read_snp_identifying_data(*bgen_file, context, &SNPID, &RSID,
                                  &chromosome, &SNP_position, &A1, &A2);
        // get the current location of bgen file, this will be used to
        // skip to current location later on
        if (check_chr(chromosome, prev_chr, chr_num, chr_error, sex_error))
        {
            if (process_snp(
                    exclusion_regions, mismatch_snp_record_name,
                    mismatch_source, SNPID,
                    std::make_unique<SNP>(RSID, chr_num, SNP_position, A1, A2,
                                          file_idx, bgen_file->tellg()),
                    processed_snps, duplicated_snps, retain_snp, genotype))
            { ++ref_target_match; }
        }

        // read in the genotype data block so that we advance the
        // ifstream pointer to the next SNP entry
        genfile::bgen::ignore_genotype_data_block(*bgen_file, context);
    }
    if (num_snp < 1000 && !m_reporter->unit_testing())
    {
        fprintf(stderr, "\r%zu SNPs processed in %s\r", num_snp, name.c_str());
    }
    else if (num_snp % 1000 == 0 && !m_reporter->unit_testing())
    {
        fprintf(stderr, "\r%zuK SNPs processed in %s   \r", num_snp / 1000,
                name.c_str());
    }
    bgen_file.reset();
    if (!m_reporter->unit_testing()) fprintf(stderr, "\n");
    return ref_target_match;
}

void BinaryGen::gen_snp_vector(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const std::string& out_prefix, Genotype* target)
{
    const std::string mismatch_snp_record_name = out_prefix + ".mismatch";
    std::unordered_set<std::string> duplicated_snps;
    std::unordered_set<std::string> processed_snps;
    auto&& genotype = (m_is_ref) ? target : this;
    std::vector<bool> retain_snp(genotype->m_existed_snps.size(), false);
    size_t total_unfiltered_snps = 0;
    size_t ref_target_match = 0;
    bool chr_error = false, sex_error = false;
    if (!m_is_ref)
    {
        auto bgen_file = misc::load_stream(
            m_genotype_file_names.front() + ".bgen", std::ios_base::binary);
        check_sample_consistent(m_context_map[0], *bgen_file);
    }
    // use as a check, just in case we can't handlex
    for (auto context : m_context_map)
    {
        if (std::numeric_limits<size_t>::max() - total_unfiltered_snps
            <= context.number_of_variants)
        {
            throw std::runtime_error(
                "Error: Too many SNPs in bgen file (> "
                + std::to_string(std::numeric_limits<size_t>::max())
                + "which is the theoretical upper bound of "
                  "number we can store");
        }
        total_unfiltered_snps += context.number_of_variants;
    }
    for (size_t file_idx = 0; file_idx < m_genotype_file_names.size();
         ++file_idx)
    {
        // now start reading each bgen file
        auto bgen_file = misc::load_stream(
            m_genotype_file_names[file_idx] + ".bgen", std::ios_base::binary);
        ref_target_match += transverse_bgen_for_snp(
            exclusion_regions, mismatch_snp_record_name, file_idx,
            std::move(bgen_file), duplicated_snps, processed_snps, retain_snp,
            chr_error, sex_error, genotype);
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
                                    const std::string& prefix,
                                    Genotype* genotype)
{
    const std::string intermediate_name = prefix + ".inter";
    std::vector<bool> retain_snps(genotype->m_existed_snps.size(), false);
    std::streampos byte_pos, tmp_byte_pos;
    size_t processed_count = 0;
    size_t cur_file_idx = 0;
    size_t retained = 0;
    uint32_t ref_count = 0;
    uint32_t het_count = 0;
    uint32_t alt_count = 0;
    uint32_t missing_count = 0;
    // we initialize the plink converter with the sample inclusion vector
    // and also the tempory genotype vector list. We also provide the hard
    // coding threshold
    // TODO: This isn't correct if there are non-founder samples in our datas
    // and if we account for ref and target, we also need to consider situation
    // where we use target as reference.
    PLINK_generator setter(m_calculate_prs.data(), m_tmp_genotype.data(),
                           m_hard_threshold, m_dose_threshold);
    // now consider if we are generating the intermediate file
    std::ofstream inter_out;
    if (m_intermediate)
    {
        auto flag = std::ios::binary;
        if (m_is_ref)
        {
            // always append with reference as worst case will be target
            // generating an empty file
            flag |= std::ios::app;
        }
        inter_out.open(intermediate_name.c_str(), flag);
    }
    // now start processing the bgen file
    double progress = 0, prev_progress = -1.0;
    const size_t total_snp = genotype->m_existed_snps.size();
    for (auto&& snp : genotype->m_existed_snps)
    {
        progress = static_cast<double>(processed_count)
                   / static_cast<double>(total_snp) * 100;
        if (!m_reporter->unit_testing() && progress - prev_progress > 0.01)
        {
            fprintf(stderr, "\rCalculating allele frequencies: %03.2f%%",
                    progress);
            prev_progress = progress;
        }
        snp->get_file_info(cur_file_idx, byte_pos, m_is_ref);
        // now read in the genotype information
        genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
            m_genotype_file, m_genotype_file_names[cur_file_idx] + ".bgen",
            m_context_map[cur_file_idx], setter, &m_buffer1, &m_buffer2,
            byte_pos);
        // no founder, much easier
        setter.get_count(ref_count, het_count, alt_count, missing_count);
        ++processed_count;
        if (filter_snp(ref_count, het_count, alt_count, ref_count, het_count,
                       alt_count, filter_info.geno, filter_info.maf,
                       missing_count))
        { continue; }

        if (setter.info_score(filter_info.info_type) < filter_info.info_score)
        {
            ++m_num_info_filter;
            continue;
        }
        // if we can reach here, it is not removed
        snp->set_counts(ref_count, het_count, alt_count, missing_count,
                        m_is_ref);
        snp->set_expected(setter.expected(), m_is_ref);
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
            inter_out.write(reinterpret_cast<char*>(m_tmp_genotype.data()),
                            m_tmp_genotype.size() * sizeof(uintptr_t));
            if (!m_is_ref)
            {
                // target file
                if (m_hard_coded)
                {
                    m_target_plink = true;
                    snp->update_file(m_genotype_file_names.size(), tmp_byte_pos,
                                     false);
                }
                if (!m_expect_reference)
                {
                    // we don't have reference, so use target as reference
                    m_ref_plink = true;
                    snp->update_file(m_genotype_file_names.size(), tmp_byte_pos,
                                     true);
                }
            }
            else
            {
                // this is the reference file
                m_ref_plink = true;
                snp->update_file(m_genotype_file_names.size(), tmp_byte_pos,
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
    if (!m_reporter->unit_testing())
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
    std::unique_ptr<PRS_Interpreter> setter = nullptr;
    if (not_first)
    {
        setter = std::make_unique<Add_PRS>(&prs_list, &m_calculate_prs,
                                           m_prs_calculation.missing_score);
    }
    else
    {
        setter = std::make_unique<First_PRS>(&prs_list, &m_calculate_prs,
                                             m_prs_calculation.missing_score);
    }
    std::vector<size_t>::const_iterator cur_idx = start_idx;
    size_t file_idx;
    std::streampos byte_pos;
    for (; cur_idx != end_idx; ++cur_idx)
    {
        auto&& snp = m_existed_snps[(*cur_idx)];
        std::tie(file_idx, byte_pos) = snp->get_file_info(m_is_ref);
        // if the file name differ, or the file isn't open, we will open it
        auto&& context = m_context_map[file_idx];
        setter->set_stat(snp->stat(), m_homcom_weight, m_het_weight,
                         m_homrar_weight, snp->is_flipped());
        // start performing the parsing
        genfile::bgen::read_and_parse_genotype_data_block<PRS_Interpreter>(
            m_genotype_file, m_genotype_file_names[file_idx] + ".bgen", context,
            *setter, &m_buffer1, &m_buffer2, byte_pos);
        if (!not_first)
        {
            setter.reset(new Add_PRS(&prs_list, &m_calculate_prs,
                                     m_prs_calculation.missing_score));
        }
        not_first = true;
    }
    setter.reset();
}


void BinaryGen::hard_code_score(
    std::vector<PRS>& prs_list,
    const std::vector<size_t>::const_iterator& start_idx,
    const std::vector<size_t>::const_iterator& end_idx, bool reset_zero)
{
    // genotype counts
    uint32_t homrar_ct = 0;
    uint32_t missing_ct = 0;
    uint32_t het_ct = 0;
    uint32_t homcom_ct = 0;
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
    genfile::bgen::Context context;
    PLINK_generator setter(m_calculate_prs.data(), m_tmp_genotype.data(),
                           m_hard_threshold, m_dose_threshold);
    std::vector<size_t>::const_iterator cur_idx = start_idx;
    uintptr_t* genotype_ptr;
    for (; cur_idx != end_idx; ++cur_idx)
    {
        auto&& cur_snp = m_existed_snps[(*cur_idx)];
        if (cur_snp->current_genotype() == nullptr)
        {
            auto [idx, byte_pos] = cur_snp->get_file_info(m_is_ref);
            if (m_intermediate)
            {
                if (!cur_snp->get_counts(homcom_ct, het_ct, homrar_ct,
                                         missing_ct,
                                         m_prs_calculation.use_ref_maf))
                {
                    throw std::logic_error(
                        "Error: Sam has a logic error in bgen");
                }
                // Have intermediate file and have the counts
                // read in the genotype information to the genotype vector
                const uintptr_t unfiltered_sample_ct4 =
                    (m_unfiltered_sample_ct + 3) / 4;
                m_genotype_file.read(
                    m_genotype_file_names[idx], byte_pos, unfiltered_sample_ct4,
                    reinterpret_cast<char*>(m_tmp_genotype.data()));
            }
            else
            {
                // now read in the genotype information
                context = m_context_map[idx];
                // start performing the parsing
                genfile::bgen::read_and_parse_genotype_data_block<
                    PLINK_generator>(
                    m_genotype_file, m_genotype_file_names[idx] + ".bgen",
                    context, setter, &m_buffer1, &m_buffer2, byte_pos);
                if (!m_prs_calculation.use_ref_maf)
                {
                    setter.get_count(homcom_ct, het_ct, homrar_ct, missing_ct);
                }
                else
                {
                    if (!cur_snp->get_counts(homcom_ct, het_ct, homrar_ct,
                                             missing_ct,
                                             m_prs_calculation.use_ref_maf))
                    {
                        throw std::runtime_error(
                            "Error: Sam forgot to load genotype count from "
                            "reference");
                    }
                }
            }
            genotype_ptr = m_tmp_genotype.data();
        }
        else
        {
            if (!cur_snp->get_counts(homcom_ct, het_ct, homrar_ct, missing_ct,
                                     m_prs_calculation.use_ref_maf))
            { throw std::runtime_error("Error: Sam has a logic error"); }
            genotype_ptr = cur_snp->current_genotype();
        }
        homcom_weight = m_homcom_weight;
        het_weight = m_het_weight;
        homrar_weight = m_homrar_weight;
        if (cur_snp->is_flipped()) { std::swap(homcom_weight, homrar_weight); }
        maf = 1.0
              - static_cast<double>(homcom_weight * homcom_ct
                                    + het_ct * het_weight
                                    + homrar_weight * homrar_ct)
                    / (static_cast<double>(homcom_ct + het_ct + homrar_ct)
                       * ploidy);

        stat = cur_snp->stat();
        adj_score = 0;
        if (is_centre) { adj_score = ploidy * stat * maf; }
        miss_score = 0;
        if (mean_impute) { miss_score = ploidy * stat * maf; }
        read_prs(genotype_ptr, prs_list, ploidy, stat, adj_score, miss_score,
                 miss_count, homcom_weight, het_weight, homrar_weight,
                 not_first);
        not_first = true;
    }
}

void BinaryGen::count_and_read_genotype(const std::unique_ptr<SNP>& snp)
{
    auto [file_idx, byte_pos] = snp->get_file_info(false);
    auto&& genotype = snp->current_genotype();
    // load into memory is useless for dosage score
    if (!m_hard_coded) return;

    if (m_intermediate)
    {
        // this is the intermediate
        uint32_t homrar_ct = 0;
        uint32_t missing_ct = 0;
        uint32_t het_ct = 0;
        uint32_t homcom_ct = 0;
        if (!snp->get_counts(homcom_ct, het_ct, homrar_ct, missing_ct,
                             m_prs_calculation.use_ref_maf))
        { throw std::logic_error("Error: Sam has a logic error in bgen"); }
        const uintptr_t unfiltered_sample_ct4 =
            (m_unfiltered_sample_ct + 3) / 4;
        m_genotype_file.read(m_genotype_file_names[file_idx], byte_pos,
                             unfiltered_sample_ct4,
                             reinterpret_cast<char*>(genotype));
    }
    else
    {

        // start performing the parsing
        PLINK_generator setter(m_calculate_prs.data(), genotype,
                               m_hard_threshold, m_dose_threshold);
        genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
            m_genotype_file, m_genotype_file_names[file_idx] + ".bgen",
            m_context_map[file_idx], setter, &m_buffer1, &m_buffer2, byte_pos);
    }
}

void BinaryGen::read_score(std::vector<PRS>& prs_list,
                           const std::vector<size_t>::const_iterator& start_idx,
                           const std::vector<size_t>::const_iterator& end_idx,
                           bool reset_zero)
{
    if (m_hard_coded)
    { hard_code_score(prs_list, start_idx, end_idx, reset_zero); }
    else
    {
        dosage_score(prs_list, start_idx, end_idx, reset_zero);
    }
}
