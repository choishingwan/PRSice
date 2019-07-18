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

BinaryGen::BinaryGen(const std::string& list_file, const std::string& file,
                     const std::string& pheno_file,
                     const std::string& out_prefix, const double hard_threshold,
                     const double dose_threshold, const size_t thread,
                     const bool use_inter, const bool use_hard_coded,
                     const bool no_regress, const bool ignore_fid,
                     const bool keep_nonfounder, const bool keep_ambig,
                     const bool is_ref, Reporter& reporter)
{
    m_intermediate = use_inter;
    m_thread = thread;
    m_ignore_fid = ignore_fid;
    m_keep_nonfounder = keep_nonfounder;
    m_keep_ambig = keep_ambig;
    m_is_ref = is_ref;
    m_hard_coded = use_hard_coded;
    m_intermediate_file = out_prefix + ".inter";
    m_hard_threshold = hard_threshold;
    m_dose_threshold = dose_threshold;
    // set the chromosome information
    // will need to add more script here if we want to support something
    // other than human
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    m_haploid_mask.resize(CHROM_MASK_WORDS, 0);
    // main use of following function is to set the max code
    init_chr();
    std::string message = "Initializing Genotype";
    std::string file_name;
    bool is_list = !(list_file.empty());
    if (!is_list) {
        file_name = file;
    }
    else
    {
        file_name = list_file;
    }
    if (is_list) {
        std::vector<std::string> token = misc::split(file_name, ",");
        bool external_sample = false;
        if (token.size() == 2) {
            m_sample_file = token[1];
            file_name = token[0];
            external_sample = true;
        }
        message.append(" info from file " + file_name + " (bgen)\n");
        if (external_sample) {
            message.append("With external sample file: " + m_sample_file
                           + "\n");
        }
        if (!external_sample && no_regress && pheno_file.empty()) {
            throw std::runtime_error("ERROR: You must provide a phenotype "
                                     "file for bgen format!\n");
        }
        m_genotype_files = load_genotype_prefix(file_name);
        if (!external_sample) {
            m_sample_file = pheno_file;
        }
    }
    else
    {
        std::vector<std::string> token = misc::split(file_name, ",");
        bool external_sample = false;
        if (token.size() == 2) {
            m_sample_file = token[1];
            file_name = token[0];
            external_sample = true;
        }
        message.append(" file: " + file_name + " (bgen)\n");
        if (external_sample) {
            message.append("With external fam file: " + m_sample_file + "\n");
        }
        if (!external_sample && no_regress && pheno_file.empty()) {
            throw std::runtime_error("ERROR: You must provide a phenotype "
                                     "file for bgen format!\n");
        }
        m_genotype_files = set_genotype_files(file_name);
        if (!external_sample) {
            m_sample_file = pheno_file;
        }
    }
    reporter.report(message);
}


std::vector<Sample_ID> BinaryGen::gen_sample_vector(const std::string& delim)
{
    // first, check if the sample file is in the sample format specified by BGEN
    bool is_sample_format = check_is_sample_format();
    std::ifstream sample_file(m_sample_file.c_str());
    if (!sample_file.is_open()) {
        std::string error_message =
            "Error: Cannot open sample file: " + m_sample_file;
        throw std::runtime_error(error_message);
    }

    std::string line;
    int sex_col = -1;
    // now check if there's a sex information
    if (is_sample_format) {
        // only do this if the file is sample format
        std::getline(sample_file, line);
        std::vector<std::string> header_names = misc::split(line);
        // might want to get the reporter involved.
        fprintf(stderr, "Detected bgen sample file format\n");
        for (size_t i = 3; i < header_names.size(); ++i) {
            // try to identify column containing SEX
            std::transform(header_names[i].begin(), header_names[i].end(),
                           header_names[i].begin(), ::toupper);
            if (header_names[i].compare("SEX") == 0) {
                sex_col = static_cast<int>(i);
                break;
            }
        }
        // now we read in the second line
        std::getline(sample_file, line);
        if (sex_col != -1) {
            // double check if the format is alright
            std::vector<std::string> header_format = misc::split(line);
            // no need range check as we know for sure sex_col must be in range
            if (header_format[static_cast<std::vector<std::string>::size_type>(
                    sex_col)]
                != "D")
            {
                // it is not the expected format, simply ignore the sex
                // information as we don't need it for any crucial
                // calculation anyway
                std::string error_message =
                    "Warning: Sex must be coded as "
                    "\"D\" in bgen sample file!\n"
                    "We will ignore the sex information.";
                std::cerr << error_message << "\n";
                sex_col = -1;
            }
        }
    }
    // now start reading the file
    int line_id = 0;
    std::vector<Sample_ID> sample_name;
    std::unordered_set<std::string> duplicated_samples;
    std::vector<std::string> duplicated_sample_id;
    std::vector<std::string> token;
    std::vector<bool> temp_inclusion_vec;
    bool inclusion = false;
    std::string FID, IID;
    while (std::getline(sample_file, line)) {
        misc::trim(line);
        line_id++;
        if (!line.empty()) {
            token = misc::split(line);
            // if it is not the sample file, check if this has a header
            // not the best way, but will do it
            if (line_id == 1) {
                std::string header_test = token[0];
                std::transform(header_test.begin(), header_test.end(),
                               header_test.begin(), ::toupper);
                if (header_test == "FID"
                    || (header_test == "IID" && m_ignore_fid))
                {
                    // this is the header, skip
                    continue;
                }
                else
                {
                    // emit a warning so people might be aware of it
                    std::cerr
                        << "We assume the following line is not a header:\n"
                        << line << "\n(first column isn't FID or IID)\n";
                }
            }
            if (static_cast<int>(token.size())
                < ((sex_col != -1) ? (sex_col) : (1 + !m_ignore_fid)))
            {
                std::string error_message =
                    "Error: Line " + std::to_string(line_id)
                    + " must have at least "
                    + std::to_string((sex_col != -1) ? (sex_col)
                                                     : (1 + !m_ignore_fid))
                    + " columns! Number of column="
                    + misc::to_string(token.size());
                throw std::runtime_error(error_message);
            }
            m_unfiltered_sample_ct++;

            if (is_sample_format || !m_ignore_fid) {
                FID = token[0];
                IID = token[1];
            }
            else
            {
                // not a sample format or used ignore fid
                FID = "";
                IID = token[0];
            }
            std::string id = (m_ignore_fid) ? IID : FID + delim + IID;
            // we assume all bgen samples are founders
            if (!m_remove_sample) {
                inclusion = (m_sample_selection_list.find(id)
                             != m_sample_selection_list.end());
            }
            else
            {
                inclusion = (m_sample_selection_list.find(id)
                             == m_sample_selection_list.end());
            }
            if (sex_col != -1) {
                try
                {
                    int sex_info = misc::convert<int>(
                        token[static_cast<std::vector<std::string>::size_type>(
                            sex_col)]);
                    m_num_male += (sex_info == 1);
                    m_num_female += (sex_info == 2);
                    m_num_ambig_sex += (sex_info != 1 && sex_info != 2);
                }
                catch (...)
                {
                    throw std::runtime_error(
                        "Error: Non-numeric sex coding!\n");
                }
            }
            else
            {
                m_num_ambig_sex++;
            }
            if (duplicated_samples.find(id) != duplicated_samples.end()) {
                duplicated_sample_id.push_back(id);
            }
            duplicated_samples.insert(id);
            temp_inclusion_vec.push_back(inclusion);
            if (!m_is_ref && inclusion) {
                // all sample must be a founder
                sample_name.emplace_back(Sample_ID(FID, IID, "", true));
            }
        }
    }
    if (!duplicated_sample_id.empty()) {
        // TODO: Produce a file containing id of all valid samples
        std::string error_message =
            "Error: A total of " + misc::to_string(duplicated_sample_id.size())
            + " duplicated samples detected! Please ensure all samples "
              "have an "
              "unique identifier";
        throw std::runtime_error(error_message);
    }


    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    m_sample_include.resize(unfiltered_sample_ctl, 0);
    m_founder_info.resize(unfiltered_sample_ctl, 0);
    m_founder_ct = 0;
    for (std::vector<bool>::size_type i = 0; i < temp_inclusion_vec.size(); ++i)
    {
        // using temp_inclusion_vec allow us to set the sample_include vector
        // and founder_info vector without generating the sample vector
        if (temp_inclusion_vec[i]) {
            ++m_sample_ct;
            SET_BIT(i, m_sample_include.data());
            // we assume all bgen samples to be founder
            SET_BIT(i, m_founder_info.data());
        }
    }
    m_founder_ct = m_sample_ct;
    sample_file.close();
    // m_prs_info.reserve(m_sample_ct);
    // initialize the PRS vector (can't reserve, otherwise seg fault(no idea
    // why))
    for (size_t i = 0; i < m_sample_ct; ++i) {
        m_prs_info.emplace_back(PRS());
    }
    // initialize regression flag
    m_in_regression.resize(m_sample_include.size(), 0);
    // I don't really like to use this but I also don't want to un-necessarily
    // pass the delimitor to gen_snp
    m_id_delim = delim;
    return sample_name;
}

bool BinaryGen::check_is_sample_format()
{
    // read the sample file
    std::ifstream sample_file(m_sample_file.c_str());
    if (!sample_file.is_open()) {
        std::string error_message =
            "Error: Cannot open sample file: " + m_sample_file;
        throw std::runtime_error(error_message);
    }
    // get the first two line of input
    std::string first_line, second_line;
    std::getline(sample_file, first_line);
    std::getline(sample_file, second_line);
    sample_file.close();
    // split the first two lines
    std::vector<std::string> first_row = misc::split(first_line);
    std::vector<std::string> second_row = misc::split(second_line);
    // each row should have the same number of column
    if (first_row.size() != second_row.size() || first_row.size() < 3) {
        return false;
    }
    // first 3 must be 0
    for (size_t i = 0; i < 3; ++i) {
        if (second_row[i] != "0") return false;
    }
    // DCPB
    for (size_t i = 4; i < second_row.size(); ++i) {
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

void BinaryGen::get_context(std::string& prefix)
{
    // get the context information for the input bgen file
    // most of these codes are copy from the bgen library
    std::string bgen_name = prefix + ".bgen";
    std::ifstream bgen_file(bgen_name.c_str(), std::ifstream::binary);
    if (!bgen_file.is_open()) {
        std::string error_message = "Error: Cannot open bgen file " + bgen_name;
        throw std::runtime_error(error_message);
    }
    // initialize the bgen context object
    genfile::bgen::Context context;
    uint32_t offset, header_size = 0, number_of_snp_blocks = 0,
                     number_of_samples = 0, flags = 0;
    char magic[4];
    std::size_t fixed_data_size = 20;
    std::vector<char> free_data;
    // read in the offset information
    genfile::bgen::read_little_endian_integer(bgen_file, &offset);
    // read in the size of the header
    genfile::bgen::read_little_endian_integer(bgen_file, &header_size);
    // the size of header should be bigger than or equal to the fixed data size
    assert(header_size >= fixed_data_size);
    genfile::bgen::read_little_endian_integer(bgen_file, &number_of_snp_blocks);
    genfile::bgen::read_little_endian_integer(bgen_file, &number_of_samples);
    // now read in the magic number
    bgen_file.read(&magic[0], 4);
    free_data.resize(header_size - fixed_data_size);
    if (free_data.size() > 0) {
        bgen_file.read(&free_data[0], free_data.size());
    }
    // now read in the flag
    genfile::bgen::read_little_endian_integer(bgen_file, &flags);
    if ((magic[0] != 'b' || magic[1] != 'g' || magic[2] != 'e'
         || magic[3] != 'n')
        && (magic[0] != 0 || magic[1] != 0 || magic[2] != 0 || magic[3] != 0))
    {
        throw std::runtime_error("Error: Incorrect magic string!\nPlease "
                                 "check you have provided a valid bgen "
                                 "file!");
    }
    if (bgen_file) {
        context.number_of_samples = number_of_samples;
        context.number_of_variants = number_of_snp_blocks;
        context.magic.assign(&magic[0], &magic[0] + 4);
        context.offset = offset;
        context.flags = flags;
        m_context_map[prefix].offset = context.offset;
        m_context_map[prefix].flags = context.flags;
        m_context_map[prefix].number_of_samples = context.number_of_samples;
        m_context_map[prefix].number_of_variants = context.number_of_variants;

        if ((flags & genfile::bgen::e_CompressedSNPBlocks)
            == genfile::bgen::e_ZstdCompression)
        {
            throw std::runtime_error(
                "Error: zstd compression currently not supported");
        }
    }
    else
    {
        throw std::runtime_error("Error: Problem reading bgen file!");
    }
}

bool BinaryGen::check_sample_consistent(const std::string& bgen_name,
                                        const std::string& delim,
                                        const genfile::bgen::Context& context)
{
    // only do this if our gen file contains the sample information
    if (context.flags & genfile::bgen::e_SampleIdentifiers) {
        std::ifstream bgen_file(bgen_name.c_str(), std::ifstream::binary);
        uint32_t tmp_offset;
        genfile::bgen::Context tmp_context;
        genfile::bgen::read_offset(bgen_file, &tmp_offset);
        genfile::bgen::read_header_block(bgen_file, &tmp_context);
        const bool has_fid = !m_sample_id.front().FID.empty();
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
        if (actual_number_of_samples != context.number_of_samples) {
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
        if (!m_is_ref) {
            size_t sample_vector_idx = 0;
            for (size_t i = 0; i < actual_number_of_samples; ++i) {
                genfile::bgen::read_length_followed_by_data(
                    bgen_file, &identifier_size, &identifier);
                if (!bgen_file)
                    throw std::runtime_error(
                        "Error: Problem reading bgen file!");
                bytes_read += sizeof(identifier_size) + identifier_size;
                // Need to double check. BGEN format might differ depends
                // if FID is provided. When FID is provided, then the ID
                // should be FID + delimitor + IID; otherwise it'd be IID
                if(IS_SET(m_sample_include.data(), i)){
                    if (m_sample_id[sample_vector_idx].IID != identifier
                            && (m_sample_id[sample_vector_idx].FID + delim +
                                m_sample_id[sample_vector_idx].IID)
                            != identifier)
                    {
                        std::string error_message =
                                "Error: Sample mismatch "
                                "between bgen and phenotype file! Name in BGEN "
                                "file is "
                                ":"
                                + identifier + " and in phentoype file is: ";
                        if (has_fid)
                            error_message.append(m_sample_id[sample_vector_idx].FID + delim
                                                 + m_sample_id[sample_vector_idx].IID);
                        else
                            error_message.append(m_sample_id[sample_vector_idx].IID);
                        error_message.append(
                                    ". Please note that PRSice require the bgen file and "
                                    "the .sample (or phenotype file if sample file is "
                                    "not provided) to have sample in the same order. (We "
                                    "might be able to losen this requirement in future "
                                    "when we have more time)");
                        throw std::runtime_error(error_message);
                    }
                    sample_vector_idx++;
                }
            }
        }
        assert(bytes_read == sample_block_size);
    }
    return true;
}

void BinaryGen::gen_snp_vector(
    const std::vector<IITree<int, int>>& exclusion_regions,
    const std::string& out_prefix, Genotype* target)
{
    std::unordered_set<std::string> duplicated_snps;
    // should only apply to SNPs that are not removed due to extract/exclude
    std::unordered_set<std::string> processed_snps;
    std::vector<std::string> alleles;
    std::vector<bool> retain_snp;
    auto&& reference = (m_is_ref) ? target : this;
    retain_snp.resize(reference->m_existed_snps.size(), false);
    std::ifstream bgen_file;
    std::ofstream mismatch_snp_record;
    std::string bgen_name;
    std::string allele;
    std::string SNPID;
    std::string RSID;
    std::string cur_id;
    std::string chromosome;
    std::string prev_chr = "";
    std::string file_name;
    std::string mismatch_snp_record_name = out_prefix + ".mismatch";
    std::string error_message = "";
    std::string A1, A2;
    std::streampos byte_pos;
    size_t total_unfiltered_snps = 0;
    size_t ref_target_match = 0;
    uint32_t SNP_position = 0;
    uint32_t offset;
    uint32_t num_snp;
    int chr_code = 0;
    bool exclude_snp = false;
    bool chr_sex_error = false;
    bool chr_error = false;
    bool prev_chr_sex_error = false;
    bool prev_chr_error = false;
    bool flipping = false;
    bool to_remove = false;
    for (auto prefix : m_genotype_files) {
        // go through each genotype file and get the context information
        get_context(prefix);
        // get the total unfiltered snp size so that we can initalize the vector
        total_unfiltered_snps += m_context_map[prefix].number_of_variants;
    }
    check_sample_consistent(std::string(m_genotype_files.front() + ".bgen"),
                            m_id_delim,
                            m_context_map[m_genotype_files.front()]);
    // we don't need to reserve the vector now, as we have already
    // read in the base file
    // it does however mean that we can get into trouble if the base
    // is rather big
    for (auto prefix : m_genotype_files) {
        // now start reading each bgen file
        bgen_name = prefix + ".bgen";
        if (bgen_file.is_open()) bgen_file.close();
        bgen_file.clear();
        bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
        if (!bgen_file.is_open()) {
            std::string error_message =
                "Error: Cannot open bgen file " + bgen_name;
            throw std::runtime_error(error_message);
        }
        // read in the offset
        offset = m_context_map[prefix].offset;
        // skip the offest
        bgen_file.seekg(offset + 4);
        num_snp = m_context_map[prefix].number_of_variants;
        // obtain the context information. We don't check out of bound as that
        // is unlikely to happen
        auto&& context = m_context_map[prefix];
        for (size_t i_snp = 0; i_snp < num_snp; ++i_snp) {
            // go through each SNP in the file
            if (i_snp % 1000 == 0) {
                fprintf(stderr, "\r%zuK SNPs processed in %s\r", i_snp / 1000,
                        bgen_name.c_str());
            }
            m_unfiltered_marker_ct++;
            // directly use the library without decompressing the genotype
            read_snp_identifying_data(
                bgen_file, context, &SNPID, &RSID, &chromosome, &SNP_position,
                [&alleles](std::size_t n) { alleles.resize(n); },
                [&alleles](std::size_t i, std::string& allele) {
                    std::transform(allele.begin(), allele.end(), allele.begin(),
                                   ::toupper);
                    alleles.at(i) = allele;
                });
            exclude_snp = false;
            if (chromosome != prev_chr) {
                chr_code = get_chrom_code_raw(chromosome.c_str());
                if (chr_code_check(chr_code, chr_sex_error, chr_error,
                                   error_message))
                {
                    if (chr_error && !prev_chr_error) {
                        std::cerr << error_message << "\n";
                        prev_chr_error = chr_error;
                    }
                    if (chr_sex_error && !prev_chr_sex_error) {
                        std::cerr << error_message << "\n";
                        prev_chr_sex_error = chr_sex_error;
                    }
                    exclude_snp = true;
                }
                if (!exclude_snp) {
                    // only update prev_chr if we want to include this SNP
                    prev_chr = chromosome;
                }
            }

            if (RSID == ".") {
                // when the rs id isn't available,
                // change it to chr:loc coding. Though I highly doubt the
                // usefulness of this feature as it will also require the
                // reference and base to have the same encoding
                RSID = std::to_string(chr_code) + ":"
                       + std::to_string(SNP_position);
            }
            // default to RS
            cur_id = RSID;
            // by this time point, we should always have the
            // m_existed_snps_index propagated with SNPs from the base. So we
            // can first check if the SNP are presented in base
            if (reference->m_existed_snps_index.find(RSID)
                    == reference->m_existed_snps_index.end()
                && reference->m_existed_snps_index.find(SNPID)
                       == reference->m_existed_snps_index.end())
            {
                // this is the reference panel, and the SNP wasn't found in the
                // target doens't matter if we use RSID or SNPID
                exclude_snp = true;
            }
            else if (reference->m_existed_snps_index.find(RSID)
                     == reference->m_existed_snps_index.end())
            {
                // we found the SNPID
                cur_id = SNPID;
            }
            else if (reference->m_existed_snps_index.find(SNPID)
                     == reference->m_existed_snps_index.end())
            {
                // we found the RSID
                cur_id = RSID;
            }


            // user exclude indicate this SNP is discard by user, not because of
            // QC matric
            if (!m_is_ref) {
                if ((!m_exclude_snp
                     && m_snp_selection_list.find(cur_id)
                            == m_snp_selection_list.end()))
                {
                    exclude_snp = true;
                }
            }
            if (processed_snps.find(cur_id) != processed_snps.end()) {
                duplicated_snps.insert(cur_id);
                exclude_snp = true;
            }
            // perform check on ambiguousity
            else if (ambiguous(alleles.back(), alleles.front()))
            {
                m_num_ambig++;
                if (!m_keep_ambig) exclude_snp = true;
            }

            to_remove = Genotype::within_region(exclusion_regions, chr_code,
                                                SNP_position);
            if (to_remove) {
                m_num_xrange++;
                exclude_snp = true;
            }
            // get the current location of bgen file, this will be used to skip
            // to current location later on
            byte_pos = bgen_file.tellg();
            // read in the genotype data block so that we advance the ifstream
            // pointer to the next SNP entry
            read_genotype_data_block(bgen_file, context, &m_buffer1);
            // if we want to exclude this SNP, we will not perform
            // decompression
            if (!exclude_snp) {
                A1 = alleles.front();
                A2 = alleles.back();
                auto&& target_index = reference->m_existed_snps_index[cur_id];
                if (!reference->m_existed_snps[target_index].matching(
                        chr_code, SNP_position, A1, A2, flipping))
                {
                    // SNP not matched
                    if (!mismatch_snp_record.is_open()) {
                        // open the file accordingly
                        if (m_mismatch_file_output) {
                            mismatch_snp_record.open(
                                mismatch_snp_record_name.c_str(),
                                std::ofstream::app);
                            if (!mismatch_snp_record.is_open()) {
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
                            if (!mismatch_snp_record.is_open()) {
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
                    if (m_is_ref) {
                        mismatch_snp_record
                            << "Reference\t" << cur_id << "\t"
                            << target->m_existed_snps[target_index].chr()
                            << "\t" << chr_code << "\t"
                            << target->m_existed_snps[target_index].loc()
                            << "\t" << SNP_position << "\t"
                            << target->m_existed_snps[target_index].ref()
                            << "\t" << A1 << "\t"
                            << target->m_existed_snps[target_index].alt()
                            << "\t" << A2 << "\n";
                    }
                    else
                    {
                        mismatch_snp_record
                            << "Base\t" << cur_id << "\t" << chr_code << "\t"
                            << m_existed_snps[target_index].chr() << "\t"
                            << SNP_position << "\t"
                            << m_existed_snps[target_index].loc() << "\t" << A1
                            << "\t" << m_existed_snps[target_index].ref()
                            << "\t" << A2 << "\t"
                            << m_existed_snps[target_index].alt() << "\n";
                    }
                    m_num_ref_target_mismatch++;
                }
                else
                {
                    processed_snps.insert(cur_id);
                    if (m_is_ref) {
                        target->m_existed_snps[target_index].add_reference(
                            prefix, byte_pos, flipping);
                    }
                    else
                    {
                        m_existed_snps[target_index].add_target(
                            prefix, byte_pos, chr_code, SNP_position, A1, A2,
                            flipping);
                    }
                    retain_snp[target_index] = true;
                    ref_target_match++;
                }
            }
        }
        bgen_file.close();
        fprintf(stderr, "\n");
    }
    if (ref_target_match != reference->m_existed_snps.size()) {
        // there are mismatch, so we need to update the snp vector
        reference->m_existed_snps.erase(
            std::remove_if(
                reference->m_existed_snps.begin(),
                reference->m_existed_snps.end(),
                [&retain_snp, &reference](const SNP& s) {
                    return !retain_snp[&s - &*begin(reference->m_existed_snps)];
                }),
            reference->m_existed_snps.end());
        reference->m_existed_snps.shrink_to_fit();
        // now update the SNP vector index
        reference->update_snp_index();
    }
    if (duplicated_snps.size() != 0) {
        // there are duplicated SNPs
        std::ofstream log_file_stream;
        std::string dup_name = out_prefix + ".valid";
        log_file_stream.open(dup_name.c_str());
        if (!log_file_stream.is_open()) {
            std::string error_message = "Error: Cannot open file: " + dup_name;
            throw std::runtime_error(error_message);
        }
        for (auto&& snp : reference->m_existed_snps) {
            if (duplicated_snps.find(snp.rs()) != duplicated_snps.end())
                continue;
            log_file_stream << snp.rs() << "\t" << snp.chr() << "\t"
                            << snp.loc() << "\t" << snp.ref() << "\t"
                            << snp.alt() << "\n";
        }

        log_file_stream.close();
        std::string error_message =
            "Error: A total of " + std::to_string(duplicated_snps.size())
            + " duplicated SNP ID detected out of "
            + std::to_string(reference->m_existed_snps.size())
            + " input SNPs!. Valid SNP ID stored at " + dup_name
            + ". You can avoid this error by using --extract " + dup_name;
        throw std::runtime_error(error_message);
    }
}


void BinaryGen::calc_freq_gen_inter(
    const double& maf_threshold, const double& geno_threshold,
    const double& info_threshold, const bool maf_filter, const bool geno_filter,
    const bool info_filter, const bool hard_coded, Genotype* target)
{
    std::vector<bool> retain_snp;
    auto&& reference = (m_is_ref) ? target : this;
    retain_snp.resize(reference->m_existed_snps.size(), false);
    if (!m_is_ref) {
        std::sort(begin(reference->m_existed_snps),
                  end(reference->m_existed_snps),
                  [](SNP const& t1, SNP const& t2) {
                      if (t1.file_name().compare(t2.file_name()) == 0) {
                          return t1.byte_pos() < t2.byte_pos();
                      }
                      else
                          return t1.file_name().compare(t2.file_name()) < 0;
                  });
    }
    else
    {
        // sortby reference positions
        std::sort(begin(reference->m_existed_snps),
                  end(reference->m_existed_snps),
                  [](SNP const& t1, SNP const& t2) {
                      if (t1.file_name().compare(t2.file_name()) == 0) {
                          return t1.ref_byte_pos() < t2.ref_byte_pos();
                      }
                      else
                          return t1.file_name().compare(t2.file_name()) < 0;
                  });
    }
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
    std::vector<bool> retain_snps(reference->m_existed_snps.size(), false);
    std::string prev_file = "";
    std::string cur_file_name = "";
    std::string bgen_name = "";
    std::ifstream bgen_file;
    double cur_maf, cur_geno;
    double sample_ct_recip =
        1.0 / (static_cast<double>(static_cast<int32_t>(m_sample_ct)));
    std::streampos byte_pos, tmp_byte_pos;
    size_t processed_count = 0;
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
    init_quaterarr_from_bitarr(m_sample_include.data(), m_unfiltered_sample_ct,
                               sample_include2.data());
    init_quaterarr_from_bitarr(m_founder_info.data(), m_unfiltered_sample_ct,
                               founder_include2.data());
    m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    // we initialize the plink converter with the sample inclusion vector and
    // also the tempory genotype vector list. We also provide the hard coding
    // threshold
    PLINK_generator setter(&m_sample_include, m_tmp_genotype.data(),
                           m_hard_threshold, m_dose_threshold);
    // now consider if we are generating the intermediate file
    std::ofstream inter_out;
    if (m_intermediate) {
        // allow generation of intermediate file
        if (m_target_plink && m_is_ref) {
            // target already generated some intermediate, now append for
            // reference
            inter_out.open(m_intermediate_file.c_str(),
                           std::ios::binary | std::ios::app);
        }
        else
        {
            // a new intermediate file
            inter_out.open(m_intermediate_file.c_str(), std::ios::binary);
        }
    }
    // now start processing the bgen file
    double progress = 0, prev_progress = -1.0;
    const size_t total_snp = reference->m_existed_snps.size();
    for (auto&& snp : reference->m_existed_snps) {
        progress = static_cast<double>(processed_count)
                   / static_cast<double>(total_snp) * 100;
        if (progress - prev_progress > 0.01) {
            fprintf(stderr, "\rCalculating allele frequencies: %03.2f%%",
                    progress);
            prev_progress = progress;
        }
        if (m_is_ref) {
            cur_file_name = snp.ref_file_name();
            byte_pos = snp.ref_byte_pos();
        }
        else
        {
            cur_file_name = snp.file_name();
            byte_pos = snp.byte_pos();
        }
        if (prev_file != cur_file_name) {
            bgen_name = cur_file_name + ".bgen";
            prev_file = cur_file_name;
            bgen_file.close();
            bgen_file.clear();
            bgen_file.open(bgen_name.c_str());
            if (!bgen_file.is_open()) {
                std::string error_message =
                    "Error: Cannot open bed file: " + bgen_name + "!\n";
                throw std::runtime_error(error_message);
            }
        }
        processed_count++;
        // bgen always seek as there are always something stored in between
        // the genotype data
        if (!bgen_file.seekg(byte_pos, std::ios_base::beg)) {
            std::string error_message =
                "Error: Cannot read the bgen file (seek): " + bgen_name;
            throw std::runtime_error(error_message);
        }
        // now read in the genotype information
        auto&& context = m_context_map[cur_file_name];
        genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
            bgen_file, context, setter, &m_buffer1, &m_buffer2);
        // there's no founder for bgen, can simply calculate the genotyping rate
        // and maf directly when we parse the probability data
        /*
        single_marker_freqs_and_hwe(
            unfiltered_sample_ctv2, m_tmp_genotype.data(),
            sample_include2.data(), founder_include2.data(), m_sample_ct,
            &ll_ct, &lh_ct, &hh_ct, m_founder_ct, &ll_ctf, &lh_ctf, &hh_ctf);
            */
        setter.get_count(ll_ct, lh_ct, hh_ct, missing);
        uii = ll_ct + lh_ct + hh_ct;
        cur_geno = 1.0 - ((static_cast<int32_t>(uii)) * sample_ct_recip);
        uii = 2 * (ll_ct + lh_ct + hh_ct);
        // tmp_total = (ll_ctf + lh_ctf + hh_ctf);
        // assert(m_founder_ct >= tmp_total);
        // missing = m_founder_ct - tmp_total;
        if (!uii) {
            cur_maf = 0.5;
        }
        else
        {
            cur_maf = (static_cast<double>(2 * hh_ct + lh_ct))
                      / (static_cast<double>(uii));
            cur_maf = (cur_maf > 0.5) ? 1-cur_maf: cur_maf;
        }
        if(snp.rs() == "RSID_149"){
            std::cerr << hh_ct << "\t" << lh_ct << "\t" << ll_ct << std::endl;
            std::cerr << "Maf is: " << cur_maf << std::endl;
        }
        // filter by genotype missingness
        if (geno_filter && geno_threshold < cur_geno) {
            m_num_geno_filter++;
            continue;
        }
        // filter by MAF
        // do not flip the MAF for now, so that we
        // are not confuse later on
        // remove SNP if maf lower than threshold
        if (maf_filter && cur_maf < maf_threshold) {
            m_num_maf_filter++;
            continue;
        }
        else if (maf_filter
                 && (misc::logically_equal(cur_maf, 0.0)
                     || misc::logically_equal(cur_maf, 1.0)))
        {
            // none of the sample contain this SNP
            // still count as MAF filtering (for now)
            m_num_maf_filter++;
            continue;
        }

        if (info_filter && setter.info_score() < info_threshold) {
            m_num_info_filter++;
            continue;
        }
        // if we can reach here, it is not removed
        if (m_is_ref) {
            snp.set_ref_counts(ll_ct, lh_ct, hh_ct, missing);
            snp.set_ref_expected(setter.expected());
        }
        else
        {
            snp.set_counts(ll_ct, lh_ct, hh_ct, missing);
            snp.set_expected(setter.expected());
        }
        retained++;
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
            inter_out.write((char*) (&m_tmp_genotype[0]),
                            m_tmp_genotype.size() * sizeof(m_tmp_genotype[0]));
            if (!m_is_ref) {
                // target file
                if (hard_coded) {
                    m_target_plink = true;
                    snp.update_target(m_intermediate_file, tmp_byte_pos);
                }
                if (!m_expect_reference) {
                    // we don't have reference
                    m_ref_plink = true;
                    snp.update_reference(m_intermediate_file, tmp_byte_pos);
                }
            }
            else
            {
                // this is the reference file
                m_ref_plink = true;
                snp.update_reference(m_intermediate_file, tmp_byte_pos);
            }
        }
    }

    fprintf(stderr, "\rCalculating allele frequencies: %03.2f%%\n", 100.0);
    // now update the vector
    if (retained != reference->m_existed_snps.size()) {
        reference->m_existed_snps.erase(
            std::remove_if(reference->m_existed_snps.begin(),
                           reference->m_existed_snps.end(),
                           [&retain_snps, &reference](const SNP& s) {
                               return !retain_snps[(
                                   &s - &*begin(reference->m_existed_snps))];
                           }),
            reference->m_existed_snps.end());
        reference->m_existed_snps.shrink_to_fit();
    }
}

BinaryGen::~BinaryGen()
{
    if (m_bgen_file.is_open()) m_bgen_file.close();
    if (m_target_plink || m_ref_plink) {
        // if we have constructed the intermediate file, we should remove it
        // to save space (plus that file isn't of any useful format and
        // can't be used by any other problem nor can it be reused)
        std::remove(m_intermediate_file.c_str());
    }
}

void BinaryGen::dosage_score(
    const std::vector<size_t>::const_iterator& start_idx,
    const std::vector<size_t>::const_iterator& end_idx, bool reset_zero,
    const bool use_ref_maf)
{
    // currently, use_ref_maf doesn't work on bgen dosage file
    // main reason is we need expected value instead of
    // the MAF (might want to disable the MAF calculation from the start?)
    m_cur_file = "";
    std::string bgen_name;
    bool not_first = !reset_zero;
    // we initialize the PRS interpretor with the required information.
    // m_prs_info is where we store the PRS information
    // and m_sample_include let us know if the sample is required.
    // m_missing_score will inform us as to how to handle the missingness
    PRS_Interpreter setter(&m_prs_info, &m_sample_include, m_missing_score);
    std::vector<size_t>::const_iterator cur_idx = start_idx;
    for (; cur_idx != end_idx; ++cur_idx) {
        auto&& snp = m_existed_snps[(*cur_idx)];
        // if the file name differ, or the file isn't open, we will open it
        if (snp.file_name() != m_cur_file || !m_bgen_file.is_open()) {
            // open the bgen file if required
            if (m_bgen_file.is_open()) m_bgen_file.close();
            bgen_name = snp.file_name() + ".bgen";
            m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
            if (!m_bgen_file.is_open()) {
                std::string error_message =
                    "Error: Cannot open bgen file: " + snp.file_name();
                throw std::runtime_error(error_message);
            }
            m_cur_file = snp.file_name();
        }
        // For bgen file, we will always perform seek as there are always
        // bunch of information between the genotype dosages
        if (!m_bgen_file.seekg(snp.byte_pos(), std::ios_base::beg)) {
            std::string error_message =
                "Error: Cannot seek within the bgen file: " + m_cur_file
                + "!\n";
            throw std::runtime_error(error_message);
        }
        auto&& context = m_context_map[m_cur_file];
        setter.set_stat(snp.stat(), m_homcom_weight, m_het_weight,
                        m_homrar_weight, snp.is_flipped(), not_first);

        // start performing the parsing
        genfile::bgen::read_and_parse_genotype_data_block<PRS_Interpreter>(
            m_bgen_file, context, setter, &m_buffer1, &m_buffer2);
        // check if this SNP has some non-missing sample, if not, invalidate
        // it
        // after reading in this SNP, we no longer need to reset the PRS
        not_first = true;
    }
}


void BinaryGen::hard_code_score(
    const std::vector<size_t>::const_iterator& start_idx,
    const std::vector<size_t>::const_iterator& end_idx, bool reset_zero,
    const bool use_ref_maf)
{
    // we need to calculate the size of possible vectors
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    uintptr_t* lbptr;
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uintptr_t ulii = 0;
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
    int ploidy = 2;
    const int miss_count =
        static_cast<int>((m_missing_score != MISSING_SCORE::SET_ZERO) * ploidy);
    const bool is_centre = (m_missing_score == MISSING_SCORE::CENTER);
    const bool mean_impute = (m_missing_score == MISSING_SCORE::MEAN_IMPUTE);
    // check if we need to reset the sample's PRS
    bool not_first = !reset_zero;
    double stat, maf, adj_score, miss_score;
    std::streampos byte_pos;
    m_cur_file = "";
    // initialize the data structure for storing the genotype
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    genfile::bgen::Context context;
    PLINK_generator setter(&m_sample_include, genotype.data(), m_hard_threshold,
                           m_dose_threshold);
    std::vector<size_t>::const_iterator cur_idx = start_idx;
    for (; cur_idx != end_idx; ++cur_idx) {
        auto&& cur_snp = m_existed_snps[(*cur_idx)];
        // read in the genotype using the modified load_and_collapse_incl
        // function. m_target_plink will inform the function wheter there's
        // an intermediate file
        // if it has the intermediate file, then we should have already
        // calculated the counts
        if (m_cur_file != cur_snp.file_name()) {
            // If we are processing a new file we will need to read it
            if (m_bgen_file.is_open()) {
                m_bgen_file.close();
            }
            m_cur_file = cur_snp.file_name();
            std::string bgen_name =
                m_cur_file + ((m_intermediate) ? "" : ".bgen");
            m_bgen_file.open(bgen_name.c_str(), std::ios::binary);
            if (!m_bgen_file.is_open()) {
                std::string error_message =
                    "Error: Cannot open bgen file: " + bgen_name
                    + ((m_intermediate) ? "(intermediate)" : "");
                throw std::runtime_error(error_message);
            }
            // reset the m_prev_loc flag to 0
            m_prev_loc = 0;
        }
        byte_pos = cur_snp.byte_pos();
        if ((m_prev_loc != byte_pos)
            && !m_bgen_file.seekg(byte_pos, std::ios_base::beg))
        {
            throw std::runtime_error(
                "Error: Cannot seek within the intermediate file!");
        }
        if (cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct,
                               use_ref_maf))
        {
            // must have intermediate file
            m_bgen_file.read((char*) genotype.data(), unfiltered_sample_ct4);
            m_prev_loc =
                static_cast<std::streampos>(unfiltered_sample_ct4) + byte_pos;
        }
        else
        {
            // now read in the genotype information
            context = m_context_map[m_cur_file];
            // start performing the parsing
            genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
                m_bgen_file, context, setter, &m_buffer1, &m_buffer2);
            setter.get_count(homcom_ct, het_ct, homrar_ct, missing_ct);
        }
        // if we haven't got the count from the genotype matrix, we will
        // need to calculate that in theory, we might not need to do the
        // counting as that is already done when we convert the dosages into
        // the binary genotypes (TODO)
        homcom_weight = m_homcom_weight;
        het_weight = m_het_weight;
        homrar_weight = m_homrar_weight;
        maf =
            static_cast<double>(homcom_weight * homcom_ct + het_ct * het_weight
                                + homrar_weight * homrar_ct)
            / (static_cast<double>(homcom_ct + het_ct + homrar_ct) * ploidy);
        if (cur_snp.is_flipped()) {
            // change the mean to reflect flipping
            maf = 1.0 - maf;
            // swap the weighting
            std::swap(homcom_weight, homrar_weight);
        }
        stat = cur_snp.stat(); // Multiply by ploidy
        // only set these value to the imputed value if we require them
        adj_score = 0;
        if (is_centre) {
            adj_score = ploidy * stat * maf;
        }
        miss_score = 0;
        if (mean_impute) {
            // again, mean_impute is stable, branch prediction should be ok
            miss_score = ploidy * stat * maf;
        }
        // start reading the genotype
        lbptr = genotype.data();
        uii = 0;
        ulii = 0;
        do
        {
            // ulii contain the numeric representation of the current genotype
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct) {
                // this is PLINK, not sure exactly what this is about
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            // ujj sample index of the current genotype block
            ujj = 0;
            while (ujj < BITCT) {
                // go through the whole genotype block
                // ukk is the current genotype
                ukk = (ulii >> ujj) & 3;
                // and the sample index can be calculated as uii+(ujj/2)
                if (uii + (ujj / 2) >= m_sample_ct) {
                    break;
                }
                auto&& sample_prs = m_prs_info[uii + (ujj / 2)];
                // now we will get all genotypes (0, 1, 2, 3)
                if (not_first) {
                    switch (ukk)
                    {
                    default:
                        // true = 1, false = 0
                        sample_prs.num_snp += ploidy;
                        sample_prs.prs += homcom_weight * stat - adj_score;
                        break;
                    case 1:
                        sample_prs.num_snp += ploidy;
                        sample_prs.prs += het_weight * stat - adj_score;
                        break;
                    case 3:
                        sample_prs.num_snp += ploidy;
                        sample_prs.prs += homrar_weight * stat - adj_score;
                        break;
                    case 2:
                        // handle missing sample
                        sample_prs.num_snp += miss_count;
                        sample_prs.prs += miss_score;
                        break;
                    }
                }
                else
                {
                    switch (ukk)
                    {
                    default:
                        // true = 1, false = 0
                        sample_prs.num_snp = ploidy;
                        sample_prs.prs = homcom_weight * stat - adj_score;
                        break;
                    case 1:
                        sample_prs.num_snp = ploidy;
                        sample_prs.prs = het_weight * stat - adj_score;
                        break;
                    case 3:
                        sample_prs.num_snp = ploidy;
                        sample_prs.prs = homrar_weight * stat - adj_score;
                        break;
                    case 2:
                        // handle missing sample
                        sample_prs.num_snp = miss_count;
                        sample_prs.prs = miss_score;
                        break;
                    }
                }
                // ulii &= ~((3 * ONELU) << ujj);
                // as each sample is represented by two byte, we will add 2 to
                // the index
                ujj += 2;
            }
            // uii is the number of samples we have finished so far
            uii += BITCT2;
        } while (uii < m_sample_ct);
        // we've finish processing the first SNP no longer need to reset the
        // PRS
        not_first = true;
    }
}


void BinaryGen::read_score(const std::vector<size_t>::const_iterator& start_idx,
                           const std::vector<size_t>::const_iterator& end_idx,
                           bool reset_zero, const bool use_ref_maf)
{
    // because I don't want to touch the code in dosage_score, we will reset
    // the sample here reset_sample_prs();
    if (m_hard_coded) {
        // for hard coded, we need to check if intermediate file is used
        // instead
        hard_code_score(start_idx, end_idx, reset_zero, use_ref_maf);
    }
    else
    {
        dosage_score(start_idx, end_idx, reset_zero, use_ref_maf);
    }
}
