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
            if (line_id == 1 && !is_sample_format) {
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
                if (IS_SET(m_sample_include.data(), i)) {
                    if (m_sample_id[sample_vector_idx].IID != identifier
                        && (m_sample_id[sample_vector_idx].FID + delim
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
                                m_sample_id[sample_vector_idx].FID + delim
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
                                                static_cast<int>(SNP_position));
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
                            prefix, byte_pos, chr_code,
                            static_cast<int>(SNP_position), A1, A2, flipping);
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
