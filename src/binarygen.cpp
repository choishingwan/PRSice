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


BinaryGen::BinaryGen(const Commander& commander, Reporter& reporter,
                     const bool is_ref)
{
    m_intermediate = commander.use_inter();
    m_thread = static_cast<uint32_t>(commander.thread());
    m_ignore_fid = commander.ignore_fid();
    m_keep_nonfounder = commander.nonfounders();
    m_keep_ambig = commander.keep_ambig();
    m_is_ref = is_ref;
    const bool no_regress = commander.no_regress();
    const std::string pheno_file = commander.pheno_file();
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
            message.append(" info from file " + reference_name + " (bgen)\n");
            if (external_sample) {
                message.append("With external sample file: " + m_sample_file
                               + "\n");
            }
            if (!external_sample && no_regress && pheno_file.empty()) {
                throw std::runtime_error("ERROR: You must provide a phenotype "
                                         "file for bgen format!\n");
            }
            m_genotype_files = load_genotype_prefix(reference_name);
            if (!external_sample) {
                m_sample_file = pheno_file;
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
            message.append(" file: " + reference_name + " (bgen)\n");
            if (external_sample) {
                message.append("With external fam file: " + m_sample_file
                               + "\n");
            }
            if (!external_sample && no_regress && pheno_file.empty()) {
                throw std::runtime_error("ERROR: You must provide a phenotype "
                                         "file for bgen format!\n");
            }
            m_genotype_files = set_genotype_files(reference_name);
            if (!external_sample) {
                m_sample_file = pheno_file;
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
            if (!external_sample && no_regress && pheno_file.empty()) {
                throw std::runtime_error("ERROR: You must provide a phenotype "
                                         "file for bgen format!\n");
            }
            m_genotype_files = load_genotype_prefix(target_name);
            if (!external_sample) {
                m_sample_file = pheno_file;
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
            if (!external_sample && no_regress && pheno_file.empty()) {
                throw std::runtime_error("ERROR: You must provide a phenotype "
                                         "file for bgen format!\n");
            }
            m_genotype_files = set_genotype_files(target_name);
            if (!external_sample) {
                m_sample_file = pheno_file;
            }
        }
    }
    reporter.report(message);
}


std::vector<Sample_ID> BinaryGen::gen_sample_vector()
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
        // read in the second line
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
                    std::cerr << "We assume the following line is not a header:"
                              << line << " (first column isn't FID or IID)\n";
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
                // not a sample format and does not use ignore fid
                FID = "";
                IID = token[0];
            }
            std::string id = (m_ignore_fid) ? IID : FID + "_" + IID;
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
                                        const genfile::bgen::Context& context)
{
    // only do this if our gen file contains the sample information
    if (context.flags & genfile::bgen::e_SampleIdentifiers) {
        std::ifstream bgen_file(bgen_name.c_str(), std::ifstream::binary);
        uint32_t tmp_offset;
        genfile::bgen::Context tmp_context;
        genfile::bgen::read_offset(bgen_file, &tmp_offset);
        genfile::bgen::read_header_block(bgen_file, &tmp_context);

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
            for (size_t i = 0; i < actual_number_of_samples; ++i) {
                genfile::bgen::read_length_followed_by_data(
                    bgen_file, &identifier_size, &identifier);
                if (!bgen_file)
                    throw std::runtime_error(
                        "Error: Problem reading bgen file!");
                bytes_read += sizeof(identifier_size) + identifier_size;
                // Only need to use IID as BGEN doesn't have the FID
                // information
                if (m_sample_id[i].IID != identifier) {
                    std::string error_message =
                        "Error: Sample mismatch "
                        "between bgen and phenotype file! Name in BGEN "
                        "file is "
                        ":"
                        + identifier
                        + " and in phentoype file is: " + m_sample_id[i].IID
                        + ". Please note that PRSice require the bgen file and "
                          "the .sample (or phenotype file if sample file is "
                          "not provided) to have sample in the same order. (We "
                          "might be able to losen this requirement in future "
                          "when we have more time)";
                    throw std::runtime_error(error_message);
                }
            }
        }
        assert(bytes_read == sample_block_size);
    }
    return true;
}
std::vector<SNP> BinaryGen::gen_snp_vector(const Commander& commander,
                                           Region& exclusion, Genotype* target)
{
    const std::string out_prefix = commander.out();
    double maf_threshold = 0.0;
    const bool maf_filter = m_is_ref ? commander.ref_maf(maf_threshold)
                                     : commander.target_maf(maf_threshold);
    double geno_threshold = 0.0;
    const bool geno_filter = m_is_ref ? commander.ref_geno(geno_threshold)
                                      : commander.target_geno(geno_threshold);
    double info_threshold = 0.0;
    const bool info_filter = m_is_ref ? commander.ref_info(info_threshold)
                                      : commander.target_info(info_threshold);
    m_hard_coded = m_is_ref
                       ? commander.ref_hard_threshold(m_hard_threshold)
                       : commander.target_hard_thresholding(m_hard_threshold);
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_sample_ct);

    std::vector<SNP> snp_res;
    std::unordered_set<std::string> duplicated_snps;
    // should only apply to SNPs that are not removed due to extract/exclude
    std::unordered_set<std::string> duplicate_check_list;
    std::vector<std::string> alleles;
    std::vector<bool> ref_retain;
    // we will need to initialize the ref_retain to have size = # SNP in target
    // as we are going to resize the target instead of the vector in ref
    if (m_is_ref) ref_retain.resize(target->m_existed_snps.size(), false);
    std::ifstream bgen_file;
    std::ofstream mismatch_snp_record;
    std::ofstream inter_out;
    std::string bgen_name;
    std::string allele;
    std::string SNPID;
    std::string RSID;
    std::string chromosome;
    std::string prev_chr = "";
    std::string file_name;
    std::string mismatch_snp_record_name = out_prefix + ".mismatch";
    std::string error_message = "";
    double cur_maf;
    std::streampos byte_pos, tmp_byte_pos;
    size_t total_unfiltered_snps = 0;
    size_t ref_target_match = 0;
    uint32_t SNP_position;
    uint32_t offset;
    uint32_t num_snp;
    uint32_t homrar_ct;
    uint32_t missing_ct;
    uint32_t het_ct;
    uint32_t homcom_ct;
    intptr_t nanal;
    int chr_code = 0;
    bool exclude_snp = false;
    bool chr_sex_error = false;
    bool chr_error = false;
    bool prev_chr_sex_error = false;
    bool prev_chr_error = false;
    bool first_bgen_file = true;
    bool user_exclude = false;
    bool dummy;
    m_intermediate_file = out_prefix + ".inter";
    m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    // we initialize the plink converter with the sample inclusion vector and
    // also the tempory genotype vector list. We also provide the hard coding
    // threshold
    PLINK_generator setter(&m_sample_include, m_tmp_genotype.data(),
                           m_hard_threshold);
    m_sample_mask.resize(pheno_nm_ctv2);
    fill_quatervec_55(static_cast<uint32_t>(m_sample_ct), m_sample_mask.data());


    for (auto prefix : m_genotype_files) {
        // go through each genotype file and get the context information
        get_context(prefix);
        if (first_bgen_file) {
            first_bgen_file = false;
            // only use the first bgen file to check sample consistency
            std::string bgen_name = prefix + ".bgen";
            check_sample_consistent(bgen_name, m_context_map[prefix]);
        }
        // get the total unfiltered snp size so that we can initalize the vector
        total_unfiltered_snps += m_context_map[prefix].number_of_variants;
    }
    // first pass to get the total SNP number such that we can speed up the
    // push back might need time to reserve large amount of memory for the
    // large number of SNPs included in bgen
    snp_res.reserve(total_unfiltered_snps);
    // to allow multiple file for one chromosome, we put these variable
    // outside the for loop
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

            // directly use the library

            read_snp_identifying_data(
                bgen_file, context, &SNPID, &RSID, &chromosome, &SNP_position,
                [&alleles](std::size_t n) { alleles.resize(n); },
                [&alleles](std::size_t i, std::string& allele) {
                    std::transform(allele.begin(), allele.end(), allele.begin(),
                                   ::toupper);
                    alleles.at(i) = allele;
                });
            exclude_snp = false;
            // but we will not process anything
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
            // user exclude indicate this SNP is discard by user, not because of
            // QC matric
            user_exclude = false;
            if (!m_is_ref) {
                if ((!m_exclude_snp
                     && m_snp_selection_list.find(RSID)
                            == m_snp_selection_list.end())
                    || (m_exclude_snp
                        && m_snp_selection_list.find(RSID)
                               != m_snp_selection_list.end()))
                {
                    user_exclude = true;
                }
            }
            else if (target->m_existed_snps_index.find(RSID)
                     == target->m_existed_snps_index.end())
            {
                // this is the reference panel, and the SNP wasn't found in the
                // target
                exclude_snp = true;
            }
            if (exclusion.check_exclusion(chr_code,
                                          static_cast<int>(SNP_position)))
            {
                // remove the SNP if it is within the exclusion region
                exclude_snp = true;
            }
            if (duplicate_check_list.find(RSID) != duplicate_check_list.end()) {
                duplicated_snps.insert(RSID);
                exclude_snp = true;
            }
            // perform check on ambiguousity
            else if (ambiguous(alleles.back(), alleles.front()))
            {
                m_num_ambig++;
                if (!m_keep_ambig) exclude_snp = true;
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
                // now filter
                file_name = prefix;
                // while we allow user excluded SNP to come into this function,
                // we will not perform filtering or generate the intermediate
                // for it
                if (!user_exclude
                    && (maf_filter || geno_filter || info_filter
                        || m_intermediate))
                {
                    genfile::bgen::uncompress_probability_data(
                        context, m_buffer1, &m_buffer2);
                    genfile::bgen::parse_probability_data<PLINK_generator>(
                        &(m_buffer2)[0], &(m_buffer2)[0] + m_buffer2.size(),
                        context, setter);
                    // this will read in the plink binary
                    if (maf_filter || geno_filter || info_filter) {
                        // only perform the QC if we need to
                        genovec_3freq(m_tmp_genotype.data(),
                                      m_sample_mask.data(), pheno_nm_ctv2,
                                      &missing_ct, &het_ct, &homcom_ct);
                        nanal = static_cast<intptr_t>(m_sample_ct) - missing_ct;
                        // seems a bit strange with this bit level conversion
                        // TODO: Try to remove the need of casting
                        homrar_ct =
                            static_cast<uint32_t>(nanal) - het_ct - homcom_ct;

                        if (nanal == 0) {
                            // non of the sample contain this SNP
                            // still count as MAF filtering (for now)
                            m_num_maf_filter++;
                            // it is now safe to use continue to break the
                            // processing as we have already advance the
                            // ifstream read pointer
                            continue;
                        }

                        if (geno_filter
                            && static_cast<double>(missing_ct)
                                       / static_cast<double>(m_sample_ct)
                                   > geno_threshold)
                        {
                            m_num_geno_filter++;
                            continue;
                        }

                        cur_maf = (static_cast<double>(het_ct + homrar_ct * 2)
                                   / (static_cast<double>(nanal) * 2.0));
                        if (cur_maf > 0.5) cur_maf = 1.0 - cur_maf;
                        // remove SNP if maf lower than threshold
                        if (maf_filter && cur_maf < maf_threshold) {
                            m_num_maf_filter++;
                            continue;
                        }
                        if (info_filter && setter.info_score() < info_threshold)
                        {
                            m_num_info_filter++;
                            continue;
                        }
                    }
                    if (m_intermediate
                        && (m_is_ref || !m_expect_reference
                            || (!m_is_ref && m_hard_coded)))
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
                                        m_tmp_genotype.size()
                                            * sizeof(m_tmp_genotype[0]));
                        file_name = m_intermediate_file;

                        if (!m_is_ref && m_hard_coded && !m_expect_reference) {
                            // not reference, use hard coding and not expecting
                            // reference
                            m_ref_plink = true;
                            m_target_plink = true;
                        }
                        else if (!m_expect_reference || m_is_ref)
                        {
                            // the reference or we are not expecting a
                            // reference. So the ref field should take up the
                            // bytepos information
                            m_ref_plink = true;
                        }
                    }
                }
                if (!m_is_ref && !user_exclude) {
                    // this is not a reference file
                    m_existed_snps_index[RSID] = snp_res.size();
                    if (m_target_plink) {
                        // we have constructed the
                        byte_pos = tmp_byte_pos;
                        file_name = m_intermediate_file;
                    }
                    else
                    {
                        file_name = prefix;
                    }
                    // we use the last alleles obtained from the BGEN as the
                    // effective allele (A1) and the first alele as the
                    // non-effective allele to match behaviour of PLINK
                    // (In PLINK we give weight to A1, in BGEN we give weight to
                    // the last allele in our diploid scenario)
                    snp_res.emplace_back(SNP(RSID, chr_code, SNP_position,
                                             alleles.back(), alleles.front(),
                                             file_name, byte_pos));
                }
                else
                {
                    auto&& target_index = target->m_existed_snps_index[RSID];
                    // this is a reference file
                    if (!user_exclude
                        && !target->m_existed_snps[target_index].matching(
                               chr_code, SNP_position, alleles.back(),
                               alleles.front(), dummy))
                    {
                        // The SNP in reference does not match the one in target
                        // and this is not a user exclude SNP
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
                            << "Reference\t" << RSID << "\t"
                            << target->m_existed_snps[target_index].chr()
                            << "\t" << chr_code << "\t"
                            << target->m_existed_snps[target_index].loc()
                            << "\t" << SNP_position << "\t"
                            << target->m_existed_snps[target_index].ref()
                            << "\t" << alleles.front() << "\t"
                            << target->m_existed_snps[target_index].alt()
                            << "\t" << alleles.back() << "\n";

                        m_num_ref_target_mismatch++;
                    }
                    else if (!user_exclude)
                    {
                        // the information matched between the reference and
                        // target and this is not a user exclude SNP
                        // therefore we allow the check of duplication
                        duplicate_check_list.insert(RSID);
                        if (m_ref_plink) {
                            byte_pos = tmp_byte_pos;
                            file_name = m_intermediate_file;
                        }
                        else
                        {
                            file_name = prefix;
                        }
                        target->m_existed_snps[target_index].add_reference(
                            file_name, byte_pos);
                        ref_retain[target_index] = true;
                        ref_target_match++;
                    }
                }
            }
        }
        bgen_file.close();
        fprintf(stderr, "\n");
    }
    // hope that this will release some memory
    // only useful for target file
    snp_res.shrink_to_fit();
    if (m_is_ref && ref_target_match != target->m_existed_snps.size()) {
        // there are mismatch between reference and target, so we need to update
        // the target snp vector
        target->m_existed_snps.erase(
            std::remove_if(
                target->m_existed_snps.begin(), target->m_existed_snps.end(),
                [&ref_retain, &target](const SNP& s) {
                    return !ref_retain[&s - &*begin(target->m_existed_snps)];
                }),
            target->m_existed_snps.end());
        target->m_existed_snps.shrink_to_fit();
        // When reading the reference file, we actually update the
        // SNP list in target file. This lead to the index search in
        // target to have the wrong index. To avoid that, we need to
        // update the SNP index accordingly
        target->update_snp_index();
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
        for (auto&& snp : (m_is_ref ? target->m_existed_snps : snp_res)) {
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
            + std::to_string(snp_res.size())
            + " input SNPs!. Valid SNP ID stored at " + dup_name
            + ". You can avoid this error by using --extract " + dup_name;
        throw std::runtime_error(error_message);
    }
    if (m_hard_coded) {
        uintptr_t unfiltered_sample_ctl =
            BITCT_TO_WORDCT(m_unfiltered_sample_ct);
        // ensure the tmp_genotype is initialized with the correct size if we
        // are going to perform hard-coding
        m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    }
    return snp_res;
}

BinaryGen::~BinaryGen()
{
    if (m_bgen_file.is_open()) m_bgen_file.close();
    if (m_target_plink || m_ref_plink) {
        // if we have constructed the intermediate file, we should remove it to
        // save space (plus that file isn't of any useful format and can't be
        // used by any other problem nor can it be reused)
        std::remove(m_intermediate_file.c_str());
    }
}


void BinaryGen::dosage_score(const size_t start_index, const size_t end_bound,
                             const size_t region_index, bool set_zero)
{
    m_cur_file = "";
    std::string bgen_name;

    bool not_first = !set_zero;
    // we initialize the PRS interpretor with the required information.
    // m_prs_info is where we store the PRS information
    // and m_sample_include let us know if the sample is required.
    // m_missing_score will inform us as to how to handle the missingness
    PRS_Interpreter setter(&m_prs_info, &m_sample_include, m_missing_score);
    for (size_t i_snp = start_index; i_snp < end_bound; ++i_snp) {
        auto&& snp = m_existed_snps[i_snp];
        // skip SNPs that are not in the region or that are invalid
        if (!snp.in(region_index) || !snp.valid()) continue;
        if (snp.file_name() == m_cur_file) {
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
        // For bgen file, we will always perform seek as there are always bunch
        // of information between the genotype dosages
        if (!m_bgen_file.seekg(snp.byte_pos(), std::ios_base::beg)) {
            throw std::runtime_error(
                "Error: Cannot seek within the bgen file!");
        }

        auto&& context = m_context_map[m_cur_file];
        setter.set_stat(snp.stat(), m_homcom_weight, m_het_weight,
                        m_homrar_weight, snp.is_flipped(), not_first);
        // start performing the parsing
        genfile::bgen::read_and_parse_genotype_data_block<PRS_Interpreter>(
            m_bgen_file, context, setter, &m_buffer1, &m_buffer2, false);
        // check if this SNP has some non-missing sample, if not, invalidate it
        if (!setter.is_valid()) snp.invalidate();
        // after reading in this SNP, we no longer need to reset the PRS
        not_first = true;
    }
}


void BinaryGen::hard_code_score(const size_t start_index,
                                const size_t end_bound,
                                const size_t region_index, bool set_zero)
{

    // we need to calculate the size of possible vectors
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    uintptr_t* lbptr;
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uintptr_t ulii = 0;
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
    const uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_sample_ct);

    const bool miss_count = (m_missing_score != MISSING_SCORE::SET_ZERO);
    const bool is_centre = (m_missing_score == MISSING_SCORE::CENTER);
    const bool mean_impute = (m_missing_score == MISSING_SCORE::MEAN_IMPUTE);
    // check if we need to reset the sample's PRS
    bool not_first = !set_zero;
    intptr_t nanal;
    double stat, maf, adj_score, miss_score;

    m_cur_file = "";
    // initialize the data structure for storing the genotype
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);

    for (size_t i_snp = start_index; i_snp < end_bound; ++i_snp) {
        // go through each SNP and skip any that does not fall into our current
        // region
        auto&& cur_snp = m_existed_snps[i_snp];
        if (!cur_snp.in(region_index) || !cur_snp.valid()) continue;

        // read in the genotype using the modified load_and_collapse_incl
        // function. m_target_plink will inform the function wheter there's an
        // intermediate file
        if (load_and_collapse_incl(cur_snp.byte_pos(), cur_snp.file_name(),
                                   genotype.data(), m_target_plink))
        {
            throw std::runtime_error("Error: Cannot read the bed file!");
        }
        // if we haven't got the count from the genotype matrix, we will need to
        // calculate that in theory, we might not need to do the counting as
        // that is already done when we convert the dosages into the binary
        // genotypes (TODO)
        if (!cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct)) {
            genovec_3freq(genotype.data(), m_sample_mask.data(), pheno_nm_ctv2,
                          &missing_ct, &het_ct, &homcom_ct);
            cur_snp.set_counts(homcom_ct, het_ct, homrar_ct, missing_ct);
        }
        nanal = static_cast<intptr_t>(m_sample_ct) - missing_ct;
        if (nanal == 0) {
            // ignore sample if all sample has missing genotyep for this SNP
            cur_snp.invalidate();
            continue;
        }
        homcom_weight = m_homcom_weight;
        het_weight = m_het_weight;
        homrar_weight = m_homrar_weight;

        maf =
            static_cast<double>(het_ct * het_weight + homrar_weight * homrar_ct)
            / static_cast<double>(nanal * 2.0);
        if (cur_snp.is_flipped()) {
            // change the mean to reflect flipping
            maf = 1.0 - maf;
            // swap the weighting
            std::swap(homcom_weight, homrar_weight);
        }
        stat = cur_snp.stat() * 2; // Multiply by ploidy

        // only set these value to the imputed value if we require them
        adj_score = 0;
        if (is_centre) adj_score = stat * maf;
        miss_score = 0;
        if (mean_impute) miss_score = stat * maf;

        lbptr = genotype.data();
        uii = 0;
        ulii = 0;
        // start reading the genotype
        do
        {
            // the genotype block we are reading
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct) {
                // PLINK not sure what it does
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            ujj = 0;
            // ujj represent the sample index in the current block
            while (ujj < BITCT) {
                // ukk is the genotype
                ukk = (ulii >> ujj) & 3;
                // uii+ujj/2 is the sample index
                auto&& sample_prs = m_prs_info[uii + (ujj / 2)];
                // now we will get all genotypes (0, 1, 2, 3)
                switch (ukk)
                {
                default:
                    if (not_first) {
                        // we will just add the information
                        ++sample_prs.num_snp;
                        sample_prs.prs +=
                            homcom_weight * stat * 0.5 - adj_score;
                    }
                    else
                    {
                        // we want to reset the prs
                        sample_prs.num_snp = 1;
                        sample_prs.prs = homcom_weight * stat * 0.5 - adj_score;
                    }
                    break;
                case 1:
                    if (not_first) {
                        // we will just add the information
                        ++sample_prs.num_snp;
                        sample_prs.prs += het_weight * stat * 0.5 - adj_score;
                    }
                    else
                    {
                        // we want to reset the prs
                        sample_prs.num_snp = 1;
                        sample_prs.prs = het_weight * stat * 0.5 - adj_score;
                    }
                    break;
                case 3:
                    if (not_first) {
                        // we will just add the information
                        ++sample_prs.num_snp;
                        sample_prs.prs +=
                            homrar_weight * stat * 0.5 - adj_score;
                    }
                    else
                    {
                        // we want to reset the prs
                        sample_prs.num_snp = 1;
                        sample_prs.prs = homrar_weight * stat * 0.5 - adj_score;
                    }
                    break;
                case 2:
                    // handle missing sample
                    if (not_first) {
                        // not first should only be false for the first SNP.
                        // Again, we might have a faster run time using if case
                        // here due to its simplicity + consistency in the
                        // true/false
                        // add 1 to the number of SNP if we are not setting
                        // missing samples to 0
                        sample_prs.num_snp += miss_count;
                        // add the current genotype weight to the score
                        sample_prs.prs += miss_score;
                    }
                    else
                    {
                        // add 1 to the number of SNP if we are not setting
                        // missing samples to 0
                        sample_prs.num_snp = miss_count;
                        // directly assign the new PRS to the storage
                        sample_prs.prs = miss_score;
                    }
                    break;
                }
                // ulii &= ~((3 * ONELU) << ujj);
                // each sample is represented by two byte, so we increment the
                // index by 2
                ujj += 2;
            }
            // we finish one block of samples.
            uii += BITCT2;
        } while (uii < m_sample_ct);
        // we've finish processing the first SNP no longer need to reset the PRS
        not_first = true;
    }
}


void BinaryGen::read_score(const size_t start_index, const size_t end_bound,
                           const size_t region_index, bool set_zero)
{
    if (m_hard_coded) {
        hard_code_score(start_index, end_bound, region_index, set_zero);
        return;
    }
    else
    {
        dosage_score(start_index, end_bound, region_index, set_zero);
        return;
    }
}

void BinaryGen::read_score(const std::vector<size_t>& index, bool reset_zero)
{
    // because I don't want to touch the code in dosage_score, we will reset
    // the sample here reset_sample_prs();
    if (m_hard_coded) {
        // for hard coded, we need to check if intermediate file is used
        // instead
        hard_code_score(index, reset_zero);
    }
    else
    {
        dosage_score(index, reset_zero);
    }
}


void BinaryGen::hard_code_score(const std::vector<size_t>& index, bool set_zero)
{
    // we need to calculate the size of possible vectors
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    uintptr_t* lbptr;
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uintptr_t ulii = 0;
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
    const uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_sample_ct);

    const bool miss_count = (m_missing_score != MISSING_SCORE::SET_ZERO);
    const bool is_centre = (m_missing_score == MISSING_SCORE::CENTER);
    const bool mean_impute = (m_missing_score == MISSING_SCORE::MEAN_IMPUTE);
    // check if we need to reset the sample's PRS
    bool not_first = !set_zero;
    intptr_t nanal;
    double stat, maf, adj_score, miss_score;

    m_cur_file = "";
    // initialize the data structure for storing the genotype
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);

    for (auto&& i_snp : index) {
        // go through each SNP and skip any that does not fall into our current
        // region
        auto&& cur_snp = m_existed_snps[i_snp];
        if (!cur_snp.valid()) continue;

        // read in the genotype using the modified load_and_collapse_incl
        // function. m_target_plink will inform the function wheter there's an
        // intermediate file
        if (load_and_collapse_incl(cur_snp.byte_pos(), cur_snp.file_name(),
                                   genotype.data(), m_target_plink))
        {
            throw std::runtime_error("Error: Cannot read the bed file!");
        }
        // if we haven't got the count from the genotype matrix, we will need to
        // calculate that in theory, we might not need to do the counting as
        // that is already done when we convert the dosages into the binary
        // genotypes (TODO)
        if (!cur_snp.get_counts(homcom_ct, het_ct, homrar_ct, missing_ct)) {
            genovec_3freq(genotype.data(), m_sample_mask.data(), pheno_nm_ctv2,
                          &missing_ct, &het_ct, &homcom_ct);
            cur_snp.set_counts(homcom_ct, het_ct, homrar_ct, missing_ct);
        }
        nanal = static_cast<intptr_t>(m_sample_ct) - missing_ct;
        if (nanal == 0) {
            // ignore sample if all sample has missing genotyep for this SNP
            cur_snp.invalidate();
            continue;
        }
        homcom_weight = m_homcom_weight;
        het_weight = m_het_weight;
        homrar_weight = m_homrar_weight;

        maf =
            static_cast<double>(het_ct * het_weight + homrar_weight * homrar_ct)
            / static_cast<double>(nanal * 2.0);
        if (cur_snp.is_flipped()) {
            // change the mean to reflect flipping
            maf = 1.0 - maf;
            // swap the weighting
            std::swap(homcom_weight, homrar_weight);
        }
        stat = cur_snp.stat() * 2; // Multiply by ploidy

        // only set these value to the imputed value if we require them
        adj_score = 0;
        if (is_centre) adj_score = stat * maf;
        miss_score = 0;
        if (mean_impute) miss_score = stat * maf;

        lbptr = genotype.data();
        uii = 0;
        ulii = 0;
        // start reading the genotype
        do
        {
            // the genotype block we are reading
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct) {
                // PLINK not sure what it does
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            ujj = 0;
            // ujj represent the sample index in the current block
            while (ujj < BITCT) {
                // ukk is the genotype
                ukk = (ulii >> ujj) & 3;
                // uii+ujj/2 is the sample index
                auto&& sample_prs = m_prs_info[uii + (ujj / 2)];
                // now we will get all genotypes (0, 1, 2, 3)
                switch (ukk)
                {
                default:
                    if (not_first) {
                        // we will just add the information
                        ++sample_prs.num_snp;
                        sample_prs.prs +=
                            homcom_weight * stat * 0.5 - adj_score;
                    }
                    else
                    {
                        // we want to reset the prs
                        sample_prs.num_snp = 1;
                        sample_prs.prs = homcom_weight * stat * 0.5 - adj_score;
                    }
                    break;
                case 1:
                    if (not_first) {
                        // we will just add the information
                        ++sample_prs.num_snp;
                        sample_prs.prs += het_weight * stat * 0.5 - adj_score;
                    }
                    else
                    {
                        // we want to reset the prs
                        sample_prs.num_snp = 1;
                        sample_prs.prs = het_weight * stat * 0.5 - adj_score;
                    }
                    break;
                case 3:
                    if (not_first) {
                        // we will just add the information
                        ++sample_prs.num_snp;
                        sample_prs.prs +=
                            homrar_weight * stat * 0.5 - adj_score;
                    }
                    else
                    {
                        // we want to reset the prs
                        sample_prs.num_snp = 1;
                        sample_prs.prs = homrar_weight * stat * 0.5 - adj_score;
                    }
                    break;
                case 2:
                    // handle missing sample
                    if (not_first) {
                        // not first should only be false for the first SNP.
                        // Again, we might have a faster run time using if case
                        // here due to its simplicity + consistency in the
                        // true/false
                        // add 1 to the number of SNP if we are not setting
                        // missing samples to 0
                        sample_prs.num_snp += miss_count;
                        // add the current genotype weight to the score
                        sample_prs.prs += miss_score;
                    }
                    else
                    {
                        // add 1 to the number of SNP if we are not setting
                        // missing samples to 0
                        sample_prs.num_snp = miss_count;
                        // directly assign the new PRS to the storage
                        sample_prs.prs = miss_score;
                    }
                    break;
                }
                // ulii &= ~((3 * ONELU) << ujj);
                // each sample is represented by two byte, so we increment the
                // index by 2
                ujj += 2;
            }
            // we finish one block of samples.
            uii += BITCT2;
        } while (uii < m_sample_ct);
        // we've finish processing the first SNP no longer need to reset the PRS
        not_first = true;
    }
}


void BinaryGen::dosage_score(const std::vector<size_t>& index, bool set_zero)
{
    m_cur_file = "";
    std::string bgen_name;

    bool not_first = !set_zero;
    // we initialize the PRS interpretor with the required information.
    // m_prs_info is where we store the PRS information
    // and m_sample_include let us know if the sample is required.
    // m_missing_score will inform us as to how to handle the missingness
    PRS_Interpreter setter(&m_prs_info, &m_sample_include, m_missing_score);
    for (auto&& i_snp : index) {
        auto&& snp = m_existed_snps[i_snp];
        // skip SNPs that are not in the region or that are invalid
        if (!snp.valid()) continue;
        if (snp.file_name() == m_cur_file) {
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
        // For bgen file, we will always perform seek as there are always bunch
        // of information between the genotype dosages
        if (!m_bgen_file.seekg(snp.byte_pos(), std::ios_base::beg)) {
            throw std::runtime_error(
                "Error: Cannot seek within the bgen file!");
        }

        auto&& context = m_context_map[m_cur_file];
        setter.set_stat(snp.stat(), m_homcom_weight, m_het_weight,
                        m_homrar_weight, snp.is_flipped(), not_first);
        // start performing the parsing
        genfile::bgen::read_and_parse_genotype_data_block<PRS_Interpreter>(
            m_bgen_file, context, setter, &m_buffer1, &m_buffer2, false);
        // check if this SNP has some non-missing sample, if not, invalidate it
        if (!setter.is_valid()) snp.invalidate();
        // after reading in this SNP, we no longer need to reset the PRS
        not_first = true;
    }
}
