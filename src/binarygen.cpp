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

BinaryGen::BinaryGen(const std::string& prefix, const std::string& sample_file,
                     const size_t thread, const bool ignore_fid,
                     const bool keep_nonfounder, const bool keep_ambig)
    : Genotype(thread, ignore_fid, keep_nonfounder, keep_ambig)
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
    m_sample_file = sample_file;
}

std::vector<Sample> BinaryGen::gen_sample_vector()
{
    bool is_sample_format = check_is_sample_format();
    std::ifstream sample_file(m_sample_file.c_str());
    if (!sample_file.is_open()) {
        std::string error_message =
            "Error: Cannot open sample file: " + m_sample_file;
        throw std::runtime_error(error_message);
    }
    std::string line;
    int sex_col = -1;
    if (is_sample_format) {
        std::getline(sample_file, line);
        std::vector<std::string> header_names = misc::split(line);
        // might want to get the reporter involved.
        fprintf(stderr, "Detected bgen sample file format\n");
        for (size_t i = 3; i < header_names.size(); ++i) {
            std::transform(header_names[i].begin(), header_names[i].end(),
                           header_names[i].begin(), ::toupper);
            if (header_names[i].compare("SEX") == 0) {
                sex_col = i;
                break;
            }
        }
        // read in the second line
        std::getline(sample_file, line);
        if (sex_col != -1) {
            // double check if the format is alright
            std::vector<std::string> header_format = misc::split(line);
            if (header_format[sex_col].compare("D") != 0) {
                std::string error_message =
                    "Error: Sex must be coded as \"D\" in bgen sample file!";
                throw std::runtime_error(error_message);
            }
        }
    }

    int line_id = 0;
    std::vector<Sample> sample_name;
    std::unordered_set<std::string> duplicated_samples;
    std::vector<std::string> duplicated_sample_id;
    std::vector<bool> temp_inclusion_vec;
    bool inclusion = false;
    while (std::getline(sample_file, line)) {
        misc::trim(line);
        line_id++;
        if (!line.empty()) {
            std::vector<std::string> token = misc::split(line);
            // if it is not the sample file, check if this has a header
            // not the best way, but will do it
            if (line_id == 1) {
                std::string header_test = token[0];
                std::transform(header_test.begin(), header_test.end(),
                               header_test.begin(), ::toupper);
                if (header_test.compare("FID") == 0
                    || (header_test.compare("IID") == 0 && m_ignore_fid))
                {
                    // this is the header, skip
                    continue;
                }
            }
            if (token.size()
                < ((sex_col != -1) ? (sex_col) : (1 + !m_ignore_fid)))
            {
                std::string error_message =
                    "Error: Line " + std::to_string(line_id)
                    + " must have at least "
                    + std::to_string((sex_col != -1) ? (sex_col)
                                                     : (1 + !m_ignore_fid))
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
                cur_sample.IID = token[0];
            }
            std::string id = (m_ignore_fid)
                                 ? cur_sample.IID
                                 : cur_sample.FID + "_" + cur_sample.IID;
            cur_sample.pheno = "";
            // true as we assume all bgen samples are founders
            cur_sample.include = true;
            cur_sample.prs = 0.0;
            cur_sample.num_snp = 0.0;
            cur_sample.in_regression = false;
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
                    int sex_info = misc::convert<int>(token[sex_col]);
                    m_num_male += (sex_info == 1);
                    m_num_female += (sex_info == 2);
                    m_num_ambig_sex += (sex_info != 1 && sex_info != 2);
                }
                catch (const std::runtime_error& er)
                {
                    throw std::runtime_error("Error: Invalid sex coding!\n");
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
                sample_name.push_back(cur_sample);
            }
        }
    }
    if (!duplicated_sample_id.empty()) {
        // TODO: Produce a file containing id of all valid samples
        std::string error_message =
            "Error: A total of " + std::to_string(duplicated_sample_id.size())
            + " duplicated samples detected! Please ensure all samples have an "
              "unique identifier";
        throw std::runtime_error(error_message);
    }

    m_founder_ct = m_unfiltered_sample_ct;
    // now set all those vectors
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    // don't bother with founder info here as we don't have this information
    m_sample_include.resize(unfiltered_sample_ctl, 0);
    m_founder_info.resize(unfiltered_sample_ctl, 0);
    m_num_male = 0, m_num_female = 0, m_num_ambig_sex = 0,
    m_num_non_founder = 0;
    for (size_t i = 0; i < temp_inclusion_vec.size(); ++i) {
        if (temp_inclusion_vec[i]) {
            SET_BIT(i, m_sample_include.data());
            // we assume all bgen samples to be founder
            SET_BIT(i, m_founder_info.data());
        }
    }
    sample_file.close();
    return sample_name;
}

bool BinaryGen::check_is_sample_format()
{
    std::ifstream sample_file(m_sample_file.c_str());
    if (!sample_file.is_open()) {
        std::string error_message =
            "Error: Cannot open sample file: " + m_sample_file;
        throw std::runtime_error(error_message);
    }
    std::string first_line, second_line;
    std::getline(sample_file, first_line);
    std::getline(sample_file, second_line);
    sample_file.close();
    std::vector<std::string> first_row = misc::split(first_line);
    std::vector<std::string> second_row = misc::split(second_line);
    if (first_row.size() != second_row.size() || first_row.size() < 3) {
        return false;
    } // first 3 must be 0
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

void BinaryGen::get_context(std::string& prefix)
{
    std::string bgen_name = prefix + ".bgen";
    std::ifstream bgen_file(bgen_name.c_str(), std::ifstream::binary);
    if (!bgen_file.is_open()) {
        std::string error_message = "Error: Cannot open bgen file " + bgen_name;
        throw std::runtime_error(error_message);
    }
    genfile::bgen::Context context;
    uint32_t offset;
    genfile::bgen::read_little_endian_integer(bgen_file, &offset);
    uint32_t header_size = 0, number_of_snp_blocks = 0, number_of_samples = 0,
             flags = 0;
    char magic[4];
    std::size_t fixed_data_size = 20;
    std::vector<char> free_data;
    genfile::bgen::read_little_endian_integer(bgen_file, &header_size);
    assert(header_size >= fixed_data_size);
    genfile::bgen::read_little_endian_integer(bgen_file, &number_of_snp_blocks);
    genfile::bgen::read_little_endian_integer(bgen_file, &number_of_samples);
    bgen_file.read(&magic[0], 4);
    free_data.resize(header_size - fixed_data_size);
    bgen_file.read(&free_data[0], free_data.size());
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
        bytes_read += 8;
        assert(actual_number_of_samples == context.number_of_samples);
        // we don't need to check the sample name when we are doing the LD
        // as we don't store any
        if(!m_is_ref){
        	for (size_t i = 0; i < actual_number_of_samples; ++i) {
				genfile::bgen::read_length_followed_by_data(
					bgen_file, &identifier_size, &identifier);
				if (!bgen_file)
					throw std::runtime_error("Error: Problem reading bgen file!");
				bytes_read += sizeof(identifier_size) + identifier_size;
				// Only need to use IID as BGEN doesn't have the FID information
				if (m_sample_names[i].IID.compare(identifier) != 0) {
					throw std::runtime_error("Error: Sample mismatch "
											 "between bgen and phenotype "
											 "file!");
				}
			}
        }
        assert(bytes_read == sample_block_size);
    }
    return true;
}
std::vector<SNP> BinaryGen::gen_snp_vector(const double geno, const double maf,
                                           const double info_score,
                                           const double hard_threshold,
                                           const bool hard_coded,
                                           const std::string& out_prefix,
                                           Genotype* target)
{
    std::vector<SNP> snp_res;
    std::unordered_set<std::string> duplicated_snps;
    std::vector<int> ref_target_overlap_index;
    // should only apply to SNPs that are not removed due to extract/exclude
    std::unordered_set<std::string> duplicate_check_list;
    m_hard_threshold = hard_threshold;
    m_hard_coded = hard_coded;
    bool chr_sex_error = false;
    bool chr_error = false;
    bool first_bgen_file = true;
    bool user_exclude = false;
    size_t chr_index = 0;
    size_t total_unfiltered_snps = 0;

    for (auto prefix : m_genotype_files) {
        get_context(prefix);
        if (first_bgen_file) {
            first_bgen_file = false;
            std::string bgen_name = prefix + ".bgen";
            check_sample_consistent(bgen_name, m_context_map[prefix]);
        }
        total_unfiltered_snps += m_context_map[prefix].number_of_variants;
    }
    // first pass to get the total SNP number such that we can speed up the push
    // back
    snp_res.reserve(total_unfiltered_snps);
    // to allow multiple file for one chromosome, we put these variable outside
    // the for loop
    std::string prev_chr = "";
    int chr_code = 0;
    for (auto prefix : m_genotype_files) {
        std::string bgen_name = prefix + ".bgen";
        std::ifstream bgen_file(bgen_name.c_str(), std::ifstream::binary);
        if (!bgen_file.is_open()) {
            std::string error_message =
                "Error: Cannot open bgen file " + bgen_name;
            throw std::runtime_error(error_message);
        }
        uint32_t offset = m_context_map[prefix].offset;
        bgen_file.seekg(offset + 4);
        uint32_t num_snp = m_context_map[prefix].number_of_variants;
        auto&& context = m_context_map[prefix];
        for (size_t i_snp = 0; i_snp < num_snp; ++i_snp) {
            if (i_snp % 1000 == 0) {
                fprintf(stderr, "\r%zuK SNPs processed in %s\r", i_snp / 1000,
                        bgen_name.c_str());
            }
            m_unfiltered_marker_ct++;
            std::string allele;
            std::string SNPID;
            std::string RSID;
            std::string chromosome;
            uint32_t SNP_position;
            std::vector<std::string> alleles;
            // directly use the libraryread_snp_identifying_data(

            read_snp_identifying_data(
                bgen_file, context, &SNPID, &RSID, &chromosome, &SNP_position,
                [&alleles](std::size_t n) { alleles.resize(n); },
                [&alleles](std::size_t i, std::string const& allele) {
                    alleles.at(i) = allele;
                });

            for (auto&& a : alleles) {

                std::transform(a.begin(), a.end(), a.begin(), ::toupper);
            }
            std::streampos byte_pos = bgen_file.tellg();
            bool exclude_snp = false;
            // but we will not process anything
            if (chromosome.compare(prev_chr) != 0) {
                prev_chr = chromosome;
                if (m_chr_order.find(chromosome) != m_chr_order.end()) {

                    throw std::runtime_error("Error: SNPs on the same "
                                             "chromosome must be clustered "
                                             "together!");
                }
                m_chr_order[chromosome] = chr_index++;
                chr_code = get_chrom_code_raw(chromosome.c_str());
                if (((const uint32_t) chr_code) > m_max_code)
                { // bigger than the maximum code, ignore it
                    if (!chr_error) {
                        fprintf(stderr,
                                "Warning: SNPs with chromosome number larger "
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
                        fprintf(stderr, "Warning: Currently not support "
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
                    exclude_snp = true;
                }
            }
            else if (target->m_existed_snps_index.find(RSID)
                     == target->m_existed_snps_index.end())
            {
                exclude_snp = true;
            }
            if (duplicate_check_list.find(RSID) != duplicate_check_list.end()) {
                duplicated_snps.insert(RSID);
            }
            else if (!exclude_snp && ambiguous(alleles.front(), alleles.back()))
            {
                m_num_ambig++;
                if (!m_keep_ambig) exclude_snp = true;
            }
            if (!user_exclude) {
                duplicate_check_list.insert(RSID);
            }
            std::vector<genfile::byte_t> buffer1;
            std::vector<genfile::byte_t>* buffer2 = nullptr;
            read_genotype_data_block(bgen_file, context, &buffer1);
            // if we want to exclude this SNP, we will not perform decompression
            if (!exclude_snp) {
                // now filter
                //
                if (!(maf <= 0.0 && geno >= 1.0 && info_score <= 0.0)) {
                    QC_Checker setter(&m_sample_include, hard_threshold,
                                      hard_coded);
                    genfile::bgen::uncompress_probability_data(context, buffer1,
                                                               buffer2);
                    genfile::bgen::parse_probability_data<QC_Checker>(
                        &(*buffer2)[0], &(*buffer2)[0] + buffer2->size(),
                        context, setter);
                    if (setter.filter_snp(maf, geno, info_score)) continue;
                }
                if (!m_is_ref) {
                    m_existed_snps_index[RSID] = snp_res.size();
                    // TODO: Update SNP constructor
                    snp_res.emplace_back(SNP(RSID, chr_code, SNP_position,
                                             alleles.front(), alleles.back(),
                                             prefix, byte_pos));
                }
                else
                {
                    auto&& target_index = target->m_existed_snps_index[RSID];
                    bool dummy;
                    if (!target->m_existed_snps[target_index].matching(
                            chr_code, SNP_position, alleles.front(),
                            alleles.back(), dummy))
                    {
                        m_num_ref_target_mismatch++;
                    }
                    target->m_existed_snps[target_index].add_reference(
                        prefix, byte_pos);
                    ref_target_overlap_index.push_back(target_index);
                }
            }
        }
        bgen_file.close();
        fprintf(stderr, "\n");
    }
    snp_res.shrink_to_fit(); // so that it will be more suitable

    if (m_is_ref) {
        target->update_snps(ref_target_overlap_index);
    }
    if (duplicated_snps.size() != 0) {
        std::ofstream log_file_stream;
        std::string dup_name = out_prefix + ".valid";
        log_file_stream.open(dup_name.c_str());
        if (!log_file_stream.is_open()) {
            std::string error_message = "Error: Cannot open file: " + dup_name;
            throw std::runtime_error(error_message);
        }
        for (auto&& snp : snp_res) {
            if (duplicated_snps.find(snp.rs()) != duplicated_snps.end())
                continue;
            log_file_stream << snp.rs() << "\t" << snp.chr() << "\t"
                            << snp.loc() << "\t" << snp.ref() << "\t"
                            << snp.alt() << std::endl;
        }
        log_file_stream.close();
        std::string error_message =
            "Error: Duplicated SNP ID detected!. Valid SNP ID stored at "
            + dup_name + ". You can avoid this error by using --extract "
            + dup_name;
        throw std::runtime_error(error_message);
    }

    return snp_res;
}

BinaryGen::~BinaryGen()
{
    if (m_bgen_file.is_open()) m_bgen_file.close();
}


void BinaryGen::dosage_score(size_t start_index, size_t end_bound,
                             const size_t region_index)
{
    m_cur_file = "";
    std::vector<genfile::byte_t> buffer1, buffer2;
    for (size_t i_snp = start_index; i_snp < end_bound; ++i_snp) {
        auto&& snp = m_existed_snps[i_snp];
        if (!snp.in(region_index)) continue;
        if (m_cur_file.empty() || snp.file_name().compare(m_cur_file) != 0) {
            if (m_bgen_file.is_open()) m_bgen_file.close();
            std::string bgen_name = snp.file_name() + ".bgen";
            m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
            if (!m_bgen_file.is_open()) {
                std::string error_message =
                    "Error: Cannot open bgen file: " + snp.file_name();
                throw std::runtime_error(error_message);
            }
            m_cur_file = snp.file_name();
        }
        m_bgen_file.seekg(snp.byte_pos(), std::ios_base::beg);

        auto&& context = m_context_map[m_cur_file];

        PRS_Interpreter setter(&m_sample_names, m_model, m_missing_score,
                               &m_sample_include, snp.stat() * 2,
                               snp.is_flipped());
        /*
        PRS_Interpreter setter(&m_sample_names, &g_prs_storage, &g_num_snps,
                               vector_pad, m_model, m_missing_score,
                               snp.stat() * 2,
                               snp.is_flipped()); // Multiple by ploidy
                               */
        // after this, m_sample contain the latest PRS score
        genfile::bgen::read_and_parse_genotype_data_block<PRS_Interpreter>(
            m_bgen_file, context, setter, &buffer1, &buffer2, false);
    }
}

void BinaryGen::hard_code_score(size_t start_index, size_t end_bound,
                                const size_t region_index)
{

    m_cur_file = "";
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uintptr_t ulii = 0;
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    uintptr_t final_mask = get_final_mask(m_sample_ct);

    size_t num_sample_read = m_sample_names.size();
    // index is w.r.t. partition, which contain all the information
    std::vector<uintptr_t> genotype(unfiltered_sample_ctl * 2, 0);
    //    uintptr_t* genotype = new uintptr_t[unfiltered_sample_ctl * 2];
    for (size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
    { // for each SNP
        auto&& cur_snp = m_existed_snps[i_snp];
        if (!cur_snp.in(region_index)) continue;
        if (load_and_collapse_incl(cur_snp.byte_pos(), cur_snp.file_name(),
                                   m_unfiltered_sample_ct, m_sample_ct,
                                   m_sample_include.data(), final_mask, false,
                                   m_tmp_genotype.data(), genotype.data()))
        {
            throw std::runtime_error("Error: Cannot read the bed file!");
        }
        uintptr_t* lbptr = genotype.data();
        uii = 0;
        std::vector<size_t> missing_samples;
        double stat = cur_snp.stat() * 2; // Multiply by ploidy
        bool flipped = cur_snp.is_flipped();
        std::vector<double> genotypes(num_sample_read);

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
                    if (sample_idx < num_sample_read) {
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
        } while (uii < num_sample_read);


        size_t i_missing = 0;

        if (num_sample_read - nmiss == 0) {
            cur_snp.invalidate();
            continue;
        }
        // due to the way the binary code works, the aa will always be 0
        // added there just for fun tbh
        aa = num_sample_read - nmiss - aA - AA;
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
                      / ((double) (num_sample_read - nmiss)
                         * 2.0)); // MAF does not count missing
        double center_score = stat * maf;
        size_t num_miss = missing_samples.size();
        size_t actual_index = 0;
        for (size_t i_sample = 0; i_sample < num_sample_read; ++i_sample) {
            auto&& sample = m_sample_names[i_sample];
            if (i_missing < num_miss
                && actual_index == missing_samples[i_missing])
            {
                if (m_missing_score == MISSING_SCORE::MEAN_IMPUTE)
                    sample.prs += center_score;
                // g_prs_storage[i_sample + vector_pad] += center_score;
                if (m_missing_score != MISSING_SCORE::SET_ZERO)
                    sample.num_snp++;
                // g_num_snps[i_sample + vector_pad]++;
                i_missing++;
            }
            else
            { // not missing sample
                if (m_missing_score == MISSING_SCORE::CENTER) {
                    // if centering, we want to keep missing at 0
                    sample.prs -= center_score;
                    // g_prs_storage[i_sample + vector_pad] -= center_score;
                }

                int g = (flipped) ? fabs(genotypes[actual_index] - 2)
                                  : genotypes[actual_index];
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
                sample.num_snp++;
                // g_prs_storage[i_sample + vector_pad] += g * stat * 0.5;
                // g_num_snps[i_sample + vector_pad]++;
            }
            actual_index++;
        }
    }
}

void BinaryGen::read_score(size_t start_index, size_t end_bound,
                           const size_t region_index)
{
    if (m_hard_coded) {
        hard_code_score(start_index, end_bound, region_index);
        return;
    }
    else
        dosage_score(start_index, end_bound, region_index);
}
