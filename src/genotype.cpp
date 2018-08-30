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

#include "genotype.hpp"


std::vector<std::string> Genotype::set_genotype_files(const std::string& prefix)
{
    std::vector<std::string> genotype_files;
    if (prefix.find("#") != std::string::npos) {
        // auto read will only include the autosomes
        for (size_t chr = 1; chr <= m_autosome_ct; ++chr) {
            std::string name = prefix;
            misc::replace_substring(name, "#", std::to_string(chr));
            genotype_files.push_back(name);
        }
    }
    else
    {
        genotype_files.push_back(prefix);
    }
    return genotype_files;
}

std::vector<std::string>
Genotype::load_genotype_prefix(const std::string& file_name)
{
    std::vector<std::string> genotype_files;
    std::ifstream multi;
    multi.open(file_name.c_str());
    if (!multi.is_open()) {
        throw std::runtime_error(
            std::string("Error: Cannot open file: " + file_name));
    }
    std::string line;
    while (std::getline(multi, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        genotype_files.push_back(line);
    }
    multi.close();
    return genotype_files;
}
void Genotype::init_chr(int num_auto, bool no_x, bool no_y, bool no_xy,
                        bool no_mt)
{
    // this initialize haploid mask as the maximum possible number

    if (num_auto < 0) {
        num_auto = -num_auto;
        m_autosome_ct = static_cast<uint32_t>(num_auto);
        m_xymt_codes[X_OFFSET] = -1;
        m_xymt_codes[Y_OFFSET] = -1;
        m_xymt_codes[XY_OFFSET] = -1;
        m_xymt_codes[MT_OFFSET] = -1;
        m_max_code = static_cast<uint32_t>(num_auto);
        fill_all_bits((static_cast<uint32_t>(num_auto)) + 1,
                      m_haploid_mask.data());
    }
    else
    {
        uint32_t unsign_num_auto = static_cast<uint32_t>(num_auto);
        m_autosome_ct = unsign_num_auto;
        m_xymt_codes[X_OFFSET] = num_auto + 1;
        m_xymt_codes[Y_OFFSET] = num_auto + 2;
        m_xymt_codes[XY_OFFSET] = num_auto + 3;
        m_xymt_codes[MT_OFFSET] = num_auto + 4;
        set_bit(unsign_num_auto + 1, m_haploid_mask.data());
        set_bit(unsign_num_auto + 2, m_haploid_mask.data());
        if (no_x) {
            m_xymt_codes[X_OFFSET] = -1;
            clear_bit(unsign_num_auto + 1, m_haploid_mask.data());
        }
        if (no_y) {
            m_xymt_codes[Y_OFFSET] = -1;
            clear_bit(unsign_num_auto + 2, m_haploid_mask.data());
        }
        if (no_xy) {
            m_xymt_codes[XY_OFFSET] = -1;
        }
        if (no_mt) {
            m_xymt_codes[MT_OFFSET] = -1;
        }
        if (m_xymt_codes[MT_OFFSET] != -1) {
            m_max_code = unsign_num_auto + 4;
        }
        else if (m_xymt_codes[XY_OFFSET] != -1)
        {
            m_max_code = unsign_num_auto + 3;
        }
        else if (m_xymt_codes[Y_OFFSET] != -1)
        {
            m_max_code = unsign_num_auto + 2;
        }
        else if (m_xymt_codes[X_OFFSET] != -1)
        {
            m_max_code = unsign_num_auto + 1;
        }
        else
        {
            m_max_code = unsign_num_auto;
        }
    }
    /*
    fill_all_bits(m_autosome_ct + 1, m_chrom_mask.data());
    for (uint32_t xymt_idx = 0; xymt_idx < XYMT_OFFSET_CT; ++xymt_idx) {
        int32_t cur_code = m_xymt_codes[xymt_idx];
        if (cur_code != -1) {
            set_bit(m_xymt_codes[xymt_idx], m_chrom_mask.data());
        }
    }
    m_chrom_start.resize(m_max_code); // 1 extra for the info
    */
}

std::unordered_set<std::string> Genotype::load_snp_list(std::string input,
                                                        Reporter& reporter)
{
    std::ifstream in;
    // first, we read in the file
    in.open(input.c_str());
    if (!in.is_open()) {
        std::string error_message = "Error: Cannot open file: " + input;
        throw std::runtime_error(error_message);
    }
    std::string line;
    // we will return the "result" variable
    std::unordered_set<std::string> result;
    std::string message;
    // allow the input to be slightly more flexible. This flexibility relies on
    // the header or the column name of the file
    std::getline(in, line);
    misc::trim(line);
    std::vector<std::string> token = misc::split(line);
    // rs_index store the location of the RS ID
    size_t rs_index = 0;
    if (token.size() != 1) {
        // if we have more than one column, we will try and see if it is some
        // format we can recognized
        bool has_snp_colname = false;
        for (auto&& name : token) {
            // we will go through the header and identify any of the followings
            // (case insensitive): 1) SNP 2) RS 3) RS_ID 4) RS.ID 5) RSID

            std::transform(name.begin(), name.end(), name.begin(), ::toupper);
            if (name == "SNP" || name == "RS" || name == "RS_ID"
                || name == "RS.ID" || name == "RSID")
            {
                /// we will assume this column to contain the SNP ID
                has_snp_colname = true;
                message = name + " assume to be column containing SNP ID";
                reporter.report(message);
                break;
            }
            // we will continue to iterate until we reaches the end or found a
            // column with the specific names
            rs_index++;
        }
        if (!has_snp_colname) {
            // if we reaches to the end of the column name list and didn't find
            // the required header, we will resort to familiar file formats
            if (token.size() == 6) {
                // with 6 column, we will assume this to be a bim file, where
                // the SNP ID is at the second column
                message = "SNP extraction/exclusion list contains 6 columns, "
                          "will assume second column contains the SNP ID";
                reporter.report(message);

                // we set the rs index to 1 (use the second column as the RS ID)
                rs_index = 1;
            }
            else
            {
                // otherwise, we will assume the first column contain the SNP ID
                message = "SNP extraction/exclusion list contains "
                          + std::to_string(token.size())
                          + " columns, "
                            "will assume first column contains the SNP ID";
                reporter.report(message);
                // we set the rs index to 0 (use the first column as the RS ID)
                rs_index = 0;
            }
        }
    }
    else
    {
        // with only one column, it is easy, just use that column as the SNP ID
        // rs_index will be 0
        message =
            "Only one column detected, will assume only SNP ID is provided";
        reporter.report(message);
    }
    in.clear();
    in.seekg(0, std::ios::beg);
    // now that we have identify the format of the file we can read in the file
    while (std::getline(in, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        // don't think the parsing . into chr:bp will be helpful in the context
        // of an extraction / exclusion list
        // as this is a set, we don't need to worry about duplicates
        result.insert(token[rs_index]);
    }
    return result;
}

std::unordered_set<std::string> Genotype::load_ref(std::string input,
                                                   bool ignore_fid)
{
    std::ifstream in;
    in.open(input.c_str());
    if (!in.is_open()) {
        std::string error_message = "Error: Cannot open file: " + input;
        throw std::runtime_error(error_message);
    }
    // read in the file
    std::string line;
    // now go through the sample file. We require the FID (if any) and IID  must
    // be the first 1/2 column of the file
    std::unordered_set<std::string> result;
    while (std::getline(in, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (ignore_fid) {
            result.insert(token[0]);
        }
        else
        {
            if (token.size() < 2)
                throw std::runtime_error(
                    "Error: Require FID and IID for extraction. "
                    "You can ignore the FID by using the --ignore-fid flag");
            result.insert(token[0] + "_" + token[1]);
        }
    }
    in.close();
    return result;
}

bool Genotype::chr_code_check(int32_t chr_code, bool& sex_error,
                              bool& chr_error, std::string& error_message)
{
    if (is_set(m_haploid_mask.data(), static_cast<uint32_t>(chr_code))
        || chr_code == m_xymt_codes[X_OFFSET]
        || chr_code == m_xymt_codes[Y_OFFSET])
    {
        // we ignore Sex chromosomes and haploid chromosome

        error_message = "Warning: Currently not support "
                        "haploid chromosome and sex "
                        "chromosomes\n";
        sex_error = true;
        return true;
    }
    // will cause slight mix up here as sex_error will = chr_error now
    if ((static_cast<uint32_t>(chr_code)) > m_max_code) {
        // bigger than the maximum code, ignore it
        error_message = "Warning: SNPs with chromosome number larger "
                        "than "
                        + std::to_string(m_max_code) + "."
                        + " They will be ignored!\n";
        chr_error = true;
        return true;
    }
    return false;
}

void Genotype::load_samples(const std::string& keep_file,
                            const std::string& remove_file, bool verbose,
                            Reporter& reporter)
{
    if (!remove_file.empty()) {
        m_sample_selection_list = load_ref(remove_file, m_ignore_fid);
    }
    else if (!keep_file.empty())
    {
        m_remove_sample = false;
        m_sample_selection_list = load_ref(keep_file, m_ignore_fid);
    }
    if (!m_is_ref) {
        // m_sample_names = gen_sample_vector();
        m_sample_id = gen_sample_vector();
    }
    else
    {
        // don't bother loading up the sample vector as it should
        // never be used for reference panel (except for the founder_info)
        gen_sample_vector();
    }
    std::string message = misc::to_string(m_unfiltered_sample_ct) + " people ("
                          + misc::to_string(m_num_male) + " male(s), "
                          + misc::to_string(m_num_female)
                          + " female(s)) observed\n";
    message.append(misc::to_string(m_founder_ct) + " founder(s) included\n");
    if (verbose) reporter.report(message);
    m_sample_selection_list.clear();
}


void Genotype::load_snps(const Commander& commander, Region& exclusion,
                         bool verbose, Reporter& reporter, Genotype* target)
{
    if (!m_is_ref) {
        if (!commander.extract_file().empty()) {
            m_exclude_snp = false;
            m_snp_selection_list =
                load_snp_list(commander.extract_file(), reporter);
        }
        if (!commander.exclude_file().empty()) {
            m_snp_selection_list =
                load_snp_list(commander.exclude_file(), reporter);
        }
    }
    m_existed_snps = gen_snp_vector(commander, exclusion, target);
    m_marker_ct = m_existed_snps.size();
    std::string message = "";
    if (m_num_ambig != 0 && !m_keep_ambig) {
        message.append(std::to_string(m_num_ambig)
                       + " ambiguous variant(s) excluded\n");
    }
    else if (m_num_ambig != 0)
    {
        message.append(std::to_string(m_num_ambig)
                       + " ambiguous variant(s) kept\n");
    }
    if (m_num_geno_filter != 0) {
        message.append(
            std::to_string(m_num_geno_filter)
            + " variant(s) excluded based on genotype missingness threshold\n");
    }
    if (m_num_maf_filter != 0) {
        message.append(std::to_string(m_num_maf_filter)
                       + " variant(s) excluded based on MAF threshold\n");
    }
    if (m_num_info_filter != 0) {
        message.append(
            std::to_string(m_num_maf_filter)
            + " variant(s) excluded based on INFO score threshold\n");
    }
    if (!m_is_ref) {
        message.append(std::to_string(m_marker_ct) + " variant(s) included\n");
    }
    else
    {
        message.append(std::to_string(target->m_existed_snps.size())
                       + " variant(s) remained\n");
    }
    if (verbose) reporter.report(message);
    m_snp_selection_list.clear();
}

Genotype::~Genotype() {}

void Genotype::read_base(const Commander& c_commander, Region& region,
                         Reporter& reporter)
{
    // can assume region is of the same order as m_existed_snp

    const std::vector<int> index = c_commander.index();
    std::vector<std::vector<std::string>::size_type> casted_index;
    for (auto&& i : index) {
        casted_index.push_back(
            static_cast<std::vector<std::string>::size_type>(i));
    }

    const std::vector<std::string>::size_type max_index =
        casted_index[+BASE_INDEX::MAX];
    const std::string input = c_commander.base_name();
    GZSTREAM_NAMESPACE::igzstream gz_snp_file;
    std::ifstream snp_file;
    std::ofstream mismatch_snp_record;
    double info_threshold;
    const bool info_filter = c_commander.base_info_score(info_threshold);
    double maf_control;
    const bool maf_control_filter = c_commander.maf_base_control(maf_control);
    double maf_case;
    const bool maf_case_filter = c_commander.maf_base_case(maf_case);
    const bool beta = c_commander.beta();
    const bool fastscore = c_commander.fastscore();
    const bool no_full = c_commander.no_full();
    // Read in threshold information
    double max_threshold = (c_commander.fastscore()) ? c_commander.bar_upper()
                                                     : c_commander.upper();
    const double bound_start = c_commander.lower();
    const double bound_end = c_commander.upper();
    const double bound_inter = c_commander.inter();
    if (!c_commander.no_full()) max_threshold = 1.0;
    // Start reading the base file. If the base file contain gz as its suffix,
    // we will read it as a gz file
    bool gz_input = false;
    if (input.substr(input.find_last_of(".") + 1).compare("gz") == 0) {
        gz_snp_file.open(input.c_str());
        if (!gz_snp_file.good()) {
            std::string error_message =
                "Error: Cannot open base file (gz) to read!\n";
            throw std::runtime_error(error_message);
        }
        gz_input = true;
    }
    // otherwise, we will assume it to be normal uncompressed file
    if (!gz_input) {
        snp_file.open(input.c_str());
        if (!snp_file.is_open()) {
            std::string error_message =
                "Error: Cannot open base file: " + input;
            throw std::runtime_error(error_message);
        }
    }

    std::vector<std::string> token;
    std::string line;
    std::string message = "Base file: " + input + "\n";
    std::string mismatch_snp_record_name = c_commander.out() + ".mismatch";
    // Some QC counts
    std::string rs_id;
    std::string ref_allele;
    std::string alt_allele;
    double maf = 1;
    double info_score = 1;
    double pvalue = 2.0;
    double stat = 0.0;
    size_t num_duplicated = 0;
    size_t num_excluded = 0;
    size_t num_ambiguous = 0;
    size_t num_haploid = 0;
    size_t num_not_found = 0;
    size_t num_mismatched = 0;
    size_t num_not_converted = 0; // this is for NA
    size_t num_negative_stat = 0;
    size_t num_line_in_base = 0;
    size_t num_info_filter = 0;
    size_t num_chr_filter = 0;
    size_t num_maf_filter = 0;
    size_t num_retained = 0;
    bool maf_filtered = false;
    bool exclude = false;
    bool flipped = false;
    bool has_chr = false;
    bool has_bp = false;
    int category = -1;
    double pthres = 0.0;
    int32_t chr_code;
    std::unordered_set<std::string> dup_index;

    // Actual reading the file, will do a bunch of QC
    std::streampos file_length = 0;
    if (gz_input) {
        // gzstream does not support seek, so we can't display progress bar
        if (!c_commander.is_index()) {
            std::getline(gz_snp_file, line);
            message.append("GZ file detected. Header of file is:\n");
            message.append(line + "\n\n");
        }
        reporter.report("Due to library restrictions, we cannot display "
                        "progress bar for gz");
    }
    else
    {
        // get the maximum file length, therefore allow us to display the
        // progress
        snp_file.seekg(0, snp_file.end);
        file_length = snp_file.tellg();
        snp_file.clear();
        snp_file.seekg(0, snp_file.beg);
        if (!c_commander.is_index()) std::getline(snp_file, line);
    }
    std::unordered_set<int> unique_thresholds;
    // vector indicate if the SNP is to be retained
    std::vector<bool> retain_snp(m_existed_snps.size(), false);
    double prev_progress = 0.0;
    while ((!gz_input && std::getline(snp_file, line))
           || (gz_input && std::getline(gz_snp_file, line)))
    {
        if (!gz_input) {
            double progress = static_cast<double>(snp_file.tellg())
                              / static_cast<double>(file_length) * 100;
            if (progress - prev_progress > 0.01) {
                fprintf(stderr, "\rReading %03.2f%%", progress);
                prev_progress = progress;
            }
        }
        misc::trim(line);
        if (line.empty()) continue;
        num_line_in_base++;
        exclude = false;
        token = misc::split(line);
        has_chr = false;
        has_bp = false;
        if (token.size() <= max_index) {
            std::string error_message = line;
            error_message.append("\nMore index than column in data");
            throw std::runtime_error(error_message);
        }

        rs_id = token[casted_index[+BASE_INDEX::RS]];
        if (m_existed_snps_index.find(rs_id) != m_existed_snps_index.end()
            && dup_index.find(rs_id) == dup_index.end())
        {
            // if this is not a duplicated SNP and it is in the target file
            dup_index.insert(rs_id);
            auto&& cur_snp = m_existed_snps[m_existed_snps_index[rs_id]];
            chr_code = -1;
            if (index[+BASE_INDEX::CHR] >= 0) {
                chr_code = get_chrom_code_raw(
                    token[casted_index[+BASE_INDEX::CHR]].c_str());
                if (static_cast<uint32_t>(chr_code) > m_max_code) {
                    if (chr_code != -1) {
                        if (chr_code >= MAX_POSSIBLE_CHROM) {
                            // it is ok as chr_code - MAX_POSSIBLE_CHROM >=0
                            chr_code =
                                m_xymt_codes[static_cast<size_t>(chr_code)
                                             - MAX_POSSIBLE_CHROM];
                            // this is the sex chromosomes
                            // we don't need to output the error as they will be
                            // filtered out before by the genotype read anyway
                            exclude = true;
                            num_haploid++;
                        }
                        else
                        {
                            exclude = true;
                            chr_code = -1;
                            num_chr_filter++;
                        }
                    }
                }
                else if (is_set(m_haploid_mask.data(),
                                static_cast<uint32_t>(chr_code))
                         || chr_code == m_xymt_codes[X_OFFSET]
                         || chr_code == m_xymt_codes[Y_OFFSET])
                {
                    exclude = true;
                    num_haploid++;
                }
            }
            has_chr = (chr_code != -1);
            ref_allele = (index[+BASE_INDEX::REF] >= 0)
                             ? token[casted_index[+BASE_INDEX::REF]]
                             : "";
            alt_allele = (index[+BASE_INDEX::ALT] >= 0)
                             ? token[casted_index[+BASE_INDEX::ALT]]
                             : "";
            std::transform(ref_allele.begin(), ref_allele.end(),
                           ref_allele.begin(), ::toupper);
            std::transform(alt_allele.begin(), alt_allele.end(),
                           alt_allele.begin(), ::toupper);
            int loc = -1;
            if (index[+BASE_INDEX::BP] >= 0) {
                // obtain the SNP coordinate
                try
                {
                    loc = misc::convert<int>(
                        token[casted_index[+BASE_INDEX::BP]].c_str());
                    if (loc < 0) {
                        std::string error_message =
                            "Error: " + rs_id + " has negative loci!\n";
                        throw std::runtime_error(error_message);
                    }
                    has_bp = true;
                }
                catch (...)
                {
                    std::string error_message =
                        "Error: Non-numeric loci for " + rs_id + "!\n";
                    throw std::runtime_error(error_message);
                }
            }
            maf = 1;
            // note, this is for reference
            maf_filtered = false;
            if (maf_control_filter && index[+BASE_INDEX::MAF] >= 0) {
                // only read in if we want to perform MAF filtering
                try
                {
                    maf = misc::convert<double>(
                        token[casted_index[+BASE_INDEX::MAF]].c_str());
                }
                catch (...)
                {
                    // exclude because we can't read the MAF, therefore assume
                    // this is problematic
                    num_maf_filter++;
                    exclude = true;
                    maf_filtered = true;
                }
                if (maf < maf_control) {
                    num_maf_filter++;
                    exclude = true;
                    maf_filtered = true;
                }
            }
            if (maf_case_filter && index[+BASE_INDEX::MAF_CASE] >= 0) {
                // only read in the case MAf if we want to perform filtering on
                // it
                try
                {
                    maf = misc::convert<double>(
                        token[casted_index[+BASE_INDEX::MAF_CASE]].c_str());
                }
                catch (...)
                {
                    // again, any problem in parsing the MAF will lead to SNP
                    // being deleted
                    num_maf_filter += !maf_filtered;
                    exclude = true;
                }
                if (maf < maf_case) {
                    num_maf_filter += !maf_filtered;
                    exclude = true;
                }
            }
            info_score = 1;
            if (info_filter && index[+BASE_INDEX::INFO] >= 0) {
                // obtain the INFO score
                try
                {
                    info_score = misc::convert<double>(
                        token[casted_index[+BASE_INDEX::INFO]].c_str());
                }
                catch (...)
                {
                    // if no info score, just assume it doesn't pass the QC
                    num_info_filter++;
                    exclude = true;
                }
                if (info_score < info_threshold) {
                    num_info_filter++;
                    exclude = true;
                }
            }
            flipped = false;
            if (!cur_snp.matching(chr_code, loc, ref_allele, alt_allele,
                                  flipped))
            {
                // Mismatched SNPs
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
                auto&& target_index = m_existed_snps_index[rs_id];
                std::string chr_out =
                    (has_chr) ? misc::to_string(chr_code) : "-";
                std::string loc_out = (has_bp) ? misc::to_string(loc) : "-";
                std::string alt_allele_out =
                    (alt_allele.empty()) ? "-" : alt_allele;
                mismatch_snp_record
                    << "Base\t" << rs_id << "\t"
                    << m_existed_snps[target_index].chr() << "\t" << chr_out
                    << "\t" << m_existed_snps[target_index].loc() << "\t"
                    << loc_out << "\t" << m_existed_snps[target_index].ref()
                    << "\t" << ref_allele << "\t"
                    << m_existed_snps[target_index].alt() << "\t"
                    << alt_allele_out << "\n";
                num_mismatched++;
                exclude = true;
            }
            pvalue = 2.0;
            try
            {
                pvalue =
                    misc::convert<double>(token[casted_index[+BASE_INDEX::P]]);
                if (pvalue < 0.0 || pvalue > 1.0) {
                    std::string error_message =
                        "Error: Invalid p-value for " + rs_id + "!\n";
                    throw std::runtime_error(error_message);
                }
                else if (pvalue > max_threshold)
                {
                    exclude = true;
                    num_excluded++;
                }
            }
            catch (...)
            {
                exclude = true;
                num_not_converted++;
            }
            stat = 0.0;
            try
            {
                stat = misc::convert<double>(
                    token[casted_index[+BASE_INDEX::STAT]]);
                if (stat < 0 && !beta) {
                    num_negative_stat++;
                    exclude = true;
                }
                else if (misc::logically_equal(stat, 0.0) && !beta)
                {
                    // we can't calculate log(0), so we will say we can't
                    // convert it
                    num_not_converted++;
                    exclude = true;
                }
                else if (!beta)
                    // if it is OR, we will perform natural log to convert it to
                    // beta
                    stat = log(stat);
            }
            catch (...)
            {
                // non-numeric statistic
                num_not_converted++;
                exclude = true;
            }

            if (!alt_allele.empty() && ambiguous(ref_allele, alt_allele)) {
                num_ambiguous++;
                exclude = !m_keep_ambig;
            }
            if (!exclude) {
                category = -1;
                pthres = 0.0;
                if (fastscore) {
                    category = c_commander.get_category(pvalue);
                    pthres = c_commander.get_threshold(category);
                }
                else
                {
                    // calculate the threshold instead
                    if (pvalue > bound_end && !no_full) {
                        category = static_cast<int>(std::ceil(
                            (bound_end + 0.1 - bound_start) / bound_inter));
                        pthres = 1.0;
                    }
                    else
                    {
                        category = static_cast<int>(
                            std::ceil((pvalue - bound_start) / bound_inter));
                        category = (category < 0) ? 0 : category;
                        pthres = category * bound_inter + bound_start;
                    }
                }
                if (flipped) cur_snp.set_flipped();
                // ignore the SE as it currently serves no purpose
                // cur_snp.set_retain();
                retain_snp[m_existed_snps_index[rs_id]] = true;
                num_retained++;
                cur_snp.set_statistic(stat, pvalue, category, pthres);
                if (unique_thresholds.find(category) == unique_thresholds.end())
                {
                    unique_thresholds.insert(category);
                    m_thresholds.push_back(pthres);
                    // m_categories.push_back(category);
                }
                // by definition (how we calculate), category cannot be
                // negative, so this comparison should be ok
                if (m_max_category < static_cast<uint32_t>(category)) {
                    m_max_category = static_cast<uint32_t>(category);
                }
            }
        }
        else if (dup_index.find(rs_id) != dup_index.end())
        {
            num_duplicated++;
        }
        else
        {
            num_not_found++;
        }
    }
    if (gz_input)
        gz_snp_file.close();
    else
        snp_file.close();

    fprintf(stderr, "\rReading %03.2f%%\n", 100.0);


    if (num_retained != m_existed_snps.size()) {
        // remove all SNPs that we don't want to retain
        m_existed_snps.erase(
            std::remove_if(m_existed_snps.begin(), m_existed_snps.end(),
                           [&retain_snp, this](const SNP& s) {
                               return !retain_snp[&s - &*begin(m_existed_snps)];
                           }),
            m_existed_snps.end());
        m_existed_snps.shrink_to_fit();
    }

    m_existed_snps_index.clear();
    // now m_existed_snps is ok and can be used directly
    // it is ok here because read_base is the last function that will alter the
    // number of SNPs included in the object before performing clumping
    // vector index is the current index on m_existed_snp
    intptr_t vector_index = 0;
    // we do it here such that the m_existed_snps is sorted correctly
    // low_bound is where the current snp should read from and last_snp is where
    // the last_snp in the vector which doesn't have the up_bound set
    intptr_t low_bound = 0, last_snp = 0;
    intptr_t prev_chr = -1, prev_loc = 0;
    intptr_t diff = 0;
    m_max_window_size = 0;
    // now we iterate thorugh all the SNPs to define the clumping window
    for (auto&& cur_snp : m_existed_snps) {
        if (prev_chr != cur_snp.chr()) {
            // if prev_chr not equal to  current chromosome, we update the
            // previous information to the current SNP (reset)
            prev_chr = cur_snp.chr();
            prev_loc = cur_snp.loc();
            low_bound = vector_index;
        }
        else if (cur_snp.loc() - prev_loc
                 > static_cast<intptr_t>(m_clump_distance))
        {
            // now the chromosome didn't change, and the distance of our current
            // SNP is further away from the previous SNP than our required
            // threshold
            while (cur_snp.loc() - prev_loc
                       > static_cast<intptr_t>(m_clump_distance)
                   && low_bound < vector_index)
            {
                ++low_bound;
                prev_loc =
                    m_existed_snps[static_cast<std::vector<SNP>::size_type>(
                                       low_bound)]
                        .loc();
            }
        }
        // now low_bound should be the first SNP where the core index SNP need
        // to read from
        cur_snp.set_low_bound(low_bound);
        // set the end of the vector as the default up bound
        cur_snp.set_up_bound(static_cast<intptr_t>(m_existed_snps.size()));
        // update all previous SNPs that are out bounud
        while (
            m_existed_snps[static_cast<std::vector<SNP>::size_type>(last_snp)]
                    .chr()
                != cur_snp.chr()
            || cur_snp.loc()
                       - m_existed_snps
                             [static_cast<std::vector<SNP>::size_type>(
                                  last_snp)]
                                 .loc()
                   > static_cast<intptr_t>(m_clump_distance))
        {
            // if the last SNP is on a differenet chromosome or it is to far
            // from the current SNP

            diff = vector_index
                   - m_existed_snps[static_cast<std::vector<SNP>::size_type>(
                                        last_snp)]
                         .low_bound();
            // we will set the up bound of that SNP to the current SNP
            m_existed_snps[static_cast<std::vector<SNP>::size_type>(last_snp)]
                .set_up_bound(vector_index);
            ++last_snp;
            if (m_max_window_size < diff) {
                m_max_window_size = diff;
            }
        }
        // assign the index
        /*
        m_existed_snps_index[cur_snp.rs()] =
            static_cast<std::vector<SNP>::size_type>(vector_index);
            */
        ++vector_index;
        // then we assign the flag for the current SNP
        cur_snp.set_flag(region);
    }
    for (int i = static_cast<int>(m_existed_snps.size()) - 1; i >= 0; i--) {
        auto&& cur_snp =
            m_existed_snps[static_cast<std::vector<SNP>::size_type>(i)];
        if (cur_snp.up_bound() != static_cast<intptr_t>(m_existed_snps.size()))
            break;
        if (m_max_window_size < cur_snp.up_bound() - cur_snp.low_bound()) {
            m_max_window_size = cur_snp.up_bound() - cur_snp.low_bound();
        }
    }
    // Suggest that we want to release memory
    // but this is only a suggest as this is non-binding request
    // Proper way of releasing memory will be to do swarp. Yet that
    // might lead to out of scope or some other error here?
    m_existed_snps.shrink_to_fit();
    m_region_size = static_cast<uint32_t>(region.size());
    message.append(std::to_string(num_line_in_base)
                   + " variant(s) observed in base file, with:\n");
    if (num_duplicated) {
        message.append(std::to_string(num_duplicated)
                       + " duplicated variant(s)\n");
    }
    if (num_excluded) {
        message.append(std::to_string(num_excluded)
                       + " variant(s) excluded due to p-value threshold\n");
    }
    if (num_chr_filter) {
        message.append(
            std::to_string(num_excluded)
            + " variant(s) excluded as they are on unknown/sex chromosome\n");
    }
    if (num_ambiguous) {
        message.append(std::to_string(num_ambiguous) + " ambiguous variant(s)");
        if (!m_keep_ambig) {
            message.append(" excluded");
        }
        message.append("\n");
    }
    if (num_haploid) {
        message.append(std::to_string(num_haploid)
                       + " variant(s) located on haploid chromosome\n");
    }
    if (num_not_found) {
        message.append(std::to_string(num_not_found)
                       + " variant(s) not found in target file\n");
    }
    if (num_mismatched) {
        message.append(std::to_string(num_mismatched)
                       + " mismatched variant(s) excluded\n");
    }
    if (num_not_converted) {
        message.append(std::to_string(num_not_converted)
                       + " NA stat/p-value observed\n");
    }
    if (num_negative_stat) {
        message.append(std::to_string(num_negative_stat)
                       + " negative statistic observed. Maybe you have "
                         "forgotten the --beta flag?\n");
    }
    if (num_info_filter) {
        message.append(std::to_string(num_info_filter)
                       + " variant(s) with INFO score less than "
                       + std::to_string(info_threshold) + "\n");
    }
    if (num_maf_filter) {
        message.append(std::to_string(num_maf_filter)
                       + " variant(s) excluded due to MAF threshold\n");
    }
    message.append(std::to_string(m_existed_snps.size())
                   + " total variant(s) included from base file\n\n");
    if (num_mismatched > 0) {
        message.append(
            "Warning: Mismatched SNPs detected between base and target!");
        message.append("You should check the files are based on the "
                       "same genome build, or that can just be InDels\n");
    }
    reporter.report(message);
    if (m_existed_snps.size() == 0) {
        throw std::runtime_error("Error: No valid variant remaining");
    }
    m_num_threshold = static_cast<uint32_t>(unique_thresholds.size());
}

void Genotype::set_info(const Commander& c_commander)
{
    m_clump_p = c_commander.clump_p();
    m_clump_r2 = c_commander.clump_r2();
    m_use_proxy = c_commander.proxy(m_clump_proxy);
    m_clump_distance = static_cast<uintptr_t>(c_commander.clump_dist());
    m_model = c_commander.model();
    m_missing_score = c_commander.get_missing_score();
    m_scoring = c_commander.get_score();
    m_seed = c_commander.seed();
    switch (m_model)
    {
    case MODEL::HETEROZYGOUS:
        m_homcom_weight = 0;
        m_het_weight = 1;
        m_homrar_weight = 0;
        break;
    case MODEL::DOMINANT:
        m_homcom_weight = 0;
        m_het_weight = 1;
        m_homrar_weight = 1;
        break;
    case MODEL::RECESSIVE:
        m_homcom_weight = 0;
        m_het_weight = 0;
        m_homrar_weight = 1;
        break;
    default: break;
    }
}

void Genotype::pearson_clump(Genotype& reference, Reporter& reporter)
{
    // All sample sizes used here must correspond to the sample size in the
    // reference panel
    // NOTE: most of the following definitions are similar to those used in
    // efficient clumping.
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(reference.m_unfiltered_sample_ct);
    const uintptr_t founder_ctv2 =
        QUATERCT_TO_ALIGNED_WORDCT(reference.m_founder_ct);
    const uintptr_t founder_ctwd = reference.m_founder_ct / BITCT2;
    const uintptr_t founder_ctwd12 = founder_ctwd / 12;
    const uintptr_t founder_ctwd12_rem = founder_ctwd - (12 * founder_ctwd12);
    const uintptr_t lshift_last =
        2 * ((0x7fffffc0 - reference.m_founder_ct) % BITCT2);
    const uintptr_t founder_ct_mld =
        (reference.m_founder_ct + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
    const uint32_t founder_ctv3 =
        BITCT_TO_ALIGNED_WORDCT(static_cast<uint32_t>(reference.m_founder_ct));
    const uint32_t founder_ctsplit = 3 * founder_ctv3;
    const uint32_t founder_ct_mld_m1 =
        (static_cast<uint32_t>(founder_ct_mld)) - 1;
#ifdef __LP64__
    // not sure what the MULTIPLEX_LD does. In the PLINK code it seems to be
    // compile time defined
    const uint32_t founder_ct_mld_rem =
        (MULTIPLEX_LD / 192)
        - (founder_ct_mld * MULTIPLEX_LD - reference.m_founder_ct) / 192;
#else
    const uint32_t founder_ct_mld_rem =
        (MULTIPLEX_LD / 48)
        - (founder_ct_mld * MULTIPLEX_LD - reference.founder_ct()) / 48;
#endif

    const uintptr_t founder_ct_192_long =
        founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2)
        + founder_ct_mld_rem * (192 / BITCT2);
    const double min_r2 =
        (m_use_proxy) ? std::min(m_clump_proxy, m_clump_r2) : m_clump_r2;
    const bool is_x = false;
    // initialize the genotype vector for storage
    std::vector<uintptr_t> genotype_vector(unfiltered_sample_ctl * 2);
    std::vector<uintptr_t> index_data(3 * founder_ctsplit + founder_ctv3);
    std::vector<uint32_t> index_tots(6);
    // pre-allocate the memory without bothering the memory pool stuff
    std::vector<uint32_t> ld_missing_count(
        static_cast<std::vector<uint32_t>::size_type>(m_max_window_size));
    uint32_t* ld_missing_ct_ptr = nullptr;
    std::unordered_set<double> used_thresholds;
    double dxx;
    double dyy;
    double cov12;
    double r2 = -1.0;
    double non_missing_ctd = 0;
    int32_t dp_result[5];
    uint32_t index_missing_ct = 0;
    uint32_t fixed_non_missing_ct = 0;
    uint32_t non_missing_ct = 0;
    std::vector<uintptr_t> founder_include2(founder_ctv2, 0);
    fill_quatervec_55(static_cast<uint32_t>(reference.m_founder_ct),
                      founder_include2.data());
    m_thresholds.clear();

// reference must have sorted

// try and get a workspace
#ifdef __APPLE__
    int32_t mib[2];
    size_t sztmp;
#endif
    unsigned char* bigstack_ua = nullptr; // ua = unaligned
    unsigned char* bigstack_initial_base;
    int64_t llxx;
    intptr_t default_alloc_mb;
    intptr_t malloc_size_mb = 0;
#ifdef __APPLE__
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    llxx = 0;

    sztmp = sizeof(int64_t);
    sysctl(mib, 2, &llxx, &sztmp, nullptr, 0);
    llxx /= 1048576;
#else
#ifdef _WIN32
    MEMORYSTATUSEX memstatus;
    memstatus.dwLength = sizeof(memstatus);
    GlobalMemoryStatusEx(&memstatus);
    llxx = memstatus.ullTotalPhys / 1048576;
#else
    llxx = ((uint64_t) sysconf(_SC_PHYS_PAGES))
           * ((size_t) sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
    if (!llxx) {
        default_alloc_mb = BIGSTACK_DEFAULT_MB;
    }
    else if (llxx < (BIGSTACK_MIN_MB * 2))
    {
        default_alloc_mb = BIGSTACK_MIN_MB;
    }
    else
    {
        default_alloc_mb = llxx / 2;
    }
    if (!malloc_size_mb) {
        malloc_size_mb = default_alloc_mb;
    }
    else if (malloc_size_mb < BIGSTACK_MIN_MB)
    {
        malloc_size_mb = BIGSTACK_MIN_MB;
    }
    std::string message = "";
#ifndef __LP64__
    if (malloc_size_mb > 2047) {
        malloc_size_mb = 2047;
    }
#endif

    size_t total_required_size = 0;
    size_t require_size = founder_ct_192_long * sizeof(intptr_t);
    require_size = round_up_pow2(require_size, CACHELINE);
    total_required_size = require_size * 2;
    require_size = static_cast<uint32_t>(m_max_window_size)
                   * founder_ct_192_long * sizeof(intptr_t);
    require_size = round_up_pow2(require_size, CACHELINE);
    total_required_size += require_size * 2;
    malloc_size_mb =
        total_required_size * sizeof(intptr_t) / (sizeof(char) * 1048576) + 1;
    if (llxx) {
        if (llxx < malloc_size_mb) {
            std::string error_message =
                "Error: Insufficient memory - require at least "
                + misc::to_string(malloc_size_mb)
                + " MB for clumping but only detected " + misc::to_string(llxx)
                + " MB";
            throw std::runtime_error(error_message);
        }
        message = std::to_string(llxx) + " MB RAM detected; reserving "
                  + std::to_string(malloc_size_mb) + " MB for clumping\n";
    }
    else
    {
        message = "Failed to calculate system memory. Attemping to reserve"
                  + std::to_string(malloc_size_mb) + " MB for clumping\n";
    }
    bigstack_ua =
        (unsigned char*) malloc(malloc_size_mb * 1048576 * sizeof(char));
    if (!bigstack_ua) {
        // if failed to allocate the memory, we will just directly error out as
        // we won't have enough memory for the analysis
        throw std::runtime_error("Failed to allocate required memory");
    }
    reporter.report(message);
    // force 64-byte align to make cache line sensitivity work
    bigstack_initial_base =
        (unsigned char*) round_up_pow2((uintptr_t) bigstack_ua, CACHELINE);

    // series of pointer use to handle the stack
    uintptr_t* index_geno = nullptr;
    uintptr_t* index_mask = nullptr;
    uintptr_t* window_data = nullptr;
    uintptr_t* geno_mask = nullptr;
    uintptr_t* geno_mask_ptr = nullptr;

    std::vector<bool> remain_core(m_existed_snps.size(), false);
    uintptr_t num_core_snps = 0;
    // index_geno points to the base of the stack
    index_geno = (uintptr_t*) bigstack_initial_base;
    require_size = round_up_pow2(require_size, CACHELINE);
    // index_mask then occupy the next stack, leaving enough memory for
    // index_geno (if my require_size calculation is wrong, then it will lead to
    // problem as each pointer might overrun and over-write each other
    index_mask = index_geno + require_size;
    // Then the window_data occupy the last bit of data
    window_data = index_mask + require_size;
    require_size = static_cast<uint32_t>(m_max_window_size)
                   * founder_ct_192_long * sizeof(intptr_t);
    require_size = round_up_pow2(require_size, CACHELINE);
    geno_mask = window_data + require_size;

    uintptr_t* window_data_ptr = nullptr;
    unsigned char* g_bigstack_end =
        &(bigstack_initial_base[(malloc_size_mb * 1048576
                                 - (uintptr_t)(bigstack_initial_base
                                               - bigstack_ua))
                                & (~(CACHELINE - ONELU))]);
    // m_max_window_size represent the maximum number of SNPs required for any
    // one window
    // and max_window_size is the number of windows we can handle in one round
    // given the memory that we have

    uintptr_t max_window_size =
        (((uintptr_t) g_bigstack_end) - ((uintptr_t) bigstack_initial_base))
        / (founder_ctv2 * sizeof(intptr_t));
    // we need two times the memory, so the window need to be divided by 2
    max_window_size /= 2;
    g_bigstack_end = nullptr;
    uintptr_t cur_window_size = 0;
    if (!max_window_size) {
        throw std::runtime_error("Error: Not enough memory for clumping!");
    }
    double prev_progress = -1.0;
    const size_t num_snp = m_existed_snps.size();
    // to improve performance, separate out pearson and non-pearson out
    for (size_t i_snp = 0; i_snp < num_snp; ++i_snp) {
        // print progress
        double progress =
            static_cast<double>(i_snp) / static_cast<double>(num_snp) * 100;
        if (progress - prev_progress > 0.01) {
            fprintf(stderr, "\rClumping Progress: %03.2f%%", progress);
            prev_progress = progress;
        }
        // get the index on the target m_existed_snp of the i_snp th most
        // significant SNP
        auto&& cur_snp_index = m_sort_by_p_index[i_snp];
        // then read in the SNP from the target genotype (this)
        auto&& cur_target_snp = m_existed_snps[cur_snp_index];
        if (cur_target_snp.clumped() || cur_target_snp.p_value() > m_clump_p)
            // skip if clumped
            continue;
        size_t start = static_cast<size_t>(cur_target_snp.low_bound());
        size_t end = static_cast<size_t>(cur_target_snp.up_bound());
        window_data_ptr = window_data;
        ld_missing_ct_ptr = ld_missing_count.data();
        geno_mask_ptr = geno_mask;
        cur_window_size = 0;
        // transversing on TARGET
        for (size_t i_pair = start; i_pair < cur_snp_index; i_pair++) {
            // first we read any SNPs the preceed the core SNP in the file (so
            // that we don't need to backtrack)
            auto&& pair_target_snp = m_existed_snps[i_pair];
            if (pair_target_snp.clumped()
                || pair_target_snp.p_value() > m_clump_p)
                // skip if it is clumped or p-value passed the threshold
                continue;
            // reset the genotype pointer content
            fill_ulong_zero(founder_ct_192_long, window_data_ptr);
            if (++cur_window_size == max_window_size) {
                throw std::runtime_error("Error: Out of memory!");
            }
            // read in the genotype
            reference.read_genotype(window_data_ptr,
                                    pair_target_snp.ref_byte_pos(),
                                    pair_target_snp.ref_file_name());
            // now reset the geno_mst and recalculate it
            fill_ulong_zero(founder_ct_192_long, geno_mask_ptr);
            ld_process_load2(window_data_ptr, geno_mask_ptr, ld_missing_ct_ptr,
                             static_cast<uint32_t>(reference.m_founder_ct),
                             is_x, nullptr);
            // not sure what this iterator is for
            ld_missing_ct_ptr++;
            // move to the next memory region
            geno_mask_ptr = &(geno_mask_ptr[founder_ct_192_long]);
            window_data_ptr = &(window_data_ptr[founder_ct_192_long]);
        }
        if (++cur_window_size == max_window_size) {
            throw std::runtime_error("Error: Out of memory!");
        }
        // not sure why we use a different initailization method here
        window_data_ptr[founder_ctv2 - 2] = 0;
        window_data_ptr[founder_ctv2 - 1] = 0;
        // reset the index genotype content
        std::fill(index_data.begin(), index_data.end(), 0);
        fill_ulong_zero(founder_ct_192_long, index_geno);
        // read in the core SNP
        reference.read_genotype(index_geno, cur_target_snp.ref_byte_pos(),
                                cur_target_snp.ref_file_name());
        // reset the mask for the index
        fill_ulong_zero(founder_ct_192_long, index_mask);
        ld_process_load2(index_geno, index_mask, &index_missing_ct,
                         static_cast<uint32_t>(reference.m_founder_ct), is_x,
                         nullptr);
        // now go back to the beginning of the vector and start readin the
        // stored content
        window_data_ptr = window_data;
        geno_mask_ptr = geno_mask;
        size_t cur_index = 0;
        for (size_t i_pair = start; i_pair < cur_snp_index; i_pair++) {
            auto&& pair_target_snp = m_existed_snps[i_pair];
            if (pair_target_snp.clumped()
                || pair_target_snp.p_value() > m_clump_p)
                // again, we skip any SNPs that are clumped
                continue;
            r2 = -1;
            // taking risk here, we don't know if this is what PLINK
            // acutally wants as they have a complete different structure
            // here (for their multi-thread)
            // Most if not all of the following codes were obtained from PLINK's
            // LD calculation and we are not completely certain what each of the
            // function does
            fixed_non_missing_ct = static_cast<uint32_t>(reference.m_founder_ct)
                                   - index_missing_ct;
            non_missing_ct = fixed_non_missing_ct - ld_missing_count[cur_index];
            if (index_missing_ct && ld_missing_count[cur_index]) {
                non_missing_ct += ld_missing_ct_intersect(
                    geno_mask_ptr, index_mask, founder_ctwd12,
                    founder_ctwd12_rem, lshift_last);
            }
            dp_result[0] = static_cast<int32_t>(reference.m_founder_ct);
            dp_result[1] = -static_cast<int32_t>(fixed_non_missing_ct);
            dp_result[2] = static_cast<int32_t>(ld_missing_count[cur_index])
                           - static_cast<int32_t>(reference.m_founder_ct);
            dp_result[3] = dp_result[1];
            dp_result[4] = dp_result[2];
            ld_dot_prod(window_data_ptr, index_geno, geno_mask_ptr, index_mask,
                        dp_result, founder_ct_mld_m1, founder_ct_mld_rem);
            non_missing_ctd =
                static_cast<double>(static_cast<int32_t>(non_missing_ct));
            dxx = dp_result[1];
            dyy = dp_result[2];
            cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
            dxx = (dp_result[3] * non_missing_ctd + dxx * dxx)
                  * (dp_result[4] * non_missing_ctd + dyy * dyy);
            // we have now obtained the R2
            r2 = (cov12 * cov12) / dxx;
            cur_index++;
            // advance the pointers
            geno_mask_ptr = &(geno_mask_ptr[founder_ct_192_long]);
            window_data_ptr = &(window_data_ptr[founder_ct_192_long]);
            if (r2 >= min_r2) {
                // perform clumping if we need to
                cur_target_snp.clump(pair_target_snp, r2, m_use_proxy,
                                     m_clump_proxy);
            }
        }
        // now we start reading SNPs that comes after the core SNP
        for (size_t i_pair = cur_snp_index + 1; i_pair < end; ++i_pair) {
            // we move to the front of the memory for each corresponding vector
            // as we no longer need to keep the information in the memory
            window_data_ptr = window_data;
            geno_mask_ptr = geno_mask;
            ld_missing_ct_ptr = ld_missing_count.data();
            auto&& pair_target_snp = m_existed_snps[i_pair];
            if (pair_target_snp.clumped()
                || pair_target_snp.p_value() > m_clump_p)
                // again Skip SNPs that we don't want
                continue;
            fill_ulong_zero(founder_ct_192_long, window_data_ptr);
            // read in the SNP information
            reference.read_genotype(window_data_ptr,
                                    pair_target_snp.ref_byte_pos(),
                                    pair_target_snp.ref_file_name());
            r2 = -1;

            fill_ulong_zero(founder_ct_192_long, geno_mask_ptr);
            ld_process_load2(window_data_ptr, geno_mask_ptr, ld_missing_ct_ptr,
                             static_cast<uint32_t>(reference.m_founder_ct),
                             is_x, nullptr);
            // taking risk here, we don't know if this is what PLINK
            // acutally wants as they have a complete different structure
            // here (for their multi-thread)
            fixed_non_missing_ct = static_cast<uint32_t>(reference.m_founder_ct)
                                   - index_missing_ct;
            non_missing_ct = fixed_non_missing_ct - ld_missing_count[cur_index];
            if (index_missing_ct && ld_missing_count[cur_index]) {
                non_missing_ct += ld_missing_ct_intersect(
                    geno_mask_ptr, index_mask, founder_ctwd12,
                    founder_ctwd12_rem, lshift_last);
            }
            dp_result[0] = static_cast<int32_t>(reference.m_founder_ct);
            dp_result[1] = -static_cast<int32_t>(fixed_non_missing_ct);
            dp_result[2] = static_cast<int32_t>(ld_missing_count[0])
                           - static_cast<int32_t>(reference.m_founder_ct);
            dp_result[3] = dp_result[1];
            dp_result[4] = dp_result[2];
            ld_dot_prod(window_data_ptr, index_geno, geno_mask_ptr, index_mask,
                        dp_result, founder_ct_mld_m1, founder_ct_mld_rem);
            non_missing_ctd =
                static_cast<double>(static_cast<int32_t>(non_missing_ct));
            dxx = dp_result[1];
            dyy = dp_result[2];
            cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
            dxx = (dp_result[3] * non_missing_ctd + dxx * dxx)
                  * (dp_result[4] * non_missing_ctd + dyy * dyy);
            // obtained the R2
            r2 = (cov12 * cov12) / dxx;
            if (r2 >= min_r2) {
                // clump the SNP if R2 > threshold
                cur_target_snp.clump(pair_target_snp, r2, m_use_proxy,
                                     m_clump_proxy);
            }
        }
        // Set the core SNP to be clumped so that it will not be clumped out by
        // other SNPs
        cur_target_snp.set_clumped();
        // remain_core_snps.push_back(cur_snp_index);
        // indicate we are keeping this SNP
        remain_core[cur_snp_index] = true;
        num_core_snps++;
        double thres = cur_target_snp.get_threshold();
        // store the p-value threshold information of this SNP so that we can
        // know which threshold contain SNPs, therefore use this information for
        // all-score output
        if (used_thresholds.find(thres) == used_thresholds.end()) {
            used_thresholds.insert(thres);
            m_thresholds.push_back(thres);
        }
    }
    fprintf(stderr, "\rClumping Progress: %03.2f%%\n\n", 100.0);
    // release the memory
    window_data = nullptr;
    window_data_ptr = nullptr;
    free(bigstack_ua);
    bigstack_ua = nullptr;
    bigstack_initial_base = nullptr;

    // we now remove the clumped SNP from our data
    if (num_core_snps != m_existed_snps.size()) {
        //  if there are clumped SNPs, we would like to remove them from the
        //  vector
        // the remain_core SNP has the same order as the m_existed_snps
        // (core_snp_index is the index correspond to the m_existed_snp of
        // target) so it is ok for us to simply perform the remove if

        // because of the algorithm, any SNP with p-value above the clump-p
        // threshold will never be process, thus the remain_core for those SNPs
        // will always be 0, and will be removed (though they will also
        // misleadingly be considered as "clumped")
        // TODO: Issue a warning if the highest p-value threshold is higher than
        // the clump-p threshold
        m_existed_snps.erase(
            std::remove_if(
                m_existed_snps.begin(), m_existed_snps.end(),
                [&remain_core, this](const SNP& s) {
                    return !remain_core[&s - &*begin(m_existed_snps)];
                }),
            m_existed_snps.end());
        m_existed_snps.shrink_to_fit();
    }
    // we no longer require the index. might as well clear it (and hope it will
    // release the memory)
    m_existed_snps_index.clear();
    // no longer require the m_existed_snps_index
    message = "";
    message.append("Number of variant(s) after clumping : "
                   + std::to_string(m_existed_snps.size()) + "\n");
    reporter.report(message);
}

void Genotype::efficient_clumping(Genotype& reference, Reporter& reporter,
                                  bool const use_pearson)
{
    // the m_existed_snp must be sorted before coming into this equation
    reporter.report("Start performing clumping");
    if (use_pearson) {
        // use person clump instead
        pearson_clump(reference, reporter);
        return;
    }
    // we want to initialize the vectors with size correspond to the sample in
    // the reference panel. Need to use the unfiltered sample size
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(reference.m_unfiltered_sample_ct);
    // again, we want to initialize the vector containing the founder
    // membership, which require us to know the number of founders in the
    // reference panel
    const uint32_t founder_ctv3 =
        BITCT_TO_ALIGNED_WORDCT(static_cast<uint32_t>(reference.m_founder_ct));
    const uint32_t founder_ctsplit = 3 * founder_ctv3;
    const uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(reference.m_founder_ct);
    const uintptr_t founder_ctv2 =
        QUATERCT_TO_ALIGNED_WORDCT(reference.m_founder_ct);

    // We only want to perform clumping if our R2 is higher than a minimum
    // threshold. Depending on whethre we do proxy clumping or not, the minimum
    // threshold can be the clump_p or clump_proxy parameter
    const double min_r2 =
        (m_use_proxy) ? std::min(m_clump_proxy, m_clump_r2) : m_clump_r2;
    // is_x is used in PLINK to indicate if the genotype is from the X
    // chromsome, as PRSice ignore any sex chromosome, we can set it as a
    // constant false
    const bool is_x = false;

    // when we read in the genotype, the genotype is stored in byte represented
    // by uintptr_t. Then we will use the PLINK 2 functions to obtain the R2
    // As this vector can be big, we will only initialize it once (memory
    // allocation can be expensive)
    std::vector<uintptr_t> genotype_vector(unfiltered_sample_ctl * 2);
    // We need a vector to indicate which SNPs are remaining after clumping so
    // that we can remove any clumped SNPs away. This is achieved by using this
    // boolean vector
    std::vector<bool> remain_core(m_existed_snps.size(), false);
    // next few parameters are used to store the intermediate output from PLINK
    // function, corresponding to the count of different genotypes
    double freq11;
    double freq11_expected;
    double freq1x;
    double freq2x;
    double freqx1;
    double freqx2;
    double dxx;
    // and this is the storage to result R2
    double r2 = -1.0;
    // The following two vectors are used for storing the intermediate output.
    // Again, put memory allocation at the beginning
    std::vector<uintptr_t> index_data(3 * founder_ctsplit + founder_ctv3);
    std::vector<uint32_t> index_tots(6);
    // pre-allocate the memory without bothering the memory pool stuff
    std::vector<uint32_t> ld_missing_count(
        static_cast<std::vector<uint32_t>::size_type>(m_max_window_size));
    // This is a data mask used by PLINK in the calculation of R2. Preallocate
    // to speed up
    std::vector<uintptr_t> founder_include2(founder_ctv2, 0);
    fill_quatervec_55(static_cast<uint32_t>(reference.m_founder_ct),
                      founder_include2.data());
    // information
    std::unordered_set<double> used_thresholds;
    // we clean out the m_thresholds information as this might change after
    // clumping (some threshold might have no more SNP in it)
    m_thresholds.clear();

// one way to speed things up as in PLINK 2 is to pre-allocate the memory space
// for what we need to do next. The following code did precisely that (borrow
// from PLINK2)
#ifdef __APPLE__
    int32_t mib[2];
    size_t sztmp;
#endif
    unsigned char* bigstack_ua = nullptr; // ua = unaligned
    unsigned char* bigstack_initial_base;
    int64_t llxx;
    intptr_t default_alloc_mb;
    intptr_t malloc_size_mb = 0;
#ifdef __APPLE__
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    llxx = 0;

    sztmp = sizeof(int64_t);
    sysctl(mib, 2, &llxx, &sztmp, nullptr, 0);
    llxx /= 1048576;
#else
#ifdef _WIN32
    MEMORYSTATUSEX memstatus;
    memstatus.dwLength = sizeof(memstatus);
    GlobalMemoryStatusEx(&memstatus);
    llxx = memstatus.ullTotalPhys / 1048576;
#else
    llxx = ((uint64_t) sysconf(_SC_PHYS_PAGES))
           * ((size_t) sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
    if (!llxx) {
        default_alloc_mb = BIGSTACK_DEFAULT_MB;
    }
    else if (llxx < (BIGSTACK_MIN_MB * 2))
    {
        default_alloc_mb = BIGSTACK_MIN_MB;
    }
    else
    {
        default_alloc_mb = llxx / 2;
    }
    if (!malloc_size_mb) {
        malloc_size_mb = default_alloc_mb;
    }
    else if (malloc_size_mb < BIGSTACK_MIN_MB)
    {
        malloc_size_mb = BIGSTACK_MIN_MB;
    }
    std::string message = "";
#ifndef __LP64__
    if (malloc_size_mb > 2047) {
        malloc_size_mb = 2047;
    }
#endif

    // m_max_window_size represent the maximum number of SNPs required for any
    // one window

    // calculate the maximum amount of memory required (the thing is, we can't
    // do the analysis if we don't have this amount as the largest region will
    // fail. So should I call this maximum amount of memory required or minimum
    // amount of memory required?
    malloc_size_mb = (static_cast<uintptr_t>(m_max_window_size) + 1)
                         * founder_ctv2 * sizeof(intptr_t) / 1048576
                     + 1;
    if (llxx) {
        // we have detected the memory, but need to check if that's enough
        if (malloc_size_mb > llxx) {
            std::string error_message =
                "Error: Insufficient memory for clumping! Require "
                + misc::to_string(malloc_size_mb) + " MB but detected only "
                + misc::to_string(llxx) + " MB";
            throw std::runtime_error(error_message);
        }
        else
        {
            message = misc::to_string(llxx) + " MB RAM detected; reserving "
                      + misc::to_string(malloc_size_mb) + " MB for clumping\n";
        }
    }
    else
    {
        message = "Failed to calculate system memory. Attemping to reserve"
                  + misc::to_string(malloc_size_mb) + " MB for clumping\n";
    }
    // now allocate the memory into a pointer
    bigstack_ua = reinterpret_cast<unsigned char*>(malloc(
        static_cast<uintptr_t>(malloc_size_mb) * 1048576 * sizeof(char)));
    // if fail, return nullptr which will then get into the while loop
    if (!bigstack_ua) {
        throw std::runtime_error("Failed to allocate required memory");
    }
    else
    {
        message.append("Allocated " + misc::to_string(malloc_size_mb)
                       + " MB successfully\n");
    }
    reporter.report(message);
    // force 64-byte align to make cache line sensitivity work (from PLINK, not
    // familiar with computer programming to know this...)
    // will stay with the old style cast to avoid trouble
    bigstack_initial_base =
        (unsigned char*) round_up_pow2((uintptr_t) bigstack_ua, CACHELINE);

    // window data is the pointer walking through the allocated memory
    uintptr_t* window_data = nullptr;
    // now we can point the window pointer to the start of the allocated memory
    window_data = (uintptr_t*) bigstack_initial_base;
    uintptr_t* window_data_ptr = nullptr;
    unsigned char* g_bigstack_end = &(
        bigstack_initial_base[(static_cast<uintptr_t>(malloc_size_mb) * 1048576
                               - static_cast<uintptr_t>(bigstack_initial_base
                                                        - bigstack_ua))
                              & (~(CACHELINE - ONELU))]);

    uintptr_t num_core_snps = 0;
    // and max_window_size is the number of windows we can handle in one round
    // given the memory that we have. 1 window = 1 SNP
    uintptr_t max_window_size =
        (((uintptr_t) g_bigstack_end) - ((uintptr_t) bigstack_initial_base))
        / (founder_ctv2 * sizeof(intptr_t));

    g_bigstack_end = nullptr;
    // a counter to count how many windows have we read
    uintptr_t cur_window_size = 0;
    // if max_window size is 0, it means we don't have enough memory for
    // clumping
    if (!max_window_size) {
        throw std::runtime_error("Error: Not enough memory for clumping!");
    }

    double prev_progress = -1.0;
    const auto num_snp = m_existed_snps.size();
    for (size_t i_snp = 0; i_snp < num_snp; ++i_snp) {
        // now start iterate through each SNP
        double progress =
            static_cast<double>(i_snp) / static_cast<double>(num_snp) * 100;
        if (progress - prev_progress > 0.01) {
            fprintf(stderr, "\rClumping Progress: %03.2f%%", progress);
            prev_progress = progress;
        }
        // get the index on target.m_existed_snp to the i_th most significant
        // SNP
        auto&& cur_snp_index = m_sort_by_p_index[i_snp];
        // read in the current SNP
        auto&& cur_target_snp = m_existed_snps[cur_snp_index];
        if (cur_target_snp.clumped() || cur_target_snp.p_value() > m_clump_p)
            // ignore any SNP that are clumped or that has a p-value higher than
            // the clump-p threshold
            continue;
        // with the new change, all SNPs are stored in target
        // so we can safely ignore finding it in the reference
        // though we still need to get the sample size information
        // from the reference panel
        // Any SNP with p-value less than clump-p will be ignored
        // because they can never be an index SNP and thus are not of our
        // interested

        // this is the first SNP we should read from
        std::vector<size_t>::size_type start =
            static_cast<std::vector<size_t>::size_type>(
                cur_target_snp.low_bound());
        // this is the first SNP we should ignore
        std::vector<size_t>::size_type end =
            static_cast<std::vector<size_t>::size_type>(
                cur_target_snp.up_bound());
        // reset our pointer to the start of the memory stack as we are working
        // on a new core SNP
        window_data_ptr = window_data;
        cur_window_size = 0;
        // transversing on TARGET
        // now we will transverse any SNP that comes before the index SNP in the
        // file
        for (size_t i_pair = start; i_pair < cur_snp_index; i_pair++) {
            // the start and end correspond to index on m_existed_snps instead
            // of m_sort_by_p, so we can skip reading from m_sort_by_p
            // the current SNP
            auto&& pair_target_snp = m_existed_snps[i_pair];
            if (pair_target_snp.clumped()
                || pair_target_snp.p_value() > m_clump_p)
                // ignore SNP that are clumped or that has higher p-value than
                // threshold
                continue;
            // Something PLINK does. I suspect this is to reset the content of
            // the pointer to 0
            window_data_ptr[founder_ctv2 - 2] = 0;
            window_data_ptr[founder_ctv2 - 1] = 0;
            if (++cur_window_size == max_window_size) {
                throw std::runtime_error("Error: Out of memory!");
            }
            // read in the genotype data from the memory
            // this depends on the type of the reference.
            // Most important information is the ref_byte_pos (reference byte
            // position for reading) and ref_file_name (which reference file
            // should we read from)
            reference.read_genotype(window_data_ptr,
                                    pair_target_snp.ref_byte_pos(),
                                    pair_target_snp.ref_file_name());
            // we then move the pointer forward to the next space in the memory
            window_data_ptr = &(window_data_ptr[founder_ctv2]);
        }

        if (++cur_window_size == max_window_size) {
            throw std::runtime_error("Error: Out of memory!");
        }
        // now we want to read in the index / core SNP
        // reset the content of the pointer again
        window_data_ptr[founder_ctv2 - 2] = 0;
        window_data_ptr[founder_ctv2 - 1] = 0;
        // reset the index_data information
        std::fill(index_data.begin(), index_data.end(), 0);
        // then we can read in the genotype from the reference panel
        // note the use of cur_target_snp
        reference.read_genotype(window_data_ptr, cur_target_snp.ref_byte_pos(),
                                cur_target_snp.ref_file_name());
        // generate the required data mask
        // Disclaimer: For the next few lines, they are from PLINK and I don't
        // fully understand what they are doing
        vec_datamask(reference.m_founder_ct, 0, window_data_ptr,
                     founder_include2.data(), index_data.data());
        // then populate the index_tots
        index_tots[0] = popcount2_longs(index_data.data(), founder_ctl2);
        vec_datamask(reference.m_founder_ct, 2, window_data_ptr,
                     founder_include2.data(), &(index_data[founder_ctv2]));
        index_tots[1] =
            popcount2_longs(&(index_data[founder_ctv2]), founder_ctl2);
        vec_datamask(reference.m_founder_ct, 3, window_data_ptr,
                     founder_include2.data(), &(index_data[2 * founder_ctv2]));
        index_tots[2] =
            popcount2_longs(&(index_data[2 * founder_ctv2]), founder_ctl2);
        // we have finished reading the index and stored the necessary
        // inforamtion, we can now calculate the R2 between the index and
        // previous SNPs
        // move back to the front of the memory (we don't need to read the SNP
        // again as the info is stored in the memory)
        window_data_ptr = window_data;
        for (size_t i_pair = start; i_pair < cur_snp_index; i_pair++) {
            auto&& pair_target_snp = m_existed_snps[i_pair];
            if (pair_target_snp.clumped()
                || pair_target_snp.p_value() > m_clump_p)
                // Again, ignore unwanted SNP
                continue;
            r2 = -1;
            uint32_t counts[18];
            // calculate the counts
            // these counts are then used for calculation of R2. However, I
            // don't fully understand the algorithm here (copy from PLINK2)
            genovec_3freq(window_data_ptr, index_data.data(), founder_ctl2,
                          &(counts[0]), &(counts[1]), &(counts[2]));
            counts[0] = index_tots[0] - counts[0] - counts[1] - counts[2];
            genovec_3freq(window_data_ptr, &(index_data[founder_ctv2]),
                          founder_ctl2, &(counts[3]), &(counts[4]),
                          &(counts[5]));
            counts[3] = index_tots[1] - counts[3] - counts[4] - counts[5];
            genovec_3freq(window_data_ptr, &(index_data[2 * founder_ctv2]),
                          founder_ctl2, &(counts[6]), &(counts[7]),
                          &(counts[8]));
            counts[6] = index_tots[2] - counts[6] - counts[7] - counts[8];
            if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x,
                                        &freqx1, &freqx2, &freq11))
            {
                // if the calculation is sucessful, we can then calculate the R2
                freq11_expected = freqx1 * freq1x;
                dxx = freq11 - freq11_expected;
                // message from PLINK:
                // if r^2 threshold is 0, let everything else through but
                // exclude the apparent zeroes.  Zeroes *are* included if
                // r2_thresh is negative,
                // though (only nans are rejected then).
                if (fabs(dxx) < SMALL_EPSILON
                    || fabs(freq11_expected * freq2x * freqx2) < SMALL_EPSILON)
                {
                    r2 = 0.0;
                }
                else
                {
                    r2 = dxx * dxx / (freq11_expected * freq2x * freqx2);
                }
            }
            if (r2 >= min_r2) {
                // if the R2 between two SNP is higher than the minim threshold,
                // we will perform clumping
                // use the core SNP to clump the pair_target_snp
                cur_target_snp.clump(pair_target_snp, r2, m_use_proxy,
                                     m_clump_proxy);
            }
            // travel to the next snp
            window_data_ptr = &(window_data_ptr[founder_ctv2]);
        }
        // now we can read the SNPs that come after the index SNP in the file

        for (size_t i_pair = cur_snp_index + 1; i_pair < end; ++i_pair) {
            // we don't need to store the SNP information, as we can process
            // each SNP immediately. Always reset the pointer to the beginning
            // of the stack
            window_data_ptr = window_data;
            // read in the SNP information from teh target
            auto&& pair_target_snp = m_existed_snps[i_pair];
            if (pair_target_snp.clumped()
                || pair_target_snp.p_value() > m_clump_p)
                // skip if not required
                continue;
            // reset data
            window_data_ptr[founder_ctv2 - 2] = 0;
            window_data_ptr[founder_ctv2 - 1] = 0;
            // read in the genotype information
            reference.read_genotype(window_data_ptr,
                                    pair_target_snp.ref_byte_pos(),
                                    pair_target_snp.ref_file_name());

            r2 = -1;
            // obtain count using PLINK magic
            uint32_t counts[18];
            genovec_3freq(window_data_ptr, index_data.data(), founder_ctl2,
                          &(counts[0]), &(counts[1]), &(counts[2]));
            counts[0] = index_tots[0] - counts[0] - counts[1] - counts[2];
            genovec_3freq(window_data_ptr, &(index_data[founder_ctv2]),
                          founder_ctl2, &(counts[3]), &(counts[4]),
                          &(counts[5]));
            counts[3] = index_tots[1] - counts[3] - counts[4] - counts[5];
            genovec_3freq(window_data_ptr, &(index_data[2 * founder_ctv2]),
                          founder_ctl2, &(counts[6]), &(counts[7]),
                          &(counts[8]));
            counts[6] = index_tots[2] - counts[6] - counts[7] - counts[8];
            if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x,
                                        &freqx1, &freqx2, &freq11))
            {
                freq11_expected = freqx1 * freq1x;
                dxx = freq11 - freq11_expected;
                // Message in PLINK:
                // if r^2 threshold is 0, let everything else through but
                // exclude the apparent zeroes.  Zeroes *are* included if
                // r2_thresh is negative,
                // though (only nans are rejected then).
                if (fabs(dxx) < SMALL_EPSILON
                    || fabs(freq11_expected * freq2x * freqx2) < SMALL_EPSILON)
                {
                    r2 = 0.0;
                }
                else
                {
                    r2 = dxx * dxx / (freq11_expected * freq2x * freqx2);
                }
            }
            // now perform clumping if required
            if (r2 >= min_r2) {
                cur_target_snp.clump(pair_target_snp, r2, m_use_proxy,
                                     m_clump_proxy);
            }
        }
        // we set the core SNP to be "clumped" so that it will no longer be
        // considered by other SNP
        cur_target_snp.set_clumped();
        // we set the remain_core to true so that we will keep it at the end
        remain_core[cur_snp_index] = true;
        num_core_snps++;
        // we also get the p-value threshold of the SNP so that we know what
        // thresholds are included in our analysis. This information is vital
        // for the generation of all score file
        double thres = cur_target_snp.get_threshold();
        if (used_thresholds.find(thres) == used_thresholds.end()) {
            used_thresholds.insert(thres);
            m_thresholds.push_back(thres);
        }
    }
    fprintf(stderr, "\rClumping Progress: %03.2f%%\n\n", 100.0);
    // now we release the memory stack
    window_data = nullptr;
    window_data_ptr = nullptr;
    free(bigstack_ua);
    bigstack_ua = nullptr;
    bigstack_initial_base = nullptr;
    m_num_threshold = static_cast<uint32_t>(m_thresholds.size());
    if (num_core_snps != m_existed_snps.size()) {
        //  if there are clumped SNPs, we would like to remove them from the
        //  vector
        // the remain_core SNP has the same order as the m_existed_snps
        // (core_snp_index is the index correspond to the m_existed_snp of
        // target) so it is ok for us to simply perform the remove if

        // because of the algorithm, any SNP with p-value above the clump-p
        // threshold will never be process, thus the remain_core for those SNPs
        // will always be 0, and will be removed (though they will also
        // misleadingly be considered as "clumped")
        // TODO: Issue a warning if the highest p-value threshold is higher than
        // the clump-p threshold
        m_existed_snps.erase(
            std::remove_if(
                m_existed_snps.begin(), m_existed_snps.end(),
                [&remain_core, this](const SNP& s) {
                    return !remain_core[&s - &*begin(m_existed_snps)];
                }),
            m_existed_snps.end());
        m_existed_snps.shrink_to_fit();
    }
    // we no longer require the index. might as well clear it (and hope it will
    // release the memory)
    m_existed_snps_index.clear();

    message = "";
    message.append("Number of variant(s) after clumping : "
                   + misc::to_string(m_existed_snps.size()) + "\n");
    reporter.report(message);
}


bool Genotype::prepare_prsice()
{
    if (m_existed_snps.size() == 0) return false;
    // first, sort the threshold
    std::sort(m_thresholds.begin(), m_thresholds.end());
    std::sort(begin(m_existed_snps), end(m_existed_snps),
              [](SNP const& t1, SNP const& t2) {
                  if (t1.category() == t2.category()) {
                      if (t1.file_name().compare(t2.file_name()) == 0) {
                          return t1.byte_pos() < t2.byte_pos();
                      }
                      else
                          return t1.file_name().compare(t2.file_name()) < 0;
                  }
                  else
                      return t1.category() < t2.category();
              });
    return true;
}

void Genotype::get_null_score(const int& set_size, const int& prev_size,
                              const std::vector<size_t>& background_list,
                              const bool first_run,
                              const bool require_statistic)
{

    if (m_existed_snps.empty()
        || static_cast<std::vector<SNP>::size_type>(set_size)
               >= m_existed_snps.size())
        return;
    // we will initailize a selected_snp_index containing the index of SNPs that
    // we'd like to add / assign to our PRS in the current round.
    // we will get anything from (prev_size , set_size]
    std::vector<size_t> selected_snp_index(background_list.begin() + prev_size,
                                           background_list.begin() + set_size);
    // background_list contain the index stored in the m_background_snp_index.
    // This index is related to the order of m_existed_snp. If user used our
    // default input (e.g. bar-levels  =1 and fastscore), this should be the
    // most optimum sorting order. Thus we shouldn't do anything (if user use
    // other threshold, then that will cause an array of other problem)
    // TODO: provide an error if user use p-value thresholding together with
    // PRSet?
    std::sort(selected_snp_index.begin(), selected_snp_index.end());
    read_score(selected_snp_index, first_run);
    if (require_statistic) {
        misc::RunningStat rs;
        size_t num_prs = m_prs_info.size();
        for (size_t i = 0; i < num_prs; ++i) {
            if (!IS_SET(m_sample_include, i)) continue;
            if (m_prs_info[i].num_snp == 0) {
                rs.push(0.0);
            }
            else
            {
                rs.push(m_prs_info[i].prs
                        / static_cast<double>(m_prs_info[i].num_snp));
            }
        }
        m_mean_score = rs.mean();
        m_score_sd = rs.sd();
    }
}


bool Genotype::get_score(int& cur_index, double& cur_threshold,
                         uint32_t& num_snp_included, const size_t region_index,
                         const bool non_cumulate, const bool require_statistic,
                         const bool first_run)
{
    if (m_existed_snps.size() == 0
        || cur_index == static_cast<int>(m_existed_snps.size()))
        return false;
    // we will reset the number of SNP count if we are not running cumulative
    // PRS
    if (non_cumulate) num_snp_included = 0;
    // we will capture the current category and find out all SNPs that fall
    // within this category

    if (cur_index == -1) // first run
    {
        cur_index = 0;
    }
    intptr_t cur_category =
        m_existed_snps[static_cast<std::vector<SNP>::size_type>(cur_index)]
            .category();
    bool ended = false;
    // obtain the current p-value threshold
    cur_threshold =
        m_existed_snps[static_cast<std::vector<SNP>::size_type>(cur_index)]
            .get_threshold();
    // existed snp should be sorted such that the SNPs should be
    // access sequentially

    std::vector<SNP>::size_type end_index =
        static_cast<std::vector<SNP>::size_type>(cur_index);
    for (; end_index < m_existed_snps.size(); ++end_index) {
        // go through each SNP
        if (m_existed_snps[end_index].category() != cur_category) {
            cur_category = m_existed_snps[end_index].category();
            // if we are no longer within the same category, we will end our
            // current iteration
            ended = true;
            break;
        }
        // we want to count the number of SNPs included in the current analysis
        if (m_existed_snps[end_index].in(region_index)) num_snp_included++;
    }
    if (!ended) {
        // if ended isn't set to true, we have reached the end of the
        // m_existed_snp (might want to avoid using ended as the flag name as
        // that is confusing)
        end_index = m_existed_snps.size();
        // we don't update teh cur_category as that will be the same as the one
        // before
    }
    // we now run the read_score function which will construct PRS based on the
    // range specified by cur_index and end_inde for the i th region and will
    // reset the PRS to 0 if  it is none cumulative PRS calculation or if this
    // is the first run
    read_score(static_cast<size_t>(cur_index), end_index, region_index,
               (non_cumulate || first_run));
    // update the current index
    cur_index = static_cast<int>(end_index);
    // if we want to perform standardization, we will need to calculate the
    // mean and SD
    if (require_statistic) {
        misc::RunningStat rs;
        size_t num_prs = m_prs_info.size();
        for (size_t i = 0; i < num_prs; ++i) {
            if (!IS_SET(m_sample_include, i)) continue;
            if (m_prs_info[i].num_snp == 0) {
                rs.push(0.0);
            }
            else
            {
                rs.push(m_prs_info[i].prs
                        / static_cast<double>(m_prs_info[i].num_snp));
            }
        }
        m_mean_score = rs.mean();
        m_score_sd = rs.sd();
    }
    return true;
}


/**
 * DON'T TOUCH AREA
 *
 */

double Genotype::calc_lnlike(double known11, double known12, double known21,
                             double known22, double center_ct_d, double freq11,
                             double freq12, double freq21, double freq22,
                             double half_hethet_share, double freq11_incr)
{
    double lnlike;
    freq11 += freq11_incr;
    freq22 += freq11_incr;
    freq12 += half_hethet_share - freq11_incr;
    freq21 += half_hethet_share - freq11_incr;
    lnlike = center_ct_d * log(freq11 * freq22 + freq12 * freq21);
    if (known11 != 0.0) {
        lnlike += known11 * log(freq11);
    }
    if (known12 != 0.0) {
        lnlike += known12 * log(freq12);
    }
    if (known21 != 0.0) {
        lnlike += known21 * log(freq21);
    }
    if (known22 != 0.0) {
        lnlike += known22 * log(freq22);
    }
    return lnlike;
}

// This is where the magic happens
uint32_t Genotype::em_phase_hethet(double known11, double known12,
                                   double known21, double known22,
                                   uint32_t center_ct, double* freq1x_ptr,
                                   double* freq2x_ptr, double* freqx1_ptr,
                                   double* freqx2_ptr, double* freq11_ptr,
                                   uint32_t* onside_sol_ct_ptr)
{
    // Returns 1 if at least one SNP is monomorphic over all valid
    // observations; returns 0 otherwise, and fills all frequencies
    // using the maximum likelihood solution to the cubic equation.
    // (We're discontinuing most use of EM phasing since better
    // algorithms have been developed, but the two marker case is
    // mathematically clean and fast enough that it'll probably remain
    // useful as an input for some of those better algorithms...)
    double center_ct_d = (int32_t) center_ct;
    double twice_tot = known11 + known12 + known21 + known22 + 2 * center_ct_d;
    uint32_t sol_start_idx = 0;
    uint32_t sol_end_idx = 1;
    double solutions[3];
    double twice_tot_recip;
    double half_hethet_share;
    double freq11;
    double freq12;
    double freq21;
    double freq22;
    double prod_1122;
    double prod_1221;
    double incr_1122;
    double best_sol;
    double best_lnlike;
    double cur_lnlike;
    double freq1x;
    double freq2x;
    double freqx1;
    double freqx2;
    double lbound;
    double dxx;
    uint32_t cur_sol_idx;
    // shouldn't have to worry about subtractive cancellation problems
    // here
    if (twice_tot == 0.0) {
        return 1;
    }
    twice_tot_recip = 1.0 / twice_tot;
    freq11 = known11 * twice_tot_recip;
    freq12 = known12 * twice_tot_recip;
    freq21 = known21 * twice_tot_recip;
    freq22 = known22 * twice_tot_recip;
    prod_1122 = freq11 * freq22;
    prod_1221 = freq12 * freq21;
    half_hethet_share = center_ct_d * twice_tot_recip;
    // the following four values should all be guaranteed nonzero except
    // in the NAN case
    freq1x = freq11 + freq12 + half_hethet_share;
    freq2x = 1.0 - freq1x;
    freqx1 = freq11 + freq21 + half_hethet_share;
    freqx2 = 1.0 - freqx1;
    if (center_ct) {
        if ((prod_1122 != 0.0) || (prod_1221 != 0.0)) {
            sol_end_idx = cubic_real_roots(
                0.5
                    * (freq11 + freq22 - freq12 - freq21
                       - 3 * half_hethet_share),
                0.5
                    * (prod_1122 + prod_1221
                       + half_hethet_share
                             * (freq12 + freq21 - freq11 - freq22
                                + half_hethet_share)),
                -0.5 * half_hethet_share * prod_1122, solutions);
            while (sol_end_idx
                   && (solutions[sol_end_idx - 1]
                       > half_hethet_share + SMALLISH_EPSILON))
            {
                sol_end_idx--;
                // assert(sol_end_idx && sol_end_idx-1 >= 0);
            }
            while ((sol_start_idx < sol_end_idx)
                   && (solutions[sol_start_idx] < -SMALLISH_EPSILON))
            {
                sol_start_idx++;
                // assert((sol_start_idx < sol_end_idx) &&sol_start_idx
                // < 3);
            }
            if (sol_start_idx == sol_end_idx) {
                // Lost a planet Master Obi-Wan has.  How
                // embarrassing... lost root must be a double root at
                // one of the boundary points, just check their
                // likelihoods
                sol_start_idx = 0;
                sol_end_idx = 2;
                solutions[0] = 0;
                solutions[1] = half_hethet_share;
            }
            else
            {
                if (solutions[sol_start_idx] < 0) {
                    solutions[sol_start_idx] = 0;
                } // checking here
                if (solutions[sol_end_idx - 1] > half_hethet_share) {
                    solutions[sol_end_idx - 1] = half_hethet_share;
                }
            }
        }
        else
        {
            solutions[0] = 0;
            if ((freq22 + SMALLISH_EPSILON < half_hethet_share + freq21)
                && (freq21 + SMALLISH_EPSILON < half_hethet_share + freq22))
            {
                sol_end_idx = 3;
                solutions[1] = (half_hethet_share + freq21 - freq22) * 0.5;
                solutions[2] = half_hethet_share;
            }
            else
            {
                sol_end_idx = 2;
                solutions[1] = half_hethet_share;
            }
        }
        best_sol = solutions[sol_start_idx];
        if (sol_end_idx > sol_start_idx + 1) {
            // select largest log likelihood
            best_lnlike = calc_lnlike(known11, known12, known21, known22,
                                      center_ct_d, freq11, freq12, freq21,
                                      freq22, half_hethet_share, best_sol);
            cur_sol_idx = sol_start_idx + 1;
            do
            {
                incr_1122 = solutions[cur_sol_idx];
                cur_lnlike = calc_lnlike(known11, known12, known21, known22,
                                         center_ct_d, freq11, freq12, freq21,
                                         freq22, half_hethet_share, incr_1122);
                if (cur_lnlike > best_lnlike) {
                    cur_lnlike = best_lnlike;
                    best_sol = incr_1122;
                }
            } while (++cur_sol_idx < sol_end_idx);
        }
        if (onside_sol_ct_ptr && (sol_end_idx > sol_start_idx + 1)) {
            if (freqx1 * freq1x >= freq11) {
                dxx = freq1x * freqx1 - freq11;
                if (dxx > half_hethet_share) {
                    dxx = half_hethet_share;
                }
            }
            else
            {
                dxx = 0.0;
            }
            // okay to NOT count suboptimal boundary points because they
            // don't permit direction changes within the main interval
            // this should exactly match haploview_blocks_classify()'s D
            // sign check
            if ((freq11 + best_sol) - freqx1 * freq1x >= 0.0) {
                if (best_sol > dxx + SMALLISH_EPSILON) {
                    lbound = dxx + SMALLISH_EPSILON;
                }
                else
                {
                    lbound = dxx;
                }
                if (best_sol < half_hethet_share - SMALLISH_EPSILON) {
                    half_hethet_share -= SMALLISH_EPSILON;
                }
            }
            else
            {
                if (best_sol > SMALLISH_EPSILON) {
                    lbound = SMALLISH_EPSILON;
                }
                else
                {
                    lbound = 0.0;
                }
                if (best_sol < dxx - SMALLISH_EPSILON) {
                    half_hethet_share = dxx - SMALLISH_EPSILON;
                }
                else
                {
                    half_hethet_share = dxx;
                }
            }
            for (cur_sol_idx = sol_start_idx; cur_sol_idx < sol_end_idx;
                 cur_sol_idx++)
            {
                if (solutions[cur_sol_idx] < lbound) {
                    sol_start_idx++;
                }
                if (solutions[cur_sol_idx] > half_hethet_share) {
                    break;
                }
            }
            if (cur_sol_idx >= sol_start_idx + 2) {
                *onside_sol_ct_ptr = cur_sol_idx - sol_start_idx;
            }
        }
        freq11 += best_sol;
    }
    else if ((prod_1122 == 0.0) && (prod_1221 == 0.0))
    {
        return 1;
    }
    *freq1x_ptr = freq1x;
    *freq2x_ptr = freq2x;
    *freqx1_ptr = freqx1;
    *freqx2_ptr = freqx2;
    *freq11_ptr = freq11;
    return 0;
}

uint32_t Genotype::em_phase_hethet_nobase(uint32_t* counts, uint32_t is_x1,
                                          uint32_t is_x2, double* freq1x_ptr,
                                          double* freq2x_ptr,
                                          double* freqx1_ptr,
                                          double* freqx2_ptr,
                                          double* freq11_ptr)
{
    // if is_x1 and/or is_x2 is set, counts[9]..[17] are male-only
    // counts.
    double known11 = (double) (2 * counts[0] + counts[1] + counts[3]);
    double known12 = (double) (2 * counts[2] + counts[1] + counts[5]);
    double known21 = (double) (2 * counts[6] + counts[3] + counts[7]);
    double known22 = (double) (2 * counts[8] + counts[5] + counts[7]);
    if (is_x1 || is_x2) {
        if (is_x1 && is_x2) {
            known11 -= (double) ((int32_t) counts[9]);
            known12 -= (double) ((int32_t) counts[11]);
            known21 -= (double) ((int32_t) counts[15]);
            known22 -= (double) ((int32_t) counts[17]);
        }
        else if (is_x1)
        {
            known11 -=
                ((double) (2 * counts[9] + counts[10])) * (1.0 - SQRT_HALF);
            known12 -=
                ((double) (2 * counts[11] + counts[10])) * (1.0 - SQRT_HALF);
            known21 -=
                ((double) (2 * counts[15] + counts[16])) * (1.0 - SQRT_HALF);
            known22 -=
                ((double) (2 * counts[17] + counts[16])) * (1.0 - SQRT_HALF);
        }
        else
        {
            known11 -=
                ((double) (2 * counts[9] + counts[12])) * (1.0 - SQRT_HALF);
            known12 -=
                ((double) (2 * counts[11] + counts[12])) * (1.0 - SQRT_HALF);
            known21 -=
                ((double) (2 * counts[15] + counts[14])) * (1.0 - SQRT_HALF);
            known22 -=
                ((double) (2 * counts[17] + counts[14])) * (1.0 - SQRT_HALF);
        }
    }
    return em_phase_hethet(known11, known12, known21, known22, counts[4],
                           freq1x_ptr, freq2x_ptr, freqx1_ptr, freqx2_ptr,
                           freq11_ptr, nullptr);
}

uint32_t Genotype::load_and_split3(uintptr_t* rawbuf,
                                   uint32_t unfiltered_sample_ct,
                                   uintptr_t* casebuf, uint32_t case_ctv,
                                   uint32_t ctrl_ctv, uint32_t do_reverse,
                                   uint32_t is_case_only,
                                   uintptr_t* nm_info_ptr)
{
    uintptr_t* rawbuf_end = &(rawbuf[unfiltered_sample_ct / BITCT2]);
    uintptr_t* ctrlbuf = &(casebuf[3 * case_ctv]);
    uintptr_t case_words[4];
    uintptr_t ctrl_words[4];
    uint32_t case_rem = 0;
    uint32_t ctrl_rem = 0;
    uint32_t read_shift_max = BITCT2;
    uint32_t sample_uidx = 0;
    uint32_t offset0_case = do_reverse * 2 * case_ctv;
    uint32_t offset2_case = (1 - do_reverse) * 2 * case_ctv;
    uint32_t offset0_ctrl = do_reverse * 2 * ctrl_ctv;
    uint32_t offset2_ctrl = (1 - do_reverse) * 2 * ctrl_ctv;
    uint32_t read_shift;
    uintptr_t read_word;
    uintptr_t ulii;

    case_words[0] = 0;
    case_words[1] = 0;
    case_words[2] = 0;
    case_words[3] = 0;
    ctrl_words[0] = 0;
    ctrl_words[1] = 0;
    ctrl_words[2] = 0;
    ctrl_words[3] = 0;
    while (1) {
        while (rawbuf < rawbuf_end) {
            read_word = *rawbuf++;
            for (read_shift = 0; read_shift < read_shift_max;
                 sample_uidx++, read_shift++)
            {
                ulii = read_word & 3; // Both is_set is always true, because
                                      // dummy_nm is set
                case_words[ulii] |= ONELU << case_rem;
                if (++case_rem == BITCT) {
                    casebuf[offset0_case] = case_words[0];
                    casebuf[case_ctv] = case_words[2];
                    casebuf[offset2_case] = case_words[3];
                    casebuf++;
                    case_words[0] = 0;
                    case_words[2] = 0;
                    case_words[3] = 0;
                    case_rem = 0;
                }
                read_word >>= 2;
            }
        }
        if (sample_uidx == unfiltered_sample_ct) {
            if (case_rem) {
                casebuf[offset0_case] = case_words[0];
                casebuf[case_ctv] = case_words[2];
                casebuf[offset2_case] = case_words[3];
            }
            if (ctrl_rem) {
                ctrlbuf[offset0_ctrl] = ctrl_words[0];
                ctrlbuf[ctrl_ctv] = ctrl_words[2];
                ctrlbuf[offset2_ctrl] = ctrl_words[3];
            }
            ulii = 3;
            if (case_words[1]) {
                ulii -= 1;
            }
            if (ctrl_words[1]) {
                ulii -= 2;
            }
            *nm_info_ptr = ulii;
            return 0;
        }
        rawbuf_end++;
        read_shift_max = unfiltered_sample_ct % BITCT2;
    }
}

void Genotype::two_locus_count_table(uintptr_t* lptr1, uintptr_t* lptr2,
                                     uint32_t* counts_3x3, uint32_t sample_ctv3,
                                     uint32_t is_zmiss2)
{
#ifdef __LP64__
    uint32_t uii;
    fill_uint_zero(9, counts_3x3);
    if (!is_zmiss2) {
        two_locus_3x3_tablev((__m128i*) lptr1, (__m128i*) lptr2, counts_3x3,
                             sample_ctv3 / 2, 3);
    }
    else
    {
        two_locus_3x3_tablev((__m128i*) lptr2, (__m128i*) lptr1, counts_3x3,
                             sample_ctv3 / 2, 2);
        uii = counts_3x3[1];
        counts_3x3[1] = counts_3x3[3];
        counts_3x3[3] = uii;
        counts_3x3[6] = counts_3x3[2];
        counts_3x3[7] = counts_3x3[5];
    }
#else
    counts_3x3[0] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
    counts_3x3[3] =
        popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
    counts_3x3[6] =
        popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
    lptr2 = &(lptr2[sample_ctv3]);
    counts_3x3[1] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
    counts_3x3[4] =
        popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
    counts_3x3[7] =
        popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
    if (!is_zmiss2) {
        lptr2 = &(lptr2[sample_ctv3]);
        counts_3x3[2] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
        counts_3x3[5] =
            popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
        counts_3x3[8] = popcount_longs_intersect(
            lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
    }
#endif
}

void Genotype::two_locus_count_table_zmiss1(uintptr_t* lptr1, uintptr_t* lptr2,
                                            uint32_t* counts_3x3,
                                            uint32_t sample_ctv3,
                                            uint32_t is_zmiss2)
{

#ifdef __LP64__
    fill_uint_zero(6, counts_3x3);
    if (is_zmiss2) {
        two_locus_3x3_zmiss_tablev((__m128i*) lptr1, (__m128i*) lptr2,
                                   counts_3x3, sample_ctv3 / 2);
    }
    else
    {
        two_locus_3x3_tablev((__m128i*) lptr1, (__m128i*) lptr2, counts_3x3,
                             sample_ctv3 / 2, 2);
    }
#else
    counts_3x3[0] = popcount_longs_intersect(lptr1, lptr2, sample_ctv3);
    counts_3x3[1] =
        popcount_longs_intersect(lptr1, &(lptr2[sample_ctv3]), sample_ctv3);
    if (!is_zmiss2) {
        counts_3x3[2] = popcount_longs_intersect(
            lptr1, &(lptr2[2 * sample_ctv3]), sample_ctv3);
        counts_3x3[5] = popcount_longs_intersect(
            &(lptr1[sample_ctv3]), &(lptr2[2 * sample_ctv3]), sample_ctv3);
    }
    lptr1 = &(lptr1[sample_ctv3]);
    counts_3x3[3] = popcount_longs_intersect(lptr1, lptr2, sample_ctv3);
    counts_3x3[4] =
        popcount_longs_intersect(lptr1, &(lptr2[sample_ctv3]), sample_ctv3);
#endif
}

#ifdef __LP64__
void Genotype::two_locus_3x3_tablev(__m128i* vec1, __m128i* vec2,
                                    uint32_t* counts_3x3, uint32_t sample_ctv6,
                                    uint32_t iter_ct)
{
    const __m128i m1 = {FIVEMASK, FIVEMASK};
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    __m128i* vec20;
    __m128i* vec21;
    __m128i* vec22;
    __m128i* vend1;
    __m128i loader1;
    __m128i loader20;
    __m128i loader21;
    __m128i loader22;
    __m128i count10;
    __m128i count11;
    __m128i count12;
    __m128i count20;
    __m128i count21;
    __m128i count22;
    __univec acc0;
    __univec acc1;
    __univec acc2;
    uint32_t ct;
    uint32_t ct2;
    while (iter_ct--) {
        ct = sample_ctv6;
        vec20 = vec2;
        vec21 = &(vec20[sample_ctv6]);
        vec22 = &(vec20[2 * sample_ctv6]);
        while (ct >= 30) {
            ct -= 30;
            vend1 = &(vec1[30]);
            acc0.vi = _mm_setzero_si128();
            acc1.vi = _mm_setzero_si128();
            acc2.vi = _mm_setzero_si128();
            do
            {
            two_locus_3x3_tablev_outer:
                loader1 = *vec1++;
                loader20 = *vec20++;
                loader21 = *vec21++;
                loader22 = *vec22++;
                count10 = _mm_and_si128(loader1, loader20);
                count11 = _mm_and_si128(loader1, loader21);
                count12 = _mm_and_si128(loader1, loader22);
                count10 = _mm_sub_epi64(
                    count10, _mm_and_si128(_mm_srli_epi64(count10, 1), m1));
                count11 = _mm_sub_epi64(
                    count11, _mm_and_si128(_mm_srli_epi64(count11, 1), m1));
                count12 = _mm_sub_epi64(
                    count12, _mm_and_si128(_mm_srli_epi64(count12, 1), m1));
            two_locus_3x3_tablev_two_left:
                // unlike the zmiss variant, this apparently does not
                // suffer from enough register spill to justify
                // shrinking the inner loop
                loader1 = *vec1++;
                loader20 = *vec20++;
                loader21 = *vec21++;
                loader22 = *vec22++;
                count20 = _mm_and_si128(loader1, loader20);
                count21 = _mm_and_si128(loader1, loader21);
                count22 = _mm_and_si128(loader1, loader22);
                count20 = _mm_sub_epi64(
                    count20, _mm_and_si128(_mm_srli_epi64(count20, 1), m1));
                count21 = _mm_sub_epi64(
                    count21, _mm_and_si128(_mm_srli_epi64(count21, 1), m1));
                count22 = _mm_sub_epi64(
                    count22, _mm_and_si128(_mm_srli_epi64(count22, 1), m1));
            two_locus_3x3_tablev_one_left:
                loader1 = *vec1++;
                loader20 = *vec20++;
                loader21 = _mm_and_si128(loader1, loader20); // half1
                loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1),
                                         m1); // half2
                count10 = _mm_add_epi64(count10, _mm_and_si128(loader21, m1));
                count20 = _mm_add_epi64(count20, loader22);
                loader20 = *vec21++;
                loader21 = _mm_and_si128(loader1, loader20);
                loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
                count11 = _mm_add_epi64(count11, _mm_and_si128(loader21, m1));
                count21 = _mm_add_epi64(count21, loader22);
                loader20 = *vec22++;
                loader21 = _mm_and_si128(loader1, loader20);
                loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
                count12 = _mm_add_epi64(count12, _mm_and_si128(loader21, m1));
                count22 = _mm_add_epi64(count22, loader22);

                count10 = _mm_add_epi64(
                    _mm_and_si128(count10, m2),
                    _mm_and_si128(_mm_srli_epi64(count10, 2), m2));
                count11 = _mm_add_epi64(
                    _mm_and_si128(count11, m2),
                    _mm_and_si128(_mm_srli_epi64(count11, 2), m2));
                count12 = _mm_add_epi64(
                    _mm_and_si128(count12, m2),
                    _mm_and_si128(_mm_srli_epi64(count12, 2), m2));
                count10 = _mm_add_epi64(
                    count10,
                    _mm_add_epi64(
                        _mm_and_si128(count20, m2),
                        _mm_and_si128(_mm_srli_epi64(count20, 2), m2)));
                count11 = _mm_add_epi64(
                    count11,
                    _mm_add_epi64(
                        _mm_and_si128(count21, m2),
                        _mm_and_si128(_mm_srli_epi64(count21, 2), m2)));
                count12 = _mm_add_epi64(
                    count12,
                    _mm_add_epi64(
                        _mm_and_si128(count22, m2),
                        _mm_and_si128(_mm_srli_epi64(count22, 2), m2)));
                acc0.vi = _mm_add_epi64(
                    acc0.vi,
                    _mm_add_epi64(
                        _mm_and_si128(count10, m4),
                        _mm_and_si128(_mm_srli_epi64(count10, 4), m4)));
                acc1.vi = _mm_add_epi64(
                    acc1.vi,
                    _mm_add_epi64(
                        _mm_and_si128(count11, m4),
                        _mm_and_si128(_mm_srli_epi64(count11, 4), m4)));
                acc2.vi = _mm_add_epi64(
                    acc2.vi,
                    _mm_add_epi64(
                        _mm_and_si128(count12, m4),
                        _mm_and_si128(_mm_srli_epi64(count12, 4), m4)));
            } while (vec1 < vend1);
            const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
            acc0.vi =
                _mm_add_epi64(_mm_and_si128(acc0.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc0.vi, 8), m8));
            acc1.vi =
                _mm_add_epi64(_mm_and_si128(acc1.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
            acc2.vi =
                _mm_add_epi64(_mm_and_si128(acc2.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
            counts_3x3[0] +=
                ((acc0.u8[0] + acc0.u8[1]) * 0x1000100010001LLU) >> 48;
            counts_3x3[1] +=
                ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
            counts_3x3[2] +=
                ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
        }
        if (ct) {
            vend1 = &(vec1[ct]);
            ct2 = ct % 3;
            acc0.vi = _mm_setzero_si128();
            acc1.vi = _mm_setzero_si128();
            acc2.vi = _mm_setzero_si128();
            ct = 0;
            if (ct2) {
                count10 = _mm_setzero_si128();
                count11 = _mm_setzero_si128();
                count12 = _mm_setzero_si128();
                if (ct2 == 2) {
                    goto two_locus_3x3_tablev_two_left;
                }
                count20 = _mm_setzero_si128();
                count21 = _mm_setzero_si128();
                count22 = _mm_setzero_si128();
                goto two_locus_3x3_tablev_one_left;
            }
            goto two_locus_3x3_tablev_outer;
        }
        counts_3x3 = &(counts_3x3[3]);
    }
}
#endif
