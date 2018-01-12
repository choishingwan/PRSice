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
        // auto read will only include the autosomes unless otherwise?
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


void Genotype::init_chr(int num_auto, bool no_x, bool no_y, bool no_xy,
                        bool no_mt)
{
    // this initialize haploid mask as the maximum possible number

    if (num_auto < 0) {
        num_auto = -num_auto;
        m_autosome_ct = num_auto;
        m_xymt_codes[X_OFFSET] = -1;
        m_xymt_codes[Y_OFFSET] = -1;
        m_xymt_codes[XY_OFFSET] = -1;
        m_xymt_codes[MT_OFFSET] = -1;
        m_max_code = num_auto;
        fill_all_bits(((uint32_t) num_auto) + 1, m_haploid_mask.data());
    }
    else
    {
        m_autosome_ct = num_auto;
        m_xymt_codes[X_OFFSET] = num_auto + 1;
        m_xymt_codes[Y_OFFSET] = num_auto + 2;
        m_xymt_codes[XY_OFFSET] = num_auto + 3;
        m_xymt_codes[MT_OFFSET] = num_auto + 4;
        set_bit(num_auto + 1, m_haploid_mask.data());
        set_bit(num_auto + 2, m_haploid_mask.data());
        if (no_x) {
            m_xymt_codes[X_OFFSET] = -1;
            clear_bit(num_auto + 1, m_haploid_mask.data());
        }
        if (no_y) {
            m_xymt_codes[Y_OFFSET] = -1;
            clear_bit(num_auto + 2, m_haploid_mask.data());
        }
        if (no_xy) {
            m_xymt_codes[XY_OFFSET] = -1;
        }
        if (no_mt) {
            m_xymt_codes[MT_OFFSET] = -1;
        }
        if (m_xymt_codes[MT_OFFSET] != -1) {
            m_max_code = num_auto + 4;
        }
        else if (m_xymt_codes[XY_OFFSET] != -1)
        {
            m_max_code = num_auto + 3;
        }
        else if (m_xymt_codes[Y_OFFSET] != -1)
        {
            m_max_code = num_auto + 2;
        }
        else if (m_xymt_codes[X_OFFSET] != -1)
        {
            m_max_code = num_auto + 1;
        }
        else
        {
            m_max_code = num_auto;
        }
    }
    fill_all_bits(m_autosome_ct + 1, m_chrom_mask.data());
    for (uint32_t xymt_idx = 0; xymt_idx < XYMT_OFFSET_CT; ++xymt_idx) {
        int32_t cur_code = m_xymt_codes[xymt_idx];
        if (cur_code != -1) {
            set_bit(m_xymt_codes[xymt_idx], m_chrom_mask.data());
        }
    }
    m_chrom_start.resize(m_max_code); // 1 extra for the info
}

std::unordered_set<std::string> Genotype::load_snp_list(std::string input,
                                                        Reporter& reporter)
{
    std::ifstream in;
    in.open(input.c_str());
    if (!in.is_open()) {
        std::string error_message = "ERROR: Cannot open file: " + input;
        throw std::runtime_error(error_message);
    }
    std::string line;
    std::unordered_set<std::string> result;
    bool error = false;
    while (std::getline(in, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (token[0].compare(".") == 0) {
            if (!error) {
                error = true;
                std::string message =
                    "WARNING: Some SNPs from the "
                    "extraction/exclusion list has rs-id of . "
                    "They will be excluded unless the file contains at least 3 "
                    "columns."
                    "When 3 columns is provided, we will assume the "
                    "second and third columns are the chromosome and "
                    "coordinates "
                    "respectively and will generate an rsid as chr:loc\n";
                reporter.report(message);
            }
            if (token.size() >= 3) {
                token[0] = token[1] + ":" + token[2];
            }
        }
        if (result.find(token[0]) == result.end()) {
            result.insert(token[0]);
        }
    }
    return result;
}

std::unordered_set<std::string> Genotype::load_ref(std::string input,
                                                   bool ignore_fid)
{
    std::ifstream in;
    in.open(input.c_str());
    if (!in.is_open()) {
        std::string error_message = "ERROR: Cannot open file: " + input;
        throw std::runtime_error(error_message);
    }
    std::string line;
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
                    "ERROR: Require FID and IID for extraction. "
                    "You can ignore the FID by using the --ignore-fid flag");
            result.insert(token[0] + "_" + token[1]);
        }
    }
    in.close();
    return result;
}

void Genotype::load_samples(const std::string& keep_file,
                            const std::string& remove_file, bool verbose,
                            Reporter& reporter)
{
    if (!remove_file.empty()) {
        m_sample_selection_list = load_ref(remove_file, m_ignore_fid);
    }
    if (!keep_file.empty()) {
        m_remove_sample = false;
        m_sample_selection_list = load_ref(keep_file, m_ignore_fid);
    }
    m_sample_names = gen_sample_vector();
    std::string message = std::to_string(m_unfiltered_sample_ct) + " people ("
                          + std::to_string(m_num_male) + " male(s), "
                          + std::to_string(m_num_female)
                          + " female(s)) observed\n";
    message.append(std::to_string(m_founder_ct) + " founder(s) included\n");
    if (verbose) reporter.report(message);
    m_sample_selection_list.clear();
}


void Genotype::load_snps(
    const std::string out_prefix,
    const std::unordered_map<std::string, size_t>& existed_snps,
    const double geno, const double maf, const double info,
    const double hard_threshold, const bool hard_coded, bool verbose,
    Reporter& reporter)
{
    // only include the valid SNPs
    for (auto&& snp : existed_snps) {
        m_snp_selection_list.insert(snp.first);
    }
    m_exclude_snp = false;
    m_existed_snps =
        gen_snp_vector(geno, maf, info, hard_threshold, hard_coded, out_prefix);
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

    message.append(std::to_string(m_marker_ct) + " variant(s) included\n");
    if (verbose) reporter.report(message);
    m_snp_selection_list.clear();
}

void Genotype::load_snps(const std::string out_prefix,
                         const std::string& extract_file,
                         const std::string& exclude_file, const double geno,
                         const double maf, const double info,
                         const double hard_threshold, const bool hard_coded,
                         bool verbose, Reporter& reporter)
{
    if (!extract_file.empty()) {
        m_exclude_snp = false;
        m_snp_selection_list = load_snp_list(extract_file, reporter);
    }
    if (!exclude_file.empty()) {
        m_snp_selection_list = load_snp_list(exclude_file, reporter);
    }

    m_existed_snps =
        gen_snp_vector(geno, maf, info, hard_threshold, hard_coded, out_prefix);
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
            + " variant(s) excluded based on genotype missingness threshold");
    }
    if (m_num_maf_filter != 0) {
        message.append(std::to_string(m_num_maf_filter)
                       + " variant(s) excluded based on MAF threshold");
    }
    if (m_num_info_filter != 0) {
        message.append(std::to_string(m_num_maf_filter)
                       + " variant(s) excluded based on INFO score threshold");
    }
    message.append(std::to_string(m_marker_ct) + " variant(s) included\n");
    if (verbose) reporter.report(message);
    m_snp_selection_list.clear();
}

Genotype::~Genotype() {}

void Genotype::read_base(const Commander& c_commander, Region& region,
                         Reporter& reporter)
{
    // can assume region is of the same order as m_existed_snp

    const std::string input = c_commander.base_name();
    const bool beta = c_commander.beta();
    const bool fastscore = c_commander.fastscore();
    const bool no_full = c_commander.no_full();
    const double info_threshold = c_commander.base_info_score();
    const double maf_control = c_commander.maf_base_control();
    const double maf_case = c_commander.maf_base_case();
    std::vector<int> index = c_commander.index();
    // now coordinates obtained from target file instead. Coordinate information
    // in base file only use for validation
    bool gz_input = false;
    GZSTREAM_NAMESPACE::igzstream gz_snp_file;
    if (input.substr(input.find_last_of(".") + 1).compare("gz") == 0) {
        gz_snp_file.open(input.c_str());
        if (!gz_snp_file.good()) {
            std::string error_message =
                "ERROR: Cannot open base file (gz) to read!\n";
            throw std::runtime_error(error_message);
        }
        gz_input = true;
    }

    std::ifstream snp_file;
    if (!gz_input) {
        snp_file.open(input.c_str());
        if (!snp_file.is_open()) {
            std::string error_message =
                "ERROR: Cannot open base file: " + input;
            throw std::runtime_error(error_message);
        }
    }
    size_t max_index = index[+BASE_INDEX::MAX];
    std::string line;
    std::string message = "Base file: " + input + "\n";
    // category related stuff
    double threshold = (c_commander.fastscore()) ? c_commander.bar_upper()
                                                 : c_commander.upper();
    const double bound_start = c_commander.lower();
    const double bound_end = c_commander.upper();
    const double bound_inter = c_commander.inter();

    threshold = (!no_full) ? 1.0 : threshold;
    std::vector<std::string> token;

    // exclude indicates we don't want this SNP
    bool exclude = false;
    // Some QC counts
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

    std::unordered_set<std::string> dup_index;
    std::vector<int> exist_index; // try to use this as quick search
    // Actual reading the file, will do a bunch of QC
    size_t file_length = 0;
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
        snp_file.seekg(0, snp_file.end);
        file_length = snp_file.tellg();
        snp_file.clear();
        snp_file.seekg(0, snp_file.beg);
        if (!c_commander.is_index()) std::getline(snp_file, line);
    }
    std::unordered_set<int> unique_thresholds;
    double prev_progress = 0.0;

    // very ugly, might want to use polymorphism (is this the right word) as an
    // alternative solution
    while ((!gz_input && std::getline(snp_file, line))
           || (gz_input && std::getline(gz_snp_file, line)))
    {
        if (!gz_input) {
            double progress =
                (double) snp_file.tellg() / (double) (file_length) *100;
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

        if (token.size() <= max_index) {
            std::string error_message = line;
            error_message.append("\nMore index than column in data");
            throw std::runtime_error(error_message);
        }

        std::string rs_id = token[index[+BASE_INDEX::RS]];
        if (m_existed_snps_index.find(rs_id) != m_existed_snps_index.end()
            && dup_index.find(rs_id) == dup_index.end())
        {
            dup_index.insert(rs_id);
            auto&& cur_snp = m_existed_snps[m_existed_snps_index[rs_id]];
            int32_t chr_code = -1;
            if (index[+BASE_INDEX::CHR] >= 0) {
                chr_code =
                    get_chrom_code_raw(token[index[+BASE_INDEX::CHR]].c_str());
                if (((const uint32_t) chr_code) > m_max_code) {
                    if (chr_code != -1) {
                        if (chr_code >= MAX_POSSIBLE_CHROM) {
                            chr_code =
                                m_xymt_codes[chr_code - MAX_POSSIBLE_CHROM];
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
                else if (is_set(m_haploid_mask.data(), chr_code)
                         || chr_code == m_xymt_codes[X_OFFSET]
                         || chr_code == m_xymt_codes[Y_OFFSET])
                {
                    exclude = true;
                    num_haploid++;
                }
            }
            std::string ref_allele = (index[+BASE_INDEX::REF] >= 0)
                                         ? token[index[+BASE_INDEX::REF]]
                                         : "";
            std::string alt_allele = (index[+BASE_INDEX::ALT] >= 0)
                                         ? token[index[+BASE_INDEX::ALT]]
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
                        token[index[+BASE_INDEX::BP]].c_str());
                    if (loc < 0) {
                        std::string error_message =
                            "ERROR: " + rs_id + " has negative loci!\n";
                        throw std::runtime_error(error_message);
                    }
                }
                catch (const std::runtime_error& error)
                {
                    std::string error_message =
                        "ERROR: Non-numeric loci for " + rs_id + "!\n";
                    throw std::runtime_error(error_message);
                }
            }
            double maf = 1;
            bool maf_filtered = false;
            if (index[+BASE_INDEX::MAF] >= 0) {
                try
                {
                    maf = misc::convert<double>(
                        token[index[+BASE_INDEX::MAF]].c_str());
                }
                catch (const std::runtime_error& error)
                {
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
            if (index[+BASE_INDEX::MAF_CASE] >= 0) {
                try
                {
                    maf = misc::convert<double>(
                        token[index[+BASE_INDEX::MAF_CASE]].c_str());
                }
                catch (const std::runtime_error& error)
                {
                    num_maf_filter += !maf_filtered;
                    exclude = true;
                }
                if (maf < maf_case) {
                    num_maf_filter += !maf_filtered;
                    exclude = true;
                }
            }
            double info_score = 1;
            if (index[+BASE_INDEX::INFO] >= 0) {
                // obtain the INFO score
                try
                {
                    info_score = misc::convert<double>(
                        token[index[+BASE_INDEX::INFO]].c_str());
                }
                catch (const std::runtime_error& error)
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
            bool flipped = false;
            if (!cur_snp.matching(chr_code, loc, ref_allele, alt_allele,
                                  flipped))
            {
                // Mismatched SNPs
                num_mismatched++;
                exclude = true;
            }
            double pvalue = 2.0;
            try
            {
                pvalue = misc::convert<double>(token[index[+BASE_INDEX::P]]);
                if (pvalue < 0.0 || pvalue > 1.0) {
                    std::string error_message =
                        "ERROR: Invalid p-value for " + rs_id + "!\n";
                    throw std::runtime_error(error_message);
                }
                else if (pvalue > threshold)
                {
                    exclude = true;
                    num_excluded++;
                }
            }
            catch (const std::runtime_error& error)
            {
                exclude = true;
                num_not_converted++;
            }
            double stat = 0.0;
            try
            {
                stat = misc::convert<double>(token[index[+BASE_INDEX::STAT]]);
                if (stat < 0 && !beta) {
                    num_negative_stat++;
                    exclude = true;
                }
                else if (!beta)
                    stat = log(stat);
            }
            catch (const std::runtime_error& error)
            {
                num_not_converted++;
                exclude = true;
            }

            if (!alt_allele.empty() && ambiguous(ref_allele, alt_allele)) {
                num_ambiguous++;
                exclude = !m_keep_ambig;
            }
            if (!exclude) {
                int category = -1;
                double pthres = 0.0;
                if (fastscore) {
                    category = c_commander.get_category(pvalue);
                    pthres = c_commander.get_threshold(category);
                }
                else
                {
                    // calculate the threshold instead
                    if (pvalue > bound_end && !no_full) {
                        category = std::ceil((bound_end + 0.1 - bound_start)
                                             / bound_inter);
                        pthres = 1.0;
                    }
                    else
                    {
                        category =
                            std::ceil((pvalue - bound_start) / bound_inter);
                        category = (category < 0) ? 0 : category;
                        pthres = category * bound_inter + bound_start;
                    }
                }
                if (flipped) cur_snp.set_flipped();
                // ignore the SE as it currently serves no purpose
                exist_index.push_back(m_existed_snps_index[rs_id]);
                cur_snp.set_statistic(stat, 0.0, pvalue, category, pthres);
                if (unique_thresholds.find(category) == unique_thresholds.end())
                {
                    unique_thresholds.insert(category);
                    m_thresholds.push_back(pthres);
                }
                m_max_category =
                    (m_max_category < category) ? category : m_max_category;
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
    if (exist_index.size() != m_existed_snps.size())
    { // only do this if we need to remove some SNPs
        // we assume exist_index doesn't have any duplicated index
        std::sort(exist_index.begin(), exist_index.end());
        int start = (exist_index.empty()) ? -1 : exist_index.front();
        int end = start;
        std::vector<SNP>::iterator last = m_existed_snps.begin();
        ;
        for (auto&& ind : exist_index) {
            if (ind == start || ind - end == 1)
                end = ind; // try to perform the copy as a block
            else
            {
                std::copy(m_existed_snps.begin() + start,
                          m_existed_snps.begin() + end + 1, last);
                last += end + 1 - start;
                start = ind;
                end = ind;
            }
        }
        if (!exist_index.empty()) {
            std::copy(m_existed_snps.begin() + start,
                      m_existed_snps.begin() + end + 1, last);
            last += end + 1 - start;
        }
        m_existed_snps.erase(last, m_existed_snps.end());
    }
    m_existed_snps_index.clear();
    // now m_existed_snps is ok and can be used directly
    size_t vector_index = 0;
    // we do it here such that the m_existed_snps is sorted correctly
    size_t low_bound = 0, last_snp = 0;
    int prev_chr = 0, prev_loc = 0;
    for (auto&& cur_snp : m_existed_snps) {
        if (prev_chr != cur_snp.chr()) {
            prev_chr = cur_snp.chr();
            prev_loc = cur_snp.loc();
            low_bound = vector_index;
        }
        else if (cur_snp.loc() - prev_loc > clump_info.distance)
        {
            while (cur_snp.loc() - prev_loc > clump_info.distance
                   && low_bound < vector_index)
            {
                low_bound++;
                prev_loc = m_existed_snps[low_bound].loc();
            }
        }
        // now low_bound should be the first SNP where the core index SNP need
        // to read from
        cur_snp.set_low_bound(low_bound);
        // set this as the default
        cur_snp.set_up_bound(m_existed_snps.size());
        // update all previous SNPs that are out bounud
        while (m_existed_snps[last_snp].chr() != cur_snp.chr()
               || cur_snp.loc() - m_existed_snps[last_snp].loc()
                      > clump_info.distance)
        {
            m_existed_snps[last_snp++].set_up_bound(vector_index);
        }
        m_existed_snps_index[cur_snp.rs()] = vector_index++;
        // cur_snp.set_flag( region.check(cur_snp.chr(), cur_snp.loc()));
        cur_snp.set_flag(region);
    }
    // Suggest that we want to release memory
    // but this is only a suggest as this is non-binding request
    // Proper way of releasing memory will be to do swarp. Yet that
    // might lead to out of scope or some other error here?
    m_existed_snps.shrink_to_fit();
    m_region_size = region.size();
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
    reporter.report(message);
    if (m_existed_snps.size() == 0) {
        throw std::runtime_error("Error: No valid variant remaining");
    }
    m_num_threshold = unique_thresholds.size();
}

void Genotype::set_info(const Commander& c_commander, const bool ld)
{
    clump_info.p_value = c_commander.clump_p();
    clump_info.r2 = c_commander.clump_r2();
    clump_info.proxy = c_commander.proxy();
    clump_info.use_proxy = c_commander.use_proxy();
    clump_info.distance = c_commander.clump_dist();
    m_model = c_commander.model();
    m_missing_score = c_commander.get_missing_score();
    m_scoring = c_commander.get_score();
}

double Genotype::get_r2(bool core_missing, std::vector<uint32_t>& index_tots,
                        std::vector<uintptr_t>& index_data,
                        std::vector<uintptr_t>& genotype_vector)
{
    bool is_x = false;
    double freq11;
    double freq11_expected;
    double freq1x;
    double freq2x;
    double freqx1;
    double freqx2;
    double dxx;
    double r2 = -1.0;
    uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(m_founder_ct);
    uintptr_t founder_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_founder_ct);
    uint32_t counts[18];
    genovec_3freq(genotype_vector.data(), index_data.data(), founder_ctl2,
                  &(counts[0]), &(counts[1]), &(counts[2]));
    counts[0] = index_tots[0] - counts[0] - counts[1] - counts[2];
    genovec_3freq(genotype_vector.data(), &(index_data[founder_ctv2]),
                  founder_ctl2, &(counts[3]), &(counts[4]), &(counts[5]));
    counts[3] = index_tots[1] - counts[3] - counts[4] - counts[5];
    genovec_3freq(genotype_vector.data(), &(index_data[2 * founder_ctv2]),
                  founder_ctl2, &(counts[6]), &(counts[7]), &(counts[8]));
    counts[6] = index_tots[2] - counts[6] - counts[7] - counts[8];

    if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x, &freqx1,
                                &freqx2, &freq11))
    {
        freq11_expected = freqx1 * freq1x;
        dxx = freq11 - freq11_expected;
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
    return r2;
}

double Genotype::get_r2(bool core_missing, bool pair_missing,
                        std::vector<uint32_t>& core_tot,
                        std::vector<uint32_t>& pair_tot,
                        std::vector<uintptr_t>& genotype_vector,
                        std::vector<uintptr_t>& pair_genotype_vector)
{

    uint32_t counts[18];
    const uint32_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
    double freq11;
    double freq11_expected;
    double freq1x;
    double freq2x;
    double freqx1;
    double freqx2;
    double dxx;
    double r2 = 0.0;

    if (core_missing) {
        two_locus_count_table_zmiss1(genotype_vector.data(),
                                     pair_genotype_vector.data(), counts,
                                     founder_ctv3, pair_missing);
        if (pair_missing) {
            counts[2] = core_tot[0] - counts[0] - counts[1];
            counts[5] = core_tot[1] - counts[3] - counts[4];
        }
        counts[6] = pair_tot[0] - counts[0] - counts[3];
        counts[7] = pair_tot[1] - counts[1] - counts[4];
        counts[8] = pair_tot[2] - counts[2] - counts[5];
    }
    else
    {
        two_locus_count_table(genotype_vector.data(),
                              pair_genotype_vector.data(), counts, founder_ctv3,
                              pair_missing);
        if (pair_missing) {
            counts[2] = core_tot[0] - counts[0] - counts[1];
            counts[5] = core_tot[1] - counts[3] - counts[4];
            counts[8] = core_tot[2] - counts[6] - counts[7];
        }
    }
    // below, the false are basically is_x1 is_x2
    if (em_phase_hethet_nobase(counts, false, false, &freq1x, &freq2x, &freqx1,
                               &freqx2, &freq11))
    {
        r2 = -1;
    }
    else
    {

        freq11_expected = freqx1 * freq1x;
        dxx = freq11 - freq11_expected;
        // also want to avoid divide by 0
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
    return r2;
}


void Genotype::efficient_clumping(Genotype& reference, Reporter& reporter)
{
    /*
     *	Threading doesn't speed things up too much
     *	but will increase the memory usage, leading to the malloc
     *	error users experiencing. So we will now use only 1 thread.
     *	Also, we will go through each SNP one by one, reading
     *	the binary from the SNPs one at a time. Then remove it
     *	right after. This should be the fastest for PRSice. The
     *	problem with this method is that the clumping might take
     *	longer for PRSet due to I/O burden in the worst case
     *	scenario (repeatingly reading and discarding the same SNP)
     */
    /*
     * we now starts with
     * m_existed_snps -> SNPs in target
     * reference.m_existed_snps -> SNPs in reference
     * reference.m_existed_snps_index -> For quick finding the index of SNPs in
     * reference m_sort_by_p_index -> vector containing order of SNPs to process
     * Also, each SNP contains the bound information, where the low bound is the
     * Index of the SNP to start process from and the up bound is the SNP to end
     */
    // sample related constants
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uint32_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
    const uint32_t founder_ctsplit = 3 * founder_ctv3;
    const int num_snp = m_existed_snps.size();
    std::vector<uintptr_t> genotype_vector(unfiltered_sample_ctl * 2);
    std::vector<int> remain_core_snps;
    const double min_r2 = (clump_info.use_proxy)
                              ? std::min(clump_info.proxy, clump_info.r2)
                              : clump_info.r2;
    bool is_x = false;
    double freq11;
    double freq11_expected;
    double freq1x;
    double freq2x;
    double freqx1;
    double freqx2;
    double dxx;
    double r2 = -1.0;
    uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(m_founder_ct);
    uintptr_t founder_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_founder_ct);
    bool mismatch_error = false;
    bool flipped = false;
    int mismatch = 0;
    std::vector<uintptr_t> index_data(3 * founder_ctsplit + founder_ctv3);
    std::vector<uint32_t> index_tots(6);
    double prev_progress = 0.0;

    std::vector<uintptr_t> founder_include2(founder_ctv2, 0);
    fill_quatervec_55(m_founder_ct, founder_include2.data());
    std::unordered_set<int> unique_threshold;
    std::unordered_set<double> used_thresholds;
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
    if (llxx) {
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
    // if fail, return nullptr which will then get into the while loop
    while (!bigstack_ua) {
        malloc_size_mb = (malloc_size_mb * 3) / 4;
        if (malloc_size_mb < BIGSTACK_MIN_MB) {
            malloc_size_mb = BIGSTACK_MIN_MB;
        }
        bigstack_ua =
            (unsigned char*) malloc(malloc_size_mb * 1048576 * sizeof(char));
        if (bigstack_ua) {
            message.append(
                "Allocated " + std::to_string(malloc_size_mb)
                + " MB successfully, after larger attempt(s) failed\n");
        }
        else if (malloc_size_mb == BIGSTACK_MIN_MB)
        {
            throw std::runtime_error("Failed to allocate required memory");
        }
    }
    // force 64-byte align to make cache line sensitivity work
    reporter.report(message);
    bigstack_initial_base =
        (unsigned char*) round_up_pow2((uintptr_t) bigstack_ua, CACHELINE);
    uintptr_t* window_data = (uintptr_t*) bigstack_initial_base;
    uintptr_t* window_data_ptr = nullptr;
    unsigned char* g_bigstack_end =
        &(bigstack_initial_base[(malloc_size_mb * 1048576
                                 - (uintptr_t)(bigstack_initial_base
                                               - bigstack_ua))
                                & (~(CACHELINE - ONELU))]);
    uintptr_t max_window_size =
        (((uintptr_t) g_bigstack_end) - ((uintptr_t) bigstack_initial_base))
        / (founder_ctv2 * sizeof(intptr_t));
    g_bigstack_end = nullptr;
    uintptr_t cur_window_size = 0;
    if (!max_window_size) {
        throw std::runtime_error("ERROR: Not enough memory for clumping!");
    }
    for (size_t i_snp = 0; i_snp < m_sort_by_p_index.size(); ++i_snp) {
        double progress = (double) i_snp / (double) num_snp * 100;
        if (progress - prev_progress > 0.01) {
            fprintf(stderr, "\rClumping Progress: %03.2f%%", progress);
            prev_progress = progress;
        }
        auto&& cur_snp_index = m_sort_by_p_index[i_snp];
        // skip any SNPs that are clumped
        auto&& cur_snp = m_existed_snps[cur_snp_index];
        if (cur_snp.clumped() || cur_snp.p_value() > clump_info.p_value)
            continue;
        auto&& target_pair = reference.m_existed_snps_index.find(cur_snp.rs());
        if (target_pair == reference.m_existed_snps_index.end()) continue;
        // Any SNP with p-value less than clump-p will be ignored
        // because they can never be an index SNP and thus are not of our
        // interested

        auto&& ref_snp = reference.m_existed_snps[target_pair->second];
        if (!cur_snp.matching(ref_snp.chr(), ref_snp.loc(), ref_snp.ref(),
                              ref_snp.alt(), flipped))
        {
            mismatch++;
            if (!mismatch_error) {
                std::string message =
                    "WARNING: Mismatched SNPs between LD reference and target!";
                message.append("Will use information from target file");
                message.append("You should check the files are based on the "
                               "same genome build\n");
                reporter.report(message);
                mismatch_error = true;
            }
        }
        // location cannot be negative
        // if this become a problem, we can also put it forward
        // and skip it just like SNPs with large p-values
        assert(cur_snp.loc() >= 0);

        // set the missing information
        // contain_missing == 3 = has missing
        size_t start = cur_snp.low_bound();
        size_t end = cur_snp.up_bound();
        window_data_ptr = window_data;
        cur_window_size = 0;
        for (size_t i_pair = start; i_pair < cur_snp_index; i_pair++) {
            auto&& pair_snp = m_existed_snps[i_pair];
            if (pair_snp.clumped() || pair_snp.p_value() > clump_info.p_value)
                continue;
            auto&& pair_index =
                reference.m_existed_snps_index.find(pair_snp.rs());
            if (pair_index == reference.m_existed_snps_index.end()) continue;
            auto&& ref_pair_snp = reference.m_existed_snps[pair_index->second];
            window_data_ptr[founder_ctv2 - 2] = 0;
            window_data_ptr[founder_ctv2 - 1] = 0;
            if (++cur_window_size == max_window_size) {
                throw std::runtime_error("ERROR: Out of memory!");
            }
            reference.read_genotype(window_data_ptr, ref_pair_snp,
                                    ref_pair_snp.file_name());
            window_data_ptr = &(window_data_ptr[founder_ctv2]);
        }

        window_data_ptr[founder_ctv2 - 2] = 0;
        window_data_ptr[founder_ctv2 - 1] = 0;
        std::fill(index_data.begin(), index_data.end(), 0);

        reference.read_genotype(window_data_ptr, cur_snp, cur_snp.file_name());
        vec_datamask(m_founder_ct, 0, window_data_ptr, founder_include2.data(),
                     index_data.data());
        index_tots[0] = popcount2_longs(index_data.data(), founder_ctl2);
        vec_datamask(m_founder_ct, 2, window_data_ptr, founder_include2.data(),
                     &(index_data[founder_ctv2]));
        index_tots[1] =
            popcount2_longs(&(index_data[founder_ctv2]), founder_ctl2);
        vec_datamask(m_founder_ct, 3, window_data_ptr, founder_include2.data(),
                     &(index_data[2 * founder_ctv2]));
        index_tots[2] =
            popcount2_longs(&(index_data[2 * founder_ctv2]), founder_ctl2);

        if (++cur_window_size == max_window_size) {
            throw std::runtime_error("ERROR: Out of memory!");
        }
        window_data_ptr = window_data;
        for (size_t i_pair = start; i_pair < cur_snp_index; i_pair++) {
            auto&& pair_snp = m_existed_snps[i_pair];
            if (pair_snp.clumped() || pair_snp.p_value() > clump_info.p_value)
                continue;
            auto&& pair_index =
                reference.m_existed_snps_index.find(pair_snp.rs());
            if (pair_index == reference.m_existed_snps_index.end()) continue;
            auto&& ref_pair_snp = reference.m_existed_snps[pair_index->second];

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
            r2 = -1;
            if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x,
                                        &freqx1, &freqx2, &freq11))
            {
                freq11_expected = freqx1 * freq1x;
                dxx = freq11 - freq11_expected;
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
                cur_snp.clump(ref_pair_snp, r2, clump_info.proxy);
            }
            window_data_ptr = &(window_data_ptr[founder_ctv2]);
        }

        for (size_t i_pair = cur_snp_index + 1; i_pair < end; ++i_pair) {
            window_data_ptr = window_data;
            auto&& pair_snp = m_existed_snps[i_pair];
            if (pair_snp.clumped() || pair_snp.p_value() > clump_info.p_value)
                continue;
            auto&& pair_index =
                reference.m_existed_snps_index.find(pair_snp.rs());
            if (pair_index == reference.m_existed_snps_index.end()) continue;
            auto&& ref_pair_snp = reference.m_existed_snps[pair_index->second];
            window_data_ptr[founder_ctv2 - 2] = 0;
            window_data_ptr[founder_ctv2 - 1] = 0;

            reference.read_genotype(window_data_ptr, ref_pair_snp,
                                    ref_pair_snp.file_name());
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
            r2 = -1;
            if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x,
                                        &freqx1, &freqx2, &freq11))
            {
                freq11_expected = freqx1 * freq1x;
                dxx = freq11 - freq11_expected;
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
                cur_snp.clump(ref_pair_snp, r2, clump_info.proxy);
            }
        }
        cur_snp.set_clumped();
        remain_core_snps.push_back(cur_snp_index);
        double thres = cur_snp.get_threshold();
        if (used_thresholds.find(thres) == used_thresholds.end()) {
            used_thresholds.insert(thres);
            m_thresholds.push_back(thres);
        }
        if (unique_threshold.find(cur_snp.category()) == unique_threshold.end())
        {
            unique_threshold.insert(cur_snp.category());
        }
    }
    fprintf(stderr, "\rClumping Progress: %03.2f%%\n\n", 100.0);

    window_data = nullptr;
    window_data_ptr = nullptr;
    delete[] bigstack_ua;
    bigstack_ua = nullptr;
    bigstack_initial_base = nullptr;

    m_existed_snps_index.clear();
    m_num_threshold = unique_threshold.size();
    if (remain_core_snps.size() != m_existed_snps.size()) {
        // only do this if we need to remove some SNPs
        // we assume exist_index doesn't have any duplicated index
        std::sort(remain_core_snps.begin(), remain_core_snps.end());
        int start = (remain_core_snps.empty()) ? -1 : remain_core_snps.front();
        int end = start;
        std::vector<SNP>::iterator last = m_existed_snps.begin();

        for (auto&& ind : remain_core_snps) {
            if (ind == start || ind - end == 1)
                end = ind; // try to perform the copy as a block
            else
            {
                std::copy(m_existed_snps.begin() + start,
                          m_existed_snps.begin() + end + 1, last);
                last += end + 1 - start;
                start = ind;
                end = ind;
            }
        }
        if (!remain_core_snps.empty()) {
            std::copy(m_existed_snps.begin() + start,
                      m_existed_snps.begin() + end + 1, last);
            last += end + 1 - start;
        }
        m_existed_snps.erase(last, m_existed_snps.end());
    }
    m_existed_snps.shrink_to_fit();
    m_existed_snps_index.clear();
    // no longer require the m_existed_snps_index
    message = "";
    if (mismatch != 0) {
        message.append("There are a total of " + std::to_string(mismatch)
                       + " mismatched variant(s) between the reference panel "
                         "and the target genotype\n");
    }
    message.append("Number of variant(s) after clumping : "
                   + std::to_string(m_existed_snps.size()) + "\n");
    reporter.report(message);
}

bool Genotype::sort_by_p()
{
    if (m_existed_snps.size() == 0) return false;
    m_sort_by_p_index = SNP::sort_by_p_chr(m_existed_snps);
    return true;
}

bool Genotype::prepare_prsice()
{
    if (m_existed_snps.size() == 0) return false;
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

bool Genotype::get_score(int& cur_index, int& cur_category,
                         double& cur_threshold, size_t& num_snp_included,
                         const size_t region_index,
                         const bool require_statistic)
{
    if (m_existed_snps.size() == 0 || cur_index == m_existed_snps.size())
        return false;
    int end_index = 0;
    bool ended = false;
    if (cur_index == -1) // first run
    {
        cur_index = 0;
        cur_category = m_existed_snps[cur_index].category();
    }
    cur_threshold = m_existed_snps[cur_index].get_threshold();
    // existed snp should be sorted such that the SNPs should be
    // access sequentially
    for (size_t i = cur_index; i < m_existed_snps.size(); ++i) {
        if (m_existed_snps[i].category() != cur_category) {
            end_index = i;
            ended = true;
            break;
        }
        //		// Use as part of the output
        if (m_existed_snps[i].in(region_index)) num_snp_included++;
    }
    if (!ended) {
        end_index = m_existed_snps.size();
        cur_category = m_existed_snps.back().category();
    }
    else
        cur_category = m_existed_snps[end_index].category();
    read_score(cur_index, end_index, region_index);
    cur_index = end_index;
    if (require_statistic) {
        misc::RunningStat rs;
        for (auto&& sample : m_sample_names) {
            if (!sample.included || !sample.has_pheno) continue;
            if (sample.num_snp == 0) {
                rs.push(0.0);
            }
            else
            {
                rs.push(sample.prs / (double) sample.num_snp);
            }
        }
        m_mean_score = rs.mean();
        m_score_sd = rs.sd();
    }
    return true;
}

void Genotype::print_snp(std::string& output, double threshold,
                         const size_t region_index)
{
    std::ofstream snp_out;
    snp_out.open(output);
    if (!snp_out.is_open()) {
        std::string error_message =
            "ERROR: Cannot open file: " + output + " to write";
        throw std::runtime_error(error_message);
    }
    for (auto&& snp : m_existed_snps) {
        snp_out << snp.rs() << "\t" << snp.chr() << "\t" << snp.loc() << "\t"
                << snp.p_value();
        if (snp.get_threshold() <= threshold && snp.in(region_index)) {
            snp_out << "\tY";
        }
        else
        {
            snp_out << "\tN";
        }
        snp_out << std::endl;
    }
    snp_out.close();
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
    // Returns 1 if at least one SNP is monomorphic over all valid observations;
    // returns 0 otherwise, and fills all frequencies using the maximum
    // likelihood solution to the cubic equation.
    // (We're discontinuing most use of EM phasing since better algorithms have
    // been developed, but the two marker case is mathematically clean and fast
    // enough that it'll probably remain useful as an input for some of those
    // better algorithms...)
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
    // shouldn't have to worry about subtractive cancellation problems here
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
    // the following four values should all be guaranteed nonzero except in the
    // NAN case
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
                // assert((sol_start_idx < sol_end_idx) &&sol_start_idx < 3);
            }
            if (sol_start_idx == sol_end_idx) {
                // Lost a planet Master Obi-Wan has.  How embarrassing...
                // lost root must be a double root at one of the boundary
                // points, just check their likelihoods
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
            // okay to NOT count suboptimal boundary points because they don't
            // permit direction changes within the main interval this should
            // exactly match haploview_blocks_classify()'s D sign check
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
    // if is_x1 and/or is_x2 is set, counts[9]..[17] are male-only counts.
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
                // unlike the zmiss variant, this apparently does not suffer
                // from enough register spill to justify shrinking the inner
                // loop
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
                loader22 =
                    _mm_and_si128(_mm_srli_epi64(loader21, 1), m1); // half2
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
