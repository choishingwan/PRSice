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

std::mutex Genotype::clump_mtx;

/* we don't want to do this as this seems to be root of some problems
void Genotype::initialize()
{
// Don't use this. For some reason, this does not work
    m_founder_ctl = BITCT_TO_WORDCT(m_founder_ct);
    m_founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
    m_founder_ctsplit = 3 * m_founder_ctv3;
    m_final_mask = get_final_mask(m_founder_ct);
    m_unfiltered_marker_ctl = BITCT_TO_WORDCT(m_unfiltered_marker_ct);
    m_marker_exclude = new uintptr_t[m_unfiltered_marker_ctl];
    std::memset(m_marker_exclude, 0x0,
m_unfiltered_marker_ctl*sizeof(uintptr_t)); m_marker_ct = m_existed_snps.size();
}
*/
void Genotype::init_chr(int num_auto, bool no_x, bool no_y, bool no_xy,
                        bool no_mt)
{
    // this initialize haploid mask as the maximum possible number

    if (num_auto < 0) {
        num_auto                = -num_auto;
        m_autosome_ct           = num_auto;
        m_xymt_codes[X_OFFSET]  = -1;
        m_xymt_codes[Y_OFFSET]  = -1;
        m_xymt_codes[XY_OFFSET] = -1;
        m_xymt_codes[MT_OFFSET] = -1;
        m_max_code              = num_auto;
        fill_all_bits(((uint32_t) num_auto) + 1, m_haploid_mask);
    }
    else
    {
        m_autosome_ct           = num_auto;
        m_xymt_codes[X_OFFSET]  = num_auto + 1;
        m_xymt_codes[Y_OFFSET]  = num_auto + 2;
        m_xymt_codes[XY_OFFSET] = num_auto + 3;
        m_xymt_codes[MT_OFFSET] = num_auto + 4;
        set_bit(num_auto + 1, m_haploid_mask);
        set_bit(num_auto + 2, m_haploid_mask);
        if (no_x) {
            m_xymt_codes[X_OFFSET] = -1;
            clear_bit(num_auto + 1, m_haploid_mask);
        }
        if (no_y) {
            m_xymt_codes[Y_OFFSET] = -1;
            clear_bit(num_auto + 2, m_haploid_mask);
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
    fill_all_bits(m_autosome_ct + 1, m_chrom_mask);
    for (uint32_t xymt_idx = 0; xymt_idx < XYMT_OFFSET_CT; ++xymt_idx) {
        int32_t cur_code = m_xymt_codes[xymt_idx];
        if (cur_code != -1) {
            set_bit(m_xymt_codes[xymt_idx], m_chrom_mask);
        }
    }
    m_chrom_start.resize(m_max_code); // 1 extra for the info
}

void Genotype::set_genotype_files(std::string prefix)
{
    if (prefix.find("#") != std::string::npos) {
        for (size_t chr = 1; chr < m_max_code; ++chr) {
            std::string name = prefix;
            misc::replace_substring(name, "#", std::to_string(chr));
            m_genotype_files.push_back(name);
        }
    }
    else
    {
        m_genotype_files.push_back(prefix);
    }
}

void Genotype::update_include(const std::vector<Sample>& inclusion)
{
    m_sample_ct                     = 0;
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    std::fill(m_sample_include, m_sample_include + unfiltered_sample_ctl, 0);
    // std::memset(m_sample_include, 0x0,
    // m_unfiltered_sample_ctl*sizeof(uintptr_t));
    for (size_t i_sample = 0; i_sample < inclusion.size(); ++i_sample) {
        if (IS_SET(m_founder_info, i_sample) && inclusion[i_sample].included) {
            SET_BIT(i_sample, m_sample_include);
            m_sample_ct++;
        }
    }
}

std::unordered_set<std::string> Genotype::load_snp_list(std::string input)
{
    std::ifstream in;
    in.open(input.c_str());
    if (!in.is_open()) {
        std::string error_message = "ERROR: Cannot open file: " + input;
        throw std::runtime_error(error_message);
    }
    std::string                     line;
    std::unordered_set<std::string> result;
    bool                            error = false;
    while (std::getline(in, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (token[0].compare(".") == 0) {
            if (!error) {
                error = true;
                fprintf(stderr, "WARNING: Some SNPs from the "
                                "extraction/exclusion list has rs-id of .\n");
                fprintf(stderr, "         They will be excluded unless the "
                                "file contains at least 3 columns\n");
                fprintf(stderr, "         When 3 columns is provided, we will "
                                "assume the second and third columns\n");
                fprintf(stderr, "         are chromosome and coordinates "
                                "respectively and will generate an rsid\n");
                fprintf(stderr, "         as chr:loc\n");
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
                                                   bool        ignore_fid)
{
    std::ifstream in;
    in.open(input.c_str());
    if (!in.is_open()) {
        std::string error_message = "ERROR: Cannot open file: " + input;
        throw std::runtime_error(error_message);
    }
    std::string                     line;
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

Genotype::Genotype(std::string prefix, std::string remove_sample,
                   std::string keep_sample, std::string extract_snp,
                   std::string exclude_snp, bool ignore_fid, int num_auto,
                   bool no_x, bool no_y, bool no_xy, bool no_mt,
                   bool keep_ambig, const size_t thread, bool verbose)
{
    m_thread = thread;
    if (!remove_sample.empty()) {
        m_remove_sample      = true;
        m_remove_sample_list = load_ref(remove_sample, ignore_fid);
    }
    if (!keep_sample.empty()) {
        m_keep_sample      = false;
        m_keep_sample_list = load_ref(keep_sample, ignore_fid);
    }
    if (!extract_snp.empty()) {
        m_extract_snp      = true;
        m_extract_snp_list = load_snp_list(extract_snp);
    }
    if (!exclude_snp.empty()) {
        m_exclude_snp      = true;
        m_exclude_snp_list = load_snp_list(exclude_snp);
    }

    /** setting the chromosome information **/
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    // we are not using the following script for now as we only support human
    m_haploid_mask = new uintptr_t[CHROM_MASK_WORDS];
    fill_ulong_zero(CHROM_MASK_WORDS, m_haploid_mask);
    m_chrom_mask = new uintptr_t[CHROM_MASK_WORDS];
    fill_ulong_zero(CHROM_MASK_WORDS, m_chrom_mask);
    // now initialize the chromosome
    init_chr(num_auto, no_x, no_y, no_xy, no_mt);

    set_genotype_files(prefix);
    m_sample_names = load_samples(ignore_fid);
    m_existed_snps = load_snps();
    if (verbose) {
        fprintf(stderr, "%zu people (%zu males, %zu females) included\n",
                m_unfiltered_sample_ct, m_num_male, m_num_female);
        if (m_num_ambig != 0 && !keep_ambig)
            fprintf(stderr, "%u ambiguous variants excluded\n", m_num_ambig);
        else if (m_num_ambig != 0)
            fprintf(stderr, "%u ambiguous variants kept\n", m_num_ambig);
        fprintf(stderr, "%zu variants included\n", m_marker_ct);
    }
}

Genotype::~Genotype()
{
    // TODO Auto-generated destructor stub
    if (m_founder_info != nullptr) delete[] m_founder_info;
    // if(m_sex_male != nullptr) delete [] m_sex_male;
    if (m_sample_include != nullptr) delete[] m_sample_include;
    // if(m_marker_exclude != nullptr) delete [] m_marker_exclude;
    if (m_haploid_mask != nullptr) delete[] m_haploid_mask;
    if (m_chrom_mask != nullptr) delete[] m_chrom_mask;
    // if(m_tmp_genotype != nullptr) delete [] m_tmp_genotype;
}

void Genotype::read_base(const Commander& c_commander, Region& region)
{
    // can assume region is of the same order as m_existed_snp
    m_scoring                   = c_commander.get_scoring();
    const std::string input     = c_commander.base_name();
    const bool        beta      = c_commander.beta();
    const bool        fastscore = c_commander.fastscore();
    const bool        full      = c_commander.full();
    std::vector<int>  index =
        c_commander.index(); // more appropriate for commander
    // now coordinates obtained from target file instead. Coordinate information
    // in base file only use for validation
    std::ifstream snp_file;
    snp_file.open(input.c_str());
    if (!snp_file.is_open()) {
        std::string error_message = "ERROR: Cannot open base file: " + input;
        throw std::runtime_error(error_message);
    }
    size_t      max_index = index[+BASE_INDEX::MAX];
    std::string line;
    if (!c_commander.has_index()) std::getline(snp_file, line);

    // category related stuff
    double threshold = (c_commander.fastscore()) ? c_commander.bar_upper()
                                                 : c_commander.upper();
    double bound_start = c_commander.lower();
    double bound_end   = c_commander.upper();
    double bound_inter = c_commander.inter();

    threshold = (full) ? 1.0 : threshold;
    std::vector<std::string> token;

    bool exclude = false;
    // Some QC countss
    size_t num_duplicated = 0;
    size_t num_excluded   = 0;
    size_t num_ambiguous  = 0;
    size_t num_haploid    = 0;

    size_t num_not_found     = 0;
    size_t num_mismatched    = 0;
    size_t num_not_converted = 0; // this is for NA
    size_t num_negative_stat = 0;

	std::unordered_set<std::string> dup_index;
	std::vector<int> exist_index; // try to use this as quick search
	// Actual reading the file, will do a bunch of QC
	snp_file.seekg (0, snp_file.end);
	size_t file_length = snp_file.tellg();
	snp_file.seekg (0, snp_file.beg);
	std::unordered_set<int> unique_thresholds;
	double prev_progress = 0.0;
	while (std::getline(snp_file, line))
	{
		double progress = (double) snp_file.tellg() / (double) (file_length)*100;
		if(progress-prev_progress > 0.01)
		{
			fprintf(stderr, "\rReading %03.2f%%", progress);
			prev_progress = progress;
		}
		misc::trim(line);
		if (line.empty()) continue;
		exclude = false;
		token = misc::split(line);

		if (token.size() <= max_index)
		{
			std::cerr << line << std::endl;
			throw std::runtime_error("More index than column in data");
		}

		std::string rs_id = token[index[+BASE_INDEX::RS]];
		if(m_existed_snps_index.find(rs_id)!=m_existed_snps_index.end() && dup_index.find(rs_id)==dup_index.end())
		{
			dup_index.insert(rs_id);
			auto &&cur_snp = m_existed_snps[m_existed_snps_index[rs_id]];
			int32_t chr_code = -1;
			if (index[+BASE_INDEX::CHR] >= 0)
			{
				chr_code = get_chrom_code_raw(token[index[+BASE_INDEX::CHR]].c_str());
				if (((const uint32_t)chr_code) > m_max_code) {
					if (chr_code != -1) {
						if (chr_code >= MAX_POSSIBLE_CHROM) {
							chr_code= m_xymt_codes[chr_code - MAX_POSSIBLE_CHROM];
							// this is the sex chromosomes
							// we don't need to output the error as they will be filtered out before by the
							// genotype read anyway
							exclude =true;
							num_haploid++;
						}
						else
						{
							std::string error_message ="ERROR: Cannot parse chromosome code: "
									+ token[index[+BASE_INDEX::CHR]];
							throw std::runtime_error(error_message);
						}
					}
				}
				else if(is_set(m_haploid_mask, chr_code) || chr_code==m_xymt_codes[X_OFFSET] ||
						chr_code==m_xymt_codes[Y_OFFSET])
				{
					// again, doesn't need to provide this message
					// the only time this will happen is when the target & base has different chromosome information
					//if(!hap_error) fprintf(stderr, "\nWARNING: Currently not supporting haploid chromosome and sex chromosomes\n");
					exclude = true;
					num_haploid++;
				}
			}
			std::string ref_allele = (index[+BASE_INDEX::REF] >= 0) ? token[index[+BASE_INDEX::REF]] : "";
			std::string alt_allele = (index[+BASE_INDEX::ALT] >= 0) ? token[index[+BASE_INDEX::ALT]] : "";
			std::transform(ref_allele.begin(), ref_allele.end(),ref_allele.begin(), ::toupper);
			std::transform(alt_allele.begin(), alt_allele.end(),alt_allele.begin(), ::toupper);
			int loc = -1;
			if (index[+BASE_INDEX::BP] >= 0)
			{
				// obtain the SNP coordinate
				try {
					loc = misc::convert<int>( token[index[+BASE_INDEX::BP]].c_str());
					if (loc < 0)
					{
						std::string error_message = "ERROR: "+rs_id+" has negative loci!\n";
						throw std::runtime_error(error_message);
					}
				} catch (const std::runtime_error &error) {
					std::string error_message = "ERROR: Non-numeric loci for "+rs_id+"!\n";
					throw std::runtime_error(error_message);
				}
			}
			bool flipped = false;
			if(!cur_snp.matching(chr_code, loc, ref_allele, alt_allele, flipped))
			{
				num_mismatched++;
				exclude = true; // hard check, as we can't tell if that is correct or not anyway
			}
			double pvalue = 2.0;
			try{
				pvalue = misc::convert<double>( token[index[+BASE_INDEX::P]]);
				if (pvalue < 0.0 || pvalue > 1.0)
				{
					std::string error_message = "ERROR: Invalid p-value for "+rs_id+"!\n";
					throw std::runtime_error(error_message);
				}
				else if (pvalue > threshold)
				{
					exclude = true;
					num_excluded++;
				}
			}catch (const std::runtime_error& error) {
				exclude = true;
				num_not_converted++;
			}
			double stat = 0.0;
			try {
				stat = misc::convert<double>( token[index[+BASE_INDEX::STAT]]);
				if(stat <0 && !beta)
				{
					num_negative_stat++;
					exclude = true;
				}
				else if (!beta) stat = log(stat);
			} catch (const std::runtime_error& error) {
				num_not_converted++;
				exclude = true;
			}

            if (!alt_allele.empty() && ambiguous(ref_allele, alt_allele)) {
                num_ambiguous++;
                exclude = !filter.keep_ambig;
            }
            if (!exclude) {
                int    category = -1;
                double pthres   = 0.0;
                if (fastscore) {
                    category = c_commander.get_category(pvalue);
                    pthres   = c_commander.get_threshold(category);
                }
                else
                {
                    // calculate the threshold instead
                    if (pvalue > bound_end && full) {
                        category = std::ceil((bound_end + 0.1 - bound_start)
                                             / bound_inter);
                        pthres   = 1.0;
                    }
                    else
                    {
                        category =
                            std::ceil((pvalue - bound_start) / bound_inter);
                        category = (category < 0) ? 0 : category;
                        pthres   = category * bound_inter + bound_start;
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
    snp_file.close();

    fprintf(stderr, "\rReading %03.2f%%\n", 100.0);
    if (exist_index.size() != m_existed_snps.size())
    { // only do this if we need to remove some SNPs
        // we assume exist_index doesn't have any duplicated index
        std::sort(exist_index.begin(), exist_index.end());
        int start = (exist_index.empty()) ? -1 : exist_index.front();
        int end   = start;
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
                end   = ind;
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
    for (auto&& cur_snp : m_existed_snps) // should be in the correct order
    {
        m_existed_snps_index[cur_snp.rs()] = vector_index++;
        // cur_snp.set_flag( region.check(cur_snp.chr(), cur_snp.loc()));
        cur_snp.set_flag(region);
    }
    m_region_size = region.size();

    if (num_duplicated)
        fprintf(stderr, "%zu duplicated variant(s) in base file\n",
                num_duplicated);
    if (num_excluded)
        fprintf(stderr, "%zu variant(s) excluded due to p-value threshold\n",
                num_excluded);
    if (num_ambiguous) {
        fprintf(stderr, "%zu ambiguous variant(s)", num_ambiguous);
        if (!filter.keep_ambig) {
            fprintf(stderr, " excluded");
        }
        fprintf(stderr, "\n");
    }
    if (num_haploid)
        fprintf(stderr, "%zu variant(s) located on haploid chromosome\n",
                num_haploid);
    if (num_not_found)
        fprintf(stderr, "%zu variant(s) not found in target file\n",
                num_not_found);
    if (num_mismatched)
        fprintf(stderr, "%zu mismatched variant(s) excluded\n", num_mismatched);
    if (num_not_converted)
        fprintf(stderr, "%zu NA stat/p-value observed\n", num_not_converted);
    if (num_negative_stat)
        fprintf(
            stderr,
            "%zu negative statistic observed. Please make sure it is really "
            "OR\n",
            num_negative_stat);
    fprintf(stderr, "%zu total SNPs included from base file\n\n",
            m_existed_snps.size());
    clump_info.p_value           = c_commander.clump_p();
    clump_info.r2                = c_commander.clump_r2();
    clump_info.proxy             = c_commander.proxy();
    clump_info.use_proxy         = c_commander.use_proxy();
    clump_info.distance          = c_commander.clump_dist();
    filter.filter_geno           = c_commander.filter_geno();
    filter.filter_maf            = c_commander.filter_maf();
    filter.filter_info           = c_commander.filter_info();
    filter.filter_hard_threshold = c_commander.filter_hard_threshold();
    filter.geno                  = c_commander.geno();
    filter.info_score            = c_commander.info_score();
    filter.maf                   = c_commander.maf();
    filter.hard_threshold        = c_commander.hard_threshold();
    filter.use_hard              = c_commander.hard_coding();
    m_num_threshold              = unique_thresholds.size();
}

void Genotype::clump(Genotype& reference)
{
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    uint32_t  founder_ctv3          = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
    uint32_t  founder_ctsplit       = 3 * founder_ctv3;
    auto&&    cur_snp               = m_existed_snps.front();
    size_t    bp_of_core            = cur_snp.loc();
    int       prev_chr              = cur_snp.chr();
    int       mismatch              = 0;
    int       total                 = 0;
    size_t    core_genotype_index   = 0;
    int       progress              = 0;
    int       num_snp               = m_existed_snps.size();
    int       begin_index           = 0;
    bool      mismatch_error        = false;
    bool      require_clump         = false;
    double    prev_progress         = 0.0;
    std::unordered_set<int> overlapped_snps;
    uintptr_t*              genotype = new uintptr_t[unfiltered_sample_ctl * 2];
    for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp) {
        auto&& snp    = m_existed_snps[i_snp];
        auto&& target = reference.m_existed_snps_index.find(snp.rs());
        if (target == reference.m_existed_snps_index.end())
            continue; // only work on SNPs that are in both
        auto&& ld_snp = reference.m_existed_snps[target->second];
        total++;
        bool flipped = false; // place holder
        if (!snp.matching(ld_snp.chr(), ld_snp.loc(), ld_snp.ref(),
                          ld_snp.alt(), flipped))
        {
            mismatch++;
            if (!mismatch_error) {
                fprintf(stderr, "WARNING: Mismatched SNPs between LD reference "
                                "and target!\n");
                fprintf(stderr,
                        "         Will use information from target file\n");
                fprintf(stderr, "         You should check the files are based "
                                "on same genome build\n");
                mismatch_error = true;
            }
        }
        assert(snp.loc() >= 0);
        if (prev_chr != snp.chr()) {
            perform_clump(core_genotype_index, begin_index, i_snp,
                          require_clump);
            prev_chr = snp.chr();
        }
        else if ((snp.loc() - bp_of_core) > clump_info.distance)
        {
            perform_clump(core_genotype_index, begin_index, i_snp,
                          require_clump);
        }
        // Now read in the current SNP

        overlapped_snps.insert(i_snp);
        std::fill(genotype, genotype + unfiltered_sample_ctl * 2, 0);
        // std::memset(genotype, 0x0,
        // m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
        reference.read_genotype(genotype, snp.snp_id(), snp.file_name());

        uintptr_t ulii = founder_ctsplit * sizeof(intptr_t)
                         + 2 * sizeof(int32_t)
                         + (m_marker_ct - 1) * 2 * sizeof(double);
        uintptr_t* geno1 = new uintptr_t[3 * founder_ctsplit + founder_ctv3];
        std::fill(geno1, geno1 + (3 * founder_ctsplit + founder_ctv3), 0);
        // std::memset(geno1, 0x0, (3*m_founder_ctsplit
        // +m_founder_ctv3)*sizeof(uintptr_t));
        load_and_split3(genotype, m_founder_ct, geno1, founder_ctv3, 0, 0, 1,
                        &ulii);

        snp.set_clump_geno(geno1, ulii);

        if (!require_clump && snp.p_value() < clump_info.p_value)
        { // Set this as the core SNP
            bp_of_core          = snp.loc();
            core_genotype_index = i_snp; // Should store the index on genotype
            require_clump       = true;
        }
        double cur_progress = (double) progress++ / (double) (num_snp) *100.0;
        if (cur_progress - prev_progress > 0.01) {
            prev_progress = cur_progress;
            fprintf(stderr, "\rClumping Progress: %03.2f%%", cur_progress);
        }
    }
    delete[] genotype;
    if (core_genotype_index != m_existed_snps.size()) {
        // this make sure this will be the last
        perform_clump(core_genotype_index, begin_index, m_existed_snps.size(),
                      require_clump);
    }

    fprintf(stderr, "\rClumping Progress: %03.2f%%\n\n", 100.0);

    std::vector<int>           remain_snps;
    std::unordered_set<double> used_thresholds;
    m_existed_snps_index.clear();
    std::vector<size_t>     p_sort_order = SNP::sort_by_p(m_existed_snps);
    bool                    proxy        = clump_info.use_proxy;
    std::unordered_set<int> unique_threshold;
    m_thresholds.clear();
    for (auto&& i_snp : p_sort_order) {
        if (overlapped_snps.find(i_snp) != overlapped_snps.end()
            && m_existed_snps[i_snp].p_value() <= clump_info.p_value)
        {
            if (proxy && !m_existed_snps[i_snp].clumped()) {
                m_existed_snps[i_snp].proxy_clump(m_existed_snps,
                                                  clump_info.proxy);
                double thres = m_existed_snps[i_snp].get_threshold();
                if (used_thresholds.find(thres) == used_thresholds.end()) {
                    used_thresholds.insert(thres);
                    m_thresholds.push_back(thres);
                }
                remain_snps.push_back(i_snp);
                if (unique_threshold.find(m_existed_snps[i_snp].category())
                    == unique_threshold.end())
                {
                    unique_threshold.insert(m_existed_snps[i_snp].category());
                }
            }
            else if (!m_existed_snps[i_snp].clumped())
            {
                m_existed_snps[i_snp].clump(m_existed_snps);
                double thres = m_existed_snps[i_snp].get_threshold();
                if (used_thresholds.find(thres) == used_thresholds.end()) {
                    used_thresholds.insert(thres);
                    m_thresholds.push_back(thres);
                }
                remain_snps.push_back(i_snp);
                if (unique_threshold.find(m_existed_snps[i_snp].category())
                    == unique_threshold.end())
                {
                    unique_threshold.insert(m_existed_snps[i_snp].category());
                }
            }
        }
        else if (m_existed_snps[i_snp].p_value() >= clump_info.p_value)
            break;
    }
    m_num_threshold = unique_threshold.size();
    if (remain_snps.size() != m_existed_snps.size())
    { // only do this if we need to remove some SNPs
        // we assume exist_index doesn't have any duplicated index
        std::sort(remain_snps.begin(), remain_snps.end());
        int start = (remain_snps.empty()) ? -1 : remain_snps.front();
        int end   = start;
        std::vector<SNP>::iterator last = m_existed_snps.begin();
        ;
        for (auto&& ind : remain_snps) {
            if (ind == start || ind - end == 1)
                end = ind; // try to perform the copy as a block
            else
            {
                std::copy(m_existed_snps.begin() + start,
                          m_existed_snps.begin() + end + 1, last);
                last += end + 1 - start;
                start = ind;
                end   = ind;
            }
        }
        if (!remain_snps.empty()) {
            std::copy(m_existed_snps.begin() + start,
                      m_existed_snps.begin() + end + 1, last);
            last += end + 1 - start;
        }
        m_existed_snps.erase(last, m_existed_snps.end());
    }
    m_existed_snps_index.clear();
    // we don't need the index any more because we don't need to match SNPs
    // anymore
    fprintf(stderr, "Number of SNPs after clumping : %zu\n\n",
            m_existed_snps.size());
}

void Genotype::perform_clump(size_t& core_genotype_index, int& begin_index,
                             int current_index, bool require_clump)
{
    if (core_genotype_index >= m_existed_snps.size()) return;
    /**
     * Previous use of SNP + genotype vector are way too complicated
     * Now do everything with m_existed_snps
     */
    int    core_chr       = m_existed_snps[core_genotype_index].chr();
    size_t core_loc       = m_existed_snps[core_genotype_index].loc();
    size_t infinite_guard = 0; // guard against infinite while loop
    size_t max_possible   = current_index - core_genotype_index;
    int    next_chr       = (current_index >= m_existed_snps.size())
                       ? m_existed_snps.back().chr() + 1
                       : m_existed_snps[current_index].chr();
    int next_loc = (current_index >= m_existed_snps.size())
                       ? m_existed_snps.back().loc() + 2 * clump_info.distance
                       : m_existed_snps[current_index].loc();
    while (require_clump
           && (core_chr != next_chr
               || (next_loc - core_loc) > clump_info.distance))
    { // as long as we still need to perform clumping
        clump_thread(core_genotype_index, begin_index, current_index);
        require_clump = false;
        core_genotype_index++;
        for (; core_genotype_index < current_index; ++core_genotype_index) {
            if (m_existed_snps[core_genotype_index].p_value()
                < clump_info.p_value)
            {
                core_chr      = m_existed_snps[core_genotype_index].chr();
                core_loc      = m_existed_snps[core_genotype_index].loc();
                require_clump = true;
                break;
            }
        }
        // New core found, need to clean things up a bit
        if (require_clump) {
            for (size_t remover = begin_index; remover < core_genotype_index;
                 ++remover)
            {
                if (core_loc - m_existed_snps[remover].loc()
                    > clump_info.distance)
                {
                    begin_index = remover + 1;
                    m_existed_snps[remover].clean_clump();
                }
                else
                {
                    begin_index = remover;
                    break;
                }
            }
        }
        infinite_guard++;
        if (infinite_guard > max_possible)
            throw std::logic_error(
                "ERROR: While loop run longer than expected");
    }
    // for this to be true, the require clump should be false or the core_snp is
    // now within reach of the new snp
    if (core_chr != next_chr)
    { // next_chr is not within m_genotype and clump_index yet
        for (size_t i = begin_index; i < current_index; ++i) {
            m_existed_snps[i].clean_clump();
        }
        begin_index = current_index;
    }
    else if (!require_clump)
    {
        // remove anything that is too far ahead
        for (size_t i = begin_index; i < current_index; ++i) {
            if (next_loc - m_existed_snps[i].loc() > clump_info.distance) {
                begin_index = i + 1;
                m_existed_snps[i].clean_clump();
            }
            else
            {
                begin_index = i;
                break;
            }
        }
    }
}

void Genotype::clump_thread(const size_t c_core_genotype_index,
                            const size_t c_begin_index,
                            const size_t c_current_index)
{

    uint32_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
    size_t   wind_size    = c_current_index - c_begin_index;
    if (wind_size <= 1) return; // nothing to do
    uint32_t tot1[6];
    tot1[0] = popcount_longs(m_existed_snps[c_core_genotype_index].clump_geno(),
                             founder_ctv3);
    tot1[1] = popcount_longs(
        &(m_existed_snps[c_core_genotype_index].clump_geno()[founder_ctv3]),
        founder_ctv3);
    tot1[2] = popcount_longs(
        &(m_existed_snps[c_core_genotype_index].clump_geno()[2 * founder_ctv3]),
        founder_ctv3);
    // in theory, std::fill, std::copy should be safer than these mem thing
    // but as a safety guard from my stupidity, let's just follow plink
    /*
            if (is_x) {
                    memcpy(geno_male, geno1, founder_ctsplit *
       sizeof(intptr_t)); bitvec_and(m_sex_male, founder_ctv3, geno_male);
                    tot1[3] = popcount_longs(geno_male, founder_ctv3);
                    bitvec_and(m_sex_male, founder_ctv3,
       &(geno_male[founder_ctv3])); tot1[4] =
       popcount_longs(&(geno_male[founder_ctv3]), founder_ctv3);
                    bitvec_and(m_sex_male, founder_ctv3, &(geno_male[2 *
       founder_ctv3])); tot1[5] = popcount_longs(&(geno_male[2 * founder_ctv3]),
       founder_ctv3);
            }
     */

    std::vector<std::thread> thread_store;
    if ((wind_size - 1) < m_thread) {
        for (size_t i_snp = c_begin_index; i_snp < c_current_index; ++i_snp) {
            if (i_snp == c_core_genotype_index) continue;
            thread_store.push_back(std::thread(
                &Genotype::compute_clump, this, c_core_genotype_index, i_snp,
                i_snp + 1,
                m_existed_snps[c_core_genotype_index].clump_missing(),
                std::ref(tot1)));
        }
    }
    else
    {
        int num_snp_per_thread =
            (int) (wind_size) / (int) m_thread; // round down
        int remain    = (int) (wind_size) % (int) m_thread;
        int cur_start = c_begin_index;
        int cur_end   = num_snp_per_thread + cur_start;
        for (size_t i_thread = 0; i_thread < m_thread; ++i_thread) {
            thread_store.push_back(std::thread(
                &Genotype::compute_clump, this, c_core_genotype_index,
                cur_start, cur_end + (remain > 0),
                m_existed_snps[c_core_genotype_index].clump_missing(),
                std::ref(tot1)));

            cur_start = cur_end + (remain > 0);
            cur_end += num_snp_per_thread + (remain > 0);
            if (cur_end > c_current_index) cur_end = c_current_index;
            remain--;
        }
    }
    for (auto&& thread_runner : thread_store) thread_runner.join();
    thread_store.clear();
    // delete [] geno_male;
}

void Genotype::compute_clump(size_t core_genotype_index, size_t i_start,
                             size_t i_end, bool nm_fixed, uint32_t* tot1)
{
    uint32_t counts[18];
    uint32_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
    double   freq11;
    double   freq11_expected;
    double   freq1x;
    double   freq2x;
    double   freqx1;
    double   freqx2;
    double   dxx;
    double   r2          = 0.0;
    bool     zmiss2      = false;
    size_t   max_size    = m_existed_snps.size();
    double   ref_p_value = m_existed_snps[core_genotype_index].p_value();
    int      ref_loc     = m_existed_snps[core_genotype_index].loc();
    auto&&   ref_geno    = m_existed_snps[core_genotype_index].clump_geno();
    std::vector<double> r2_store;
    std::vector<size_t>
        target_index_store; // index we want to push into the current index
    for (size_t i_snp = i_start; i_snp < i_end && i_snp < max_size; ++i_snp) {
        zmiss2 = false;
        if (i_snp != core_genotype_index
            && ((m_existed_snps[i_snp].p_value() > ref_p_value)
                || (m_existed_snps[i_snp].p_value() == ref_p_value
                    && (m_existed_snps[i_snp].loc() > ref_loc))))
        {
            // if the target is not as significant as the reference SNP or if it
            // is the same significance but with appear later in the genome

            auto&&    target_geno = m_existed_snps[i_snp].clump_geno();
            uintptr_t uiptr[3];
            uiptr[0] = popcount_longs(target_geno, founder_ctv3);
            uiptr[1] =
                popcount_longs(&(target_geno[founder_ctv3]), founder_ctv3);
            uiptr[2] =
                popcount_longs(&(target_geno[2 * founder_ctv3]), founder_ctv3);
            if (m_existed_snps[i_snp].clump_missing()) {
                zmiss2 = true;
            }

            if (nm_fixed) {
                two_locus_count_table_zmiss1(ref_geno, target_geno, counts,
                                             founder_ctv3, zmiss2);
                if (zmiss2) {
                    counts[2] = tot1[0] - counts[0] - counts[1];
                    counts[5] = tot1[1] - counts[3] - counts[4];
                }
                counts[6] = uiptr[0] - counts[0] - counts[3];
                counts[7] = uiptr[1] - counts[1] - counts[4];
                counts[8] = uiptr[2] - counts[2] - counts[5];
            }
            else
            {
                two_locus_count_table(ref_geno, target_geno, counts,
                                      founder_ctv3, zmiss2);
                if (zmiss2) {
                    counts[2] = tot1[0] - counts[0] - counts[1];
                    counts[5] = tot1[1] - counts[3] - counts[4];
                    counts[8] = tot1[2] - counts[6] - counts[7];
                }
            }

            /*
            // good thing is that the x1 and x2 must always be the same
            if (is_x1 || is_x2) {
                    two_locus_count_table(geno_male, ulptr, &(counts[9]),
            founder_ctv3, zmiss2); if (zmiss2) { counts[11] = tot1[3] -
            counts[9] - counts[10]; counts[14] = tot1[4] - counts[12] -
            counts[13]; counts[17] = tot1[5] - counts[15] - counts[16];
                    }
            }
             */
            // below, the false are basically is_x1 is_x2

            if (em_phase_hethet_nobase(counts, false, false, &freq1x, &freq2x,
                                       &freqx1, &freqx2, &freq11))
            {
                r2 = -1;
            }
            else
            {
                freq11_expected = freqx1 * freq1x;
                dxx             = freq11 - freq11_expected;
                if (fabs(dxx) < SMALL_EPSILON) {
                    r2 = 0.0;
                }
                else
                {
                    r2 = dxx * dxx / (freq11_expected * freq2x * freqx2);
                }
            }
            if (r2 >= clump_info.r2) {
                target_index_store.push_back(i_snp);
                r2_store.push_back(r2);
            }
        }
    }

    Genotype::clump_mtx.lock();
    m_existed_snps[core_genotype_index].add_clump(target_index_store);
    m_existed_snps[core_genotype_index].add_clump_r2(r2_store);
    Genotype::clump_mtx.unlock();
}

bool Genotype::prepare_prsice()
{
    if (m_existed_snps.size() == 0) return false;
    std::sort(begin(m_existed_snps), end(m_existed_snps),
              [](SNP const& t1, SNP const& t2) {
                  if (t1.category() == t2.category()) {
                      if (t1.file_name().compare(t2.file_name()) == 0) {
                          return t1.snp_id() < t2.snp_id();
                      }
                      else
                          return t1.file_name().compare(t2.file_name()) < 0;
                  }
                  else
                      return t1.category() < t2.category();
              });
    return true;
}

bool Genotype::get_score(misc::vec2d<Sample_lite>& prs_score, int& cur_index,
                         int& cur_category, double& cur_threshold,
                         std::vector<size_t>& num_snp_included)
{
    if (m_existed_snps.size() == 0 || cur_index == m_existed_snps.size())
        return false;
    int  end_index = 0;
    bool ended     = false;
    if (cur_index == -1) // first run
    {
        cur_index    = 0;
        cur_category = m_existed_snps[cur_index].category();
    }
    cur_threshold = m_existed_snps[cur_index].get_threshold();
    // existed snp should be sorted such that the SNPs should be
    // access sequentially
    for (size_t i = cur_index; i < m_existed_snps.size(); ++i) {
        if (m_existed_snps[i].category() != cur_category) {
            end_index = i;
            ended     = true;
            break;
        }
        //		// Use as part of the output
        for (size_t i_region = 0; i_region < m_region_size; ++i_region) {
            if (m_existed_snps[i].in(i_region)) num_snp_included[i_region]++;
        }
    }
    std::ofstream debug;
    debug.open("DEBUG");
    for (size_t i = 0; i < num_snp_included.size(); ++i) {
        debug << i << "\t" << num_snp_included[i]
              << std::endl; // need to match with region info
    }
    debug.close();
    // exit(0);
    if (!ended) {
        end_index    = m_existed_snps.size();
        cur_category = m_existed_snps.back().category();
    }
    else
        cur_category = m_existed_snps[end_index].category();
    read_score(prs_score, cur_index, end_index);

    cur_index = end_index;
    return true;
}

void Genotype::print_snp(std::string& output, double threshold)
{
    std::ofstream snp_out;
    snp_out.open(output);
    if (!snp_out.is_open()) {
        std::string error_message =
            "ERROR: Cannot open file: " + output + " to write";
        throw std::runtime_error(error_message);
    }
    for (auto&& snp : m_existed_snps) {
        if (snp.get_threshold() <= threshold) {
            snp_out << snp.rs() << "\t" << snp.chr() << "\t" << snp.loc()
                    << "\t" << snp.p_value() << std::endl;
        }
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
    uint32_t sol_end_idx   = 1;
    double   solutions[3];
    double   twice_tot_recip;
    double   half_hethet_share;
    double   freq11;
    double   freq12;
    double   freq21;
    double   freq22;
    double   prod_1122;
    double   prod_1221;
    double   incr_1122;
    double   best_sol;
    double   best_lnlike;
    double   cur_lnlike;
    double   freq1x;
    double   freq2x;
    double   freqx1;
    double   freqx2;
    double   lbound;
    double   dxx;
    uint32_t cur_sol_idx;
    // shouldn't have to worry about subtractive cancellation problems here
    if (twice_tot == 0.0) {
        return 1;
    }
    twice_tot_recip   = 1.0 / twice_tot;
    freq11            = known11 * twice_tot_recip;
    freq12            = known12 * twice_tot_recip;
    freq21            = known21 * twice_tot_recip;
    freq22            = known22 * twice_tot_recip;
    prod_1122         = freq11 * freq22;
    prod_1221         = freq12 * freq21;
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
                sol_end_idx   = 2;
                solutions[0]  = 0;
                solutions[1]  = half_hethet_share;
            }
            else
            {
                if (solutions[sol_start_idx] < 0) {
                    solutions[sol_start_idx] = 0;
                }
                // checking here
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
                sol_end_idx  = 3;
                solutions[1] = (half_hethet_share + freq21 - freq22) * 0.5;
                solutions[2] = half_hethet_share;
            }
            else
            {
                sol_end_idx  = 2;
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
                incr_1122  = solutions[cur_sol_idx];
                cur_lnlike = calc_lnlike(known11, known12, known21, known22,
                                         center_ct_d, freq11, freq12, freq21,
                                         freq22, half_hethet_share, incr_1122);
                if (cur_lnlike > best_lnlike) {
                    cur_lnlike = best_lnlike;
                    best_sol   = incr_1122;
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
                                   uint32_t   unfiltered_sample_ct,
                                   uintptr_t* casebuf, uint32_t case_ctv,
                                   uint32_t ctrl_ctv, uint32_t do_reverse,
                                   uint32_t   is_case_only,
                                   uintptr_t* nm_info_ptr)
{
    uintptr_t* rawbuf_end = &(rawbuf[unfiltered_sample_ct / BITCT2]);
    uintptr_t* ctrlbuf    = &(casebuf[3 * case_ctv]);
    uintptr_t  case_words[4];
    uintptr_t  ctrl_words[4];
    uint32_t   case_rem       = 0;
    uint32_t   ctrl_rem       = 0;
    uint32_t   read_shift_max = BITCT2;
    uint32_t   sample_uidx    = 0;
    uint32_t   offset0_case   = do_reverse * 2 * case_ctv;
    uint32_t   offset2_case   = (1 - do_reverse) * 2 * case_ctv;
    uint32_t   offset0_ctrl   = do_reverse * 2 * ctrl_ctv;
    uint32_t   offset2_ctrl   = (1 - do_reverse) * 2 * ctrl_ctv;
    uint32_t   read_shift;
    uintptr_t  read_word;
    uintptr_t  ulii;

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
                    casebuf[case_ctv]     = case_words[2];
                    casebuf[offset2_case] = case_words[3];
                    casebuf++;
                    case_words[0] = 0;
                    case_words[2] = 0;
                    case_words[3] = 0;
                    case_rem      = 0;
                }
                read_word >>= 2;
            }
        }
        if (sample_uidx == unfiltered_sample_ct) {
            if (case_rem) {
                casebuf[offset0_case] = case_words[0];
                casebuf[case_ctv]     = case_words[2];
                casebuf[offset2_case] = case_words[3];
            }
            if (ctrl_rem) {
                ctrlbuf[offset0_ctrl] = ctrl_words[0];
                ctrlbuf[ctrl_ctv]     = ctrl_words[2];
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
        uii           = counts_3x3[1];
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
    lptr2         = &(lptr2[sample_ctv3]);
    counts_3x3[1] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
    counts_3x3[4] =
        popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
    counts_3x3[7] =
        popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
    if (!is_zmiss2) {
        lptr2         = &(lptr2[sample_ctv3]);
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
                                            uint32_t  sample_ctv3,
                                            uint32_t  is_zmiss2)
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
    lptr1         = &(lptr1[sample_ctv3]);
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
    __m128i*      vec20;
    __m128i*      vec21;
    __m128i*      vec22;
    __m128i*      vend1;
    __m128i       loader1;
    __m128i       loader20;
    __m128i       loader21;
    __m128i       loader22;
    __m128i       count10;
    __m128i       count11;
    __m128i       count12;
    __m128i       count20;
    __m128i       count21;
    __m128i       count22;
    __univec      acc0;
    __univec      acc1;
    __univec      acc2;
    uint32_t      ct;
    uint32_t      ct2;
    while (iter_ct--) {
        ct    = sample_ctv6;
        vec20 = vec2;
        vec21 = &(vec20[sample_ctv6]);
        vec22 = &(vec20[2 * sample_ctv6]);
        while (ct >= 30) {
            ct -= 30;
            vend1   = &(vec1[30]);
            acc0.vi = _mm_setzero_si128();
            acc1.vi = _mm_setzero_si128();
            acc2.vi = _mm_setzero_si128();
            do
            {
            two_locus_3x3_tablev_outer:
                loader1  = *vec1++;
                loader20 = *vec20++;
                loader21 = *vec21++;
                loader22 = *vec22++;
                count10  = _mm_and_si128(loader1, loader20);
                count11  = _mm_and_si128(loader1, loader21);
                count12  = _mm_and_si128(loader1, loader22);
                count10  = _mm_sub_epi64(
                    count10, _mm_and_si128(_mm_srli_epi64(count10, 1), m1));
                count11 = _mm_sub_epi64(
                    count11, _mm_and_si128(_mm_srli_epi64(count11, 1), m1));
                count12 = _mm_sub_epi64(
                    count12, _mm_and_si128(_mm_srli_epi64(count12, 1), m1));
            two_locus_3x3_tablev_two_left:
                // unlike the zmiss variant, this apparently does not suffer
                // from enough register spill to justify shrinking the inner
                // loop
                loader1  = *vec1++;
                loader20 = *vec20++;
                loader21 = *vec21++;
                loader22 = *vec22++;
                count20  = _mm_and_si128(loader1, loader20);
                count21  = _mm_and_si128(loader1, loader21);
                count22  = _mm_and_si128(loader1, loader22);
                count20  = _mm_sub_epi64(
                    count20, _mm_and_si128(_mm_srli_epi64(count20, 1), m1));
                count21 = _mm_sub_epi64(
                    count21, _mm_and_si128(_mm_srli_epi64(count21, 1), m1));
                count22 = _mm_sub_epi64(
                    count22, _mm_and_si128(_mm_srli_epi64(count22, 1), m1));
            two_locus_3x3_tablev_one_left:
                loader1  = *vec1++;
                loader20 = *vec20++;
                loader21 = _mm_and_si128(loader1, loader20); // half1
                loader22 =
                    _mm_and_si128(_mm_srli_epi64(loader21, 1), m1); // half2
                count10  = _mm_add_epi64(count10, _mm_and_si128(loader21, m1));
                count20  = _mm_add_epi64(count20, loader22);
                loader20 = *vec21++;
                loader21 = _mm_and_si128(loader1, loader20);
                loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
                count11  = _mm_add_epi64(count11, _mm_and_si128(loader21, m1));
                count21  = _mm_add_epi64(count21, loader22);
                loader20 = *vec22++;
                loader21 = _mm_and_si128(loader1, loader20);
                loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
                count12  = _mm_add_epi64(count12, _mm_and_si128(loader21, m1));
                count22  = _mm_add_epi64(count22, loader22);

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
            vend1   = &(vec1[ct]);
            ct2     = ct % 3;
            acc0.vi = _mm_setzero_si128();
            acc1.vi = _mm_setzero_si128();
            acc2.vi = _mm_setzero_si128();
            ct      = 0;
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
