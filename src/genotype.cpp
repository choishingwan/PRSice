/*
 * genotype.cpp
 *
 *  Created on: 27 Mar 2017
 *      Author: shingwanchoi
 */

#include "genotype.hpp"

void Genotype::init_chr(int num_auto, bool no_x, bool no_y, bool no_xy, bool no_mt)
{
	// this initialize haploid mask as the maximum possible number
	chrom_info.haploid_mask = new uintptr_t[CHROM_MASK_WORDS];
	fill_ulong_zero(CHROM_MASK_WORDS, chrom_info.haploid_mask);

	if(num_auto < 0)
	{
		num_auto = -num_auto;
		chrom_info.autosome_ct = num_auto;
		chrom_info.xymt_codes[X_OFFSET] = -1;
		chrom_info.xymt_codes[Y_OFFSET] = -1;
		chrom_info.xymt_codes[XY_OFFSET] = -1;
		chrom_info.xymt_codes[MT_OFFSET] = -1;
		chrom_info.max_code = num_auto;
		fill_all_bits(((uint32_t)num_auto) + 1, chrom_info.haploid_mask);
	}
	else
	{
		chrom_info.autosome_ct = num_auto;
		chrom_info.xymt_codes[X_OFFSET] = num_auto+1;
		chrom_info.xymt_codes[Y_OFFSET] = num_auto+2;
		chrom_info.xymt_codes[XY_OFFSET] = num_auto+3;
		chrom_info.xymt_codes[MT_OFFSET] = num_auto+4;
		set_bit(num_auto + 1, chrom_info.haploid_mask);
		set_bit(num_auto + 2, chrom_info.haploid_mask);
		if(no_x){
			chrom_info.xymt_codes[X_OFFSET] = -1;
			clear_bit(num_auto + 1, chrom_info.haploid_mask);
		}
		if(no_y)
		{
			chrom_info.xymt_codes[Y_OFFSET] = -1;
			clear_bit(num_auto + 2, chrom_info.haploid_mask);
		}
		if(no_xy)
		{
			chrom_info.xymt_codes[XY_OFFSET] = -1;
		}
		if(no_mt)
		{
			chrom_info.xymt_codes[MT_OFFSET] = -1;
		}
		if (chrom_info.xymt_codes[MT_OFFSET] != -1) {
			chrom_info.max_code = num_auto + 4;
		} else if (chrom_info.xymt_codes[XY_OFFSET] != -1) {
			chrom_info.max_code = num_auto + 3;
		} else if (chrom_info.xymt_codes[Y_OFFSET] != -1) {
			chrom_info.max_code = num_auto + 2;
		} else if (chrom_info.xymt_codes[X_OFFSET] != -1) {
			chrom_info.max_code = num_auto + 1;
		} else {
			chrom_info.max_code = num_auto;
		}
	}
	chrom_info.chrom_mask = new uintptr_t[CHROM_MASK_WORDS];
	fill_ulong_zero(CHROM_MASK_WORDS, chrom_info.chrom_mask);
	fill_all_bits(chrom_info.autosome_ct + 1, chrom_info.chrom_mask);
	for (uint32_t xymt_idx = 0; xymt_idx < XYMT_OFFSET_CT; ++xymt_idx) {
		int32_t cur_code = chrom_info.xymt_codes[xymt_idx];
		if (cur_code != -1) {
			set_bit(chrom_info.xymt_codes[xymt_idx], chrom_info.chrom_mask);
		}
	}
	chrom_info.chrom_start.resize(chrom_info.max_code);// 1 extra for the info
}

void Genotype::set_genotype_files(std::string prefix)
{
	if(prefix.find("#")!=std::string::npos)
	{
		for(size_t chr = 1; chr < chrom_info.max_code; ++chr)
		{
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

Genotype::Genotype(std::string prefix, int num_auto,
		bool no_x, bool no_y, bool no_xy, bool no_mt, const size_t thread, bool verbose)
{
	m_xymt_codes.resize(XYMT_OFFSET_CT);
	init_chr(num_auto, no_x, no_y, no_xy, no_mt);
	set_genotype_files(prefix);
	m_sample_names = load_samples();
	m_existed_snps = load_snps();
	if(verbose)
	{
		fprintf(stderr, "%zu people (%zu males, %zu females) included\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
		if(m_num_ambig!=0) fprintf(stderr, "%zu ambiguous variants excluded\n", m_num_ambig);
		fprintf(stderr, "%zu variants included\n", m_marker_ct);
	}
}

void Genotype::load(std::vector<SNP> &snp_info, std::unordered_map<std::string, int> &snp_index,
		const Commander &c_commander, bool verbose)
{
	// load sample should always be straightforward
	load_sample();
	// now load snps, this is complicated because we want
	// 1. hh_exist (Most important!)
	// 2. filtering by mind if required (On all genes, before extract?)
	// 3. filtering by maf and geno if required
	// need to remember what SNPs causes the hh_exist.
	// maybe after filtering, hh_exist is no longer required
	// sequence is mind -> geno -> maf
	// let the user worries about the sequence of filtering
	load_snps(snp_info, snp_index, c_commander);
	if(verbose)
	{
		fprintf(stderr, "%zu people (%zu males, %zu females) loaded from .fam\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
		fprintf(stderr, "%zu variants included\n", m_marker_ct);
	}
}

Genotype::~Genotype() {
	// TODO Auto-generated destructor stub
	if(m_founder_info!=nullptr) delete [] m_founder_info;
	if(m_sex_male != nullptr) delete [] m_sex_male;
	if(m_sample_exclude != nullptr) delete [] m_sample_exclude;
	if(m_marker_exclude != nullptr) delete [] m_marker_exclude;
	if(chrom_info.haploid_mask != nullptr) delete [] chrom_info.haploid_mask;
	if(chrom_info.chrom_mask != nullptr) delete [] chrom_info.chrom_mask;
}

void Genotype::read_base(const Commander &c_commander, Region &region)
{
	// can assume region is of the same order as m_existed_snp
	const std::string input = c_commander.base_name();
	const bool beta = c_commander.beta();
	const bool fastscore = c_commander.fastscore();
	const bool full = c_commander.full();
	std::vector<int> index = c_commander.index(); // more appropriate for commander
	// now coordinates obtained from target file instead. Coordinate information
	// in base file only use for validation
	std::ifstream snp_file;
	snp_file.open(input.c_str());
	if(!snp_file.is_open())
	{
		std::string error_message = "ERROR: Cannot open base file: " +input;
		throw std::runtime_error(error_message);
	}
	int max_index = index[+BASE_INDEX::MAX];
	std::string line;
	if (!c_commander.has_index()) std::getline(snp_file, line);

	// category related stuff
	double threshold = (c_commander.fastscore())? c_commander.bar_upper() : c_commander.upper();
	double bound_start = c_commander.lower();
	double bound_end = c_commander.upper();
	double bound_inter = c_commander.inter();

	threshold = (full)? 1.0 : threshold;
	std::vector < std::string > token;

	bool exclude = false;
	bool hap_error = false;
	// Some QC countss
	size_t num_duplicated = 0;
	size_t num_excluded = 0;
	size_t num_not_found = 0;
	size_t num_mismatched = 0;
	size_t num_not_converted = 0; // this is for NA
	size_t num_negative_stat = 0;
	std::unordered_set<std::string> dup_index;
	std::vector<int> exist_index; // try to use this as quick search
	// Actual reading the file, will do a bunch of QC
	while (std::getline(snp_file, line))
	{
		misc::trim(line);
		if (line.empty()) continue;
		exclude = false;
		token = misc::split(line);
		if (token.size() <= max_index)
			throw std::runtime_error("More index than column in data");
		else
		{
			std::string rs_id = token[index[+BASE_INDEX::RS]];
			auto target = m_existed_snps_index.find(rs_id);
			if(target!=m_existed_snps_index.end() && dup_index.find(rs_id)==dup_index.end())
			{
				dup_index.insert(rs_id);
				auto &&cur_snp = m_existed_snps[target->second];
				if(!cur_snp.is_required()) num_excluded++;
				else
				{
					int32_t chr_code = -1;
					if (index[+BASE_INDEX::CHR] >= 0)
					{
						chr_code = get_chrom_code_raw(token[index[+BASE_INDEX::CHR]].c_str());
						if (((const uint32_t)chr_code) > m_max_code) {
							if (chr_code != -1) {
								if (chr_code >= MAX_POSSIBLE_CHROM) {
									chr_code= m_xymt_codes[chr_code - MAX_POSSIBLE_CHROM];
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
							fprintf(stderr, "WARNING: Currently not support haploid chromosome and sex chromosomes\n");
							exclude = true;
							num_excluded++;
						}
					}
					std::string ref_allele = (index[+BASE_INDEX::REF] >= 0) ? token[index[+BASE_INDEX::REF]] : "";
					std::string alt_allele = (index[+BASE_INDEX::ALT] >= 0) ? token[index[+BASE_INDEX::ALT]] : "";
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
						misc::convert<double>( token[index[+BASE_INDEX::P]]);
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
						num_not_converted = true;
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


					if(!alt_allele.empty() && SNP::ambiguous(ref_allele, alt_allele)){
						num_excluded++;
						exclude= true;
					}
					if(!exclude)
					{
						int category = -1;
						double pthres = 0.0;
						if (fastscore)
						{
							category = c_commander.get_category(pvalue);
							pthres = c_commander.get_threshold(category);
						}
						else
						{
							// calculate the threshold instead
							if (pvalue > bound_end && full)
							{
								category = std::ceil((bound_end + 0.1 - bound_start) / bound_inter);
								pthres = 1.0;
							}
							else
							{
								category = std::ceil((pvalue - bound_start) / bound_inter);
								category = (category < 0) ? 0 : category;
								pthres = category * bound_inter + bound_start;
							}
						}
						if(flipped) cur_snp.set_flipped();
						// ignore the SE as it currently serves no purpose
						exist_index.push_back(target->second);
						cur_snp.set_statistic(stat, 0.0, pvalue, category, pthres);
					}
				}
			}
			else if(dup_index.find(rs_id)!=dup_index.end())
			{
				num_duplicated++;
			}
			else
			{
				num_not_found++;
			}
		}
	}
	snp_file.close();

	if(!num_duplicated) fprintf(stderr, "%zu duplicated variant(s) in base file\n", num_duplicated);
	if(!num_excluded) fprintf(stderr, "%zu variant(s) excluded\n", num_excluded);
	if(!num_not_found) fprintf(stderr, "%zu variant(s) not found in target file\n", num_not_found);
	if(!num_mismatched) fprintf(stderr ,"%zu mismatched variant(s) excluded\n", num_mismatched);
	if(!num_not_converted) fprintf(stderr, "%zu NA statistic observed\n", num_not_converted);
	if(!num_negative_stat) fprintf(stderr, "%zu negative statistic observed. Please make sure it is really OR\n", num_negative_stat);
	fprintf(stderr, "%zu total SNPs included from base file\n");
	if(exist_index.size() != m_existed_snps.size())
	{ // only do this if we need to remove some SNPs
		// we assume exist_index doesn't have any duplicated index
		std::sort(exist_index.begin(), exist_index.end());
		int start = (exist_index.empty())? -1:exist_index.front();
		int end = start;
		std::vector<SNP>::iterator last = m_existed_snps.begin();;
		for(auto && ind : exist_index)
		{
			if(ind==start||ind-end==1) end=ind; // try to perform the copy as a block
			else{
				std::copy(m_existed_snps.begin()+start, m_existed_snps.begin()+end+1,last);
				last += end+1-start;
				start =ind;
				end = ind;
			}
		}
		if(!exist_index.empty())
		{
			std::copy(m_existed_snps.begin()+start, m_existed_snps.begin()+end+1, last);
			last+= end+1-start;
		}
		m_existed_snps.erase(last, m_existed_snps.end());
	}
	m_existed_snps_index.clear();
	// now m_existed_snps is ok and can be used directly
	size_t vector_index = 0;
	for(auto &&cur_snp : m_existed_snps) // should be in the correct order
	{
		m_existed_snps_index[cur_snp.rs()] = vector_index++;
		cur_snp.set_flag( region.check(cur_snp.chr(), cur_snp.loc()));
	}
	clump_info.p_value = c_commander.clump_p();
	clump_info.r2 =  c_commander.clump_r2();
	clump_info.proxy = c_commander.proxy();
	clump_info.use_proxy = c_commander.use_proxy();
	filter.filter_geno = c_commander.filter_geno();
	filter.filter_maf = c_commander.filter_maf();
	filter.filter_info = c_commander.filter_info();
	filter.geno = c_commander.geno();
	filter.info_score = c_commander.info_score();
	filter.maf = c_commander.maf();
}


void Genotype::clump(Genotype &reference)
{
	uintptr_t unfiltered_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_unfiltered_sample_ct);
	uintptr_t final_mask = get_final_mask(m_founder_ct);
	// use the old algorithm to speed up development time. Current focus should be
	// to get bgen working asap
	auto &&cur_snp = m_existed_snps.front();
	size_t snp_id_in_list = cur_snp.snp_id();
	std::string prev_file=cur_snp.file_name();
	int prev_chr= cur_snp.chr();
	size_t bp_of_core =cur_snp.loc();
	for(auto &&snp = m_existed_snps)
	{
		auto &&target = reference.m_existed_snps_index.find(snp.rs());
		if(target==reference.m_existed_snps_index.end()) continue; // only work on SNPs that are in both


	}
}



