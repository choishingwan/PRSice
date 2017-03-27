/*
 * genotype.cpp
 *
 *  Created on: 27 Mar 2017
 *      Author: shingwanchoi
 */

#include "genotype.hpp"

void Genotype::init_chr(int num_auto, bool x, bool y, bool xy, bool mt)
{
	// this initialize haploid mask as the maximum possible number
	m_haploid_mask = new uintptr_t[CHROM_MASK_WORDS];
	fill_ulong_zero(CHROM_MASK_WORDS, m_haploid_mask);

	if(num_auto < 0)
	{
		num_auto = -num_auto;
		m_autosome_ct = num_auto;
		m_xymt_codes[X_OFFSET] = -1;
		m_xymt_codes[Y_OFFSET] = -1;
		m_xymt_codes[XY_OFFSET] = -1;
		m_xymt_codes[MT_OFFSET] = -1;
		m_max_code = num_auto;
		fill_all_bits(((uint32_t)num_auto) + 1, m_haploid_mask);
	}
	else
	{
		m_autosome_ct = num_auto;
		m_xymt_codes[X_OFFSET] = num_auto+1;
		m_xymt_codes[Y_OFFSET] = num_auto+2;
		m_xymt_codes[XY_OFFSET] = num_auto+3;
		m_xymt_codes[MT_OFFSET] = num_auto+4;
		set_bit(num_auto + 1, m_haploid_mask);
		set_bit(num_auto + 2, m_haploid_mask);
		if(!x){
			m_xymt_codes[X_OFFSET] = -1;
			clear_bit(num_auto + 1, m_haploid_mask);
		}
		if(!y)
		{
			m_xymt_codes[Y_OFFSET] = -1;
			clear_bit(num_auto + 2, m_haploid_mask);
		}
		if(!xy)
		{
			m_xymt_codes[XY_OFFSET] = -1;
		}
		if(!mt)
		{
			m_xymt_codes[MT_OFFSET] = -1;
		}
		if (m_xymt_codes[MT_OFFSET] != -1) {
			m_max_code = num_auto + 4;
		} else if (m_xymt_codes[XY_OFFSET] != -1) {
			m_max_code = num_auto + 3;
		} else if (m_xymt_codes[Y_OFFSET] != -1) {
			m_max_code = num_auto + 2;
		} else if (m_xymt_codes[X_OFFSET] != -1) {
			m_max_code = num_auto + 1;
		} else {
			m_max_code = num_auto;
		}
	}
}

void Genotype::set_genotype_files(std::string prefix)
{
	if(prefix.find("#")!=std::string::npos)
	{
		for(size_t chr = 1; chr < m_max_code; ++chr)
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
		bool x, bool y, bool xy, bool mt, const size_t thread, bool verbose) {
	// TODO Auto-generated constructor stub
	init_chr(num_auto, x, y, xy, mt);
	// obtain files
	set_genotype_files(prefix);
	load_sample();
	load_snps();
	if(verbose)
	{
		fprintf(stderr, "%zu people (%zu males, %zu females) loaded from .fam\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
		fprintf(stderr, "%zu variants included\n", m_marker_ct);
	}
}

Genotype::~Genotype() {
	// TODO Auto-generated destructor stub
	if(!m_founder_info==nullptr) delete [] m_founder_info;
	if(!m_sex_male == nullptr) delete [] m_sex_male;
	if(!m_sample_exclude == nullptr) delete [] m_sample_exclude;
	if(!m_marker_exclude == nullptr) delete [] m_marker_exclude;
	if(!m_haploid_mask == nullptr) delete [] m_haploid_mask;
}


