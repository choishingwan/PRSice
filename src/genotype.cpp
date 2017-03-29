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
		if(no_x){
			m_xymt_codes[X_OFFSET] = -1;
			clear_bit(num_auto + 1, m_haploid_mask);
		}
		if(no_y)
		{
			m_xymt_codes[Y_OFFSET] = -1;
			clear_bit(num_auto + 2, m_haploid_mask);
		}
		if(no_xy)
		{
			m_xymt_codes[XY_OFFSET] = -1;
		}
		if(no_mt)
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
		bool no_x, bool no_y, bool no_xy, bool no_mt, const size_t thread, bool verbose) {
	// TODO Auto-generated constructor stub
	init_chr(num_auto, no_x, no_y, no_xy, no_mt);
	// obtain files
	set_genotype_files(prefix);
	load_sample();
	m_existed_snps=load_snps();
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
	if(m_haploid_mask != nullptr) delete [] m_haploid_mask;
}



double Genotype::update_existed( Genotype &reference)
{
	int miss_match = 0;
	int matched = 0;
	for(auto &&snp : m_existed_snps)
	{
		auto target =reference.m_existed_snps_index.find(snp.get_rs());
		if(target==reference.m_existed_snps_index.end())
		{
			snp.not_required();
		}
		else
		{
			matched ++;
			//soft check, no biggy if they are not the same
			miss_match+=(snp==m_existed_snps[target->second])? 0 : 1;
		}
	}
	return (matched==0)? -1 : (double)miss_match/matched;
}

double Genotype::update_existed(const std::unordered_map<std::string, int> &ref_index,
		const std::vector<SNP> &reference)
{
	int miss_match = 0;
	int matched = 0;
	for(auto &&snp : m_existed_snps)
	{
		auto target = ref_index.find(snp.get_rs());
		if(target==ref_index.end())
		{
			snp.not_required();
		}
		else
		{
			matched ++;
			//soft check, no biggy if they are not the same
			miss_match+=(snp==m_existed_snps[target->second])? 0 : 1;
		}
	}
	return (matched==0)? -1 : (double)miss_match/matched;
}

void Genotype::reset_existed()
{
	for(auto &&snp : m_existed_snps)
	{
		snp.required();
	}
}


void Genotype::read_snps(const Commander &commander, const Region &region)
{
	// want to also sort out the partitioning here
	const std::string input = c_commander.base_name()
	const bool beta = c_commander.beta();
	if (beta && c_commander.statistic().compare("OR") == 0)
		fprintf(stderr, "WARNING: OR detected but user suggest the input is beta!\n");
	std::vector<int> index = SNP::get_index(c_commander, input); // more appropriate for commander
	// now coordinates obtained from target file instead. Coordinate information
	// in base file only use for validation

		// Open the file
		std::ifstream snp_file;
		snp_file.open(input.c_str());
		if (!snp_file.is_open())
		{
			std::string error_message = "Cannot open base file: " + input;
			throw std::runtime_error(error_message);
		}
		// Some QC countss
		int num_duplicated = 0;
		size_t num_stat_not_convertible = 0;
		size_t num_p_not_convertible = 0;
		size_t num_indel = 0;
		size_t num_se_not_convertible = 0;
		size_t num_exclude = 0;
		int max_index = index[+SNP_Index::MAX];
		std::string line;
		// remove header if index is not provided
		if (!c_commander.index()) std::getline(snp_file, line);
		bool read_error = false;
		bool not_converted = false;
		bool exclude = false;
		bool chr_error =false;
		// category related stuff
		double threshold = (c_commander.fastscore())? c_commander.get_bar_upper() : c_commander.get_upper();
		threshold = (c_commander.full())? 1.0 : threshold;
		std::vector < std::string > token;
		// Actual reading the file, will do a bunch of QC
		while (std::getline(snp_file, line))
		{
			misc::trim(line);
			if (!line.empty())
			{
				not_converted = false;
				exclude = false;
				token = misc::split(line);
				if (token.size() <= max_index)
					throw std::runtime_error("More index than column in data");
				else
				{
					std::string rs_id = token[index[+SNP_Index::RS]];
					std::string chr = (index[+SNP_Index::CHR] >= 0)? token[index[+SNP_Index::CHR]] : "";
					int32_t chr_code = get_chrom_code_raw(chr.c_str());
					if (((const uint32_t)chr_code) > PLINK::g_max_code) {
						if (chr_code != -1) {
							if (chr_code >= MAX_POSSIBLE_CHROM) {
								chr_code= PLINK::g_xymt_codes[chr_code - MAX_POSSIBLE_CHROM];
							}
							else
							{
								std::string error_message ="ERROR: Cannot parse chromosome code: " + chr;
								throw std::runtime_error(error_message);
							}
						}
					}
					if(chr_code == (uint32_t)PLINK::g_xymt_codes[X_OFFSET] ||
							chr_code == (uint32_t)PLINK::g_xymt_codes[Y_OFFSET])
					{
						exclude= true;
						if(!chr_error)
						{
							fprintf(stderr, "WARNING: Sex chromosome detected. They will be ignored\n");
							chr_error = true;
						}
					}
					std::string ref_allele = (index[+SNP_Index::REF] >= 0) ? token[index[+SNP_Index::REF]] : "";
					std::string alt_allele = (index[+SNP_Index::ALT] >= 0) ? token[index[+SNP_Index::ALT]] : "";
					double pvalue = 0.0;
					if (index[+SNP_Index::P] >= 0)
					{
						try {
							// obtain the p-value and calculate the corresponding threshold & category
							pvalue = misc::convert<double>(
									token[index[+SNP_Index::P]]);
							if (pvalue < 0.0 || pvalue > 1.0)
							{
								read_error = true;
								fprintf(stderr, "ERROR: %s's p-value is %f\n", rs_id.c_str(), pvalue);
							}
							else if (pvalue > threshold)
							{
								exclude = true;
								num_exclude++;
							}
						} catch (const std::runtime_error &error) {
							num_p_not_convertible++;
							not_converted = true;
						}
					}
					double stat = 0.0;
					if (index[+SNP_Index::STAT] >= 0)
					{
						//Obtain the test statistic
						try {
							stat = misc::convert<double>( token[index[+SNP_Index::STAT]]);
							if (!beta) stat = log(stat);
						} catch (const std::runtime_error& error) {
							num_stat_not_convertible++;
							not_converted = true;
						}
					}
					double se = 0.0;
					if (index[+SNP_Index::SE] >= 0)
					{
						// obtain the standard error (though it is currently useless)
						try {
							se = misc::convert<double>( token[index[+SNP_Index::SE]]);
						} catch (const std::runtime_error &error) {
							num_se_not_convertible++;
						}
					}
					int loc = -1;
					if (index[+SNP_Index::BP] >= 0)
					{
						// obtain the SNP coordinate
						try {
							int temp = misc::convert<int>( token[index[+SNP_Index::BP]].c_str());
							if (temp < 0)
							{
								read_error = true;
								fprintf(stderr, "ERROR: %s has negative loci\n", rs_id.c_str());
							}
							else loc = temp;
						} catch (const std::runtime_error &error) { }
					}
					if (!SNP::valid_snp(ref_allele))
					{
						num_indel++;
					}
					else if (!alt_allele.empty() && !SNP::valid_snp(alt_allele)) num_indel++;
					else if(!alt_allele.empty() && SNP::ambiguous(ref_allele, alt_allele)){
						num_exclude++;
					}
					else if (!not_converted && !exclude)
					{
						m_snp_list.push_back( new SNP(rs_id, chr, loc, ref_allele, alt_allele, stat, se, pvalue));
					} else {
	//		    			We skip any SNPs with non-convertible stat and p-value as we don't know how to
	//		    			handle them. Most likely those will be NA, which should be ignored anyway
					}
				}
			}
			if (read_error) throw std::runtime_error( "Please check if you have the correct input");
		}
		snp_file.close();

		m_snp_list.sort();

		// The problem of this filtering is that if the same SNP is duplicated but with different
		// p-value etc, they will not be removed
		size_t before = m_snp_list.size();
		m_snp_list.erase(std::unique(m_snp_list.begin(), m_snp_list.end()), m_snp_list.end());
		size_t after = m_snp_list.size();

		std::unordered_set<std::string> unique_chr;
		std::unordered_set<std::string> unique_rsid;

		for (size_t i_snp = 0; i_snp < m_snp_list.size(); ++i_snp)
		{
			m_snp_index[m_snp_list[i_snp].get_rs_id()] = i_snp;
			std::string cur_chr = m_snp_list[i_snp].get_chr();
			if (unique_chr.find(cur_chr) == unique_chr.end())
			{
				unique_chr.insert(cur_chr);
				m_chr_list.push_back(cur_chr);
			}
			std::string rs = m_snp_list[i_snp].get_rs_id();
			if(unique_rsid.find(rs)==unique_rsid.end())
			{
				unique_rsid.insert(rs);
			}
			else
			{
				std::string error_message = "WARNING: Duplicated SNP ID: " + rs;
				throw std::runtime_error(error_message);
			}
			if (index[+SNP_Index::CHR] >= 0 && index[+SNP_Index::BP] >= 0) // if we have chr and bp information
			{
				m_snp_list[i_snp].set_flag( region.check(cur_chr, m_snp_list[i_snp].get_loc()));
			}
		}

		num_duplicated = (int) before - (int) after;
		if (num_indel != 0) fprintf(stderr, "Number of invalid SNPs    : %zu\n", num_indel);
		if (num_exclude != 0) fprintf(stderr, "Number of SNPs excluded: %zu\n", num_exclude);
		if (num_duplicated != 0) fprintf(stderr, "Number of duplicated SNPs : %d\n", num_duplicated);
		if (num_stat_not_convertible != 0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_convertible);
		if (num_p_not_convertible != 0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_convertible);
		if (num_se_not_convertible != 0) fprintf(stderr, "Failed to convert %zu SE\n", num_se_not_convertible);
		fprintf(stderr, "Final Number of SNPs from base  : %zu\n", m_snp_list.size());
		PLINK::set_chromosome(m_chr_list); // update the chromosome information for PLINK
	}

}
