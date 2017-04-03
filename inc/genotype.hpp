/*
 * genotype.hpp
 *
 *  Created on: 27 Mar 2017
 *      Author: shingwanchoi
 */

#ifndef SRC_GENOTYPE_HPP_
#define SRC_GENOTYPE_HPP_

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cstring>
#include <cctype>
#include <algorithm>
#include "snp.hpp"
#include "commander.hpp"
#include "region.hpp"
#include "plink_common.hpp"
#include "misc.hpp"
#include "storage.hpp"

class Genotype {
public:
	Genotype(std::string prefix, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
			bool no_mt=false, const size_t thread=1, bool verbose=false);
	virtual ~Genotype();
	double update_existed( Genotype &reference);
	double update_existed(const std::unordered_map<std::string, int> &ref_index,
			const std::vector<SNP> &reference);

	inline bool existed (const std::string &rs) const
	{
		return m_existed_snps_index.find(rs)!=m_existed_snps_index.end();
	};

	inline bool matched (const SNP &target) const
	{
		if(existed(target.get_rs())){
			return m_existed_snps.at(m_existed_snps_index.at(target.get_rs()))==target;
		}
		return false;
	};
	void read_base(const Commander &commander, Region &region);

protected:

	void finalize_snps(Region &region, const int distance);

	void set_genotype_files(std::string prefix);
	std::vector<std::string> m_genotype_files;


	void init_chr(int num_auto, bool no_x, bool no_y, bool no_xy, bool no_mt);
	uint32_t m_autosome_ct;
	int32_t m_xymt_codes[XYMT_OFFSET_CT];
	uint32_t m_max_code;
	uintptr_t* m_haploid_mask;


	virtual void load_sample(){};
	uintptr_t* m_founder_info = nullptr;
	uintptr_t* m_sex_male = nullptr;
	uintptr_t* m_sample_exclude = nullptr;
	size_t m_num_male=0;
	size_t m_num_female=0;
	size_t m_num_ambig_sex=0;
	uintptr_t m_founder_ct = 0;
	uintptr_t m_unfiltered_sample_ct = 0;
	uintptr_t m_unfiltered_sample_ctl = 0;
	uintptr_t m_unfiltered_sample_ct4 = 0;

	// normally, the vector might go our of scope. Key is, use this as return value
	// then use c++11 which move instead of copy the variable if allowed
	virtual std::vector<SNP> load_snps(){ return std::vector<SNP>(0); };
	uintptr_t m_unfiltered_marker_ct = 0;
	uintptr_t m_unfiltered_marker_ctl = 0;
	uintptr_t m_marker_ct = 0;
	uintptr_t m_marker_exclude_ct = 0;
	uintptr_t* m_marker_exclude = nullptr;
	std::unordered_map<std::string, size_t> m_existed_snps_index;
	std::vector<SNP> m_existed_snps;
	std::unordered_map<std::string, int> m_chr_order;

	std::unordered_map<std::string, int> get_chr_order() const { return m_chr_order; };
	virtual void read_genotype(){};
	//hh_exists
	inline bool ambiguous(std::string ref_allele, std::string alt_allele)
	{
		return (ref_allele == "A" && alt_allele == "T")
				|| (ref_allele == "a" && alt_allele == "t")
				|| (ref_allele == "G" && alt_allele == "C")
				|| (ref_allele == "g" && alt_allele == "c");
	}
};

class Plink: public Genotype{
public:
	Plink(std::string prefix, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
			bool no_mt=false, const size_t thread=1, bool verbose=false);
private:
	void load_sample();
	std::vector<SNP> load_snps();
};


class BinaryPlink: public Genotype{
public:
	BinaryPlink(std::string prefix, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
			bool no_mt=false, const size_t thread=1, bool verbose=false);
private:
	uintptr_t m_bed_offset = 3;
	void load_sample();
	std::vector<SNP> load_snps();
	void load_bed();
	FILE* m_bedfile = nullptr;
};


class GenomeFactory {
  public:
	std::unique_ptr<Genotype> createGenotype(const Commander &commander, const std::string &prefix, bool verbose)
	{
		return std::unique_ptr<Genotype>(new BinaryPlink(prefix, commander.num_auto(),
				commander.no_x(), commander.no_y(), commander.no_xy(), commander.no_mt(),
				commander.thread(), verbose));
	}

};

#endif /* SRC_GENOTYPE_HPP_ */
