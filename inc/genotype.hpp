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
#include <cstdio>
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
			bool no_mt=false, const size_t thread=1);
	virtual ~Genotype();

	std::unordered_map<std::string, int> get_chr_order() const { return m_chr_order; };
	void clump(Genotype &reference);

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

	void filter_mind(double mind);
	void clump(Genotype &reference);

protected:
	struct{
		double r2;
		double proxy;
		double p_value;
		bool use_proxy;
	} clump_info;

	struct{
		double maf;
		double geno;
		double info_score;
		bool filter_maf;
		bool filter_geno;
		bool filter_info;
	} filter;

	struct sample_id{
		std::string FID;
		std::string IID;
	};

	void finalize_snps(Region &region, const int distance);

	void set_genotype_files(std::string prefix);
	std::vector<std::string> m_genotype_files;

	void init_chr(int num_auto, bool no_x, bool no_y, bool no_xy, bool no_mt);
	uint32_t m_autosome_ct;
	std::vector<int32_t> m_xymt_codes;
	uint32_t m_max_code;
	uintptr_t* m_haploid_mask;


	virtual std::vector<Sample> load_samples(){ return std::vector<Sample>(0); };
	uintptr_t m_unfiltered_sample_ct = 0;
	uintptr_t m_unfiltered_sample_ctl = 0;
	uintptr_t m_unfiltered_sample_ct4 = 0;
	std::vector<Sample> m_sample_names;

	void init_sample_vectors(){};
	uintptr_t* m_founder_info = nullptr;
	uintptr_t* m_sex_male = nullptr;
	uintptr_t* m_sample_exclude = nullptr;
	size_t m_num_male=0;
	size_t m_num_female=0;
	size_t m_num_ambig_sex=0;
	uintptr_t m_founder_ct = 0;

	virtual void load_snps(std::vector<SNP> &snp_info, std::unordered_map<std::string, int> &snp_index,
			const Commander &c_commander);
	uintptr_t m_unfiltered_marker_ct = 0;
	uintptr_t m_unfiltered_marker_ctl = 0;
	uintptr_t m_marker_ct = 0;
	uintptr_t m_marker_exclude_ct = 0;
	uintptr_t* m_marker_exclude = nullptr;
	std::unordered_map<std::string, size_t> m_existed_snps_index;
	std::vector<SNP> m_existed_snps;
	std::unordered_map<std::string, int> m_chr_order;
	uint32_t m_hh_exists; // might be a bit harsh, but should also read in maf when loading SNPs
	uint32_t m_num_ambig;
	struct{
		double r2;
		double proxy;
		double p_value;
		bool use_proxy;
	} clump_info;


	virtual void read_genotype(){};
	//hh_exists
	inline bool ambiguous(std::string ref_allele, std::string alt_allele)
	{
		return (ref_allele == "A" && alt_allele == "T")
				|| (ref_allele == "a" && alt_allele == "t")
				|| (ref_allele == "G" && alt_allele == "C")
				|| (ref_allele == "g" && alt_allele == "c");
	}

	void filter_mind(double mind);
};

/*
class Plink: public Genotype{
public:
	Plink(std::string prefix, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
			bool no_mt=false, const size_t thread=1);
private:
	void load_sample();
	std::vector<SNP> load_snps();
};
*/

class BinaryPlink: public Genotype{
public:
	BinaryPlink(std::string prefix, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
			bool no_mt=false, const size_t thread=1);
private:
	uintptr_t m_bed_offset = 3;
	std::vector<Sample> load_samples();
	std::vector<SNP> load_snps();
	void load_bed();
	FILE* m_bedfile = nullptr;
};


class GenomeFactory {
private:
	std::unordered_map<std::string, int> file_type {
		{ "bed", 0 },
		{ "ped", 1 },
		{ "bgen", 2}
	};
public:
	std::unique_ptr<Genotype> createGenotype(const Commander &commander, const std::string &prefix,
			const std::string &type)
	{
		fprintf(stderr, "Loading Genotype file: %s ", prefix.c_str());
		int code = (file_type.find(type)!=file_type.end())? file_type[type]: 0;
		switch(code)
		{
		case 1:
			/*
			return std::unique_ptr<Genotype>(new Plink(prefix, commander.num_auto(),
							commander.no_x(), commander.no_y(), commander.no_xy(), commander.no_mt(),
							commander.thread(), verbose));
							*/
		case 2:
		default:
		case 0:
			fprintf(stderr, "(bed)\n");
			return std::unique_ptr<Genotype>(new BinaryPlink(prefix, commander.num_auto(),
							commander.no_x(), commander.no_y(), commander.no_xy(), commander.no_mt(),
							commander.thread()));

		}
	}

};

#endif /* SRC_GENOTYPE_HPP_ */
