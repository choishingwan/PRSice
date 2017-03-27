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
#include "plink_common.hpp"

class Genotype {
public:
	Genotype(std::string prefix, int num_auto=22, bool x=true, bool y=true, bool xy=true, bool mt=true,
			const size_t thread=1, bool verbose=false);
	virtual ~Genotype();
protected:
	struct SNP{
		inline SNP(std::string chr, size_t loc, std::string ref, std::string alt, std::string file,
				size_t line, bool include): chr(chr), ref(ref), alt(alt), loc(loc), file(file),
						line(line), include(include);
		std::string chr;
		std::string ref;
		std::string alt;
		size_t loc;
		std::string file;
		size_t line;
		bool include;
	};

	void set_genotype_files(std::string prefix);
	std::vector<std::string> m_genotype_files;


	void init_chr(int num_auto, bool x, bool y, bool xy, bool mt);
	uint32_t m_autosome_ct;
	int32_t m_xymt_codes[XYMT_OFFSET_CT];
	uint32_t m_max_code;
	uintptr_t* m_haploid_mask;


	virtual void load_sample();
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

	virtual void load_snps();
	uintptr_t m_unfiltered_marker_ct = 0;
	uintptr_t m_unfiltered_marker_ctl = 0;
	uintptr_t m_marker_ct = 0;
	uintptr_t m_marker_exclude_ct = 0;
	uintptr_t* m_marker_exclude = nullptr;
	std::unordered_map<std::string, SNP> m_existed_snps;

	virtual void read_genotype();
	//hh_exists


};

class Plink: public Genotype{
public:
	Plink(std::string prefix, int num_auto=22, bool x=true, bool y=true, bool xy=true, bool mt=true,
			const size_t thread=1, bool verbose=false);
private:
	void load_sample();
	void load_snps();
};


class BinaryPlink: public Genotype{
public:
	BinaryPlink(std::string prefix, int num_auto=22, bool x=true, bool y=true, bool xy=true, bool mt=true,
			const size_t thread=1, bool verbose=false);
private:
	uintptr_t m_bed_offset = 3;
	void load_sample();
	void load_snps();
};


class GenomeFactory {
  public:
	 std::unique_ptr<BinaryPlink> createBinaryPlink(std::string prefix, int num_auto=22, bool x=true,
			 bool y=true, bool xy=true, bool mt=true, const size_t thread=1, bool verbose=false)
	 {
		 return std::unique_ptr<BinaryPlink>(new BinaryPlink(prefix, num_auto, x, y, xy, mt
				 thread, verbose));
	 };

	 std::unique_ptr<Plink> createPlink(std::string prefix, int num_auto=22, bool x=true,
			 bool y=true, bool xy=true, bool mt=true, const size_t thread=1, bool verbose=false)
	{
		 return std::unique_ptr<Plink>(new BinaryPlink(prefix, num_auto, x, y, xy, mt,
				 thread, verbose));
	}

};

#endif /* SRC_GENOTYPE_HPP_ */
