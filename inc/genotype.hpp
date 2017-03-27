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
#include "plink_common.hpp"

class Genotype {
public:
	Genotype(std::string prefix, int num_auto=22, bool x=true, bool y=true, bool xy=true, bool mt=true,
			const size_t thread=1, bool verbose=false);
	virtual ~Genotype();
protected:
	void init_chr(int num_auto, bool x, bool y, bool xy, bool mt);
	virtual void read_genotype();
	uint32_t m_autosome_ct;
	int32_t m_xymt_codes[XYMT_OFFSET_CT];
	uint32_t m_max_code;
	uintptr_t* m_haploid_mask;
	//hh_exists
	/*
	uintptr_t* m_haploid_mask;
	uintptr_t* m_xymt_codes;
	uint32_t m_autosome_ct;
	uint32_t m_max_code;
	uint32_t m_hh_exists;
	*/
};

class BinaryPlink: public Genotype{
public:
	BinaryPlink(std::string prefix, const size_t thread=1, uint32_t species_code=SPECIES_HUMAN, bool verbose=false);
};

class Plink: public Genotype{
public:
	Plink(std::string prefix, const size_t thread=1, uint32_t species_code=SPECIES_HUMAN, bool verbose=false);
};

class GenomeFactory {
  public:
	 std::unique_ptr<BinaryPlink> createBinaryPlink(std::string prefix, const size_t thread=1,
			 uint32_t species_code=SPECIES_HUMAN, bool verbose=false)
	 {
		 return std::unique_ptr<BinaryPlink>(new BinaryPlink(prefix, thread, species_code, verbose));
	 };

	 std::unique_ptr<BinaryPlink> createPlink(std::string prefix, const size_t thread=1,
			 uint32_t species_code=SPECIES_HUMAN, bool verbose=false)
	{
		 return std::unique_ptr<BinaryPlink>(new BinaryPlink(prefix, thread, species_code, verbose));
	}

};

#endif /* SRC_GENOTYPE_HPP_ */
