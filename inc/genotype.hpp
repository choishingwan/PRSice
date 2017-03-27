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

#define SPECIES_HUMAN 0
#define SPECIES_COW 1
#define SPECIES_DOG 2
#define SPECIES_HORSE 3
#define SPECIES_MOUSE 4
#define SPECIES_RICE 5
#define SPECIES_SHEEP 6
#define SPECIES_UNKNOWN 7
#define SPECIES_DEFAULT SPECIES_HUMAN

class Genotype {
public:
	Genotype(std::string prefix, const size_t thread=1, uint32_t species_code=SPECIES_HUMAN, bool verbose=false);
	virtual ~Genotype();
protected:
	virtual void read_genotype();
	Chrom_info m_chrom_info;
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
