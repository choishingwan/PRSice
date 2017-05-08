/*
 * storage.hpp
 *
 *  Created on: 14 Oct 2016
 *      Author: shingwanchoi
 */

#ifndef PRSICE_INC_STORAGE_HPP_
#define PRSICE_INC_STORAGE_HPP_
#include <memory>
#include <cstdint>
#include <string>
	// From http://stackoverflow.com/a/12927952/1441789
	struct Sample{
		std::string FID;
		std::string IID;
		std::string pheno;
		double prs;
		int num_snp;
		bool included;
	};

	struct Sample_lite{
		double prs;
		int num_snp;
	};

    struct prsice_result{
		double threshold;
		double r2;
		double r2_adj;
		double coefficient;
		double p;
		double emp_p;
		int num_snp;
    };


	template< typename e >
	struct enumeration_traits;

	struct enumeration_trait_indexing {
		static constexpr bool does_index = true;
	};

	template< typename e >
	constexpr
	typename std::enable_if< enumeration_traits< e >::does_index,
    		typename std::underlying_type< e >::type >::type
	operator + ( e val )
    		{ return static_cast< typename std::underlying_type< e >::type >( val ); }

	template< typename e >
	typename std::enable_if< enumeration_traits< e >::does_index, e & >::type
	operator ++ ( e &val )
    		{ return val = static_cast< e >( + val + 1 ); }
	// END

	enum class BIM{CHR, RS, CM, BP, A1, A2};
	enum class BASE_INDEX {CHR, REF, ALT, STAT, RS, BP, SE, P, MAX };
	enum class FAM {FID, IID, FATHER, MOTHER, SEX, PHENOTYPE};
    enum class GTF{CHR, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME, ATTRIBUTE};
    	// Mean imputed, no-mean imputed, centering is currently too complicated based on our algorithm
    enum class SCORING{MEAN_IMPUTE, SET_ZERO, CENTER};

	template<> struct enumeration_traits< BASE_INDEX > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< GTF > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< FAM > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< BIM > : enumeration_trait_indexing {};
#if defined(__LP64__) || defined(_WIN64)
    typedef std::uint64_t long_type;
	#define ONE  0x1LLU

#else
    typedef std::uint32_t long_type;
	#define ONE  0x1LU
#endif
#endif /* PRSICE_INC_STORAGE_HPP_ */
