/*
 * storage.hpp
 *
 *  Created on: 14 Oct 2016
 *      Author: shingwanchoi
 */

#ifndef PRSICE_INC_STORAGE_HPP_
#define PRSICE_INC_STORAGE_HPP_
#include <memory>
	// From http://stackoverflow.com/a/12927952/1441789
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

	enum class SNP_Index {CHR, REF, ALT, STAT, RS, BP, SE, P, MAX };
	enum class FAM {FID, IID, FATHER, MOTHER, SEX, PHENOTYPE};
	enum class BIM{CHR, RS, CM, BP, A1, A2};
    enum class FILE_INFO { FILE, LINE, INDEX  }; // This is for clumping in PLINK
    enum class PRS{IID=0, PRS, RS=0, LINE, CATEGORY, INDEX, FILENAME, THRESHOLD=0, R2, NSNP, COEFF, P, R2ADJ};
    enum help_index{CATEGORY, SHORT, LONG, DESCRIPTION};
	template<> struct enumeration_traits< SNP_Index > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< FAM > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< BIM > : enumeration_trait_indexing {};
    template<> struct enumeration_traits< FILE_INFO > : enumeration_trait_indexing {};
    template<> struct enumeration_traits< PRS > : enumeration_trait_indexing {};
	//List of const for use with the GET
	// IID PRS
	typedef std::pair<std::string, double> prs_score;
	// rsid, line number in bim, category, snp_list index
	typedef std::tuple<std::string, size_t, int, size_t, std::string> p_partition;
	// threshold, r2,  num_snps, p, coefficient r2 adjust
	typedef std::tuple<double, double, size_t, double, double, double> PRSice_result;
	//  threshold, r2, num_snps, coefficient (The r2 is for determining if it is the best)
	typedef std::tuple<double, double, size_t, double> PRSice_best;
    typedef std::tuple<std::string, char, std::string, std::string> help;

#endif /* PRSICE_INC_STORAGE_HPP_ */
