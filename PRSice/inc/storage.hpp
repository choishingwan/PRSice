/*
 * storage.hpp
 *
 *  Created on: 14 Oct 2016
 *      Author: shingwanchoi
 */

#ifndef PRSICE_INC_STORAGE_HPP_
#define PRSICE_INC_STORAGE_HPP_
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
	template<> struct enumeration_traits< SNP_Index > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< FAM > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< BIM > : enumeration_trait_indexing {};

	//List of const for use with the GET
	const size_t IID=0;
	const size_t PRS=1;
	const size_t RS=0;
	const size_t LINE=1;
	const size_t CATEGORY=2;
	const size_t INDEX=3;
	const size_t FILENAME=4;
	const size_t THRESHOLD=0;
	const size_t R2 = 1;
	const size_t P = 3;
	const size_t NSNP=2;
	const size_t R2ADJ=4;
	// IID PRS
	typedef std::pair<std::string, double> prs_score;
	// rsid, line number in bim, category, snp_list index
	typedef std::tuple<std::string, size_t, int, size_t, std::string> p_partition;
	// threshold, r2,  num_snps, r2 adjust, p
	typedef std::tuple<double, double, size_t, double, double> PRSice_result;
	//  threshold, r2, num_snps (The r2 is for determining if it is the best)
	typedef std::tuple<double, double, size_t> PRSice_best;

#endif /* PRSICE_INC_STORAGE_HPP_ */
