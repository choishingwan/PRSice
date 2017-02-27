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
    enum class FILE_INFO { FILE, LINE, INDEX, HAPLOID, X, Y  }; // This is for clumping in PLINK
    enum class PRS{FID=0, IID, PRS, NNMISS, RS=0, LINE, CATEGORY, INDEX, FILENAME, P_THRES, THRESHOLD=0, R2,
    	NSNP, COEFF, P, EMPIRICAL_P, R2ADJ};
    enum class BOUNDARY{CHR, START, END};
    enum class GTF{CHR, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME, ATTRIBUTE};
    	// Mean imputed, no-mean imputed, centering is currently too complicated based on our algorithm
    enum class SCORING{MEAN_IMPUTE, SET_ZERO, CENTER};
    enum help_index{CATEGORY, SHORT, LONG, DESCRIPTION};
    template<> struct enumeration_traits< BOUNDARY > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< GTF > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< SNP_Index > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< FAM > : enumeration_trait_indexing {};
	template<> struct enumeration_traits< BIM > : enumeration_trait_indexing {};
    template<> struct enumeration_traits< FILE_INFO > : enumeration_trait_indexing {};
    template<> struct enumeration_traits< PRS > : enumeration_trait_indexing {};
    template<> struct enumeration_traits< SCORING > : enumeration_trait_indexing {};
	//List of const for use with the GET
	// FID IID PRS Number of non-missing SNP
	typedef std::tuple<std::string, std::string, double, size_t> prs_score;
	// rsid, line number in bim, category, snp_list index
	typedef std::tuple<std::string, size_t, int, size_t, std::string, double> p_partition;
	// threshold, r2,  num_snps, p, coefficient r2 adjust, number of better
	typedef std::tuple<double, double, size_t, double, double, size_t, double> PRSice_result;
	//  threshold, r2, num_snps, coefficient, pvalue (The r2 is for determining if it is the best), number of better
	typedef std::tuple<double, double, size_t, double, double, int> PRSice_best;
    typedef std::tuple<std::string, char, std::string, std::string> help;
    typedef std::tuple<std::string, size_t, size_t> boundary;
    typedef std::tuple<std::string, int, size_t, bool, bool, bool> snp_link;


#if defined(__LP64__) || defined(_WIN64)
    typedef std::uint64_t long_type;
	#define ONE  0x1LLU

#else
    typedef std::uint32_t long_type;
	#define ONE  0x1LU
#endif
#endif /* PRSICE_INC_STORAGE_HPP_ */
