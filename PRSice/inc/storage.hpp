/*
 * storage.hpp
 *
 *  Created on: 14 Oct 2016
 *      Author: shingwanchoi
 */

#ifndef PRSICE_INC_STORAGE_HPP_
#define PRSICE_INC_STORAGE_HPP_
	// From http://stackoverflow.com/a/12927952/1441789
	enum class SNP_Index {CHR, REF, ALT, STAT, RS, BP, SE, P, MAX };
	enum class FAM {FID, IID, FATHER, MOTHER, SEX, PHENOTYPE};
	enum class BIM{CHR, RS, CM, BP, A1, A2};
    enum class FILE_INFO { FILE, LINE, INDEX  }; // This is for clumping in PLINK
    enum class PRS{IID=0, PRS, RS=0, LINE, CATEGORY, INDEX, FILENAME, THRESHOLD=0, R2, NSNP, P, R2ADJ};
    enum help_index{CATEGORY, SHORT, LONG, DESCRIPTION};
    int operator + ( SNP_Index val ) { return static_cast< int >( val ); }
    int operator + ( FAM val ) { return static_cast< int >( val ); }
    int operator + ( BIM val ) { return static_cast< int >( val ); }
    int operator + ( FILE_INFO val ) { return static_cast< int >( val ); }
    int operator + ( PRS val ) { return static_cast< int >( val ); }
    //END
	//List of const for use with the GET
	// IID PRS
	typedef std::pair<std::string, double> prs_score;
	// rsid, line number in bim, category, snp_list index
	typedef std::tuple<std::string, size_t, int, size_t, std::string> p_partition;
	// threshold, r2,  num_snps, p, r2 adjust
	typedef std::tuple<double, double, size_t, double, double> PRSice_result;
	//  threshold, r2, num_snps (The r2 is for determining if it is the best)
	typedef std::tuple<double, double, size_t> PRSice_best;
    typedef std::tuple<std::string, char, std::string, std::string> help;

#endif /* PRSICE_INC_STORAGE_HPP_ */
