/*
 * storage.hpp
 *
 *  Created on: 14 Oct 2016
 *      Author: shingwanchoi
 */

#ifndef PRSICE_INC_STORAGE_HPP_
#define PRSICE_INC_STORAGE_HPP_

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
