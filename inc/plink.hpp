/*
 * plink.hpp
 *
 *  Created on: 19 Feb 2017
 *      Author: shingwan
 */

#ifndef PLINK_HPP_
#define PLINK_HPP_

#include <stdexcept>
#include <stdio.h>
#include <cstring>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <deque>
#include <unistd.h> // for _SC_PHYS_PAGES
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
#include "snp.hpp"
#include "plink_common.hpp"
#include "plink_set.hpp"
#include "storage.hpp"
#include "misc.hpp"

class PLINK {
public:
	typedef std::unordered_map<std::string, size_t> catelog;
	PLINK(std::string prefix, const size_t thread=1, const catelog &inclusion=catelog());
	static void initialize();
	virtual ~PLINK();
	static void set_chromosome(std::vector<std::string> chr)
	{
		g_chr_list = chr;
	};
	void start_clumping(boost::ptr_vector<SNP> &snp_list, double p_threshold, double r2_threhsold,
			size_t kb_threshold, double proxy);
private:
	/*
	 * As I am unfamiliar with the alien language used by Chris, I might run into problem
	 * wrt memory control. So might be better for me to refer to the ld_report_dprime instead
	 * of the clump report part
	 * neat thing is, he also got the threading sorted there XP
	 */

	// bigstack_double_reset(bigstack_mark, bigstack_end_mark); <- how chris clean the memory I guess...
	int32_t load_bim(const catelog &inclusion=catelog());
	int32_t load_fam();
	int32_t load_bed();
	void lerase(int num);
	void perform_clump(std::deque<size_t> &clump_snp_index, boost::ptr_vector<SNP> &snp_list,
			size_t &core_snp_index, bool &require_clump, double p_threshold, double r2_threshold,
			size_t kb_threshold, std::string next_chr, size_t next_loc);
	size_t m_thread=1;
	void void PLINK::clump_thread(const size_t c_core_index, const std::deque<size_t> &c_clump_snp_index,
			boost::ptr_vector<SNP> &snp_list, const double c_r2_threshold);
	static std::vector<std::string> g_chr_list;
	FILE* m_bedfile = nullptr;
	std::vector<std::string> m_prefix;
	std::vector<snp_link> m_snp_link;
	std::vector<uintptr_t*> m_genotype;
	size_t m_num_male;
	size_t m_num_female;
	size_t m_num_ambig_sex;
	uintptr_t m_bed_offset = 3;
	uintptr_t m_unfiltered_marker_ct = 0;
	uintptr_t m_marker_ct = 0;
	uintptr_t m_marker_exclude_ct = 0;
	uintptr_t m_unfiltered_sample_ct = 0;
	uintptr_t m_unfiltered_sample_ct4 = 0;
	uintptr_t m_unfiltered_sample_ctl = 0;
	uintptr_t m_founder_ct = 0;
	uintptr_t* m_founder_info = nullptr;
	uintptr_t* m_sex_male = nullptr;
	uintptr_t* m_sample_exclude = nullptr;
	uintptr_t* m_marker_exclude = nullptr;
	uintptr_t* m_marker_reverse = nullptr;

	uint32_t em_phase_hethet(double known11, double known12, double known21, double known22, uint32_t center_ct,
			double* freq1x_ptr, double* freq2x_ptr, double* freqx1_ptr, double* freqx2_ptr, double* freq11_ptr,
			uint32_t* onside_sol_ct_ptr);
	uint32_t em_phase_hethet_nobase(uint32_t* counts, uint32_t is_x1, uint32_t is_x2, double* freq1x_ptr,
			double* freq2x_ptr, double* freqx1_ptr, double* freqx2_ptr, double* freq11_ptr);
	double calc_lnlike(double known11, double known12, double known21, double known22, double center_ct_d,
			double freq11, double freq12, double freq21, double freq22, double half_hethet_share, double freq11_incr);

};

#endif /* PLINK_HPP_ */
