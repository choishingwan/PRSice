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
#include <mutex>
#include <cstring>
#include <string>
#include <cassert>
#include <fstream>
#include <thread>
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
	PLINK(std::string prefix, bool verbose=false, const size_t thread=1, const catelog &inclusion=catelog());
	static void initialize();
	virtual ~PLINK();
	static void set_chromosome(std::vector<std::string> chr)
	{
		g_chr_list = chr;
	};
	void start_clumping(catelog& inclusion, boost::ptr_vector<SNP> &snp_list, double p_threshold, double r2_threhsold,
			size_t kb_threshold, double proxyy_threshold);
	void clear()
	{
		fclose(m_bedfile);
		m_bedfile = nullptr;
	}
	void get_score(const std::vector<p_partition> &partition,
            const boost::ptr_vector<SNP> &snp_list, std::vector< std::vector<prs_score> > &prs_score,
            size_t start_index, size_t end_bound, size_t num_region, SCORING scoring);
private:
	/*
	 * As I am unfamiliar with the alien language used by Chris, I might run into problem
	 * wrt memory control. So might be better for me to refer to the ld_report_dprime instead
	 * of the clump report part
	 * neat thing is, he also got the threading sorted there XP
	 */

	// bigstack_double_reset(bigstack_mark, bigstack_end_mark); <- how chris clean the memory I guess...
	static std::mutex clump_mtx;
	int32_t load_bim(const catelog &inclusion=catelog());
	int32_t load_fam();
	int32_t load_bed();
	int32_t load_bed(const std::string &bedname);
	void lerase(int num);
	void perform_clump(std::deque<size_t> &clump_snp_index, boost::ptr_vector<SNP> &snp_list,
			size_t &core_snp_index, bool &require_clump, double p_threshold, double r2_threshold,
			size_t kb_threshold, std::string next_chr, size_t next_loc);
	size_t m_thread=1;
	void clump_thread(const size_t c_core_index, const std::deque<size_t> &c_clump_snp_index,
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
	uint32_t load_and_split3(uintptr_t* rawbuf, uint32_t unfiltered_sample_ct, uintptr_t* casebuf,
			uintptr_t* pheno_nm, uintptr_t* pheno_c, uint32_t case_ctv, uint32_t ctrl_ctv, uint32_t do_reverse,
			uint32_t is_case_only, uintptr_t* nm_info_ptr);
	void two_locus_count_table(uintptr_t* lptr1, uintptr_t* lptr2, uint32_t* counts_3x3, uint32_t sample_ctv3,
			uint32_t is_zmiss2);
	void two_locus_count_table_zmiss1(uintptr_t* lptr1, uintptr_t* lptr2, uint32_t* counts_3x3,
			uint32_t sample_ctv3, uint32_t is_zmiss2);
	void compute_clump( size_t core_snp_index, size_t i_start, size_t i_end, boost::ptr_vector<SNP> &snp_list,
			const std::deque<size_t> &clump_snp_index, const double r2_threshold, uintptr_t* geno1,
			bool nm_fixed, uint32_t* tot1);
#ifdef __LP64__
	void two_locus_3x3_tablev(__m128i* vec1, __m128i* vec2, uint32_t* counts_3x3, uint32_t sample_ctv6,
			uint32_t iter_ct);

	inline void two_locus_3x3_zmiss_tablev(__m128i* veca0, __m128i* vecb0, uint32_t* counts_3x3, uint32_t sample_ctv6) {
	  const __m128i m1 = {FIVEMASK, FIVEMASK};
	  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
	  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
	  __m128i* vecb1 = &(vecb0[sample_ctv6]);
	  __m128i* veca1 = &(veca0[sample_ctv6]);
	  __m128i* vend;
	  __m128i loadera0;
	  __m128i loaderb0;
	  __m128i loaderb1;
	  __m128i loadera1;
	  __m128i countx00;
	  __m128i countx01;
	  __m128i countx11;
	  __m128i countx10;
	  __m128i county00;
	  __m128i county01;
	  __m128i county11;
	  __m128i county10;
	  __univec acc00;
	  __univec acc01;
	  __univec acc11;
	  __univec acc10;
	  uint32_t ct2;
	  while (sample_ctv6 >= 30) {
	    sample_ctv6 -= 30;
	    vend = &(veca0[30]);
	    acc00.vi = _mm_setzero_si128();
	    acc01.vi = _mm_setzero_si128();
	    acc11.vi = _mm_setzero_si128();
	    acc10.vi = _mm_setzero_si128();
	    do {
	    two_locus_3x3_zmiss_tablev_outer:
	      loadera0 = *veca0++;
	      loaderb0 = *vecb0++;
	      loaderb1 = *vecb1++;
	      loadera1 = *veca1++;
	      countx00 = _mm_and_si128(loadera0, loaderb0);
	      countx01 = _mm_and_si128(loadera0, loaderb1);
	      countx11 = _mm_and_si128(loadera1, loaderb1);
	      countx10 = _mm_and_si128(loadera1, loaderb0);
	      countx00 = _mm_sub_epi64(countx00, _mm_and_si128(_mm_srli_epi64(countx00, 1), m1));
	      countx01 = _mm_sub_epi64(countx01, _mm_and_si128(_mm_srli_epi64(countx01, 1), m1));
	      countx11 = _mm_sub_epi64(countx11, _mm_and_si128(_mm_srli_epi64(countx11, 1), m1));
	      countx10 = _mm_sub_epi64(countx10, _mm_and_si128(_mm_srli_epi64(countx10, 1), m1));
	      countx00 = _mm_add_epi64(_mm_and_si128(countx00, m2), _mm_and_si128(_mm_srli_epi64(countx00, 2), m2));
	      countx01 = _mm_add_epi64(_mm_and_si128(countx01, m2), _mm_and_si128(_mm_srli_epi64(countx01, 2), m2));
	      countx11 = _mm_add_epi64(_mm_and_si128(countx11, m2), _mm_and_si128(_mm_srli_epi64(countx11, 2), m2));
	      countx10 = _mm_add_epi64(_mm_and_si128(countx10, m2), _mm_and_si128(_mm_srli_epi64(countx10, 2), m2));
	    two_locus_3x3_zmiss_tablev_one_left:
	      loadera0 = *veca0++;
	      loaderb0 = *vecb0++;
	      loaderb1 = *vecb1++;
	      loadera1 = *veca1++;
	      county00 = _mm_and_si128(loadera0, loaderb0);
	      county01 = _mm_and_si128(loadera0, loaderb1);
	      county11 = _mm_and_si128(loadera1, loaderb1);
	      county10 = _mm_and_si128(loadera1, loaderb0);
	      county00 = _mm_sub_epi64(county00, _mm_and_si128(_mm_srli_epi64(county00, 1), m1));
	      county01 = _mm_sub_epi64(county01, _mm_and_si128(_mm_srli_epi64(county01, 1), m1));
	      county11 = _mm_sub_epi64(county11, _mm_and_si128(_mm_srli_epi64(county11, 1), m1));
	      county10 = _mm_sub_epi64(county10, _mm_and_si128(_mm_srli_epi64(county10, 1), m1));
	      countx00 = _mm_add_epi64(countx00, _mm_add_epi64(_mm_and_si128(county00, m2), _mm_and_si128(_mm_srli_epi64(county00, 2), m2)));
	      countx01 = _mm_add_epi64(countx01, _mm_add_epi64(_mm_and_si128(county01, m2), _mm_and_si128(_mm_srli_epi64(county01, 2), m2)));
	      countx11 = _mm_add_epi64(countx11, _mm_add_epi64(_mm_and_si128(county11, m2), _mm_and_si128(_mm_srli_epi64(county11, 2), m2)));
	      countx10 = _mm_add_epi64(countx10, _mm_add_epi64(_mm_and_si128(county10, m2), _mm_and_si128(_mm_srli_epi64(county10, 2), m2)));
	      acc00.vi = _mm_add_epi64(acc00.vi, _mm_add_epi64(_mm_and_si128(countx00, m4), _mm_and_si128(_mm_srli_epi64(countx00, 4), m4)));
	      acc01.vi = _mm_add_epi64(acc01.vi, _mm_add_epi64(_mm_and_si128(countx01, m4), _mm_and_si128(_mm_srli_epi64(countx01, 4), m4)));
	      acc11.vi = _mm_add_epi64(acc11.vi, _mm_add_epi64(_mm_and_si128(countx11, m4), _mm_and_si128(_mm_srli_epi64(countx11, 4), m4)));
	      acc10.vi = _mm_add_epi64(acc10.vi, _mm_add_epi64(_mm_and_si128(countx10, m4), _mm_and_si128(_mm_srli_epi64(countx10, 4), m4)));
	    } while (veca0 < vend);
	    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
	    acc00.vi = _mm_add_epi64(_mm_and_si128(acc00.vi, m8), _mm_and_si128(_mm_srli_epi64(acc00.vi, 8), m8));
	    acc01.vi = _mm_add_epi64(_mm_and_si128(acc01.vi, m8), _mm_and_si128(_mm_srli_epi64(acc01.vi, 8), m8));
	    acc11.vi = _mm_add_epi64(_mm_and_si128(acc11.vi, m8), _mm_and_si128(_mm_srli_epi64(acc11.vi, 8), m8));
	    acc10.vi = _mm_add_epi64(_mm_and_si128(acc10.vi, m8), _mm_and_si128(_mm_srli_epi64(acc10.vi, 8), m8));
	    counts_3x3[0] += ((acc00.u8[0] + acc00.u8[1]) * 0x1000100010001LLU) >> 48;
	    counts_3x3[1] += ((acc01.u8[0] + acc01.u8[1]) * 0x1000100010001LLU) >> 48;
	    counts_3x3[4] += ((acc11.u8[0] + acc11.u8[1]) * 0x1000100010001LLU) >> 48;
	    counts_3x3[3] += ((acc10.u8[0] + acc10.u8[1]) * 0x1000100010001LLU) >> 48;
	  }
	  if (sample_ctv6) {
	    vend = &(veca0[sample_ctv6]);
	    ct2 = sample_ctv6 % 2;
	    sample_ctv6 = 0;
	    acc00.vi = _mm_setzero_si128();
	    acc01.vi = _mm_setzero_si128();
	    acc11.vi = _mm_setzero_si128();
	    acc10.vi = _mm_setzero_si128();
	    if (ct2) {
	      countx00 = _mm_setzero_si128();
	      countx01 = _mm_setzero_si128();
	      countx11 = _mm_setzero_si128();
	      countx10 = _mm_setzero_si128();
	      goto two_locus_3x3_zmiss_tablev_one_left;
	    }
	    goto two_locus_3x3_zmiss_tablev_outer;
	  }
	}
#endif
};

#endif /* PLINK_HPP_ */
