/*
 * genotype.hpp
 *
 *  Created on: 27 Mar 2017
 *      Author: shingwanchoi
 */

#ifndef SRC_GENOTYPE_HPP_
#define SRC_GENOTYPE_HPP_


#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <deque>
#include <Eigen/Dense>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include "commander.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "region.hpp"
#include "snp.hpp"
#include "storage.hpp"

class Genotype {
public:
	Genotype(std::string prefix, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
			bool no_mt=false, const size_t thread=1, bool verbose=false);
	virtual ~Genotype();
	std::unordered_map<std::string, int> get_chr_order() const { return m_chr_order; };
	void read_base(const Commander &c_commander, Region &region);
	void clump(Genotype &reference);

	inline bool existed (const std::string &rs) const
	{
		return m_existed_snps_index.find(rs)!=m_existed_snps_index.end();
	};

	inline bool matched (const SNP &target) const
	{
		if(existed(target.get_rs())){
			return m_existed_snps.at(m_existed_snps_index.at(target.get_rs()))==target;
		}
		return false;
	};
	std::vector<Sample> sample_names() const { return m_sample_names;};
	size_t max_category() const { return m_max_category; };
	bool get_score(Eigen::MatrixXd &current_prs_score, int &cur_index, int &cur_category,
			std::vector<size_t> &num_snp_included);
	bool prepare_prsice();

protected:
	static std::mutex clump_mtx;
	size_t m_max_category = 0;
	size_t m_region_size = 1;
	SCORING m_scoring;
	void lerase(int num);
	std::deque<uintptr_t*> m_genotype;
	struct{
		double r2;
		double proxy;
		double p_value;
		int distance;
		bool use_proxy;
	} clump_info;

	struct{
		double maf;
		double geno;
		double info_score;
		bool filter_maf;
		bool filter_geno;
		bool filter_info;
	} filter;

	void finalize_snps(Region &region, const int distance);

	void set_genotype_files(std::string prefix);
	std::vector<std::string> m_genotype_files;

	void init_chr(int num_auto, bool no_x, bool no_y, bool no_xy, bool no_mt);
	uint32_t m_autosome_ct;
	std::vector<int32_t> m_xymt_codes;
	std::vector<int32_t> m_chrom_start;
	uint32_t m_max_code;
	uintptr_t* m_haploid_mask;
	uintptr_t* m_chrom_mask;


	virtual std::vector<Sample> load_samples(){ return std::vector<Sample>(0); };
	uintptr_t m_unfiltered_sample_ct = 0;
	uintptr_t m_unfiltered_sample_ctl = 0;
	uintptr_t m_unfiltered_sample_ct4 = 0;
	std::vector<Sample> m_sample_names;

	void init_sample_vectors(){};
	uintptr_t* m_founder_info = nullptr;
	uintptr_t* m_sex_male = nullptr;
	uintptr_t* m_sample_exclude = nullptr;
	size_t m_num_male=0;
	size_t m_num_female=0;
	size_t m_num_ambig_sex=0;
	uintptr_t m_founder_ct = 0;

	virtual std::vector<SNP> load_snps();
	uintptr_t m_unfiltered_marker_ct = 0;
	uintptr_t m_unfiltered_marker_ctl = 0;
	uintptr_t m_marker_ct = 0;
	uintptr_t m_marker_exclude_ct = 0;
	uintptr_t* m_marker_exclude = nullptr;
	std::unordered_map<std::string, size_t> m_existed_snps_index;
	std::vector<SNP> m_existed_snps;
	std::unordered_map<std::string, int> m_chr_order;
	uint32_t m_hh_exists; // might be a bit harsh, but should also read in maf when loading SNPs
	uint32_t m_num_ambig;

	uint32_t m_thread;
	virtual void read_genotype(uintptr_t* genotype, const uint32_t snp_index, const std::string &file_name){genotype=nullptr;};
	void  virtual read_score(std::vector< std::vector<Sample> > &current_prs_score, size_t start_index, size_t end_bound);
	//hh_exists
	inline bool ambiguous(std::string ref_allele, std::string alt_allele)
	{
		return (ref_allele == "A" && alt_allele == "T")
				|| (ref_allele == "a" && alt_allele == "t")
				|| (ref_allele == "G" && alt_allele == "C")
				|| (ref_allele == "g" && alt_allele == "c");
	}

	void perform_clump(int core_genotype_index, bool require_clump, int chr, int loc,
			std::deque<int> &clump_index);
	void clump_thread(const size_t c_core_genotype_index, std::deque<int> &clump_index);
	void compute_clump( size_t core_genotype_index, size_t i_start, size_t i_end,
			uintptr_t* geno1, bool nm_fixed, uint32_t* tot1, const std::deque<int> &clump_index);




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

/*
class Plink: public Genotype{
public:
	Plink(std::string prefix, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
			bool no_mt=false, const size_t thread=1);
private:
	void load_sample();
	std::vector<SNP> load_snps();
};
*/

class BinaryPlink: public Genotype{
public:
	BinaryPlink(std::string prefix, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
			bool no_mt=false, const size_t thread=1, bool verbose=false);
	 ~BinaryPlink();
private:
	uintptr_t m_bed_offset = 3;
	std::vector<Sample> load_samples();
	std::vector<SNP> load_snps();
	void check_bed();
	void read_genotype(uintptr_t* genotype, const uint32_t snp_index, const std::string &file_name);
	void read_score(Eigen::MatrixXd &current_prs_score, size_t start_index, size_t end_bound);
	FILE* m_bedfile = nullptr;
	std::string m_cur_file;
	uintptr_t m_final_mask;
	uintptr_t *m_tmp_genotype;
};


class GenomeFactory {
private:
	std::unordered_map<std::string, int> file_type {
		{ "bed", 0 },
		{ "ped", 1 },
		{ "bgen", 2}
	};
public:
	std::unique_ptr<Genotype> createGenotype(const Commander &commander, const std::string &prefix,
			const std::string &type, bool verbose)
	{
		fprintf(stderr, "Loading Genotype file: %s ", prefix.c_str());
		int code = (file_type.find(type)!=file_type.end())? file_type[type]: 0;
		switch(code)
		{
		case 1:
			/*
			return std::unique_ptr<Genotype>(new Plink(prefix, commander.num_auto(),
							commander.no_x(), commander.no_y(), commander.no_xy(), commander.no_mt(),
							commander.thread(), verbose));
							*/
		case 2:
		default:
		case 0:
			fprintf(stderr, "(bed)\n");
			return std::unique_ptr<Genotype>(new BinaryPlink(prefix, commander.num_auto(),
							commander.no_x(), commander.no_y(), commander.no_xy(), commander.no_mt(),
							commander.thread(), verbose));

		}
	}

};

#endif /* SRC_GENOTYPE_HPP_ */
