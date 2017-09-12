// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef GENOTYPE_H
#define GENOTYPE_H

#include "commander.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "region.hpp"
#include "snp.hpp"
#include "storage.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <deque>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Genotype
{
public:
    Genotype(){};
    Genotype(std::string prefix, std::string remove_sample,
             std::string keep_sample, std::string extract_snp,
             std::string exclude_snp, bool ignore_fid, int num_auto = 22,
             bool no_x = false, bool no_y = false, bool no_xy = false,
             bool no_mt = false, bool keep_ambig = false,
             const size_t thread = 1, bool verbose = false);

    virtual ~Genotype();
    std::unordered_map<std::string, int> get_chr_order() const
    {
        return m_chr_order;
    };

    std::vector<double> get_thresholds() const { return m_thresholds; };
    std::vector<Sample> sample_names() const { return m_sample_names; };
    size_t max_category() const { return m_max_category; };
    bool get_score(misc::vec2d<Sample_lite>& current_prs_score, int& cur_index,
                   int& cur_category, double& cur_threshold,
                   std::vector<size_t>& num_snp_included);
    bool prepare_prsice();
    void print_snp(std::string& output, double threshold);
    size_t num_threshold() const { return m_num_threshold; };
    void update_include(const std::vector<Sample>& inclusion);
    void read_base(const Commander& c_commander, Region& region);
    void clump(Genotype& reference);
    void set_clump_info(const Commander& c_commander);

protected:
    // for loading the sample inclusion / exclusion set
    std::unordered_set<std::string> load_ref(std::string input,
                                             bool ignore_fid);
    // for loading the SNP inclusion / exclusion set
    std::unordered_set<std::string> load_snp_list(std::string input);

    // for processing the genotype file (mainly for multiple chromosome input)
    void set_genotype_files(std::string prefix);
    // storage of genotype file names
    std::vector<std::string> m_genotype_files;


    /** chromosome information **/

    void init_chr(int num_auto, bool no_x, bool no_y, bool no_xy, bool no_mt);
    std::vector<int32_t> m_xymt_codes;
    std::vector<int32_t> m_chrom_start;
    uint32_t m_autosome_ct;
    uint32_t m_max_code;
    std::vector<uintptr_t> m_haploid_mask;
    std::vector<uintptr_t> m_chrom_mask;

    // for easy calculatable items, remove them from the class to keep things
    // more organized and easier to debug

    /** sample information **/
    virtual std::vector<Sample> load_samples(bool ignore_fid)
    {
        return std::vector<Sample>(0);
    };
    uintptr_t m_unfiltered_sample_ct = 0; //number of unfiltered samples
    uintptr_t m_sample_ct = 0; // number of final samples
    // storage of sample information
    std::vector<Sample> m_sample_names; // vector storing the sample information

    /** SNP/Sample selection **/
    // default is remove (if not provided, the list is empty,
    // thus nothing will be removed)
    bool m_remove_sample = true;
    std::unordered_set<std::string> m_sample_selection_list;

    bool m_exclude_snp = true;
    std::unordered_set<std::string> m_snp_selection_list;

    // tempory storage for genotype information
    // use for plink read
    std::vector<uintptr_t> m_tmp_genotype;

    static std::mutex clump_mtx;

    /** Misc information **/
    size_t m_max_category = 0;
    size_t m_region_size = 1;
    size_t m_num_threshold = 0;
    SCORING m_scoring;

    struct
    {
        double r2;
        double proxy;
        double p_value;
        int distance;
        bool use_proxy;
    } clump_info;

    struct
    {
        double maf;
        double geno;
        double info_score;
        double hard_threshold;
        bool filter_hard_threshold;
        bool filter_maf;
        bool filter_geno;
        bool filter_info;
        bool keep_ambig;
        bool use_hard;
    } filter;


    std::vector<double> m_thresholds;


    // founder_info stores the samples that are included from
    // the genotype file. All non-founder and removed samples
    // will have bit set to 0
    std::vector<uintptr_t> m_founder_info;

    std::vector<uintptr_t> m_sex_male;
    std::vector<uintptr_t> m_sample_include;
    size_t m_num_male = 0;
    size_t m_num_female = 0;
    size_t m_num_ambig_sex = 0;
    size_t m_num_non_founder = 0;
    uintptr_t m_founder_ct = 0;

    virtual std::vector<SNP> load_snps() { return std::vector<SNP>(0); };

    uintptr_t m_unfiltered_marker_ct = 0;
    uintptr_t m_marker_ct = 0;
    uint32_t m_num_ambig =0;

    std::unordered_map<std::string, size_t> m_existed_snps_index;
    std::vector<SNP> m_existed_snps;
    std::unordered_map<std::string, int> m_chr_order;
    // uint32_t m_hh_exists;



    uint32_t m_thread;
    virtual inline void read_genotype(uintptr_t* genotype,
                                      const uint32_t snp_index,
                                      const std::string& file_name)
    {
        genotype = nullptr;
    };
    virtual void read_score(misc::vec2d<Sample_lite>& current_prs_score,
                            size_t start_index, size_t end_bound){};

    // hh_exists
    inline bool ambiguous(std::string ref_allele, std::string alt_allele)
    {
        return (ref_allele == "A" && alt_allele == "T")
               || (ref_allele == "a" && alt_allele == "t")
               || (ref_allele == "G" && alt_allele == "C")
               || (ref_allele == "g" && alt_allele == "c")
               || (alt_allele == "A" && ref_allele == "T")
               || (alt_allele == "a" && ref_allele == "t")
               || (alt_allele == "G" && ref_allele == "C")
               || (alt_allele == "g" && ref_allele == "c");
    };

    void perform_clump(size_t& core_genotype_index, int& begin_index,
                       int current_index, bool require_clump);
    void clump_thread(const size_t c_core_genotype_index,
                      const size_t c_begin_index, const size_t c_current_index);
    void compute_clump(size_t core_genotype_index, size_t i_start, size_t i_end,
                       bool nm_fixed, uint32_t* tot1);


    // no touchy area (PLINK Code)
    uint32_t em_phase_hethet(double known11, double known12, double known21,
                             double known22, uint32_t center_ct,
                             double* freq1x_ptr, double* freq2x_ptr,
                             double* freqx1_ptr, double* freqx2_ptr,
                             double* freq11_ptr, uint32_t* onside_sol_ct_ptr);
    uint32_t em_phase_hethet_nobase(uint32_t* counts, uint32_t is_x1,
                                    uint32_t is_x2, double* freq1x_ptr,
                                    double* freq2x_ptr, double* freqx1_ptr,
                                    double* freqx2_ptr, double* freq11_ptr);
    double calc_lnlike(double known11, double known12, double known21,
                       double known22, double center_ct_d, double freq11,
                       double freq12, double freq21, double freq22,
                       double half_hethet_share, double freq11_incr);
    uint32_t load_and_split3(uintptr_t* rawbuf, uint32_t unfiltered_sample_ct,
                             uintptr_t* casebuf, uint32_t case_ctv,
                             uint32_t ctrl_ctv, uint32_t do_reverse,
                             uint32_t is_case_only, uintptr_t* nm_info_ptr);
    void two_locus_count_table(uintptr_t* lptr1, uintptr_t* lptr2,
                               uint32_t* counts_3x3, uint32_t sample_ctv3,
                               uint32_t is_zmiss2);
    void two_locus_count_table_zmiss1(uintptr_t* lptr1, uintptr_t* lptr2,
                                      uint32_t* counts_3x3,
                                      uint32_t sample_ctv3, uint32_t is_zmiss2);
#ifdef __LP64__
    void two_locus_3x3_tablev(__m128i* vec1, __m128i* vec2,
                              uint32_t* counts_3x3, uint32_t sample_ctv6,
                              uint32_t iter_ct);

    inline void two_locus_3x3_zmiss_tablev(__m128i* veca0, __m128i* vecb0,
                                           uint32_t* counts_3x3,
                                           uint32_t sample_ctv6)
    {
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
            do
            {
            two_locus_3x3_zmiss_tablev_outer:
                loadera0 = *veca0++;
                loaderb0 = *vecb0++;
                loaderb1 = *vecb1++;
                loadera1 = *veca1++;
                countx00 = _mm_and_si128(loadera0, loaderb0);
                countx01 = _mm_and_si128(loadera0, loaderb1);
                countx11 = _mm_and_si128(loadera1, loaderb1);
                countx10 = _mm_and_si128(loadera1, loaderb0);
                countx00 = _mm_sub_epi64(
                    countx00, _mm_and_si128(_mm_srli_epi64(countx00, 1), m1));
                countx01 = _mm_sub_epi64(
                    countx01, _mm_and_si128(_mm_srli_epi64(countx01, 1), m1));
                countx11 = _mm_sub_epi64(
                    countx11, _mm_and_si128(_mm_srli_epi64(countx11, 1), m1));
                countx10 = _mm_sub_epi64(
                    countx10, _mm_and_si128(_mm_srli_epi64(countx10, 1), m1));
                countx00 = _mm_add_epi64(
                    _mm_and_si128(countx00, m2),
                    _mm_and_si128(_mm_srli_epi64(countx00, 2), m2));
                countx01 = _mm_add_epi64(
                    _mm_and_si128(countx01, m2),
                    _mm_and_si128(_mm_srli_epi64(countx01, 2), m2));
                countx11 = _mm_add_epi64(
                    _mm_and_si128(countx11, m2),
                    _mm_and_si128(_mm_srli_epi64(countx11, 2), m2));
                countx10 = _mm_add_epi64(
                    _mm_and_si128(countx10, m2),
                    _mm_and_si128(_mm_srli_epi64(countx10, 2), m2));
            two_locus_3x3_zmiss_tablev_one_left:
                loadera0 = *veca0++;
                loaderb0 = *vecb0++;
                loaderb1 = *vecb1++;
                loadera1 = *veca1++;
                county00 = _mm_and_si128(loadera0, loaderb0);
                county01 = _mm_and_si128(loadera0, loaderb1);
                county11 = _mm_and_si128(loadera1, loaderb1);
                county10 = _mm_and_si128(loadera1, loaderb0);
                county00 = _mm_sub_epi64(
                    county00, _mm_and_si128(_mm_srli_epi64(county00, 1), m1));
                county01 = _mm_sub_epi64(
                    county01, _mm_and_si128(_mm_srli_epi64(county01, 1), m1));
                county11 = _mm_sub_epi64(
                    county11, _mm_and_si128(_mm_srli_epi64(county11, 1), m1));
                county10 = _mm_sub_epi64(
                    county10, _mm_and_si128(_mm_srli_epi64(county10, 1), m1));
                countx00 = _mm_add_epi64(
                    countx00,
                    _mm_add_epi64(
                        _mm_and_si128(county00, m2),
                        _mm_and_si128(_mm_srli_epi64(county00, 2), m2)));
                countx01 = _mm_add_epi64(
                    countx01,
                    _mm_add_epi64(
                        _mm_and_si128(county01, m2),
                        _mm_and_si128(_mm_srli_epi64(county01, 2), m2)));
                countx11 = _mm_add_epi64(
                    countx11,
                    _mm_add_epi64(
                        _mm_and_si128(county11, m2),
                        _mm_and_si128(_mm_srli_epi64(county11, 2), m2)));
                countx10 = _mm_add_epi64(
                    countx10,
                    _mm_add_epi64(
                        _mm_and_si128(county10, m2),
                        _mm_and_si128(_mm_srli_epi64(county10, 2), m2)));
                acc00.vi = _mm_add_epi64(
                    acc00.vi,
                    _mm_add_epi64(
                        _mm_and_si128(countx00, m4),
                        _mm_and_si128(_mm_srli_epi64(countx00, 4), m4)));
                acc01.vi = _mm_add_epi64(
                    acc01.vi,
                    _mm_add_epi64(
                        _mm_and_si128(countx01, m4),
                        _mm_and_si128(_mm_srli_epi64(countx01, 4), m4)));
                acc11.vi = _mm_add_epi64(
                    acc11.vi,
                    _mm_add_epi64(
                        _mm_and_si128(countx11, m4),
                        _mm_and_si128(_mm_srli_epi64(countx11, 4), m4)));
                acc10.vi = _mm_add_epi64(
                    acc10.vi,
                    _mm_add_epi64(
                        _mm_and_si128(countx10, m4),
                        _mm_and_si128(_mm_srli_epi64(countx10, 4), m4)));
            } while (veca0 < vend);
            const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
            acc00.vi =
                _mm_add_epi64(_mm_and_si128(acc00.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc00.vi, 8), m8));
            acc01.vi =
                _mm_add_epi64(_mm_and_si128(acc01.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc01.vi, 8), m8));
            acc11.vi =
                _mm_add_epi64(_mm_and_si128(acc11.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc11.vi, 8), m8));
            acc10.vi =
                _mm_add_epi64(_mm_and_si128(acc10.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc10.vi, 8), m8));
            counts_3x3[0] +=
                ((acc00.u8[0] + acc00.u8[1]) * 0x1000100010001LLU) >> 48;
            counts_3x3[1] +=
                ((acc01.u8[0] + acc01.u8[1]) * 0x1000100010001LLU) >> 48;
            counts_3x3[4] +=
                ((acc11.u8[0] + acc11.u8[1]) * 0x1000100010001LLU) >> 48;
            counts_3x3[3] +=
                ((acc10.u8[0] + acc10.u8[1]) * 0x1000100010001LLU) >> 48;
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
 Plink(std::string prefix, int num_auto=22, bool no_x=false, bool no_y=false,
 bool no_xy=false, bool no_mt=false, const size_t thread=1); private: void
 load_sample(); std::vector<SNP> load_snps();
 };
 */

#endif /* SRC_GENOTYPE_HPP_ */
