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
#include "reporter.hpp"
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
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef __APPLE__
#include <sys/sysctl.h> // sysctl()
#endif

class Genotype
{
public:
    Genotype(){};
    Genotype(const size_t thread, const bool ignore_fid,
             const bool keep_nonfounder, const bool keep_ambig)
        : m_thread(thread)
        , m_keep_nonfounder(keep_nonfounder)
        , m_ignore_fid(ignore_fid)
        , m_keep_ambig(keep_ambig){};
    virtual ~Genotype();


    void load_samples(const std::string& keep_file,
                      const std::string& remove_file, bool verbose,
                      Reporter& reporter);
    void load_snps(const std::string out_prefix,
                   const std::string& extract_file,
                   const std::string& exclude_file, const double geno,
                   const double maf, const double info,
                   const double hard_threshold, const bool hard_coded,
                   bool verbose, Reporter& reporter);

    void load_snps(const std::string out_prefix,
                   const std::unordered_map<std::string, size_t>& existed_snps,
                   const double geno, const double maf, const double info,
                   const double hard_threshold, const bool hard_coded,
                   bool verbose, Reporter& reporter);
    std::unordered_map<std::string, int> get_chr_order() const
    {
        return m_chr_order;
    };

    std::unordered_map<std::string, size_t> index() const
    {
        return m_existed_snps_index;
    };
    std::vector<double> get_thresholds() const { return m_thresholds; };
    std::vector<Sample> sample_names() const { return m_sample_names; };
    size_t max_category() const { return m_max_category; };
    size_t num_sample() const { return m_sample_names.size(); }
    bool get_score(int& cur_index, int& cur_category, double& cur_threshold,
                   size_t& num_snp_included, const size_t region_index,
                   const bool require_statistic);
    bool sort_by_p();
    void print_snp(std::string& output, double threshold,
                   const size_t region_index);
    size_t num_threshold() const { return m_num_threshold; };
    void read_base(const Commander& c_commander, Region& region,
                   Reporter& reporter);
    // void clump(Genotype& reference);
    void efficient_clumping(Genotype& reference, Reporter& reporter);
    void set_info(const Commander& c_commander, const bool ld = false);
    void reset_sample()
    {
        for (auto&& sample : m_sample_names) {
            sample.num_snp = 0;
            sample.prs = 0.0;
            sample.has_pheno = false;
        }
    };
    void reset_prs()
    {
        for (auto&& sample : m_sample_names) {
            sample.num_snp = 0;
            sample.prs = 0.0;
        }
    };
    bool prepare_prsice();
    std::string sample_id(size_t i) const
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        if (m_ignore_fid)
            return m_sample_names[i].IID;
        else
            return m_sample_names[i].FID + "_" + m_sample_names[i].IID;
    }
    bool sample_included(size_t i) const
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        return m_sample_names[i].included;
    };
    void got_pheno(size_t i)
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        m_sample_names[i].has_pheno = true;
    }
    void invalid_pheno(size_t i)
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        m_sample_names[i].has_pheno = false;
    }
    bool has_pheno(size_t i) const
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        return m_sample_names[i].has_pheno;
    }
    std::string pheno(size_t i) const
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        return m_sample_names[i].pheno;
    }
    bool pheno_is_na(size_t i) const
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        return m_sample_names[i].pheno.compare("NA") == 0;
    }
    std::string fid(size_t i) const
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        return m_sample_names[i].FID;
    }
    std::string iid(size_t i) const
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        return m_sample_names[i].IID;
    }
    double calculate_score(SCORING score_type, size_t i) const
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        double score = 0;
        switch (score_type)
        {
        case SCORING::AVERAGE:
            if (m_sample_names[i].num_snp == 0) return 0;
            return m_sample_names[i].prs / (double) m_sample_names[i].num_snp;
            break;
        case SCORING::SUM: return m_sample_names[i].prs; break;
        case SCORING::STANDARDIZE:
            if (m_sample_names[i].num_snp == 0)
                return -m_mean_score / m_score_sd;
            else
                return ((m_sample_names[i].prs
                         / (double) m_sample_names[i].num_snp)
                        - m_mean_score)
                       / m_score_sd;
            break;
        }
        return score;
    }

protected:
    // variable storages
    // vector storing all the genotype files
    std::vector<std::string> m_genotype_files;
    // sample file name. Fam for plink
    std::string m_sample_file;
    /** chromosome information **/
    // here such that we can extend PRSice to other organism if we have time
    // also use for handling sex chromosome (which is something that require
    // further research)
    std::vector<int32_t> m_xymt_codes;
    std::vector<int32_t> m_chrom_start;
    uint32_t m_autosome_ct;
    uint32_t m_max_code;
    std::vector<uintptr_t> m_haploid_mask;
    std::vector<uintptr_t> m_chrom_mask;
    // misc variables
    uint32_t m_thread = 1;
    bool m_keep_nonfounder = false;
    bool m_ignore_fid = false;
    bool m_keep_ambig = false;
    // sample information
    std::vector<Sample> m_sample_names;
    bool m_remove_sample = true;
    std::unordered_set<std::string> m_sample_selection_list;
    uintptr_t m_unfiltered_sample_ct = 0; // number of unfiltered samples
    uintptr_t m_sample_ct = 0;            // number of final samples
    std::vector<uintptr_t> m_founder_info;
    std::vector<uintptr_t> m_sample_include;
    std::vector<uintptr_t> m_sex_male;
    size_t m_num_male = 0;
    size_t m_num_female = 0;
    size_t m_num_ambig_sex = 0;
    size_t m_num_non_founder = 0;
    uintptr_t m_founder_ct = 0;
    // snp information
    bool m_exclude_snp = true;
    std::unordered_set<std::string> m_snp_selection_list;
    uintptr_t m_marker_ct = 0;
    uint32_t m_num_ambig = 0;
    uint32_t m_num_maf_filter = 0;
    uint32_t m_num_geno_filter = 0;
    uint32_t m_num_info_filter = 0;
    double m_hard_threshold;
    bool m_hard_coded = false;
    std::unordered_map<std::string, size_t> m_existed_snps_index;
    std::vector<size_t> m_sort_by_p_index;
    std::vector<SNP> m_existed_snps;
    std::unordered_map<std::string, int> m_chr_order;
    std::vector<uintptr_t> m_tmp_genotype;
    // functions
    // function to substitute the # in the sample name
    std::vector<std::string> set_genotype_files(const std::string& prefix);
    void init_chr(int num_auto = 22, bool no_x = false, bool no_y = false,
                  bool no_xy = false, bool no_mt = false);
    // responsible for reading in the sample
    virtual std::vector<Sample> gen_sample_vector()
    {
        return std::vector<Sample>(0);
    };
    virtual std::vector<SNP> gen_snp_vector(const double geno, const double maf,
                                            const double info,
                                            const double hard_threshold,
                                            const bool hard_coded,
                                            const std::string& out_prefix)
    {
        return std::vector<SNP>(0);
    };
    // for loading the sample inclusion / exclusion set
    std::unordered_set<std::string> load_ref(std::string input,
                                             bool ignore_fid);
    // for loading the SNP inclusion / exclusion set
    std::unordered_set<std::string> load_snp_list(std::string input,
                                                  Reporter& reporter);
    double get_r2(bool core_missing, bool pair_missing,
                  std::vector<uint32_t>& core_tot,
                  std::vector<uint32_t>& pair_tot,
                  std::vector<uintptr_t>& genotype_vector,
                  std::vector<uintptr_t>& pair_genotype_vector);
    double get_r2(bool core_missing, std::vector<uint32_t>& index_tots,
                  std::vector<uintptr_t>& index_data,
                  std::vector<uintptr_t>& genotype_vector);
    /** Misc information **/
    size_t m_max_category = 0;
    size_t m_region_size = 1;
    size_t m_num_threshold = 0;
    MODEL m_model = MODEL::ADDITIVE;
    MISSING_SCORE m_missing_score = MISSING_SCORE::MEAN_IMPUTE;
    SCORING m_scoring = SCORING::AVERAGE;
    double m_mean_score = 0.0;
    double m_score_sd = 0.0;
    struct
    {
        double r2;
        double proxy;
        double p_value;
        int distance;
        bool use_proxy;
    } clump_info;


    std::vector<double> m_thresholds;


    uintptr_t m_unfiltered_marker_ct = 0;
    // uint32_t m_hh_exists;


    virtual inline void read_genotype(uintptr_t* genotype, const SNP& snp,
                                      const std::string& file_name){};
    virtual void read_score(size_t start_index, size_t end_bound,
                            const size_t region_index){};


    // hh_exists
    inline bool ambiguous(const std::string& ref_allele,
                          const std::string& alt_allele) const
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

    /*
    void perform_clump(size_t& core_genotype_index, int& begin_index,
                       int current_index, bool require_clump);
    void clump_thread(const size_t c_core_genotype_index,
                      const size_t c_begin_index, const size_t c_current_index);
    void compute_clump(size_t core_genotype_index, size_t i_start, size_t i_end,
                       bool nm_fixed, uint32_t* tot1);
*/

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
