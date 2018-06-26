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
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef __APPLE__
#include <sys/sysctl.h> // sysctl()
#endif

// class BinaryPlink;
// class BinaryGen;
#define MULTIPLEX_LD 1920
#define MULTIPLEX_2LD (MULTIPLEX_LD * 2)

class Genotype
{
public:
    Genotype(){};
    Genotype(const size_t thread, const bool ignore_fid,
             const bool keep_nonfounder, const bool keep_ambig,
             const bool is_ref = false)
        : m_is_ref(is_ref)
        , m_thread(thread)
        , m_keep_nonfounder(keep_nonfounder)
        , m_ignore_fid(ignore_fid)
        , m_keep_ambig(keep_ambig){};
    virtual ~Genotype();

    // after load samples, samples that we don't want to include
    // will all have their FID, IID and phenotype masked by empty string
    void load_samples(const std::string& keep_file,
                      const std::string& remove_file, bool verbose,
                      Reporter& reporter);
    // we use pointer for the target such that we can have a nullptr as default
    // therefore not providing any input argument for target file processing
    void load_snps(const std::string out_prefix,
                   const std::string& extract_file,
                   const std::string& exclude_file, const double geno,
                   const double maf, const double info,
                   const double hard_threshold, const bool hard_coded,
                   Region& exclusion, bool verbose, Reporter& reporter,
                   Genotype* target = nullptr);

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
                   size_t& num_snp_included, const size_t region_index, const bool cumulate,
                   const bool require_statistic);
    bool sort_by_p();
    void print_snp(std::string& output, double threshold,
                   const size_t region_index);
    size_t num_threshold() const { return m_num_threshold; };
    void read_base(const Commander& c_commander, Region& region,
                   Reporter& reporter);
    // void clump(Genotype& reference);
    void efficient_clumping(Genotype& reference, Reporter& reporter,
                            bool const pearson);
    void set_info(const Commander& c_commander, const bool ld = false);

    void reset_sample_pheno()
    {
        for (auto&& sample : m_sample_names) {
            sample.num_snp = 0;
            sample.prs = 0.0;
            sample.in_regression = false;
        }
    };

    void reset_sample_prs()
    {
        for (auto&& sample : m_sample_names) {
            sample.prs = 0.0;
            sample.num_snp = 0.0;
        }
    };
    bool prepare_prsice(Reporter& reporter);
    std::string sample_id(size_t i) const
    {
        if (i > m_sample_names.size())
            throw std::out_of_range("Sample name vector out of range");
        if (m_ignore_fid)
            return m_sample_names[i].IID;
        else
            return m_sample_names[i].FID + "_" + m_sample_names[i].IID;
    }


    bool sample_in_regression(size_t i) const
    {
        return m_sample_names.at(i).in_regression;
    }
    bool is_include(size_t i) const { return m_sample_names.at(i).include; }
    void set_in_regression(size_t i, bool within)
    {
        m_sample_names.at(i).in_regression = within;
    }
    bool include_for_regression(size_t i) const
    {
        return m_sample_names.at(i).in_regression;
    };
    std::string pheno(size_t i) const { return m_sample_names.at(i).pheno; }
    bool pheno_is_na(size_t i) const
    {
        return m_sample_names.at(i).pheno.compare("NA") == 0;
    }
    std::string fid(size_t i) const { return m_sample_names.at(i).FID; }
    std::string iid(size_t i) const { return m_sample_names.at(i).IID; }
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
    void update_index(size_t i)
    {
        if (i == m_existed_snps.size())
            m_cur_category_index = m_categories.size();
        else
            m_cur_category_index =
                m_categories_index[m_existed_snps[i].category()];
    }
    uintptr_t founder_ct() const { return m_founder_ct; }
    uintptr_t unfiltered_sample_ct() const { return m_unfiltered_sample_ct; }
    void count_snp_in_region(Region& region, const std::string& name,
                             const bool print_snp) const
    {
        std::string snp_file_name = name + ".snp";
        std::ofstream print_file;
        if (print_snp) {
            print_file.open(snp_file_name.c_str());
            if (!print_file.is_open()) {
                throw std::runtime_error(
                    "Error: Cannot open snp file to write: " + snp_file_name);
            }
            print_file << "SNP\tCHR\tBP\tP";
            for (size_t i_region = 0; i_region < region.size(); ++i_region) {
                print_file << "\t" << region.get_name(i_region);
            }
            print_file << std::endl;
        }
        std::vector<int> result(region.size(), 0);
        for (auto&& snp : m_existed_snps) {
            if (print_snp)
                print_file << snp.rs() << "\t" << snp.chr() << "\t" << snp.loc()
                           << "\t" << snp.p_value();
            for (size_t i_region = 0; i_region < region.size(); ++i_region) {
                result[i_region] += snp.in(i_region);
                if (print_snp)
                    print_file << "\t" << (snp.in(i_region) ? "Y" : "N");
            }
            if (print_snp) print_file << std::endl;
        }
        if (print_snp) print_file.close();
        region.post_clump_count(result);
    };

    void get_null_score(const size_t& set_size, const size_t& num_selected_snps,
                        const std::vector<size_t>& selection_list,
                        const bool require_standardize);
    void init_background_index(const size_t& background_index)
    {
        m_background_snp_index.clear();
        for (size_t i = 0; i < m_existed_snps.size(); ++i)
            if (m_existed_snps[i].in(background_index))
                m_background_snp_index.push_back(i);
    }
    size_t num_background() const { return m_background_snp_index.size(); };

protected:
    friend class BinaryPlink;
    friend class BinaryGen;
    // variable storages
    // vector storing all the genotype files
    std::vector<std::string> m_genotype_files;
    // sample file name. Fam for plink
    std::string m_sample_file;

    bool m_is_ref = false;
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
    uint32_t m_num_ref_target_mismatch = 0;
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
    std::vector<size_t> m_background_snp_index;
    // functions
    // function to substitute the # in the sample name
    std::vector<std::string> set_genotype_files(const std::string& prefix);
    std::vector<std::string> load_genotype_prefix(const std::string& file_name);
    void init_chr(int num_auto = 22, bool no_x = false, bool no_y = false,
                  bool no_xy = false, bool no_mt = false);
    // responsible for reading in the sample
    virtual std::vector<Sample> gen_sample_vector()
    {
        return std::vector<Sample>(0);
    };
    virtual std::vector<SNP>
    gen_snp_vector(const double geno, const double maf, const double info,
                   const double hard_threshold, const bool hard_coded,
                   Region& exclusion, const std::string& out_prefix,
                   Genotype* target = nullptr)
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
    unsigned int m_seed;
    size_t m_max_category = 0;
    size_t m_region_size = 1;
    size_t m_num_threshold = 0;
    size_t m_max_window_size = 0;
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
    std::vector<int> m_categories;
    std::unordered_map<int, int> m_categories_index;
    int m_cur_category_index = 0;

    uintptr_t m_unfiltered_marker_ct = 0;
    // uint32_t m_hh_exists;

    void update_snps(std::vector<int>& retained_index);
    virtual inline void read_genotype(uintptr_t* genotype,
                                      const std::streampos byte_pos,
                                      const std::string& file_name){};
    virtual void read_score(size_t start_index, size_t end_bound,
                            const size_t region_index){};
    virtual void read_score(std::vector<size_t>& index){};


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

    // mask_buf has the same size as the geno_buf and is used in later situation
    // missing_ct_ptr seems to be a int for number of missing
    // in plink, because they use multi-threading, they use a vector for this
    // but as we are clumping, we can't do multi-threading in the way
    // plink does (though might be possible if we start to incorporate the
    // highly complicated memory pool of plink)
    void ld_process_load2(uintptr_t* geno_buf, uintptr_t* mask_buf,
                          uint32_t* missing_ct_ptr, uint32_t founder_ct,
                          uint32_t is_x, uintptr_t* founder_male_include2)
    {
        // if is_x is always false, then we can safely ignore
        // founder_male_include2
        uintptr_t* geno_ptr = geno_buf;
        uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(founder_ct);
        uintptr_t* geno_end = &(geno_buf[founder_ctl2]);
        uintptr_t* mask_buf_ptr = mask_buf;
        uintptr_t cur_geno;
        uintptr_t shifted_masked_geno;
        uintptr_t new_geno;
        uintptr_t new_mask;
        do
        {
            cur_geno = *geno_ptr;
            shifted_masked_geno = (cur_geno >> 1) & FIVEMASK;
            new_geno = cur_geno - shifted_masked_geno;
            *geno_ptr++ = new_geno;
            new_mask = (((~cur_geno) & FIVEMASK) | shifted_masked_geno) * 3;
            *mask_buf_ptr++ = new_mask;
        } while (geno_ptr < geno_end);
        if (is_x) {
            geno_ptr = geno_buf;
            do
            {
                new_geno = *geno_ptr;
                *geno_ptr++ = new_geno
                              + ((~(new_geno | (new_geno >> 1)))
                                 & (*founder_male_include2++));
            } while (geno_ptr < geno_end);
        }
        if (founder_ct % BITCT2) {
            mask_buf[founder_ct / BITCT2] &=
                (ONELU << (2 * (founder_ct % BITCT2))) - ONELU;
        }
        *missing_ct_ptr =
            founder_ct - (popcount_longs(mask_buf, founder_ctl2) / 2);
    }


// amazing bitwise R2 calculation function from PLINK
#ifdef __LP64__
    static inline void ld_dot_prod_batch(__m128i* vec1, __m128i* vec2,
                                         __m128i* mask1, __m128i* mask2,
                                         int32_t* return_vals, uint32_t iters)
    {
        // Main routine for computation of \sum_i^M (x_i - \mu_x)(y_i - \mu_y),
        // where x_i, y_i \in \{-1, 0, 1\}, but there are missing values.
        //
        //
        // We decompose this sum into
        //   \sum_i x_iy_i - \mu_y\sum_i x_i - \mu_x\sum_i y_i +
        //   (M - # missing)\mu_x\mu_y.
        // *Without* missing values, this can be handled very cleanly.  The last
        // three terms can all be precomputed, and \sum_i x_iy_i can be handled
        // in a manner very similar to bitwise Hamming distance.  This is
        // several times as fast as the lookup tables used for relationship
        // matrices.
        //
        // Unfortunately, when missing values are present,
        // \mu_y\sum_{i: nonmissing from y} x_i and
        // \mu_x\sum_{i: nonmissing from x} y_i must also be evaluated (and, in
        // practice, \mu_y\sum_{i: nonmissing from y} x_i^2 and
        // \mu_x\sum_{i: nonmissing from x} y_i^2 should be determined here as
        // well); this removes much of the speed advantage, and the best
        // applications of the underlying ternary dot product algorithm used
        // here lie elsewhere. Nevertheless, it is still faster, so we use it.
        // (possible todo: accelerated function when there really are no missing
        // values, similar to what is now done for --fast-epistasis)
        //
        //
        // Input:
        // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
        // * mask1 and mask2 mask out missing values (i.e. 00 for missing, 11
        // for
        //   nonmissing).
        // * return_vals provides space for return values.
        // * iters is the number of 48-byte windows to process, anywhere from 1
        // to 10
        //   inclusive.
        //
        // This function performs the update
        //   return_vals[0] += (-N) + \sum_i x_iy_i
        //   return_vals[1] += N_y + \sum_{i: nonmissing from y} x_i
        //   return_vals[2] += N_x + \sum_{i: nonmissing from x} y_i
        //   return_vals[3] += N_y - \sum_{i: nonmissing from y} x_i^2
        //   return_vals[4] += N_x - \sum_{i: nonmissing from x} y_i^2
        // where N is the number of samples processed after applying the
        // missingness masks indicated by the subscripts.
        //
        // Computation of terms [1]-[4] is based on the identity
        //   N_y + \sum_{i: nonmissing from y} x_i = popcount2(vec1 & mask2)
        // where "popcount2" refers to starting with two-bit integers instead of
        // one-bit integers in our summing process (this allows us to skip a few
        // operations).  (Once we can assume the presence of hardware popcount,
        // a slightly different implementation may be better.)
        //
        // The trickier [0] computation currently proceeds as follows:
        //
        // 1. zcheck := (vec1 | vec2) & 0x5555...
        // Detects whether at least one member of the pair has a 0/missing
        // value.
        //
        // 2. popcount2(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
        // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i
        // dot product.
        //
        // MULTIPLEX_LD sets of values are usually handled per function call. If
        // fewer values are present, the ends of all input vectors should be
        // zeroed out.

        const __m128i m1 = {FIVEMASK, FIVEMASK};
        const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
        const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
        __m128i loader1;
        __m128i loader2;
        __m128i sum1;
        __m128i sum2;
        __m128i sum11;
        __m128i sum22;
        __m128i sum12;
        __m128i tmp_sum1;
        __m128i tmp_sum2;
        __m128i tmp_sum12;
        __univec acc;
        __univec acc1;
        __univec acc2;
        __univec acc11;
        __univec acc22;
        acc.vi = _mm_setzero_si128();
        acc1.vi = _mm_setzero_si128();
        acc2.vi = _mm_setzero_si128();
        acc11.vi = _mm_setzero_si128();
        acc22.vi = _mm_setzero_si128();
        do
        {
            loader1 = *vec1++;
            loader2 = *vec2++;
            sum1 = *mask2++;
            sum2 = *mask1++;
            sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
            // sum11 = _mm_and_si128(_mm_and_si128(_mm_xor_si128(sum1, m1), m1),
            // loader1); sum22 = _mm_and_si128(_mm_and_si128(_mm_xor_si128(sum2,
            // m1), m1), loader2);
            sum1 = _mm_and_si128(sum1, loader1);
            sum2 = _mm_and_si128(sum2, loader2);
            sum11 = _mm_and_si128(sum1, m1);
            sum22 = _mm_and_si128(sum2, m1);
            // use andnot to eliminate need for 0xaaaa... to occupy an xmm
            // register
            loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12),
                                       _mm_xor_si128(loader1, loader2));
            sum12 = _mm_or_si128(sum12, loader1);

            // sum1, sum2, and sum12 now store the (biased) two-bit sums of
            // interest; merge to 4 bits to prevent overflow.  this merge can be
            // postponed for sum11 and sum22 because the individual terms are
            // 0/1 instead of 0/1/2.
            sum1 = _mm_add_epi64(_mm_and_si128(sum1, m2),
                                 _mm_and_si128(_mm_srli_epi64(sum1, 2), m2));
            sum2 = _mm_add_epi64(_mm_and_si128(sum2, m2),
                                 _mm_and_si128(_mm_srli_epi64(sum2, 2), m2));
            sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2),
                                  _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

            loader1 = *vec1++;
            loader2 = *vec2++;
            tmp_sum1 = *mask2++;
            tmp_sum2 = *mask1++;
            tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
            tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
            tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
            sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
            sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
            loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12),
                                       _mm_xor_si128(loader1, loader2));
            tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

            sum1 = _mm_add_epi64(
                sum1,
                _mm_add_epi64(_mm_and_si128(tmp_sum1, m2),
                              _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
            sum2 = _mm_add_epi64(
                sum2,
                _mm_add_epi64(_mm_and_si128(tmp_sum2, m2),
                              _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
            sum12 = _mm_add_epi64(
                sum12,
                _mm_add_epi64(_mm_and_si128(tmp_sum12, m2),
                              _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

            loader1 = *vec1++;
            loader2 = *vec2++;
            tmp_sum1 = *mask2++;
            tmp_sum2 = *mask1++;
            tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
            tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
            tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
            sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
            sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
            loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12),
                                       _mm_xor_si128(loader1, loader2));
            tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

            sum1 = _mm_add_epi64(
                sum1,
                _mm_add_epi64(_mm_and_si128(tmp_sum1, m2),
                              _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
            sum2 = _mm_add_epi64(
                sum2,
                _mm_add_epi64(_mm_and_si128(tmp_sum2, m2),
                              _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
            sum11 = _mm_add_epi64(_mm_and_si128(sum11, m2),
                                  _mm_and_si128(_mm_srli_epi64(sum11, 2), m2));
            sum22 = _mm_add_epi64(_mm_and_si128(sum22, m2),
                                  _mm_and_si128(_mm_srli_epi64(sum22, 2), m2));
            sum12 = _mm_add_epi64(
                sum12,
                _mm_add_epi64(_mm_and_si128(tmp_sum12, m2),
                              _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

            acc1.vi = _mm_add_epi64(
                acc1.vi,
                _mm_add_epi64(_mm_and_si128(sum1, m4),
                              _mm_and_si128(_mm_srli_epi64(sum1, 4), m4)));
            acc2.vi = _mm_add_epi64(
                acc2.vi,
                _mm_add_epi64(_mm_and_si128(sum2, m4),
                              _mm_and_si128(_mm_srli_epi64(sum2, 4), m4)));
            acc11.vi = _mm_add_epi64(
                acc11.vi,
                _mm_add_epi64(_mm_and_si128(sum11, m4),
                              _mm_and_si128(_mm_srli_epi64(sum11, 4), m4)));
            acc22.vi = _mm_add_epi64(
                acc22.vi,
                _mm_add_epi64(_mm_and_si128(sum22, m4),
                              _mm_and_si128(_mm_srli_epi64(sum22, 4), m4)));
            acc.vi = _mm_add_epi64(
                acc.vi,
                _mm_add_epi64(_mm_and_si128(sum12, m4),
                              _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
        } while (--iters);
        // moved down because we've almost certainly run out of xmm registers
        const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
#if MULTIPLEX_LD > 960
        acc1.vi = _mm_add_epi64(_mm_and_si128(acc1.vi, m8),
                                _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
        acc2.vi = _mm_add_epi64(_mm_and_si128(acc2.vi, m8),
                                _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
        acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                               _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
        acc1.vi = _mm_and_si128(
            _mm_add_epi64(acc1.vi, _mm_srli_epi64(acc1.vi, 8)), m8);
        acc2.vi = _mm_and_si128(
            _mm_add_epi64(acc2.vi, _mm_srli_epi64(acc2.vi, 8)), m8);
        acc.vi =
            _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
        acc11.vi = _mm_and_si128(
            _mm_add_epi64(acc11.vi, _mm_srli_epi64(acc11.vi, 8)), m8);
        acc22.vi = _mm_and_si128(
            _mm_add_epi64(acc22.vi, _mm_srli_epi64(acc22.vi, 8)), m8);

        return_vals[0] -= ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
        return_vals[1] +=
            ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
        return_vals[2] +=
            ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
        return_vals[3] +=
            ((acc11.u8[0] + acc11.u8[1]) * 0x1000100010001LLU) >> 48;
        return_vals[4] +=
            ((acc22.u8[0] + acc22.u8[1]) * 0x1000100010001LLU) >> 48;
    }

    void ld_dot_prod(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1,
                     uintptr_t* mask2, int32_t* return_vals,
                     uint32_t batch_ct_m1, uint32_t last_batch_size)
    {
        while (batch_ct_m1--) {
            ld_dot_prod_batch((__m128i*) vec1, (__m128i*) vec2,
                              (__m128i*) mask1, (__m128i*) mask2, return_vals,
                              MULTIPLEX_LD / 192);
            vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
            vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
            mask1 = &(mask1[MULTIPLEX_LD / BITCT2]);
            mask2 = &(mask2[MULTIPLEX_LD / BITCT2]);
        }
        ld_dot_prod_batch((__m128i*) vec1, (__m128i*) vec2, (__m128i*) mask1,
                          (__m128i*) mask2, return_vals, last_batch_size);
    }

    static inline int32_t ld_dot_prod_nm_batch(__m128i* vec1, __m128i* vec2,
                                               uint32_t iters)
    {
        // faster ld_dot_prod_batch() for no-missing-calls case.
        const __m128i m1 = {FIVEMASK, FIVEMASK};
        const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
        const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
        const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
        __m128i loader1;
        __m128i loader2;
        __m128i sum12;
        __m128i tmp_sum12;
        __univec acc;
        acc.vi = _mm_setzero_si128();
        do
        {
            loader1 = *vec1++;
            loader2 = *vec2++;
            sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
            loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12),
                                       _mm_xor_si128(loader1, loader2));
            sum12 = _mm_or_si128(sum12, loader1);
            sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2),
                                  _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

            loader1 = *vec1++;
            loader2 = *vec2++;
            tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
            loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12),
                                       _mm_xor_si128(loader1, loader2));
            tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);
            sum12 = _mm_add_epi64(
                sum12,
                _mm_add_epi64(_mm_and_si128(tmp_sum12, m2),
                              _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

            loader1 = *vec1++;
            loader2 = *vec2++;
            tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
            loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12),
                                       _mm_xor_si128(loader1, loader2));
            tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);
            sum12 = _mm_add_epi64(
                sum12,
                _mm_add_epi64(_mm_and_si128(tmp_sum12, m2),
                              _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

            acc.vi = _mm_add_epi64(
                acc.vi,
                _mm_add_epi64(_mm_and_si128(sum12, m4),
                              _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
        } while (--iters);
#if MULTIPLEX_LD > 960
        acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                               _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
        acc.vi =
            _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
        return (uint32_t)(((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48);
    }

    int32_t ld_dot_prod_nm(uintptr_t* vec1, uintptr_t* vec2,
                           uint32_t founder_ct, uint32_t batch_ct_m1,
                           uint32_t last_batch_size)
    {
        // accelerated implementation for no-missing-loci case
        int32_t result = (int32_t) founder_ct;
        while (batch_ct_m1--) {
            result -= ld_dot_prod_nm_batch((__m128i*) vec1, (__m128i*) vec2,
                                           MULTIPLEX_LD / 192);
            vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
            vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
        }
        result -= ld_dot_prod_nm_batch((__m128i*) vec1, (__m128i*) vec2,
                                       last_batch_size);
        return result;
    }
#else
    static inline void ld_dot_prod_batch(uintptr_t* vec1, uintptr_t* vec2,
                                         uintptr_t* mask1, uintptr_t* mask2,
                                         int32_t* return_vals, uint32_t iters)
    {
        uint32_t final_sum1 = 0;
        uint32_t final_sum2 = 0;
        uint32_t final_sum11 = 0;
        uint32_t final_sum22 = 0;
        uint32_t final_sum12 = 0;
        uintptr_t loader1;
        uintptr_t loader2;
        uintptr_t sum1;
        uintptr_t sum2;
        uintptr_t sum11;
        uintptr_t sum22;
        uintptr_t sum12;
        uintptr_t tmp_sum1;
        uintptr_t tmp_sum2;
        uintptr_t tmp_sum12;
        do
        {
            // (The important part of the header comment on the 64-bit version
            // is copied below.)
            //
            // Input:
            // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
            // * mask1 and mask2 mask out missing values (i.e. 00 for missing,
            // 11 for
            //   nonmissing).
            // * return_vals provides space for return values.
            // * iters is the number of 12-byte windows to process, anywhere
            // from 1 to
            //   40 inclusive.  (No, this is not the interface you'd use for a
            //   general-purpose library.)  [32- and 64-bit differ here.]
            //
            // This function performs the update
            //   return_vals[0] += (-N) + \sum_i x_iy_i
            //   return_vals[1] += N_y + \sum_{i: nonmissing from y} x_i
            //   return_vals[2] += N_x + \sum_{i: nonmissing from x} y_i
            //   return_vals[3] += N_y - \sum_{i: nonmissing from y} x_i^2
            //   return_vals[4] += N_x - \sum_{i: nonmissing from x} y_i^2
            // where N is the number of samples processed after applying the
            // missingness masks indicated by the subscripts.
            //
            // Computation of terms [1]-[4] is based on the identity
            //   N_y + \sum_{i: nonmissing from y} x_i = popcount2(vec1 & mask2)
            // where "popcount2" refers to starting with two-bit integers
            // instead of one-bit integers in our summing process (this allows
            // us to skip a few operations).  (Once we can assume the presence
            // of hardware popcount, a slightly different implementation may be
            // better.)
            //
            // The trickier [0] computation currently proceeds as follows:
            //
            // 1. zcheck := (vec1 | vec2) & 0x5555...
            // Detects whether at least one member of the pair has a 0/missing
            // value.
            //
            // 2. popcount2(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
            // Subtracting this *from* a bias will give us our desired \sum_i
            // x_iy_i dot product.

            loader1 = *vec1++;
            loader2 = *vec2++;
            sum1 = *mask2++;
            sum2 = *mask1++;
            sum12 = (loader1 | loader2) & FIVEMASK;

            sum1 = sum1 & loader1;
            sum2 = sum2 & loader2;
            loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
            sum12 = sum12 | loader1;
            sum11 = sum1 & FIVEMASK;
            sum22 = sum2 & FIVEMASK;

            sum1 = (sum1 & 0x33333333) + ((sum1 >> 2) & 0x33333333);
            sum2 = (sum2 & 0x33333333) + ((sum2 >> 2) & 0x33333333);
            sum12 = (sum12 & 0x33333333) + ((sum12 >> 2) & 0x33333333);

            loader1 = *vec1++;
            loader2 = *vec2++;
            tmp_sum1 = *mask2++;
            tmp_sum2 = *mask1++;
            tmp_sum12 = (loader1 | loader2) & FIVEMASK;
            tmp_sum1 = tmp_sum1 & loader1;
            tmp_sum2 = tmp_sum2 & loader2;

            loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
            tmp_sum12 = tmp_sum12 | loader1;
            sum11 += tmp_sum1 & FIVEMASK;
            sum22 += tmp_sum2 & FIVEMASK;

            sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
            sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
            sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

            loader1 = *vec1++;
            loader2 = *vec2++;
            tmp_sum1 = *mask2++;
            tmp_sum2 = *mask1++;
            tmp_sum12 = (loader1 | loader2) & FIVEMASK;

            tmp_sum1 = tmp_sum1 & loader1;
            tmp_sum2 = tmp_sum2 & loader2;
            loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
            tmp_sum12 = tmp_sum12 | loader1;
            sum11 += tmp_sum1 & FIVEMASK;
            sum22 += tmp_sum2 & FIVEMASK;

            sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
            sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
            sum11 = (sum11 & 0x33333333) + ((sum11 >> 2) & 0x33333333);
            sum22 = (sum22 & 0x33333333) + ((sum22 >> 2) & 0x33333333);
            sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

            sum1 = (sum1 & 0x0f0f0f0f) + ((sum1 >> 4) & 0x0f0f0f0f);
            sum2 = (sum2 & 0x0f0f0f0f) + ((sum2 >> 4) & 0x0f0f0f0f);
            sum11 = (sum11 & 0x0f0f0f0f) + ((sum11 >> 4) & 0x0f0f0f0f);
            sum22 = (sum22 & 0x0f0f0f0f) + ((sum22 >> 4) & 0x0f0f0f0f);
            sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);

            // technically could do the multiply-and-shift only once every two
            // rounds
            final_sum1 += (sum1 * 0x01010101) >> 24;
            final_sum2 += (sum2 * 0x01010101) >> 24;
            final_sum11 += (sum11 * 0x01010101) >> 24;
            final_sum22 += (sum22 * 0x01010101) >> 24;
            final_sum12 += (sum12 * 0x01010101) >> 24;
        } while (--iters);
        return_vals[0] -= final_sum12;
        return_vals[1] += final_sum1;
        return_vals[2] += final_sum2;
        return_vals[3] += final_sum11;
        return_vals[4] += final_sum22;
    }

    void ld_dot_prod(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1,
                     uintptr_t* mask2, int32_t* return_vals,
                     uint32_t batch_ct_m1, uint32_t last_batch_size)
    {
        while (batch_ct_m1--) {
            ld_dot_prod_batch(vec1, vec2, mask1, mask2, return_vals,
                              MULTIPLEX_LD / 48);
            vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
            vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
            mask1 = &(mask1[MULTIPLEX_LD / BITCT2]);
            mask2 = &(mask2[MULTIPLEX_LD / BITCT2]);
        }
        ld_dot_prod_batch(vec1, vec2, mask1, mask2, return_vals,
                          last_batch_size);
    }

    static inline int32_t ld_dot_prod_nm_batch(uintptr_t* vec1, uintptr_t* vec2,
                                               uint32_t iters)
    {
        uint32_t final_sum12 = 0;
        uintptr_t loader1;
        uintptr_t loader2;
        uintptr_t sum12;
        uintptr_t tmp_sum12;
        do
        {
            loader1 = *vec1++;
            loader2 = *vec2++;
            sum12 = (loader1 | loader2) & FIVEMASK;
            loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
            sum12 = sum12 | loader1;
            sum12 = (sum12 & 0x33333333) + ((sum12 >> 2) & 0x33333333);

            loader1 = *vec1++;
            loader2 = *vec2++;
            tmp_sum12 = (loader1 | loader2) & FIVEMASK;
            loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
            tmp_sum12 = tmp_sum12 | loader1;
            sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

            loader1 = *vec1++;
            loader2 = *vec2++;
            tmp_sum12 = (loader1 | loader2) & FIVEMASK;
            loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
            tmp_sum12 = tmp_sum12 | loader1;
            sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);
            sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);

            final_sum12 += (sum12 * 0x01010101) >> 24;
        } while (--iters);
        return final_sum12;
    }

    int32_t ld_dot_prod_nm(uintptr_t* vec1, uintptr_t* vec2,
                           uint32_t founder_ct, uint32_t batch_ct_m1,
                           uint32_t last_batch_size)
    {
        int32_t result = (int32_t) founder_ct;
        while (batch_ct_m1--) {
            result -= ld_dot_prod_nm_batch(vec1, vec2, MULTIPLEX_LD / 48);
            vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
            vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
        }
        result -= ld_dot_prod_nm_batch(vec1, vec2, last_batch_size);
        return result;
    }
#endif // __LP64__


    uint32_t ld_missing_ct_intersect(uintptr_t* lptr1, uintptr_t* lptr2,
                                     uintptr_t word12_ct, uintptr_t word12_rem,
                                     uintptr_t lshift_last)
    {
        // variant of popcount_longs_intersect()
        uintptr_t tot = 0;
        uintptr_t* lptr1_end2;
#ifdef __LP64__
        const __m128i m1 = {FIVEMASK, FIVEMASK};
        const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
        const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
        const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
        __m128i* vptr1 = (__m128i*) lptr1;
        __m128i* vptr2 = (__m128i*) lptr2;
        __m128i* vend1;
        __m128i loader1;
        __m128i loader2;
        __univec acc;

        while (word12_ct >= 10) {
            word12_ct -= 10;
            vend1 = &(vptr1[60]);
        ld_missing_ct_intersect_main_loop:
            acc.vi = _mm_setzero_si128();
            do
            {
                loader1 =
                    _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
                loader2 =
                    _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
                loader1 = _mm_add_epi64(
                    loader1,
                    _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
                loader2 = _mm_add_epi64(
                    loader2,
                    _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
                loader1 = _mm_add_epi64(
                    loader1,
                    _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
                loader2 = _mm_add_epi64(
                    loader2,
                    _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
                loader1 = _mm_add_epi64(
                    _mm_and_si128(loader1, m2),
                    _mm_and_si128(_mm_srli_epi64(loader1, 2), m2));
                loader1 = _mm_add_epi64(
                    loader1,
                    _mm_add_epi64(
                        _mm_and_si128(loader2, m2),
                        _mm_and_si128(_mm_srli_epi64(loader2, 2), m2)));
                acc.vi = _mm_add_epi64(
                    acc.vi, _mm_add_epi64(
                                _mm_and_si128(loader1, m4),
                                _mm_and_si128(_mm_srli_epi64(loader1, 4), m4)));
            } while (vptr1 < vend1);
            acc.vi =
                _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                              _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
            tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
        }
        if (word12_ct) {
            vend1 = &(vptr1[word12_ct * 6]);
            word12_ct = 0;
            goto ld_missing_ct_intersect_main_loop;
        }
        lptr1 = (uintptr_t*) vptr1;
        lptr2 = (uintptr_t*) vptr2;
#else
        uintptr_t* lptr1_end = &(lptr1[word12_ct * 12]);
        uintptr_t tmp_stor;
        uintptr_t loader1;
        uintptr_t loader2;
        while (lptr1 < lptr1_end) {
            loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
            loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
            tmp_stor = (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);

            loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
            loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
            loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
            tmp_stor += (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);
            tot += (tmp_stor * 0x01010101) >> 24;
        }
#endif
        lptr1_end2 = &(lptr1[word12_rem]);
        while (lptr1 < lptr1_end2) {
            tot += popcount2_long((~((*lptr1++) | (*lptr2++))) & FIVEMASK);
        }
        if (lshift_last) {
            tot += popcount2_long(((~((*lptr1) | (*lptr2))) & FIVEMASK)
                                  << lshift_last);
        }
        return tot;
    }
};

#endif /* SRC_GENOTYPE_HPP_ */
