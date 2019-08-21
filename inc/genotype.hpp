// This file is part of PRSice-2, copyright (C) 2016-2019
// Shing Wan Choi, Paul F. O’Reilly
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

#include "IITree.h"
#include "commander.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
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
#include <functional>
#include <memory>
#include <mio.hpp>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef __APPLE__
#include <sys/sysctl.h> // sysctl()
#endif


#ifdef _WIN32
#include <mingw.mutex.h>
#include <mingw.thread.h>
#else
#include <mutex>
#include <thread>
#endif
// class BinaryPlink;
// class BinaryGen;
#define MULTIPLEX_LD 1920
#define MULTIPLEX_2LD (MULTIPLEX_LD * 2)
class Genotype
{
public:
    /*!
     * \brief default constructor of Genotype
     */
    Genotype() {}
    /*!
     * \brief Genotype constructor
     * \param thread = number of thread
     * \param ignore_fid = if FID should be ignored
     * \param keep_nonfounder = if we should keep non-founders
     * \param keep_ambig = if we should keep ambiguous SNPs
     * \param is_ref = if this is a reference genome
     */
    Genotype(const uint32_t thread, const bool ignore_fid,
             const bool keep_nonfounder, const bool keep_ambig,
             const bool is_ref = false)
        : m_thread(thread)
        , m_ignore_fid(ignore_fid)
        , m_is_ref(is_ref)
        , m_keep_nonfounder(keep_nonfounder)
        , m_keep_ambig(keep_ambig)
    {
    }
    virtual ~Genotype();

    // after load samples, samples that we don't want to include
    // will all have their FID, IID and phenotype masked by empty string
    /*!
     * \brief Function to load samples into the genotype class object. Will call
     *        gen_sample_vector. Will only store the sample info if this is the
     *        target.
     * \param keep_file is the file containing samples to be included
     * \param remove_file is the file containing samples to be excluded
     * \param verbose boolean indicate whether we should output the count
     * \param reporter is the logger
     */
    void load_samples(const std::string& keep_file,
                      const std::string& remove_file, const std::string& delim,
                      bool verbose);

    // We need the exclusion_region parameter because when we read in the base
    // we do allow users to provide a base without the CHR and LOC, which forbid
    // us to do the regional filtering. However, as exclusion and extractions
    // are based on the SNP ID, we can be confident that they were already done
    // in read_base
    void load_snps(const std::string& out,
                   const std::vector<IITree<size_t, size_t>>& exclusion_regions,
                   bool verbose, Genotype* target = nullptr);

    void calc_freqs_and_intermediate(
        const double& maf_threshold, const double& geno_threshold,
        const double& info_threshold, const bool maf_filter,
        const bool geno_filter, const bool info_filter, const bool hard_coded,
        bool verbose, Genotype* target = nullptr);
    /*!
     * \brief Return the number of SNPs, use for unit test
     * \return reuturn the number of SNPs included
     */
    std::vector<SNP>::size_type num_snps() const
    {
        return m_existed_snps.size();
    }
    /*!
     * \brief Function to re-propagate the m_existed_snps_index after reading
     * the reference panel
     */
    void update_snp_index()
    {
        m_existed_snps_index.clear();
        for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp)
        { m_existed_snps_index[m_existed_snps[i_snp].rs()] = i_snp; }
    }
    /*!
     * \brief Return the number of sample we wish to perform PRS on
     * \return the number of sample
     */
    size_t num_sample() const { return m_sample_id.size(); }

    /*!
     * \brief Function to prepare clumping. Should sort all the SNPs by their
     * chromosome number and then by their p-value
     * \return true if there are SNPs to sort
     */
    bool sort_by_p()
    {
        if (m_existed_snps.size() == 0) return false;
        m_sort_by_p_index = SNP::sort_by_p_chr(m_existed_snps);
        return true;
    }
    /*!
     * \brief Function to carry out clumping. One of the most complicated
     * function.
     *
     * \param reference is the reference genotype object. Will need to use it
     * for clumping. When target is used for clumping, reference=this
     *
     * \param reporter the logger
     * \param use_pearson indicate if we want to perform pearson correlation
     * calculation instead of the haplotype likelihood clumping. If this is
     * true, efficient_clumping will call the pearson_clumping algorithm
     */
    void efficient_clumping(Genotype& reference);
    /*!
     * \brief This function helps to load all command line dependencies into the
     * object so that we don't need to pass along the commander any more
     *
     * \param c_commander the container containing the required information
     */
    void set_info(const Commander& c_commander)
    {
        m_clump_p = c_commander.clump_p();
        m_clump_r2 = c_commander.clump_r2();
        m_use_proxy = c_commander.proxy(m_clump_proxy);
        m_clump_distance = static_cast<uintptr_t>(c_commander.clump_dist());
        m_model = c_commander.model();
        m_missing_score = c_commander.get_missing_score();
        m_scoring = c_commander.get_score();
        m_seed = c_commander.seed();
        switch (m_model)
        {
        case MODEL::HETEROZYGOUS:
            m_homcom_weight = 0;
            m_het_weight = 1;
            m_homrar_weight = 0;
            break;
        case MODEL::DOMINANT:
            m_homcom_weight = 0;
            m_het_weight = 1;
            m_homrar_weight = 1;
            break;
        case MODEL::RECESSIVE:
            m_homcom_weight = 0;
            m_het_weight = 0;
            m_homrar_weight = 1;
            break;
        default:
            m_homcom_weight = 0;
            m_het_weight = 1;
            m_homrar_weight = 2;
            break;
        }
    }

    /*!
     * \brief Before each run of PRSice, we need to reset the in regression flag
     * to false and propagate it later on to indicate if the sample is used in
     * the regression model
     */
    void reset_in_regression_flag()
    {
        std::fill(m_in_regression.begin(), m_in_regression.end(), 0);
    }
    /*!
     * \brief Function to prepare the object for PRSice. Will sort the
     * m_existed_snp vector according to their p-value.
     * \return True if there are SNPs to process
     */
    bool prepare_prsice();
    /*!
     * \brief This function will return the sample ID
     * \param i is the index of the sample
     * \return the sample ID
     */
    std::string sample_id(const size_t i, const std::string& delim) const
    {
        if (i > m_sample_id.size())
            throw std::out_of_range("Sample name vector out of range");
        if (m_ignore_fid)
            return m_sample_id[i].IID;
        else
            return m_sample_id[i].FID + delim + m_sample_id[i].IID;
    }

    /*!
     * \brief Funtion return whether sample is founder (whether sample should be
     * included in regression)
     *
     * \param i is the sample index
     * \return true if sample is to be included
     */
    bool is_founder(size_t i) const { return m_sample_id.at(i).founder; }

    // as we don't check i is within range, the following 3 functions
    // have a potential of array out of bound
    /*!
     * \brief Check if the i th sample has been included in regression
     * \param i is the index of the sample
     * \return true if it is included
     */
    bool sample_in_regression(size_t i) const
    {
        return IS_SET(m_in_regression.data(), i);
    }
    /*!
     * \brief Inform PRSice that the i th sample has been included in the
     * regression
     *
     * \param i is the sample index
     */
    void set_in_regression(size_t i) { SET_BIT(i, m_in_regression.data()); }
    /*!
     * \brief Return the phenotype stored in the fam file of the i th sample
     * \param i is the index of the  sample
     * \return the phenotype of the sample
     */
    std::string pheno(size_t i) const { return m_sample_id[i].pheno; }
    /*!
     * \brief This function return if the i th sample has NA as phenotype
     * \param i is the sample ID
     * \return true if the phenotype is NA
     */
    bool pheno_is_na(size_t i) const { return m_sample_id[i].pheno == "NA"; }
    /*!
     * \brief Return the fid of the i th sample
     * \param i is the sample index
     * \return FID of the i th sample
     */
    std::string fid(size_t i) const { return m_sample_id.at(i).FID; }
    /*!
     * \brief Return the iid fo the i th sample
     * \param i is the sample index
     * \return IID of the i th sample
     */
    std::string iid(size_t i) const { return m_sample_id.at(i).IID; }
    /*!
     * \brief This function will calculate the required PRS for the i th sample
     * \param score_type is the type of score user want to calculate
     * \param i is the sample index
     * \return the PRS
     */
    inline double calculate_score(SCORING score_type, size_t i) const
    {
        if (i >= m_prs_info.size())
            throw std::out_of_range("Sample name vector out of range");
        double prs = m_prs_info[i].prs;
        int num_snp = m_prs_info[i].num_snp;
        double avg = prs;
        if (num_snp == 0) { avg = 0.0; }
        else
        {
            avg = prs / static_cast<double>(num_snp);
        }

        switch (score_type)
        {
        case SCORING::SUM: return m_prs_info[i].prs;
        case SCORING::STANDARDIZE: return (avg - m_mean_score) / m_score_sd;
        default:
            // default is avg
            return avg;
        }
    }
    /*!
     * \brief Function for calculating the PRS from the null set
     * \param set_size is the size of the set
     * \param prev_size is the amount of SNPs we have already processed
     * \param background_list is the vector containing the permuted background
     * index
     * \param first_run is a boolean representing if we need to reset the PRS to
     * 0
     * \param require_standardize is a boolean representing if we need to
     * calculate the mean and SD
     */
    void get_null_score(const size_t& set_size, const size_t& prev_size,
                        std::vector<size_t>& background_list,
                        const bool first_run, const bool require_standardize,
                        const bool use_ref_maf);
    /*!
     * \brief return the largest chromosome allowed
     * \return  the largest chromosome
     */
    uint32_t max_chr() const { return m_max_code; }
    /*!
     * \brief Indicate we will be using a reference file. Use for bgen
     * intermediate output generation
     */
    void expect_reference() { m_expect_reference = true; }
    void snp_extraction(const std::string& extract_snps,
                        const std::string& exclude_snps);
    void read_base(const std::string& base_file,
                   const std::vector<size_t>& col_index,
                   const std::vector<bool>& has_col,
                   const std::vector<double>& barlevels,
                   const double& bound_start, const double& bound_inter,
                   const double& bound_end,
                   const std::vector<IITree<size_t, size_t>>& exclusion_region,
                   const double& maf_control, const double& maf_case,
                   const double& info_threshold, const bool fastscore,
                   const bool no_full, const bool is_beta, const bool is_index,
                   const bool keep_ambig);
    void build_clump_windows();
    void build_membership_matrix(std::vector<size_t>& region_membership,
                                 std::vector<size_t>& region_start_idx,
                                 const size_t num_sets, const std::string& out,
                                 const std::vector<std::string>& region_name,
                                 const bool print_snps);
    size_t num_threshold() const { return m_num_thresholds; }
    std::vector<double> get_thresholds() const { return m_thresholds; }
    bool get_score(std::vector<size_t>::const_iterator& start_index,
                   const std::vector<size_t>::const_iterator& end_index,
                   double& cur_threshold, uint32_t& num_snp_included,
                   const bool non_cumulate, const bool require_statistic,
                   const bool first_run, const bool use_ref_maf);
    static bool within_region(const std::vector<IITree<size_t, size_t>>& cr,
                              const size_t chr, const size_t loc)
    {
        std::vector<size_t> output;
        if (chr >= cr.size()) return false;
        size_t start = (loc == 0) ? 1 : loc;
        cr[chr].overlap(start, loc + 1, output);
        return !output.empty();
    }


    static void construct_flag(
        const std::string& rs,
        const std::vector<IITree<size_t, size_t>>& gene_sets,
        const std::unordered_map<std::string, std::vector<size_t>>& snp_in_sets,
        std::vector<uintptr_t>& flag, const size_t required_size,
        const size_t chr, const size_t bp, const bool genome_wide_background)
    {
        if (flag.size() != required_size) { flag.resize(required_size); }
        std::fill(flag.begin(), flag.end(), 0);
        SET_BIT(0, flag.data());
        if (genome_wide_background) { SET_BIT(1, flag.data()); }
        // because the chromosome number is undefined. It will not be presented
        // in any of the region (we filter out any region with undefined chr)
        if (!gene_sets.empty())
        {
            std::vector<size_t> out;
            if (chr >= gene_sets.size()) return;
            gene_sets[chr].overlap(bp - 1, bp + 1, out);
            size_t idx;
            for (auto&& j : out)
            {
                idx = gene_sets[chr].data(j);
                // idx= cr_label(gene_sets, b[j]);
                SET_BIT(idx, flag.data());
            }
        }
        if (snp_in_sets.empty() || rs.empty()) return;
        auto&& snp_idx = snp_in_sets.find(rs);
        if (snp_idx != snp_in_sets.end())
        {
            for (auto&& i : snp_idx->second) { SET_BIT(i, flag.data()); }
        }
        return;
    }
    void add_flags(
        const std::vector<IITree<size_t, size_t>>& cr,
        const std::unordered_map<std::string, std::vector<size_t>>& snp_in_sets,
        const size_t num_sets, const bool genome_wide_background);

    virtual void init_mmap() {}

protected:
    // friend with all child class so that they can also access the
    // protected elements
    friend class BinaryPlink;
    friend class BinaryGen;
    // vector storing all the genotype files
    // std::vector<Sample> m_sample_names;
    std::vector<SNP> m_existed_snps;
    std::unordered_map<std::string, size_t> m_existed_snps_index;
    std::unordered_set<std::string> m_sample_selection_list;
    std::unordered_set<std::string> m_snp_selection_list;
    std::vector<Sample_ID> m_sample_id;
    std::vector<PRS> m_prs_info;
    std::vector<std::string> m_genotype_file_names;
    std::vector<mio::mmap_source> m_genotype_files;
    std::vector<double> m_thresholds;
    std::vector<uintptr_t> m_tmp_genotype;
    // std::vector<uintptr_t> m_chrom_mask;
    std::vector<uintptr_t> m_founder_info;
    std::vector<uintptr_t> m_founder_include2;
    std::vector<uintptr_t> m_sample_include;
    std::vector<uintptr_t> m_sample_include2;
    std::vector<uintptr_t> m_in_regression;
    std::vector<uintptr_t> m_haploid_mask;
    std::vector<size_t> m_sort_by_p_index;
    // std::vector<uintptr_t> m_sex_male;
    std::vector<int32_t> m_xymt_codes;
    // std::vector<int32_t> m_chrom_start;
    // sample file name. Fam for plink
    std::string m_sample_file;
    double m_mean_score = 0.0;
    double m_score_sd = 0.0;
    double m_hard_threshold = 0.0;
    double m_dose_threshold = 0.0;
    double m_clump_r2 = 0.0;
    double m_clump_proxy = 0.0;
    double m_clump_p = 0.0;
    double m_homcom_weight = 0;
    double m_het_weight = 1;
    double m_homrar_weight = 2;
    size_t m_num_thresholds = 0;
    size_t m_thread = 1; // number of final samples
    size_t m_max_window_size = 0;
    size_t m_num_ambig = 0;
    size_t m_num_maf_filter = 0;
    size_t m_num_geno_filter = 0;
    size_t m_num_info_filter = 0;
    size_t m_num_xrange = 0;
    size_t m_base_missed = 0;
    uintptr_t m_unfiltered_sample_ct = 0; // number of unfiltered samples
    uintptr_t m_unfiltered_marker_ct = 0;
    uintptr_t m_clump_distance = 0;
    uintptr_t m_sample_ct = 0;
    uintptr_t m_founder_ct = 0;
    uintptr_t m_marker_ct = 0;
    uint32_t m_max_category = 0;
    uint32_t m_autosome_ct = 0;
    uint32_t m_max_code = 0;
    std::random_device::result_type m_seed = 0;
    uint32_t m_num_ref_target_mismatch = 0;
    uint32_t m_num_male = 0;
    uint32_t m_num_female = 0;
    uint32_t m_num_ambig_sex = 0;
    uint32_t m_num_non_founder = 0;
    bool m_use_proxy = false;
    bool m_ignore_fid = false;
    bool m_is_ref = false;
    bool m_keep_nonfounder = false;
    bool m_keep_ambig = false;
    bool m_remove_sample = true;
    bool m_exclude_snp = true;
    bool m_hard_coded = false;
    bool m_expect_reference = false;
    bool m_mismatch_file_output = false;
    Reporter* m_reporter;
    MODEL m_model = MODEL::ADDITIVE;
    MISSING_SCORE m_missing_score = MISSING_SCORE::MEAN_IMPUTE;
    SCORING m_scoring = SCORING::AVERAGE;

    /*!
     * \brief Calculate the threshold bin based on the p-value and bound
     * info \param pvalue the input p-value \param bound_start is the start
     * of p-value threshold \param bound_inter is the step size of p-value
     * threshold \param bound_end is the end of p-value threshold \param
     * pthres return the name of p-value threshold this SNP belongs to
     * \param no_full indicate if we want the p=1 threshold
     * \return the category where this SNP belongs to
     */
    int calculate_category(const double& pvalue, const double& bound_start,
                           const double& bound_inter, const double& bound_end,
                           double& pthres, const bool no_full)
    {
        // NOTE: Threshold is x < p <= end and minimum category is 0
        int category = 0;
        if (pvalue > bound_end && !no_full)
        {
            category = static_cast<int>(
                std::ceil((bound_end + 0.1 - bound_start) / bound_inter));
            pthres = 1.0;
        }
        else
        {
            category = static_cast<int>(
                std::ceil((pvalue - bound_start) / bound_inter));
            category = (category < 0) ? 0 : category;
            pthres = category * bound_inter + bound_start;
        }
        return category;
    }

    /*!
     * \brief Calculate the category based on the input pvalue and barlevels
     * \param pvalue the input pvalue
     * \param barlevels the bar levels
     * \param pthres the p-value threshold this SNP belong to
     * \return the category of this SNP
     */
    int calculate_category(const double& pvalue,
                           const std::vector<double>& barlevels, double& pthres)
    {
        for (size_t i = 0; i < barlevels.size(); ++i)
        {
            if (pvalue < barlevels[i]
                || misc::logically_equal(pvalue, barlevels[i]))
            {
                pthres = barlevels[i];
                return static_cast<int>(i);
            }
        }
        pthres = 1.0;
        return static_cast<int>(barlevels.size());
    }

    /*!
     * \brief Replace # in the name of the genotype file and generate list
     * of file for subsequent analysis \param prefix contains the name of
     * the genotype file \return a vector of string containing names of the
     * genotype files
     */
    std::vector<std::string> set_genotype_files(const std::string& prefix);
    /*!
     * \brief Read in the genotype list file and add the genotype file names to
     *        the vector
     * \param file_name is the name of the list file
     * \return a vector of string containing names of all genotype files
     */
    std::vector<std::string> load_genotype_prefix(const std::string& file_name);
    /*!
     * \brief Initialize vector related to chromosome information e.g.haplotype.
     *        Currently not really useful except for setting the max_code which
     *        is later used to transform chromosome strings to chromosome code
     * \param num_auto is the number of autosome, we fix it to 22 for human
     * \param no_x indicate if chrX is missing for this organism
     * \param no_y indicate if chrY is missing for this organism
     * \param no_xy indicate if chrXY is missing for this organism
     * \param no_mt indicate if chrMT is missing for this organism
     */
    void init_chr(int num_auto = 22, bool no_x = false, bool no_y = false,
                  bool no_xy = false, bool no_mt = false);
    /*!
     * \brief For a given error code, check if this chromsome should be kept
     * \param chr_code is the chromosome code
     * \param sex_error indicate if the error is related to sex chromosome
     * \param chr_error indicate if the error is related to chromosome number
     * too big
     * \param error_message contain the error message
     * \return true if we want to skip this chromosome
     */
    bool chr_code_check(int32_t chr_code, bool& sex_error, bool& chr_error,
                        std::string& error_message);
    /*!
     * \brief Function to read in the sample. Any subclass must implement this
     * function. They \b must initialize the \b m_sample_info \b m_founder_info
     * \b m_founder_ct \b m_sample_ct \b m_prs_info \b m_in_regression and \b
     * m_tmp_genotype (optional)
     * \return vector containing the sample information
     */
    virtual std::vector<Sample_ID> gen_sample_vector(const std::string& delim)
    {
        return std::vector<Sample_ID>(0);
    }

    virtual void gen_snp_vector(
        const std::vector<IITree<size_t, size_t>>& /*exclusion_regions*/,
        const std::string& /*out_prefix*/, Genotype* /*target*/)
    {
    }
    virtual void calc_freq_gen_inter(
        const double& /*maf_threshold*/, const double& /*geno_threshold*/,
        const double& /*info_threshold*/, const bool /*maf_filter*/,
        const bool /*geno_filter*/, const bool /*info_filter*/,
        const bool /*hard_coded*/, Genotype* /*target=nullptr*/)
    {
    }
    void read_prs(std::vector<uintptr_t>& genotype, const size_t ploidy,
                  const double stat, const double adj_score,
                  const double miss_score, const size_t miss_count,
                  const double homcom_weight, const double het_weight,
                  const double homrar_weight, const bool not_first)
    {
        uintptr_t* lbptr = genotype.data();
        uintptr_t ulii;
        uint32_t uii;
        uint32_t ujj;
        uint32_t ukk;
        uii = 0;
        ulii = 0;
        do
        {
            // ulii contain the numeric representation of the current genotype
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct)
            {
                // this is PLINK, not sure exactly what this is about
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            // ujj sample index of the current genotype block
            ujj = 0;
            while (ujj < BITCT)
            {
                // go through the whole genotype block
                // ukk is the current genotype
                ukk = (ulii >> ujj) & 3;
                // and the sample index can be calculated as uii+(ujj/2)
                if (uii + (ujj / 2) >= m_sample_ct) { break; }
                auto&& sample_prs = m_prs_info[uii + (ujj / 2)];
                // now we will get all genotypes (0, 1, 2, 3)
                if (not_first)
                {
                    switch (ukk)
                    {
                    default:
                        // true = 1, false = 0
                        sample_prs.num_snp += ploidy;
                        sample_prs.prs += homcom_weight * stat - adj_score;
                        break;
                    case 1:
                        sample_prs.num_snp += ploidy;
                        sample_prs.prs += het_weight * stat - adj_score;
                        break;
                    case 3:
                        sample_prs.num_snp += ploidy;
                        sample_prs.prs += homrar_weight * stat - adj_score;
                        break;
                    case 2:
                        // handle missing sample
                        sample_prs.num_snp += miss_count;
                        sample_prs.prs += miss_score;
                        break;
                    }
                }
                else
                {
                    switch (ukk)
                    {
                    default:
                        // true = 1, false = 0
                        sample_prs.num_snp = ploidy;
                        sample_prs.prs = homcom_weight * stat - adj_score;
                        break;
                    case 1:
                        sample_prs.num_snp = ploidy;
                        sample_prs.prs = het_weight * stat - adj_score;
                        break;
                    case 3:
                        sample_prs.num_snp = ploidy;
                        sample_prs.prs = homrar_weight * stat - adj_score;
                        break;
                    case 2:
                        // handle missing sample
                        sample_prs.num_snp = miss_count;
                        sample_prs.prs = miss_score;
                        break;
                    }
                }
                // ulii &= ~((3 * ONELU) << ujj);
                // as each sample is represented by two byte, we will add 2 to
                // the index
                ujj += 2;
            }
            // uii is the number of samples we have finished so far
            uii += BITCT2;
        } while (uii < m_sample_ct);
    }

    /*!
     * \brief Function to read in the genotype in PLINK binary format. Any
     * subclass must implement this function to assist the processing of their
     * specific file type. The first argument is the genotype vector use to
     * store the PLINK binary whereas the second parameter is the streampos,
     * allowing us to seekg to the specific location of the file as indicate in
     * the thrid parameter
     */
    virtual inline void read_genotype(uintptr_t* /*genotype*/,
                                      const unsigned long long /*byte_pos*/,
                                      const size_t& /*file_index*/)
    {
    }
    virtual void
    read_score(const std::vector<size_t>::const_iterator& /*start*/,
               const std::vector<size_t>::const_iterator& /*end*/,
               bool /*reset_zero*/, const bool /*use_ref_maf*/)
    {
    }

    // for loading the sample inclusion / exclusion set
    /*!
     * \brief Function to load in the sample extraction exclusion list
     * \param input the file name
     * \param ignore_fid whether we should ignore the FID (use 2 column or 1)
     * \return an unordered_set use for checking if the sample is in the file
     */
    std::unordered_set<std::string> load_ref(const std::string& input,
                                             const std::string& delim,
                                             bool ignore_fid);
    /*!
     * \brief Function to load in SNP extraction exclusion list
     * \param input the file name of the SNP list
     * \param reporter the logger
     * \returnan unordered_set use for checking if the SNP is in the file
     */
    std::unordered_set<std::string> load_snp_list(const std::string& input);
    void shrink_snp_vector(const std::vector<bool>& retain)
    {
        m_existed_snps.erase(
            std::remove_if(m_existed_snps.begin(), m_existed_snps.end(),
                           [&retain, this](const SNP& s) {
                               return !retain[(&s - &*begin(m_existed_snps))];
                           }),
            m_existed_snps.end());
        m_existed_snps.shrink_to_fit();
    }
    /** Misc information **/
    // uint32_t m_hh_exists;
    /*!
     * \brief Function to perform the pearson correlation calculation
     * \param reference is the reference genotype
     * \param reporter is the logger
     */
    void pearson_clump(Genotype& reference, Reporter& reporter);

    /*!
     * \brief Function to check if the two alleles are ambiguous
     * \param ref_allele the reference allele
     * \param alt_allele the alternative allele
     * \return true if it is ambiguous
     */
    inline bool ambiguous(const std::string& ref_allele,
                          const std::string& alt_allele) const
    {
        // the allele should all be in upper case but whatever
        // true if equal
        if (ref_allele.empty() || alt_allele.empty()) return false;
        if (ref_allele == alt_allele) return true;
        return (ref_allele == "A" && alt_allele == "T")
               || (alt_allele == "A" && ref_allele == "T")
               || (ref_allele == "G" && alt_allele == "C")
               || (alt_allele == "G" && ref_allele == "C");
    }


    // modified version of the
    // single_marker_freqs_and_hwe function from PLINK (plink_filter.c)
    // we remove the HWE calculation (we don't want to include that yet,
    // as that'd require us to implement a version for bgen)
    inline void single_marker_freqs_and_hwe(
        uintptr_t unfiltered_sample_ctl2, uintptr_t* lptr,
        uintptr_t* sample_include2, uintptr_t* founder_include2,
        uintptr_t sample_ct, uint32_t* ll_ctp, uint32_t* lh_ctp,
        uint32_t* hh_ctp, uintptr_t sample_f_ct, uint32_t* ll_ctfp,
        uint32_t* lh_ctfp, uint32_t* hh_ctfp)
    {
        uint32_t tot_a = 0;
        uint32_t tot_b = 0;
        uint32_t tot_c = 0;
        uint32_t tot_a_f = 0;
        uint32_t tot_b_f = 0;
        uint32_t tot_c_f = 0;
        uintptr_t* lptr_end = &(lptr[unfiltered_sample_ctl2]);
        uintptr_t loader;
        uintptr_t loader2;
        uintptr_t loader3;
#ifdef __LP64__
        uintptr_t cur_decr = 120;
        uintptr_t* lptr_12x_end;
        unfiltered_sample_ctl2 -= unfiltered_sample_ctl2 % 12;
        while (unfiltered_sample_ctl2 >= 120)
        {
        single_marker_freqs_and_hwe_loop:
            lptr_12x_end = &(lptr[cur_decr]);
            count_3freq_1920b((__m128i*) lptr, (__m128i*) lptr_12x_end,
                              (__m128i*) sample_include2, &tot_a, &tot_b,
                              &tot_c);
            count_3freq_1920b((__m128i*) lptr, (__m128i*) lptr_12x_end,
                              (__m128i*) founder_include2, &tot_a_f, &tot_b_f,
                              &tot_c_f);
            lptr = lptr_12x_end;
            sample_include2 = &(sample_include2[cur_decr]);
            founder_include2 = &(founder_include2[cur_decr]);
            unfiltered_sample_ctl2 -= cur_decr;
        }
        if (unfiltered_sample_ctl2)
        {
            cur_decr = unfiltered_sample_ctl2;
            goto single_marker_freqs_and_hwe_loop;
        }
#else
        uintptr_t* lptr_twelve_end =
            &(lptr[unfiltered_sample_ctl2 - unfiltered_sample_ctl2 % 12]);
        while (lptr < lptr_twelve_end)
        {
            count_3freq_48b(lptr, sample_include2, &tot_a, &tot_b, &tot_c);
            count_3freq_48b(lptr, founder_include2, &tot_a_f, &tot_b_f,
                            &tot_c_f);
            lptr = &(lptr[12]);
            sample_include2 = &(sample_include2[12]);
            founder_include2 = &(founder_include2[12]);
        }
#endif
        while (lptr < lptr_end)
        {
            loader = *lptr++;
            loader2 = *sample_include2++;
            loader3 = (loader >> 1) & loader2;
            loader2 &= loader;
            // N.B. because of the construction of sample_include2, only
            // even-numbered bits can be present here.  So popcount2_long is
            // safe.
            tot_a += popcount2_long(loader2);
            tot_b += popcount2_long(loader3);
            tot_c += popcount2_long(loader & loader3);
            loader2 = *founder_include2++;
            loader3 = (loader >> 1) & loader2;
            loader2 &= loader;
            tot_a_f += popcount2_long(loader2);
            tot_b_f += popcount2_long(loader3);
            tot_c_f += popcount2_long(loader & loader3);
        }
        *hh_ctp = tot_c;
        *lh_ctp = tot_b - tot_c;
        *ll_ctp = sample_ct - tot_a - *lh_ctp;
        *hh_ctfp = tot_c_f;
        *lh_ctfp = tot_b_f - tot_c_f;
        *ll_ctfp = sample_f_ct - tot_a_f - *lh_ctfp;
    }
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
};

#endif /* SRC_GENOTYPE_HPP_ */
