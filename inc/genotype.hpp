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

#include "cgranges.h"
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
                      const std::string& remove_file, bool verbose,
                      Reporter& reporter);
    /*!
     * \brief Function to load SNPs into the genotype class object. Will call
     *        gen_snp_vector function. Will only generate the vector if this is
     *        the target
     * * \param maf_threshold is the maf threshold
     * \param maf_filter is the boolean indicate if we want to perform maf
     * filtering
     * \param geno_threshold is the geno threshold
     * \param geno_filter is the boolean indicate if we want to perform geno
     * filtering
     * \param hard_threshold is the hard coding threshold
     * \param hard_coded is the boolean indicate if hard coding should be
     * performed
     * \param info_threshold is the INFO score threshold
     * \param info_filter is the boolean indicate if we want to perform INFO
     * score filtering
     * \param exclusion is the region of exclusion
     * \param verbose is the boolean indicate if we want to print the message
     * \param reporter is the logger
     * \param target is the target genotype. Equal to nullptr if this is the
     * target
     */
    void load_snps(const std::string& out, const std::string& exclude,
                   const std::string& extract, const double& maf_threshold,
                   const double& geno_threshold, const double& info_threshold,
                   const double& hard_threshold, const bool maf_filter,
                   const bool geno_filter, const bool info_filter,
                   const bool hard_coded, cgranges_t* exclusion_region,
                   bool verbose, Reporter& reporter,
                   Genotype* target = nullptr);

    // do a quick filtering before we actually read in and process the genotypes
    void load_snps(const std::string& out, const std::string& exclude,
                   const std::string& extract, cgranges_t* exclusion_region,
                   bool verbose, Reporter& reporter,
                   Genotype* target = nullptr);

    void calc_freqs_and_intermediate(
        const double& maf_threshold, const double& geno_threshold,
        const double& info_threshold, const double& hard_threshold,
        const bool maf_filter, const bool geno_filter, const bool info_filter,
        const bool hard_coded, bool verbose, Reporter& reporter,
        Genotype* target = nullptr);
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
        for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp) {
            m_existed_snps_index[m_existed_snps[i_snp].rs()] = i_snp;
        }
    }

    size_t max_category() const { return m_max_category; }
    /*!
     * \brief Return the number of sample we wish to perform PRS on
     * \return the number of sample
     */
    size_t num_sample() const { return m_sample_id.size(); }
    /*!
     * \brief Function to obtain PRS score from the genotype file. Will assign
     * the result information to our PRS vector
     *
     * \param cur_index is the index of the starting SNP that we want to read
     * from
     *
     * \param cur_threshold return the p-value threshold of the current
     * processing chunck. This information is use for generating the output
     *
     * \param num_snp_included return the total number of SNP included in this
     * chunck
     *
     * \param region_index is the index of the region of interest
     * \param non_cumulate indicate if we want to perform cumulated PRS
     * calculation or not
     *
     * \param require_statistic if we want to standardize the PRS
     * \param first_run if we want to add or reset the PRS input
     * \return true if we can run, false otherwise
     */
    bool get_score(int& cur_index, double& cur_threshold,
                   uint32_t& num_snp_included, const size_t region_index,
                   const bool non_cumulate, const bool require_statistic,
                   const bool first_run);
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
    void efficient_clumping(Genotype& reference, Reporter& reporter,
                            bool const use_pearson);

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
     * \brief This function will provide the snp coordinate of a single snp,
     * return true when the SNP is found and false when it is not
     * \param rs_id is the rs id of the SNP
     * \param chr is chr return value
     * \param loc is the bp of the SNP
     * \return true if found, false if not
     */
    bool get_snp_loc(const std::string& rs_id, intptr_t& chr,
                     intptr_t& loc) const
    {
        auto&& snp_index = m_existed_snps_index.find(rs_id);
        if (snp_index == m_existed_snps_index.end()) return false;
        // TODO: Might want to unify the int type usage
        chr = m_existed_snps[snp_index->second].chr();
        loc = m_existed_snps[snp_index->second].loc();
        return true;
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
    std::string sample_id(size_t i) const
    {
        if (i > m_sample_id.size())
            throw std::out_of_range("Sample name vector out of range");
        if (m_ignore_fid)
            return m_sample_id[i].IID;
        else
            return m_sample_id[i].FID + "_" + m_sample_id[i].IID;
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
    double calculate_score(SCORING score_type, size_t i) const
    {
        if (i >= m_prs_info.size())
            throw std::out_of_range("Sample name vector out of range");
        double prs = m_prs_info[i].prs;
        int num_snp = m_prs_info[i].num_snp;
        double avg = prs;
        if (num_snp == 0) {
            avg = 0.0;
        }
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
    void get_null_score(const int& set_size, const int& prev_size,
                        const std::vector<size_t>& background_list,
                        const bool first_run, const bool require_standardize);
    /*!
     * \brief Return the number of SNPs included in the background
     * \return  the number of background SNPs
     */
    size_t num_background() const { return m_background_snp_index.size(); }
    /*!
     * \brief Get the index of background SNPs on the m_existed_snp vector
     * \return the index of background SNPs on the m_existed_snp vector
     */
    std::vector<size_t> background_index() const
    {
        return m_background_snp_index;
    }
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
    /*!
     * \brief Return the i th SNP, only use for unit testing
     * \param i is the index to the SNP
     * \return the ith SNP object
     */
    SNP get_snp(size_t i) const { return m_existed_snps.at(i); }
    void read_base(const std::string& base_file,
                   const std::vector<size_t>& col_index,
                   const std::vector<bool>& has_col,
                   const std::vector<double>& barlevels,
                   const double& bound_start, const double& bound_inter,
                   const double& bound_end, const double& maf_control,
                   const double& maf_case, const double& info_threshold,
                   const bool maf_control_filter, const bool maf_case_filter,
                   const bool info_filter, const bool fastscore,
                   const bool no_full, const bool is_beta, const bool is_index,
                   const bool keep_ambig, Reporter& reporter);
    void build_clump_windows();
    void build_membership_matrix(std::vector<size_t>& region_membership,
                                 std::vector<size_t>& region_start_idx,
                                 const size_t num_sets, const std::string& out,
                                 const std::vector<std::string>& region_name,
                                 const bool print_snps);
    size_t num_threshold() const { return m_num_thresholds; }
    std::vector<double> get_thresholds() const { return m_thresholds; }
    bool get_score(const std::vector<size_t>& region_membership,
                   size_t & start_index, const size_t &end_index,  double& cur_threshold,
                   uint32_t& num_snp_included,
                   const bool non_cumulate, const bool require_statistic,
                   const bool first_run);
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
    std::vector<std::string> m_genotype_files;
    std::vector<double> m_thresholds;
    std::vector<uintptr_t> m_tmp_genotype;
    // std::vector<uintptr_t> m_chrom_mask;
    std::vector<uintptr_t> m_founder_info;
    std::vector<uintptr_t> m_sample_include;
    std::vector<uintptr_t> m_in_regression;
    std::vector<uintptr_t> m_haploid_mask;
    std::vector<size_t> m_sort_by_p_index;
    std::vector<size_t> m_background_snp_index;
    // std::vector<uintptr_t> m_sex_male;
    std::vector<int32_t> m_xymt_codes;
    // std::vector<int32_t> m_chrom_start;
    // sample file name. Fam for plink
    std::string m_sample_file;
    std::mutex m_get_maf_mutex;
    std::mutex m_update_beta_mutex;
    double m_mean_score = 0.0;
    double m_score_sd = 0.0;
    double m_hard_threshold = 0.0;
    double m_clump_r2 = 0.0;
    double m_clump_proxy = 0.0;
    double m_clump_p = 0.0;
    double m_homcom_weight = 0;
    double m_het_weight = 1;
    double m_homrar_weight = 2;
    size_t m_num_thresholds = 0;
    uintptr_t m_unfiltered_sample_ct = 0; // number of unfiltered samples
    uintptr_t m_unfiltered_marker_ct = 0;
    uintptr_t m_clump_distance = 0;
    uintptr_t m_sample_ct = 0;
    uintptr_t m_founder_ct = 0;
    uintptr_t m_marker_ct = 0;
    intptr_t m_max_window_size = 0;
    uint32_t m_max_category = 0;
    uint32_t m_thread = 1; // number of final samples
    uint32_t m_autosome_ct = 0;
    uint32_t m_max_code = 0;
    std::random_device::result_type m_seed = 0;
    uint32_t m_num_ambig = 0;
    uint32_t m_num_ref_target_mismatch = 0;
    uint32_t m_num_maf_filter = 0;
    uint32_t m_num_geno_filter = 0;
    uint32_t m_num_info_filter = 0;
    uint32_t m_num_male = 0;
    uint32_t m_num_female = 0;
    uint32_t m_num_ambig_sex = 0;
    uint32_t m_num_non_founder = 0;
    uint32_t m_base_missed = 0;
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
        if (pvalue > bound_end && !no_full) {
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
        for (std::vector<double>::size_type i = 0; i < barlevels.size(); ++i) {
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
    virtual std::vector<Sample_ID> gen_sample_vector()
    {
        return std::vector<Sample_ID>(0);
    }
    /*!
     * \brief Function to read in the SNP information. Any subclass must
     * implement this function \return a vector containing the SNP information
     */
    virtual std::vector<SNP>
    gen_snp_vector(const std::string& /*out_prefix*/,
                   const double& /*maf_threshold*/, const bool /*maf_filter*/,
                   const double& /*geno_threshold*/, const bool /*geno_filter*/,
                   const double& /*hard_threshold*/, const bool /*hard_coded*/,
                   const double& /*info_threshold*/, const bool /*info_filter*/,
                   cgranges_t* /*exclusion*/, Genotype* /*target*/)
    {
        return std::vector<SNP>(0);
    }

    virtual void gen_snp_vector(const std::string& /*out_prefix*/,
                                cgranges_t* /*exclusion*/, Genotype* /*target*/)
    {
    }
    virtual void calc_freq_gen_inter(
        const double& /*maf_threshold*/, const double& /*geno_threshold*/,
        const double& /*info_threshold*/, const double& /*hard_threshold*/,
        const bool /*maf_filter*/, const bool /*geno_filter*/,
        const bool /*info_filter*/, const bool /*hard_coded*/,
        Genotype* /*target=nullptr*/)
    {
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
                                      const std::streampos /*byte_pos*/,
                                      const std::string& /*file_name*/)
    {
    }
    /*!
     * \brief read_score is the master function for performing the score
     * reading. All subclass must implement this function to assist the
     * calculation of PRS in the corresponding file format
     */
    virtual void read_score(const size_t /*start_index*/,
                            const size_t /*end_bound*/,
                            const size_t /*region_index*/, bool /*reset_zero*/)
    {
    }
    /*!
     * \brief similar to read_score function, but instead of using the index as
     * an input, we use a vector containing the required index as an input
     */
    virtual void read_score(const std::vector<size_t>& /*index*/,
                            bool /*reset_zero*/)
    {
    }


    // for loading the sample inclusion / exclusion set
    /*!
     * \brief Function to load in the sample extraction exclusion list
     * \param input the file name
     * \param ignore_fid whether we should ignore the FID (use 2 column or 1)
     * \return an unordered_set use for checking if the sample is in the file
     */
    std::unordered_set<std::string> load_ref(std::string input,
                                             bool ignore_fid);
    /*!
     * \brief Function to load in SNP extraction exclusion list
     * \param input the file name of the SNP list
     * \param reporter the logger
     * \returnan unordered_set use for checking if the SNP is in the file
     */
    std::unordered_set<std::string> load_snp_list(std::string input,
                                                  Reporter& reporter);

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
        while (unfiltered_sample_ctl2 >= 120) {
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
        if (unfiltered_sample_ctl2) {
            cur_decr = unfiltered_sample_ctl2;
            goto single_marker_freqs_and_hwe_loop;
        }
#else
        uintptr_t* lptr_twelve_end =
            &(lptr[unfiltered_sample_ctl2 - unfiltered_sample_ctl2 % 12]);
        while (lptr < lptr_twelve_end) {
            count_3freq_48b(lptr, sample_include2, &tot_a, &tot_b, &tot_c);
            count_3freq_48b(lptr, founder_include2, &tot_a_f, &tot_b_f,
                            &tot_c_f);
            lptr = &(lptr[12]);
            sample_include2 = &(sample_include2[12]);
            founder_include2 = &(founder_include2[12]);
        }
#endif
        while (lptr < lptr_end) {
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
