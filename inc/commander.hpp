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

#ifndef COMMANDER_H
#define COMMANDER_H

#include "gzstream.h"
#include "misc.hpp"
#include "storage.hpp"
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <random>
#include <reporter.hpp>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>
#ifdef _WIN32
#include <windows.h>
#endif

const std::string version = "2.2.7";
const std::string date = "2019-08-21";
class Commander
{
public:
    Commander();
    virtual ~Commander();
    bool init(int argc, char* argv[], Reporter& reporter);

    // base
    /*!
     * \brief Return the column index of the base file
     * \return The column index of the base
     */
    std::vector<size_t> index() const { return m_base_col_index; }
    std::vector<bool> has_col() const { return m_base_has_col; }
    /*!
     * \brief Get the base file name
     * \return Base file name
     */
    std::string base_name() const { return m_base_file; }

    /*!
     * \brief Get the base info score threshold
     * \return return the base info score filtering threshold, 0.0 if not used
     */
    double base_info_score() const { return m_base_info_threshold; }
    /*!
     * \brief Get the base maf threshold for controls
     * \return return the maf base control filtering threshold, 0.0 if not used
     */
    double maf_base_control() const { return m_control_maf_threshold; }

    /*!
     * \brief Get the base maf threshold for cases
     * \return return the maf base case filtering threshold, 0.0 if not used
     */
    double maf_base_case() const { return m_case_maf_threshold; }
    /*!
     * \brief Check if column input is in index form
     * \return true if input is in index
     */
    bool is_index() const { return m_input_is_index; }
    /*!
     * \brief Check if input statistic is beta
     * \return true if it is beta
     */
    bool beta() const { return m_stat_is_beta; }

    // clump
    /*!
     * \brief Get proxy clump threshold
     * \param storage of threshold
     * \return true if proxy clump is to be performed
     */
    bool proxy(double& threshold) const
    {
        threshold = m_proxy_threshold;
        return m_use_proxy_clump;
    }
    /*!
     * \brief Get p-value threshold for clumping
     * \return p-value threshold
     */
    double clump_p() const { return m_clump_p; }
    /*!
     * \brief Get r2 threshold for clumping
     * \return R2 threshold
     */
    double clump_r2() const { return m_clump_r2; }
    /*!
     * \brief Get clumping distance in \b bp
     * \return Clumping distance in bp
     */
    int clump_dist() const { return m_clump_distance; }
    /*!
     * \brief Check if user want to skip clumping
     * \return true if user want to skip clumping
     */
    bool no_clump() const { return m_no_clump; }
    /*!
     * \brief Check if user want to perform proxy clumping
     * \return true if user want to perform proxy clumping
     */
    bool use_proxy() const { return m_use_proxy_clump; }

    // covariate
    /*!
     * \brief Return the covariate file name
     * \return The covariate file name. An empty string if user did not
     *         provide a covariate file
     */
    std::string get_cov_file() const { return m_cov_file; }
    /*!
     * \brief Get the name of the covariates
     * \return The name of the covariates in form of vector
     */
    std::vector<std::string> get_cov_name() const { return m_cov_colname; }
    /*!
     * \brief Get the index of the covariates
     * \return The column index of the covariates
     */
    std::vector<size_t> get_cov_index() const { return m_col_index_of_cov; }
    /*!
     * \brief Get the index of the covariate that are factor
     * \return The column index of the factor covariates
     */
    std::vector<size_t> get_factor_cov_index() const
    {
        return m_col_index_of_factor_cov;
    }
    // reference panel
    /*!
     * \brief Function to return the reference file name
     * \param String to store the reference file name to
     * \return True if reference file is provided, false otherwise
     */
    bool ref_name(std::string& ref_name) const
    {
        ref_name = m_ref_file;
        return m_use_reference;
    }
    /*!
     * \brief Function to directly return the reference file name
     * \return empty if reference file name is not provided
     */
    std::string ref_name() const
    {
        if (m_use_reference) return m_ref_file;
        return "";
    }
    /*!
     * \brief Function to return the file name containing names to the reference
     *        files
     * \param String to store the file name
     * \return true if file is provided, false otherwise
     */
    bool ref_list(std::string& ref_list) const
    {
        ref_list = m_ref_list;
        return m_ref_list_provided;
    }
    /*!
     * \brief function to directly return the file name contain names to the
     * reference files
     * \return empty string if reference list file is not provided
     */
    std::string ref_list() const
    {
        if (m_ref_list_provided) return m_ref_list;
        return "";
    }
    /*!
     * \brief Return the type of the reference file
     * \return the type of the reference file
     */
    std::string ref_type() const { return m_ref_type; }
    /*!
     * \brief Return the file containing samples to keep in reference
     * \return The file name of the keep file
     */
    std::string ref_keep_file() const { return m_ref_keep; }
    /*!
     * \brief Return the file containing samples to remove in reference
     * \return The file name of the remove file
     */
    std::string ref_remove_file() const { return m_ref_remove; }
    /*!
     * \brief Check if reference panel is used
     * \return true if reference panel is used
     */
    bool use_ref() const { return m_use_reference; }
    /*!
     * \brief Check if intermediate file is allowed
     * \return true if intermeidate file is allowed
     */
    bool use_inter() const { return m_allow_inter; }
    // reference filtering
    /*!
     * \brief Get genotype filtering threshold for reference
     * \param geno is the storage of threshold
     * \return true if user want to perform genotype missingness filtering
     */
    bool ref_geno(double& geno) const
    {
        geno = m_ref_geno;
        return m_perform_ref_geno_filter;
    }
    bool get_ref_hard_threshold() const { return m_ref_hard_threshold; }
    bool target_hard_threshold() const { return m_target_hard_thresholding; }
    bool ref_hard_threshold() const { return m_perform_ref_hard_thresholding; }
    /*!
     * \brief Get the MAF filtering threshold for reference
     * \param maf is the storage of the threshold
     * \return true if maf filtering is required
     */
    bool ref_maf(double& maf) const
    {
        maf = m_ref_maf;
        return m_perform_ref_maf_filter;
    }
    /*!
     * \brief Get the info score threshold for reference
     * \param info is the storage of threshold
     * \return true if info score filtering is required
     */
    bool ref_info(double& info) const
    {
        info = m_ref_info_score;
        return m_perform_ref_info_filter;
    }
    // misc
    std::string delim() const { return m_id_delim; }
    /*!
     * \brief Get output prefix
     * \return the output prefix
     */
    std::string out() const { return m_out_prefix; }
    /*!
     * \brief Get the exclusion_range
     * \return The string containign the exclusion range
     */
    std::string exclusion_range() const { return m_exclusion_range; }
    /*!
     * \brief Return if all score should be printed
     * \return True if all score should be printed
     */
    bool all_scores() const { return m_print_all_scores; }
    /*!
     * \brief Return if non-cumulative PRS should be calculated
     * \return True for non-cumulative PRS calculation
     */
    bool non_cumulate() const { return m_non_cumulate_prs; }
    /*!
     * \brief Return if we should ignore the FID field for all input
     * \return True if we should
     */
    bool ignore_fid() const { return m_ignore_fid; }
    /*!
     * \brief Return if we should perform logistic regression for permutation
     * \return true if we should
     */
    bool logit_perm() const { return m_logit_perm; }
    /*!
     * \brief Return if we should generate the .snp file
     * \return true if we should
     */
    bool print_snp() const { return m_print_snp; }
    /*!
     * \brief Return if pearson correlation should be used instead
     * \return True if pearson correlation should be used
     */
    bool pearson() const { return m_pearson; }
    /*!
     * \brief Get the number of permutation to be performed
     * \param perm stores the number of permutation to be performed
     * \return true if we want to perform permutation
     */
    bool num_perm(size_t& perm) const
    {
        perm = m_permutation;
        return m_perform_permutation;
    }
    /*!
     * \brief Return the seed used for the run
     * \return the seed
     */
    std::random_device::result_type seed() const { return m_seed; }
    /*!
     * \brief Return the number of thread to be used
     * \return the number of thread
     */
    size_t thread() const { return static_cast<size_t>(m_thread); }
    /*!
     * \brief Return the maximum memory allowed to use
     * \param detected should be the available memory
     * \return Maximum memory we can use
     */
    unsigned long long max_memory(const unsigned long long detected) const
    {
        return (!m_provided_memory || m_memory > detected) ? detected
                                                           : m_memory;
    }
    unsigned long long memory() const { return m_memory; }
    /*!
     * \brief Return the full bar level
     * \return the vector containing the sorted bar levels
     */
    std::vector<double> bar_levels() const { return m_barlevel; }
    /*!
     * \brief Return the p-value threshold corresponding to the category.
     *        For fastscoring only
     * \param i is the category id
     * \return the p-value threshold corresponding to the category id
     */
    double get_threshold(int i) const
    {
        // Any category less than 0 is assume to be 0
        if (i < 0) return m_barlevel.front();
        // Any category larger than the bar-level size will be assumed to have p
        // threshold of 1
        // use intermediate to avoid implicite conversion
        auto&& index = static_cast<size_t>(i);
        if (index > m_barlevel.size()) return 1.0;
        return m_barlevel.at(index);
    }
    /*!
     * \brief Get the starting threhsold for high resolution scoring
     * \return The starting threshold
     */
    double lower() const { return m_lower_threshold; }
    /*!
     * \brief Get the last threhsold for high resolution scoring
     * \return The last threshold
     */
    double upper() const { return m_upper_threshold; }
    /*!
     * \brief Get the step size for high resolution scoring
     * \return The step size
     */
    double inter() const { return m_inter_threshold; }


    /*!
     * \brief Check if user do not want to include p-value threshold of 1
     * \return True if user do not want to include p-value threshold of 1
     */
    bool no_full() const { return m_no_full; }
    /*!
     * \brief Check if user wants to perform fastscoring
     * \return true if user wants to perform fastscoring
     */
    bool fastscore() const { return m_fastscore; }

    // prs_calculation
    /*!
     * \brief Get the missing score handling method
     * \return return the missing score enum
     */
    MISSING_SCORE get_missing_score() const { return m_missing_score; }
    /*!
     * \brief Get the PRS score calculation method
     * \return  return the scoring enum
     */
    SCORING get_score() const { return m_scoring_method; }
    /*!
     * \brief Get the genetic model use for calculation
     * \return return the model enum
     */
    MODEL model() const { return m_genetic_model; }
    /*!
     * \brief Indicate if we want to skip performing regression
     * \return true if regression is to be skipped
     */
    bool no_regress() const { return m_no_regress; }

    // prs_snp_filtering
    /*!
     * \brief Return the file name of the exclusion file
     * \return the name of the exclusion file
     */
    std::string exclude_file() const { return m_exclude_file; }
    /*!
     * \brief Return the file name of the extraction file
     * \return the name of the extract file
     */
    std::string extract_file() const { return m_extract_file; }
    /*!
     * \brief Get the genotype filtering threshold of target
     * \param The threshold is stored in geno
     * \return True if genotype filtering is required
     */
    bool target_geno(double& geno) const
    {
        geno = m_target_geno;
        return m_target_geno_filter;
    }
    double get_target_hard_threshold() const { return m_target_hard_threshold; }
    /*!
     * \brief Get the maf filtering threshold of target
     * \param The threshold is stored in maf
     * \return True if maf filtering is required
     */
    bool target_maf(double& maf) const
    {
        maf = m_target_maf;
        return m_target_maf_filter;
    }

    /*!
     * \brief Get the INFO filtering threshold of target
     * \param The threshold is stored in info
     * \return True if INFO filtering is required
     */
    bool target_info(double& info) const
    {
        info = m_target_info_score;
        return m_target_info_filter;
    }

    /*!
     * \brief check if target should be hard coded
     * \return true if target should be hard coded
     */
    bool hard_coded() const { return m_target_is_hard_coded; }
    /*!
     * \brief Check if we should retain ambiguous SNPs
     * \return true if ambiguous SNPs should be retained
     */
    bool keep_ambig() const { return m_keep_ambig; }

    // prset
    /*!
     * \brief Return the list of bed files
     * \return vector contain bed file names
     */
    std::vector<std::string> bed() const { return m_bed_files; }
    /*!
     * \brief Return the features to be included
     * \return vector containing the features
     */
    std::vector<std::string> feature() const { return m_feature; }
    /*!
     * \brief Return the gtf file name
     * \return GTF file name
     */
    std::string gtf() const { return m_gtf; }
    /*!
     * \brief Return the msigdb file name
     * \return MSigDB file name
     */
    std::vector<std::string> msigdb() const { return m_msigdb; }
    /*!
     * \brief snp_set return file name containing SNP set(s)
     * \return File name containing the SNP set(s)
     */
    std::vector<std::string> snp_set() const { return m_snp_set; }

    /*!
     * \brief Return file name containing background regions
     * \return File name containing backgroudn regions
     */
    std::string background() const { return m_background; }

    /*!
     * \brief Get the number of permutation to be performed for each set
     * \param perm stores the number of permutation
     * \return Return true if we want to perform set based permutation
     */
    bool set_perm(size_t& perm) const
    {
        perm = m_set_perm;
        return m_perform_set_perm;
    }
    /*!
     * \brief Check if we want to perform set based permutation
     * \return Return true if we want to perform set based permutation
     */
    bool perform_set_perm() const { return m_perform_set_perm; }
    /*!
     * \brief Check if we want the whole genome to be act as the background
     * \return True if we want the whole genome to be used as the background
     */
    bool genome_wide_background() const { return m_full_background; }
    size_t window_5() const { return m_window_5; }
    /*!
     * \brief Return the 3' extension of regions in \b bp
     * \return The length of extension to 3' regions
     */
    size_t window_3() const { return m_window_3; }
    // target
    /*!
     * \brief Return the target file name
     * \return the target file name. If target list is set, will return an empty
     * string
     */
    std::string target_name() const { return m_target_file; }
    /*!
     * \brief Return the file containing list of target files
     * \param target_list will be used to store the file name
     * \return true if target list is used
     */
    bool target_list(std::string& target_list) const
    {
        target_list = m_target_file_list;
        return m_use_target_list;
    }
    /*!
     * \brief Directly return the target file list
     * \return Return empty string if not provided
     */
    std::string target_list() const
    {
        if (m_use_target_list) return m_target_file_list;
        return "";
    }
    /*!
     * \brief Return the type of target
     * \return Type of target
     */
    std::string target_type() const { return m_target_type; }
    /*!
     * \brief Return the name of the phenotype file
     * \return name of phenotype file, empty if not provided
     */
    std::string pheno_file() const { return m_pheno_file; }
    /*!
     * \brief name of phenotype given the index
     * \param index is the column index of the phenotype file
     * \return The name of the phenotype. Error if out of bound
     */
    std::string pheno_col(size_t index) const { return m_pheno_col.at(index); }
    /*!
     * \brief File containing sample ID to be kept in target
     * \return name of the keep file
     */
    std::string keep_sample_file() const { return m_target_keep; }
    /*!
     * \brief File containing sample ID to be removed in target
     * \return name of the remove file
     */
    std::string remove_sample_file() const { return m_target_remove; }
    /*!
     * \brief Return the vector containing the phenotype names
     * \return Vector containing the phenotype names
     */
    std::vector<std::string> pheno_col() const { return m_pheno_col; }
    /*!
     * \brief Return vector containing if target phenotype is binary or not
     * \return Vector containing if teh target phenotype is binary or not
     */
    std::vector<bool> is_binary() const { return m_is_binary; }
    /*!
     * \brief Return vector containing the prevalence of target binary
     * phenotypes
     * \return vector containing the prevalence
     */
    std::vector<double> prevalence() const { return m_prevalence; }
    /*!
     * \brief return whether prevalence is provided
     * \return true if prevalence is provided
     */
    bool has_prevalence() const { return !m_prevalence.empty(); }
    /*!
     * \brief Check if we have phenotype names
     * \return true if we have phenotype names
     */
    bool has_pheno_col() const { return !m_pheno_col.empty(); }
    /*!
     * \brief Check if the ith phenotype is binary
     * \param index of phenotype
     * \return true if phenotype should be binary
     */
    bool is_binary(size_t index) const { return m_is_binary.at(index); }
    /*!
     * \brief Check if we want to include non-founders in our analysis
     * \return True if non-founders should be included
     */
    bool nonfounders() const { return m_include_nonfounders; }
    bool disable_mmap() const { return m_disable_mmap; }
    bool use_ref_maf() const { return m_use_ref_maf; }
    double ref_dose_thres() const { return m_ref_dose_thres; }
    double target_dose_thres() const { return m_target_dose_thres; }

protected:
private:
    const std::vector<std::string> supported_types = {"bed", "ped", "bgen"};
    std::vector<std::string> m_bed_files;
    std::vector<std::string> m_cov_colname;
    std::vector<std::string> m_factor_cov;
    std::vector<std::string> m_msigdb;
    std::vector<std::string> m_feature;
    std::vector<std::string> m_pheno_col;
    std::vector<std::string> m_snp_set;
    std::vector<double> m_barlevel;
    std::vector<double> m_prevalence;
    std::vector<size_t> m_col_index_of_cov;
    std::vector<size_t> m_col_index_of_factor_cov;
    std::vector<size_t> m_base_col_index;
    std::vector<bool> m_base_has_col;
    std::vector<bool> m_is_binary;
    std::string m_target_file = "";
    std::string m_target_file_list = "";
    std::string m_target_keep = "";
    std::string m_target_remove = "";
    std::string m_target_type = "bed";
    std::string m_pheno_file = "";
    std::string m_gtf = "";
    std::string m_background = "";
    std::string m_cov_file;
    std::string m_base_file = "";
    std::string m_chr = "CHR";
    std::string m_effect_allele = "A1";
    std::string m_non_effect_allele = "A2";
    std::string m_statistic = "";
    std::string m_snp = "SNP";
    std::string m_bp = "BP";
    std::string m_p_value = "P";
    std::string m_id_delim = " ";
    std::string m_info_col = "INFO,0.9";
    std::string m_maf_col;
    std::string m_out_prefix = "PRSice";
    std::string m_exclusion_range = "";
    std::string m_ref_file = "";
    std::string m_ref_list = "";
    std::string m_ref_type = "bed";
    std::string m_ref_keep = "";
    std::string m_ref_remove = "";
    std::string m_exclude_file = "";
    std::string m_extract_file = "";
    std::string m_help_message;
    double m_base_info_threshold = 0.0;
    // control maf threhsold is also used for quantitative trait
    double m_control_maf_threshold = 0.0;
    double m_case_maf_threshold = 0.0;
    double m_proxy_threshold = -1.0;
    double m_clump_p = 1.0;
    double m_clump_r2 = 0.1;
    double m_ref_dose_thres = 0.0;
    double m_ref_geno = 1.0;
    double m_ref_hard_threshold = 0.1;
    double m_ref_maf = 0.0;
    double m_ref_info_score = 0.0;
    // TODO: might consider using 1e-8 instead
    double m_lower_threshold = 5e-8;
    double m_inter_threshold = 0.00005;
    double m_upper_threshold = 0.5;
    double m_target_dose_thres = 0.0;
    double m_target_geno = 1.0;
    double m_target_hard_threshold = 0.1;
    double m_target_maf = 0.0;
    double m_target_info_score = 0.0;
    unsigned long long m_memory = 1e10;
    size_t m_clump_distance = 250000;
    size_t m_permutation = 0;
    size_t m_set_perm = 0;
    size_t m_window_5 = 0;
    size_t m_window_3 = 0;
    std::random_device::result_type m_seed = std::random_device()();
    MISSING_SCORE m_missing_score = MISSING_SCORE::MEAN_IMPUTE;
    SCORING m_scoring_method = SCORING::AVERAGE;
    MODEL m_genetic_model = MODEL::ADDITIVE;

    int m_thread = 1;
    int m_allow_inter = false;
    int m_disable_mmap = false;
    int m_fastscore = false;
    int m_full_background = false;
    int m_ignore_fid = false;
    int m_include_nonfounders = false;
    int m_input_is_index = false;
    int m_keep_ambig = false;
    int m_logit_perm = false;
    int m_no_clump = false;
    int m_no_full = false;
    int m_no_regress = false;
    int m_non_cumulate_prs = false;
    int m_pearson = false;
    int m_print_all_scores = false;
    int m_print_snp = false;
    int m_stat_is_beta = false;
    int m_stat_is_or = false;
    int m_target_is_hard_coded = false;
    int m_user_no_default = false;
    int m_use_ref_maf = false;
    bool m_use_reference = false;
    bool m_ref_list_provided = false;
    bool m_provided_seed = false;
    bool m_provided_memory = false;
    bool m_perform_permutation = false;
    bool m_provided_chr_col = false;
    bool m_provided_clump_dist = false;
    bool m_provided_effect_allele = false;
    bool m_provided_non_effect_allele = false;
    bool m_provided_statistic = false;
    bool m_provided_snp_id = false;
    bool m_provided_bp = false;
    bool m_provided_p_value = false;
    bool m_provided_info_threshold = false;
    bool m_perform_base_maf_filter = false;
    bool m_use_proxy_clump = false;
    bool m_perform_ref_geno_filter = false;
    bool m_perform_ref_hard_thresholding = false;
    bool m_perform_ref_maf_filter = false;
    bool m_perform_ref_info_filter = false;
    bool m_set_use_thresholds = false;
    bool m_set_delim = false;
    bool m_target_geno_filter = false;
    bool m_target_hard_thresholding = false;
    bool m_target_maf_filter = false;
    bool m_target_info_filter = false;
    bool m_perform_prset = false;
    bool m_perform_set_perm = false;
    bool m_use_target_list = false;

    ////////////////////////////////////////////
    //
    // end of struct definition
    //
    ////////////////////////////////////////////

    /*!
     * \brief Parsing all command line arguments
     * \param argc
     * \param argv
     * \param optString is the string containing the short flags
     * \param longOpts is the struct containing all the options
     * \param reporter is the object to report all messages
     * \return true if we want to continue the program
     */
    bool parse_command(int argc, char* argv[], const char* optString,
                       const struct option longOpts[], Reporter& reporter);
    /*!
     * \brief Print the usage information
     */
    void usage();
    /*!
     * \brief Set the help message
     */
    void set_help_message();
    /*!
     * \brief Function to check if parameters related to base are correct
     * \param message the parameter storage
     * \param error_message the storage for error messages
     * \return true if parameters are alright
     */
    bool base_check(std::map<std::string, std::string>& message,
                    std::string& error_message);
    /*!
     * \brief Function to check if parameters for clumping
     * \param message is the parameter storage
     * \param
     * error_message is the storage for error messages
     * \return  return true if
     * clumping parameters are alright
     */
    bool clump_check(std::map<std::string, std::string>& message,
                     std::string& error_message);
    /*!
     * \brief Function to check if parameters for reference panel are correct
     * \param message is the parameter storage
     * \param
     * error_message is the storage for error messages
     * \return  return true if
     * reference panel parameters are alright
     */
    bool ref_check(std::map<std::string, std::string>& message,
                   std::string& error_message);
    /*!
     * \brief Function to check if parameters for covariate are correct
     * \param error_message is the storage for error messages
     * \return true if parameters are alright
     */
    bool covariate_check(std::string& error_message);
    /*!
     * \brief Function to check if parameters for filtering are correct
     * \param error_message is the storage for error messages
     * \return true if parameters are alright
     */
    bool filter_check(std::string& error_message);
    /*!
     * \brief Function to check if parameters for misc options are correct
     * \param error_message is the storage for error messages
     * \return true if parameters are alright
     */
    bool misc_check(std::map<std::string, std::string>& message,
                    std::string& error_message);
    /*!
     * \brief Function to check if parameters for prset are correct
     * \param error_message is the storage for error messages
     * \return true if parameters are alright
     */
    bool prset_check(std::map<std::string, std::string>& message,
                     std::string& error_message);
    /*!
     * \brief Function to check if parameters for prsice thresholding are
     * correct \param error_message is the storage for error messages \return
     * true if parameters are alright
     */
    bool prsice_check(std::map<std::string, std::string>& message,
                      std::string& error_message);
    /*!
     * \brief Function to check if parameters for target are correct
     * \param error_message is the storage for error messages
     * \return true if parameters are alright
     */
    bool target_check(std::map<std::string, std::string>& message,
                      std::string& error_message);

    /*!
     * \brief Function to parse the binary vector input. Support short form
     * (e.g. 4T)
     * \param input is the input string
     * \param message is the parameter storage
     * \param error_message is the error message storage
     * \param target is the target variable
     * \param c is the flag
     * \return true if conversion is successful
     */
    inline bool parse_binary_vector(const std::string& input,
                                    std::map<std::string, std::string>& message,
                                    std::string& error_message,
                                    std::vector<bool>& target,
                                    const std::string& c)
    {
        // allow more complex regrex
        // should not come in without something
        if (input.empty()) return false;
        if (message.find(c) == message.end()) { message[c] = input; }
        else
        {
            message[c] = "," + input;
        }
        if (!input.empty() && input.back() == ',')
        {
            error_message.append("Warning: , detected at end of input: " + input
                                 + ". Have you accidentally included space in "
                                   "your input? (Space is not allowed)");
        }
        std::vector<std::string> token = misc::split(input, ",");
        try
        {
            for (auto&& bin : token)
            {
                // check if this is true or false, if, not, try parsing
                std::transform(bin.begin(), bin.end(), bin.begin(), ::toupper);
                if (bin == "T" || bin == "TRUE") { target.push_back(true); }
                else if (bin == "F" || bin == "FALSE")
                {
                    target.push_back(false);
                }
                // we likely have numeric input
                else
                {
                    std::string value_str = bin.substr(bin.length() - 1);
                    try
                    {
                        size_t repeat =
                            misc::string_to_size_t(value_str.c_str());
                        bool value;
                        if (bin.back() == 'T') { value = true; }
                        else if (bin.back() == 'F')
                        {
                            value = false;
                        }
                        else
                        {
                            throw std::runtime_error("");
                        }
                        for (size_t i = 0; i < repeat; ++i)
                        { target.push_back(value); }
                    }
                    catch (const std::runtime_error&)
                    {
                        throw std::runtime_error("");
                    }
                }
            }
        }
        catch (...)
        {
            error_message.append(
                "Error: Invalid argument passed to " + c + ": " + input
                + "! Require binary arguments e.g. T/F, True/False\n");
            return false;
        }
        return true;
    }
    /*!
     * \brief Load the string input into a vector
     * \param input the input
     * \param message the parameter storage
     * \param target where we store the output
     * \param c the flag of command
     * \param error_message the storage of the error messages
     */
    inline void load_string_vector(const std::string& input,
                                   std::map<std::string, std::string>& message,
                                   std::vector<std::string>& target,
                                   const std::string& c,
                                   std::string& error_message)
    {
        if (input.empty()) return;
        if (message.find(c) == message.end()) { message[c] = input; }
        else
        {
            message[c] = "," + input;
        }
        if (!input.empty() && input.back() == ',')
        {
            error_message.append("Warning: , detected at end of input: " + input
                                 + ". Have you accidentally included space in "
                                   "your input? (Space is not allowed)\n");
        }
        std::vector<std::string> token = misc::split(input, ",");
        target.insert(target.end(), token.begin(), token.end());
    }

    /*!
     * \brief Convert input into a vector of numeric values
     * \param input is the input value
     * \param message is the parameter storage
     * \param error_message is the error message storage
     * \param target is where we store the result
     * \param c the command flag
     * \return true if this is sucess
     */
    template <typename T>
    inline bool load_numeric_vector(const std::string& input,
                                    std::map<std::string, std::string>& message,
                                    std::string& error_message,
                                    std::vector<T>& target,
                                    const std::string& c)
    {
        // should always have an input
        if (input.empty()) return false;
        if (message.find(c) == message.end()) { message[c] = input; }
        else
        {
            // we will append instead of overwrite
            message[c] = "," + input;
        }
        if (!input.empty() && input.back() == ',')
        {
            error_message.append("Warning: , detected at end of input: " + input
                                 + ". Have you accidentally included space in "
                                   "your input? (Space is not allowed)\n");
        }
        std::vector<std::string> token = misc::split(input, ",");
        try
        {
            for (auto&& bar : token) target.push_back(misc::convert<T>(bar));
        }
        catch (...)
        {
            error_message.append("Error: Non numeric argument passed to " + c
                                 + ": " + input + "!\n");
            return false;
        }
        return true;
    }
    /*!
     * \brief Convert input into numeric value
     * \param input is the input value
     * \param message is the current output message related to command parse
     * \param error_message is the current error emssage store
     * \param target is where we are going to store the value
     * \param target_boolean is a boolean to indicate if the value is set
     * \param c is the command flag
     * \return true if we have sucessfully set the value
     */
    template <typename Type>
    inline bool set_numeric(const std::string& input,
                            std::map<std::string, std::string>& message,
                            std::string& error_message, Type& target,
                            bool& target_boolean, const std::string& c)
    {
        if (message.find(c) != message.end())
        {
            error_message.append("Warning: Duplicated argument --" + c + "\n");
        }
        message[c] = input;
        try
        {
            target = misc::convert<Type>(input);
            target_boolean = true;
        }
        catch (...)
        {
            error_message.append("Error: Non numeric argument passed to " + c
                                 + ": " + input + "!\n");
            return false;
        }
        return true;
    }

    /*!
     * \brief Function parsing string into genetic model used
     * \param in the input
     * \param message the parameter storage
     * \param error_message the error message storage
     * \return true if successfully parse the input into genetic model enum
     */
    inline bool set_model(const std::string& in,
                          std::map<std::string, std::string>& message,
                          std::string& error_message)
    {
        std::string input = in;
        if (input.empty())
        {
            error_message.append("Error: Model cannot be empty!\n");
            return false;
        }
        std::transform(input.begin(), input.end(), input.begin(), ::toupper);
        if (input.at(0) == 'A')
        {
            input = "add";
            m_genetic_model = MODEL::ADDITIVE;
        }
        else if (input.at(0) == 'D')
        {
            input = "dom";
            m_genetic_model = MODEL::DOMINANT;
        }
        else if (input.at(0) == 'R')
        {
            input = "rec";
            m_genetic_model = MODEL::RECESSIVE;
        }
        else if (input.at(0) == 'H')
        {
            input = "het";
            m_genetic_model = MODEL::HETEROZYGOUS;
        }
        else
        {
            error_message.append("Error: Unrecognized model: " + input + "!\n");
            return false;
        }
        if (message.find("model") != message.end())
        { error_message.append("Warning: Duplicated argument --model\n"); }
        message["model"] = input;
        return true;
    }

    /*!
     * \brief Function parsing string into MISSING enum
     * \param in the input
     * \param message the parameter storage
     * \param error_message the error message storage
     * \return true if successfully parse the input into MISSING enum
     */
    inline bool set_score(const std::string& in,
                          std::map<std::string, std::string>& message,
                          std::string& error_message)
    {
        std::string input = in;
        if (input.empty())
        {
            error_message.append(
                "Error: Score method string cannot be empty!\n");
            return false;
        }
        std::transform(input.begin(), input.end(), input.begin(), ::toupper);
        if (input.at(0) == 'A')
        {
            input = "avg";
            m_scoring_method = SCORING::AVERAGE;
        }
        else if (input.at(0) == 'S')
        {
            if (input.length() < 2)
            {
                error_message.append("Error: Ambiguous scoring method: " + input
                                     + "!\n");
                return false;
            }
            if (input.at(1) == 'T')
            {
                input = "standard";
                m_scoring_method = SCORING::STANDARDIZE;
            }
            else
            {
                input = "sum";
                m_scoring_method = SCORING::SUM;
            }
        }
        else
        {
            error_message.append("Error: Unrecognized scoring method: " + input
                                 + "!\n");
            return false;
        }
        if (message.find("score") != message.end())
        { error_message.append("Warning: Duplicated argument --score\n"); }
        message["score"] = input;
        return true;
    }

    /*!
     * \brief Function parsing string into PRS calculation enum
     * \param in the input
     * \param message the parameter storage
     * \param error_message the error message storage
     * \return true if successfully parse the input into PRS scoring enum
     */
    inline bool set_missing(const std::string& in,
                            std::map<std::string, std::string>& message,
                            std::string& error_message)
    {
        std::string input = in;
        if (input.empty())
        {
            error_message.append(
                "Error: Missing handling method string cannot be empty!\n");
            return false;
        }
        std::transform(input.begin(), input.end(), input.begin(), ::toupper);
        if (input.at(0) == 'C')
        {
            input = "CENTER";
            m_missing_score = MISSING_SCORE::CENTER;
        }
        else if (input.at(0) == 'M')
        {
            input = "MEAN_IMPUTE";
            m_missing_score = MISSING_SCORE::MEAN_IMPUTE;
        }
        else if (input.at(0) == 'S')
        {
            input = "SET_ZERO";
            m_missing_score = MISSING_SCORE::SET_ZERO;
        }
        else
        {
            error_message.append("Error: Unrecognized Missing handling method: "
                                 + input + "!\n");
            return false;
        }
        if (message.find("missing") != message.end())
        { error_message.append("Warning: Duplicated argument --score\n"); }
        message["missing"] = input;
        return true;
    }
    std::vector<std::string> transform_covariate(const std::string& cov_in);
    /*!
     * \brief Simply store the string input into the target variable
     * \param input The input string
     * \param message the parameter storage
     * \param target the target variable
     * \param target_boolean boolean indicate if the target is set
     * \param c the command flag
     * \param error_message error message storage
     */
    inline void set_string(const std::string& input,
                           std::map<std::string, std::string>& message,
                           std::string& target, bool& target_boolean,
                           const std::string& c, std::string& error_message,
                           bool add_quote = false)
    {
        if (message.find(c) != message.end())
        {
            error_message.append("Warning: Duplicated argument --" + c + "\n");
        }
        if (add_quote) { message[c] = "\"" + input + "\""; }
        else
        {
            message[c] = input;
        }
        target = input;
        target_boolean = true;
    }

    // currently, we don't anticipate anything more than TB
    /*!
     * \brief Convert user input into valid memory format
     * \param input the input string
     * \param message the parameter storage
     * \param error_messages the string containing all error messages
     * \return true if successfully set the memory info
     */
    inline bool set_memory(const std::string& input,
                           std::map<std::string, std::string>& message,
                           std::string& error_messages)
    {
        m_provided_memory = true;
        std::string in = input;
        if (message.find("memory") != message.end())
        { error_messages.append("Warning: Duplicated argument --memory\n"); }
        try
        {
            size_t memory = misc::convert<unsigned long long>(input);
            m_memory = memory;
        }
        catch (...)
        {
            // contain MB KB or B here
            if (input.length() >= 2)
            {
                try
                {
                    std::transform(in.begin(), in.end(), in.begin(), ::toupper);
                    std::string unit = in.substr(in.length() - 2);
                    std::string value = in.substr(0, in.length() - 2);
                    if (unit == "KB")
                    {
                        m_memory =
                            misc::convert<unsigned long long>(value) * 1024;
                        message["memory"] = misc::to_string(value) + "KB";
                    }
                    else if (unit == "MB")
                    {
                        m_memory = misc::convert<unsigned long long>(value)
                                   * 1024 * 1024;
                        message["memory"] = misc::to_string(value) + "MB";
                    }
                    else if (unit == "GB")
                    {
                        m_memory = misc::convert<unsigned long long>(value)
                                   * 1024 * 1024 * 1024;
                        message["memory"] = misc::to_string(value) + "GB";
                    }
                    else if (unit == "TB")
                    {
                        m_memory = misc::convert<unsigned long long>(value)
                                   * 1024 * 1024 * 1024 * 1024;
                        message["memory"] = misc::to_string(value) + "TB";
                    }
                    else
                    {
                        // maybe only one input?
                        unit = in.substr(in.length() - 1);
                        value = in.substr(0, in.length() - 1);
                        if (unit == "B")
                        {
                            m_memory = misc::convert<unsigned long long>(value);
                            message["memory"] = misc::to_string(value) + "B";
                        }
                        else if (unit == "K")
                        {
                            m_memory =
                                misc::convert<unsigned long long>(value) * 1024;
                            message["memory"] = misc::to_string(value) + "KB";
                        }
                        else if (unit == "M")
                        {
                            m_memory = misc::convert<unsigned long long>(value)
                                       * 1024 * 1024;
                            message["memory"] = misc::to_string(value) + "MB";
                        }
                        else if (unit == "G")
                        {
                            m_memory = misc::convert<unsigned long long>(value)
                                       * 1024 * 1024 * 1024;
                            message["memory"] = misc::to_string(value) + "GB";
                        }
                        else if (unit == "T")
                        {
                            m_memory = misc::convert<unsigned long long>(value)
                                       * 1024 * 1024 * 1024 * 1024;
                            message["memory"] = misc::to_string(value) + "TB";
                        }
                        std::cerr << "New memory: " << m_memory << std::endl;
                    }
                }
                catch (...)
                {
                    error_messages.append("Error: Undefined memory input: "
                                          + in);
                    return false;
                }
            }
            else
            {
                error_messages.append("Error: Undefined memory input: " + in);
                return false;
            }
        }

        message["memory"] = in;
        return true;
    }

    inline bool valid_distance(const std::string& str, const size_t& unit,
                               size_t& res, std::string& error_messages)
    {
        double cur_dist = misc::convert<double>(str) * unit;
        res = static_cast<size_t>(cur_dist);
        if (trunc(cur_dist) != cur_dist)
        {
            error_messages.append("Error: Non-integer distance obtained: "
                                  + misc::to_string(str) + " x "
                                  + misc::to_string(unit) + "\n");
        }
        return (trunc(cur_dist) == cur_dist);
    }
    inline size_t set_distance(const std::string& input,
                               const std::string& command, size_t default_unit,
                               std::map<std::string, std::string>& message,
                               bool& error, std::string& error_messages)
    {
        if (message.find(command) != message.end())
        {
            error_messages.append("Warning: Duplicated argument --" + command
                                  + "\n");
        }
        std::string in = input;
        std::transform(in.begin(), in.end(), in.begin(), ::tolower);
        const size_t u = 1000;
        message[command] = in;
        size_t dist = 0;
        bool cur_error = false;
        try
        {
            // when no unit is provided, we multiply based on default
            cur_error =
                !valid_distance(input, default_unit, dist, error_messages);
            if (cur_error)
            {
                error = true;
                return ~size_t(0);
            }
            if (default_unit == u)
                message[command] = input + "kb";
            else
                message[command] = input + "bp";
            return dist;
        }
        catch (...)
        {
            if (input.length() >= 2)
            {
                try
                {
                    std::string unit = in.substr(in.length() - 2);
                    size_t multiplication = 1;
                    std::string value = in.substr(0, in.length() - 2);
                    if (unit == "bp") {}
                    else if (unit == "kb")
                        multiplication = u;
                    else if (unit == "mb")
                        multiplication = u * u;
                    else if (unit == "gb")
                        multiplication = u * u * u;
                    else if (unit == "tb")
                        multiplication = u * u * u * u;
                    else
                    {
                        unit = in.substr(in.length() - 1);
                        value = in.substr(0, in.length() - 1);
                        if (unit == "b") {}
                        else if (unit == "k")
                            multiplication = u;
                        else if (unit == "m")
                            multiplication = u * u;
                        else if (unit == "g")
                            multiplication = u * u * u;
                        else if (unit == "t")
                            multiplication = u * u * u * u;
                        else
                        {
                            throw std::runtime_error("");
                        }
                    }
                    error |= !valid_distance(value, multiplication, dist,
                                             error_messages);
                    message[command] = value + unit;
                    return dist;
                }
                catch (...)
                {
                    error_messages.append("Error: Undefined distance input: "
                                          + in);
                    message.erase(command);
                    error = true;
                    return ~size_t(0);
                }
            }
            else
            {
                error_messages.append("Error: Undefined distance input: " + in);
                message.erase(command);
                error = true;
                return ~size_t(0);
            }
        }
        error = true;
        return ~size_t(0);
    }
    /*!
     * \brief Get the column index based on file header and the input string
     * \param target is the input string
     * \param ref is the vector containing the column header
     * \return the index to the column containing the column name. -1 if it is
     * not found
     */
    inline bool index_check(const std::string& target,
                            const std::vector<std::string>& ref, size_t& index,
                            bool case_sensitive = true) const
    {
        std::string tmp = "";
        for (size_t i = 0; i < ref.size(); ++i)
        {
            tmp = ref[i];
            if (!case_sensitive)
            {
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            }
            if (target == tmp)
            {
                index = i;
                return true;
            }
        }
        return false;
    }

    // return true if ok
    inline bool set_base_info_threshold(const std::vector<std::string>& ref,
                                        std::string& error_message)
    {
        std::vector<std::string> info = misc::split(m_info_col, ",");
        size_t index;
        bool contain = index_check(info[0], ref, index);
        m_base_has_col[+BASE_INDEX::INFO] = contain;
        if (contain) m_base_col_index[+BASE_INDEX::INFO] = index;
        if (!contain)
        {
            error_message.append("Warning: INFO field not found in base file,");
            error_message.append("will ignore INFO filtering\n");
            return true;
        }
        if (info.size() != 2)
        {
            error_message.append("Error: Invalid format of "
                                 "--base-info. Should be "
                                 "ColName,Threshold.\n");
            return false;
        }
        try
        {
            m_base_info_threshold = misc::convert<double>(info[1]);
            if (m_base_info_threshold < 0 || m_base_info_threshold > 1)
            {
                error_message.append("Error: Base INFO threshold "
                                     "must be within 0 and 1!\n");
                return false;
            }
        }
        catch (...)
        {
            error_message.append(
                "Error: Invalid argument passed to --base-info: " + m_info_col
                + "! Second argument must be numeric\n");
            return false;
        }
        return true;
    }
    // return true if valid
    inline bool set_base_maf_filter(const std::vector<std::string>& ref,
                                    std::string& error_message)
    {
        std::string maf_error = "Error: Invalid format of --base-maf. "
                                "Should be ColName,Threshold."
                                "or ColName,Threshold:ColName,Threshold.\n";
        std::vector<std::string> case_control = misc::split(m_maf_col, ":");
        if (case_control.size() > 2)
        {
            error_message.append(maf_error);
            return false;
        }
        bool control = true;
        size_t index = 0;
        bool contain = false;
        for (auto&& maf : case_control)
        {
            std::vector<std::string> detail = misc::split(maf, ",");
            if (control)
            {
                contain = index_check(detail[0], ref, index);
                m_base_has_col[+BASE_INDEX::MAF] = contain;
                if (contain) m_base_col_index[+BASE_INDEX::MAF] = index;
                if (!contain)
                {
                    error_message.append(
                        "Warning: MAF field not found in base file. "
                        "Will not perform MAF filtering on the base file\n");
                    return true;
                }
            }
            else
            {
                contain = index_check(detail[0], ref, index);
                if (contain) m_base_col_index[+BASE_INDEX::MAF_CASE] = index;
                m_base_has_col[+BASE_INDEX::MAF_CASE] = contain;
                if (!contain)
                {
                    error_message.append(
                        "Warning: Case MAF field not found in base file"
                        "Will not perform MAF filtering on the base file\n");
                    return true;
                }
            }
            double cur_maf;
            try
            {
                cur_maf = misc::convert<double>(detail[1]);
                if (cur_maf < 0 || cur_maf > 1)
                {
                    error_message.append("Error: Base MAF threshold must "
                                         "be within 0 and 1!\n");
                    return false;
                }
            }
            catch (...)
            {
                error_message.append(
                    "Error: Invalid argument passed to --base-maf: " + m_maf_col
                    + "! Threshold must be numeric\n");
                return false;
            }
            if (control) { m_control_maf_threshold = cur_maf; }
            else
            {
                m_case_maf_threshold = cur_maf;
            }
            control = false;
        }
        return true;
    }

    bool in_file(const std::vector<std::string>& column_names,
                 const std::string& field, const std::string& field_name,
                 std::map<std::string, std::string>& message, BASE_INDEX index,
                 bool case_sensitive = true)
    {
        size_t col_index;
        bool has_col =
            index_check(field, column_names, col_index, case_sensitive);
        if (has_col)
        {
            m_base_col_index[+index] = col_index;
            message[field_name] = field;
        }
        m_base_has_col[+index] = has_col;
        return has_col;
    }
};

#endif // COMMANDER_H
