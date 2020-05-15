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


#ifndef PRSICE_H
#define PRSICE_H

#include "commander.hpp"
#include "genotype.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "regression.hpp"
#include "reporter.hpp"
#include "snp.hpp"
#include "storage.hpp"
#include "thread_queue.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <errno.h>
#include <fstream>
#include <iomanip>
#include <map>
#include <math.h>
#include <mutex>
#include <random>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _WIN32
#include <process.h>
#include <windows.h>
//#define pthread_t HANDLE
//#define THREAD_RET_TYPE unsigned __stdcall
//#define THREAD_RETURN return 0
// we give an extra space for window just in case
#define NEXT_LENGTH 1LL
#else
#define NEXT_LENGTH 0LL
//#include <pthread.h>
#endif
#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/mach_host.h>
#include <mach/mach_init.h>
#include <mach/mach_types.h>
#include <mach/vm_statistics.h>
#endif
// This should be the class to handle all the procedures

class PRSice
{
public:
    PRSice() {}
    PRSice(const CalculatePRS& prs_info, const PThresholding& p_info,
           const Permutations& perm, const std::string& output,
           const size_t max_fid, const size_t max_iid, const bool binary,
           Reporter* reporter)
        : m_prefix(output)
        , m_max_fid_length(max_fid)
        , m_max_iid_length(max_iid)
        , m_binary_trait(binary)
        , m_reporter(reporter)
        , m_prs_info(prs_info)
        , m_p_info(p_info)
        , m_perm_info(perm)
    {
    }

    virtual ~PRSice();
    /*!
     * \brief This function will read in the phenotype information and determine
     * which phenotype to include
     * \param file_name is the name of the phenotype file
     * \param is the column name of the desired phenotypes
     * \param reporter is the logger
     */
    static void pheno_check(const bool no_regress, Phenotype& pheno,
                            Reporter& reporter);

    // init_matrix whenever phenotype changes
    /*!
     * \brief init_matrix will initialize the independent and dependent matrix
     * \param c_commander contains all user input
     * \param pheno_index is the index of the phenotype we are currently
     * processing
     * \param target is the target genotype. Provide the sample ID information
     * \param reporter is the logger
     */
    void init_matrix(const size_t pheno_index, const std::string& delim,
                     Genotype& target);
    void new_phenotype(Genotype& target);
    bool run_prsice(const size_t pheno_index, const size_t region_index,
                    const std::vector<std::vector<size_t>>& region_membership,
                    const bool all_scores, Genotype& target);
    /*!
     * \brief Before calling this function, the target should have loaded the
     * PRS. Then this function will fill in the m_independent_variable matrix
     * and call the required regression algorithms. It will then check if we
     * encounter a more significant result
     * \param target is the target genotype file containing the PRS information
     * \param threshold is the current p-value threshold, use for output
     * \param thread is the number of thread allowed
     * \param pheno_index is the index of the current phenotype
     * \param iter_threshold is the index of the current threshold
     */
    void regress_score(Genotype& target, const double threshold,
                       const int thread, const size_t pheno_index,
                       const size_t prs_result_idx);

    /*!
     * \brief Function responsible for generating the .prsice file
     * \param c_commander contains the user inputs
     * \param region contains the region information
     * \param pheno_index is the index of the current phenotype
     * \param region_index is teh index of the current region
     */
    void output(const std::vector<std::string>& region_names,
                const size_t pheno_index, const size_t region_index);
    /*!
     * \brief Function that prepare the output files by writing out white
     * spaces, this allow us to generate a nice vertical file
     * \param out is the output prefix
     * \param all_score indicate if we want to generate the all score file
     * \param has_prev indicate if prevalence is provided. If it is
     * provided, we need to also provide the adjusted R2 \param target is
     * the target genotype object containing the sample names \param
     * region_name is the name of regions involved \param pheno_index is the
     * index of the current phenotype
     */
    void prep_output(const Genotype& target,
                     const std::vector<std::string>& region_name,
                     const size_t pheno_index, const bool all_score);
    /*!
     * \brief This function will summarize all PRSice / PRSet results and
     * generate the .summary file
     * \param c_commander contains all user input
     * \param reporter is the logger
     */
    void summarize();
    /*!
     * \brief Calculate the number of processes required
     * \param commander the user input, provide information on the number of
     * permutation
     * \param num_region the number of region to process
     * \param num_thresholds the number of thresholds to process
     */
    void init_progress_count(const std::vector<std::set<double>>& thresholds)
    {
        const size_t num_region = thresholds.size();
        const bool set_perm = m_perm_info.run_set_perm;
        const bool perm = m_perm_info.run_perm;
        // competitive p-value calculation will have its own progress counts
        const size_t num_pheno_perm =
            set_perm ? 1 : (perm) ? m_perm_info.num_permutation + 1 : 1;
        const size_t num_set_perm = set_perm ? m_perm_info.num_permutation : 0;
        m_total_process = 0;
        for (size_t i = 0; i < thresholds.size(); ++i)
        {
            if (i == 1) continue; // ignore background
            m_total_process += thresholds[i].size() * num_pheno_perm;
        }
        // now calculate the number of competitive permutation we need
        // we only use the best threshold for set based permutation
        m_total_competitive_process = (num_region - 2) * num_set_perm;
        m_analysis_done = 0;
        m_total_competitive_perm_done = 0;
    }

    PRSice(const PRSice&) = delete;            // disable copying
    PRSice& operator=(const PRSice&) = delete; // disable assignment
    /*!
     * \brief responsible for printing the progress bar
     * \param completed is use for cheating. To print the 100% at the end
     */
    void print_progress(bool completed = false)
    {
        double cur_progress = (static_cast<double>(m_analysis_done)
                               / static_cast<double>(m_total_process))
                              * 100.0;
        // progress bar can be slow when permutation + thresholding is used due
        // to the huge amount of processing required
        if (cur_progress - m_previous_percentage > 0.01)
        {
            fprintf(stderr, "\rProcessing %03.2f%%", cur_progress);
            m_previous_percentage = cur_progress;
        }

        if (m_previous_percentage >= 100.0 || completed)
        { fprintf(stderr, "\rProcessing %03.2f%%", 100.0); }
    }
    /*!
     * \brief The master function for performing the competitive analysis
     * \param target is the target genotype object
     * \param commander contains all user inputs
     * \param pheno_index is the index of the current phenotype
     */
    void
    run_competitive(Genotype& target,
                    const std::vector<size_t>::const_iterator& bk_start_idx,
                    const std::vector<size_t>::const_iterator& bk_end_idx,
                    const size_t pheno_index);

    /*!
     * \brief Function responsible to generate the best score file
     * \param target is the target genotype, mainly for ID and in_regression
     * flag
     * \param pheno_index  the index of the current phenotype
     * \param  commander is the container of all user inputs
     */
    void print_best(Genotype& target,
                    const std::vector<std::string>& region_name,
                    const size_t pheno_index);

protected:
    static void parse_pheno_header(std::unique_ptr<std::istream> pheno_file,
                                   Phenotype& pheno_info, Reporter& reporter);
    static std::tuple<size_t, bool>
    get_pheno_idx(const std::vector<std::string_view>& column,
                  const Phenotype& pheno_info, const std::string& pheno);

    struct prsice_result
    {
        double threshold;
        double r2;
        double r2_adj;
        double coefficient;
        double p;
        double emp_p;
        double se;
        double competitive_p;
        size_t num_snp; // num snp should always be positive
    };
    struct prsice_summary
    {
        prsice_result result;
        std::string pheno;
        std::string set;
        double r2_null;
        double top;
        double bottom;
        double prevalence;
        bool has_competitive;
    };
    struct column_file_info
    {
        long long header_length;
        long long skip_column_length;
        long long line_width;
        long long processed_threshold;
        column_file_info()
        {
            header_length = 0;
            skip_column_length = 0;
            line_width = 0;
            processed_threshold = 0;
        }
    };

    // store the number of non-sig, margin sig, and sig pathway & phenotype
    static std::mutex lock_guard;
    // As R has a default precision of 7, we will go a bit
    // higher to ensure we use up all precision
    const std::string m_prefix;
    const long long m_precision = 9;
    // the 7 are:
    // 1 for sign
    // 1 for dot
    // 2 for e- (scientific)
    // 3 for exponent (max precision is somewhere around +-e297, so 3 is enough
    const long long m_numeric_width = m_precision + 7;
    const size_t m_max_fid_length = 3;
    const size_t m_max_iid_length = 3;
    const bool m_binary_trait = true;
    Eigen::MatrixXd m_independent_variables;
    // TODO: Use other method for faster best output
    Eigen::MatrixXd m_fast_best_output;
    Eigen::VectorXd m_phenotype;
    std::unordered_map<std::string, size_t> m_sample_with_phenotypes;
    std::vector<prsice_result> m_prs_results;
    std::vector<prsice_summary> m_prs_summary; // for multiple traits
    std::vector<double> m_perm_result;
    std::vector<double> m_permuted_pheno;
    std::vector<double> m_best_sample_score;
    std::vector<size_t> m_matrix_index;
    std::vector<size_t> m_significant_store {0, 0, 0};
    std::ofstream m_all_out, m_best_out, m_prsice_out;
    column_file_info m_all_file, m_best_file;
    std::string m_out;
    std::mutex m_thread_mutex;
    double m_previous_percentage = -1.0;
    double m_null_r2 = 0.0;
    double m_null_p = 1.0;
    double m_null_se = 0.0;
    double m_null_coeff = 0.0;
    size_t m_total_process = 0;
    size_t m_total_competitive_process = 0;
    size_t m_total_competitive_perm_done = 0;
    uint32_t m_num_snp_included = 0;
    uint32_t m_analysis_done = 0;


    int m_best_index = -1;
    bool m_quick_best = true;
    bool m_printed_warning = false;
    Reporter* m_reporter;
    CalculatePRS m_prs_info;
    PThresholding m_p_info;
    Permutations m_perm_info;
    Phenotype m_pheno_info;

    // Functions
    static bool empty_name(const std::string& in) { return in.empty(); }
    /*!
     * \brief permutation is the master function to call the subfunctions
     * responsible for calculating the permuted t-value
     * \param n_thread indicate the number of threads allowed
     * \param is_binary indicate if the current phenotype is binary
     */
    void permutation(const int n_thread, const bool is_binary);
    /*!
     * \brief This function will calculate the maximum length of the FID and
     * IID, generate the matrix index for quicker search and also set the in
     * regression flag for each sample
     * \param target is the target genotype object
     */
    void update_sample_included(const std::string& delim, const bool binary,
                                Genotype& target);
    void slow_print_best(Genotype& target, const size_t pheno_index);
    /*!
     * \brief gen_pheno_vec is the function responsible for generating the
     * phenotype vector
     * \param target is the target genotype object, providing information on
     * FID, IID and also phenotype (e.g. from fam file)
     * \param pheno_file_name contains the name of the phenotype file
     * \param pheno_index contain the index of the current phenotype
     * \param reporter is the logger
     */
    void gen_pheno_vec(Genotype& target, const size_t pheno_index,
                       const std::string& delim);
    std::tuple<std::vector<double>, size_t, size_t, int> process_phenotype_file(
        const std::string& file_name, const std::string& delim,
        const std::size_t pheno_idx, const bool ignore_fid, Genotype& target);
    std::tuple<std::vector<double>, size_t, int>
    process_phenotype_info(const std::string& delim, Genotype& target);
    std::tuple<bool, size_t, size_t>
    binary_pheno_is_valid(const int max_pheno_code,
                          std::vector<double>& pheno_store);
    bool quantitative_pheno_is_valid(const std::vector<double>& pheno_store);
    void print_pheno_log(const std::string& name, const size_t sample_ct,
                         const size_t num_not_found, const size_t invalid_pheno,
                         const int max_pheno_code, const bool ignore_fid,
                         std::vector<double>& pheno_store);
    /*!
     * \brief Function to generate the m_independent_variable matrix
     * \param c_cov_file is the name of the covariate file
     * \param cov_header_name is a string vector containing the name of each
     * covaraites
     * \param cov_header_index is a vector of int containing the column index of
     * each covariates
     * \param factor_cov_index is a vector of int containing the column index of
     * each factor covariates
     * \param reporter is the logger
     */
    void gen_cov_matrix(const std::string& delim);
    void adjustment_factor(const double prevalence, double& top,
                           double& bottom);
    void print_na(const std::string& region_name, const double threshold,
                  const size_t num_snp, const bool has_prevalence);
    void store_best(const std::string& pheno_name,
                    const std::string& region_name, const double top,
                    const double bottom, const double prevalence,
                    const bool is_base);
    /*!
     * \brief Function use to process the covariate file, should be able to
     * determine the level of factors
     * \param cov_file is the name of the covariate file
     * \param factor_cov_index is the column index for factor covariates
     * \param cov_start_index is the starting position of each column (return)
     * \param cov_index is the column index for all covariates
     * \param cov_name is the name of each covariates
     * \param factor_levels is a structure to store the factor levels. It's size
     * should equal to the number of factor level. The nested map should contain
     * the variable to level matching information
     * \param num_column is the number of column required (return)
     * \param reporter is the logger
     */
    void process_cov_file(
        std::vector<size_t>& cov_start_index,
        std::vector<std::unordered_map<std::string, size_t>>& factor_levels,
        Eigen::Index& num_column, const std::string& delim);
    bool validate_covariate(const std::string& covariate,
                            const size_t num_factors, const size_t idx,
                            size_t& factor_level_idx,
                            std::vector<size_t>& missing_count);
    void update_sample_matrix(
        const std::vector<size_t>& missing_count,
        std::vector<std::pair<std::string, size_t>>& valid_sample_index);
    void get_se_matrix(
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& PQR,
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType&
            Pmat,
        const Eigen::MatrixXd& Rinv, const Eigen::Index p,
        const Eigen::Index rank, Eigen::VectorXd& se_base);
    void observe_set_perm(Thread_Queue<size_t>& progress_observer,
                          size_t total_perm);
    double get_coeff_resid_norm(const Regress& decomposed,
                                const Eigen::VectorXd& prs,
                                Eigen::VectorXd& beta, Eigen::VectorXd effects);
    template <typename T>
    void subject_set_perm(T& progress_observer, Genotype& target,
                          std::vector<size_t> background,
                          std::map<size_t, std::vector<size_t>>& set_index,
                          std::vector<std::atomic<size_t>>& set_perm_res,
                          const std::vector<double>& obs_t_value,
                          const std::random_device::result_type seed,
                          const Regress& decomposed, const size_t num_perm,
                          const bool is_binary);
    /*!
     * \brief Once PRS analysis and permutation has been performed for all
     * p-value thresholds we will run this function to calculate the
     * empirical p-value
     */
    void process_permutations();


    /*!
     * \brief Function to generate PRS for null set when multiple threading is
     * used
     * \param q is teh queue used to communicate with the consumer
     * \param target is the target genotype, responsible for the generation of
     * PRS
     * \param num_consumer is the number of consumer. use for restricting the
     * number of PRS read in at one time
     * \param set_index is the dictionary containing the sizes of sets
     * \param num_perm is the number of permutation to erpfrom
     * \param require_standardize is a boolean, indicating if we want a
     * standardized PRS
     */
    void
    produce_null_prs(Thread_Queue<std::pair<std::vector<double>, size_t>>& q,
                     Genotype& target, std::vector<size_t> background,
                     size_t num_consumer,
                     std::map<size_t, std::vector<size_t>>& set_index);
    /*!
     * \brief This is the "consumer" function responsible for reading in the PRS
     * and perform the regression analysis
     * \param q is the queue used for communication between the producer and
     * consumer
     * \param set_index is the dictionary containing index to ori_t_value for
     * sets with size specified in the key
     * \param ori_t_value contain the observed t-statistic for the  sets
     * \param set_perm_res is the vector storing the result of permutation.
     * Counting the number of time the permuted T is bigger than the observed T
     * for a specific set
     * \param is_binary indicate if the phenotype is binary or not
     */
    void consume_prs(Thread_Queue<std::pair<std::vector<double>, size_t>>& q,
                     const Regress& decomposed,
                     std::map<size_t, std::vector<size_t>>& set_index,
                     const std::vector<double>& obs_t_value,
                     std::vector<std::atomic<size_t>>& set_perm_res,
                     const bool is_binary);

    void null_set_no_thread(
        Genotype& target, const size_t num_background,
        std::vector<size_t> background,
        const std::map<size_t, std::vector<size_t>>& set_index,
        const Regress& decomposed, std::vector<double>& obs_t_value,
        std::vector<std::atomic<size_t>>& set_perm_res, const bool is_binary);
    double get_t_value(const Regress& decomposed, const Eigen::VectorXd& prs,
                       Eigen::VectorXd& beta, Eigen::VectorXd& effects,
                       double& coefficient, double& standard_error);
    /*!
     * \brief The "producer" for generating the permuted phenotypes
     * \param q is the queue for contacting the consumers
     * \param num_consumer is the number of consumer
     */
    void gen_null_pheno(Thread_Queue<std::pair<Eigen::VectorXd, size_t>>& q,
                        size_t num_consumer);
    /*!
     * \brief The "consumer" for calculating the T-value on permuted phenotypes
     * \param q is the queue where the producer generated the permuted phenotype
     * \param decomposed is the pre-computed decomposition
     * \param rank is the pre-computed rank
     * \param pre_se is the pre-computed SE matrix
     * \param run_glm is a boolean indicate if we want to run logistic
     * regression
     */
    void consume_null_pheno(
        Thread_Queue<std::pair<Eigen::VectorXd, size_t>>& q,
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& PQR,
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType&
            Pmat,
        const Eigen::MatrixXd& R, bool run_glm);
    /*!
     * \brief Funtion to perform single threaded permutation
     * \param decomposed is the pre-decomposed independent matrix. If run glm is
     * true, this will be ignored
     * \param rank is the rank of the decomposition
     * \param pre_se is the pre-computed SE matrix required for calculating the
     * final SE
     * \param run_glm indicate if we want to run GLM instead of using
     * precomputed matrix
     */
    void run_null_perm_no_thread(
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& PQR,
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType&
            Pmat,
        const Eigen::MatrixXd& R, const bool run_glm);

    void parse_pheno(const std::string& pheno, std::vector<double>& pheno_store,
                     int& max_pheno_code);

    std::unordered_map<std::string, std::string>
    load_pheno_map(const std::string& delim, const size_t idx,
                   const bool ignore_fid,
                   std::unique_ptr<std::istream> pheno_file);

    void reset_result_containers(const Genotype& target,
                                 const size_t region_idx);

    void fisher_yates(std::vector<size_t>& idx, std::mt19937& g, size_t n);
    template <typename T>
    class dummy_reporter
    {
        PRSice& m_parent;
        bool m_completed = false;

    public:
        dummy_reporter(PRSice& p) : m_parent(p) {}
        void emplace(T&& /*item*/)
        {
            ++m_parent.m_analysis_done;
            m_parent.print_progress();
        }
        void completed() { m_completed = true; }
    };
};

#endif // PRSICE_H
