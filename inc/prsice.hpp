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

#ifndef PRSICE_H
#define PRSICE_H

#include "commander.hpp"
#include "genotype.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "region.hpp"
#include "regression.hpp"
#include "reporter.hpp"
#include "snp.hpp"
#include "storage.hpp"
#include "thread_queue.hpp"
#include <Eigen/Dense>
#include <algorithm>
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
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _WIN32
#include <mingw.mutex.h>
#include <mingw.thread.h>
#include <process.h>
#include <windows.h>
//#define pthread_t HANDLE
//#define THREAD_RET_TYPE unsigned __stdcall
//#define THREAD_RETURN return 0
// we give an extra space for window just in case
#define NEXT_LENGTH 1
#else
#include <thread>
#define NEXT_LENGTH 0
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
    /*!
     * \brief PRSice constructor
     * \param commander contains all user input
     * \param reporter is the logger
     */
    PRSice(const Commander& commander, const bool prset, Reporter& reporter)
        : m_target_binary(commander.is_binary())
        , m_target(commander.target_name())
        , m_out(commander.out())
        , m_score(commander.get_score())
        , m_missing_score(commander.get_missing_score())
        , m_ignore_fid(commander.ignore_fid())
        , m_perform_prset(prset)
    {
        m_perform_perm = commander.num_perm(m_num_perm);
        m_logit_perm = commander.logit_perm();
        m_seed = commander.seed();
        bool has_binary = false;
        for (auto b : m_target_binary) {
            if (b) {
                has_binary = true;
                break;
            }
        }
        if (m_perform_perm) {
            if (has_binary) {
                if (!m_logit_perm) {
                    std::string message =
                        "Warning: To speed up the permutation, "
                        "we perform linear regression instead of logistic "
                        "regression within the permutation and uses the "
                        "p-value to rank the thresholds. Our assumptions "
                        "are as follow:\n";
                    message.append("1) Linear Regression & Logistic "
                                   "Regression produce similar p-values\n");
                    message.append("2) P-value is correlated with R2\n\n");
                    message.append("If you must, you can run logistic "
                                   "regression instead by setting the "
                                   "--logit-perm flag\n\n");
                    reporter.report(message);
                }
                else
                {
                    std::string message = "Warning: Using --logit-perm can be "
                                          "ridiculously slow\n";
                    reporter.report(message);
                }
            }
        }
    }

    virtual ~PRSice();
    /*!
     * \brief This function will read in the phenotype information and determine
     * which phenotype to include
     * \param file_name is the name of the phenotype file
     * \param is the column name of the desired phenotypes
     * \param reporter is the logger
     */
    void pheno_check(const std::string& file_name,
                     const std::vector<std::string>& col_names,
                     const std::vector<bool>& is_binary, Reporter& reporter);
    // init_matrix whenever phenotype changes
    /*!
     * \brief init_matrix will initialize the independent and dependent matrix
     * \param c_commander contains all user input
     * \param pheno_index is the index of the phenotype we are currently
     * processing
     * \param target is the target genotype. Provide the sample ID information
     * \param reporter is the logger
     */
    void init_matrix(const Commander& c_commander, const size_t pheno_index,
                     Genotype& target, Reporter& reporter);
    /*!
     * \brief Return the total number of phenotype involved
     * \return the total number of phenotype to process
     */
    size_t num_phenotype() const
    {
        return (pheno_info.use_pheno) ? pheno_info.name.size() : 1;
    }
    void run_prsice(const Commander& c_commander, const size_t pheno_index,
                    const size_t region_index,
                    const std::vector<size_t>& region_membership,
                    const std::vector<size_t>& region_start_idx,
                    Genotype& target);
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
    void regress_score(Genotype& target, const double threshold, size_t thread,
                       const size_t pheno_index, const size_t prs_result_idx);

    /*!
     * \brief Function responsible for generating the .prsice file
     * \param c_commander contains the user inputs
     * \param region contains the region information
     * \param pheno_index is the index of the current phenotype
     * \param region_index is teh index of the current region
     */
    void output(const Commander& c_commander,
                const std::vector<std::string>& region_names,
                const size_t pheno_index, const size_t region_index);
    /*!
     * \brief Function that prepare the output files by writing out white
     * spaces, this allow us to generate a nice vertical file
     * \param out is the output prefix
     * \param all_score indicate if we want to generate the all score file
     * \param has_prev indicate if prevalence is provided. If it is provided, we
     * need to also provide the adjusted R2
     * \param target is the target genotype object containing the sample names
     * \param region_name is the name of regions involved
     * \param pheno_index is the index of the current phenotype
     */
    void prep_output(const std::string& out, const bool all_score,
                     const bool has_prev, const Genotype& target,
                     const std::vector<std::string>& region_name,
                     const size_t pheno_index);
    /*!
     * \brief This function will summarize all PRSice / PRSet results and
     * generate the .summary file
     * \param c_commander contains all user input
     * \param reporter is the logger
     */
    void summarize(const Commander& c_commander, Reporter& reporter);
    /*!
     * \brief Calculate the number of processes required
     * \param commander the user input, provide information on the number of
     * permutation
     * \param num_region the number of region to process
     * \param num_thresholds the number of thresholds to process
     */
    void init_process_count(const Commander& commander, size_t num_region,
                            size_t num_thresholds)
    {
        size_t num_perm;
        const bool perm = commander.num_perm(num_perm);
        size_t num_set_perm;
        const bool set_perm = commander.set_perm(num_set_perm);
        // the number of raw PRSice run
        m_total_process = num_thresholds * num_phenotype()
                          * ((num_region > 2) ? num_region - 1 : 1);
        if (perm) {
            m_total_process *= (num_perm + 1);
        }
        if (set_perm) {
            // the additional permutation we've got to run, num_region -2 as we
            // don't perform permutation on the background set nor the base set
            m_total_process += num_phenotype()
                               * (num_region - 2)
                               * num_set_perm;
        }
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
        if (cur_progress - m_previous_percentage > 0.01) {
            fprintf(stderr, "\rProcessing %03.2f%%", cur_progress);
            m_previous_percentage = cur_progress;
        }

        if (m_previous_percentage >= 100.0 || completed) {
            fprintf(stderr, "\rProcessing %03.2f%%", 100.0);
        }
    }
    /*!
     * \brief The master function for performing the competitive analysis
     * \param target is the target genotype object
     * \param commander contains all user inputs
     * \param pheno_index is the index of the current phenotype
     */
    void run_competitive(Genotype& target,
                         const std::vector<size_t>::const_iterator & bk_start_idx,
                         const std::vector<size_t>::const_iterator& bk_end_idx, const Commander& commander,
                         const size_t pheno_index, Reporter &reporter);


protected:
private:
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
        int header_length;
        int skip_column_length;
        int line_width;
        int processed_threshold;
    };
    struct Pheno_Info
    {
        std::vector<int> col;
        std::vector<std::string> name;
        std::vector<int> order;
        std::vector<bool> binary;
        bool use_pheno;
    } pheno_info;

    // store the number of non-sig, margin sig, and sig pathway & phenotype
    static std::mutex lock_guard;

    Eigen::MatrixXd m_independent_variables;
    Eigen::VectorXd m_phenotype;
    std::unordered_map<std::string, size_t> m_sample_with_phenotypes;
    std::vector<prsice_result> m_prs_results;
    std::vector<prsice_summary> m_prs_summary; // for multiple traits
    std::vector<double> m_perm_result;
    std::vector<double> m_permuted_pheno;
    std::vector<double> m_best_sample_score;
    std::vector<size_t> m_matrix_index;
    std::vector<size_t> m_significant_store{0, 0, 0};
    std::vector<bool> m_target_binary;
    column_file_info m_all_file, m_best_file;
    std::string m_target;
    std::string m_out;
    std::mutex m_thread_mutex;
    double m_previous_percentage = -1.0;
    double m_null_r2 = 0.0;
    double m_null_p = 1.0;
    double m_null_se = 0.0;
    double m_null_coeff = 0.0;
    std::random_device::result_type m_seed = 0;
    size_t m_total_process = 0;
    uint32_t m_num_snp_included = 0;
    uint32_t m_analysis_done = 0;
    // As R has a default precision of 7, we will go a bit
    // higher to ensure we use up all precision
    int32_t m_precision = 9;
    // the 7 are:
    // 1 for sign
    // 1 for dot
    // 2 for e- (scientific)
    // 3 for exponent (max precision is somewhere around +-e297, so 3 is enough
    int32_t m_numeric_width = m_precision + 7;
    int32_t m_max_fid_length = 3;
    int32_t m_max_iid_length = 3;
    int m_best_index = -1;
    size_t m_num_perm = 0;
    SCORING m_score = SCORING::AVERAGE;
    MISSING_SCORE m_missing_score = MISSING_SCORE::MEAN_IMPUTE;
    bool m_ignore_fid = false;
    bool m_perform_prset = false;
    bool m_perform_competitive = false;
    bool m_perform_perm = false;
    bool m_logit_perm = false;
    // Functions

    /*!
     * \brief permutation is the master function to call the subfunctions
     * responsible for calculating the permuted t-value
     * \param n_thread indicate the number of threads allowed
     * \param is_binary indicate if the current phenotype is binary
     */
    void permutation(const size_t n_thread, bool is_binary);
    /*!
     * \brief This function will calculate the maximum length of the FID and
     * IID, generate the matrix index for quicker search and also set the in
     * regression flag for each sample
     * \param target is the target genotype object
     */
    void update_sample_included(Genotype& target);
    /*!
     * \brief gen_pheno_vec is the function responsible for generating the
     * phenotype vector
     * \param target is the target genotype object, providing information on
     * FID, IID and also phenotype (e.g. from fam file)
     * \param pheno_file_name contains the name of the phenotype file
     * \param pheno_index contain the index of the current phenotype
     * \param reporter is the logger
     */
    void gen_pheno_vec(Genotype& target, const std::string& pheno_file_name,
                       const size_t pheno_index, Reporter& reporter);
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
    void gen_cov_matrix(const std::string& c_cov_file,
                        const std::vector<std::string> cov_header_name,
                        const std::vector<uint32_t> cov_header_index,
                        const std::vector<uint32_t> factor_cov_index,
                        Reporter& reporter);
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
        const std::string& cov_file,
        const std::vector<uint32_t>& factor_cov_index,
        std::vector<size_t>& cov_start_index,
        const std::vector<uint32_t>& cov_index,
        const std::vector<std::string>& cov_name,
        std::vector<std::unordered_map<std::string, size_t>>& factor_levels,
        size_t& num_column, Reporter& reporter);

    /*!
     * \brief Once PRS analysis and permutation has been performed for all
     * p-value thresholds we will run this function to calculate the empirical
     * p-value
     */
    void process_permutations();

    /*!
     * \brief Function responsible to generate the best score file
     * \param target is the target genotype, mainly for ID and in_regression
     * flag
     * \param pheno_index  the index of the current phenotype
     * \param  commander is the container of all user inputs
     */
    void print_best(Genotype& target, const size_t pheno_index,
                    const Commander& commander);

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
    produce_null_prs(Thread_Queue<std::pair<std::vector<double>, size_t> > &q,
                     Genotype& target,
                     const std::vector<size_t>::const_iterator &bk_start_idx,
                     const std::vector<size_t>::const_iterator &bk_end_idx, size_t num_consumer,
                     std::map<size_t, std::vector<size_t> > &set_index,
                     const size_t num_perm, const bool require_standardize, const bool use_ref_maf);
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
    void consume_prs(Thread_Queue<std::pair<std::vector<double>, size_t> > &q,
                     std::map<size_t, std::vector<size_t> > &set_index,
                     std::vector<double>& obs_t_value,
                     std::vector<size_t> &set_perm_res, const bool is_binary);
    /*!
     * \brief Function responsible for running the permutation required for
     * computing the competitive p-value
     * \param target is the target genotype object
     * \param set_index  is the dictionary for index of sets with the same size
     * \param ori_t_value contain the observed T-value for each set
     * \param set_perm_res is the vector storing the result of permutation.
     * Counting the number of time the permuted T is bigger than the observed T
     * for a specific set
     * \param num_perm is the number of permutation to perform
     * \param is_binary indicate if the phenotype is binary
     * \param require_standardize indicate if we require the standardization of
     * the genotype
     */
    void null_set_no_thread(Genotype& target,
                            const std::vector<size_t>::const_iterator &bk_start_idx,
                            const std::vector<size_t>::const_iterator &bk_end_idx,
                            const std::map<size_t, std::vector<size_t>>& set_index,
                            std::vector<double>& obs_t_value,
                            std::vector<size_t>& set_perm_res,
                            const size_t num_perm, const bool is_binary,
                            const bool require_standardize, const bool use_ref_maf);

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
    void
    consume_null_pheno(Thread_Queue<std::pair<Eigen::VectorXd, size_t>>& q,
                       const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> &decomposed,
                       int rank, const Eigen::VectorXd& pre_se, bool run_glm);
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
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& decomposed,
        const Eigen::Index rank, const Eigen::VectorXd& pre_se,
        const bool run_glm);
};

#endif // PRSICE_H
