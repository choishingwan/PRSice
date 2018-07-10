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
    PRSice(const std::string& base_name, const Commander& commander,
           const bool prset, const size_t sample_ct, Reporter& reporter)
        : m_ignore_fid(commander.ignore_fid())
        , m_prset(prset)
        , m_num_perm(commander.permutation())
        , m_score(commander.get_score())
        , m_missing_score(commander.get_missing_score())
        , m_base_name(base_name)
        , m_target(commander.target_name())
        , m_out(commander.out())
        , m_target_binary(commander.is_binary())
    {

        m_logit_perm = commander.logit_perm();
        // we calculate the number of permutation we can run at one time
        bool perm = (commander.permutation() > 0);
        m_seed = commander.seed();
        // check if there's any binary phenotype
        bool has_binary = false;
        for (auto&& b : m_target_binary) {
            if (b) {
                has_binary = true;
                break;
            }
        }
        if (perm) {
            // instead of reserving the required memory, we should just adjust
            // the thread number, though most likely it will be light on memory
            // now

            // gen_perm_memory(commander, sample_ct, reporter);

            // Additional slice to keep
            // DEBUG here
            // m_remain_slice = m_num_perm % m_perm_per_slice;
            if (has_binary) {
                if (!m_logit_perm) {
                    std::string message =
                        "Warning: To speed up the permutation, "
                        "we perform  linear regression instead of logistic "
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
    };

    virtual ~PRSice();
    void pheno_check(const Commander& c_commander, Reporter& reporter);
    // init_matrix whenever phenotype changes
    void init_matrix(const Commander& c_commander, const size_t pheno_index,
                     Genotype& target, Reporter& reporter,
                     const bool prslice = false);
    size_t num_phenotype() const
    {
        return (pheno_info.use_pheno) ? pheno_info.name.size() : 1;
    };
    void run_prsice(const Commander& c_commander, const Region& region,
                    const size_t pheno_index, const size_t region_index,
                    Genotype& target);
    void regress_score(Genotype& target, const double threshold, size_t thread,
                       const size_t pheno_index, const size_t iter_threshold);

    void prsice(const Commander& c_commander, const Region& c_region,
                const size_t c_pheno_index, bool prslice = false);
    void output(const Commander& c_commander, const Region& region,
                const size_t pheno_index, const size_t region_index,
                Genotype& target);
    void prep_output(const Commander& commander, Genotype& target,
                     std::vector<std::string> region_name,
                     const size_t pheno_index);
    void transpose_all(const Commander& c_commander, const Region& c_region,
                       size_t pheno_index) const;
    void summarize(const Commander& c_commander, Reporter& reporter);
    void init_process_count(const Commander& commander, size_t num_region,
                            size_t num_thresholds)
    {
        const bool perm = (commander.permutation() > 0);
        const bool set_perm = (commander.perform_set_perm());
        // the number of raw PRSice run
        m_total_process = num_thresholds * num_phenotype()
                          * ((num_region > 1) ? num_region - 1 : 1);
        if (perm) {
            m_total_process *= (commander.permutation() + 1);
        }
        if (set_perm) {
            // the additional permutation we've got to run, num_region -2 as we
            // don't perform permutation on the background set nor the base set
            m_total_process +=
                num_phenotype() * (num_region - 2) * (commander.set_perm());
        }
    }
    PRSice(const PRSice&) = delete;            // disable copying
    PRSice& operator=(const PRSice&) = delete; // disable assignment
    void print_progress(bool completed = false)
    {
        double cur_progress =
            ((double) m_analysis_done / (double) m_total_process) * 100.0;
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
        int num_snp;
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
    };

    struct Pheno_Info
    {
        std::vector<int> col;
        std::vector<std::string> name;
        std::vector<int> order;
        std::vector<bool> binary;
        bool use_pheno;
    } pheno_info;


    bool m_ignore_fid = false;
    bool m_prset = false;
    bool m_logit_perm = false;

    double m_null_r2 = 0.0;
    double m_null_p = 1.0;
    double m_null_se = 0.0;
    double m_null_coeff = 0.0;
    int m_best_index = -1;
    int m_perm_per_slice = 0;
    int m_remain_slice = 0;
    unsigned int m_seed = 0;
    size_t m_num_perm = 0;
    size_t m_num_snp_included = 0;
    size_t m_region_index = 0;
    size_t m_all_thresholds = 0;
    size_t m_max_fid_length = 3;
    size_t m_max_iid_length = 3;
    size_t m_total_process = 0;
    size_t m_analysis_done = 0;
    double m_previous_percentage = -1.0;
    SCORING m_score = SCORING::AVERAGE;
    MISSING_SCORE m_missing_score = MISSING_SCORE::MEAN_IMPUTE;
    std::string m_log_file;
    std::string m_base_name;
    std::string m_target;
    std::string m_out;
    std::vector<bool> m_target_binary;
    std::vector<prsice_result> m_prs_results;
    std::vector<prsice_summary> m_prs_summary; // for multiple traits
    std::vector<double> m_best_sample_score;
    std::vector<size_t> m_significant_store{0, 0, 0}; // store the number of
                                                      // non-sig, margin sig,
                                                      // and sig pathway &
                                                      // phenotype
    std::unordered_map<int, std::vector<double>> m_null_store;
    std::unordered_map<std::string, size_t> m_sample_with_phenotypes;
    static std::mutex lock_guard;

    struct column_file_info
    {
        size_t header_length;
        size_t skip_column_length;
        size_t line_width;
        size_t processed_threshold;
    };

    column_file_info m_all_file, m_best_file, m_snp_output;
    // As R has a default precision of 7, we will go a bit
    // higher to ensure we use up all precision
    size_t m_precision = 9;
    // the 7 are:
    // 1 for sign
    // 1 for dot
    // 2 for e- (scientific)
    // 3 for exponent (max precision is somewhere around +-e297, so 3 is enough
    size_t m_numeric_width = m_precision + 7;
    Eigen::MatrixXd m_independent_variables;
    Eigen::VectorXd m_phenotype;
    std::vector<double> m_perm_result;
    std::vector<double> m_permuted_pheno;
    std::mutex m_thread_mutex;
    // Functions
    void thread_score(size_t region_start, size_t region_end, double threshold,
                      size_t thread, const size_t c_pheno_index,
                      const size_t iter_threshold);
    void thread_perm(Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& decomposed,
                     size_t start, size_t end, int rank,
                     const Eigen::VectorXd& pre_se, size_t processed);
    void permutation(Genotype& target, const size_t n_thread, bool is_binary);
    void update_sample_included(Genotype& target);
    void gen_pheno_vec(Genotype& target, const std::string& pheno_file_name,
                       const int pheno_index, bool regress, Reporter& reporter);
    std::vector<size_t> get_cov_index(const std::string& c_cov_file,
                                      std::vector<std::string>& cov_header,
                                      Reporter& reporter);
    void gen_cov_matrix(const std::string& c_cov_file,
                        std::vector<std::string>& cov_header,
                        Reporter& reporter);
    void check_factor_cov(
        const std::string& c_cov_file,
        const std::vector<std::string>& c_cov_header,
        const std::vector<size_t>& cov_index,
        std::vector<std::unordered_map<std::string, int>>& factor_levels);
    // This should help us to update the m_prs_results
    void process_permutations();
    void summary();
    void gen_perm_memory(const Commander& commander, const size_t sample_ct,
                         Reporter& reporter);
    void print_best(Genotype& target, const size_t pheno_index,
                    const Commander& commander);
    void run_competitive(Genotype& target, const Commander& commander,
                         const size_t num_snp, const bool store_null,
                         const bool binary);
    void produce_null_prs(Thread_Queue<std::vector<double>>& q,
                          Genotype& target, std::vector<int>& sample_index,
                          size_t num_consumer, size_t num_perm, size_t set_size,
                          size_t num_selected_snps, double original_p,
                          bool require_standardize);
    void consume_prs(Thread_Queue<std::vector<double>>& q, double original_p,
                     int& num_significant, bool is_binary, bool store_p);
    void null_set_no_thread(Genotype& target, std::vector<int>& sample_index,
                            int& num_significant, size_t num_perm,
                            size_t set_size, size_t num_selected_snps,
                            double original_p, bool require_standardize,
                            bool is_binary, bool store_p);
    void gen_null_pheno(Thread_Queue<std::pair<Eigen::VectorXd, size_t>>& q,
                        size_t num_consumer);

    void
    consume_null_pheno(Thread_Queue<std::pair<Eigen::VectorXd, size_t>>& q,
                       Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& decomposed,
                       int rank, const Eigen::VectorXd& pre_se, bool run_glm);
    void run_null_perm_no_thread(
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& decomposed, int rank,
        const Eigen::VectorXd& pre_se, bool run_glm);
};

#endif // PRSICE_H
