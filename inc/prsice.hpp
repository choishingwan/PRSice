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
#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
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
// This should be the class to handle all the procedures
class PRSice
{
public:
    PRSice(const std::string& base_name, const Commander& commander,
           const bool prset, const size_t sample_ct, Reporter& reporter)
        : m_ignore_fid(commander.ignore_fid())
        , m_prset(prset)
        , m_logit_perm(commander.logit_perm())
        , m_num_perm(commander.permutation())
        , m_score(commander.get_score())
		, m_missing_score(commander.get_missing_score())
        , m_base_name(base_name)
        , m_target(commander.target_name())
        , m_out(commander.out())
        , m_target_binary(commander.is_binary())
    {

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
            // first check for ridiculously large sample size
            // allow 10 GB here
            if (CHAR_BIT * sample_ct > 1000000000) {
                m_perm_per_slice = 1;
            }
            else
            {
                // in theory, most of the time, perm_per_slice should be
                // equal to c_commander.num_permutation();

                int sample_memory = CHAR_BIT * sample_ct;
                m_perm_per_slice = 1000000000 / sample_memory;
                m_perm_per_slice = (m_perm_per_slice > m_num_perm)
                                       ? m_num_perm
                                       : m_perm_per_slice;
                // Additional slice to keep
                m_remain_slice = m_num_perm % m_perm_per_slice;
            }
            if (has_binary) {
                if (!m_logit_perm) {
                    std::string message =
                        "Warning: To speed up the permutation, "
                        "we perform  linear regression instead of logistic "
                        "regression within the permutation and uses the "
                        "p-value to rank the thresholds. Our assumptions "
                        "are as follow:\n";
                    message.append("         1. Linear Regression & Logistic "
                                   "Regression produce similar p-values\n");
                    message.append(
                        "         2. P-value is correlated with R2\n\n");
                    message.append("         If you must, you can run logistic "
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
    void run_prsice(const Commander& c_commander,
                    const std::string& region_name, const size_t pheno_index,
                    const size_t region_index, Genotype& target);
    void regress_score(Genotype &target, const double threshold, size_t thread,
                       const size_t pheno_index, const size_t iter_threshold);

    void prsice(const Commander& c_commander, const Region& c_region,
                const size_t c_pheno_index, bool prslice = false);
    void output(const Commander& c_commander, const Region& region,
                const size_t pheno_index, const size_t region_index,
                Genotype& target);
    void transpose_all(const Commander& c_commander, const Region& c_region,
                       size_t pheno_index) const;
    void summarize(const Commander& c_commander, Reporter& reporter);

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

    struct
    {
        std::vector<int> col;
        std::vector<std::string> name;
        std::vector<int> order;
        std::vector<bool> binary;
        bool use_pheno;
    } pheno_info;

    static std::mutex score_mutex;


    bool m_ignore_fid = false;
    bool m_prset = false;
    bool m_logit_perm = false;
    Eigen::VectorXd m_phenotype;
    Eigen::MatrixXd m_independent_variables;
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
    SCORING m_score = SCORING::AVERAGE;
    MISSING_SCORE m_missing_score = MISSING_SCORE::MEAN_IMPUTE;
    std::string m_log_file;
    std::string m_base_name;
    std::string m_target;
    std::string m_out;
    std::vector<bool> m_target_binary;
    std::vector<double> m_perm_result;
    std::vector<prsice_result> m_prs_results;
    std::vector<prsice_summary> m_prs_summary; // for multiple traits
    std::vector<double> m_best_sample_score;
    std::vector<size_t> m_sample_index;
    std::vector<size_t> m_significant_store{0, 0, 0}; // store the number of
                                                      // non-sig, margin sig,
                                                      // and sig pathway &
                                                      // phenotype

    std::unordered_map<std::string, size_t> m_sample_with_phenotypes;


    void thread_score(size_t region_start, size_t region_end, double threshold,
                      size_t thread, const size_t c_pheno_index,
                      const size_t iter_threshold);
    void thread_perm(Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& decomposed,
                     std::vector<Eigen::MatrixXd>& pheno_perm, size_t start,
                     size_t end, int rank, const Eigen::VectorXd& pre_se,
                     size_t processed, bool logit_perm);

    void permutation(Genotype &target, const size_t n_thread, bool logit_perm);
    void update_sample_included(Genotype &target);
    void gen_pheno_vec(Genotype &target, const std::string& pheno_file_name,
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
};

#endif // PRSICE_H
