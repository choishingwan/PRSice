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
    PRSice(std::string base_name, const Commander commander, bool prset)
        : m_base_name(base_name)
    {

        m_target = commander.target_name();
        m_target_binary = commander.is_binary();
        m_score = commander.get_scoring();
        m_ignore_fid = commander.ignore_fid();
        m_prset = prset;
        m_out = commander.out();
        bool perm = commander.permute();
        m_num_perm = commander.num_permutation();
        m_logit_perm = commander.logit_perm();
        m_seed = std::random_device()(); // cerr valgrind doesn't like this
        if (commander.seeded()) m_seed = commander.seed();
        fprintf(stderr, "Seed: %zu\n", m_seed);
        std::ofstream log_file_stream;
        log_file_stream.open(m_log_file.c_str(), std::ofstream::app);
        if (!log_file_stream.is_open()) {
        	std::string error_message =
        			"ERROR: Cannot open log file: " + m_log_file;
        	throw std::runtime_error(error_message);
        }
        log_file_stream << "Seed: "<< m_seed << std::endl;
        log_file_stream.close();
        bool has_binary = false;
        for (auto&& b : m_target_binary) {
            if (b) {
                has_binary = true;
                break;
            }
        }
        m_log_file = m_out + ".log";
        if (perm) {
            // first check for ridiculously large sample size
            // allow 10 GB here
            if (CHAR_BIT * m_sample_included.size() > 1000000000) {
                m_perm_per_slice = 1;
            }
            else
            {
                // in theory, most of the time, perm_per_slice should be
                // equal to c_commander.num_permutation();
                int sample_memory = CHAR_BIT * m_sample_included.size();
                m_perm_per_slice = 1000000000 / sample_memory;
                m_perm_per_slice = (m_perm_per_slice > m_num_perm)
                                       ? m_num_perm
                                       : m_perm_per_slice;
                // Additional slice to keep
                m_remain_slice = m_num_perm % m_perm_per_slice;
            }
            if (has_binary) {
                fprintf(stderr,
                        "\nWARNING: To speed up the permutation, we perform\n");
                fprintf(stderr,
                        "         linear regression instead of logistic\n");
                fprintf(
                    stderr,
                    "         regression within the permutation and uses\n");
                fprintf(stderr,
                        "         the p-value to rank the thresholds. Our "
                        "assumptions\n");
                fprintf(stderr, "         are as follow:\n");
                fprintf(stderr, "         1. Linear Regression & Logistic "
                                "Regression produce\n");
                fprintf(stderr, "            similar p-values\n");
                if (!m_logit_perm) {
                    fprintf(stderr,
                            "         2. P-value is correlated with R2\n\n");
                    fprintf(stderr,
                            "         If you must, you can run logistic "
                            "regression instead\n");
                    fprintf(stderr,
                            "         by setting the --logit-perm flag\n\n");
                }
                else
                {
                    fprintf(stderr, "         Using --logit-perm can be "
                                    "ridiculously slow\n");
                }
            }
        }
    };

    virtual ~PRSice();
    void pheno_check(const Commander& c_commander);
    // init_matrix whenever phenotype changes
    void init_matrix(const Commander& c_commander, const size_t pheno_index,
                     Genotype& target, const bool prslice = false);
    size_t num_phenotype() const
    {
        return (pheno_info.use_pheno) ? pheno_info.name.size() : 1;
    };
    void run_prsice(const Commander& c_commander, const std::string region_name,
                    const size_t pheno_index, const size_t region_index,
                    Genotype& target);
    void regress_score(const double threshold, size_t thread,
                       const size_t pheno_index, const size_t iter_threshold);

    void prsice(const Commander& c_commander, const Region& c_region,
                const size_t c_pheno_index, bool prslice = false);
    void output(const Commander& c_commander, const Region& region,
                const size_t pheno_index, const size_t region_index,
                Genotype& target);
    void transpose_all(const Commander& c_commander, const Region& c_region,
                       size_t pheno_index) const;
    void summarize(const Commander& c_commander);

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
    bool m_average_score = true;
    Eigen::VectorXd m_phenotype;
    Eigen::MatrixXd m_independent_variables;
    double m_null_r2 = 0.0;
    int m_best_index = -1;
    int m_perm_per_slice = 0;
    int m_remain_slice = 0;
    int m_num_perm = 0;
    unsigned int m_seed = 0;
    size_t m_num_snp_included = 0;
    size_t m_region_index = 0;
    size_t m_all_thresholds = 0;
    size_t m_max_fid_length = 3;
    size_t m_max_iid_length = 3;
    SCORING m_score = SCORING::MEAN_IMPUTE;
    std::string m_log_file;
    std::string m_base_name;
    std::string m_target;
    std::string m_out;
    std::vector<bool> m_target_binary;
    std::vector<double> m_perm_result;
    std::vector<prsice_result> m_prs_results;
    std::vector<prsice_summary> m_prs_summary; // for multiple traits
    std::vector<Sample_lite> m_current_sample_score;
    std::vector<Sample_lite> m_best_sample_score;
    std::vector<std::string> m_sample_included;
    std::vector<size_t> m_sample_index;
    std::vector<size_t> m_significant_store{0, 0, 0}; // store the number of
                                                      // non-sig, margin sig,
                                                      // and sig pathway &
                                                      // phenotype
    std::vector<Sample> m_sample_names; // might want to not storing it here
    std::unordered_map<std::string, size_t> m_sample_with_phenotypes;


    void thread_score(size_t region_start, size_t region_end, double threshold,
                      size_t thread, const size_t c_pheno_index,
                      const size_t iter_threshold);
    void thread_perm(Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& decomposed,
                     std::vector<Eigen::MatrixXd>& pheno_perm, size_t start,
                     size_t end, int rank, const Eigen::VectorXd& pre_se,
                     size_t processed, bool logit_perm);
    void permutation(unsigned int seed, int perm_per_slice, int remain_slice,
                     int total_permutation, int n_thread, bool logit_perm);
    void permutation(const int n_thread, bool logit_perm);
    void update_sample_included();
    void gen_pheno_vec(const std::string& pheno_file_name,
                       const int pheno_index, bool regress);
    std::vector<size_t> get_cov_index(const std::string& c_cov_file,
                                      std::vector<std::string>& cov_header);
    void gen_cov_matrix(const std::string& c_cov_file,
                        std::vector<std::string>& cov_header);
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
