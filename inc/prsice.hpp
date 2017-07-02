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

#include <algorithm>
#include <chrono>
#include <Eigen/Dense>
#include <fstream>
#include <math.h>
#include <map>
#include <mutex>
#include <random>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <thread>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "commander.hpp"
#include "genotype.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "region.hpp"
#include "regression.hpp"
#include "snp.hpp"
#include "storage.hpp"
//This should be the class to handle all the procedures
class PRSice
{
public:

    PRSice(std::string base_name, std::string target, std::vector<bool> target_binary,
    		SCORING score, size_t num_region, bool ignore_fid):
    			m_base_name(base_name), m_target(target), m_target_binary(target_binary),
				m_score(score), m_region_size(num_region), m_ignore_fid(ignore_fid){ };
    virtual ~PRSice();
    void pheno_check(const Commander &c_commander);
    void init_matrix(const Commander &c_commander, const size_t pheno_index, Genotype &target,
    		const bool prslice=false);
    size_t num_phenotype() const { return (pheno_info.use_pheno)? pheno_info.name.size() : 1; };
    void prsice(const Commander &c_commander, const std::vector<std::string> &region_name, const size_t c_pheno_index,
    		Genotype &target, bool prslice=false);

    //working in progress
    void prsice(const Commander &c_commander, const Region &c_region, const size_t c_pheno_index,  bool prslice=false);
    void output(const Commander &c_commander, const Region &c_region, size_t pheno_index,
    		Genotype &target) const;
    //void output(const Commander &c_commander, size_t pheno_index) const;
    // PRSlice related stuff
    //void prslice_windows(const Commander &c_commander, const Region &c_region);
    //void prslice(const Commander &c_commander, const Region &c_region, const size_t c_pheno_index);
protected:
private:

    inline double lee_adjust(double r2, double top, double bottom) const
    {
        return top*r2/(bottom*r2);
    }
    void set_lee(double prevalence, double case_ratio, double &top, double &bottom) const;
    struct{
            std::vector<int> col;
            std::vector<std::string> name;
            std::vector<int> order;
            std::vector<bool> binary;
            bool use_pheno;
    }pheno_info;

    //slowly update the class
    //input related
    std::string m_base_name;
    std::string m_target;
    std::vector<bool> m_target_binary;
	SCORING m_score = SCORING::MEAN_IMPUTE;
    size_t m_region_size=1;
    bool m_ignore_fid = false;

    misc::vec2d<double> m_region_perm_result;
    std::vector<Sample> m_sample_names;
	std::unordered_map<std::string,size_t> m_sample_with_phenotypes;

	Eigen::VectorXd m_phenotype;
	Eigen::MatrixXd m_independent_variables;

	// Use struct for more elegant coding
	// struct stored in storage so that Genotype class can also use it
	misc::vec2d<prsice_result> m_prs_results;
	//std::vector< std::vector<prsice_result> > m_prs_results; // 1d = region, 2d=results
	std::vector<size_t> m_best_index; // only need to store the index for the best region
	// we are safe to assume that the order of samples follow the read in from fam
	// due to the way we initialize the m_pheno and m_independent_variable
	//misc::vec2d<Sample_lite> m_best_score;
	misc::vec2d<Sample_lite> m_current_sample_score;
	misc::vec2d<Sample_lite> m_best_sample_score;
	std::vector<std::string> m_sample_included;
	std::vector<int> m_sample_index;
	//std::vector< std::vector<Sample_lite> > m_best_score; // PRS for the best threshold
	//std::vector< std::vector<Sample_lite> > m_current_score; // PRS for the current threshold
	std::vector<size_t> m_num_snp_included;
	/**
	 * function area
	 */
	void gen_pheno_vec(const std::string &pheno_file_name, const int pheno_index, bool regress);
	std::vector<size_t> get_cov_index(const std::string &c_cov_file,
			const std::vector<std::string> &c_cov_header);
	void gen_cov_matrix(const std::string &c_cov_file, const std::vector<std::string> &c_cov_header);
	//This should help us to update the m_prs_results
	void process_permutations();
    // valid sample information

	// matrix for regression
	// important guide for all operation
	//std::vector<p_partition> m_partition;
	// snp information

	// Null information
	double m_null_r2 = 0.0;
	// others
	// For thread safety
    static std::mutex score_mutex;

    // Region info
    // PRSlice related storages
    //enum prslice_wind{WIND,SNPS, R2, P, NSNP, COEFF};
    //typedef std::tuple<std::string, std::vector<p_partition>, double, double, double, double > windows;
    //std::vector<windows> m_best_snps;




    /**
     * Check whether if the SNP is included in the target file. This should update
     * the m_include_snp
     * @param c_target_bim_name The target bim file name
     * @param num_ambig Number of ambiguous SNPs
     * @param not_found Number of base SNPs not found in the target file
     * @param num_duplicate Number of duplicated SNPs found in the target file
     */

    //void update_line(std::unordered_map<std::string, size_t> &partition_index);
    void thread_score(size_t region_start, size_t region_end, double threshold,
            size_t thread, const size_t c_pheno_index, const size_t iter_threshold);
    void thread_perm(Eigen::ColPivHouseholderQR<Eigen::MatrixXd> &decomposed,
            std::vector<Eigen::MatrixXd> &pheno_perm, size_t start, size_t end,
            size_t i_region, int rank, const Eigen::VectorXd &pre_se, size_t processed,
            bool logit_perm);
    void permutation(unsigned int seed, int perm_per_slice, int remain_slice,
                int total_permutation, int n_thread, bool logit_perm);
};

#endif // PRSICE_H
