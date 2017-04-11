#ifndef PRSICE_H
#define PRSICE_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <Eigen/Dense>
#include <mutex>
#include <boost/ptr_container/ptr_vector.hpp>
#include <vector>
#include <algorithm>
#include <map>
#include <stdio.h>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include "genotype.hpp"
#include "commander.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "snp.hpp"
#include "region.hpp"
#include "regression.hpp"
#include "storage.hpp"
//This should be the class to handle all the procedures
class PRSice
{
public:

    PRSice(std::string base_name, std::string target, std::vector<bool> target_binary,
    		size_t permutation, SCORING score, size_t num_region, bool ignore_fid):
    			m_base_name(base_name), m_target(target), m_target_binary(target_binary),
				m_perm(permutation), m_score(score), m_region_size(num_region),
				m_ignore_fid(ignore_fid){ };
    virtual ~PRSice();
    void pheno_check(const Commander &c_commander);
    void init_matrix(const Commander &c_commander, const size_t pheno_index, Genotype &target,
    		const bool prslice=false);
    size_t num_phenotype() const { return m_pheno_names.size(); };
    void prsice(const Commander &c_commander, const std::vector<std::string> &region_name, const size_t c_pheno_index,
    		Genotype &target, bool prslice=false);

    //working in progress
    void prsice(const Commander &c_commander, const Region &c_region, const size_t c_pheno_index,  bool prslice=false);
    void output(const Commander &c_commander, const Region &c_region, size_t pheno_index) const;
    void output(const Commander &c_commander, size_t pheno_index) const;
    // PRSlice related stuff
    void prslice_windows(const Commander &c_commander, const Region &c_region);
    void prslice(const Commander &c_commander, const Region &c_region, const size_t c_pheno_index);
protected:
private:
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
	size_t m_perm = 0;
	SCORING m_score = SCORING::MEAN_IMPUTE;
    size_t m_region_size=1;
    bool m_ignore_fid = false;

    std::vector<Sample> m_sample_names;
	std::unordered_map<std::string,size_t> m_sample_with_phenotypes;

	Eigen::VectorXd m_phenotype;
	Eigen::MatrixXd m_independent_variables;

	// Use struct for more elegant coding
	// struct stored in storage so that Genotype class can also use it
	std::vector< std::vector<prsice_result> > m_prs_results; // 1d = region, 2d=results
	std::vector<size_t> m_best_index; // only need to store the index for the best region
	// we are safe to assume that the order of samples follow the read in from fam
	// due to the way we initialize the m_pheno and m_independent_variable
	Eigen::MatrixXd m_best_score; // PRS for the best threshold
	Eigen::MatrixXd m_current_score; // PRS for the current threshold

	/**
	 * function area
	 */
	void gen_pheno_vec(const std::string &pheno_file_name, const int pheno_index, bool regress);
	std::vector<size_t> get_cov_index(const std::string &c_cov_file,
			const std::vector<std::string> &c_cov_header);
	void gen_cov_matrix(const std::string &c_cov_file, const std::vector<std::string> &c_cov_header);

    // valid sample information

	// matrix for regression
	// important guide for all operation
	std::vector<p_partition> m_partition;
	// snp information
    boost::ptr_vector<SNP> m_snp_list;
    std::unordered_map<std::string, size_t> m_snp_index; // only use for reading in information
    std::vector<std::string> m_chr_list; // chromosome information
    std::unordered_map<std::string, size_t> m_include_snp; // the information provider until categorize

    // Phenotype storages
    enum pheno_store{FILE_NAME, INDEX, NAME, ORDER};
    typedef std::tuple<std::string, size_t, std::string, size_t> pheno_storage;
	std::vector<pheno_storage> m_pheno_names;
	// Null information
	double m_null_r2 = 0.0;
	// others
	// For thread safety
    static std::mutex score_mutex;

    // Holder vector containing the sample names in the target file
    std::vector<size_t> m_num_snp_included;
    std::vector<std::vector<prs_score> > m_current_prs;
    // Region info
    // PRSlice related storages
    enum prslice_wind{WIND,SNPS, R2, P, NSNP, COEFF};
    typedef std::tuple<std::string, std::vector<p_partition>, double, double, double, double > windows;
    std::vector<windows> m_best_snps;




    /**
     * Check whether if the SNP is included in the target file. This should update
     * the m_include_snp
     * @param c_target_bim_name The target bim file name
     * @param num_ambig Number of ambiguous SNPs
     * @param not_found Number of base SNPs not found in the target file
     * @param num_duplicate Number of duplicated SNPs found in the target file
     */

    void update_line(std::unordered_map<std::string, size_t> &partition_index);
    bool get_prs_score(size_t &cur_index, PLINK &score_plink);
    void thread_score(size_t region_start, size_t region_end, double threshold, size_t thread, const size_t c_pheno_index);

};

#endif // PRSICE_H
