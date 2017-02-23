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
#include "commander.hpp"
#include "misc.hpp"
#include "plink.hpp"
#include "snp.hpp"
#include "region.hpp"
#include "regression.hpp"
#include "storage.hpp"
//This should be the class to handle all the procedures
class PRSice
{
public:
    PRSice(std::string base_name, int index, std::string target, std::vector<bool> target_binary,
    		size_t permutation, SCORING score, size_t num_region, bool ignore_fid):
    			m_base_name(base_name), m_base_index(index), m_target(target),
				m_target_binary(target_binary), m_perm(permutation), m_score(score),
				m_region_size(num_region), m_ignore_fid(ignore_fid){
        if(index < 0)
        {
            throw std::out_of_range("Index cannot be less than 0");
        }
    };
    virtual ~PRSice();
    void get_snp(const Commander &c_commander, Region &region);

    // Under construction
    void perform_clump(const Commander &c_commander);
    void pheno_check(const Commander &c_commander);

    void init_matrix(const Commander &c_commander, const size_t c_pheno_index, const bool prslice);
    size_t num_phenotype() const { return m_pheno_names.size(); };
    void categorize(const Commander &c_commander);
    void prsice(const Commander &c_commander, const Region &c_region, const size_t c_pheno_index,  bool prslice=false);
    void output(const Commander &c_commander, const Region &c_region, size_t pheno_index) const;
    void output(const Commander &c_commander, size_t pheno_index) const;
    // PRSlice related stuff
    void prslice_windows(const Commander &c_commander, const Region &c_region);
    void prslice(const Commander &c_commander, const Region &c_region, const size_t c_pheno_index);
protected:
private:
    //slowly update the class
    //input related
    std::string m_base_name;
    int m_base_index;
    std::string m_target;
    std::vector<bool> m_target_binary;
	size_t m_perm = 0;
	SCORING m_score = SCORING::MEAN_IMPUTE;
    size_t m_region_size=1;
    bool m_ignore_fid = false;

    // valid sample information
	std::unordered_map<std::string,size_t> m_sample_with_phenotypes;
	std::vector<prs_score> m_sample_names;
	// matrix for regression
	Eigen::VectorXd m_phenotype;
	Eigen::MatrixXd m_independent_variables;
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
    std::vector<PRSice_best> m_best_threshold;
    std::vector<std::vector<prs_score> > m_best_score;
    std::vector<std::vector<prs_score> > m_current_prs;
    std::vector<std::vector<PRSice_result> > m_prs_results;
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
    void check_inclusion(const std::string &c_target_bim_name,
    		size_t &num_ambig, size_t &not_found, size_t &num_duplicate);

    void update_line(std::unordered_map<std::string, size_t> &partition_index);
    void gen_pheno_vec(const std::string c_pheno, const int pheno_index, const int col_index, bool regress);
    void gen_cov_matrix(const std::string &c_cov_file, const std::vector<std::string> &c_cov_header);
    bool get_prs_score(size_t &cur_index, PLINK &score_plink);
    void thread_score(size_t region_start, size_t region_end, double threshold, size_t thread, const size_t c_pheno_index);

};

#endif // PRSICE_H
