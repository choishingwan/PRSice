#ifndef PRSICE_H
#define PRSICE_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <Eigen/Dense>
#include <boost/ptr_container/ptr_vector.hpp>
#include <vector>
#include <algorithm>
#include <map>
#include <stdio.h>
#include <thread>
#include <unordered_map>
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
    PRSice(std::string base_name, int index, std::string target, std::vector<bool> target_binary, size_t permutation): m_base_name(base_name), m_base_index(index),
    		m_target(target), m_target_binary(target_binary), m_perm(permutation)
    {
        if(index < 0)
        {
            throw std::out_of_range("Index cannot be less than 0");
        }
    };
    virtual ~PRSice();
    void get_snp(const Commander &c_commander, Region &region);



    /**
     * This function will perform clumping based on SNPs found in both the target and LD files
     * All required parameters are obtained from c_commander. Index SNPs will be stored in
     * m_include_snp. As different base file might contain different SNPs, this has to be
     * done separately for each base files
     * @param c_commander Contain all the parameters e.g. clumping threshold.
     */
    void clump(const Commander &c_commander);
    /**
     * This function will compute the required phenotype that will be used for PRS.
     * Result phenotype will be stored in m_pheno_names
     * @param c_commander Again, c_commander serves as the container of all parameters
     */
    void init_pheno(const Commander &c_commander);
    /**
     * This function compute the required matrix and the null r2
     * @param c_commander The paramter container
     * @param c_pheno_index Current phenotype
     * @param prslice Whether if this is for PRSlice. For PRSlice, we will not allow output of all
     */
    void init_matrix(const Commander &c_commander, const size_t c_pheno_index, const bool prslice);
    /**
     * This function return the number of valid phenotype for the analysis
     * @return Number of valid phenotype for the analysis
     */
    size_t num_phenotype() const { return m_pheno_names.size(); };
    /**
     * Essentially calculate the index for PRSice run. Categorize SNPs based on
     * their p-values. Will update the m_partition vector
     * @param c_commander Parameter container
     */
    void categorize(const Commander &c_commander);
    void prsice(const Commander &c_commander, const Region &c_region, const size_t c_pheno_index,  bool prslice=false);
    /**
     * Output the results
     * @param c_commander List of parameters
     * @param c_region List of regions
     * @param pheno_index Indication of which phenotype we are working with
     */
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


    // Phenotype storages
    enum pheno_store{FILE_NAME, INDEX, NAME, ORDER};
    typedef std::tuple<std::string, size_t, std::string, size_t> pheno_storage;
	std::vector<pheno_storage> m_pheno_names;
	size_t m_pheno_index=0;
	size_t m_perm = 0;
	// Regression related storages
	double m_null_r2 = 0.0;
	Eigen::VectorXd m_phenotype;
	Eigen::MatrixXd m_independent_variables;
	std::unordered_map<std::string,size_t> m_sample_with_phenotypes;
	// For thread safety
    static std::mutex score_mutex;

    // Holder vector containing the sample names in the target file
	std::vector<prs_score> m_sample_names;
	// PRS storages
    std::vector<p_partition> m_partition; //
    std::vector<size_t> m_num_snp_included; //
    std::vector<PRSice_best> m_best_threshold; //
    std::vector<std::vector<prs_score> > m_best_score; //
    std::vector<std::vector<prs_score> > m_current_prs; //
    std::vector<std::vector<PRSice_result> > m_prs_results; //
    // SNPs found in the base file
    boost::ptr_vector<SNP> m_snp_list; //
    std::unordered_map<std::string, size_t> m_snp_index; //
    // PRSlice related storages
    enum prslice_wind{WIND,SNPS, R2, P, NSNP, COEFF};
    typedef std::tuple<std::string, std::vector<p_partition>, double, double, double, double > windows;
    std::vector<windows> m_best_snps;
    /**
     * This vector contain the chromosome included in the base file
     * Should be used to guide the plink class to read the multi-chromosome
     * files
     */
    std::vector<std::string> m_chr_list;
    std::unordered_map<std::string, size_t> m_include_snp; //


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
    bool get_prs_score(size_t &cur_index);
    void thread_score(size_t region_start, size_t region_end, double threshold, size_t thread, const size_t c_pheno_index);

};

#endif // PRSICE_H
