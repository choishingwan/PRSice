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
#include "commander.hpp"
#include "misc.hpp"
#include "plink.hpp"
#include "snp.hpp"
#include "region.hpp"
#include "regression.hpp"

//This should be the class to handle all the procedures
class PRSice
{
    public:
        PRSice();
        virtual ~PRSice();
        void run(const Commander &c_commander, Region &region);
        void process(const std::string &c_input, bool binary, const Commander &c_commander, Region &region);
        // IID PRS
        typedef std::pair<std::string, double> prs_score;
    protected:
    private:
//      This is used for thead safety
        static std::mutex score_mutex;
        // rsid, line number in bim, category, snp_list index
        typedef std::tuple<std::string, size_t, int, size_t> p_partition;
        // threshold, r2, r2 adjust, p, num_snps
        typedef std::tuple<double, double, double, double, size_t> PRSice_result;
        // r2, threshold, num_snps (The r2 is for determining if it is the best)
        typedef std::tuple<double,double,  size_t> PRSice_best;
        std::string m_current_base;
        // This function is responsible for obtaining the PRS from the genotype file
        // quick_ref is used as the classifier to indicate whether if the snp is required
        // snp_list contain the SNP information
        // target is the plink prefix of the target file
        // score is the vector containing the PRS score (before dividing by num_snp)
        // num_snp_included is the number of SNP included so far
        // cur_index indicate which part of quick_ref is at. Will be updated after the
        //           function call
        bool get_prs_score(const std::vector<PRSice::p_partition> &quick_ref,
        			const boost::ptr_vector<SNP> &snp_list, const std::string &target,
				std::vector< std::vector<std::pair<std::string, double> > > &score,
				std::vector<size_t> &num_snp_included, size_t &cur_index);
        // If an independent plink file is provided for LD calculation, we would like to
        // make sure only SNPs also observed in the target file are used for clumping
        // inclusion is the map of SNPs that should be included and their corresponding
        //           index in the snp_list vector
        // c_target_bim_name is the name of the target file
        // snp_list contain information of SNPs
        // c_snp_index contain the index for each SNP
        void update_inclusion(std::map<std::string, size_t> &inclusion, const std::string &c_target_bim_name,
        			boost::ptr_vector<SNP> &snp_list, const std::map<std::string, size_t> &c_snp_index);
        void get_snp(boost::ptr_vector<SNP> &snp_list, std::map<std::string, size_t> &snp_index,
        			const std::string &c_input, bool beta, const Commander &c_commander, Region &region);
        Eigen::MatrixXd gen_cov_matrix(const std::string &c_target, const std::string &c_cov_file,
        			const std::vector<std::string> &c_cov_header, std::map<std::string, size_t> &fam_index);
        Eigen::VectorXd gen_pheno_vec(const std::string &c_target, const std::string c_pheno,
        			const std::vector<std::string> pheno_col, bool target_binary,
					std::map<std::string, size_t> &fam_index, std::vector<std::pair<std::string, double> > &prs_score);
        void calculate_score(const Commander &c_commander, const std::map<std::string, size_t> &inclusion,
        		const boost::ptr_vector<SNP> &snp_list, const Region &c_region);
        void target_process(const Commander &c_commander, const std::map<std::string, size_t> &inclusion,
        		const boost::ptr_vector<SNP> &snp_list, const Region &c_region);
        void thread_score( Eigen::MatrixXd &independent_variables, const Eigen::VectorXd &c_pheno,
        		const std::vector<std::vector<PRSice::prs_score > > &c_prs_region_score,
        		const std::vector<size_t> &c_num_snp_included, const std::map<std::string, size_t> &c_fam_index,
        		std::vector<PRSice::PRSice_best > &prs_best_info,
        		std::vector<std::vector<PRSice::prs_score> > & prs_best_score,
			std::vector<std::vector<PRSice::PRSice_result> > &region_result,
        		size_t region_start, size_t region_end, bool target_binary, double threshold, size_t thread);
};

#endif // PRSICE_H
