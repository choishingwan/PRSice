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
#include "commander.hpp"
#include "misc.hpp"
#include "plink.hpp"
#include "snp.hpp"
#include "regression.h"
#include "region.hpp"

//This should be the class to handle all the procedures
class PRSice
{
    public:
        PRSice();
        virtual ~PRSice();
        void run(const Commander &c_commander, Region &region);
        void process(const std::string &c_input, bool binary, const Commander &c_commander, Region &region);
    protected:
    private:
        typedef std::tuple<std::string, size_t, size_t, size_t> p_partition;
        bool get_prs_score(const std::vector<PRSice::p_partition> &quick_ref,
        			const boost::ptr_vector<SNP> &snp_list, const std::string &target,
				std::vector<p_partition> &score, size_t &num_snp_included, size_t &cur_index);
        void update_inclusion(std::map<std::string, size_t> &inclusion, const std::string &c_target_bim_name,
        			boost::ptr_vector<SNP> &snp_list, const std::map<std::string, size_t> &c_snp_index);
        void get_snp(boost::ptr_vector<SNP> &snp_list, std::map<std::string, size_t> &snp_index,
        			const std::string &c_input, bool beta, const Commander &c_commander, Region &region);
        Eigen::MatrixXd gen_cov_matrix(const std::string &c_target, const std::string &c_cov_file,
        			const std::vector<std::string> &c_cov_header, std::map<std::string, size_t> &fam_index);
        Eigen::VectorXd gen_pheno_vec(const std::string &c_target, const std::string c_pheno,
        			bool target_binary, std::map<std::string, size_t> &fam_index);
        void calculate_score(const Commander &c_commander, bool target_binary,
        			const std::string c_target, const std::map<std::string, size_t> &inclusion,
				const boost::ptr_vector<SNP> &snp_list);
};

#endif // PRSICE_H
