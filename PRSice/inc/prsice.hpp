#ifndef PRSICE_H
#define PRSICE_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <Eigen/Dense>
#include <vector>
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
        bool score(const std::map<std::string, size_t> &inclusion, std::vector<SNP> &snp_list,
        		const std::string &target, std::vector<double> &score,
			double threshold_lower, double threshold_upper, size_t &num_snp);
        void update_inclusion(std::map<std::string, size_t> &inclusion, const std::string &target_bim_name,
                       std::vector<SNP> &snp_list, const std::map<std::string, size_t> &snp_index);
        Eigen::MatrixXd gen_cov_matrix(const std::string &c_target, const std::string &c_cov_file, const std::vector<std::string> &c_cov_header, std::map<std::string, int> &pheno_missing, size_t num_sample);
        Eigen::VectorXd gen_pheno_vec(const std::string &c_target, const std::string &c_pheno, bool target_binary, double &num_sample, std::map<std::string, int> &pheno_missing, std::vector<bool> &pheno_index);
};

#endif // PRSICE_H
