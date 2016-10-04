#ifndef COMMANDER_H
#define COMMANDER_H

#include <string>
#include <getopt.h>
#include <unistd.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "misc.hpp"

class Commander
{
    public:
        Commander();
        bool initialize(int argc, char *argv[]);
        virtual ~Commander();
        std::vector<std::string> get_base() const{ return m_base; };
        bool index() const{ return m_index;};
        bool gen_bed() const{ return m_gen_bed; };
        std::string chr() const{ return m_chr; };
        std::string ref() const{ return m_ref_allele; };
        std::string alt() const{ return m_alt_allele; };
        std::string statistic() const{ return m_statistic; };
        std::string snp() const{ return m_snp; };
        std::string bp() const{ return m_bp; };
        std::string se() const{ return m_standard_error; };
        std::string p() const{ return m_p_value; };
        std::string ld_prefix() const{ return m_ld_prefix;};
        std::string get_gtf() const { return m_gtf; };
        std::vector<std::string> get_target() const{ return m_target; };
        std::vector<std::string> get_bed() const { return m_bed_list; };
        std::vector<std::string> get_cov_header() const { return m_covariates; };
        std::vector<std::string> get_pheno() const { return m_pheno_file; };
        std::string get_pheno(size_t i) const { return (m_pheno_file.size()==0)? "": m_pheno_file.at(i); };
        std::string get_target(size_t i) const{ return m_target.at(i); };
        std::string get_msigdb() const{ return m_msigdb; };
        std::string get_out() const { return m_out; };
        std::string get_cov_file() const { return m_covariate_file; };
        size_t get_thread() const { return m_thread; };
        size_t get_clump_kb() const { return m_clump_kb; };
        double get_clump_p() const { return m_clump; };
        double get_clump_r2() const { return m_clump_r2; };
        double get_lower() const { return m_lower; };
        double get_upper() const { return m_upper; };
        double get_inter() const { return m_inter; };
        double get_proxy() const { return m_proxy; };
        int get_category(double p) const {
        		for(int i = 0; i < m_barlevel.size(); ++i){
        			if(p < m_barlevel[i]){
        				return i-1;
        			}
        		}
        		return -2;
        };
        bool get_target_binary(size_t i) const { return m_target_is_binary.at(i); };
        bool get_base_binary(size_t i) const { return m_use_beta.at(i); };
        bool no_regression() const { return m_no_regress; };
        bool fastscore() const { return m_fastscore; };
//        bool proxy() const { return m_proxy; };
    protected:
    private:
        void usage();
        std::vector<std::string> m_target;
        std::vector<bool> m_target_is_binary;
        std::vector<bool> m_use_beta;
        std::vector<std::string> m_base;
        std::vector<std::string> m_covariates; //Should be mutrally exclusive with m_covariate_files
        std::vector<std::string> m_bed_list;
        std::vector<std::string> m_pheno_file;
        std::vector<double> m_barlevel;
        std::string m_covariate_file;
        std::string m_ancestry_dim;
        std::string m_chr = "CHR";
        std::string m_ref_allele="A1";
        std::string m_alt_allele="A2";
        std::string m_statistic = "OR";
        std::string m_snp="SNP";
        std::string m_bp="BP";
        std::string m_standard_error = "SE";
        std::string m_p_value = "P";
        std::string m_ld_prefix=""; //Allow people to use external reference for LD estimation
        std::string m_gtf="";
        std::string m_msigdb="";
        std::string m_out = "PRSice";
        bool m_clumping=true;
        bool m_pruning = false;
        bool m_binary_target=true;
        bool m_covary=true;
        bool m_multiple_base_phenotypes=false;
        bool m_remove_mhc=false;
        bool m_fastscore =false;
        bool m_index =false; //Indicate that instead of a header name, the index is given
        bool m_gen_bed = false;
        bool m_no_regress = false;
//        bool m_proxy = false;
        double m_proxy = -1;
        double m_clump = 1;
        //double m_clump_p2 = 1;
        double m_clump_r2 = 0.1;
        size_t m_clump_kb = 250000;
        //We disable prunning

        double m_lower = 0.0001;
        double m_upper = 0.5;
        double m_inter = 0.00005;

        size_t m_thread=1;
};

#endif // COMMANDER_H
