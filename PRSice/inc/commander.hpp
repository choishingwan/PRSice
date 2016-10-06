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
        std::vector<std::string> get_bed() const { return m_bed_list; };
        std::vector<std::string> get_cov_header() const { return m_covariates; };
        std::vector<std::string> get_pheno_col() const { return m_pheno_col; };
        std::string get_pheno() const { return m_pheno_file; };
        std::string get_target() const{ return m_target; };
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
        bool target_is_binary() const { return m_target_is_binary; };
        bool get_base_binary(size_t i) const { return m_use_beta.at(i); };
        bool no_regression() const { return m_no_regress; };
        bool fastscore() const { return m_fastscore; };
        bool full() const { return m_full; };
        bool all() const { return m_all; };
//        bool proxy() const { return m_proxy; };
    protected:
    private:
        void usage();
        std::vector<bool> m_use_beta;
        std::vector<std::string> m_base;
        std::vector<std::string> m_covariates; //Should be mutrally exclusive with m_covariate_files
        std::vector<std::string> m_bed_list;
        std::vector<std::string> m_pheno_col;
        std::vector<double> m_barlevel;
        std::string m_target;
        std::string m_pheno_file;
        std::string m_covariate_file;
        std::string m_ancestry_dim;
        std::string m_chr;
        std::string m_ref_allele;
        std::string m_alt_allele;
        std::string m_statistic ;
        std::string m_snp;
        std::string m_bp;
        std::string m_standard_error;
        std::string m_p_value;
        std::string m_ld_prefix; //Allow people to use external reference for LD estimation
        std::string m_gtf;
        std::string m_msigdb;
        std::string m_out;
        bool m_target_is_binary;
        bool m_fastscore;
        bool m_index;//Indicate that instead of a header name, the index is given
        bool m_gen_bed;
        bool m_no_regress;
        bool m_all;
        bool m_full;
        double m_proxy;
        double m_clump;
        double m_clump_r2;
        size_t m_clump_kb;
        double m_lower;
        double m_upper;
        double m_inter;
        size_t m_thread;
};

#endif // COMMANDER_H
