#ifndef COMMANDER_H
#define COMMANDER_H

#include <string>
#include <getopt.h>
#include <unistd.h>
#include <vector>
#include <stdexcept>
#include "misc.hpp"

class Commander
{
    public:
        Commander();
        bool initialize(int argc, char *argv[]);
        virtual ~Commander();
        std::vector<std::string> get_base() const{ return m_base;};
        bool index() const{ return m_index;};
        std::string chr() const{ return m_chr;};
        std::string ref() const{ return m_ref_allele;};
        std::string alt() const{ return m_alt_allele;};
        std::string statistic() const{ return m_statistic;};
        std::string snp() const{ return m_snp;};
        std::string bp() const{ return m_bp;};
        std::string se() const{ return m_standard_error;};
        std::string p() const{ return m_p_value;};
        std::string ld_prefix() const{return m_ld_prefix;};
        std::vector<std::string> get_target() const{return m_target; };
    protected:
    private:
        void usage();
        std::vector<std::string> m_target;
        std::vector<bool> m_target_is_binary;
        std::vector<std::string> m_base;
        std::vector<std::string> m_covariates; //Should be mutrally exclusive with m_covariate_files
        std::vector<std::string> m_covariate_files;
        std::string m_ancestry_dim;
        std::string m_pheno_file;
        std::string m_chr = "CHR";
        std::string m_ref_allele="A1";
        std::string m_alt_allele="A2";
        std::string m_statistic = "OR";
        std::string m_snp="SNP";
        std::string m_bp="BP";
        std::string m_standard_error = "SE";
        std::string m_p_value = "P";
        std::string m_ld_prefix=""; //Allow people to use external reference for LD estimation
    
        bool m_clumping=true;
        bool m_pruning = false;
        bool m_binary_target=true;
        bool m_covary=true;
        bool m_multiple_base_phenotypes=false;
        bool m_remove_mhc=false;
        bool m_use_beta=false;
        bool m_fastscore =false;
        bool m_index =false; //Indicate that instead of a header name, the index is given

        double m_clump_p1 = 1;
        double m_clump_p2 = 1;
        double m_clump_r2 = 0.1;
        size_t m_clump_kb = 250;
        //We disable prunning

        double m_lower = 0.0001;
        double m_upper = 0.5;
        double m_inter = 0.00005;

        size_t m_thread=1;
};

#endif // COMMANDER_H
