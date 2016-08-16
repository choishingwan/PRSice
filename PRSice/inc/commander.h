#ifndef COMMANDER_H
#define COMMANDER_H

#include <string>
#include <getopt.h>
#include <unistd.h>

class Commander
{
    public:
        Commander();
        Commander(const int argc, const char *argv[]);
        virtual ~Commander();
        std::string get_target() const{ return m_target; };
    protected:
    private:
        void usage();
        std::vector<std::string> m_target;
        std::vector<std::string> m_base;
        std::vector<std::string> m_covariates; //Should be mutrally exclusive with m_covariate_files
        std::vector<std::string> m_covariate_files;
        std::string m_ancestry_dim;
        std::string m_pheno_file;
        std::string m_chr = "CHR";
        std::string m_ref_allele="A1;
        std::string m_alt_allele="A2";
        std::string m_statistic = "OR";
        std::string m_snp="SNP";
        std::string m_bp="BP";
        std::string m_standard_error = "SE";
        std::string m_p_value = "P";

        bool m_clumping=true;
        bool m_pruning = false;
        bool m_binary_target=true;
        bool m_covary=true;
        bool m_multiple_base_phenotypes=false;
        bool m_remove_mhc=false;
        bool m_use_beta=false;
        bool m_fastscore =false;

        double m_clump_p1 = 1;
        double m_clump_p2 = 1;
        double m_clump_r2 = 0.1;
        double m_clump_kb = 250;

        //Might want to disable the pruning part
        double m_prune_kb_wind = 50;
        double m_prune_kb_step = 2;
        double m_prune_kb_r2 = 0.8;

        double m_slower = 0.0001;
        double m_supper = 0.5;
        double m_sinc = 0.00005;

        size_t m_thread=1;
};

#endif // COMMANDER_H
