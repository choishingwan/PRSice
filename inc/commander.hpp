#ifndef COMMANDER_H
#define COMMANDER_H

#include "storage.hpp"
#include <string>
#include <getopt.h>
#include <unistd.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "misc.hpp"

class Commander
{
public:
	// so that we don't need to include plink_common.hpp here
    Commander();
    virtual ~Commander();
    bool initialize(int argc, char *argv[]);

    //base
    std::vector<int> index() const { return base.col_index; };
    bool has_index() const { return base.index; };
    bool beta() const { return base.beta; };
    std::string base_name() const { return base.name; };

    //clump
    std::string ld_prefix() const { return clumping.ld; };
    std::string ld_type() const { return clumping.type; };
    bool no_clump() const { return clumping.no_clump; };
    bool use_proxy() const { return clumping.provide_proxy; };
    double proxy() const { return clumping.proxy; };
    double clump_p() const { return clumping.p_value; };
    double clump_r2() const { return clumping.r2; };
    double clump_dist() const { return clumping.distance; };

    //covariate
    std::string get_cov_file() const { return covariate.name; };
    std::vector<std::string> get_cov_header() const { return covariate.covariates; };

    //misc
    bool all() const { return misc.all; };
    bool print_snp() const { return misc.print_snp; };
    bool ignore_fid() const { return misc.ignore_fid; };
    std::string out() const { return misc.out; };
    int permutation() const { return misc.permutation; };
    int thread() const { return misc.thread; };

    // prset
    std::vector<std::string> bed() const { return prset.bed; };
    std::vector<std::string> feature() const { return prset.feature; };
    std::string gtf() const { return prset.gtf; };
    std::string msigdb() const { return prset.msigdb; };

    // prsice
    std::string missing_score() const { return prsice.missing_score; };
    SCORING get_scoring() const{
    		if(prsice.missing_score == "no_mean_imputation") return SCORING::SET_ZERO;
    		else if(prsice.missing_score == "center") return SCORING::CENTER;
    		else return SCORING::MEAN_IMPUTE;
    }

    double lower() const { return prsice.lower; };
    double upper() const { return prsice.upper; };
    double inter() const { return prsice.inter; };
    bool no_regress() const { return prsice.no_regress; };
    bool full() const { return prsice.full; };
    bool fastscore() const { return prsice.fastscore; };
    double bar_upper() const { return prsice.barlevel.front(); };// we have sorted it
    double bar_lower() const { return prsice.barlevel.back(); };

    //prslice
    bool perform_prslice() const { return prslice.provided; };
    int prslice_size() const { return prslice.size; };

    //species
    int num_auto() const { return species.num_auto; };
    bool no_x() const { return species.no_x; };
    bool no_y() const { return species.no_y; };
    bool no_xy() const { return species.no_xy; };
    bool no_mt() const { return species.no_mt; };

    //target
    std::string target_name() const { return target.name; };
    std::string target_type() const { return target.type; };
    std::string pheno_file() const { return target.pheno_file; };
    std::string pheno_col(size_t index) const { return target.pheno_col.at(index); };
    std::vector<std::string> pheno_col() const { return target.pheno_col; };
    std::vector<bool> is_binary() const { return target.is_binary; };
    bool is_binary(size_t index) const { return target.is_binary.at(index); };

    int get_category(double p) const
    {
        for(int i = 0; i < prsice.barlevel.size(); ++i)
        {
            if(p <= prsice.barlevel[i])
            {
                return i;
            }
        }
        if(p > prsice.barlevel.back()) return prsice.barlevel.size();
        return -2;
    };
    double get_threshold(int i) const
    {
    	return(i < 0)? prsice.barlevel.at(0) :
    			((i>=prsice.barlevel.size())? 1.0 :prsice.barlevel.at(i));
    };

    void user_input() const;
protected:
private:
    std::string version ="2.0beta";
    std::string date = "29 March 2017";
    std::vector<std::string> supported_types = {"bed", "ped", "bgen"};
    struct{
    	std::string name;
        std::string chr;
        std::string ref_allele;
        std::string alt_allele;
        std::string statistic;
        std::string snp;
        std::string bp;
        std::string standard_error;
        std::string p_value;
        std::vector<int> col_index;
    	int beta;
        int index;
        bool provided_chr;
        bool provided_ref;
        bool provided_alt;
        bool provided_stat;
        bool provided_snp;
        bool provided_bp;
        bool provided_se;
        bool provided_p;
        bool provided_beta;
    } base;

    struct{
        std::string ld;
        std::string type;
        int no_clump;
        double proxy;
        double p_value;
        double r2;
        int distance;
        bool provide_proxy;
    } clumping;

    struct{
        std::string name;
        std::string ancestry_dim;
        std::vector<std::string> covariates;
    } covariate;

    struct{
    		double maf;
    		double mind;
    		double geno;
    		double info_score;
    		int use_maf;
    		int use_mind;
    		int use_geno;
    		int use_info;
    } filter;

    struct{
        int all;
        std::string out;
        int print_snp;
        int ignore_fid;
        int permutation;
        int thread;
    } misc;

    struct{
        std::vector<std::string> bed;
        std::vector<std::string> feature;
        std::string gtf;
        std::string msigdb;
    } prset;

    struct{
        std::string missing_score;
        std::vector<double> barlevel;
        double lower;
        double upper;
        double inter;
        int fastscore;
        int no_regress;
        int full;
    } prsice;

    struct{
        int size;
        bool provided;
    } prslice;

    struct{
        int num_auto;
        int no_x;
        int no_y;
        int no_xy;
        int no_mt;
        bool double_set;
    } species;

    struct{
    	std::string name;
        std::string pheno_file;
        std::string type;
        std::vector<std::string> pheno_col;
    	std::vector<bool> is_binary;
    } target;


    void usage();
    void info();
    void program_info();
    void base_check(bool &error, std::string &error_message);
    void clump_check(bool &error, std::string &error_message);
    void misc_check(bool &error, std::string &error_message);
    void prset_check(bool &error, std::string &error_message);
    void prslice_check(bool &error, std::string &error_message);
    void target_check(bool &error, std::string &error_message);
    void prsice_check(bool &error, std::string &error_message);

    inline void set_species(int num_auto, bool no_x, bool no_y, bool no_xy, bool no_mt, bool &error,
    		std::string &error_message, bool &species_error)
    {
    	if(species.double_set && !species_error){
    		species_error = true;
    		error = true;
    		error_message.append("ERROR: Can only specify one species\n");
    	}
    	species.num_auto = num_auto;
    	species.no_x = no_x;
    	species.no_y = no_y;
    	species.no_xy = no_xy;
    	species.no_mt = no_mt;
    	species.double_set = true;
    };

    inline void set_stat(std::string stat)
    {
    	base.statistic = stat;
    	base.provided_stat = true;
    };

    inline void set_beta(bool beta)
    {
    	base.beta = beta;
    	base.provided_beta = true;
    };

    inline void set_chr(std::string chr)
    {
    	base.chr = chr;
    	base.provided_chr = true;
    };

    inline void set_ref(std::string ref)
    {
    	base.ref_allele = ref;
    	base.provided_ref = true;
    };

    inline void set_alt(std::string alt)
    {
    	base.alt_allele = alt;
    	base.provided_alt = true;
    };

    inline void set_snp(std::string snp)
    {
    	base.snp = snp;
    	base.provided_snp = true;
    };

    inline void set_bp(std::string bp)
    {
    	base.bp = bp;
    	base.provided_bp = true;
    };

    inline void set_se(std::string se)
    {
    	base.standard_error = se;
    	base.provided_se = true;
    };

    inline void set_p(std::string p)
    {
    	base.p_value = p;
    	base.provided_p = true;
    };

    inline void set_prslice(int size)
    {
    	prslice.size = size;
    	prslice.provided = true;
    }

    inline int index_check(const std::string &target, const std::vector<std::string> ref) const
    {
    	for(size_t i = 0; i < ref.size(); ++i)
    	{
    		if(target.compare(ref[i])==0)
    		{
    			return i;
    		}
    	}
    	return -1;
    };

    inline int index_check(const std::string &target, const int max, bool &error,
    		std::string &error_message, std::string name)
    {
    	try{
    		int index = misc::convert<int>(optarg);
    		if(index >= max){
    			error =true;
    			error_message.append("ERROR: "+name+" index out of bound!\n");
    			return -1;
    		}
    		if(index < 0){
    			error =true;
    			error_message.append("ERROR: Negative "+name+" index!\n");
    			return -1;
    		}
    		return index;
    	}
    	catch(const std::runtime_error &er)
    	{
			error =true;
			error_message.append("ERROR: "+name+" index is not numeric!\n");
			return -1;
    	}
    }
};

#endif // COMMANDER_H
