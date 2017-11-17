// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef COMMANDER_H
#define COMMANDER_H

#include "misc.hpp"
#include "storage.hpp"
#include <chrono>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <random>
#include <reporter.hpp>
#include <stdexcept>
#include <string>
#include <thread>
#include <unistd.h>
#include <unordered_set>
#include <vector>
const std::string version = "2.0.15.beta";
const std::string date = "27 October 2017";
class Commander
{
public:
    Commander();
    virtual ~Commander();
    bool init(int argc, char* argv[], Reporter& reporter);

    // base
    std::vector<int> index() const { return base.col_index; };
    bool has_index() const { return base.index; };
    bool beta() const { return base.beta; };
    bool filter_base_info() const { return base.use_info; };
    bool filter_base_maf() const { return base.provided_maf; };
    double base_info_score() const { return base.info_score; };
    double maf_base() const { return base.maf_threshold; };
    std::string base_name() const { return base.name; };

    // clump
    std::string ld_prefix() const { return clumping.ld; };
    std::string ld_type() const { return clumping.type; };
    bool no_clump() const { return clumping.no_clump; };
    bool use_proxy() const { return clumping.provide_proxy; };
    double proxy() const { return clumping.proxy; };
    double clump_p() const { return clumping.p_value; };
    double clump_r2() const { return clumping.r2; };
    double clump_dist() const { return clumping.distance * 1000; };

    // covariate
    std::string get_cov_file() const { return covariate.name; };
    std::vector<std::string> get_cov_header() const
    {
        return covariate.covariates;
    };

    // filtering
    bool filter_maf() const { return filter.use_maf; };
    bool filter_geno() const { return filter.use_geno; };
    bool filter_mind() const { return filter.use_mind; };
    bool filter_info() const { return filter.info_filtering; };
    bool filter_ld_maf() const { return filter.use_ld_maf; };
    bool filter_ld_geno() const { return filter.use_ld_geno; };
    bool filter_ld_info() const { return filter.use_ld_info; };
    bool filter_hard_threshold() const { return filter.use_hard_thres; };
    bool hard_coding() const { return filter.hard_coding; };
    bool keep_ambig() const { return filter.keep_ambig; };
    std::string extract_snp_file() const
    {
        return filter.extract ? filter.extract_file : "";
    };
    std::string exclude_snp_file() const
    {
        return filter.exclude ? filter.exclude_file : "";
    };
    double maf() const { return filter.maf; };
    double geno() const { return filter.geno; };
    double mind() const { return filter.mind; };
    double info() const { return filter.info_score; };
    double ld_maf() const { return filter.ld_maf; };
    double ld_geno() const { return filter.ld_geno; };
    double ld_info() const { return filter.ld_info; };
    double hard_threshold() const { return filter.hard_threshold; };

    // misc
    bool all() const { return misc.all; };
    bool ignore_fid() const { return misc.ignore_fid; };
    bool logit_perm() const { return misc.logit_perm; };
    bool permute() const { return misc.provided_permutation; };
    bool print_snp() const { return misc.print_snp; };
    bool seeded() const { return misc.provided_seed; };
    bool provided_memory() const { return misc.provided_memory; };
    std::string out() const { return misc.out; };
    int num_permutation() const { return misc.permutation; };
    int seed() const { return misc.seed; };
    int thread() const { return misc.thread; };
    int memory() const { return misc.memory; };


    // prset
    std::vector<std::string> bed() const { return prset.bed; };
    std::vector<std::string> feature() const { return prset.feature; };
    std::string gtf() const { return prset.gtf; };
    std::string msigdb() const { return prset.msigdb; };

    // prsice
    std::string missing_score() const { return prsice.missing_score; };
    SCORING get_scoring() const
    {
        if (prsice.missing_score == "no_mean_imputation")
            return SCORING::SET_ZERO;
        else if (prsice.missing_score == "center")
            return SCORING::CENTER;
        else
            return SCORING::MEAN_IMPUTE;
    }

    double lower() const { return prsice.lower; };
    double upper() const { return prsice.upper; };
    double inter() const { return prsice.inter; };
    bool no_regress() const { return prsice.no_regress; };
    bool full() const { return prsice.full; };
    bool fastscore() const { return prsice.fastscore; };
    double bar_upper() const
    {
        return prsice.barlevel.back();
    }; // we have sorted it
    int model() const { return prsice.model; };

    // prslice
    bool perform_prslice() const { return prslice.provided; };
    int prslice_size() const { return prslice.size; };

    // species
    int num_auto() const { return species.num_auto; };
    bool no_x() const { return species.no_x; };
    bool no_y() const { return species.no_y; };
    bool no_xy() const { return species.no_xy; };
    bool no_mt() const { return species.no_mt; };

    // target
    std::string target_name() const { return target.name; };
    std::string target_type() const { return target.type; };
    std::string pheno_file() const { return target.pheno_file; };
    std::string pheno_col(size_t index) const
    {
        return target.pheno_col.at(index);
    };
    std::string keep_sample_file() const
    {
        return target.keep_sample ? target.keep_file : "";
    };
    std::string remove_sample_file() const
    {
        return target.remove_sample ? target.remove_file : "";
    };
    std::vector<std::string> pheno_col() const { return target.pheno_col; };
    std::vector<bool> is_binary() const { return target.is_binary; };
    std::vector<double> prevalence() const { return target.prevalence; };
    bool has_pheno_col() const { return !target.pheno_col.empty(); };
    bool is_binary(size_t index) const { return target.is_binary.at(index); };
    bool nonfounders() const { return target.nonfounders; };
    bool keep_sample() const { return target.keep_sample; };
    bool remove_sample() const { return target.remove_sample; };

    int get_category(double p) const
    {
        for (size_t i = 0; i < prsice.barlevel.size(); ++i) {
            if (p <= prsice.barlevel[i]) {
                return i;
            }
        }
        if (p > prsice.barlevel.back()) return prsice.barlevel.size();
        return -2; // this is impossible because something will always either be
                   // bigger than the end or smaller than the front
    };
    double get_threshold(int i) const
    {
        return (i < 0) ? prsice.barlevel.at(0)
                       : ((i >= (int) prsice.barlevel.size())
                              ? 1.0
                              : prsice.barlevel.at(i));
    };

protected:
private:
    bool process(int argc, char* argv[], const char* optString,
                 const struct option longOpts[], Reporter& reporter);
    std::vector<std::string> supported_types = {"bed", "ped", "bgen"};
    struct
    {
        std::string name;
        std::string chr;
        std::string ref_allele;
        std::string alt_allele;
        std::string statistic;
        std::string snp;
        std::string bp;
        std::string standard_error;
        std::string p_value;
        std::string info_col;
        std::string maf;
        std::vector<int> col_index;
        int beta;
        int index;
        bool provided_maf;
        bool provided_chr;
        bool provided_ref;
        bool provided_alt;
        bool provided_stat;
        bool provided_snp;
        bool provided_bp;
        bool provided_se;
        bool provided_p;
        bool use_info;
        double info_score;
        double maf_threshold;
    } base;

    struct
    {
        std::string ld;
        std::string type;
        std::string keep_file;
        std::string remove_file;
        int no_clump;
        double proxy;
        double p_value;
        double r2;
        int distance;
        bool provide_r2;
        bool provide_p;
        bool provide_distance;
        bool provide_proxy;
        bool keep_sample;
        bool remove_sample;
        bool use_type;
    } clumping;

    struct
    {
        std::string name;
        std::string ancestry_dim;
        std::vector<std::string> covariates;
    } covariate;

    struct
    {
        std::string exclude_file;
        std::string extract_file;
        double geno;
        double hard_threshold;
        double maf;
        double mind;
        double info_score;
        double ld_maf;
        double ld_geno;
        double ld_info;
        int hard_coding;
        int use_prob;
        int keep_ambig;
        bool extract;
        bool exclude;
        bool use_hard_thres;
        bool info_filtering;
        bool use_geno;
        bool use_maf;
        bool use_ld_maf;
        bool use_ld_geno;
        bool use_ld_info;
        bool use_mind;
    } filter;

    struct
    {
        std::string out;
        int all;
        int ignore_fid;
        int logit_perm;
        int memory;
        int permutation;
        int print_snp;
        int print_all_samples;
        int thread;
        size_t seed;
        bool provided_seed;
        bool provided_permutation;
        bool provide_thread;
        bool provided_memory;
        bool provided_output;
        // I want to include cross-validation here
        // do it after writing up the paper. A useful resource is here
        // https://stats.stackexchange.com/questions/103459/how-do-i-know-which-method-of-cross-validation-is-best
    } misc;

    struct
    {
        std::vector<std::string> bed;
        std::vector<std::string> feature;
        std::string gtf;
        std::string msigdb;
        bool perform_prset;
    } prset;

    struct
    {
        std::string missing_score;
        std::vector<double> barlevel;
        double lower;
        double upper;
        double inter;
        bool provide_lower;
        bool provide_upper;
        bool provide_inter;
        bool provided_model;
        int model; // Use MODEL enum
        int fastscore;
        int no_regress;
        int full;
    } prsice;

    struct
    {
        int size;
        bool provided;
    } prslice;

    struct
    {
        int num_auto;
        int no_x;
        int no_y;
        int no_xy;
        int no_mt;
        bool double_set;
    } species;

    struct
    {
        std::string name;
        std::string keep_file;
        std::string remove_file;
        std::string pheno_file;
        std::string type;
        std::vector<std::string> pheno_col;
        std::vector<double>
            prevalence; // should equal to number of binary target
        int nonfounders;
        bool keep_sample;
        bool remove_sample;
        bool use_type;
        std::vector<bool> is_binary;
    } target;

    std::string help_message;
    void usage();
    void info();
    void base_check(std::map<std::string, std::string>& message, bool& error,
                    std::string& error_message);
    void clump_check(std::map<std::string, std::string>& message, bool& error,
                     std::string& error_message);
    void covariate_check(bool& error, std::string& error_message);
    void filter_check(bool& error, std::string& error_message);
    void misc_check(std::map<std::string, std::string>& message, bool& error,
                    std::string& error_message);
    void prset_check(std::map<std::string, std::string>& message, bool& error,
                     std::string& error_message);
    void prslice_check(bool& error, std::string& error_message);
    void prsice_check(std::map<std::string, std::string>& message, bool& error,
                      std::string& error_message);
    void target_check(std::map<std::string, std::string>& message, bool& error,
                      std::string& error_message);

    inline void set_species(int num_auto, bool no_x, bool no_y, bool no_xy,
                            bool no_mt, bool& error, std::string& error_message,
                            bool& species_error)
    {
        if (species.double_set && !species_error) {
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

    inline void load_binary_vector(const std::string& input,
                                   std::map<std::string, std::string>& message,
                                   std::string& error_message,
                                   std::vector<bool>& target, bool& error,
                                   const std::string& c)
    {

        message[c] = input;
        std::vector<std::string> token = misc::split(input, ",");
        try
        {
            for (auto&& bin : token) target.push_back(misc::to_bool(bin));
        }
        catch (const std::runtime_error& er)
        {
            error_message.append("ERROR: Invalid argument passed to " + c + ": "
                                 + input + "!\n");
            error_message.append(
                "       Require binary arguments e.g. T/F, True/False\n");
        }
    }

    inline void load_string_vector(const std::string& input,
                                   std::map<std::string, std::string>& message,
                                   std::vector<std::string>& target,
                                   const std::string& c,
                                   std::string& error_message)
    {

        message[c] = input;
        std::vector<std::string> token = misc::split(input, ",");
        target.insert(target.end(), token.begin(), token.end());
    }

    template <typename T>
    inline void load_numeric_vector(const std::string& input,
                                    std::map<std::string, std::string>& message,
                                    std::string& error_message,
                                    std::vector<T>& target, bool& error,
                                    const std::string& c)
    {
        message[c] = input;
        std::vector<std::string> token = misc::split(optarg, ",");
        try
        {
            for (auto&& bar : token) target.push_back(misc::convert<T>(bar));
        }
        catch (const std::runtime_error& er)
        {
            error_message.append("ERROR: Non numeric argument passed to " + c
                                 + ": " + input + "!\n");
            error = true;
        }
    }

    template <typename Type>
    inline void set_numeric(const std::string& input,
                            std::map<std::string, std::string>& message,
                            std::string& error_message, Type& target,
                            bool& target_boolean, bool& error,
                            const std::string& c)
    {
        if (message.find(c) != message.end()) {
            error_message.append("Warning: Duplicated argument --" + c + "\n");
        }
        message[c] = input;
        try
        {
            target = misc::convert<Type>(input);
            target_boolean = true;
        }
        catch (const std::runtime_error& er)
        {
            error_message.append("ERROR: Non numeric argument passed to " + c
                                 + ": " + input + "!\n");
        }
    }

    inline void set_model(const std::string& in,
                          std::map<std::string, std::string>& message,
                          std::string& error_message, bool& error)
    {
        std::string input = in;
        if (input.empty()) {
            error_message.append("ERROR: Model cannot be empty!\n");
            error = true;
        }
        prsice.provided_model = true;
        std::transform(input.begin(), input.end(), input.begin(), ::toupper);
        if (input.at(0) == 'A') {
            input = "add";
            prsice.model = +MODEL::ADDITIVE;
        }
        else if (input.at(0) == 'D')
        {
            input = "dom";
            prsice.model = +MODEL::DOMINANT;
        }
        else if (input.at(0) == 'R')
        {
            input = "rec";
            prsice.model = +MODEL::RECESSIVE;
        }
        else if (input.at(0) == 'H')
        {
            input = "het";
            prsice.model = +MODEL::HETEROZYGOUS;
        }
        else
        {
            error = true;
            error_message.append("ERROR: Unrecognized model: " + input + "!\n");
        }
        if (message.find("model") != message.end()) {
            error_message.append("Warning: Duplicated argument --model\n");
        }
        message["model"] = input;
    }
    inline void set_string(const std::string& input,
                           std::map<std::string, std::string>& message,
                           std::string& target, bool& target_boolean,
                           const std::string& c, std::string& error_message)
    {

        if (message.find(c) != message.end()) {
            error_message.append("Warning: Duplicated argument --" + c + "\n");
        }
        message[c] = input;
        target = input;
        target_boolean = true;
    }

    inline int index_check(const std::string& target,
                           const std::vector<std::string>& ref) const
    {
        for (size_t i = 0; i < ref.size(); ++i) {
            if (target.compare(ref[i]) == 0) {
                return i;
            }
        }
        return -1;
    };

    inline int index_check(const std::string& target, const int max,
                           bool& error, std::string& error_message,
                           const std::string& name)
    {
        try
        {
            int index = misc::convert<int>(optarg);
            if (index >= max) {
                error = true;
                error_message.append("ERROR: " + name
                                     + " index out of bound!\n");
                return -1;
            }
            if (index < 0) {
                error = true;
                error_message.append("ERROR: Negative " + name + " index!\n");
                return -1;
            }
            return index;
        }
        catch (const std::runtime_error& er)
        {
            error = true;
            error_message.append("ERROR: " + name + " index is not numeric!\n");
            return -1;
        }
    }
};

#endif // COMMANDER_H
