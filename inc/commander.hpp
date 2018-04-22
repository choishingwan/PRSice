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

#include "gzstream.h"
#include "misc.hpp"
#include "storage.hpp"
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <random>
#include <reporter.hpp>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <unordered_set>
#include <vector>
#include <zlib.h>
#ifdef _WIN32
#include <windows.h>
#endif
const std::string version = "2.1.1.beta";
const std::string date = "10 April 2018";
class Commander
{
public:
    Commander();
    virtual ~Commander();
    bool init(int argc, char* argv[], Reporter& reporter);

    // base
    std::vector<int> index() const { return base.col_index; };
    bool is_index() const { return base.is_index; };
    bool beta() const { return base.is_beta; };
    double base_info_score() const { return base.info_score_threshold; };
    double maf_base_control() const { return base.maf_control_threshold; };
    std::string base_name() const { return base.name; };
    double maf_base_case() const { return base.maf_case_threshold; };

    // clump
    bool no_clump() const { return clumping.no_clump; };
    bool use_proxy() const { return clumping.provided_proxy; };
    double proxy() const { return clumping.proxy; };
    double clump_p() const { return clumping.p_value; };
    double clump_r2() const { return clumping.r2; };
    double clump_dist() const { return clumping.distance * 1000; };

    // covariate
    std::string get_cov_file() const { return covariate.file_name; };
    std::vector<std::string> get_cov_header() const
    {
        return covariate.covariates;
    };
    // reference panel
    std::string ref_name() const { return reference_panel.file_name; };
    std::string ref_type() const { return reference_panel.type; };
    std::string ld_keep_file() const { return reference_panel.keep_file; };
    std::string ld_remove_file() const { return reference_panel.remove_file; };
    // reference filtering
    double ld_geno() const { return reference_snp_filtering.geno; };
    double ld_hard_threshold() const
    {
        return reference_snp_filtering.hard_threshold;
    };
    double ld_maf() const { return reference_snp_filtering.maf; };
    double ld_info() const { return reference_snp_filtering.info_score; };

    // misc
    std::string out() const { return misc.out; };
    bool all_scores() const { return misc.print_all_scores; };
    bool ignore_fid() const { return misc.ignore_fid; };
    bool logit_perm() const { return misc.logit_perm; };
    bool print_snp() const { return misc.print_snp; };
    bool pearson() const { return misc.pearson; };
    int permutation() const { return misc.permutation; };
    int seed() const { return misc.seed; };
    int thread() const { return misc.thread; };
    size_t max_memory(const size_t detected) const
    {
        if (!misc.provided_memory)
            return detected;
        else
            return (misc.memory > detected) ? detected : misc.memory;
    }
    // p_thresholds
    double bar_upper() const { return p_thresholds.barlevel.back(); };
    double get_threshold(int i) const
    {
        return (i < 0) ? p_thresholds.barlevel.at(0)
                       : ((i >= (int) p_thresholds.barlevel.size())
                              ? 1.0
                              : p_thresholds.barlevel.at(i));
    };
    double lower() const { return p_thresholds.lower; };
    double upper() const { return p_thresholds.upper; };
    double inter() const { return p_thresholds.inter; };

    int get_category(double p) const
    {
        for (size_t i = 0; i < p_thresholds.barlevel.size(); ++i) {
            if (p <= p_thresholds.barlevel[i]) {
                return i;
            }
        }
        if (p > p_thresholds.barlevel.back())
            return p_thresholds.barlevel.size();
        return -2; // this is impossible because something will always either be
                   // bigger than the end or smaller than the front
    };
    bool no_full() const { return p_thresholds.no_full; };
    bool fastscore() const { return p_thresholds.fastscore; };

    // prs_calculation
    MISSING_SCORE get_missing_score() const
    {
        std::string s = prs_calculation.missing_score;
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        if (s == "NO_MEAN_IMPUTATION")
            return MISSING_SCORE::SET_ZERO;
        else if (s == "CENTER")
            return MISSING_SCORE::CENTER;
        else
            return MISSING_SCORE::MEAN_IMPUTE;
    }

    SCORING get_score() const
    {
        std::string s = prs_calculation.score_calculation;
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        if (s == "STD")
            return SCORING::STANDARDIZE;
        else if (s == "SUM")
            return SCORING::SUM;
        else
            return SCORING::AVERAGE;
    }

    MODEL model() const { return prs_calculation.model; };
    bool no_regress() const { return prs_calculation.no_regress; };

    // prs_snp_filtering
    std::string exclude_file() const { return prs_snp_filtering.exclude_file; };
    std::string extract_file() const { return prs_snp_filtering.extract_file; };
    double geno() const { return prs_snp_filtering.geno; };
    double hard_threshold() const { return prs_snp_filtering.hard_threshold; };
    double maf() const { return prs_snp_filtering.maf; };
    double info() const { return prs_snp_filtering.info_score; };
    bool hard_coded() const { return prs_snp_filtering.is_hard_coded; };
    bool keep_ambig() const { return prs_snp_filtering.keep_ambig; };

    // prset
    std::vector<std::string> bed() const { return prset.bed; };
    std::vector<std::string> feature() const { return prset.feature; };
    std::string gtf() const { return prset.gtf; };
    std::string msigdb() const { return prset.msigdb; };
    std::string background() const { return prset.background; };
    int set_perm() const { return prset.set_perm; };
    // prslice
    bool perform_prslice() const { return prslice.provided; };
    int prslice_size() const { return prslice.size; };
    int window_5() const { return prset.window_5; };
    int window_3() const { return prset.window_3; };
    // target
    std::string target_name() const { return target.name; };
    std::string target_type() const { return target.type; };
    std::string pheno_file() const { return target.pheno_file; };
    std::string pheno_col(size_t index) const
    {
        return target.pheno_col.at(index);
    };
    std::string keep_sample_file() const { return target.keep_file; };
    std::string remove_sample_file() const { return target.remove_file; };
    std::vector<std::string> pheno_col() const { return target.pheno_col; };
    std::vector<bool> is_binary() const { return target.is_binary; };
    std::vector<double> prevalence() const { return target.prevalence; };
    bool has_pheno_col() const { return !target.pheno_col.empty(); };
    bool is_binary(size_t index) const { return target.is_binary.at(index); };
    bool nonfounders() const { return target.include_nonfounders; };


protected:
private:
    bool process(int argc, char* argv[], const char* optString,
                 const struct option longOpts[], Reporter& reporter);
    std::vector<std::string> supported_types = {"bed", "ped", "bgen"};
    struct Base
    {
        std::string name;
        std::string chr;
        std::string effect_allele;
        std::string non_effect_allele;
        std::string statistic;
        std::string snp;
        std::string bp;
        std::string standard_error;
        std::string p_value;
        std::string info_col;
        std::string maf_col;
        std::vector<int> col_index;
        int is_beta;
        int is_index;
        int no_default;
        double info_score_threshold;
        double maf_control_threshold;
        double maf_case_threshold;
        // determine if we are going to use default
        bool provided_chr;
        bool provided_effect_allele;
        bool provided_non_effect_allele;
        bool provided_statistic;
        bool provided_snp;
        bool provided_bp;
        bool provided_standard_error;
        bool provided_p_value;
        bool provided_info;
    } base;

    struct Clump
    {
        int no_clump;
        double proxy;
        double p_value;
        double r2;
        int distance;
        bool provided_proxy;
    } clumping;

    struct Covar
    {
        std::string file_name;
        std::vector<std::string> covariates;
        // Numeric factors should be defined with ""
    } covariate;

    struct Misc
    {
        std::string out;
        int print_all_scores;
        int ignore_fid;
        int logit_perm;
        int pearson;
        int permutation;
        int print_snp;
        int thread;
        size_t memory;
        size_t seed;
        bool provided_seed;
        bool provided_memory;
    } misc;

    struct Reference
    {
        std::string file_name;
        std::string type;
        std::string keep_file;
        std::string remove_file;
    } reference_panel;

    struct Ref_filtering
    {
        double geno;
        double hard_threshold;
        double maf;
        double info_score;
        // Need to have the same snp as the target
        // so not necessary to have another set of exclude/extract files
        // must be hard coded for reference
    } reference_snp_filtering;

    struct Thresholding
    {
        std::vector<double> barlevel;
        double lower;
        double inter;
        double upper;
        int fastscore;
        int no_full;
        bool set_use_thresholds;
    } p_thresholds;

    struct Calculation
    {
        std::string missing_score;
        std::string model_name;
        std::string score_calculation;
        MODEL model; // use model enum
        int no_regress;
    } prs_calculation;

    struct Filtering
    {
        std::string exclude_file;
        std::string extract_file;
        double geno;
        double hard_threshold;
        double maf;
        double info_score;
        int is_hard_coded;
        int keep_ambig;
        int predict_ambig; // PRSoS stuff?
    } prs_snp_filtering;

    struct PRSet
    {
        std::vector<std::string> bed;
        std::vector<std::string> feature;
        std::string gtf;
        std::string msigdb;
        std::string background;
        int set_perm;
        int window_5;
        int window_3;
        bool perform_prset;
    } prset;

    struct PRSlice
    {
        int size;
        bool provided;
    } prslice;

    struct Target
    {
        std::string name;
        std::string keep_file;
        std::string remove_file;
        std::string pheno_file;
        std::string type;
        std::vector<std::string> pheno_col;
        // should equal to number of binary target
        std::vector<double> prevalence;
        std::vector<bool> is_binary;
        int include_nonfounders;
    } target;

    ////////////////////////////////////////////
    //
    // end of struct definition
    //
    ////////////////////////////////////////////


    std::string help_message;
    void usage();
    void set_help_message();
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

    inline void load_binary_vector(const std::string& input,
                                   std::map<std::string, std::string>& message,
                                   std::string& error_message,
                                   std::vector<bool>& target, bool& error,
                                   const std::string& c)
    {
        if (input.empty()) return;
        message[c] = message[c] + input;
        // check if the input is ended with , which is usually the case
        // when someone mixed in the space
        if (!input.empty() && input.back() == ',') {
            error_message.append("Warning: , detected at end of input: " + input
                                 + ". Have you accidentally included space in "
                                   "your input? (Space is not allowed)");
        }
        std::vector<std::string> token = misc::split(input, ",");
        try
        {
            for (auto&& bin : token) target.push_back(misc::to_bool(bin));
        }
        catch (const std::runtime_error& er)
        {
            error = true;
            error_message.append(
                "Error: Invalid argument passed to " + c + ": " + input
                + "! Require binary arguments e.g. T/F, True/False\n");
        }
    }

    inline void load_string_vector(const std::string& input,
                                   std::map<std::string, std::string>& message,
                                   std::vector<std::string>& target,
                                   const std::string& c,
                                   std::string& error_message)
    {

        if (input.empty()) return;
        message[c] = message[c] + input;
        if (!input.empty() && input.back() == ',') {
            error_message.append("Warning: , detected at end of input: " + input
                                 + ". Have you accidentally included space in "
                                   "your input? (Space is not allowed)\n");
        }
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
        if (input.empty()) return;
        message[c] = message[c] + input;
        if (!input.empty() && input.back() == ',') {
            error_message.append("Warning: , detected at end of input: " + input
                                 + ". Have you accidentally included space in "
                                   "your input? (Space is not allowed)\n");
        }
        std::vector<std::string> token = misc::split(optarg, ",");
        try
        {
            for (auto&& bar : token) target.push_back(misc::convert<T>(bar));
        }
        catch (const std::runtime_error& er)
        {
            error_message.append("Error: Non numeric argument passed to " + c
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
            error = true;
            error_message.append("Error: Non numeric argument passed to " + c
                                 + ": " + input + "!\n");
        }
    }

    inline void set_model(const std::string& in,
                          std::map<std::string, std::string>& message,
                          std::string& error_message, bool& error)
    {
        std::string input = in;
        if (input.empty()) {
            error_message.append("Error: Model cannot be empty!\n");
            error = true;
        }
        std::transform(input.begin(), input.end(), input.begin(), ::toupper);
        if (input.at(0) == 'A') {
            input = "add";
            prs_calculation.model = MODEL::ADDITIVE;
        }
        else if (input.at(0) == 'D')
        {
            input = "dom";
            prs_calculation.model = MODEL::DOMINANT;
        }
        else if (input.at(0) == 'R')
        {
            input = "rec";
            prs_calculation.model = MODEL::RECESSIVE;
        }
        else if (input.at(0) == 'H')
        {
            input = "het";
            prs_calculation.model = MODEL::HETEROZYGOUS;
        }
        else
        {
            error = true;
            error_message.append("Error: Unrecognized model: " + input + "!\n");
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

    inline void set_memory(const std::string& input,
                           std::map<std::string, std::string>& message,
                           std::string& error_messages, bool& error)
    {
        misc.provided_memory = true;
        if (message.find("memory") != message.end()) {
            error_messages.append("Warning: Duplicated argument --memory\n");
        }
        try
        {
            int memory = misc::convert<int>(input);
            misc.memory = memory;
        }
        catch (const std::runtime_error& er)
        {
            // contain MB KB or B here
            if (input.length() >= 2) {
                try
                {
                    std::string unit = input.substr(input.length() - 2);
                    std::string value = input.substr(0, input.length() - 2);
                    std::transform(unit.begin(), unit.end(), unit.begin(),
                                   ::toupper);
                    if (unit.compare("KB") == 0) {
                        misc.memory = misc::convert<size_t>(value) * 1024;
                    }
                    else if (unit.compare("MB") == 0)
                    {
                        misc.memory =
                            misc::convert<size_t>(value) * 1024 * 1024;
                    }
                    else if (unit.compare("GB") == 0)
                    {
                        misc.memory =
                            misc::convert<size_t>(value) * 1024 * 1024 * 1024;
                    }
                    else if (unit.compare("TB") == 0)
                    {
                        misc.memory = misc::convert<size_t>(value) * 1024 * 1024
                                      * 1024 * 1024;
                    }
                    else
                    {
                        // maybe only one input?
                        unit = input.substr(input.length() - 1);
                        value = input.substr(0, input.length() - 1);
                        std::transform(unit.begin(), unit.end(), unit.begin(),
                                       ::toupper);
                        if (unit.compare("B") == 0) {
                            misc.memory = misc::convert<size_t>(value);
                            std::cerr << "Update to: " << misc.memory
                                      << std::endl;
                        }
                    }
                }
                catch (const std::runtime_error& er_info)
                {
                    error = true;
                    error_messages.append("Error: Undefined memory input: "
                                          + input);
                }
            }
            else
            {
                error = true;
                error_messages.append("Error: Undefined memory input: "
                                      + input);
            }
        }

        message["memory"] = input;
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
                error_message.append("Error: " + name
                                     + " index out of bound!\n");
                return -1;
            }
            if (index < 0) {
                error = true;
                error_message.append("Error: Negative " + name + " index!\n");
                return -1;
            }
            return index;
        }
        catch (const std::runtime_error& er)
        {
            error = true;
            error_message.append("Error: " + name + " index is not numeric!\n");
            return -1;
        }
    }
};

#endif // COMMANDER_H
