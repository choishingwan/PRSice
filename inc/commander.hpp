// This file is part of PRSice-2, copyright (C) 2016-2019
// Shing Wan Choi, Paul F. O’Reilly
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

#include "enumerators.h"
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
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>
#ifdef _WIN32
#include <windows.h>
#endif

const std::string version = "2.2.13";
const std::string date = "2020-02-27";
class Commander
{
public:
    Commander();
    virtual ~Commander();
    bool init(int argc, char* argv[], Reporter& reporter);
    std::string get_base_name() const
    {
        return misc::remove_extension<std::string>(
            misc::base_name<std::string>(m_base_info.file_name));
    }
    const BaseFile& get_base() const { return m_base_info; }
    const CalculatePRS& get_prs_instruction() const { return m_prs_info; }
    const Clumping& get_clump_info() const { return m_clump_info; }
    const GeneSets& get_set() const { return m_prset; }
    const GenoFile& get_target() const { return m_target; }
    const GenoFile& get_reference() const { return m_reference; }
    const Permutations& get_perm() const { return m_perm_info; }
    const Phenotype& get_pheno() const { return m_pheno_info; }
    const PThresholding& get_p_threshold() const { return m_p_thresholds; }
    const QCFiltering& get_base_qc() const { return m_base_filter; }
    const QCFiltering& get_target_qc() const { return m_target_filter; }
    const QCFiltering& get_ref_qc() const { return m_ref_filter; }
    bool use_ref() const
    {
        return !(m_reference.file_list.empty()
                 && m_reference.file_name.empty());
    }
    bool need_ref() const
    {
        return !m_clump_info.no_clump || m_prs_info.use_ref_maf;
    }
    bool use_inter() const { return m_allow_inter; }
    std::string delim() const { return m_id_delim; }
    std::string out() const { return m_out_prefix; }
    std::string exclusion_range() const { return m_exclusion_range; }
    bool all_scores() const { return m_print_all_scores; }
    bool print_snp() const { return m_print_snp; }
    unsigned long long max_memory(const unsigned long long detected) const
    {
        return (!m_provided_memory || m_memory > detected) ? detected
                                                           : m_memory;
    }
    unsigned long long memory() const { return m_memory; }
    std::string exclude_file() const { return m_exclude_file; }
    std::string extract_file() const { return m_extract_file; }
    bool keep_ambig() const { return m_keep_ambig; }
    bool nonfounders() const { return m_include_nonfounders; }
    bool ultra_aggressive() const { return m_ultra_aggressive; }
    bool ambig_no_flip() const { return m_ambig_no_flip; }

protected:
    const std::vector<std::string> supported_types = {"bed", "ped", "bgen"};
    std::string m_id_delim = " ";
    std::string m_out_prefix = "PRSice";
    std::string m_exclusion_range = "";
    std::string m_exclude_file = "";
    std::string m_extract_file = "";
    std::string m_help_message;
    size_t m_memory = 1e10;
    int m_allow_inter = false;
    int m_ambig_no_flip = false;
    int m_include_nonfounders = false;
    int m_keep_ambig = false;
    int m_print_all_scores = false;
    int m_print_snp = false;
    int m_ultra_aggressive = false;
    int m_user_no_default = false;
    bool m_provided_memory = false;
    bool m_set_delim = false;

    BaseFile m_base_info;
    CalculatePRS m_prs_info;
    Clumping m_clump_info;
    GeneSets m_prset;
    GenoFile m_target;
    GenoFile m_reference;
    Permutations m_perm_info;
    QCFiltering m_target_filter;
    QCFiltering m_ref_filter;
    QCFiltering m_base_filter;
    Phenotype m_pheno_info;
    PThresholding m_p_thresholds;

    std::map<std::string, std::string> m_parameter_log;
    std::string m_error_message = "";
    ////////////////////////////////////////////
    //
    // end of struct definition
    //
    ////////////////////////////////////////////

    /*!
     * \brief Parsing all command line arguments
     * \param argc
     * \param argv
     * \param optString is the string containing the short flags
     * \param longOpts is the struct containing all the options
     * \param reporter is the object to report all messages
     * \return true if we want to continue the program
     */
    bool parse_command(int argc, char* argv[], const char* optString,
                       const struct option longOpts[], Reporter& reporter);
    void usage();
    void set_help_message();

    bool clump_check();
    bool ref_check();
    bool covariate_check();
    bool filter_check();
    bool misc_check();
    bool prset_check();
    bool prsice_check();
    bool target_check();
    bool pheno_check();


    static std::vector<std::string>
    transform_covariate(const std::string& cov_in);


    /////////////////////////////////////////////////
    /// REFACTORED FUNCTIONS
    /////////////////////////////////////////////////

    int32_t maximum_thread()
    {
        int32_t max_threads = 1;
#if defined(WIN32) || defined(_WIN32) \
    || defined(__WIN32) && !defined(__CYGWIN__)
        // max thread estimation using windows
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        max_threads = sysinfo.dwNumberOfProcessors;
        int32_t known_procs = max_threads;
#else
        int32_t known_procs =
            static_cast<int32_t>(sysconf(_SC_NPROCESSORS_ONLN));
        max_threads = (known_procs == -1) ? 1 : known_procs;
#endif
        return max_threads;
    }
    inline void set_string(const std::string& input, const std::string& c,
                           size_t base_index)
    {
        set_string<int>(input, c, m_base_info.column_name[base_index],
                        m_base_info.has_column[base_index]);
    }
    inline void set_string(const std::string& input, const std::string& c,
                           std::string& target)
    {
        bool dumy;
        set_string<bool>(input, c, target, dumy);
    }
    template <typename T>
    inline void set_string(const std::string& input, const std::string& c,
                           std::string& target, T& target_boolean,
                           bool add_quote = false)
    {
        check_duplicate(c);
        std::string in = input;
        if (add_quote) { in = "\"" + input + "\""; }
        m_parameter_log[c] = in;
        target = input;
        target_boolean = true;
    }
    template <typename T>
    inline bool convert_to_numeric_vector(const std::vector<std::string>& token,
                                          std::vector<T>& target)
    {
        try
        {
            for (auto&& bar : token) target.push_back(misc::convert<T>(bar));
        }
        catch (...)
        {
            return false;
        }
        return true;
    }
    template <typename T>
    inline bool load_numeric_vector(const std::string& input,
                                    const std::string& c,
                                    std::vector<T>& target)
    {
        // should always have an input
        if (input.empty()) return false;
        std::string comma = "";
        if (m_parameter_log.find(c) != m_parameter_log.end()) { comma = ","; }
        m_parameter_log[c].append(comma + input);
        if (!input.empty() && input.back() == ',')
        {
            m_error_message.append(
                "Warning: , detected at end of input: " + input
                + ". Have you accidentally included space in "
                  "your input? (Space is not allowed)\n");
        }
        std::vector<std::string> token = misc::split(input, ",");
        if (!convert_to_numeric_vector<T>(token, target))
        {
            m_error_message.append("Error: Non numeric argument passed to " + c
                                   + ": " + input + "!\n");
            return false;
        }
        return true;
    }

    inline bool check_true_false_ending(const std::string& in, size_t& length)
    {
        length = 1;
        if (misc::hasEnding(in, "F")) { return false; }
        else if (misc::hasEnding(in, "T"))
        {
            return true;
        }
        else if (misc::hasEnding(in, "TRUE"))
        {
            length = 4;
            return true;
        }
        else if (misc::hasEnding(in, "FALSE"))
        {
            length = 5;
            return false;
        }
        throw std::runtime_error("Error: Undefined input");
    }
    size_t number_boolean(const std::string& input, bool& result)
    {
        size_t bool_length = 0;
        try
        {
            result = check_true_false_ending(input, bool_length);
        }
        catch (const std::runtime_error& er)
        {
            throw er;
        }
        try
        {
            if (bool_length != input.length())
            {
                return misc::string_to_size_t(
                    input.substr(0, input.length() - bool_length).c_str());
            }
            else
            {
                return 1;
            }
        }
        catch (const std::runtime_error&)
        {
            throw std::runtime_error("Error: None Numeric Pattern");
        }
    }
    inline bool parse_binary_vector(const std::string& input,
                                    const std::string& c,
                                    std::vector<bool>& target)
    {
        if (input.empty()) return false;
        std::string comma = "";
        if (m_parameter_log.find(c) != m_parameter_log.end()) { comma = ","; }
        m_parameter_log[c].append(comma + input);
        if (!input.empty() && input.back() == ',')
        {
            m_error_message.append(
                "Warning: , detected at end of input: " + input
                + ". Have you accidentally included space in "
                  "your input? (Space is not allowed)");
        }
        std::vector<std::string> token = misc::split(input, ",");
        try
        {
            for (auto&& bin : token)
            {
                // check if this is true or false, if, not, try parsing
                std::transform(bin.begin(), bin.end(), bin.begin(), ::toupper);
                try
                {
                    bool value = false;
                    size_t num_repeat = number_boolean(bin, value);
                    for (size_t i = 0; i < num_repeat; ++i)
                    { target.push_back(value); }
                }
                catch (const std::runtime_error&)
                {
                    m_error_message.append(
                        "Error: Invalid argument passed to " + c + ": " + input
                        + "! Require binary arguments e.g. T/F, "
                          "True/False, or numerical combination e.g. "
                          "10T,3F etc\n");
                    return false;
                }
            }
        }
        catch (...)
        {
            m_error_message.append(
                "Error: Invalid argument passed to " + c + ": " + input
                + "! Require binary arguments e.g. T/F, True/False\n");
            return false;
        }
        return true;
    }

    inline bool extract_unit(const std::string& input, double& value,
                             std::string& unit)
    {
        bool valid = false;
        size_t unit_length = 0;
        for (; unit_length <= 2; ++unit_length)
        {
            try
            {
                value = misc::convert<double>(
                    input.substr(0, input.length() - unit_length));
                valid = true;
                break;
            }
            catch (const std::runtime_error&)
            {
            }
        }
        if (!valid) return false;
        if (unit_length == 0) { unit = "b"; }
        else
        {
            unit = input.substr(input.length() - unit_length);
        }
        return true;
    }
    inline size_t unit_power(const std::string& unit)
    {
        const std::unordered_map<std::string, size_t> unit_map = {
            {"b", 0},  {"bp", 0}, {"k", 1},  {"kb", 1}, {"m", 2},
            {"mb", 2}, {"g", 3},  {"gb", 3}, {"t", 4},  {"tb", 4}};
        return unit_map.at(unit);
    }

    inline bool parse_unit_value(const std::string& input, const std::string& c,
                                 const size_t default_power, size_t& target,
                                 bool memory = false)
    {
        check_duplicate(c);
        std::string in = input;
        std::transform(in.begin(), in.end(), in.begin(), ::tolower);
        m_parameter_log[c] = in;
        const size_t weight = memory ? 1024 : 1000;
        double value;
        std::string unit;
        if (!extract_unit(in, value, unit))
        {
            m_error_message.append("Error: Invalid input: " + in + "\n");
            return false;
        }
        size_t unit_power_level;
        try
        {
            unit_power_level = unit_power(unit) + default_power;
        }
        catch (...)
        {
            m_error_message.append("Error: Invalid input: " + in + "\n");
            return false;
        }
        value *= pow(weight, unit_power_level);
        if (trunc(value) != value && value < 0)
        {
            m_error_message.append("Error: Non-integer value obtained: "
                                   + misc::to_string(target) + "\n");
            return false;
        }
        target = static_cast<size_t>(value);
        return true;
    }


    void check_duplicate(const std::string& c)
    {
        if (m_parameter_log.find(c) != m_parameter_log.end())
        {
            m_error_message.append("Warning: Duplicated argument --" + c
                                   + "\n");
        }
    }

    template <typename Type>
    inline bool set_numeric(const std::string& input, const std::string& c,
                            Type& target)
    {
        bool dumy;
        return set_numeric<Type>(input, c, target, dumy);
    }
    template <typename Type>
    inline bool set_numeric(const std::string& input, const std::string& c,
                            Type& target, bool& target_boolean)
    {
        check_duplicate(c);
        m_parameter_log[c] = input;
        try
        {
            target = misc::convert<Type>(input);
            target_boolean = true;
            return true;
        }
        catch (...)
        {
            m_error_message.append("Error: Invalid numeric argument passed to "
                                   + c + ": " + input + "!\n");
            return false;
        }
    }

    inline void load_string_vector(const std::string& input,
                                   const std::string& c,
                                   std::vector<std::string>& target)
    {
        if (input.empty()) return;
        if (m_parameter_log.find(c) == m_parameter_log.end())
        { m_parameter_log[c] = input; }
        else
        {
            m_parameter_log[c] = "," + input;
        }
        if (!input.empty() && input.back() == ',')
        {
            m_error_message.append(
                "Warning: , detected at end of input: " + input
                + ". Have you accidentally included space in "
                  "your input? (Space is not allowed)\n");
        }
        std::vector<std::string> token = misc::split(input, ",");
        target.insert(target.end(), token.begin(), token.end());
    }

    inline bool set_memory(const std::string& input)
    {
        bool error = parse_unit_value(input, "memory", 2, m_memory, true);
        return !error;
    }

    inline bool set_missing(const std::string& in)
    {
        std::string input = in;
        check_duplicate("missing");
        std::transform(input.begin(), input.end(), input.begin(), ::tolower);
        switch (input.at(0))
        {
        case 'c':
            input = "centre";
            m_prs_info.missing_score = MISSING_SCORE::CENTER;
            break;
        case 'm':
            input = "mean_impute";
            m_prs_info.missing_score = MISSING_SCORE::MEAN_IMPUTE;
            break;
        case 's':
            input = "set_zero";
            m_prs_info.missing_score = MISSING_SCORE::SET_ZERO;
            break;
        case 'i':
            input = "impute_control";
            m_prs_info.missing_score = MISSING_SCORE::IMPUTE_CONTROL;
            break;
        default:
            m_error_message.append(
                "Error: Unrecognized Missing handling method: " + in + "!\n");
            return false;
        }
        m_parameter_log["missing"] = input;
        return true;
    }
    inline bool set_model(const std::string& in)
    {
        std::string input = in;
        std::transform(input.begin(), input.end(), input.begin(), ::tolower);
        check_duplicate("model");
        switch (input.at(0))
        {
        case 'a':
            input = "add";
            m_prs_info.genetic_model = MODEL::ADDITIVE;
            break;
        case 'd':
            input = "dom";
            m_prs_info.genetic_model = MODEL::DOMINANT;
            break;
        case 'r':
            input = "rec";
            m_prs_info.genetic_model = MODEL::RECESSIVE;
            break;
        case 'h':
            input = "het";
            m_prs_info.genetic_model = MODEL::HETEROZYGOUS;
            break;
        default:
            m_error_message.append("Error: Unrecognized model: " + in + "!\n");
            return false;
        }
        m_parameter_log["model"] = input;
        return true;
    }
    inline bool set_score(const std::string& in)
    {
        std::string input = in;
        std::transform(input.begin(), input.end(), input.begin(), ::tolower);
        check_duplicate("score");
        if (input == "avg") { m_prs_info.scoring_method = SCORING::AVERAGE; }
        else if (input == "std")
        {
            m_prs_info.scoring_method = SCORING::STANDARDIZE;
        }
        else if (input == "sum")
        {
            m_prs_info.scoring_method = SCORING::SUM;
        }
        else if (input == "con-std")
        {
            m_prs_info.scoring_method = SCORING::CONTROL_STD;
        }
        else
        {
            m_error_message.append("Error: Unrecognized scoring method: " + in
                                   + "!\n");
            return false;
        }
        m_parameter_log["score"] = input;
        return true;
    }

    bool in_file(const std::vector<std::string>& column_names,
                 const size_t index, const std::string& warning,
                 bool no_default, bool case_sensitive = true,
                 bool print_error = true)
    {
        if ((no_default && !static_cast<bool>(m_base_info.has_column[index]))
            || m_base_info.column_name[index].empty())
        {
            m_base_info.has_column[index] = false;
            return false;
        }
        size_t col_index;
        bool has_col = index_check(m_base_info.column_name[index], column_names,
                                   col_index, case_sensitive);
        if (has_col) { m_base_info.column_index[index] = col_index; }
        else if (m_base_info.has_column[index] && print_error)
        {
            m_error_message.append(warning + ": "
                                   + m_base_info.column_name[index]
                                   + " not found in base file\n");
        }
        m_base_info.has_column[index] = has_col;
        return has_col;
    }
    /**
     * @return True if we don't need to error out. False otherwise
     */
    inline bool set_base_info_threshold(const std::vector<std::string>& ref)
    {
        const std::vector<std::string> info =
            misc::split(m_base_info.column_name[+BASE_INDEX::INFO], ",");
        const bool has_input = m_base_info.has_column[+BASE_INDEX::INFO];

        size_t index;
        const bool found = index_check(info[0], ref, index);
        if (found) m_base_info.column_index[+BASE_INDEX::INFO] = index;
        m_base_info.has_column[+BASE_INDEX::INFO] = found;
        if (!found)
        {
            if (has_input)
            {
                m_error_message.append("Warning: INFO field not found in base "
                                       "file, will ignore INFO filtering\n");
            }
        }
        else if (info.size() != 2) // assume default always valid
        {
            m_error_message.append("Error: Invalid format of "
                                   "--base-info. Should be "
                                   "ColName,Threshold.\n");
            return false;
        }
        try
        {
            m_base_filter.info_score = misc::convert<double>(info[1]);
            if (!misc::within_bound<double>(m_base_filter.info_score, 0.0, 1.0))
            {
                if (has_input)
                {
                    m_error_message.append("Error: Base INFO threshold "
                                           "must be within 0 and 1!\n");
                    return false;
                }
            }
        }
        catch (...)
        {
            if (has_input)
            {
                m_error_message.append(
                    "Error: Invalid argument passed to --base-info: "
                    + m_base_info.column_name[+BASE_INDEX::INFO]
                    + "! Second argument must be numeric\n");
                return false;
            }
        }
        return true;
    }
    /**
     * @return True if we don't need to error out. False otherwise
     */
    bool process_maf(const std::vector<std::string>& ref,
                     const std::vector<std::string>& detail,
                     size_t& column_index, int& has_column, double& maf)
    {
        size_t index = 0;
        bool found = index_check(detail[0], ref, index);
        has_column = found;
        if (found) column_index = index;
        if (!found)
        {
            m_error_message.append(
                "Warning: MAF field not found in base file. "
                "Will not perform MAF filtering on the base file\n");
            return true;
        }
        double cur_maf;
        try
        {
            cur_maf = misc::convert<double>(detail[1]);
            if (!misc::within_bound<double>(cur_maf, 0.0, 1.0))
            {
                m_error_message.append("Error: Base MAF threshold must "
                                       "be within 0 and 1!\n");
                return false;
            }
        }
        catch (...)
        {
            m_error_message.append(
                "Error: Invalid argument passed to --base-maf: "
                + m_base_info.column_name[+BASE_INDEX::MAF]
                + "! Threshold must be numeric\n");
            return false;
        }
        maf = cur_maf;
        return true;
    }
    // return true if valid
    inline bool set_base_maf_filter(const std::vector<std::string>& ref)
    {
        const std::string maf_error =
            "Error: Invalid format of --base-maf. "
            "Should be ColName,Threshold."
            "or ColName,Threshold:ColName,Threshold.\n";
        std::vector<std::string> case_control =
            misc::split(m_base_info.column_name[+BASE_INDEX::MAF], ":");
        const bool print_error = m_base_info.has_column[+BASE_INDEX::MAF];
        // only process the maf filter if it is provided
        if (!print_error) return true;
        if (case_control.size() > 2)
        {
            if (print_error) { m_error_message.append(maf_error); }
            return false;
        }
        std::vector<std::string> detail;
        if (case_control.size() > 0)
        {
            detail = misc::split(case_control.front(), ",");
            return process_maf(
                ref, detail, m_base_info.column_index[+BASE_INDEX::MAF],
                m_base_info.has_column[+BASE_INDEX::MAF], m_base_filter.maf);
        }
        if (case_control.size() == 2)
        {
            detail = misc::split(case_control.back(), ",");
            return process_maf(ref, detail,
                               m_base_info.column_index[+BASE_INDEX::MAF_CASE],
                               m_base_info.has_column[+BASE_INDEX::MAF_CASE],
                               m_base_filter.maf_case);
        }
        return true;
    }
    /*!
     * \brief Get the column index based on file header and the input string
     * \param target is the input string
     * \param ref is the vector containing the column header
     * \return the index to the column containing the column name. -1 if it is
     * not found
     */
    inline bool index_check(const std::string& target,
                            const std::vector<std::string>& ref, size_t& index,
                            bool case_sensitive = true) const
    {
        std::string tmp = "";
        for (size_t i = 0; i < ref.size(); ++i)
        {
            tmp = ref[i];
            if (!case_sensitive)
            {
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            }
            if (target == tmp)
            {
                index = i;
                return true;
            }
        }
        return false;
    }

    bool get_statistic_column(const std::vector<std::string>& column_names);
    bool base_check();
    bool get_statistic_flag();
    std::string get_program_header(const std::string& name);
};

#endif // COMMANDER_H
